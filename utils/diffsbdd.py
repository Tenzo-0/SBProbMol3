import subprocess
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
import xml.etree.ElementTree as ET


class AutoDockGPUPipelineError(RuntimeError):
    pass


@dataclass
class PreparedReceptor:
    protein_pdb: Path
    out_prefix: str
    receptor_pdbqt: Path
    receptor_gpf: Path
    maps_fld: Path


def run_command(cmd: List[str], cwd: Optional[Path] = None) -> str:
    """Run a shell command and return stdout, raise if non-zero exit."""
    print(f"[CMD] {' '.join(shlex.quote(c) for c in cmd)}")
    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    print(result.stdout)
    if result.returncode != 0:
        raise AutoDockGPUPipelineError(
            f"Command failed with exit code {result.returncode}: {' '.join(cmd)}\n"
            f"Output:\n{result.stdout}"
        )
    return result.stdout


def prepare_receptor(
    protein_pdb: Path,
    ligand_sdf: Path,
    out_prefix: str,
    padding: float = 6.0,
    workdir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """
    Run mk_prepare_receptor.py to produce receptor PDBQT and GPF.

    mk_prepare_receptor.py -i <protein>.pdb -o <name> -p -v -g \
        --box_enveloping <ligand>.sdf --padding 6
    """
    workdir = workdir or protein_pdb.parent
    receptor_pdbqt = workdir / f"{out_prefix}.pdbqt"
    receptor_gpf = workdir / f"{out_prefix}.gpf"

    cmd = [
        "mk_prepare_receptor.py",
        "-i", str(protein_pdb),
        "-o", out_prefix,
        "-p",
        "-v",
        "-g",
        "--box_enveloping", str(ligand_sdf),
        "--padding", str(padding),
    ]
    run_command(cmd, cwd=workdir)

    if not receptor_pdbqt.exists():
        raise AutoDockGPUPipelineError(f"Expected receptor PDBQT not found: {receptor_pdbqt}")
    if not receptor_gpf.exists():
        raise AutoDockGPUPipelineError(f"Expected receptor GPF not found: {receptor_gpf}")
    return receptor_pdbqt, receptor_gpf


def prepare_ligand_with_obabel(
    ligand_sdf: Path,
    ligand_pdbqt: Optional[Path] = None,
) -> Path:
    """
    Convert SDF -> PDBQT using Open Babel:

        obabel <ligand>.sdf -O <ligand>.pdbqt
    """
    if ligand_pdbqt is None:
        ligand_pdbqt = ligand_sdf.with_suffix(".pdbqt")
    cmd = [
        "obabel",
        str(ligand_sdf),
        "-O", str(ligand_pdbqt),
    ]
    run_command(cmd, cwd=ligand_sdf.parent)

    if not ligand_pdbqt.exists():
        raise AutoDockGPUPipelineError(f"Expected ligand PDBQT not found: {ligand_pdbqt}")
    return ligand_pdbqt


def run_autogrid(
    receptor_gpf: Path,
    workdir: Optional[Path] = None,
) -> Path:
    """
    Run AutoGrid:

        autogrid4 -p <protein>.gpf -l <protein>.glg

    Returns path to <protein>.maps.fld
    """
    workdir = workdir or receptor_gpf.parent
    prefix = receptor_gpf.with_suffix("")  # remove .gpf
    fld = Path(str(prefix) + ".maps.fld")

    cmd = [
        "autogrid4",
        "-p", str(receptor_gpf),
        "-l", str(prefix) + ".glg",
    ]
    run_command(cmd, cwd=workdir)

    if not fld.exists():
        raise AutoDockGPUPipelineError(f"Expected maps.fld not found: {fld}")
    return fld


def run_autodock_gpu(
    maps_fld: Path,
    ligand_pdbqt: Path,
    nrun: int = 50,
) -> Path:
    """
    Run docking with AutoDock-GPU (NO -O; expect XML):

        autodock_gpu_128wi --ffile <protein>.maps.fld \
                           --lfile <ligand>.pdbqt \
                           --nrun 50

    AutoDock-GPU will create an XML file, typically <ligand_stem>.xml
    """
    workdir = maps_fld.parent
    cmd = [
        "autodock_gpu_128wi",
        "--ffile", str(maps_fld),
        "--lfile", str(ligand_pdbqt),
        "--nrun", str(nrun),
    ]
    run_command(cmd, cwd=workdir)

    # Assume XML has same stem as ligand PDBQT: e.g. 3rfm_B_CFF.pdbqt -> 3rfm_B_CFF.xml
    xml_file = ligand_pdbqt.with_suffix(".xml")
    if not xml_file.exists():
        raise AutoDockGPUPipelineError(f"Expected XML output not found: {xml_file}")
    return xml_file


def parse_autodock_xml(xml_file: Path) -> Dict[str, Any]:
    """
    Parse AutoDock-GPU XML to extract energies.

    - Reads all <free_NRG_binding> under <runs><run>...</run>
    - Reads best cluster (cluster_rank="1") from <clustering_histogram>
      to get:
        lowest_binding_energy (our 'vina-like' score)
        mean_binding_energy

    Returns:
        {
          "run_energies": [float, ...],
          "best_cluster_lowest": float or None,
          "best_cluster_mean": float or None,
        }
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # 1) Energies for each run
    run_energies: List[float] = []
    runs_elem = root.find("runs")
    if runs_elem is not None:
        for run in runs_elem.findall("run"):
            e_elem = run.find("free_NRG_binding")
            if e_elem is not None and e_elem.text is not None:
                try:
                    run_energies.append(float(e_elem.text))
                except ValueError:
                    pass

    # 2) Best cluster (cluster_rank="1") from clustering_histogram
    best_lowest = None
    best_mean = None
    result_elem = root.find("result")
    if result_elem is not None:
        hist_elem = result_elem.find("clustering_histogram")
        if hist_elem is not None:
            for clus in hist_elem.findall("cluster"):
                if clus.attrib.get("cluster_rank") == "1":
                    lb = clus.attrib.get("lowest_binding_energy")
                    mb = clus.attrib.get("mean_binding_energy")
                    if lb is not None:
                        try:
                            best_lowest = float(lb)
                        except ValueError:
                            pass
                    if mb is not None:
                        try:
                            best_mean = float(mb)
                        except ValueError:
                            pass
                    break

    return {
        "run_energies": run_energies,
        "best_cluster_lowest": best_lowest,
        "best_cluster_mean": best_mean,
    }


def prepare_receptor_maps(
    protein_pdb: str,
    ligand_sdf: str,
    out_prefix: Optional[str] = None,
    padding: float = 6.0,
) -> PreparedReceptor:
    """Prepare receptor once (PDBQT + GPF + maps) using the provided ligand box."""
    protein_pdb_path = Path(protein_pdb).resolve()
    ligand_sdf_path = Path(ligand_sdf).resolve()
    prefix = out_prefix or protein_pdb_path.stem
    receptor_pdbqt, receptor_gpf = prepare_receptor(
        protein_pdb=protein_pdb_path,
        ligand_sdf=ligand_sdf_path,
        out_prefix=prefix,
        padding=padding,
    )
    maps_fld = run_autogrid(receptor_gpf)
    return PreparedReceptor(
        protein_pdb=protein_pdb_path,
        out_prefix=prefix,
        receptor_pdbqt=receptor_pdbqt,
        receptor_gpf=receptor_gpf,
        maps_fld=maps_fld,
    )


def score_ligand_with_prepared(
    prepared: PreparedReceptor,
    ligand_sdf: str,
    nrun: int = 50,
) -> Tuple[Optional[float], Path, Path]:
    """Score a ligand using a precomputed receptor/mapping."""
    ligand_sdf_path = Path(ligand_sdf).resolve()
    ligand_pdbqt = prepare_ligand_with_obabel(ligand_sdf_path)
    xml_file = run_autodock_gpu(
        maps_fld=prepared.maps_fld,
        ligand_pdbqt=ligand_pdbqt,
        nrun=nrun,
    )
    parsed = parse_autodock_xml(xml_file)
    best_score = parsed["best_cluster_lowest"]
    if best_score is None and parsed["run_energies"]:
        best_score = min(parsed["run_energies"])
    return best_score, xml_file, ligand_pdbqt


def run_full_pipeline(
    protein_pdb: str,
    ligand_sdf: str,
    out_prefix: Optional[str] = None,
    padding: float = 6.0,
    nrun: int = 50,
) -> Tuple[Optional[float], Path]:
    """
    Full pipeline:

      1) mk_prepare_receptor.py -i <protein>.pdb -o <name> -p -v -g \
            --box_enveloping <ligand>.sdf --padding 6
      2) obabel <ligand>.sdf -O <ligand>.pdbqt
      3) autogrid4 -p <protein>.gpf -l <protein>.glg
      4) autodock_gpu_128wi --ffile <protein>.maps.fld \
                            --lfile <ligand>.pdbqt \
                            --nrun <nrun>

    Then parse XML and return 'vina-like' best score:

      - Prefer best_cluster_lowest (cluster_rank=1)
      - Fallback to min(run_energies) if needed

    Returns:
        (best_score, xml_file)
    """
    prepared = prepare_receptor_maps(
        protein_pdb=protein_pdb,
        ligand_sdf=ligand_sdf,
        out_prefix=out_prefix,
        padding=padding,
    )
    best_score, xml_file, _ = score_ligand_with_prepared(
        prepared=prepared,
        ligand_sdf=ligand_sdf,
        nrun=nrun,
    )
    return best_score, xml_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run AutoDock-GPU pipeline and extract Vina-like score from XML."
    )
    parser.add_argument("--protein", "-r", required=True, help="Receptor PDB file")
    parser.add_argument("--ligand", "-l", required=True, help="Ligand SDF file")
    parser.add_argument("--out_prefix", "-o", help="Base name for receptor/gpf/maps files")
    parser.add_argument("--padding", "-p", type=float, default=6.0, help="Box padding in Å")
    parser.add_argument("--nrun", "-n", type=int, default=50, help="Number of docking runs")
    args = parser.parse_args()

    try:
        best, xml_path = run_full_pipeline(
            protein_pdb=args.protein,
            ligand_sdf=args.ligand,
            out_prefix=args.out_prefix,
            padding=args.padding,
            nrun=args.nrun,
        )
        print(f"XML output: {xml_path}")
        if best is not None:
            print(f"Best Vina-like score (kcal/mol): {best:.3f}")
        else:
            print("Could not determine best score from XML.")
    except AutoDockGPUPipelineError as e:
        print(f"[ERROR] {e}")
        raise SystemExit(1)
