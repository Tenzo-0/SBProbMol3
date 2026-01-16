import argparse
from ase import units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory

from fairchem.core import pretrained_mlip, FAIRChemCalculator


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Structure file: xyz/pdb/sdf/cif/...")
    ap.add_argument("--model", default="uma-s-1p1", help="UMA model name")
    ap.add_argument("--task", default="omol", help="Task name: omol/omat/oc20/odac/omc")
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--dt_fs", type=float, default=0.5)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--friction_fs", type=float, default=0.01, help="Langevin friction in 1/fs")
    ap.add_argument("--traj", default="md.traj")
    ap.add_argument("--log", default="md.log")
    args = ap.parse_args()

    atoms = read(args.input)

    predictor = pretrained_mlip.get_predict_unit(args.model, device=args.device)
    atoms.calc = FAIRChemCalculator(predictor, task_name=args.task)

    # Initialize velocities (NVT start)
    MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)

    dyn = Langevin(
        atoms,
        timestep=args.dt_fs * units.fs,
        temperature_K=args.temp,
        friction=args.friction_fs / units.fs,
        logfile=args.log,
    )

    traj = Trajectory(args.traj, "w", atoms)
    dyn.attach(traj.write, interval=10)

    dyn.run(args.steps)

    # Optional: also dump last frame as xyz
    write("final.xyz", atoms)
    print(f"Done. Wrote {args.traj}, {args.log}, final.xyz")


if __name__ == "__main__":
    main()
