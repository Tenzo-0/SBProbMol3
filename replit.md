# MoleculeAI - Drug Design Platform

## Overview

MoleculeAI is a web-based drug design platform that provides researchers with tools for computational drug discovery. The application focuses on protein target selection, binding site analysis, and molecular sampling for drug development workflows. It features an interactive interface for browsing protein structures, configuring molecular parameters, and executing drug design pipelines.

### Recent Enhancements (August 2025)
- **Add New Protein**: Users can now add custom protein targets to the system
- **Advanced Pocket Definition**: Auto-detection of binding pockets with manual override options
- **Comprehensive Sampling Parameters**: Detailed configuration interface matching research-grade tools
- **Enhanced Workflow**: Multi-step guided process from target selection to molecule generation

## User Preferences

Preferred communication style: Simple, everyday language.

## System Architecture

### Frontend Architecture
- **Template Engine**: Jinja2 templates with Flask for server-side rendering
- **UI Framework**: Bootstrap 5 for responsive design and component styling
- **Icon System**: Feather icons for consistent iconography
- **JavaScript Architecture**: Vanilla ES6 classes for application logic with modular design
- **Styling**: CSS custom properties for theming with a dark, professional color scheme

### Backend Architecture
- **Web Framework**: Flask with a simple, modular structure
- **Route Organization**: Separated route definitions in dedicated modules to avoid circular imports
- **Session Management**: Flask sessions with configurable secret keys for security
- **Data Layer**: Mock data implementation using Python dictionaries for protein information
- **Logging**: Python logging module configured for debugging and monitoring

### Application Structure
- **Entry Point**: `app.py` creates the Flask application instance with environment-based configuration
- **Route Handling**: `routes.py` defines URL endpoints for protein browsing, searching, detailed views, and new functionality
- **Data Management**: `data/proteins.py` provides comprehensive protein data with PDB IDs, molecular properties, and dynamic addition capabilities
- **Template Structure**: Base template system with inheritance for consistent layout and styling
- **Enhanced Templates**: Specialized interfaces for protein addition, pocket definition, and sampling parameters

### Workflow Design
- **Multi-step Process**: Four-stage drug design workflow (Target Selection, Pocket Definition, Sampling Parameters, Results)
- **Interactive Navigation**: Step-based progression with visual workflow indicators and page transitions
- **State Management**: Client-side JavaScript manages workflow state and user selections with session persistence
- **Pocket Detection**: Auto-detection algorithms with manual override options for binding site definition
- **Parameter Configuration**: Research-grade sampling controls with real-time validation

## External Dependencies

### Frontend Dependencies
- **Bootstrap 5**: CSS framework for responsive design and UI components
- **Feather Icons**: Lightweight icon library for user interface elements

### Backend Dependencies
- **Flask**: Python web framework for application structure and request handling

### Development Dependencies
- **Python Logging**: Built-in logging for debugging and application monitoring

### Data Sources
- **Protein Data Bank (PDB)**: References to PDB IDs for protein structures (mock implementation currently)
- **Mock Data**: Sample protein structures including HIV-1 Protease, SARS-CoV-2 Mpro, Human CDK2, and EGFR Kinase

### Potential Integration Points
- **Molecular Visualization**: Framework ready for integration with molecular viewers (Three.js, NGL, etc.)
- **Computational Backends**: Architecture supports integration with molecular dynamics and docking engines
- **Database Systems**: Current mock data layer can be replaced with relational or document databases
- **Authentication Systems**: Session framework ready for user authentication implementation