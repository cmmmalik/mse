# @cmmmalik/mse - Material Simulation Environment

## Overview

**@cmmmalik/mse** (Material Simulation Environment) is a Python package designed to streamline, automate, and manage computational materials science workflows. It provides tools to set up, run, analyze, and optimize high-throughput simulations, primarily using electronic structure codes such as GPAW and VASP. The package is tailored for researchers working on atomistic simulations, materials discovery, and computational chemistry.

## Features

- **Workflow Automation:** Easily define and execute complex simulation workflows for various materials and calculators (GPAW, VASP).
- **Job Management:** Classes to initialize, monitor, reset, and rerun computational jobs locally or on HPC clusters.
- **Server Tools:** SSH and SFTP handlers for seamless file transfer, job submission, and directory management on remote servers.
- **Optimization:** Interfaces for automated structure optimization and property convergence (k-points, energy cutoffs, etc.).
- **Database Utilities:** Integration with ASE databases and Materials Project for fetching and storing simulation data.
- **Structure Analysis:** Tools for symmetry analysis, energy evaluation, and formula extraction from simulated structures.
- **Directory Utilities:** Robust helpers for local and remote directory management, including rsync-based file synchronization.

## Repository Structure

- `mse/workflows/` – Workflow definitions and base workflow classes.
- `mse/Jobs/` – Core job classes and specific job implementations for GPAW and other calculators.
- `mse/optimizer/` – Optimization routines and interfaces.
- `mse/servertools/` – SSH, SFTP, and Slurm interfaces for remote job management.
- `mse/automatization/` – Tools for automatizing workflows and integration with external databases.
- `mse/formation_analysis/` – Tools for analyzing and extracting structural and energetic data.
- `mse/system/` – Directory and file operation utilities.

## Example Usage

```python
from mse.workflows.base import Baseworkflow
from ase.atoms import Atoms

atoms = Atoms('H2O')

# Initialize a workflow for GPAW
workflow = Baseworkflow(atoms=atoms, working_directory='h2o_sim', calculator_type='gpaw')
workflow.initialize_job(name="h2o_relaxation")
workflow.make_ready()
workflow.submit()
```

## Getting Started

1. **Installation:**
   - Clone this repository:
     ```
     git clone https://github.com/cmmmalik/mse.git
     ```
   - Install dependencies:
     ```
     pip install -r requirements.txt
     ```

2. **Configuration:**
   - Configure remote servers and calculators as required in your environment.

3. **Running Workflows:**
   - Use the workflow and job classes to define and submit your simulations.

## Requirements

- Python 3.6+
- ASE (Atomic Simulation Environment)
- GPAW, VASP (for calculations)
- Paramiko (for SSH/SFTP)
- pymatgen, monty, and other scientific libraries (see `requirements.txt`)

## Contributing

Contributions are welcome! Please submit pull requests or open issues for bugs and feature requests.

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- ASE, GPAW, VASP, pymatgen, and other open-source projects that make this package possible.
