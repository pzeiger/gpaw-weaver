# gpaw-weaver

Small helpers for running and organizing [GPAW](https://gpaw.readthedocs.io) calculations with [ASE](https://wiki.fysik.dtu.dk/ase/).

- Build typed parameter dicts for plane-wave and finite-difference calculations
- Run a GPAW calculation and store the initial geometry, converged geometry, and full SCF history in an ASE database in one call
- Reload a stored calculation (atoms + calculator) from the database

## Requirements

- Python ≥ 3.10
- [GPAW](https://gpaw.readthedocs.io/install.html) and its dependencies (OpenBLAS / ScaLAPACK, libxc, …)
- [ASE](https://wiki.fysik.dtu.dk/ase/install.html)

A pre-configured devcontainer is included (see [`.devcontainer/`](.devcontainer/)) that builds GPAW with MPI, ScaLAPACK, libxc, libvdwxc, and FFTW on Ubuntu.

## Installation

Install directly from GitHub (no PyPI release yet):

```bash
pip install git+https://github.com/pzeiger/gpaw-weaver.git
```

For development, clone and install in editable mode:

```bash
git clone https://github.com/pzeiger/gpaw-weaver.git
cd gpaw-weaver
pip install -e ".[dev]"
```

## Usage

### Building calculation parameters

```python
from gpaw_weaver import make_pw_params, make_fd_params, get_mode_filestr

# Plane-wave calculation
params = make_pw_params(
    ecut=500,
    kpts={"size": (8, 8, 8), "gamma": True},
    xc="PBE",
    convergence={"density": 1e-6},
)

# Finite-difference calculation — specify either h or gpts
params = make_fd_params(h=0.18, kpts={"size": (4, 4, 4), "gamma": True})
params = make_fd_params(gpts=(48, 48, 48), kpts={"size": (4, 4, 4)})

# Human-readable string for labelling files or directories
label = get_mode_filestr(params)  # e.g. 'pw_ecut500_8x8x8_gamma'
```

Any extra keyword arguments are forwarded to GPAW and stored in the database:

```python
params = make_pw_params(500, kpts, nbands=40, charge=0, setups={"Fe": ":d,4.0"})
```

### Running and storing a calculation

```python
from ase import Atoms
from ase.db import connect
from gpaw_weaver import make_pw_params, run_and_store_gpaw_calculation

atoms = Atoms("Fe2", positions=[(0, 0, 0), (1.435, 1.435, 1.435)],
              cell=[2.87, 2.87, 2.87], pbc=True)
atoms.set_initial_magnetic_moments([2.2, 2.2])

params = make_pw_params(500, {"size": (8, 8, 8), "gamma": True},
                        convergence={"density": 1e-6})

db = connect("calculations.db")

atoms_converged, initial_id, converged_id = run_and_store_gpaw_calculation(
    atoms, params, system="bcc-Fe", db=db,
    save_gpw=True,          # write a .gpw restart file
    legacy_gpaw=True,       # False for the new GPAW implementation
)
```

This writes two rows to the database:
- **initial** — starting geometry plus all calculation parameters
- **converged** — relaxed geometry plus full SCF history (`scf_energies`, `scf_log10_eigst`, `scf_log10_dens`, `scf_magmoms`)

The initial row gains a `converged_id` pointer and the converged row gains an `initial_id` pointer so the pair is always findable from either end.

Log files are written to `gpw_logs/<id>.txt` and restart files to `gpw_files/<id>.gpw`, where `<id>` is the database row ID of the converged entry. Both directories are created automatically.

### Loading a stored calculation

```python
from gpaw_weaver import load_gpaw_calculation

atoms, calc = load_gpaw_calculation(
    params,               # calc_params used when the calculation was run
    db,
    system="bcc-Fe",
)
```

### Querying the database

The database is a standard ASE SQLite database and can be queried with the usual ASE tools:

```bash
ase db calculations.db
ase db calculations.db system=bcc-Fe
```

```python
for row in db.select(system="bcc-Fe"):
    print(row.id, row.energy, row.data["scf_energies"])
```

## Running tests

```bash
pip install -e ".[dev]"
pytest
```

All calculation tests run against both `legacy_gpaw=True` and `legacy_gpaw=False` automatically.

## License

MIT — see [LICENSE](LICENSE).
