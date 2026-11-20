"""Shared test helpers: fake log generators and FakeGPAW factory."""
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


def legacy_row(i, energy="", eigst="", dens="", custom=(), magmom=None):
    """Return one SCF row as legacy GPAW writes it (gpaw/scf.py), with blank
    cells for empty strings."""
    line = f"iter:{i:4d} 10:30:{i:02d} {energy:>12s}  {eigst:>6s}  {dens:>5s}  "
    for value in custom:
        line += f"{value:>5s}  "
    if magmom is not None:
        line += f"  {magmom}" if "," not in magmom else f" {magmom}"
    return line.rstrip() + "\n"


def make_log(n_iters=3, magmom=None):
    """Return SCF log text in legacy GPAW's layout.

    magmom=None          → non-magnetic (no magmom column)
    magmom='+1.0000'     → collinear
    magmom='+1.000,+0.000,-1.000' → non-collinear
    """
    lines = [legacy_row(i, f"{-100.0 - i * 0.01:.4f}", "-2.00", "-3.00",
                        magmom=magmom)
             for i in range(1, n_iters + 1)]
    return "".join(lines) + "\nConverged after 3 iterations.\n"


def make_fake_gpaw_class(n_spins=1, log_content=None):
    """Return a FakeGPAW Calculator class with the given spin config.

    The class writes *log_content* to the txt path on construction (matching
    real GPAW behaviour) and exposes the minimal interface used by
    run_and_store_gpaw_calculation / load_gpaw_calculation.
    """

    class FakeGPAW(Calculator):
        implemented_properties = ["energy", "forces", "magmoms"]
        default_parameters = {}
        name = "fakegpaw"

        def __init__(self, *args, txt=None, **kwargs):
            super().__init__()
            # Record construction arguments so tests can assert what
            # run_and_store forwarded to the GPAW constructor.
            self.init_args = args
            self.init_kwargs = kwargs
            if txt and log_content is not None:
                Path(txt).parent.mkdir(parents=True, exist_ok=True)
                Path(txt).write_text(log_content)

        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            n = len(atoms) if atoms is not None else 1
            self.results = {
                "energy": -100.0,
                "forces": np.zeros((n, 3)),
            }
            if n_spins > 1:
                self.results["magmoms"] = np.ones(n)

        def get_number_of_spins(self):
            return n_spins

        def write(self, path, mode="calculation"):
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).touch()

    return FakeGPAW
