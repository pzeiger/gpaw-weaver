"""Shared test helpers: fake log generators and FakeGPAW factory."""
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


def make_log(n_iters=3, magmom=None):
    """Return fake GPAW SCF log text.

    magmom=None          → non-magnetic (no magmom column)
    magmom='+1.0000'     → collinear
    magmom='+1.000,+0.000,-1.000' → non-collinear
    """
    lines = [" iter  time      total     log10-change:\n"]
    for i in range(1, n_iters + 1):
        e = -100.0 - i * 0.01
        mag = f"  {magmom}" if magmom is not None else ""
        lines.append(
            f"iter:  {i:2d}  10:30:{i:02d}  {e:.4f}   -2.000   -3.000{mag}\n"
        )
    lines.append("Converged\n")
    return "".join(lines)


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
