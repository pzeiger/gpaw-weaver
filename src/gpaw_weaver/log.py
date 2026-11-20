import re

# A row of GPAW's SCF table: ``iter:   2 ...`` (legacy implementation) or
# ``|iter:   2| ...`` (new implementation, cells separated by ``|``).
_ROW = re.compile(r'\|?iter:\s*(\d+)\s*\|?\s*\d{2}:\d{2}:\d{2}\s*\|?(.*)')
# GPAW writes the collinear magmom as ``%+.4f`` and the noncollinear one as a
# ``%+.1f`` triplet; eigst, dens and custom criteria have fewer decimals.
_MAGMOM = re.compile(
    r'[-+]\d+\.\d{4}|[-+]?\d+\.\d+,[-+]?\d+\.\d+,[-+]?\d+\.\d+')


def _number(cell):
    cell = cell.strip().rstrip('c')
    return float(cell) if cell else None


def _magmom(cell):
    if ',' in cell:
        return [float(v) for v in cell.split(',')]
    return float(cell)


def _parse_row(rest):
    """Return (energy, log10_eigst, log10_dens, magmom) from the part of a row
    after its time stamp."""
    if '|' in rest:
        # New layout: fixed cells, an empty cell is a missing value.
        cells = rest.rstrip().rstrip('|').split('|')
        energy, eigst, dens = (_number(c) for c in (cells + ['', ''])[:3])
        magmoms = [c.strip() for c in cells[3:]
                   if _MAGMOM.fullmatch(c.strip())]
    else:
        # Legacy layout: blank columns leave no token, so the magmom is
        # recognised by its format and the remaining numbers are taken in
        # column order (eigst, then dens).
        tokens = rest.split()
        magmoms = [t for t in tokens[1:] if _MAGMOM.fullmatch(t)]
        others = [_number(t) for t in tokens[1:] if t not in magmoms]
        energy = _number(tokens[0])
        eigst, dens = (others + [None, None])[:2]
    magmom = _magmom(magmoms[-1]) if magmoms else None
    return energy, eigst, dens, magmom


def extract_scf_convergence(logfile):
    """Extract SCF iteration data from a GPAW txt log file.

    Reads the ``iter:`` rows of the SCF table written by GPAW's legacy
    implementation (whitespace-separated, no header) and by its new one
    (``|``-separated cells), from the first row to the next blank line or
    ``Converged`` line. A value GPAW leaves blank (eigst and dens in the first
    iteration) is ``None``. ``'magmom'`` is a float for collinear
    spin-polarized output, a list of three floats for noncollinear output,
    and ``None`` for non-magnetic calculations.
    """
    iterations = []
    with open(logfile, 'r') as f:
        in_scf = False
        for line in f:
            if in_scf and (line.strip() == '' or 'Converged' in line):
                break
            match = _ROW.match(line.strip())
            if not match:
                continue
            in_scf = True
            energy, eigst, dens, magmom = _parse_row(match.group(2))
            iterations.append({
                'iter': int(match.group(1)),
                'energy': energy,
                'log10_eigst': eigst,
                'log10_dens': dens,
                'magmom': magmom,
            })
    return iterations
