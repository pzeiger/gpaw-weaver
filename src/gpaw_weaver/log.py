import re


def extract_scf_convergence(logfile):
    """Extract SCF iteration data from a GPAW txt log file.

    Handles collinear output (magmom: one float, e.g. ``+0.0000``),
    noncollinear output (magmom: three comma-separated floats), and
    non-magnetic calculations (no magmom column → ``'magmom'`` is ``None``).
    """
    iterations = []
    with open(logfile, 'r') as f:
        in_scf = False
        for line in f:
            if re.match(r'\s*iter\s+time\s+total', line):
                in_scf = True
                continue
            if in_scf and (line.strip() == '' or 'Converged' in line):
                break
            if in_scf:
                match = re.match(
                    r'iter:\s+(\d+)\s+\d{2}:\d{2}:\d{2}\s+'
                    r'([-+]?\d+\.\d+)c?\s+'
                    r'([-+]?\d+\.?\d*)?c?\s+'
                    r'([-+]?\d+\.?\d*)?c?'
                    r'(?:\s+([-+]?\d+\.\d+(?:,[-+]?\d+\.\d+,[-+]?\d+\.\d+)?))?',
                    line
                )
                if match:
                    magmom_str = match.group(5)
                    if magmom_str is None:
                        magmom = None
                    elif ',' in magmom_str:
                        magmom = [float(v) for v in magmom_str.split(',')]
                    else:
                        magmom = float(magmom_str)
                    iterations.append({
                        'iter': int(match.group(1)),
                        'energy': float(match.group(2)),
                        'log10_eigst': float(match.group(3)) if match.group(3) else None,
                        'log10_dens': float(match.group(4)) if match.group(4) else None,
                        'magmom': magmom,
                    })
    return iterations
