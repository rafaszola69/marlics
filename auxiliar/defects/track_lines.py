#!/usr/bin/env python3
"""Disclination lines of three-dimensional Marlics runs.

Usage:
    python3 track_lines.py RUN_DIR [options]

For each snapshot director_field_<k>.csv, the director is carried around every lattice plaquette
(xy, yz and zx faces), choosing at each corner the sign that continues the previous one; if it comes
back reversed, a disclination line crosses the plaquette (the Z2 charge of 3D nematics). Crossed
plaquettes of the same lattice cube are joined, which builds the lines (loops in periodic boxes).

Along each line the local profile is Omega.t, the cosine between the axis about which the director
turns around the line and the line tangent: about +1 for a +1/2 wedge profile, -1 for a -1/2 wedge,
0 for a twist profile. Lines often change profile along their length.

Outputs (--out PREFIX, default RUN_DIR/lines):
    PREFIX_points.csv    snapshot, time, line, x, y, z, profile   (crossed plaquette centres, sites)
    PREFIX_summary.csv   snapshot, time, lines, length, wedge_plus, wedge_minus, twist
                         (length in lattice spacings; fractions of the length by profile:
                          Omega.t > 0.5, < -0.5, otherwise)
    --plot FILE.png      lines of the last snapshot and the length against time (needs matplotlib)

Times come from the "Snapshot k: t" lines of the Marlics output (--log, or RUN_DIR/marlics.log).
"""
import argparse, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from track_defects import snapshot_times


def load_director(path):
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    Nx, Ny, Nz = (int(d[:, c].max()) + 1 for c in range(3))
    n = np.zeros((Nz, Ny, Nx, 3)); pt = np.zeros((Nz, Ny, Nx), int)
    i, j, k = (d[:, c].astype(int) for c in range(3))
    n[k, j, i] = d[:, 3:6]
    pt[k, j, i] = d[:, 11].astype(int)
    return n, pt


def shift(a, axis, periodic):
    """a at the next site along axis (0: x, 1: y, 2: z); wraps around (masked later if not periodic)."""
    return np.roll(a, -1, axis=2 - axis)


def crossed_plaquettes(n, pt, periodic):
    """For each orientation (normal along x, y, z): boolean array of crossed plaquettes, their profile
    Omega.t, and validity. Plaquette with normal a at site s has corners s, s+e_b, s+e_b+e_c, s+e_c
    with (a, b, c) cyclic, i.e. counterclockwise seen from +a."""
    out = []
    ok_site = pt > 0
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        n0 = n
        n1 = shift(n, b, periodic)
        n2 = shift(n1, c, periodic)
        n3 = shift(n, c, periodic)
        ok = ok_site & shift(ok_site, b, periodic) & shift(shift(ok_site, b, periodic), c, periodic) & shift(ok_site, c, periodic)
        for ax in (b, c):
            if not periodic[ax]:
                idx = [slice(None)] * 3; idx[2 - ax] = -1; ok[tuple(idx)] = False
        m = [n0]
        for nk in (n1, n2, n3):
            s = np.sign(np.sum(m[-1] * nk, axis=-1, keepdims=True)); s[s == 0] = 1
            m.append(nk * s)
        back = np.sum(m[-1] * n0, axis=-1)
        crossed = ok & (back < 0)
        # Rotation axis: sum of the turns between consecutive directors, the last one towards -n0.
        m.append(-n0)
        omega = sum(np.cross(m[q], m[q + 1]) for q in range(4))
        norm = np.linalg.norm(omega, axis=-1); norm[norm == 0] = 1
        profile = omega[..., a] / norm
        out.append((crossed, profile))
    return out


def build_lines(planes, shape, periodic):
    """Joins crossed plaquettes sharing a lattice cube. Returns a list of lines, each a list of
    (x, y, z, profile) points at the plaquette centres."""
    Nz, Ny, Nx = shape
    faces = []                               # (a, k, j, i)
    for a, (crossed, _) in enumerate(planes):
        for k, j, i in zip(*np.nonzero(crossed)):
            faces.append((a, k, j, i))
    index = {f: q for q, f in enumerate(faces)}
    parent = list(range(len(faces)))

    def find(q):
        while parent[q] != q:
            parent[q] = parent[parent[q]]; q = parent[q]
        return q

    # The cube with corner (i, j, k) has, for each orientation a, the faces at site s and s+e_a.
    cubes = {}
    for q, (a, k, j, i) in enumerate(faces):
        for back in (0, 1):
            c = [i, j, k]; c[a] -= back
            if c[a] < 0:
                if not periodic[a]: continue
                c[a] += (Nx, Ny, Nz)[a]
            cubes.setdefault(tuple(c), []).append(q)
    for members in cubes.values():
        for q in members[1:]:
            parent[find(q)] = find(members[0])

    lines = {}
    for q, (a, k, j, i) in enumerate(faces):
        x, y, z = i + 0.5, j + 0.5, k + 0.5
        if a == 0: x = i
        elif a == 1: y = j
        else: z = k
        lines.setdefault(find(q), []).append((x, y, z, planes[a][1][k, j, i]))
    return sorted(lines.values(), key=len, reverse=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('run_dir')
    p.add_argument('--log', help='Marlics output with the "Snapshot k: t" lines (default RUN_DIR/marlics.log)')
    p.add_argument('--out', help='output prefix (default RUN_DIR/lines)')
    p.add_argument('--periodic', default='x,y,z', help='periodic axes (default x,y,z: bulk geometry)')
    p.add_argument('--start', type=int, default=0, help='first snapshot number to analyse (0)')
    p.add_argument('--min-length', type=int, default=3, help='ignore lines shorter than this, in plaquettes (3)')
    p.add_argument('--plot', help='figure of the lines of the last snapshot and of the length (PNG)')
    a = p.parse_args(argv)

    files = sorted(glob.glob(os.path.join(a.run_dir, 'director_field_*.csv')),
                   key=lambda f: int(re.findall(r'(\d+)\.csv$', f)[0]))
    files = [f for f in files if int(re.findall(r'(\d+)\.csv$', f)[0]) >= a.start]
    if not files: sys.exit('No director_field_*.csv in ' + a.run_dir)
    numbers = [int(re.findall(r'(\d+)\.csv$', f)[0]) for f in files]
    times, from_log = snapshot_times(a.run_dir, a.log, numbers)
    periodic = tuple(ax in a.periodic for ax in 'xyz')

    out = a.out or os.path.join(a.run_dir, 'lines')
    summary, last = [], None
    with open(out + '_points.csv', 'w') as fp:
        fp.write('snapshot,time,line,x,y,z,profile\n')
        for f, num, t in zip(files, numbers, times):
            n, pt = load_director(f)
            lines = [l for l in build_lines(crossed_plaquettes(n, pt, periodic), pt.shape, periodic) if len(l) >= a.min_length]
            prof = np.array([q[3] for l in lines for q in l]) if lines else np.zeros(0)
            L = len(prof)
            summary.append((num, t, len(lines), L,
                            (prof > 0.5).sum() / max(L, 1), (prof < -0.5).sum() / max(L, 1), (abs(prof) <= 0.5).sum() / max(L, 1)))
            for li, l in enumerate(lines):
                for x, y, z, pr in l:
                    fp.write('%d,%g,%d,%.1f,%.1f,%.1f,%.3f\n' % (num, t, li, x, y, z, pr))
            last = lines
    with open(out + '_summary.csv', 'w') as fs:
        fs.write('snapshot,time,lines,length,wedge_plus,wedge_minus,twist\n')
        for r in summary:
            fs.write('%d,%g,%d,%d,%.3f,%.3f,%.3f\n' % r)

    print('%d snapshots%s -> %s_points.csv, %s_summary.csv' % (len(files), '' if from_log else ' (no log: time = snapshot number)', out, out))
    for r in summary[:: max(1, len(summary) // 10)] + ([summary[-1]] if len(summary) > 1 else []):
        print('  t=%g: %d lines, length %d, profile +1/2 %.0f%%, -1/2 %.0f%%, twist %.0f%%'
              % (r[1], r[2], r[3], 100 * r[4], 100 * r[5], 100 * r[6]))
    if a.plot:
        _plot(a.plot, last, summary, pt.shape)


def _plot(path, lines, summary, shape):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    for l in lines:
        q = np.array(l)
        ax.scatter(q[:, 0], q[:, 1], q[:, 2], c=q[:, 3], cmap='coolwarm', vmin=-1, vmax=1, s=4)
    Nz, Ny, Nx = shape
    ax.set_xlim(0, Nx); ax.set_ylim(0, Ny); ax.set_zlim(0, Nz)
    ax.set_title('Lines, last snapshot (colour: profile, red +1/2 wedge, blue -1/2, white twist)', fontsize=9)
    bx = fig.add_subplot(1, 2, 2)
    s = np.array(summary, float)
    bx.plot(s[:, 1], s[:, 3], 'k-')
    bx.set_yscale('log'); bx.set_xlabel('t'); bx.set_ylabel('total length (lattice spacings)')
    bx2 = bx.twinx(); bx2.plot(s[:, 1], s[:, 2], 'g--'); bx2.set_ylabel('number of lines', color='g')
    fig.savefig(path, dpi=110, bbox_inches='tight')


if __name__ == '__main__':
    main()
