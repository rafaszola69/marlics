#!/usr/bin/env python3
"""Point defects of two-dimensional Marlics runs: detection, tracking and events.

Usage:
    python3 track_defects.py RUN_DIR [options]

RUN_DIR holds the snapshots director_field_<k>.csv. For each snapshot the Q tensor is rebuilt from
the CSV (n, l, S, P) and its in-plane part, Q_xx-Q_yy and Q_xy, gives the in-plane orientation
angle. The charge of each lattice plaquette is the winding of that angle around it. Using Q instead
of the director keeps the detection valid when the director leaves the plane: where it points
along z the in-plane part vanishes and the plaquette is skipped (--min-order), so escaped regions
do not produce false defects. Charged plaquettes that touch are merged into one defect; its position
is refined below the lattice spacing from the drop of the in-plane order around the core.

Defects are linked between consecutive snapshots by charge and distance (--max-jump). Tracks that
end or start together and close to each other are reported as annihilations or creations; tracks
ending or starting next to the boundary as leaving or entering through it.

Outputs (--out PREFIX, default RUN_DIR/defects):
    PREFIX_tracks.csv   track, snapshot, time, x, y, charge   (x, y in lattice units, or nm with --dx)
    PREFIX_events.csv   snapshot, time, event, x, y, tracks
    --plot FILE.png     trajectories over the last snapshot (needs matplotlib)

The snapshot times come from the Marlics output ("Snapshot k: t" lines): --log FILE, or
RUN_DIR/marlics.log when present; otherwise the snapshot number is used.
"""
import argparse, glob, os, re, sys
import numpy as np


def load_snapshot(path, layer=0):
    """Returns the in-plane order (a, b) = ((Q_xx-Q_yy)/2, Q_xy) and the point types of one layer."""
    d = np.loadtxt(path, delimiter=',', skiprows=1)
    Nx, Ny = int(d[:, 0].max()) + 1, int(d[:, 1].max()) + 1
    d = d[d[:, 2] == layer]
    i, j = d[:, 0].astype(int), d[:, 1].astype(int)
    n, l = d[:, 3:6], d[:, 6:9]
    S, P = d[:, 9], d[:, 10]
    # Eigenvalues: S (director n), (P-S)/2 (l) and the rest (m = n x l).
    e0, e1 = S, 0.5 * (P - S)
    e2 = -e0 - e1
    m = np.cross(n, l)
    Q = lambda a, b: e0 * n[:, a] * n[:, b] + e1 * l[:, a] * l[:, b] + e2 * m[:, a] * m[:, b]
    a = np.zeros((Ny, Nx)); b = np.zeros((Ny, Nx)); pt = np.zeros((Ny, Nx), int)
    a[j, i] = 0.5 * (Q(0, 0) - Q(1, 1))
    b[j, i] = Q(0, 1)
    pt[j, i] = d[:, 11].astype(int)
    return a, b, pt


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def detect(a, b, pt, periodic=(False, False), min_order=0.05):
    """Defects of one snapshot: list of (x, y, charge)."""
    Ny, Nx = a.shape
    amp = np.hypot(a, b)
    ang = np.arctan2(b, a)                      # twice the in-plane orientation angle
    ok_site = (pt > 0) & (amp > min_order * np.median(amp[pt > 0]))

    # Plaquette (i, j) has corners (i, j), (i+1, j), (i+1, j+1), (i, j+1), counterclockwise.
    ip = (np.arange(Nx) + 1) % Nx; jp = (np.arange(Ny) + 1) % Ny
    A = ang; B = ang[:, ip]; C = ang[jp][:, ip]; D = ang[jp]
    w = (wrap(B - A) + wrap(C - B) + wrap(D - C) + wrap(A - D)) / (4 * np.pi)
    ok = ok_site & ok_site[:, ip] & ok_site[jp][:, ip] & ok_site[jp]
    if not periodic[0]: ok[:, -1] = False
    if not periodic[1]: ok[-1, :] = False
    charged = ok & (np.abs(w) > 0.25)
    q = np.where(charged, np.round(2 * w) / 2, 0.)

    # Merge touching charged plaquettes (8-neighbourhood, periodic where the domain is).
    seen = np.zeros_like(charged)
    defects = []
    for j0, i0 in zip(*np.nonzero(charged)):
        if seen[j0, i0]: continue
        stack, cells = [(j0, i0)], []
        seen[j0, i0] = True
        while stack:
            j, i = stack.pop(); cells.append((j, i))
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    jj, ii = j + dj, i + di
                    if periodic[1]: jj %= Ny
                    if periodic[0]: ii %= Nx
                    if 0 <= jj < Ny and 0 <= ii < Nx and charged[jj, ii] and not seen[jj, ii]:
                        seen[jj, ii] = True; stack.append((jj, ii))
        charge = sum(q[c] for c in cells)
        if charge == 0: continue                # a dipole within one lattice spacing
        x, y = _centre(cells, Nx, Ny, periodic)
        x, y = _refine(x, y, amp, pt, periodic)
        defects.append((x, y, charge))
    return defects


def _centre(cells, Nx, Ny, periodic):
    """Mean of the plaquette centres, unwrapped around periodic edges."""
    xs = np.array([c[1] for c in cells], float); ys = np.array([c[0] for c in cells], float)
    if periodic[0]: xs = xs[0] + wrap((xs - xs[0]) * 2 * np.pi / Nx) * Nx / (2 * np.pi)
    if periodic[1]: ys = ys[0] + wrap((ys - ys[0]) * 2 * np.pi / Ny) * Ny / (2 * np.pi)
    return xs.mean() + 0.5, ys.mean() + 0.5


def _refine(x, y, amp, pt, periodic):
    """Core position from the drop of the in-plane order in the 4x4 sites around (x, y)."""
    Ny, Nx = amp.shape
    i0, j0 = int(np.floor(x)) - 1, int(np.floor(y)) - 1
    ws, xs, ys = [], [], []
    ref = np.median(amp[pt > 0])
    for j in range(j0, j0 + 4):
        for i in range(i0, i0 + 4):
            ii, jj = (i % Nx if periodic[0] else i), (j % Ny if periodic[1] else j)
            if not (0 <= ii < Nx and 0 <= jj < Ny) or pt[jj, ii] == 0: continue
            wgt = max(ref - amp[jj, ii], 0.) ** 2
            ws.append(wgt); xs.append(i); ys.append(j)
    if not ws or sum(ws) == 0: return x, y
    ws = np.array(ws)
    xr, yr = np.dot(ws, xs) / ws.sum(), np.dot(ws, ys) / ws.sum()
    # Stay within the charged plaquettes found by the winding.
    return float(np.clip(xr, x - 0.5, x + 0.5)), float(np.clip(yr, y - 0.5, y + 0.5))


def distance(p, q, shape, periodic):
    Ny, Nx = shape
    dx, dy = p[0] - q[0], p[1] - q[1]
    if periodic[0]: dx = (dx + Nx / 2) % Nx - Nx / 2
    if periodic[1]: dy = (dy + Ny / 2) % Ny - Ny / 2
    return np.hypot(dx, dy)


def track(frames, shape, periodic, max_jump):
    """Links the defects of consecutive snapshots. Returns rows (track, frame, x, y, charge)."""
    rows, active, next_id = [], {}, 0          # active: track id -> (x, y, charge)
    for f, defects in enumerate(frames):
        pairs = sorted((distance(active[t][:2], d[:2], shape, periodic), t, k)
                       for t in active for k, d in enumerate(defects) if active[t][2] == d[2])
        used_t, used_k, new_active = set(), set(), {}
        for dist, t, k in pairs:
            if dist > max_jump or t in used_t or k in used_k: continue
            used_t.add(t); used_k.add(k); new_active[t] = defects[k]
        for k, d in enumerate(defects):
            if k not in used_k:
                new_active[next_id] = d; next_id += 1
        for t, d in new_active.items():
            rows.append((t, f, d[0], d[1], d[2]))
        active = new_active
    return rows


def events(rows, n_frames, pts, shape, periodic, max_jump):
    """Annihilations/creations (opposite tracks ending/starting together and close) and
    tracks leaving/entering next to the boundary."""
    first, last = {}, {}
    for t, f, x, y, q in rows:
        if t not in first: first[t] = (f, x, y, q)
        last[t] = (f, x, y, q)
    out = []
    for kind, ends, step in (('annihilation', last, 1), ('creation', first, 0)):
        for f in range(n_frames):
            group = [t for t, e in ends.items() if e[0] == f and (f < n_frames - 1 if step else f > 0)]
            used = set()
            for t in sorted(group, key=lambda t: ends[t][3]):
                if t in used: continue
                # Partners: the opposite charges closest to it, until the charges add up to zero.
                partners, charge = [t], ends[t][3]
                for u in sorted(group, key=lambda u: distance(ends[t][1:3], ends[u][1:3], shape, periodic)):
                    if u in used or u in partners or charge == 0: continue
                    if np.sign(ends[u][3]) != np.sign(charge) and \
                       distance(ends[t][1:3], ends[u][1:3], shape, periodic) <= 2 * max_jump:
                        partners.append(u); charge += ends[u][3]
                if charge == 0 and len(partners) > 1:
                    used.update(partners)
                    x = np.mean([ends[u][1] for u in partners]); y = np.mean([ends[u][2] for u in partners])
                    out.append((f + step, kind, x, y, partners))
            for t in group:
                if t in used: continue
                f0, x, y, q = ends[t]
                at_wall = _near_boundary(x, y, pts[f0], periodic, 3)
                if step: out.append((f + 1, 'left through the boundary' if at_wall else 'lost', x, y, [t]))
                else:    out.append((f, 'entered through the boundary' if at_wall else 'appeared', x, y, [t]))
    return sorted(out, key=lambda e: e[0])


def _near_boundary(x, y, pt, periodic, r):
    Ny, Nx = pt.shape
    for j in range(int(y) - r, int(y) + r + 2):
        for i in range(int(x) - r, int(x) + r + 2):
            ii, jj = (i % Nx if periodic[0] else i), (j % Ny if periodic[1] else j)
            if not (0 <= ii < Nx and 0 <= jj < Ny) or pt[jj, ii] != 1: return True
    return False


def snapshot_times(run_dir, log, numbers):
    log = log or os.path.join(run_dir, 'marlics.log')
    times = {}
    if os.path.exists(log):
        for line in open(log, errors='replace'):
            m = re.match(r'Snapshot (\d+): (\S+)', line)
            if m: times[int(m.group(1))] = float(m.group(2))
    return [times.get(k, float(k)) for k in numbers], bool(times)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('run_dir')
    p.add_argument('--log', help='Marlics output with the "Snapshot k: t" lines (default RUN_DIR/marlics.log)')
    p.add_argument('--out', help='output prefix (default RUN_DIR/defects)')
    p.add_argument('--periodic', default='', help='periodic axes, e.g. x,y (bulk geometry) or x (channel)')
    p.add_argument('--max-jump', type=float, default=5., help='largest displacement between snapshots, sites (5)')
    p.add_argument('--min-order', type=float, default=0.05,
                   help='skip sites whose in-plane order is below this fraction of the median (0.05)')
    p.add_argument('--layer', type=int, default=0, help='z layer to analyse (0)')
    p.add_argument('--start', type=int, default=0,
                   help='first snapshot number to analyse, e.g. 1 to skip a random initial state (0)')
    p.add_argument('--dx', type=float, default=1., help='lattice spacing for the positions written (1: sites)')
    p.add_argument('--plot', help='figure of the trajectories (PNG)')
    a = p.parse_args(argv)

    files = sorted(glob.glob(os.path.join(a.run_dir, 'director_field_*.csv')),
                   key=lambda f: int(re.findall(r'(\d+)\.csv$', f)[0]))
    files = [f for f in files if int(re.findall(r'(\d+)\.csv$', f)[0]) >= a.start]
    if not files: sys.exit('No director_field_*.csv in ' + a.run_dir)
    numbers = [int(re.findall(r'(\d+)\.csv$', f)[0]) for f in files]
    times, from_log = snapshot_times(a.run_dir, a.log, numbers)
    periodic = ('x' in a.periodic, 'y' in a.periodic)

    frames, pts = [], []
    for f in files:
        qa, qb, pt = load_snapshot(f, a.layer)
        frames.append(detect(qa, qb, pt, periodic, a.min_order)); pts.append(pt)
    shape = pts[0].shape
    rows = track(frames, shape, periodic, a.max_jump)
    evs = events(rows, len(frames), pts, shape, periodic, a.max_jump)

    out = a.out or os.path.join(a.run_dir, 'defects')
    with open(out + '_tracks.csv', 'w') as fh:
        fh.write('track,snapshot,time,x,y,charge\n')
        for t, f, x, y, q in rows:
            fh.write('%d,%d,%g,%.3f,%.3f,%g\n' % (t, numbers[f], times[f], x * a.dx, y * a.dx, q))
    with open(out + '_events.csv', 'w') as fh:
        fh.write('snapshot,time,event,x,y,tracks\n')
        for f, kind, x, y, tr in evs:
            fh.write('%d,%g,%s,%.3f,%.3f,%s\n' % (numbers[f], times[f], kind, x * a.dx, y * a.dx, ' '.join(map(str, tr))))

    n_tracks = len({r[0] for r in rows})
    print('%d snapshots%s, %d tracks, %d events -> %s_tracks.csv, %s_events.csv'
          % (len(files), '' if from_log else ' (no log: time = snapshot number)', n_tracks, len(evs), out, out))
    if len(evs) <= 20:
        for f, kind, x, y, tr in evs:
            print('  t=%g: %s at (%.1f, %.1f), tracks %s' % (times[f], kind, x * a.dx, y * a.dx, tr))
    else:
        kinds = sorted({e[1] for e in evs})
        print('  ' + ', '.join('%d %s' % (sum(e[1] == k for e in evs), k) for k in kinds) + ' (see the events file)')

    if a.plot:
        _plot(a.plot, rows, files[-1], a.layer, times)


def _plot(path, rows, last_file, layer, times):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    qa, qb, pt = load_snapshot(last_file, layer)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(np.where(pt > 0, np.hypot(qa, qb), np.nan), origin='lower', cmap='Greys_r', alpha=0.6)
    for t in sorted({r[0] for r in rows}):
        tr = [r for r in rows if r[0] == t]
        q = tr[0][4]
        c = 'red' if q > 0 else 'deepskyblue'
        ax.plot([r[2] for r in tr], [r[3] for r in tr], '-', color=c, lw=1.5)
        ax.plot(tr[0][2], tr[0][3], 'o' if q > 0 else 's', mfc='white', mec=c, ms=6)
        ax.plot(tr[-1][2], tr[-1][3], 'o' if q > 0 else 's', color=c, mec='black', ms=7)
    ax.set_title('Defect trajectories (red: +, blue: -; open: start, filled: end)', fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(path, dpi=110, bbox_inches='tight')


if __name__ == '__main__':
    main()
