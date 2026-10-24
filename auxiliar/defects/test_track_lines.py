#!/usr/bin/env python3
"""Checks track_lines.py on synthetic 3D director fields with known disclination lines:
straight +1/2 and -1/2 wedge lines, a pair of twist lines, a closed loop and a uniform field.
Run: python3 auxiliar/defects/test_track_lines.py
"""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track_lines as tl

failed = 0


def check(what, ok, value=''):
    global failed
    failed += not ok
    print('  %-58s %-16s %s' % (what, value, 'ok' if ok else 'FAILED'))


def lines_of(n, periodic):
    pt = np.ones(n.shape[:3], int)
    return tl.build_lines(tl.crossed_plaquettes(n, pt, periodic), pt.shape, periodic)


def grid(Nx, Ny, Nz):
    return np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz), indexing='xy')


def field(nx, ny, nz):
    """(Nz, Ny, Nx, 3) array from meshgrid-ordered (Ny, Nx, Nz) components."""
    n = np.stack([nx, ny, nz], axis=-1).transpose(2, 0, 1, 3)
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


Nx, Ny, Nz = 24, 16, 10
X, Y, Z = grid(Nx, Ny, Nz)
th1, th2 = np.arctan2(Y - 8.4, X - 6.3), np.arctan2(Y - 8.4, X - 17.3)

# Wedge pair along z: +1/2 at x=6.3, -1/2 at x=17.3.
phi = 0.5 * th1 - 0.5 * th2
L = lines_of(field(np.cos(phi), np.sin(phi), 0 * phi), (True, True, True))
check('wedge pair: two lines', len(L) == 2, str(len(L)))
check('wedge pair: each through the box (length Nz=10)', [len(l) for l in L] == [10, 10], str([len(l) for l in L]))
for l in L:
    x = np.mean([q[0] for q in l]); pr = np.mean([q[3] for q in l])
    k = '+1/2' if x < 12 else '-1/2'
    check('wedge %s at x=%.1f: profile Omega.t' % (k, x), (pr > 0.9) if k == '+1/2' else (pr < -0.9), '%.3f' % pr)

# Twist pair along z: the director turns about x, perpendicular to the lines.
phi = 0.5 * th1 - 0.5 * th2
L = lines_of(field(0 * phi, np.cos(phi), np.sin(phi)), (True, True, True))
prof = np.array([q[3] for l in L for q in l])
check('twist pair: two lines', len(L) == 2, str(len(L)))
check('twist pair: profile |Omega.t| (max)', np.abs(prof).max() < 0.2, '%.3f' % np.abs(prof).max())

# Loop of radius R=6 in the plane z=12.4, centred at (12.3, 12.6): n turns about the azimuthal
# direction (the loop tangent) around the loop core. Not periodic.
N = 25; R = 6.
X, Y, Z = grid(N, N, N)
rho = np.hypot(X - 12.3, Y - 12.6); az = np.arctan2(Y - 12.6, X - 12.3)
phi = 0.5 * np.arctan2(Z - 12.4, rho - R)
n = field(np.cos(phi) * np.cos(az), np.cos(phi) * np.sin(az), np.sin(phi))
L = lines_of(n, (False, False, False))
check('loop: one closed line', len(L) == 1, str(len(L)))
r = np.array([np.hypot(q[0] - 12.3, q[1] - 12.6) for q in L[0]])
check('loop: radius of its points (expected 6)', abs(r.mean() - R) < 0.5, '%.2f' % r.mean())
check('loop: plaquettes crossed ~ 8R (L1 length of a circle)', 6 * R < len(L[0]) < 9 * R, str(len(L[0])))
check('loop: wedge profile |Omega.t| (mean)', np.mean([abs(q[3]) for q in L[0]]) > 0.8, '%.2f' % np.mean([abs(q[3]) for q in L[0]]))

# Uniform field: no line.
check('uniform field: no line', len(lines_of(field(1 + 0 * X, 0 * X, 0 * X), (True, True, True))) == 0)

# Whole tool on snapshots written in the Marlics format.
with tempfile.TemporaryDirectory() as tmp:
    X, Y, Z = grid(Nx, Ny, Nz)
    with open(os.path.join(tmp, 'director_field_0.csv'), 'w') as f:
        f.write('x,y,z,nx,ny,nz,lx,ly,lz,S,P,pt\n')
        ph = 0.5 * np.arctan2(Y - 8.4, X - 6.3) - 0.5 * np.arctan2(Y - 8.4, X - 17.3)
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    a = ph[j, i, k]
                    f.write('%d,%d,%d,%.3f,%.3f,0,%.3f,%.3f,0,0.5,0,1\n' % (i, j, k, np.cos(a), np.sin(a), -np.sin(a), np.cos(a)))
    tl.main([tmp, '--out', os.path.join(tmp, 'l')])
    s = open(os.path.join(tmp, 'l_summary.csv')).read().split('\n')[1].split(',')
    check('tool on a CSV snapshot: lines, length', s[2:4] == ['2', '20'], ','.join(s[2:4]))

sys.exit(1 if failed else 0)
