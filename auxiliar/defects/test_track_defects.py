#!/usr/bin/env python3
"""Checks track_defects.py on synthetic snapshots whose defects are known:
  1. a box with a +1/2/-1/2 pair approaching and annihilating, a +1/2 leaving through the top wall,
     a fixed +1 defect and a region where the director escapes along z (no false defects there);
  2. a periodic domain whose pair crosses the x edge: it must stay two tracks.
Run: python3 auxiliar/defects/test_track_defects.py
"""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track_defects as td

failed = 0


def check(what, ok, value=''):
    global failed
    failed += not ok
    print('  %-60s %-14s %s' % (what, value, 'ok' if ok else 'FAILED'))


def write_snapshot(path, nx, ny, nz, S, pt):
    """CSV in the Marlics format for a uniaxial field (P=0, l normal to n)."""
    Ny, Nx = S.shape
    with open(path, 'w') as f:
        f.write('x,y,z,nx,ny,nz,lx,ly,lz,S,P,pt\n')
        for j in range(Ny):
            for i in range(Nx):
                n = np.array([nx[j, i], ny[j, i], nz[j, i]])
                l = np.cross(n, [0, 0, 1.]) if abs(n[2]) < 0.9 else np.cross(n, [1., 0, 0])
                l /= np.linalg.norm(l)
                f.write('%d,%d,0,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d\n'
                        % (i, j, n[0], n[1], n[2], l[0], l[1], l[2], S[j, i], 0., pt[j, i]))


def field(Nx, Ny, defects, escape=None):
    X, Y = np.meshgrid(np.arange(Nx), np.arange(Ny))
    phi = np.zeros((Ny, Nx)); S = 0.5 * np.ones((Ny, Nx))
    for x0, y0, k in defects:
        phi += k * np.arctan2(Y - y0, X - x0)
        S *= 1 - np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / 2.)
    tilt = np.zeros((Ny, Nx))
    if escape:                                  # director turning to z inside a disk
        x0, y0, R = escape
        r = np.hypot(X - x0, Y - y0)
        tilt = np.where(r < R, 0.5 * np.pi * (1 - r / R), 0.)
    return np.cos(phi) * np.cos(tilt), np.sin(phi) * np.cos(tilt), np.sin(tilt), S


def box_case(tmp):
    Nx, Ny = 61, 41
    pt = np.ones((Ny, Nx), int); pt[0, :] = pt[-1, :] = pt[:, 0] = pt[:, -1] = 2
    truth = []
    for f in range(10):
        d = [(10.4, 30.6, 1.0)]                                    # fixed +1
        if f < 7: d += [(15.3 + 2 * f, 15.6, 0.5), (45.3 - 2 * f, 15.6, -0.5)]   # pair
        if f < 6: d += [(52.6, 28.2 + 2 * f, 0.5)]                 # leaves through y=Ny-1 at f=6
        truth.append(d)
        write_snapshot(os.path.join(tmp, 'director_field_%d.csv' % f), *field(Nx, Ny, d, escape=(32, 32, 6)), pt)
    with open(os.path.join(tmp, 'marlics.log'), 'w') as fh:
        for f in range(10): fh.write('Snapshot %d: %g\n' % (f, 10. * f))

    td.main([tmp, '--out', os.path.join(tmp, 'box')])
    rows = np.genfromtxt(os.path.join(tmp, 'box_tracks.csv'), delimiter=',', names=True)
    err, wrong = 0., 0
    for f, d in enumerate(truth):
        found = rows[rows['snapshot'] == f]
        wrong += len(found) != len(d)
        for x0, y0, k in d:
            dist = [np.hypot(r['x'] - x0, r['y'] - y0) if r['charge'] == k else 1e9 for r in found]
            err = max(err, min(dist) if dist else 1e9)
    check('box: number of defects in every snapshot', wrong == 0, '%d wrong' % wrong)
    check('box: largest position error (sites)', err < 0.5, '%.2f' % err)
    check('box: no false defect in the escaped region',
          not any(np.hypot(r['x'] - 32, r['y'] - 32) < 7 for r in rows), '')
    check('box: 4 tracks (fixed +1, pair, leaving +1/2)', len(set(rows['track'])) == 4, '%d' % len(set(rows['track'])))
    times = sorted(set(rows['time']))
    check('box: times read from marlics.log', times[:3] == [0., 10., 20.], str(times[:3]))
    ev = open(os.path.join(tmp, 'box_events.csv')).read()
    check('box: annihilation of the pair at t=70', ',70,annihilation,' in ev)
    check('box: +1/2 leaving through the boundary at t=60', ',60,left through the boundary,' in ev)


def periodic_case(tmp):
    Nx, Ny = 80, 40
    pt = np.ones((Ny, Nx), int)
    nx, ny, nz, S = field(Nx, Ny, [(30.3, 20.6, 0.5), (50.3, 20.6, -0.5)])
    for f in range(10):
        s = 7 * f
        write_snapshot(os.path.join(tmp, 'director_field_%d.csv' % f),
                       np.roll(nx, s, 1), np.roll(ny, s, 1), np.roll(nz, s, 1), np.roll(S, s, 1), pt)
    td.main([tmp, '--out', os.path.join(tmp, 'per'), '--periodic', 'x,y', '--max-jump', '8'])
    rows = np.genfromtxt(os.path.join(tmp, 'per_tracks.csv'), delimiter=',', names=True)
    check('periodic: pair crossing the x edge stays two tracks', len(set(rows['track'])) == 2, '%d' % len(set(rows['track'])))
    plus = rows[rows['charge'] > 0]
    x_expected = (30.3 + 7 * plus['snapshot']) % Nx
    err = np.max(np.abs((plus['x'] - x_expected + Nx / 2) % Nx - Nx / 2))
    check('periodic: +1/2 position across the edge (sites)', err < 0.5, '%.2f' % err)


print('test_track_defects')
with tempfile.TemporaryDirectory() as t1, tempfile.TemporaryDirectory() as t2:
    box_case(t1)
    periodic_case(t2)
sys.exit(1 if failed else 0)
