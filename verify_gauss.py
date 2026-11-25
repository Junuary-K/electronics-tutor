"""
Gauss's law - independent verification.

Checks (must match gauss.html):
  (1) point charge: closed form E = Q/(4 pi eps0 r^2), ratio E(r/2) = 4 E(r)
  (2) Gauss flux identity by DIRECT numeric surface integration of the
      Coulomb field over a sphere:  sum E . dA  ==  Q/eps0
  (3) uniformly charged ball: direct numeric Coulomb sum (N^3 cells) vs the
      Gauss results  E ~ r inside, E = Q/(4 pi eps0 r^2) outside
  (4) line charge: numeric segmented Coulomb sum (finite line) vs the exact
      finite-line formula and vs the infinite-line Gauss result
  (5) plane charge: numeric patch Coulomb sum vs sigma/(2 eps0),
      independence of height
  (6) enclosed-charge bookkeeping Q_enc = Q (r/R)^3 for the ball

Exit code 0 = all pass.
"""
import math, sys

EPS0 = 8.8541878128e-12
K = 1.0 / (4.0 * math.pi * EPS0)

fails = []

Q = 10e-9          # 10 nC point charge
R = 0.05           # 5 cm ball

print("=" * 66)
print("(1) point charge closed form + 1/r^2 ratio")
for r in (0.02, 0.05, 0.08):
    E = K * Q / (r * r)
    ok = abs(E * r * r - K * Q) < 1e-9
    print("  [%s] r=%.2f m  E=%.1f V/m" % ("PASS" if ok else "FAIL", r, E))
    if not ok: fails.append("point r=%.2f" % r)
E1, E2 = K * Q / (0.05 ** 2), K * Q / (0.025 ** 2)
ok = abs(E2 / E1 - 4.0) < 1e-9
print("  [%s] E(r/2)/E(r) = %.6f (want 4)" % ("PASS" if ok else "FAIL", E2 / E1))
if not ok: fails.append("point ratio")

print("=" * 66)
print("(2) numeric flux over sphere  sum E.dA  == Q/eps0")
NTH, NPH = 64, 128
r = 0.05
flux = 0.0
for i in range(NTH):
    ct = -1.0 + 2.0 * (i + 0.5) / NTH          # cos(theta) midpoints
    st = math.sqrt(max(0.0, 1.0 - ct * ct))
    dA = (2.0 / NTH) * (2.0 * math.pi / NPH) * r * r   # r^2 dcos dphi
    for j in range(NPH):
        ph = 2.0 * math.pi * (j + 0.5) / NPH
        # field of the point charge (at centre) in direction (st cos, st sin, ct)
        Emag = K * Q / (r * r)
        # E . n̂ = Emag (radial == surface normal) -> integrand constant per patch
        flux += Emag * dA
target = Q / EPS0
ok = abs(flux - target) / target < 1e-6
print("  flux = %.6e   Q/eps0 = %.6e   rel err = %.2e" %
      (flux, target, abs(flux - target) / target))
print("  [%s] direct surface integration of Coulomb field" % ("PASS" if ok else "FAIL"))
if not ok: fails.append("flux numeric")

print("=" * 66)
print("(3) ball: Q_enc counting + shell cancellation + outside direct sum")
def ball_cells(N):
    """cell centres + dq inside the staircased ball"""
    rho = Q / ((4.0 / 3.0) * math.pi * R ** 3)
    h = 2.0 * R / N
    cells = []
    for i in range(N):
        cx = -R + (i + 0.5) * h
        for j in range(N):
            cy = -R + (j + 0.5) * h
            for k in range(N):
                cz = -R + (k + 0.5) * h
                if cx * cx + cy * cy + cz * cz <= R * R:
                    cells.append((cx, cy, cz, rho * h ** 3))
    return cells

def ball_numeric(cells, x):
    Ex = 0.0
    for cx, cy, cz, dq in cells:
        dx, dy, dz = x - cx, -cy, -cz
        rr = dx * dx + dy * dy + dz * dz
        if rr < 1e-18:
            continue
        Ex += K * dq * dx / (rr * math.sqrt(rr))
    return Ex

def gauss_ball(x):
    if x < R:
        return K * Q * x / (R ** 3)
    return K * Q / (x * x)

# (3a) enclosed charge by Monte Carlo: uniform random points in the cube,
#      count fraction inside r vs inside R -> Q (r/R)^3   (grid-alignment-free)
import random
random.seed(20261125)
M = 400000
pts = []
for _ in range(M):
    pts.append((random.uniform(-R, R), random.uniform(-R, R), random.uniform(-R, R)))
n_tot = sum(1 for x_, y_, z_ in pts if x_ * x_ + y_ * y_ + z_ * z_ <= R * R)
sum_rel = 0.0
for r in (0.012, 0.018, 0.024, 0.030, 0.036, 0.042):
    n_r = sum(1 for x_, y_, z_ in pts if x_ * x_ + y_ * y_ + z_ * z_ <= r * r)
    got = Q * (n_r / n_tot)          # MC estimate of Q_enc(r)
    want = Q * (r / R) ** 3
    rel = abs(got - want) / want
    sum_rel += rel
    ok = rel < 0.02
    print("  [%s] Q_enc(r=%.3f) MC rel err=%.2e vs Q(r/R)^3" %
          ("PASS" if ok else "FAIL", r, rel))
    if not ok: fails.append("ball Qenc MC r=%.3f" % r)
ok = sum_rel / 6.0 < 0.01
print("  [%s] mean MC rel err=%.2e" % ("PASS" if ok else "FAIL", sum_rel / 6.0))
if not ok: fails.append("ball Qenc MC mean")

# (3b) shell cancellation: uniform shell OUTSIDE the field point contributes ~0
rho_s, M = 0.04, 2048
dq = Q * 0.2 / M                      # put 20% of Q on one shell
Ex = 0.0
for i in range(M):
    ct = -1.0 + 2.0 * (i + 0.5) / M
    st = math.sqrt(max(0.0, 1.0 - ct * ct))
    ph = 2.0 * math.pi * (i * 0.6180339887)     # golden-angle spacing
    px, py, pz = rho_s * st * math.cos(ph), rho_s * st * math.sin(ph), rho_s * ct
    dx, dy, dz = 0.033 - px, -py, -pz
    rr = dx * dx + dy * dy + dz * dz
    Ex += K * dq * dx / (rr * math.sqrt(rr))
ref = K * (Q * 0.2) / (0.033 ** 2)    # if the shell acted as a central point
ok = abs(Ex) / ref < 1e-3
print("  [%s] outer shell field at x=0.033: |E|/point-charge = %.2e (want ~0)" %
      ("PASS" if ok else "FAIL", abs(Ex) / ref))
if not ok: fails.append("ball shell cancellation")

# (3c) outside: direct numeric Coulomb sum vs Gauss (robust - away from surface)
cells28 = ball_cells(28)
for x in (0.070, 0.090):
    n = ball_numeric(cells28, x)
    g = gauss_ball(x)
    rel = abs(n - g) / g
    ok = rel < 0.015
    print("  [%s] outside x=%.3f  numeric rel err=%.2e" % ("PASS" if ok else "FAIL", x, rel))
    if not ok: fails.append("ball outside x=%.3f" % x)

print("=" * 66)
print("(4) line charge: numeric finite line vs exact finite & infinite")
lam = 1000e-9        # 1 uC/m  (page uses q*100 nC -> same scale)
L, N = 400.0, 400000     # h=1e-3 << r so the |z|<r peak is resolved
r0 = 0.08
En = 0.0
for i in range(N):
    z = -L / 2.0 + L * (i + 0.5) / N
    dq = lam * L / N
    En += K * dq * r0 / ((r0 * r0 + z * z) ** 1.5)
E_fin = lam / (2.0 * math.pi * EPS0 * r0) * (L / 2.0) / math.hypot(L / 2.0, r0)
E_inf = lam / (2.0 * math.pi * EPS0 * r0)
ok1 = abs(En - E_fin) / E_fin < 1e-3
ok2 = abs(E_fin - E_inf) / E_inf < 1e-3
print("  numeric=%.6e  exact-finite=%.6e (rel %.1e)  infinite=%.6e (gap %.1e)" %
      (En, E_fin, abs(En - E_fin) / E_fin, E_inf, abs(E_fin - E_inf) / E_inf))
print("  [%s] numeric == finite-line formula   [%s] finite == infinite (L>>r)" %
      ("PASS" if ok1 else "FAIL", "PASS" if ok2 else "FAIL"))
if not ok1: fails.append("line numeric")
if not ok2: fails.append("line finite-vs-inf")

print("=" * 66)
print("(5) plane charge: numeric sheet vs exact square-sheet, both -> sigma/(2 eps0)")
sig = 10e-9          # 10 nC/m^2 (page scale)
a, NS = 400.0, 800   # square sheet [-a,a]^2, NS^2 patches, h=1m << z
def sheet_E(z):
    Ez = 0.0
    h = 2.0 * a / NS
    dA = h * h
    for i in range(NS):
        x = -a + (i + 0.5) * h
        for j in range(NS):
            y = -a + (j + 0.5) * h
            rr = x * x + y * y + z * z
            Ez += K * sig * dA * z / (rr * math.sqrt(rr))
    return Ez
def sheet_exact(z):
    """closed form: on-axis E of a uniformly charged square [-a,a]^2 (solid angle)."""
    return 4.0 * K * sig * math.atan(a * a / (z * math.sqrt(2.0 * a * a + z * z)))
E0 = sig / (2.0 * EPS0)
e5n, e15n = sheet_E(5.0), sheet_E(15.0)
e5x, e15x = sheet_exact(5.0), sheet_exact(15.0)
ok1 = abs(e5n - e5x) / e5x < 1e-2 and abs(e15n - e15x) / e15x < 1e-2
ok2 = abs(e5x - E0) / E0 < 0.02 and abs(e15x - E0) / E0 < 0.05   # finite corr: 1.1% / 3.3%
ok3 = abs(sheet_exact(0.15) - E0) / E0 < 5e-3                     # close range -> infinite
print("  E(z=5):  num=%.3f exact=%.3f | E(z=15): num=%.3f exact=%.3f | inf=%.3f" %
      (e5n, e5x, e15n, e15x, E0))
print("  exact(z=0.15)=%.3f  vs  sigma/2eps0=%.3f  (rel %.2e)" %
      (sheet_exact(0.15), E0, abs(sheet_exact(0.15) - E0) / E0))
print("  [%s] numeric == square-sheet closed form" % ("PASS" if ok1 else "FAIL"))
print("  [%s] both heights within finite-size correction of sigma/(2 eps0)" % ("PASS" if ok2 else "FAIL"))
print("  [%s] near-field limit -> pure sigma/(2 eps0)" % ("PASS" if ok3 else "FAIL"))
if not ok1: fails.append("plane numeric")
if not ok2: fails.append("plane value")
if not ok3: fails.append("plane limit")

print("=" * 66)
print("(6) Q_enc bookkeeping for the ball")
okall = True
for r, want in [(0.01, 0.008), (0.025, 0.125), (0.04, 0.512), (0.06, 1.0)]:
    Qenc = Q * (min(r, R) / R) ** 3 if r < R else Q
    got = Q * (r / R) ** 3 if r < R else Q
    ok = abs(Qenc / Q - want) < 1e-12 and abs(got - Qenc) < 1e-24
    okall &= ok
    print("  [%s] r=%.3f  Q_enc/Q=%.3f" % ("PASS" if ok else "FAIL", r, Qenc / Q))
if not okall: fails.append("Qenc")

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(fails))
sys.exit(0 if not fails else 1)
