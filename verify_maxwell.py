"""
Maxwell's equations demos - independent verification.

Mirrors the four numeric demos in maxwell.html:
  (1) Gauss E: patch-sum flux through a sphere around a point charge == Q/eps0,
      and the same value at two different radii (r-independence).
  (2) Gauss B: two-monopole bar-magnet model; total flux through an enclosing
      sphere == 0 (in and out cancel). Also: a single off-center monopole's
      flux through the sphere is still 4*pi*k*q_m (solid-angle invariance).
  (3) Faraday: closed-form EMF = -A*w*B0*cos(wt) vs central-difference
      -d(A*B)/dt at several phases; Lenz sign; amplitude A*w*B0.
  (4) Ampere-Maxwell: conduction current C*dV/dt vs displacement current
      eps0*A*dE/dt (C = eps0*A/d, V = k*t), matched numerically via a
      finite difference of the E-flux.
  (5) Wave speed implied by the four equations: c = 1/sqrt(mu0*eps0).

Exit code 0 = all pass.
"""
import math, sys

EPS0 = 8.8541878128e-12
MU0 = 4.0 * math.pi * 1e-7
K = 1.0 / (4.0 * math.pi * EPS0)

fails = []

print("=" * 66)
print("(1) Gauss E: numeric surface flux == Q/eps0, r-independent")
Q = 10e-9
def flux_point_clean(Q, r, nth=64, nph=128):
    f = 0.0
    for i in range(nth):
        dcos = 2.0 / nth
        dphi = 2.0 * math.pi / nph
        dA = dcos * dphi * r * r
        Emag = K * Q / (r * r)
        for j in range(nph):
            f += Emag * dA
    return f
target = Q / EPS0
for r in (0.03, 0.07, 0.12):
    f = flux_point_clean(Q, r)
    rel = abs(f - target) / target
    ok = rel < 1e-6
    print("  [%s] r=%.2f m  flux=%.4f  rel err=%.2e" % ("PASS" if ok else "FAIL", r, f, rel))
    if not ok: fails.append("gaussE r=%.2f" % r)
f1, f2 = flux_point_clean(Q, 0.05), flux_point_clean(Q, 0.02)
ok = abs(f1 - f2) < 1e-6 * abs(target)
print("  [%s] r-independence: flux(0.05) - flux(0.02) = %.2e" % ("PASS" if ok else "FAIL", abs(f1 - f2)))
if not ok: fails.append("gaussE r-indep")

print("=" * 66)
print("(2) Gauss B: two-monopole magnet -> enclosing flux = 0")
def flux_two_poles(sep, rs, nth=64, nph=128):
    """+q_m at z=+sep/2, -q_m at z=-sep/2 (magnetic Coulomb, unit q_m).
    Integrand is the RADIAL component B.n̂ = K*(rs²-s z)*rs/|P-S|³ (not B_z)."""
    f = 0.0
    dcos = 2.0 / nth
    dphi = 2.0 * math.pi / nph
    for i in range(nth):
        ct = -1.0 + 2.0 * (i + 0.5) / nth
        st = math.sqrt(max(0.0, 1.0 - ct * ct))
        for j in range(nph):
            ph = 2.0 * math.pi * (j + 0.5) / nph
            x = rs * st * math.cos(ph); y = rs * st * math.sin(ph); z = rs * ct
            fz = 0.0
            for s, sg in ((sep / 2.0, 1.0), (-sep / 2.0, -1.0)):
                dx, dy, dz = x, y, z - s
                rr = dx * dx + dy * dy + dz * dz
                dot = rs * rs - s * z          # (P-S).P
                fz += sg * K * dot * rs / (rr * math.sqrt(rr))
            f += fz * dcos * dphi
    return f
def flux_one_pole(s, rs, nth=64, nph=128):
    f = 0.0
    dcos = 2.0 / nth
    dphi = 2.0 * math.pi / nph
    for i in range(nth):
        ct = -1.0 + 2.0 * (i + 0.5) / nth
        st = math.sqrt(max(0.0, 1.0 - ct * ct))
        for j in range(nph):
            ph = 2.0 * math.pi * (j + 0.5) / nph
            x = rs * st * math.cos(ph); y = rs * st * math.sin(ph); z = rs * ct
            dx, dy, dz = x, y, z - s
            rr = dx * dx + dy * dy + dz * dz
            dot = rs * rs - s * z
            f += K * dot * rs / (rr * math.sqrt(rr)) * dcos * dphi
    return f
for sep, rs in [(0.04, 0.12), (0.10, 0.12), (0.12, 0.20)]:
    f = flux_two_poles(sep, rs)
    scale = 4.0 * math.pi * K          # one pole's total flux (unit q_m)
    ok = abs(f) / scale < 1e-6
    print("  [%s] sep=%.2f rs=%.2f  total flux=%.3e  (|.|/single=%.1e)" %
          ("PASS" if ok else "FAIL", sep, rs, f, abs(f) / scale))
    if not ok: fails.append("gaussB %.2f/%.2f" % (sep, rs))
# solid-angle invariance: one off-center pole's TOTAL flux is still 4*pi*K,
# but the patch sum only approximates it (integrand no longer ring-constant),
# so the tolerance reflects the midpoint-rule discretization, not the physics
f_off = flux_one_pole(0.05, 0.12)
rel = abs(f_off - 4 * math.pi * K) / (4 * math.pi * K)
ok = rel < 2e-3
print("  [%s] single off-center pole total flux (rel %.1e, discretization-bound)" % ("PASS" if ok else "FAIL", rel))
if not ok: fails.append("pole invariance")

print("=" * 66)
print("(3) Faraday: EMF closed form vs numeric -dPhi/dt")
def emf_closed(A, w, B0, t):
    return -A * w * B0 * math.cos(w * t)
worst = 0.0
for A, w, B0 in [(25e-4, 6.0, 20e-3), (50e-4, 2.0, 10e-3), (5e-4, 20.0, 5e-3)]:
    for t in (0.0, 0.13, 0.37, 0.62, 0.95):
        h = 1e-5
        num = -(A * (B0 * math.sin(w * (t + h)) - B0 * math.sin(w * (t - h))) / (2 * h))
        cf = emf_closed(A, w, B0, t)
        rel = abs(num - cf) / max(abs(cf), 1e-15)
        worst = max(worst, rel)
ok = worst < 1e-6
print("  [%s] central difference vs closed form, 15 points, worst rel=%.2e" % ("PASS" if ok else "FAIL", worst))
if not ok: fails.append("faraday numeric")
# Lenz sign: EMF > 0 exactly when dB/dt < 0
okall = True
for w, B0 in [(6.0, 20e-3)]:
    for t in [i * 0.1 for i in range(63)]:
        dB = B0 * w * math.cos(w * t)
        emf = -25e-4 * dB
        if not ((emf > 0) == (dB < 0) or abs(dB) < 1e-12):
            okall = False
print("  [%s] Lenz: sign(EMF) = -sign(dB/dt) across one period" % ("PASS" if okall else "FAIL"))
if not okall: fails.append("lenz sign")
# amplitude
A, w, B0 = 25e-4, 6.0, 20e-3
amp = max(abs(emf_closed(A, w, B0, i * 0.01)) for i in range(1571))
ok = abs(amp - A * w * B0) < 1e-15
print("  [%s] amplitude = A*w*B0 = %.4f mV" % ("PASS" if ok else "FAIL", amp * 1e3))
if not ok: fails.append("faraday amplitude")

print("=" * 66)
print("(4) Ampere-Maxwell: I_cond == I_d, both by two formulas + finite diff")
for A, d, kk in [(50e-4, 1e-3, 50.0), (100e-4, 2e-3, 10.0), (10e-4, 5e-3, 100.0)]:
    C = EPS0 * A / d
    I_cond = C * kk
    I_d = EPS0 * A * (kk / d)
    # numeric: Phi_E(t) = E(t)*A = (k t / d) * A ; I_d = eps0 * dPhi_E/dt
    h = 1e-6
    dPhi = ((kk * (0.5 + h) / d) * A - (kk * (0.5 - h) / d) * A) / (2 * h)
    I_d_num = EPS0 * dPhi
    ok = abs(I_cond - I_d) < 1e-15 and abs(I_d - I_d_num) < 1e-15 and abs(I_cond - I_d_num) < 1e-12
    print("  [%s] A=%dcm2 d=%.1fmm k=%d: I=%.4e A  I_d=%.4e A  I_d(num)=%.4e" %
          ("PASS" if ok else "FAIL", A * 1e4, d * 1e3, kk, I_cond, I_d, I_d_num))
    if not ok: fails.append("ampere A=%.0e" % A)

print("=" * 66)
print("(5) the four equations imply light: c = 1/sqrt(mu0*eps0)")
c = 1.0 / math.sqrt(MU0 * EPS0)
ok = abs(c - 2.99792458e8) < 1.0
print("  c = %.4f m/s  (definition-exact 299792458)" % c)
print("  [%s] matches the defined light speed within 1 m/s" % ("PASS" if ok else "FAIL"))
if not ok: fails.append("lightspeed")

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(sorted(set(fails))))
sys.exit(0 if not fails else 1)
