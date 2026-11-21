"""
MOSFET (n-channel, 1st-order square law + channel-length modulation) reference.

Checks the numbers our HTML app must reproduce, and cross-checks the
region boundary continuity and the small-signal parameters two ways.

Exit code 0 = all pass.
"""
import sys

# device under test
VT0 = 0.7     # V
gamma = 0.5   # sqrt(V)
phi2 = 0.6    # 2*phi_F, V
kn = 120e-6   # A/V^2
WL = 10.0
lam = 0.02    # 1/V

fail = []


def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print("  [%s] %-44s got %.10g  want %.10g" % ("PASS" if ok else "FAIL", name, got, want))
    if not ok:
        fail.append(name)


def vth(vbs):
    return VT0 + gamma * ((phi2 + max(0.0, -vbs)) ** 0.5 - phi2 ** 0.5)


def idrain(vgs, vds, vbs=0.0):
    vt = vth(vbs)
    vov = vgs - vt
    beta = kn * WL
    if vov <= 0:
        return 0.0, 'cutoff', vt, vov
    if vds < vov:
        i = beta * (vov * vds - 0.5 * vds * vds) * (1 + lam * vds)
        return i, 'triode', vt, vov
    i = 0.5 * beta * vov * vov * (1 + lam * vds)
    return i, 'sat', vt, vov


print("=" * 74)
print("n-MOSFET  VT0=%.2f V, gamma=%.2f, 2phiF=%.1f, kn'=%.0f uA/V^2, W/L=%.1f, lambda=%.3f"
      % (VT0, gamma, phi2, kn * 1e6, WL, lam))
print("=" * 74)

cases = [
    (0.5, 1.0, "V_GS below V_T -> cutoff"),
    (2.0, 0.5, "linear region"),
    (2.0, 2.0, "saturation"),
    (1.2, 0.2, "near threshold, small V_DS"),
]
for vgs, vds, note in cases:
    i, reg, vt, vov = idrain(vgs, vds)
    print("  V_GS=%.2f V_DS=%.2f -> %-7s  V_T=%.4f V_ov=%+.4f  I_D=%.6e A   (%s)"
          % (vgs, vds, reg, vt, vov, i, note))

print()
print("  boundary continuity at V_DS = V_ov (V_GS=2.0):")
vt = vth(0.0)
vov = 2.0 - vt
b = kn * WL
i_lin = b * (vov * vov - 0.5 * vov * vov) * (1 + lam * vov)
i_sat = 0.5 * b * vov * vov * (1 + lam * vov)
check("I_D(linear) == I_D(sat) at V_DS=V_ov", i_lin, i_sat)

print()
print("  small-signal, V_GS=2.0 V, V_DS=2.0 V (saturation):")
i0, _, _, vov0 = idrain(2.0, 2.0)
gm_analytic = b * vov0 * (1 + lam * 2.0)
h = 1e-6
gm_num = (idrain(2.0 + h, 2.0)[0] - idrain(2.0 - h, 2.0)[0]) / (2 * h)
check("g_m analytic vs numeric", gm_analytic, gm_num, 1e-6)
ro_analytic = (1 + lam * 2.0) / (lam * i0)   # exact: r_o=(1+lam*V_DS)/(lam*I_D)
ro_num = 1.0 / ((idrain(2.0, 2.0 + 1e-4)[0] - idrain(2.0, 2.0 - 1e-4)[0]) / 2e-4)
check("r_o analytic vs numeric", ro_analytic, ro_num, 1e-4)
print("     g_m = %.4f mA/V,  r_o = %.2f kohm" % (gm_analytic * 1e3, ro_analytic / 1e3))

print()
print("  body effect (V_BS = -1 V):")
check("V_T(V_BS=-1) from gamma law", vth(-1.0),
      VT0 + gamma * ((0.6 + 1.0) ** 0.5 - 0.6 ** 0.5))
check("V_T(V_BS=0) == V_T0", vth(0.0), VT0)

print()
print("=" * 74)
print("RESULT:", "ALL PASS" if not fail else "FAIL: " + ", ".join(fail))
sys.exit(0 if not fail else 1)
