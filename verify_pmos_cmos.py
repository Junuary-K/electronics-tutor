"""
p-channel MOSFET and CMOS inverter reference values.

Model: 1st-order square law + channel-length modulation, magnitudes for pMOS.
Exit code 0 = all checks pass.
"""
import sys

VDD = 3.3
VTN = 0.7
VTP = -0.7      # negative for pMOS
KN = 120e-6     # A/V^2
KP = 60e-6      # A/V^2  (per-unit k_p', CMOS section uses it directly as beta)
LAM = 0.02

fail = []


def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print("  [%s] %-46s got %.10g  want %.10g" % ("PASS" if ok else "FAIL", name, got, want))
    if not ok:
        fail.append(name)


def idn(vgs, vds):
    vov = vgs - VTN
    if vov <= 0:
        return 0.0
    if vds < vov:
        return KN * (vov * vds - 0.5 * vds * vds) * (1 + LAM * vds)
    return 0.5 * KN * vov * vov * (1 + LAM * vds)


def idp(vsg, vsd):
    vov = vsg - abs(VTP)
    if vov <= 0:
        return 0.0
    if vsd < vov:
        return KP * (vov * vsd - 0.5 * vsd * vsd) * (1 + LAM * vsd)
    return 0.5 * KP * vov * vov * (1 + LAM * vsd)


print("=" * 76)
print("p-MOSFET  V_TP=%.2f V, beta_p = k_p'*W/L = %.0f uA/V^2, lambda=%.3f" % (VTP, KP * 10 * 1e6, LAM))
print("=" * 76)
for vsg, vsd, note in [(0.5, 1.0, "below |V_TP| -> cutoff"),
                       (1.5, 0.3, "linear"),
                       (1.5, 1.5, "saturation")]:
    i = idp(vsg, vsd) * 10   # app uses W/L = 10 in the pMOS panel
    vov = vsg - abs(VTP)
    reg = "cutoff" if vov <= 0 else ("linear" if vsd < vov else "saturation")
    print("  V_SG=%.1f V_SD=%.1f -> %-10s V_ov=%.3f  I_SD=%.6e A   (%s)"
          % (vsg, vsd, reg, vov, i, note))

print()
print("=" * 76)
print("CMOS inverter  V_DD=%.1f, V_TN=%.2f, V_TP=%.2f, k_n=%.0f, k_p=%.0f"
      % (VDD, VTN, VTP, KN * 1e6, KP * 1e6))
print("=" * 76)


def vout_of(vin):
    """Solve I_n(vin,vout) = I_p(vin,vout) by bisection."""
    if vin <= VTN:
        return VDD          # nMOS off
    if (VDD - vin) <= abs(VTP):
        return 0.0          # pMOS off
    f = lambda v: idn(vin, v) - idp(VDD - vin, VDD - v)
    lo, hi = 0.0, VDD
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print("  V_in -> V_out")
for vin in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.3):
    print("    %.2f V -> %.4f V" % (vin, vout_of(vin)))

check("V_out = V_DD when V_in = 0", vout_of(0.0), VDD)
check("V_out = 0 when V_in = V_DD", vout_of(VDD), 0.0)

# switching threshold
lo, hi = 0.0, VDD
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if vout_of(mid) - mid > 0:
        lo = mid
    else:
        hi = mid
VM = 0.5 * (lo + hi)
print("  V_M (V_in = V_out) = %.4f V" % VM)
check("V_M satisfies V_out(V_M) = V_M", vout_of(VM), VM, 1e-6)

# current balance at V_M
i_n = idn(VM, VM)
i_p = idp(VDD - VM, VDD - VM)
print("  at V_M: I_n = %.4e A, I_p = %.4e A, P_static = %.4e W" % (i_n, i_p, VDD * i_n))
check("I_n = I_p at V_M", i_n, i_p, 1e-6)

# rails: essentially no current
check("I_DD ~ 0 at V_in = 0", idp(VDD - 0.0, VDD - VDD), 0.0, 1e-12)
check("I_DD ~ 0 at V_in = V_DD", idn(VDD, 0.0), 0.0, 1e-12)

print()
print("=" * 76)
print("RESULT:", "ALL PASS" if not fail else "FAIL: " + ", ".join(fail))
sys.exit(0 if not fail else 1)
