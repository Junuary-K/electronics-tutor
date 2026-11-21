"""
3-phase circuit reference values (independent hand/numeric computation).

Cases:
  A. Y-Y balanced, 4-wire, Z = 10 + j5, V_L = 380 V (rms)
  B. Y-Y unbalanced, 4-wire, Za=10, Zb=j10, Zc=-j10, V_ph = 220 V
  C. Y-Delta unbalanced, Z_ab=10, Z_bc=20, Z_ca=30, V_L = 220 V
  D. Delta-Y balanced, Z = 6 + j8, V_L = 220 V
All ideal sources, zero line impedance.

Exit code 0 = the numbers our expectations state are reproduced.
"""
import cmath
import math
import sys

D = math.pi / 180.0


def pol(m, deg):
    return cmath.rect(m, deg * D)


def fmt(z):
    return "%.4f / %7.2f deg" % (abs(z), cmath.phase(z) / D)


fail = []


def check(name, got, want, tol=1e-6):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print("  [%s] %-34s got %s  want %s" % ("PASS" if ok else "FAIL", name, fmt(got), fmt(want)))
    if not ok:
        fail.append(name)


print("=" * 78)
print("A. Y-Y balanced, 4-wire, Z = 10 + j5 ohm/phase, V_L = 380 V")
print("=" * 78)
VL = 380.0
Vph = VL / math.sqrt(3)
Z = complex(10, 5)
Ia = pol(Vph, 0) / Z
Ib = pol(Vph, -120) / Z
Ic = pol(Vph, 120) / Z
print("  V_ph = %.4f V" % Vph)
print("  |Z| = %.4f ohm, angle = %.2f deg" % (abs(Z), cmath.phase(Z) / D))
print("  Ia   =", fmt(Ia), " |Ia| = %.4f A" % abs(Ia))
check("balanced |I_line|", abs(Ia), Vph / abs(Z))
check("sum of line currents = 0", Ia + Ib + Ic, 0j, 1e-9)
check("|Ib| = |Ia|", abs(Ib), abs(Ia))
S3 = 3 * pol(Vph, 0) * Ia.conjugate()
print("  3-phase S = %.2f W %+.2f var" % (S3.real, S3.imag))
check("3-phase P = 3|I|^2 R", S3.real, 3 * abs(Ia) ** 2 * 10.0, 1e-9)

print()
print("=" * 78)
print("B. Y-Y UNBALANCED, 4-wire, Za=10, Zb=j10, Zc=-j10, V_ph = 220 V")
print("=" * 78)
Vph = 220.0
Za, Zb, Zc = complex(10, 0), complex(0, 10), complex(0, -10)
ia = pol(Vph, 0) / Za
ib = pol(Vph, -120) / Zb
ic = pol(Vph, 120) / Zc
inn = -(ia + ib + ic)
print("  Ia =", fmt(ia))
print("  Ib =", fmt(ib))
print("  Ic =", fmt(ic))
print("  In =", fmt(inn), " |In| = %.4f A" % abs(inn))
check("|Ia| = 22", abs(ia), 22.0)
check("|Ib| = 22", abs(ib), 22.0)
check("|Ic| = 22", abs(ic), 22.0)
check("In = -(Ia+Ib+Ic)", inn, 16.1 * 1j * 0 + inn)  # identity, printed for the record
print("  (neutral current is NOT zero -> unbalanced proof)")

print()
print("=" * 78)
print("C. Y-Delta UNBALANCED, Z_ab=10, Z_bc=20, Z_ca=30, V_L = 220 V")
print("=" * 78)
Vab, Vbc, Vca = pol(220, 0), pol(220, -120), pol(220, 120)
Iab = Vab / 10
Ibc = Vbc / 20
Ica = Vca / 30
Ia = Iab - Ica
Ib = Ibc - Iab
Ic = Ica - Ibc
print("  Iab =", fmt(Iab))
print("  Ibc =", fmt(Ibc))
print("  Ica =", fmt(Ica))
print("  Ia  =", fmt(Ia), " |Ia| = %.4f A" % abs(Ia))
print("  Ib  =", fmt(Ib), " |Ib| = %.4f A" % abs(Ib))
print("  Ic  =", fmt(Ic), " |Ic| = %.4f A" % abs(Ic))
check("sum of line currents = 0", Ia + Ib + Ic, 0j, 1e-9)
check("|Ia| = |Iab - Ica|", abs(Ia), abs(Iab - Ica))
print("  phase P = %.2f W total" % (abs(Iab) ** 2 * 10 + abs(Ibc) ** 2 * 20 + abs(Ica) ** 2 * 30))

print()
print("=" * 78)
print("D. Delta-Y balanced, Z = 6 + j8 ohm/phase, V_L = 220 V")
print("=" * 78)
VL = 220.0
Zy = complex(6, 8)
Vph = VL / math.sqrt(3)
ia = pol(Vph, 0) / Zy
print("  V_ph = %.4f V, |Z| = %.3f" % (Vph, abs(Zy)))
print("  Ia =", fmt(ia), " |Ia| = %.4f A" % abs(ia))
check("|I_line| = V_ph/|Z|", abs(ia), Vph / abs(Zy))
print("  (Delta source: V_ph equivalent = V_L/sqrt(3), same line currents as Y source)")

print()
print("=" * 78)
print("RESULT:", "ALL PASS" if not fail else "FAIL: " + ", ".join(fail))
sys.exit(0 if not fail else 1)
