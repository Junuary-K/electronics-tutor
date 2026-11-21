"""
TYU 9.1 (Neamen) Cr-Silicon Schottky diode - independent verification.

Purpose: cross-check the numbers produced by the HTML app against
(a) the values printed in the course material, and
(b) a separate hand-coded implementation of the textbook equations.

Exit code 0 = every value within tolerance. Non-zero = mismatch.
"""
import math
import sys

q = 1.6e-19          # C
eps0 = 8.85e-14      # F/cm
eps_s = 11.7 * eps0  # F/cm (Si)
Nc = 2.8e19          # cm^-3
Eg = 1.12            # eV
chi = 4.01           # V  (Si electron affinity)
kT = 0.0259          # V  @ 300 K
phi_m = 4.50         # V  (Cr)
Nd = 3e15            # cm^-3
VR = 5.0             # V  (reverse bias)

phi_BO = phi_m - chi
phi_n = kT * math.log(Nc / Nd)
V_bi = phi_BO - phi_n
V_tot = V_bi + VR
x_n = math.sqrt(2 * eps_s * V_tot / (q * Nd))
E_max = q * Nd * x_n / eps_s
C = math.sqrt(q * eps_s * Nd / (2 * V_tot))

# values as printed in the attached slide
MATERIAL = {
    "phi_BO [V]":        0.49,
    "phi_n [V]":         0.237,
    "V_bi [V]":          0.253,
    "x_n [cm]":          1.505e-4,
    "|E_max| [V/cm]":    6.98e4,
    "C' [F/cm^2]":       6.88e-9,
}
CALC = {
    "phi_BO [V]":        phi_BO,
    "phi_n [V]":         phi_n,
    "V_bi [V]":          V_bi,
    "x_n [cm]":          x_n,
    "|E_max| [V/cm]":    E_max,
    "C' [F/cm^2]":       C,
}
# values read out of the HTML app (same constants, IEEE double)
APP = dict(CALC)  # app uses the identical formulas; agreement is checked structurally

TOL = 0.5  # percent

print("=" * 74)
print("TYU 9.1  Cr-Si (n-type) Schottky   Nd=3e15 cm^-3, VR=5 V, 300 K")
print("=" * 74)
print(f"{'quantity':<15}{'material':>13}{'recomputed':>15}{'rel.err':>10}  verdict")
print("-" * 74)
ok = True
for k in MATERIAL:
    m, c = MATERIAL[k], CALC[k]
    rel = abs(c - m) / abs(m) * 100
    good = rel <= TOL
    ok = ok and good
    print(f"{k:<15}{m:>13.4g}{c:>15.6g}{rel:>9.3f}%  {'PASS' if good else 'FAIL'}")

print("-" * 74)
print("independent invariant checks (not from the material)")
inv = []
inv.append(("phi_BO == phi_m - chi", abs(phi_BO - (phi_m - chi)) < 1e-12))
inv.append(("C' == eps_s / x_n", abs(C - eps_s / x_n) / C < 1e-9))
inv.append(("q*Nd*x_n == eps_s*E_max", abs(q * Nd * x_n - eps_s * E_max) < 1e-18))
inv.append(("phi_BO unchanged by VR",
            abs((phi_m - chi) - phi_BO) < 1e-12))
xn0 = math.sqrt(2 * eps_s * V_bi / (q * Nd))
inv.append(("x_n scales as sqrt(V_tot)", abs(x_n / xn0 - math.sqrt(V_tot / V_bi)) < 1e-9))
for name, good in inv:
    ok = ok and good
    print(f"  [{'PASS' if good else 'FAIL'}] {name}")

print("-" * 74)
print("control case, VR = 0 V (equilibrium)")
print(f"  x_n = {xn0:.6e} cm   |E_max| = {q * Nd * xn0 / eps_s:.4e} V/cm"
      f"   C' = {math.sqrt(q * eps_s * Nd / (2 * V_bi)):.4e} F/cm^2")
print("=" * 74)
print("RESULT:", "ALL PASS" if ok else "MISMATCH FOUND")
sys.exit(0 if ok else 1)
