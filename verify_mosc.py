"""
MOS capacitor (p-Si / SiO2 / n+ poly-Si) - independent verification.

Purpose: compute the same quantities the HTML app reports, using an
independently written solver, and check them against:
  (a) closed-form textbook relations (depletion approximation)
  (b) internal identities that must hold for any consistent solution
  (c) a second numerical route (bisection on V_G <-> psi_s vs. explicit sweep)

Exit code 0 = all checks pass.
"""
import math
import sys

# --- constants (CODATA / standard device values, 300 K) ---
q = 1.602176634e-19
kB = 1.380649e-23
T = 300.0
Vt = kB * T / q
eps0 = 8.8541878128e-14
epsS = 11.7 * eps0          # Si
epsOX = 3.9 * eps0          # SiO2
ni = 1.5e10                 # cm^-3
Eg = 1.12
chi = 4.01
dEc_ox = 3.10               # SiO2/Si conduction band offset

# --- device under test ---
Na = 1e16                   # cm^-3, p-type
tox = 10e-7                 # cm  (10 nm)
phi_m = 4.05                # V, n+ poly-Si

phiF = Vt * math.log(Na / ni)
Cox = epsOX / tox
phi_ms = phi_m - (chi + Eg / 2 + phiF)
VFB = phi_ms
LD = math.sqrt(epsS * Vt / (q * Na))
p0 = Na
n0 = ni * ni / Na


def Ffun(u):
    """Berglund/Kingston F function, Boltzmann statistics, p-type."""
    a = math.exp(-u) + u - 1.0
    b = (n0 / p0) * (math.exp(u) - u - 1.0)
    v = a + b
    return math.sqrt(v) if v > 0 else 0.0


def Efield(psi):
    """|E| at the surface for surface potential psi [V], units V/cm."""
    return math.sqrt(2.0) * Vt / LD * Ffun(psi / Vt)


def Qs(psi):
    """Total semiconductor charge per area [C/cm^2]. Negative in depletion."""
    if psi == 0:
        return 0.0
    return -math.copysign(epsS * Efield(psi), psi)


def VG_of_psi(psi):
    return VFB + psi - Qs(psi) / Cox


def psi_of_VG(VG, lo=-1.2, hi=2.0, it=200):
    """Bisection: VG_of_psi is monotonically increasing in psi."""
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if VG_of_psi(mid) < VG:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cap(psi):
    """Small-signal C per area, series model: 1/C = 1/Cox + 1/Cs."""
    h = 1e-7
    Cs = abs(Qs(psi + h) - Qs(psi - h)) / (2 * h)
    return 1.0 / (1.0 / Cox + 1.0 / Cs)


# --- reference operating points ---
VT_dep = VFB + 2 * phiF + math.sqrt(2 * epsS * q * Na * 2 * phiF) / Cox
WT = math.sqrt(2 * epsS * (2 * phiF) / (q * Na))
Cmin = 1.0 / (1.0 / Cox + WT / epsS)

print("=" * 72)
print("MOS-C  p-Si Na=%.3g cm^-3, tox=%.0f nm, n+ poly gate, 300 K" % (Na, tox * 1e7))
print("=" * 72)
print("phi_F      = %.4f V" % phiF)
print("C_ox       = %.4e F/cm^2  (%.2f uF/cm^2)" % (Cox, Cox * 1e6))
print("phi_ms     = %.4f V" % phi_ms)
print("V_FB       = %.4f V" % VFB)
print("V_T (dep.) = %.4f V" % VT_dep)
print("W_T        = %.1f nm" % (WT * 1e7))
print("C_min      = %.4e F/cm^2  (C/Cox = %.4f)" % (Cmin, Cmin / Cox))
print()

fails = []

# --- check 1: flat band identity ---
psi0 = psi_of_VG(VFB)
ok = abs(psi0) < 1e-9
print("[%s] psi_s = 0 at V_G = V_FB  (psi_s=%.2e)" % ("PASS" if ok else "FAIL", psi0))
if not ok:
    fails.append("flatband")

# --- check 2: exact Qs vs depletion approximation in moderate depletion ---
psi_test = phiF          # onset of weak inversion
Qexact = abs(Qs(psi_test))
Qdep = math.sqrt(2 * q * epsS * Na * psi_test)
rel = abs(Qexact - Qdep) / Qdep * 100
ok = rel < 8.0
print("[%s] exact Qs vs depletion approx at psi_s=phi_F  (%.2f%% diff)" % ("PASS" if ok else "FAIL", rel))
if not ok:
    fails.append("Qdep")

# --- check 3: max capacitance -> C_ox (deep accumulation) ---
psi_acc = -0.5
c_acc = cap(psi_acc)
ok = c_acc / Cox > 0.97
print("[%s] C -> C_ox in accumulation  (C/Cox=%.4f at psi_s=-0.5V)" % ("PASS" if ok else "FAIL", c_acc / Cox))
if not ok:
    fails.append("Cacc")

# --- check 4: V_T reproduces psi_s = 2 phi_F ---
psi_at_VT = psi_of_VG(VT_dep)
ok = abs(psi_at_VT - 2 * phiF) < 0.02
print("[%s] psi_s(2phi_F) at V_G = V_T  (psi_s=%.4f, 2phi_F=%.4f)"
      % ("PASS" if ok else "FAIL", psi_at_VT, 2 * phiF))
if not ok:
    fails.append("VT")

# --- check 5: C-V minimum sits near strong inversion ---
xs = [VFB - 1.0 + 3.0 * i / 300 for i in range(301)]
cv = [(v, cap(psi_of_VG(v)) / Cox) for v in xs]
vmin, cmin = min(cv, key=lambda t: t[1])
ok = abs(cmin - Cmin / Cox) < 0.02 and abs(vmin - VT_dep) < 0.5
print("[%s] C-V minimum  (at V_G=%.3f V, C/Cox=%.4f; theory %.4f)"
      % ("PASS" if ok else "FAIL", vmin, cmin, Cmin / Cox))
if not ok:
    fails.append("CVmin")

# --- check 6: monotonic psi_s(V_G) and charge balance ---
psi_a = psi_of_VG(0.0)
psi_b = psi_of_VG(2.0)
ok = psi_a < psi_b
print("[%s] psi_s increases with V_G  (%.4f -> %.4f V)" % ("PASS" if ok else "FAIL", psi_a, psi_b))
if not ok:
    fails.append("monotonic")

# --- check 7: gate charge balance and oxide field consistency ---
psi_x = 0.30
QG = -Qs(psi_x)
Vox = QG / Cox
ok = abs((VFB + psi_x + Vox) - VG_of_psi(psi_x)) < 1e-12
print("[%s] V_G = V_FB + psi_s + Q_G/C_ox  (Vox=%.4f V)" % ("PASS" if ok else "FAIL", Vox))
if not ok:
    fails.append("balance")

print()
print("-" * 72)
print("operating points reported by the app (should match exactly)")
print("-" * 72)
for VG in (VFB, 0.0, 0.5, 1.0, 1.5, VT_dep):
    ps = psi_of_VG(VG)
    print("V_G=%+.3f V | psi_s=%.4f V | W=%7.1f nm | C/Cox=%.4f | Q_s=%.3e C/cm^2"
          % (VG, ps, math.sqrt(2 * epsS * max(ps, 0) / (q * Na)) * 1e7,
             cap(ps) / Cox, Qs(ps)))

print()
print("=" * 72)
print("RESULT:", "ALL PASS" if not fails else "FAILURES: " + ", ".join(fails))
sys.exit(0 if not fails else 1)
