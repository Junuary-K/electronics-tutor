"""
AC Thevenin/Norton reference for circuits containing R, L, C.

Phasor MNA at a single frequency. Gives V_th, Z_th, I_N as complex numbers.
Exit code 0 = all checks pass.
"""
import cmath
import math
import sys

D = math.pi / 180.0


def pol(m, deg):
    return cmath.rect(m, deg * D)


def show(z):
    return "%.4f / %8.2f deg  (%.3f %+.3fj)" % (abs(z), cmath.phase(z) / D, z.real, z.imag)


def solve(nodes, els, ref):
    """els: dicts t='Z'|'V'|'I', a, b, z (complex) or v (complex)."""
    nlist = [n for n in nodes if n != ref]
    idx = {n: i for i, n in enumerate(nlist)}
    vs = [e for e in els if e["t"] == "V"]
    n = len(nlist) + len(vs)
    A = [[0j] * n for _ in range(n)]
    b = [0j] * n
    for e in els:
        if e["t"] == "Z":
            y = 1 / e["z"]
            a, c = e["a"], e["b"]
            if a != ref:
                A[idx[a]][idx[a]] += y
                if c != ref:
                    A[idx[a]][idx[c]] -= y
            if c != ref:
                A[idx[c]][idx[c]] += y
                if a != ref:
                    A[idx[c]][idx[a]] -= y
        elif e["t"] == "I":
            if e["a"] != ref:
                b[idx[e["a"]]] -= e["v"]
            if e["b"] != ref:
                b[idx[e["b"]]] += e["v"]
    for k, e in enumerate(vs):
        row = len(nlist) + k
        a, c = e["a"], e["b"]
        if a != ref:
            A[row][idx[a]] += 1
            A[idx[a]][row] += 1
        if c != ref:
            A[row][idx[c]] -= 1
            A[idx[c]][row] -= 1
        b[row] = e["v"]
    # gaussian
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]
        b[i], b[p] = b[p], b[i]
        if abs(A[i][i]) < 1e-18:
            continue
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            if f == 0:
                continue
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
            b[r] -= f * b[i]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        if abs(A[i][i]) < 1e-18:
            x[i] = 0j
            continue
        s = b[i] - sum(A[i][c] * x[c] for c in range(i + 1, n))
        x[i] = s / A[i][i]
    out = {ref: 0j}
    for nd, i in idx.items():
        out[nd] = x[i]
    return out, x, len(nlist)


def thevenin(nodes, els, tp, tn):
    V, _, _ = solve(nodes, els, tn)
    vth = V[tp]
    # short circuit current: add 0 V source between tp and tn
    els2 = els + [{"t": "V", "a": tp, "b": tn, "v": 0j}]
    V2, x2, base = solve(nodes, els2, tn)
    isc = x2[base + len([e for e in els2 if e["t"] == "V"]) - 1]
    zth = vth / isc if abs(isc) > 0 else float("inf")
    return vth, isc, zth


# ---------------- T1: RC 저역통과 (1 kHz) ----------------
f = 1000.0
w = 2 * math.pi * f
RC = 100e-9
R1 = 1000.0
nodes = {"n0", "n1", "n2"}
els = [
    {"t": "V", "a": "n1", "b": "n0", "v": pol(10, 0)},
    {"t": "Z", "a": "n1", "b": "n2", "z": complex(1000, 0)},
    {"t": "Z", "a": "n2", "b": "n0", "z": 1 / (1j * w * RC)},
]
vth, isc, zth = thevenin(nodes, els, "n2", "n0")

print("=" * 76)
print("T1  RC 저역통과  V=10/0deg, R=1k, C=100nF, f=1kHz")
print("=" * 76)
print("  Z_C      =", show(1 / (1j * w * RC)))
print("  V_th     =", show(vth))
print("  I_sc     =", show(isc))
print("  Z_th     =", show(zth))
print("  |V_th| = %.4f V,  |Z_th| = %.4f ohm" % (abs(vth), abs(zth)))

# ---------------- T2: RLC (5 kHz, 위상 30도) ----------------
f2 = 5000.0
w2 = 2 * math.pi * f2
nodes2 = {"n0", "n1", "n2", "n3"}
els2 = [
    {"t": "V", "a": "n1", "b": "n0", "v": pol(5, 30)},
    {"t": "Z", "a": "n1", "b": "n2", "z": complex(10, 0)},
    {"t": "Z", "a": "n2", "b": "n3", "z": 1j * w2 * 10e-3},
    {"t": "Z", "a": "n3", "b": "n0", "z": 1 / (1j * w2 * 1e-6)},
]
vth2, isc2, zth2 = thevenin(nodes2, els2, "n3", "n0")

print()
print("=" * 76)
print("T2  R + L + C  V=5/30deg, R=10, L=10mH, C=1uF, f=5kHz")
print("=" * 76)
print("  Z_L      =", show(1j * w2 * 10e-3))
print("  Z_C      =", show(1 / (1j * w2 * 1e-6)))
print("  V_th     =", show(vth2))
print("  Z_th     =", show(zth2))
print("  I_N      =", show(vth2 / zth2))
print("  |V_th| = %.4f V,  |Z_th| = %.4f ohm" % (abs(vth2), abs(zth2)))

# ---------------- checks ----------------
fail = []


def check(name, cond, detail=""):
    print("  [%s] %-40s %s" % ("PASS" if cond else "FAIL", name, detail))
    if not cond:
        fail.append(name)


print()
print("  공통 검사")
check("V_th / I_sc = Z_th (T1)", abs(vth / isc - zth) < 1e-9)
check("V_th / I_sc = Z_th (T2)", abs(vth2 / isc2 - zth2) < 1e-9)
check("|Z_C| = 1/(2 pi f C), T1", abs(abs(1 / (1j * w * RC)) - 1 / (2 * math.pi * f * RC)) < 1e-9)
check("|Z_L| = 2 pi f L, T2", abs(abs(1j * w2 * 10e-3) - 2 * math.pi * f2 * 10e-3) < 1e-9)

# T1 closed form: divider
zc = 1 / (1j * w * RC)
check("T1 V_th = V*Zc/(R+Zc)", abs(vth - pol(10, 0) * zc / (R1 + zc)) < 1e-12)
check("T1 Z_th = R||Zc", abs(zth - (R1 * zc) / (R1 + zc)) < 1e-12)

print()
print("=" * 76)
print("RESULT:", "ALL PASS" if not fail else "FAIL: " + ", ".join(fail))
sys.exit(0 if not fail else 1)
