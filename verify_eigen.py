"""
Eigenvalues / eigenvectors of 2x2 matrices - independent verification.

Checks (must match eigen.html):
  (1) closed-form roots of lam^2 - tr*lam + det substituted back into the
      characteristic polynomial
  (2) eigenvector residual ||Av - lam v|| / ||v|| using the same row-choice
      rule as the page (pick the longer row of A - lam*I)
  (3) Vieta: lam1+lam2 = tr, lam1*lam2 = det
  (4) power iteration converges to the dominant eigenvalue (when one exists)
  (5) complex case: |lam|^2 = det and the displayed rotation angle
      atan2(c-b, a+d) equals the construction angle

Exit code 0 = all pass.
"""
import math, sys

def eig(A):
    a, b, c, d = A
    tr = a + d
    det = a * d - b * c
    D = tr * tr - 4.0 * det
    if D < -1e-12:
        re, im = tr / 2.0, math.sqrt(-D) / 2.0
        ang = math.atan2(c - b, a + d)
        return {"type": "complex", "re": re, "im": im, "abs": math.hypot(re, im),
                "ang": ang, "tr": tr, "det": det}
    s = math.sqrt(max(0.0, D))
    l1, l2 = (tr + s) / 2.0, (tr - s) / 2.0

    def vec(l):
        r1 = (a - l, b)
        r2 = (c, d - l)
        n1 = math.hypot(*r1); n2 = math.hypot(*r2)
        if max(n1, n2) < 1e-9:
            return None                      # A = lam*I
        p, qq = r1 if n1 >= n2 else r2       # choose the longer row of A - lam*I
        return (-qq, p)                       # eigenvector = its null space (perpendicular)

    return {"type": "real", "l1": l1, "l2": l2, "v1": vec(l1), "v2": vec(l2),
            "tr": tr, "det": det}

def apply(A, v):
    return (A[0] * v[0] + A[1] * v[1], A[2] * v[0] + A[3] * v[1])

TESTS = {
    "shear":  (1.0, 0.6, 0.0, 1.0),
    "scale":  (2.0, 0.0, 0.0, 0.5),
    "refl":   (1.0, 0.0, 0.0, -1.0),
    "rot30":  (math.cos(math.pi / 6), -math.sin(math.pi / 6),
               math.sin(math.pi / 6),  math.cos(math.pi / 6)),
    "rotsh":  (0.8 * math.cos(0.7), -0.8 * math.sin(0.7),
               0.8 * math.sin(0.7),  0.8 * math.cos(0.7)),
    "rotm90": (0.0, 1.0, -1.0, 0.0),
    "jordan": (2.0, 1.0, 0.0, 2.0),
    "id":     (1.0, 0.0, 0.0, 1.0),
    "gen":    (3.0, 1.0, 2.0, 4.0),
}

fails = []

print("=" * 66)
print("(1)+(2)+(3) per matrix: char-poly residual, Av=lam v residual, Vieta")
for name, A in TESTS.items():
    E = eig(A)
    ok = True
    msgs = []
    if E["type"] == "real":
        for l in (E["l1"], E["l2"]):
            res = abs(l * l - E["tr"] * l + E["det"])
            if res > 1e-9: ok = False; msgs.append("poly res %.3e" % res)
        worst = 0.0
        for l, v in ((E["l1"], E["v1"]), (E["l2"], E["v2"])):
            if v is None: continue          # A = lam*I: every vector works
            Av = apply(A, v)
            r = math.hypot(Av[0] - l * v[0], Av[1] - l * v[1]) / math.hypot(*v)
            worst = max(worst, r)
        if worst > 1e-9: ok = False; msgs.append("vec res %.3e" % worst)
        if abs((E["l1"] + E["l2"]) - E["tr"]) > 1e-12: ok = False; msgs.append("Vieta sum")
        if abs(E["l1"] * E["l2"] - E["det"]) > 1e-12: ok = False; msgs.append("Vieta prod")
        print("  [%s] %-6s lam=(%+.4f, %+.4f)  vec_res=%.2e  %s" %
              ("PASS" if ok else "FAIL", name, E["l1"], E["l2"], worst,
               "; ".join(msgs) or ""))
    else:
        ok &= abs(E["abs"] ** 2 - E["det"]) < 1e-9
        # residual of lam in char poly (complex substitution)
        re, im = E["re"], E["im"]
        res = math.hypot(re * re - im * im - E["tr"] * re + E["det"],
                         (2 * re - E["tr"]) * im)
        if res > 1e-9: ok = False; msgs.append("poly res %.3e" % res)
        print("  [%s] %-6s lam=%+.3f±%.3fi |lam|=%.4f  poly_res=%.2e  %s" %
              ("PASS" if ok else "FAIL", name, re, im, E["abs"], res,
               "; ".join(msgs) or ""))
    if not ok: fails.append(name)

print("=" * 66)
print("(4) power iteration -> dominant |lam| where it exists")
def power(A, iters=2000):
    v = (0.6, 0.8)
    for _ in range(iters):
        w = apply(A, v)
        n = math.hypot(*w)
        if n < 1e-15: return 0.0
        v = (w[0] / n, w[1] / n)
    Av = apply(A, v)
    return (v[0] * Av[0] + v[1] * Av[1]) / (v[0] ** 2 + v[1] ** 2)

# (jordan excluded: power iteration converges O(1/k) there - known slow, not a bug)
for name, want in [("scale", 2.0), ("gen", 5.0)]:
    A = TESTS[name]
    got = power(A)
    ok = abs(abs(got) - abs(want)) < 1e-6
    print("  [%s] %-6s power-iteration lam=%.6f (dominant |lam|=%.3f)" %
          ("PASS" if ok else "FAIL", name, got, abs(want)))
    if not ok: fails.append("power " + name)

print("=" * 66)
print("(5) complex rotation angle = construction angle")
for name, theta in [("rot30", math.pi / 6), ("rotsh", 0.7), ("rotm90", -math.pi / 2)]:
    A = TESTS[name]
    E = eig(A)
    diff = (E["ang"] - theta + math.pi) % (2 * math.pi) - math.pi
    ok = abs(diff) < 1e-9 and abs(E["abs"] - math.hypot(*[A[0], A[2]])) < 1e-9
    print("  [%s] %-6s angle=%.4f rad (want %.4f)  |lam|=%.4f" %
          ("PASS" if ok else "FAIL", name, E["ang"], theta, E["abs"]))
    if not ok: fails.append("angle " + name)

print("=" * 66)
print("(6) determinant identity det(ad-bc) vs lam1*lam2 / |lam|^2")
for name, A in TESTS.items():
    E = eig(A)
    det = A[0] * A[3] - A[1] * A[2]
    prod = E["l1"] * E["l2"] if E["type"] == "real" else E["abs"] ** 2
    ok = abs(prod - det) < 1e-9
    if not ok:
        fails.append("det " + name)
        print("  [FAIL] %-6s det=%.6f vs lam-product=%.6f" % (name, det, prod))
print("  [%s] all 9 matrices det == lam-product" % ("PASS" if not any(f.startswith("det") for f in fails) else "FAIL"))

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(fails))
sys.exit(0 if not fails else 1)
