"""
Differential equations solving methods - independent verification.

Mirrors ode.html's two representative methods:
  (1) 1st-order linear  y' + a y = b e^{kt}, integrating factor
      (non-resonant k != -a and resonant k == -a)
  (2) 2nd-order constant-coefficient  a y'' + b y' + c y = f(t),
      characteristic equation (over/crit/under-damped) + undetermined
      coefficients for f = const and f = F cos(wt) (incl. resonance)

Independent route: the closed-form algebra here is checked against
  * a fresh RK4 integrator (h = T/2000) on the trajectory,
  * direct substitution residuals into the original ODE,
  * textbook exact coefficient values for the classic damping cases.

Exit code 0 = all pass.
"""
import math, sys

fails = []

# ---------------- fresh RK4 ----------------
def rk4(f, t0, s0, T, N):
    h = T / N
    t, s = t0, list(s0)
    out = []
    for i in range(N + 1):
        out.append((t, s[0]))
        k1 = f(t, s)
        k2 = f(t + h / 2, [s[j] + h / 2 * k1[j] for j in range(2)])
        k3 = f(t + h / 2, [s[j] + h / 2 * k2[j] for j in range(2)])
        k4 = f(t + h, [s[j] + h * k3[j] for j in range(2)])
        s = [s[j] + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(2)]
        t += h
    return out

# ---------------- 1st order ----------------
def solve1(a, b, k, y0, t):
    if k != -a:
        A = b / (k + a)
        C = y0 - A
        return C * math.exp(-a * t) + A * math.exp(k * t)
    C = y0
    return (b * t + C) * math.exp(-a * t)

def solve1_d(a, b, k, y0, t):
    if k != -a:
        A = b / (k + a); C = y0 - A
        return -a * C * math.exp(-a * t) + A * k * math.exp(k * t)
    C = y0
    return b * math.exp(-a * t) - a * (b * t + C) * math.exp(-a * t)

print("=" * 66)
print("(1) 1st-order: closed form vs RK4 + plug residual")
CASES1 = [(2.0, 3.0, 1.0, 0.0), (1.0, 2.0, -1.0, 0.5), (0.5, -3.0, 2.0, 1.0), (3.0, 1.0, -3.0, -2.0)]
for a, b, k, y0 in CASES1:
    tag = "resonant " if k == -a else "normal    "
    def f(t, s):
        return [s[1], b * k * math.exp(k * t) - a * s[1]]   # v' = d/dt(b e^{kt} - a y)
    traj = rk4(f, 0, [y0, b - a * y0], 8.0, 2000)
    emax = 0.0
    for i in range(0, 2001, 25):
        t, v = traj[i]
        emax = max(emax, abs(v - solve1(a, b, k, y0, t)))
    emax = emax / (1.0 + abs(solve1(a, b, k, y0, 8.0)))     # scale by trajectory end
    plug = 0.0
    for t in (0.3, 1.2, 2.7, 5.5):
        plug = max(plug, abs(solve1_d(a, b, k, y0, t) + a * solve1(a, b, k, y0, t) - b * math.exp(k * t)))
    ok = emax < 1e-8 and plug < 1e-12 and abs(solve1(a, b, k, y0, 0) - y0) < 1e-12
    print("  [%s] a=%.1f b=%.1f k=%+.1f y0=%+.1f  RK4 rel err=%.2e  plug=%.2e" %
          ("PASS" if ok else "FAIL", a, b, k, y0, emax, plug))
    if not ok: fails.append("1st a=%s" % a)

# ---------------- 2nd order ----------------
def solve2(a, b, c, ftype, F, w, y0, y1):
    """returns (y(t), yp(t), yp_d(t), yp_dd(t), case, extras) via closed algebra."""
    if ftype == "const":
        yp  = lambda t: F / c
        ypd = lambda t: 0.0
        ypdd = lambda t: 0.0
    elif ftype == "cos":
        den = (c - a * w * w) ** 2 + (b * w) ** 2
        if abs(den) < 1e-15:
            A = F / (2 * a * w)
            yp   = lambda t: A * t * math.sin(w * t)
            ypd  = lambda t: A * (math.sin(w * t) + w * t * math.cos(w * t))
            ypdd = lambda t: A * (2 * w * math.cos(w * t) - w * w * t * math.sin(w * t))
        else:
            M = F * (c - a * w * w) / den
            N = F * b * w / den
            yp   = lambda t: M * math.cos(w * t) + N * math.sin(w * t)
            ypd  = lambda t: -M * w * math.sin(w * t) + N * w * math.cos(w * t)
            ypdd = lambda t: -w * w * (M * math.cos(w * t) + N * math.sin(w * t))
    else:
        yp = ypd = ypdd = lambda t: 0.0

    D = b * b - 4 * a * c
    e0, e1 = y0 - yp(0), y1 - ypd(0)
    if D > 1e-12:
        r1 = (-b + math.sqrt(D)) / (2 * a); r2 = (-b - math.sqrt(D)) / (2 * a)
        C2 = (e1 - r1 * e0) / (r2 - r1); C1 = e0 - C2
        y = lambda t: C1 * math.exp(r1 * t) + C2 * math.exp(r2 * t) + yp(t)
        return y, yp, ypd, ypdd, "over", dict(r1=r1, r2=r2, C1=C1, C2=C2, D=D)
    if D > -1e-12:
        r = -b / (2 * a)
        C1 = e0; C2 = e1 - r * C1
        y = lambda t: (C1 + C2 * t) * math.exp(r * t) + yp(t)
        return y, yp, ypd, ypdd, "crit", dict(r=r, C1=C1, C2=C2, D=D)
    al = -b / (2 * a); be = math.sqrt(-D) / (2 * a)
    C1 = e0; C2 = (e1 - al * C1) / be
    y = lambda t: math.exp(al * t) * (C1 * math.cos(be * t) + C2 * math.sin(be * t)) + yp(t)
    return y, yp, ypd, ypdd, "under", dict(al=al, be=be, C1=C1, C2=C2, D=D)

print("=" * 66)
print("(2) 2nd-order: closed form vs RK4 + particular residual + char roots")
CASES2 = [
    ("over  ", 1.0, 3.0, 2.0, "0",    0.0, 0.0, 1.0, 0.0),
    ("crit  ", 1.0, 2.0, 1.0, "0",    0.0, 0.0, 0.0, 1.0),
    ("under ", 1.0, 2.0, 5.0, "0",    0.0, 0.0, 1.0, 0.0),
    ("const ", 1.0, 3.0, 2.0, "const", 4.0, 0.0, 1.0, 0.0),
    ("cos   ", 1.0, 2.0, 2.0, "cos",  4.0, 1.0, 1.0, 0.0),
    ("rescos", 1.0, 0.0, 4.0, "cos",  4.0, 2.0, 0.0, 1.0),
    ("over2 ", 2.0, 5.0, 2.0, "cos",  3.0, 1.5, -1.0, 2.0),
]
for name, a, b, c, ftype, F, w, y0, y1 in CASES2:
    y, yp, ypd, ypdd, case, ex = solve2(a, b, c, ftype, F, w, y0, y1)
    def f(t, s):
        ft = F if ftype == "const" else (F * math.cos(w * t) if ftype == "cos" else 0.0)
        return [s[1], (ft - b * s[1] - c * s[0]) / a]
    traj = rk4(f, 0, [y0, y1], 8.0, 2000)
    emax = 0.0
    for i in range(0, 2001, 25):
        t, v = traj[i]
        emax = max(emax, abs(v - y(t)))
    pres = 0.0
    for t in (0.4, 1.3, 3.1, 6.2):
        ft = F if ftype == "const" else (F * math.cos(w * t) if ftype == "cos" else 0.0)
        pres = max(pres, abs(a * ypdd(t) + b * ypd(t) + c * yp(t) - ft))
    cres = 0.0
    if case == "over":
        cres = max(abs(a * ex["r1"] ** 2 + b * ex["r1"] + c), abs(a * ex["r2"] ** 2 + b * ex["r2"] + c))
    elif case == "crit":
        cres = abs(a * ex["r"] ** 2 + b * ex["r"] + c)
    else:
        re = a * (ex["al"] ** 2 - ex["be"] ** 2) + b * ex["al"] + c
        im = b * ex["be"] + 2 * a * ex["al"] * ex["be"]
        cres = max(abs(re), abs(im))
    ic = abs(y(0) - y0)
    ok = emax < 1e-6 and pres < 1e-9 and cres < 1e-9 and ic < 1e-12
    print("  [%s] %-6s RK4 err=%.2e  yp residual=%.2e  char=%.2e  ic=%.1e" %
          ("PASS" if ok else "FAIL", name, emax, pres, cres, ic))
    if not ok: fails.append("2nd %s" % name)

print("=" * 66)
print("(3) textbook exact values (damping cases + resonance)")
def coeffs(a, b, c, ftype, F, w, y0, y1):
    _, _, _, _, case, ex = solve2(a, b, c, ftype, F, w, y0, y1)
    return case, ex["C1"], ex["C2"]

tests = [
    ((1, 3, 2, "0", 0, 0, 1, 0), "over", 2.0, -1.0),          # y = 2e^-t - e^-2t
    ((1, 2, 1, "0", 0, 0, 0, 1), "crit", 0.0, 1.0),            # y = t e^-t
    ((1, 2, 5, "0", 0, 0, 1, 0), "under", 1.0, 0.5),          # e^-t(cos2t + .5 sin2t)
    ((1, 3, 2, "const", 4, 0, 1, 0), "over", -2.0, 1.0),      # -2e^-t + e^-2t + 2
]
okall = True
for args, want_case, want_c1, want_c2 in tests:
    case, c1, c2 = coeffs(*args)
    ok = case == want_case and abs(c1 - want_c1) < 1e-12 and abs(c2 - want_c2) < 1e-12
    okall &= ok
    print("  [%s] %-6s C1=%+.4f C2=%+.4f (want %s %+.1f %+.1f)" %
          ("PASS" if ok else "FAIL", case, c1, c2, want_case, want_c1, want_c2))
    if not ok: fails.append("textbook %s" % want_case)

# 1st-order resonance exact: y' + 2y = 5e^{-2t}, y(0)=3 -> (5t+3)e^{-2t}
v = solve1(2.0, 5.0, -2.0, 3.0, 1.7)
want = (5 * 1.7 + 3) * math.exp(-3.4)
ok = abs(v - want) < 1e-12
okall &= ok
print("  [%s] 1st resonance y(1.7) = %.6f (want %.6f)" % ("PASS" if ok else "FAIL", v, want))
if not ok: fails.append("1st resonance exact")

print("=" * 66)
print("(4) non-resonant cos particular: M,N solve the algebraic system")
a, b, c, F, w = 1.0, 2.0, 2.0, 4.0, 1.0
_, yp, ypd, ypdd, _, _ = solve2(a, b, c, "cos", F, w, 0, 0)
# M cos + N sin coefficients from the closed form at t=0, t=pi/2w
M = yp(0.0); N = yp(math.pi / (2 * w))
lhs1 = (c - a * w * w) * M + b * w * N      # should equal F
lhs2 = -b * w * M + (c - a * w * w) * N     # should equal 0
ok = abs(lhs1 - F) < 1e-12 and abs(lhs2) < 1e-12
print("  [%s] (c-aw^2)M + bwN = %.6f (want %.1f), -bwM + (c-aw^2)N = %.2e (want 0)" %
      ("PASS" if ok else "FAIL", lhs1, F, abs(lhs2)))
if not ok: fails.append("cos MN system")

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(sorted(set(fails))))
sys.exit(0 if not fails else 1)
