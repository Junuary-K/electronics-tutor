"""
Convolution y(t) = x(t)*h(t) - independent verification.

Convention (must match convolution.html):
  * signals are piecewise-defined with inclusive edges,
  * every signal argument is quantized to 10 decimals (q) so both convolution
    orders make identical boundary decisions,
  * the discrete sum y(tau) = sum_k x(k) h(tau-k) dt therefore carries an
    O(dt) edge-cell term: tolerance for closed-form comparison is dt itself,
    and refinement must converge at first order (error halves with dt).

Checks:
  (1) discrete sum vs closed forms (rect*rect trapezoids, rect*exp piecewise),
      including first-order scaling between dt=0.04 and dt=0.02
  (2) commutativity  x*h == h*x  (exact, by construction of q)
  (3) delta identity x*delta == x (half-open unit-area pulse)
  (4) refinement: interior error halves with dt and stays <= dt
  (5) hand-checkable spot values within dt

Exit code 0 = all pass.
"""
import math, sys

T0, T1 = -3.0, 9.0
def q(v): return round(v, 10)

def grid(dt):
    n = int(round((T1 - T0) / dt))
    return [q(T0 + i * dt) for i in range(n + 1)]

SIG = {
    "rect1":  lambda t: 1.0 if 0.0 <= t <= 1.0 else 0.0,
    "rect2":  lambda t: 1.0 if 0.0 <= t <= 2.0 else 0.0,
    "rect3":  lambda t: 1.0 if 0.0 <= t <= 3.0 else 0.0,
    "rect05": lambda t: 1.0 if 0.0 <= t <= 0.5 else 0.0,
    "exp":    lambda t: math.exp(-t) if t >= 0.0 else 0.0,
    "cos":    lambda t: math.cos(3.0 * math.pi * t) if 0.0 <= t <= 2.0 else 0.0,
}

def conv(xf, hf, dt):
    g = grid(dt)
    out = {}
    for tau in g:
        s = 0.0
        for k in g:
            hv = hf(q(tau - k))
            if hv != 0.0:
                xv = xf(q(k))
                if xv != 0.0:
                    s += xv * hv * dt
        out[tau] = s
    return out

def overlap_len(a, b, w, tau):
    return max(0.0, min(b, tau) - max(a, tau - w))

def analytic(xk, hk, tau):
    if xk == "rect2" and hk == "rect3": return overlap_len(0, 2, 3, tau)
    if xk == "rect1" and hk == "rect3": return overlap_len(0, 1, 3, tau)
    if xk == "rect2" and hk == "rect05": return overlap_len(0, 2, 0.5, tau)
    if xk == "rect1" and hk == "rect05": return overlap_len(0, 1, 0.5, tau)
    if xk == "rect2" and hk == "exp":
        if tau <= 0: return 0.0
        if tau <= 2: return 1.0 - math.exp(-tau)
        return math.exp(-(tau - 2)) * (1.0 - math.exp(-2.0))
    if xk == "rect1" and hk == "exp":
        if tau <= 0: return 0.0
        if tau <= 1: return 1.0 - math.exp(-tau)
        return math.exp(-(tau - 1)) * (1.0 - math.exp(-1.0))
    return None

fails = []

# ---------- (1) closed forms + first-order scaling ----------
print("=" * 66)
print("(1) discrete sum vs closed form  (edge-cell error must be <= dt)")
COMBOS = [("rect2","rect3"), ("rect1","rect3"), ("rect2","rect05"),
          ("rect1","rect05"), ("rect2","exp"), ("rect1","exp")]
err_by_dt = {}
for dt in (0.04, 0.02):
    err_by_dt[dt] = {}
    for xk, hk in COMBOS:
        y = conv(SIG[xk], SIG[hk], dt)
        emax = 0.0
        for tau, v in y.items():
            a = analytic(xk, hk, tau)
            if a is None: continue
            emax = max(emax, abs(v - a))
        err_by_dt[dt][(xk, hk)] = emax
for xk, hk in COMBOS:
    e02 = err_by_dt[0.02][(xk, hk)]
    e04 = err_by_dt[0.04][(xk, hk)]
    ok = e02 <= 0.02 + 1e-9 and e02 <= e04 / 1.9 + 1e-12
    print("  [%s] %-6s*%-6s  err(dt=.02)=%.3e  err(dt=.04)=%.3e  halving=%s" %
          ("PASS" if ok else "FAIL", xk, hk, e02, e04, "yes" if e02 <= e04 / 1.9 else "NO"))
    if not ok: fails.append("analytic %s*%s" % (xk, hk))

# ---------- (2) commutativity ----------
print("=" * 66)
print("(2) commutativity x*h == h*x (quantized -> exact)")
for xk, hk in [("rect2","exp"), ("cos","rect05"), ("rect1","rect3")]:
    a = conv(SIG[xk], SIG[hk], 0.05)
    b = conv(SIG[hk], SIG[xk], 0.05)
    em = max(abs(va - b[tau]) for tau, va in a.items())
    ok = em < 1e-9
    print("  [%s] %-6s*%-6s  max|err|=%.3e" % ("PASS" if ok else "FAIL", xk, hk, em))
    if not ok: fails.append("commutative %s*%s" % (xk, hk))

# ---------- (3) delta identity ----------
print("=" * 66)
print("(3) x * delta == x   (half-open unit-area pulse, grid-aligned)")
dt = 0.02
delta = lambda t: (1.0 / dt) if 0.0 <= t < dt else 0.0
for xk in ["rect2", "cos", "exp"]:
    y = conv(SIG[xk], delta, dt)
    em = max(abs(v - SIG[xk](tau)) for tau, v in y.items())
    ok = em < 1e-9
    print("  [%s] x=%-6s  max|err|=%.3e" % ("PASS" if ok else "FAIL", xk, em))
    if not ok: fails.append("delta %s" % xk)

# ---------- (4) refinement on smooth interior ----------
print("=" * 66)
print("(4) refinement: interior error halves with dt and stays <= dt")
def interior_err(dt):
    y = conv(SIG["rect2"], SIG["exp"], dt)
    e = 0.0
    for tau, v in y.items():
        if 0.35 < tau < 1.6:
            e = max(e, abs(v - analytic("rect2", "exp", tau)))
    return e
e1 = interior_err(0.04); e2 = interior_err(0.02)
ok = e2 <= e1 / 1.9 and e2 <= 0.02 + 1e-9
print("  err(dt=.04)=%.3e  err(dt=.02)=%.3e  ratio=%.2f" % (e1, e2, e1 / e2 if e2 else float("inf")))
print("  [%s] first-order refinement" % ("PASS" if ok else "FAIL"))
if not ok: fails.append("refinement")

# ---------- (5) spot values ----------
print("=" * 66)
print("(5) spot values: rect2*rect3 within dt (tau=1 -> 1, 2.5 -> 2, 4 -> 1)")
y = conv(SIG["rect2"], SIG["rect3"], 0.02)
okall = True
for tau, want in [(1.0, 1.0), (2.5, 2.0), (4.0, 1.0)]:
    got = y[q(tau)]
    ok = abs(got - want) <= 0.02 + 1e-9
    okall &= ok
    print("  [%s] y(%.1f) = %.4f (want %.1f +-dt)" % ("PASS" if ok else "FAIL", tau, got, want))
if not okall: fails.append("spot")

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(fails))
sys.exit(0 if not fails else 1)
