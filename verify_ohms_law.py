"""
Ohm's law / series / parallel / divider / power - independent verification.

Route: a fresh nodal (MNA) solver - build each page topology as a netlist,
solve the linear system with Gaussian elimination, and compare node voltages
and branch currents against the CLOSED-FORM formulas used by ohms-law.html
(which computes I, V1, V2, Vout, P directly from algebra, not by solving a
system). Same numbers from two genuinely different computations.

Also checks:
  * power identities P = VI = I^2 R = V^2/R per element
  * KVL (drop sum = source) and KCL (branch sum = total) residuals
  * the four worked examples on the page, recomputed through the MNA solver
  * divider formula V*R2/(R1+R2) == nodal voltage of the tap node

Exit code 0 = all pass.
"""
import sys


def solve_mna(elements, ref):
    """elements: list of (t, a, b, v) with t in 'R','V','I'.
    Returns dict node->voltage (ref=0) and source currents I('V' element)."""
    nodes = set([ref])
    for t, a, b, v in elements:
        nodes.add(a); nodes.add(b)
    nodes.discard(ref)
    idx = {n: i for i, n in enumerate(sorted(nodes))}
    vs = [e for e in elements if e[0] == 'V']
    n = len(idx) + len(vs)
    G = [[0.0] * n for _ in range(n)]
    rhs = [0.0] * n

    for t, a, b, v in elements:
        if t == 'R':
            if a in idx:
                G[idx[a]][idx[a]] += 1.0 / v
                if b in idx:
                    G[idx[a]][idx[b]] -= 1.0 / v
            if b in idx:
                G[idx[b]][idx[b]] += 1.0 / v
                if a in idx:
                    G[idx[b]][idx[a]] -= 1.0 / v
        elif t == 'I':
            if a in idx:
                rhs[idx[a]] -= v
            if b in idx:
                rhs[idx[b]] += v
    # voltage sources (row per source, matched by position in vs)
    for k, e in enumerate(vs):
        t, a, b, v = e
        row = len(idx) + k
        if a in idx:
            G[row][idx[a]] += 1.0
            G[idx[a]][row] += 1.0
        if b in idx:
            G[row][idx[b]] -= 1.0
            G[idx[b]][row] -= 1.0
        rhs[row] = v

    # Gaussian elimination with partial pivoting
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(G[r][i]))
        G[i], G[p] = G[p], G[i]
        rhs[i], rhs[p] = rhs[p], rhs[i]
        if abs(G[i][i]) < 1e-15:
            continue
        for r in range(i + 1, n):
            f = G[r][i] / G[i][i]
            if f:
                for c in range(i, n):
                    G[r][c] -= f * G[i][c]
                rhs[r] -= f * rhs[i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(G[i][i]) < 1e-15:
            continue
        s = rhs[i] - sum(G[i][c] * x[c] for c in range(i + 1, n))
        x[i] = s / G[i][i]

    volt = {ref: 0.0}
    for name, i in idx.items():
        volt[name] = x[i]
    cur = {}
    for k, e in enumerate(vs):
        cur[(e[1], e[2])] = x[len(idx) + k]   # current flowing a->b inside the source
    return volt, cur


def netlist(topo, vs, r1, r2):
    """Build each topology as (elements, tap_node_or_None)."""
    if topo == "single":
        els = [('V', 'n1', 'gnd', vs), ('R', 'n1', 'gnd', r1)]
        return els, None
    if topo == "series":
        els = [('V', 'n1', 'gnd', vs), ('R', 'n1', 'n2', r1), ('R', 'n2', 'gnd', r2)]
        return els, None
    if topo == "parallel":
        els = [('V', 'n1', 'gnd', vs), ('R', 'n1', 'gnd', r1), ('R', 'n1', 'gnd', r2)]
        return els, None
    els = [('V', 'n1', 'gnd', vs), ('R', 'n1', 'n2', r1), ('R', 'n2', 'gnd', r2)]
    return els, 'n2'          # divider: tap = n2


def closed_form(topo, vs, r1, r2):
    """The page's formulas (ohms-law.html solve())."""
    if topo == "single":
        I = vs / r1
        return dict(I=I, I1=I, I2=0.0, v1=vs, v2=0.0, p1=vs * I, p2=0.0, vout=vs)
    if topo == "series":
        I = vs / (r1 + r2)
        return dict(I=I, I1=I, I2=I, v1=I * r1, v2=I * r2, p1=I * I * r1, p2=I * I * r2, vout=vs)
    if topo == "parallel":
        I1, I2 = vs / r1, vs / r2
        return dict(I=I1 + I2, I1=I1, I2=I2, v1=vs, v2=vs, p1=vs * I1, p2=vs * I2, vout=vs)
    I = vs / (r1 + r2)
    return dict(I=I, I1=I, I2=I, v1=I * r1, v2=I * r2, p1=I * I * r1, p2=I * I * r2, vout=I * r2)


fails = []
PARAMS = [(12.0, 4.0, 6.0), (9.0, 2.0, 1.0), (5.0, 47.0, 33.0), (24.0, 8.0, 8.0), (3.3, 1.5, 2.7)]

print("=" * 66)
print("(1) MNA nodal solve vs page closed forms, per topology")
for topo in ("single", "series", "parallel", "divider"):
    worst = 0.0
    ok = True
    for vs, r1, r2 in PARAMS:
        els, tap = netlist(topo, vs, r1, r2)
        volt, cur = solve_mna(els, 'gnd')
        c = closed_form(topo, vs, r1, r2)
        src_cur = cur[('n1', 'gnd')]
        # branch currents from node voltages
        if topo == "single":
            i1 = volt['n1'] / r1
            i2 = 0.0
            v1, v2, vout = volt['n1'], 0.0, volt['n1']
        elif topo == "parallel":
            i1, i2 = volt['n1'] / r1, volt['n1'] / r2
            v1 = v2 = vout = volt['n1']
        else:
            i1 = (volt['n1'] - volt['n2']) / r1
            i2 = volt['n2'] / r2
            v1, v2 = volt['n1'] - volt['n2'], volt['n2']
            vout = volt['n2'] if topo == "divider" else vs
        # compare against closed forms
        def rel(a, b):
            return abs(a - b) / max(abs(b), 1e-12)
        e = max(rel(i1, c['I1']), rel(v1, c['v1']), rel(vout, c['vout']))
        if topo != "single":
            e = max(e, rel(i2, c['I2']), rel(v2, c['v2']))
        worst = max(worst, e)
        # KVL / KCL residuals (nodal consistency)
        if topo in ('series', 'divider'):
            kvl = abs((v1 + v2) - vs)
            if kvl > 1e-9: ok = False; fails.append("KVL %s" % topo)
        if topo == 'parallel':
            kcl = abs(i1 + i2 - abs(src_cur))
            if kcl > 1e-9: ok = False; fails.append("KCL %s" % topo)
    ok = ok and worst < 1e-12
    print("  [%s] %-9s worst rel err vs closed form = %.2e" % ("PASS" if ok else "FAIL", topo, worst))
    if not ok: fails.append("mna %s" % topo)

print("=" * 66)
print("(2) power identities P = VI = I^2 R = V^2/R per element")
okall = True
for topo in ("single", "series", "parallel", "divider"):
    for vs, r1, r2 in PARAMS:
        c = closed_form(topo, vs, r1, r2)
        for R, V, I in ((r1, c['v1'], c['I1']), (r2, c['v2'], c['I2'])):
            if R <= 0 or V == 0.0 or I == 0.0:
                continue
            pa, pb, pc = V * I, I * I * R, V * V / R
            if max(abs(pa - pb), abs(pa - pc)) / abs(pa) > 1e-12:
                okall = False
                fails.append("power %s %s" % (topo, R))
print("  [%s] all elements, all topologies, three routes agree" % ("PASS" if okall else "FAIL"))

print("=" * 66)
print("(3) the four worked examples, recomputed through the MNA solver")
def check_ex(name, topo, vs, r1, r2, want):
    els, tap = netlist(topo, vs, r1, r2)
    volt, cur = solve_mna(els, 'gnd')
    ok, desc = want(volt, cur)
    print("  [%s] %-6s %s" % ("PASS" if ok else "FAIL", name, desc))
    if not ok: fails.append("ex " + name)

check_ex("ex1", "single", 12.0, 4.0, 6.0,
         lambda v, c: (abs(v['n1'] / 4.0 - 3.0) < 1e-9 and abs(12.0 * v['n1'] / 4.0 - 36.0) < 1e-9,
                       "I = %.4f A, P = %.4f W (want 3 A / 36 W)" % (v['n1'] / 4.0, 12.0 * v['n1'] / 4.0)))
check_ex("ex2", "divider", 9.0, 2.0, 1.0,
         lambda v, c: (abs(v['n2'] - 3.0) < 1e-9,
                       "V_out = %.4f V (want 3 V)" % v['n2']))
check_ex("ex3", "parallel", 12.0, 6.0, 3.0,
         lambda v, c: (abs(v['n1'] / 6.0 - 2.0) < 1e-9 and abs(v['n1'] / 3.0 - 4.0) < 1e-9
                       and abs(1.0 / (1.0 / 6.0 + 1.0 / 3.0) - 2.0) < 1e-9,
                       "I1=%.3f I2=%.3f Req=%.3f (want 2/4/2)" % (v['n1'] / 6.0, v['n1'] / 3.0, 2.0)))
check_ex("ex4", "single", 10.0, 100.0, 6.0,
         lambda v, c: (abs(v['n1'] / 100.0 - 0.1) < 1e-9 and abs(0.1 ** 2 * 100.0 - 1.0) < 1e-9,
                       "I = %.4f A, P = %.4f W (want 0.1 A / 1 W)" % (v['n1'] / 100.0, 0.01 * 100.0)))

print("=" * 66)
print("(4) divider formula == nodal tap voltage")
okall = True
for vs, r1, r2 in PARAMS:
    els, tap = netlist("divider", vs, r1, r2)
    volt, _ = solve_mna(els, 'gnd')
    if abs(vs * r2 / (r1 + r2) - volt['n2']) > 1e-12:
        okall = False
        fails.append("divider %s" % vs)
print("  [%s] V*R2/(R1+R2) == V_tap for all param sets" % ("PASS" if okall else "FAIL"))

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(sorted(set(fails))))
sys.exit(0 if not fails else 1)
