"""
Thevenin / Norton equivalent - independent verification.

Two independent routes for R_th:
  (1) test-source method on the deactivated network (1 A injected, read V)
  (2) short-circuit current  I_sc, then R_th = V_oc / I_sc

and the series/parallel reduction must land on the same R_th.

Exit code 0 = all pass.
"""
import sys

# ---------------- tiny MNA solver ----------------
class Circuit:
    def __init__(self, elements, nodes, ref):
        self.e = elements        # list of dicts: t='R'|'V'|'I', a, b, v
        self.nodes = list(nodes)
        self.ref = ref

    def solve(self):
        nodes = [n for n in self.nodes if n != self.ref]
        idx = {n: i for i, n in enumerate(nodes)}
        vs = [e for e in self.e if e['t'] == 'V']
        n = len(nodes) + len(vs)
        A = [[0.0] * n for _ in range(n)]
        b = [0.0] * n

        for e in self.e:
            if e['t'] == 'R':
                R = e['v']
                a, bb = e['a'], e['b']
                if a != self.ref:
                    A[idx[a]][idx[a]] += 1.0 / R
                    if bb != self.ref:
                        A[idx[a]][idx[bb]] -= 1.0 / R
                if bb != self.ref:
                    A[idx[bb]][idx[bb]] += 1.0 / R
                    if a != self.ref:
                        A[idx[bb]][idx[a]] -= 1.0 / R
            elif e['t'] == 'I':
                a, bb, I = e['a'], e['b'], e['v']
                if a != self.ref:
                    b[idx[a]] -= I      # current leaves a
                if bb != self.ref:
                    b[idx[bb]] += I
        for k, e in enumerate(vs):
            row = len(nodes) + k
            a, bb = e['a'], e['b']
            if a != self.ref:
                A[row][idx[a]] += 1.0
                A[idx[a]][row] += 1.0
            if bb != self.ref:
                A[row][idx[bb]] -= 1.0
                A[idx[bb]][row] -= 1.0
            b[row] = e['v']

        # gaussian elimination with partial pivoting
        for i in range(n):
            p = max(range(i, n), key=lambda r: abs(A[r][i]))
            A[i], A[p] = A[p], A[i]
            b[i], b[p] = b[p], b[i]
            if abs(A[i][i]) < 1e-15:
                continue
            for r in range(i + 1, n):
                f = A[r][i] / A[i][i]
                if f:
                    for c in range(i, n):
                        A[r][c] -= f * A[i][c]
                    b[r] -= f * b[i]
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            if abs(A[i][i]) < 1e-15:
                x[i] = 0.0
                continue
            s = b[i] - sum(A[i][c] * x[c] for c in range(i + 1, n))
            x[i] = s / A[i][i]
        out = {self.ref: 0.0}
        for nn, i in idx.items():
            out[nn] = x[i]
        return out


def deactivate(elements):
    """V -> short (node union), I -> open (removed)."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for e in elements:
        find(e['a'])
        find(e['b'])
        if e['t'] == 'V':
            union(e['a'], e['b'])
    res = []
    for e in elements:
        if e['t'] == 'R':
            res.append({'t': 'R', 'a': find(e['a']), 'b': find(e['b']), 'v': e['v']})
    return res, find


def rth_test(elements, tp, tn):
    """Inject 1 A from tp to tn on the deactivated network; R = V(tp)-V(tn)."""
    net, find = deactivate(elements)
    nodes = set([tp, tn])
    for e in net:
        nodes.add(e['a'])
        nodes.add(e['b'])
    R = find(tn)
    els = list(net) + [{'t': 'I', 'a': R, 'b': find(tp), 'v': 1.0}]
    v = Circuit(els, nodes | {R}, R).solve()
    return v[find(tp)] - v[R]


def voc(elements, tp, tn, nodes, ref):
    v = Circuit(elements, nodes, ref).solve()
    return v[tp] - v[tn]


def isc(elements, tp, tn, nodes, ref):
    """Short the terminals with a 0 V source and read its current."""
    els = [dict(e) for e in elements] + [{'t': 'V', 'a': tp, 'b': tn, 'v': 0.0}]
    c = Circuit(els, nodes, ref)
    v = c.solve()
    # current through the 0 V source: solve again reading the source unknown
    nodes2 = [n for n in c.nodes if n != ref]
    idx = {n: i for i, n in enumerate(nodes2)}
    vs = [e for e in els if e['t'] == 'V']
    n = len(nodes2) + len(vs)
    A = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for e in els:
        if e['t'] == 'R':
            R, a, bb = e['v'], e['a'], e['b']
            if a != ref:
                A[idx[a]][idx[a]] += 1 / R
                if bb != ref:
                    A[idx[a]][idx[bb]] -= 1 / R
            if bb != ref:
                A[idx[bb]][idx[bb]] += 1 / R
                if a != ref:
                    A[idx[bb]][idx[a]] -= 1 / R
        elif e['t'] == 'I':
            a, bb, I = e['a'], e['b'], e['v']
            if a != ref:
                b[idx[a]] -= I
            if bb != ref:
                b[idx[bb]] += I
    for k, e in enumerate(vs):
        row = len(nodes2) + k
        a, bb = e['a'], e['b']
        if a != ref:
            A[row][idx[a]] += 1
            A[idx[a]][row] += 1
        if bb != ref:
            A[row][idx[bb]] -= 1
            A[idx[bb]][row] -= 1
        b[row] = e['v']
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]
        b[i], b[p] = b[p], b[i]
        if abs(A[i][i]) < 1e-15:
            continue
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            if f:
                for cc in range(i, n):
                    A[r][cc] -= f * A[i][cc]
                b[r] -= f * b[i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(A[i][i]) < 1e-15:
            x[i] = 0.0
            continue
        s = b[i] - sum(A[i][cc] * x[cc] for cc in range(i + 1, n))
        x[i] = s / A[i][i]
    # last V source is the short; current flows a->b internally
    return x[len(nodes2) + len(vs) - 1]


# ---------------- series / parallel reduction ----------------
def reduce_r(net, tp, tn, track):
    """net: list of {a,b,v}. Returns final R and appends steps to track."""
    net = [dict(e) for e in net]
    steps = []
    while True:
        # parallel
        done = False
        for i in range(len(net)):
            for j in range(i + 1, len(net)):
                same = (net[i]['a'] == net[j]['a'] and net[i]['b'] == net[j]['b']) or \
                       (net[i]['a'] == net[j]['b'] and net[i]['b'] == net[j]['a'])
                if same:
                    R = 1 / (1 / net[i]['v'] + 1 / net[j]['v'])
                    a, b = net[i]['a'], net[i]['b']
                    desc = "R%d || R%d = %.4g Ω" % (i + 1, j + 1, R)
                    net = [e for k, e in enumerate(net) if k not in (i, j)]
                    net.append({'a': a, 'b': b, 'v': R})
                    steps.append({'type': 'parallel', 'desc': desc, 'net': [dict(e) for e in net]})
                    done = True
                    break
            if done:
                break
        if done:
            continue
        # series
        deg = {}
        for e in net:
            deg[e['a']] = deg.get(e['a'], 0) + 1
            deg[e['b']] = deg.get(e['b'], 0) + 1
        for k in list(deg.keys()):
            if k in (tp, tn) or deg[k] != 2:
                continue
            inc = [i for i, e in enumerate(net) if e['a'] == k or e['b'] == k]
            i, j = inc
            other = []
            for idx in (i, j):
                other.append(net[idx]['b'] if net[idx]['a'] == k else net[idx]['a'])
            R = net[i]['v'] + net[j]['v']
            net = [e for m, e in enumerate(net) if m not in (i, j)]
            net.append({'a': other[0], 'b': other[1], 'v': R})
            steps.append({'type': 'series',
                          'desc': "R%d + R%d = %.4g Ω (node %s 제거)" % (i + 1, j + 1, R, k),
                          'net': [dict(e) for e in net]})
            done = True
            break
        if done:
            continue
        break
    return (net[0]['v'] if len(net) == 1 else None), steps


# ---------------- circuits under test ----------------
def preset1(v1=12.0, r1=4.0, r2=6.0, r3=3.0):
    nodes = ['n0', 'n1', 'n2', 't1']
    els = [{'t': 'V', 'a': 'n1', 'b': 'n0', 'v': v1},
           {'t': 'R', 'a': 'n1', 'b': 'n2', 'v': r1},
           {'t': 'R', 'a': 'n2', 'b': 'n0', 'v': r2},
           {'t': 'R', 'a': 'n2', 'b': 't1', 'v': r3}]
    return nodes, els, 't1', 'n0'


def preset2(i1=2.0, r1=5.0, r2=10.0, r3=10.0, r4=4.0):
    nodes = ['n0', 'n1', 'n2', 't1']
    els = [{'t': 'I', 'a': 'n0', 'b': 'n1', 'v': i1},
           {'t': 'R', 'a': 'n1', 'b': 'n0', 'v': r1},
           {'t': 'R', 'a': 'n1', 'b': 'n2', 'v': r2},
           {'t': 'R', 'a': 'n2', 'b': 'n0', 'v': r3},
           {'t': 'R', 'a': 'n2', 'b': 't1', 'v': r4}]
    return nodes, els, 't1', 'n0'


fails = []
for name, (nodes, els, tp, tn) in [("preset1", preset1()), ("preset2", preset2())]:
    ref = 'n0'
    V = voc(els, tp, tn, nodes, ref)
    R1 = rth_test(els, tp, tn)
    I = isc(els, tp, tn, nodes, ref)
    R2 = V / I if abs(I) > 1e-12 else float('inf')
    net, find = deactivate(els)
    R3, steps = reduce_r(net, find(tp), find(tn), None)

    print("=" * 66)
    print(name)
    print("  V_oc            = %.6f V" % V)
    print("  R_th (test src) = %.6f ohm" % R1)
    print("  I_sc            = %.6f A" % I)
    print("  R_th (V_oc/I_sc)= %.6f ohm" % R2)
    print("  R_th (series/parallel steps) = %.6f ohm" % R3 if R3 else "  reduce failed")
    print("  steps: %d" % len(steps))
    for s in steps:
        print("    -", s['desc'])

    ok = R3 is not None and abs(R1 - R2) / max(abs(R1), 1e-12) < 1e-6 and \
         abs(R1 - R3) / max(abs(R1), 1e-12) < 1e-6
    print("  [%s] three independent R_th agree" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append(name)

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ",".join(fails))
sys.exit(0 if not fails else 1)
