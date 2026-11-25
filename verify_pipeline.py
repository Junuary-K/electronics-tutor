"""
5-stage pipeline RAW hazards / stalls / forwarding - independent verification.

Rules (must match pipeline.html, documented on the page):
  * stages IF ID EX MEM WB, one cycle each, in-order issue
  * forwarding ON : ALU producer -> consumer EX allowed at EX_j+1 or later,
                    lw  producer -> consumer EX allowed at MEM_j+1 or later
    (a value is also "forwarded" when the regfile would still be stale:
     consumer IDnat <= producer WB)
  * forwarding OFF: consumer ID must wait until WB_j + 1 (write, then read)
  * stalls freeze IF: IF_{i+1} = IF_i + 1 + stalls_i
  * $0 is never a producer; sw stores but its address base is a source;
    the *nearest* previous writer of a register is the producer

Checked against hand-derived values:
  load-use pair + 2 independent     fwd: 9 cycles / 1 stall   (textbook)
                                    off: 14 cycles / 6 stalls
  ALU chain (3 deps on one add)     fwd: 8 / 0
                                    off: 11 / 3
  independent x4                    both: 8 / 0
  distance-2 dependency             fwd: 7 / 0   off: 9 / 3

Plus structural identities for every case:
  cycles == N + 4 + stalls, last WB == cycles,
  stalls_off >= stalls_fwd, IF chaining invariant.

Exit code 0 = all pass.
"""
import re, sys

def parse(line):
    t = line.strip().replace(",", " ").split()
    op = t[0].lower()
    def R(s):
        m = re.match(r"\$(\d+)", s)
        return "$" + m.group(1) if m else None
    if op == "nop":
        return {"op": op, "dest": None, "srcs": [], "raw": line}
    if op == "lw":
        m = re.search(r"\((\$\d+)\)", t[2])
        return {"op": op, "dest": R(t[1]), "srcs": [m.group(1) if m else None], "raw": line}
    if op == "sw":
        m = re.search(r"\((\$\d+)\)", t[2])
        return {"op": op, "dest": None, "srcs": [R(t[1]), m.group(1) if m else None], "raw": line}
    return {"op": op, "dest": R(t[1]), "srcs": [R(t[2]), R(t[3])], "raw": line}

def schedule(instrs, fwd):
    cur = 1
    rows, fwd_events = [], []
    for i, ins in enumerate(instrs):
        IF, IDnat = cur, cur + 1
        req_ex, req_id = IDnat + 1, IDnat
        for r in ins["srcs"]:
            if not r or r == "$0":
                continue
            for j in range(i - 1, -1, -1):
                if instrs[j]["dest"] == r and instrs[j]["dest"] != "$0":
                    p = rows[j]
                    if fwd:
                        need = p["MEM"] + 1 if instrs[j]["op"] == "lw" else p["EX"] + 1
                        req_ex = max(req_ex, need)
                        if need >= IDnat + 1 or IDnat <= p["WB"]:
                            fwd_events.append((j, i, r))
                    else:
                        req_id = max(req_id, p["WB"] + 1)
                    break
        ID = max(IDnat, req_id)
        EX = max(ID + 1, req_ex)
        st = EX - (IDnat + 1)
        rows.append({"IF": IF, "ID": ID, "IDnat": IDnat, "EX": EX,
                     "MEM": EX + 1, "WB": EX + 2, "st": st})
        cur = IF + 1 + st
    stalls = sum(r["st"] for r in rows)
    total = rows[-1]["WB"] if rows else 0
    return rows, stalls, fwd_events, total

fails = []

CODE = {
    "loaduse": ["lw $1, 0($2)", "add $3, $1, $4", "sub $5, $3, $6", "or $7, $8, $9"],
    "alu":     ["add $1, $2, $3", "sub $4, $1, $5", "or $6, $7, $1", "and $8, $1, $9"],
    "clean":   ["add $1, $2, $3", "sub $4, $5, $6", "or $7, $8, $9", "and $10, $11, $12"],
    "far":     ["add $1, $2, $3", "or $6, $7, $8", "sub $4, $1, $5"],
}

EXPECT = {   # name -> {fwd: (cycles, stalls), off: (cycles, stalls)}
    "loaduse": {True: (9, 1),  False: (14, 6)},
    "alu":     {True: (8, 0),  False: (11, 3)},
    "clean":   {True: (8, 0),  False: (8, 0)},
    "far":     {True: (7, 0),  False: (9, 2)},
}

print("=" * 66)
print("(1) hand-derived schedule values")
for name, code in CODE.items():
    instrs = [parse(l) for l in code]
    for fwd in (True, False):
        rows, stalls, ev, total = schedule(instrs, fwd)
        want_c, want_s = EXPECT[name][fwd]
        ok = (total, stalls) == (want_c, want_s)
        print("  [%s] %-8s fwd=%-5s  cycles=%2d stalls=%d  (want %d/%d)" %
              ("PASS" if ok else "FAIL", name, fwd, total, stalls, want_c, want_s))
        if not ok:
            fails.append("%s fwd=%s" % (name, fwd))
            for r in rows:
                print("      ", r)
        # structural identities
        n = len(instrs)
        okA = total == n + 4 + stalls
        okB = rows[-1]["WB"] == total
        okC = all(rows[i + 1]["IF"] == rows[i]["IF"] + 1 + rows[i]["st"]
                  for i in range(n - 1))
        if not (okA and okB and okC):
            fails.append("identity %s fwd=%s" % (name, fwd))
            print("      identity FAIL A=%s B=%s C=%s" % (okA, okB, okC))
print("  [%s] cycles==N+4+stalls, last WB==cycles, IF chaining (all 8 runs)" %
      "PASS" if not any(f.startswith("identity") for f in fails) else "FAIL")

print("=" * 66)
print("(2) stalls_off >= stalls_fwd, per case")
for name, code in CODE.items():
    instrs = [parse(l) for l in code]
    _, s_on, _, _ = schedule(instrs, True)
    _, s_off, _, _ = schedule(instrs, False)
    ok = s_off >= s_on
    print("  [%s] %-8s  off=%d >= on=%d" % ("PASS" if ok else "FAIL", name, s_off, s_on))
    if not ok: fails.append("order " + name)

print("=" * 66)
print("(3) hazard bookkeeping details")
# nearest producer wins: second write of $1 shadows the first
instrs = [parse(l) for l in ["add $1, $2, $3", "add $1, $4, $5", "or $6, $1, $7"]]
rows, stalls, ev, total = schedule(instrs, True)
# consumer at i=2 must wait for i=1 (EX@4+1=5), not i=0 (EX@3+1=4)
ok = rows[2]["EX"] == 5 and stalls == 0
print("  [%s] nearest-producer shadowing: EX=%d (want 5)" % ("PASS" if ok else "FAIL", rows[2]["EX"]))
if not ok: fails.append("shadowing")

# $0 never stalls a consumer
instrs = [parse(l) for l in ["add $0, $1, $2", "or $3, $0, $4"]]
rows, stalls, ev, total = schedule(instrs, True)
ok = stalls == 0 and total == 6
print("  [%s] $0 as source causes no stall: stalls=%d total=%d" % ("PASS" if ok else "FAIL", stalls, total))
if not ok: fails.append("$0")

# sw address base is a source (dependency), sw itself is not a producer
instrs = [parse(l) for l in ["lw $1, 0($2)", "sw $1, 4($1)", "or $5, $6, $7"]]
rows, stalls, ev, total = schedule(instrs, True)
# sw depends on $1 from lw at distance 1 -> MEM->EX ... sw uses $1 as address,
# needed in EX: 1 stall (load-use into sw counts)
ok = stalls == 1 and total == 3 + 4 + 1
print("  [%s] lw->sw base dep: stalls=%d total=%d (want 1/%d)" %
      ("PASS" if ok else "FAIL", stalls, total, 3 + 4 + 1))
if not ok: fails.append("sw base")

instrs = [parse(l) for l in ["sw $1, 4($2)", "or $5, $1, $6"]]
rows, stalls, ev, total = schedule(instrs, True)
ok = stalls == 0
print("  [%s] sw is not a producer: stalls=%d" % ("PASS" if ok else "FAIL", stalls))
if not ok: fails.append("sw producer")

# nop contributes cycles but never hazards
instrs = [parse(l) for l in ["add $1, $2, $3", "nop", "or $4, $1, $5"]]
rows, stalls, ev, total = schedule(instrs, True)
ok = stalls == 0 and total == 7
print("  [%s] nop fills a slot: stalls=%d total=%d (want 0/7)" % ("PASS" if ok else "FAIL", stalls, total))
if not ok: fails.append("nop")

print("=" * 66)
print("RESULT:", "ALL PASS" if not fails else "FAIL: " + ", ".join(fails))
sys.exit(0 if not fails else 1)
