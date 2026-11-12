"""Independent check of a characteristic printed by STP for a gen_model.py model.

Recomputes everything from the S-box tables and the rotation definition of MixWord
(not from gen_model's matrix or CVC):
  diff: every cell transition has nonzero DDT entry; L_{i+1} = MW(U_i, V_i) xor R_i; R_{i+1} = L_i
  lin : every cell has nonzero LAT entry; masks satisfy  <L_{i+1}, MW(y1,y2)> = <v1,y1> xor <v2,y2>
        for the S-box output masks (v1, v2) on basis vectors; R_i = L_{i+1};
        L_i = R_{i+1} xor (U_i >>> b1) xor (V_i >>> b2)
and returns the total weight.

Use: from check_trail import parse, check
"""
import re

from gen_model import S1, S2, mw, rol, parity


def parse(text):
    vals = {}
    for name, v in re.findall(r"ASSERT\(\s*([A-Za-z]+\d*)\s*=\s*0x([0-9A-Fa-f]+)\s*\);", text):
        vals[name] = int(v, 16)
    return vals


def ddt(s, a, b):
    return sum(1 for x in range(16) if s[x] ^ s[x ^ a] == b)


def lat(s, a, b):
    return sum(1 for x in range(16) if parity(a & x) == parity(b & s[x])) - 8


def cells(x, n):
    return [(x >> (4 * c)) & 0xF for c in range(n // 4)]


def check(mode, vals, n, b1, b2, rounds):
    w = 0
    L = lambda i: vals.get(f"L{i}", 0)
    R = lambda i: vals.get(f"R{i}", 0)
    if L(0) == 0 and R(0) == 0:
        raise AssertionError("zero input")
    cl = n // 4
    for i in range(rounds):
        U, V = vals.get(f"U{i}", 0), vals.get(f"V{i}", 0)
        if mode == "jdiff":
            # nibble j of L feeds S1 cell j+b1/4 and S2 cell j+b2/4; joint count over the same x
            for j, a in enumerate(cells(L(i), n)):
                y1, y2 = cells(U, n)[(j + b1 // 4) % cl], cells(V, n)[(j + b2 // 4) % cl]
                c = sum(1 for x in range(16) if S1[x] ^ S1[x ^ a] == y1 and S2[x] ^ S2[x ^ a] == y2)
                assert c > 0, f"round {i}: impossible joint transition {a:x}->({y1:x},{y2:x})"
                w += {16: 0, 8: 1, 4: 2, 2: 3}[c]
            assert L(i + 1) == mw(U, V, n, b1, b2) ^ R(i), f"round {i}: Feistel/MixWord mismatch"
            assert R(i + 1) == L(i), f"round {i}: R mismatch"
            continue
        if mode == "jlin":
            g = L(i + 1)
            v1 = sum(parity(g & mw(1 << k, 0, n, b1, b2)) << k for k in range(n))
            v2 = sum(parity(g & mw(0, 1 << k, n, b1, b2)) << k for k in range(n))
            Wm = vals.get(f"W{i}", 0)
            for j, u in enumerate(cells(Wm, n)):
                y1, y2 = cells(v1, n)[(j + b1 // 4) % cl], cells(v2, n)[(j + b2 // 4) % cl]
                s = abs(sum(1 - 2 * (parity(u & x) ^ parity(y1 & S1[x]) ^ parity(y2 & S2[x])) for x in range(16)))
                assert s > 0, f"round {i}: zero joint correlation {u:x}->({y1:x},{y2:x})"
                w += {16: 0, 8: 1, 4: 2}[s]
            assert R(i) == g, f"round {i}: R mask mismatch"
            assert L(i) == R(i + 1) ^ Wm, f"round {i}: L mask mismatch"
            continue
        if mode == "diff":
            for s, b, o in ((S1, b1, U), (S2, b2, V)):
                for a, y in zip(cells(rol(L(i), b, n), n), cells(o, n)):
                    c = ddt(s, a, y)
                    assert c > 0, f"round {i}: impossible S-box transition {a:x}->{y:x}"
                    w += {16: 0, 8: 1, 4: 2, 2: 3}[c]
            assert L(i + 1) == mw(U, V, n, b1, b2) ^ R(i), f"round {i}: Feistel/MixWord mismatch"
            assert R(i + 1) == L(i), f"round {i}: R mismatch"
        else:
            g = L(i + 1)
            v1 = sum(parity(g & mw(1 << k, 0, n, b1, b2)) << k for k in range(n))
            v2 = sum(parity(g & mw(0, 1 << k, n, b1, b2)) << k for k in range(n))
            for s, u, v in ((S1, U, v1), (S2, V, v2)):
                for a, y in zip(cells(u, n), cells(v, n)):
                    l = abs(lat(s, a, y))
                    assert l > 0, f"round {i}: zero-correlation S-box mask {a:x}->{y:x}"
                    w += {8: 0, 4: 1, 2: 2}[l]
            assert R(i) == g, f"round {i}: R mask mismatch"
            ror = lambda x, r: rol(x, n - r, n)
            assert L(i) == R(i + 1) ^ ror(U, b1) ^ ror(V, b2), f"round {i}: L mask mismatch"
    return w
