"""STP (CVC) models of differential and linear characteristics for the LLBC-type cipher of the
official implementation, on either the 32-bit part cipher or the full 128-bit cipher.

Cipher (mini_llbc.py / split_test.py): halves of n bits, cells of 4 bits,
    y1 = S1(L <<< b1),  y2 = S2(L <<< b2),  L' = MW(y1, y2) xor R xor key,  R' = L xor key
    full official code: n = 64, (b1, b2) = (32, 16);   one part: n = 16, (b1, b2) = (8, 4).
MW is used through its GF(2) matrix, computed from the rotation definition below.

Weights: differential  w = -log2(DDT[a][b] / 16)   (DDT entries 16, 4, 2 -> 0, 2, 3)
         linear        w = -log2(|LAT[a][b]| / 8)  (correlation weight; |LAT| 8, 4, 2 -> 0, 1, 2)
Characteristic weight = sum over all S-box cells of all rounds.  Keys/constants do not matter.

Run: python gen_model.py {diff|lin} n b1 b2 rounds W  > model.cvc
     The model is satisfiable iff a characteristic of weight <= W with nonzero input exists.
"""
import sys

S1 = [0x4, 0x7, 0x5, 0x1, 0xc, 0xb, 0xd, 0x8, 0xe, 0xf, 0x6, 0x3, 0xa, 0x9, 0x2, 0x0]
S2 = [0xe, 0xc, 0x8, 0x4, 0xb, 0xd, 0x1, 0x5, 0xf, 0x3, 0xa, 0x2, 0x9, 0x7, 0x0, 0x6]
INVALID = 0xF


def rol(x, r, n):
    r %= n
    return ((x << r) | (x >> (n - r))) & ((1 << n) - 1) if r else x


def mw(y1, y2, n, b1, b2):
    t = y1 ^ y2
    m0 = rol(t, b1, n) ^ y1
    m1 = rol(m0, b2, n) ^ t
    m2 = rol(m1, b1, n) ^ m0
    m3 = rol(m2, b2, n) ^ m1
    return m2 ^ m3


def mw_matrix(n, b1, b2):
    """cols[i] = image of input bit i (i < n: bit i of y1, else bit i-n of y2)."""
    return [mw(1 << i, 0, n, b1, b2) if i < n else mw(0, 1 << (i - n), n, b1, b2) for i in range(2 * n)]


def ddt_weights(s):
    w = {}
    for a in range(16):
        cnt = [0] * 16
        for x in range(16):
            cnt[s[x] ^ s[x ^ a]] += 1
        for b in range(16):
            w[(a, b)] = {16: 0, 8: 1, 4: 2, 2: 3, 0: INVALID}[cnt[b]]
    return w


def parity(x):
    return bin(x).count("1") & 1


def lat_weights(s):
    w = {}
    for a in range(16):
        for b in range(16):
            l = sum(1 for x in range(16) if parity(a & x) == parity(b & s[x])) - 8
            w[(a, b)] = {8: 0, 4: 1, 2: 2, 0: INVALID}[abs(l)]
    return w


def bit(var, i):
    return f"{var}[{i}:{i}]"


def xor_bits(terms):
    if not terms:
        return "0bin0"
    e = terms[0]
    for t in terms[1:]:
        e = f"BVXOR({e},{t})"
    return e


def cat_bits(exprs_msb_first):
    return "@".join(exprs_msb_first)


def joint_ddt_weights():
    """x -> (S1(x), S2(x)) as one 4-to-8-bit map (both S-boxes read the same unkeyed nibble)."""
    w = {}
    for a in range(16):
        for b1 in range(16):
            for b2 in range(16):
                c = sum(1 for x in range(16) if S1[x] ^ S1[x ^ a] == b1 and S2[x] ^ S2[x ^ a] == b2)
                w[(a, b1, b2)] = {16: 0, 8: 1, 4: 2, 2: 3, 0: INVALID}[c]
    return w


def joint_lat_weights():
    w = {}
    for u in range(16):
        for v1 in range(16):
            for v2 in range(16):
                s = sum(1 - 2 * (parity(u & x) ^ parity(v1 & S1[x]) ^ parity(v2 & S2[x])) for x in range(16))
                w[(u, v1, v2)] = {16: 0, 8: 1, 4: 2, 0: INVALID}[abs(s)]
    return w


def generate(mode, n, b1, b2, rounds, W, fixed=None):
    cells = n // 4
    cols = mw_matrix(n, b1, b2)
    out = []
    if mode in ("jdiff", "jlin"):
        # S1 cell (j + b1/4) and S2 cell (j + b2/4) both read nibble j of L (needs 4 | b1, b2)
        assert b1 % 4 == 0 and b2 % 4 == 0
        out.append("J : ARRAY BITVECTOR(12) OF BITVECTOR(4);")
        for (a, x, y), w in (joint_ddt_weights() if mode == "jdiff" else joint_lat_weights()).items():
            out.append(f"ASSERT(J[0hex{a:x}{x:x}{y:x}] = 0hex{w:x});")
    else:
        tabs = (ddt_weights if mode == "diff" else lat_weights)
        for name, s in (("T1", S1), ("T2", S2)):
            out.append(f"{name} : ARRAY BITVECTOR(8) OF BITVECTOR(4);")
            for (a, b), w in tabs(s).items():
                out.append(f"ASSERT({name}[0hex{a:x}{b:x}] = 0hex{w:x});")
    for i in range(rounds + 1):
        out.append(f"L{i}, R{i} : BITVECTOR({n});")
    for i in range(rounds):
        out.append(f"U{i}, V{i}, W{i} : BITVECTOR({n});")  # S1 / S2 input (lin) or output (diff); W: joint lin mask
    terms = []
    nib = lambda var, c: f"{var}[{4 * c + 3}:{4 * c}]"
    for i in range(rounds):
        if mode == "jdiff":
            for j in range(cells):
                c1, c2 = (j + b1 // 4) % cells, (j + b2 // 4) % cells
                e = f"J[{nib(f'L{i}', j)}@{nib(f'U{i}', c1)}@{nib(f'V{i}', c2)}]"
                out.append(f"ASSERT(NOT({e} = 0hex{INVALID:x}));")
                terms.append(e)
            for j in range(n):
                srcs = [bit(f"U{i}", k) for k in range(n) if (cols[k] >> j) & 1] + \
                       [bit(f"V{i}", k) for k in range(n) if (cols[n + k] >> j) & 1]
                out.append(f"ASSERT({bit(f'L{i+1}', j)} = {xor_bits(srcs + [bit(f'R{i}', j)])});")
            out.append(f"ASSERT(R{i+1} = L{i});")
            continue
        if mode == "jlin":
            out.append(f"ASSERT(R{i} = L{i+1});")
            vmask = [xor_bits([bit(f"L{i+1}", j) for j in range(n) if (cols[k] >> j) & 1]) for k in range(2 * n)]
            for j in range(cells):
                c1, c2 = (j + b1 // 4) % cells, (j + b2 // 4) % cells
                v1 = cat_bits([vmask[4 * c1 + t] for t in (3, 2, 1, 0)])
                v2 = cat_bits([vmask[n + 4 * c2 + t] for t in (3, 2, 1, 0)])
                e = f"J[{nib(f'W{i}', j)}@{v1}@{v2}]"
                out.append(f"ASSERT(NOT({e} = 0hex{INVALID:x}));")
                terms.append(e)
            out.append(f"ASSERT(L{i} = BVXOR(R{i+1}, W{i}));")
            continue
        if mode == "diff":
            # cell c of (L_i <<< b): bits 4c..4c+3 = L_i bits (4c+t-b) mod n; outputs U (S1), V (S2)
            for tname, b, o in (("T1", b1, f"U{i}"), ("T2", b2, f"V{i}")):
                for c in range(cells):
                    a = cat_bits([bit(f"L{i}", (4 * c + t - b) % n) for t in (3, 2, 1, 0)])
                    e = f"{tname}[{a}@{o}[{4 * c + 3}:{4 * c}]]"
                    out.append(f"ASSERT(NOT({e} = 0hex{INVALID:x}));")
                    terms.append(e)
            # L_{i+1} = MW(U, V) xor R_i ; R_{i+1} = L_i
            for j in range(n):
                srcs = [bit(f"U{i}", k) for k in range(n) if (cols[k] >> j) & 1] + \
                       [bit(f"V{i}", k) for k in range(n) if (cols[n + k] >> j) & 1]
                out.append(f"ASSERT({bit(f'L{i+1}', j)} = {xor_bits(srcs + [bit(f'R{i}', j)])});")
            out.append(f"ASSERT(R{i+1} = L{i});")
        else:
            # masks: state (L_i, R_i) before round i.  F output mask = L_{i+1}; R_i = L_{i+1};
            # S-box output masks (v1, v2) = MW^T(L_{i+1}); input masks U (S1), V (S2) in rotated coordinates;
            # L_i = R_{i+1} xor (U >>> b1) xor (V >>> b2)
            out.append(f"ASSERT(R{i} = L{i+1});")
            vmask = []
            for k in range(2 * n):
                vmask.append(xor_bits([bit(f"L{i+1}", j) for j in range(n) if (cols[k] >> j) & 1]))
            for tname, base, o in (("T1", 0, f"U{i}"), ("T2", n, f"V{i}")):
                for c in range(cells):
                    vb = cat_bits([vmask[base + 4 * c + t] for t in (3, 2, 1, 0)])
                    e = f"{tname}[{o}[{4 * c + 3}:{4 * c}]@{vb}]"
                    out.append(f"ASSERT(NOT({e} = 0hex{INVALID:x}));")
                    terms.append(e)
            for j in range(n):
                out.append(f"ASSERT({bit(f'L{i}', j)} = "
                           f"{xor_bits([bit(f'R{i+1}', j), bit(f'U{i}', (j + b1) % n), bit(f'V{i}', (j + b2) % n)])});")
    out.append(f"ASSERT(NOT(L0@R0 = 0bin{'0' * (2 * n)}));")
    if fixed is not None:
        # impossible-differential / zero-correlation check: input and output fixed, any weight.
        # Unsatisfiable = no characteristic connects them (every S-box transition must be possible).
        hx = n // 4
        for name, v in zip(("L0", "R0", f"L{rounds}", f"R{rounds}"), fixed):
            out.append(f"ASSERT({name} = 0hex{v:0{hx}x});")
    if W is not None:
        out.append("total : BITVECTOR(12);")
        out.append(f"ASSERT(total = BVPLUS(12, {', '.join('0bin00000000@' + t for t in terms)}));")
        out.append(f"ASSERT(BVLE(total, 0hex{W:03x}));")
    out.append("QUERY(FALSE);")
    out.append("COUNTEREXAMPLE;")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    mode, n, b1, b2, rounds, W = sys.argv[1], *map(int, sys.argv[2:7])
    sys.stdout.write(generate(mode, n, b1, b2, rounds, W))
