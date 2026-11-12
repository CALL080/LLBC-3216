"""Which rotation pairs (b1, b2) make the LLBC round structure split, and into how many parts?

Graph on the 128 state bit positions (L: 0..63, R: 64..127, bit i of a 64-bit half = weight 2^i).
Edges = the dependencies of one round: the Feistel moves (L -> R, R -> L) and, for the new L,
the dependencies of MixWord(S1(L<<<b1), S2(L<<<b2)) on the bits of L.  Bitwise XOR of keys and
constants adds no edge.  The connected components are the finest partition of the bit positions
that every round preserves, i.e. the independent sub-ciphers.

Two graphs are built:
  cell model : every output bit of an S-box cell is taken to depend on all four input bits
               (an over-approximation of the dependencies);
  exact      : new L bit j depends on L bit i iff the Boolean function really changes when L bit i
               is flipped.  Bit i enters one S1 cell and one S2 cell; the derivative of new L bit j
               with respect to L bit i is <a1, S1(x1) xor S1(x1 xor e)> xor <a2, S2(x2) xor S2(x2 xor e')>,
               where a1, a2 are the MixWord coefficients of bit j on the two cells.  The edge exists
               iff this derivative is not identically zero (cancellations between S1 and S2 and
               inside MixWord are taken into account).

Prediction: g = gcd(64, b1, b2); if 4 | g the cipher splits into g/4 parts (a real split needs
8 | g), otherwise it does not split.
"""
from collections import Counter
from math import gcd

N = 64
S1 = [0x4, 0x7, 0x5, 0x1, 0xc, 0xb, 0xd, 0x8, 0xe, 0xf, 0x6, 0x3, 0xa, 0x9, 0x2, 0x0]
S2 = [0xe, 0xc, 0x8, 0x4, 0xb, 0xd, 0x1, 0x5, 0xf, 0x3, 0xa, 0x2, 0x9, 0x7, 0x0, 0x6]


def find(p, x):
    while p[x] != x:
        p[x] = p[p[x]]
        x = p[x]
    return x


def rol(x, r):
    r %= N
    return ((x << r) | (x >> (N - r))) & ((1 << N) - 1) if r else x


def mixword(ml, mr, b1, b2):
    t0 = rol(ml ^ mr, b1) ^ ml
    t1 = rol(t0, b2) ^ (ml ^ mr)
    u = rol(t1, b1) ^ t0
    return u ^ (rol(u, b2) ^ t1)


def parity(v):
    return bin(v).count("1") & 1


_dep = {}


def depends(p1, p2, a1, a2):
    """Does <a1, S1(x1)> xor <a2, S2(x2)> depend on the L bit that is bit p1 of cell x1 and bit p2 of cell x2?
    The other bits of x1 and x2 are L bits too; bit q of x1 is bit q - p1 + p2 of x2 when that is in 0..3."""
    key = (p1, p2, a1, a2)
    if key not in _dep:
        d = p2 - p1                     # x2 bit index = x1 bit index + d for the shared L bits
        free2 = [q for q in range(4) if not 0 <= q - d <= 3]      # bits of x2 that are not in x1
        ans = False
        for x1 in range(16):
            for f in range(1 << len(free2)):
                x2 = 0
                for q in range(4):
                    if 0 <= q - d <= 3:
                        x2 |= ((x1 >> (q - d)) & 1) << q
                for n, q in enumerate(free2):
                    x2 |= ((f >> n) & 1) << q
                v = parity(a1 & (S1[x1] ^ S1[x1 ^ (1 << p1)])) ^ parity(a2 & (S2[x2] ^ S2[x2 ^ (1 << p2)]))
                if v:
                    ans = True
                    break
            if ans:
                break
        _dep[key] = ans
    return _dep[key]


def components(b1, b2, exact):
    p = list(range(2 * N))
    union = lambda a, b: p.__setitem__(find(p, a), find(p, b))
    # Feistel: new L bit i depends on R bit i (xor); new R bit i = L bit i
    for i in range(N):
        union(i, N + i)
    # MixWord is linear in the S-box outputs: column v of the two matrices
    col1 = [mixword(1 << v, 0, b1, b2) for v in range(N)]
    col2 = [mixword(0, 1 << v, b1, b2) for v in range(N)]
    if not exact:
        for cols, r in ((col1, b1), (col2, b2)):
            for c in range(16):
                ins = [((4 * c + t) - r) % N for t in range(4)]        # L bits feeding cell c
                outs = 0
                for t in range(4):
                    outs |= cols[4 * c + t]
                for i in ins[1:]:
                    union(ins[0], i)
                for j in range(N):
                    if (outs >> j) & 1:
                        union(ins[0], j)
        return len({find(p, x) for x in range(2 * N)})
    for i in range(N):
        v1, v2 = (i + b1) % N, (i + b2) % N        # position of L bit i in the inputs of S1 and S2
        c1, p1, c2, p2 = v1 // 4, v1 % 4, v2 // 4, v2 % 4
        for j in range(N):
            a1 = sum(((col1[4 * c1 + t] >> j) & 1) << t for t in range(4))
            a2 = sum(((col2[4 * c2 + t] >> j) & 1) << t for t in range(4))
            if (a1 or a2) and depends(p1, p2, a1, a2):
                union(i, j)
    return len({find(p, x) for x in range(2 * N)})


if __name__ == "__main__":
    for name, exact in (("cell model", False), ("exact dependencies", True)):
        split, mismatch = [], 0
        for b1 in range(1, N):
            for b2 in range(1, N):
                c = components(b1, b2, exact)
                g = gcd(gcd(N, b1), b2)
                pred = g // 4 if g % 4 == 0 else 1
                mismatch += c != pred
                if c > 1:
                    split.append((b1, b2, c))
        print(f"[{name}] pairs checked:", (N - 1) ** 2, " prediction mismatches:", mismatch)
        print(f"[{name}] pairs that split:", len(split),
              " number of parts -> count:", dict(Counter(c for _, _, c in split)))
        print(f"[{name}] (3,19):", components(3, 19, exact), " (32,16):", components(32, 16, exact),
              " (5,19) [256-bit version]:", components(5, 19, exact))
    for name, S in (("S1", S1), ("S2", S2)):
        miss = [(i, o) for i in range(4) for o in range(4)
                if not any(((S[x] ^ S[x ^ (1 << i)]) >> o) & 1 for x in range(16))]
        print(f"{name}: (input bit, output bit) pairs without dependency:", miss)
