"""Joint S-box model of the reference implementation (Section 5.1, Table 4).

In the reference implementation S1 and S2 read the same nibble x of L without any key, so the
nibble map x -> (S1(x), S2(x)) is one function from 4 bits to 8 bits.  This script prints
  - the joint difference distribution: for every nonzero input difference, the number of
    possible output-difference pairs and their counts;
  - the comparison with the independent model for input difference 4 (Table 4);
  - the nonzero linear masks (input, S1 output, S2 output) with correlation +-1;
  - the key-averaged one-round differential probability when a uniform key nibble is XORed to
    the S2 input (as in the specification): it equals the product of the two S-box probabilities.
"""
from fractions import Fraction

S1 = [0x4, 0x7, 0x5, 0x1, 0xc, 0xb, 0xd, 0x8, 0xe, 0xf, 0x6, 0x3, 0xa, 0x9, 0x2, 0x0]
S2 = [0xe, 0xc, 0x8, 0x4, 0xb, 0xd, 0x1, 0x5, 0xf, 0x3, 0xa, 0x2, 0x9, 0x7, 0x0, 0x6]
par = lambda v: bin(v).count("1") & 1


def ddt(S):
    t = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for x in range(16):
            t[a][S[x] ^ S[x ^ a]] += 1
    return t


D1, D2 = ddt(S1), ddt(S2)
joint = {a: {} for a in range(16)}
for a in range(16):
    for x in range(16):
        k = (S1[x] ^ S1[x ^ a], S2[x] ^ S2[x ^ a])
        joint[a][k] = joint[a].get(k, 0) + 1

print("joint DDT: input difference -> number of output pairs, set of counts (out of 16)")
for a in range(1, 16):
    print(f"  {a:x}: {len(joint[a])} pairs, counts {sorted(set(joint[a].values()))}, "
          f"both output differences nonzero: {all(b1 and b2 for b1, b2 in joint[a])}")
print("every nonzero input difference: exactly 8 pairs, each with count 2 (probability 2/16, weight 3):",
      all(len(joint[a]) == 8 and set(joint[a].values()) == {2} for a in range(1, 16)))

a = 4
indep = [(b1, b2) for b1 in range(16) for b2 in range(16) if D1[a][b1] and D2[a][b2]]
print(f"\nTable 4 (input difference {a}): independent model allows {len(indep)} pairs, joint model {len(joint[a])}")
for b in ((4, 1), (8, 4), (3, 4), (4, 6)):
    print(f"  output differences {b}: independent {Fraction(D1[a][b[0]] * D2[a][b[1]], 256)}, "
          f"actual {Fraction(joint[a].get(b, 0), 16)}")
print("  joint pairs:", sorted(joint[a]))

print("\nlinear masks (input, S1 output, S2 output) with correlation +-1, not all zero:")
n = 0
for g in range(16):
    for h1 in range(16):
        for h2 in range(16):
            if g | h1 | h2:
                w = sum(1 - 2 * (par(g & x) ^ par(h1 & S1[x]) ^ par(h2 & S2[x])) for x in range(16))
                if abs(w) == 16:
                    n += 1
                    print(f"  ({g:x}, {h1:x}, {h2:x}): Walsh sum {w:+d}")
print("number of such masks:", n)

ok = True
for a in range(1, 16):
    for b1 in range(16):
        for b2 in range(16):
            cnt = sum(1 for k in range(16) for x in range(16)
                      if S1[x] ^ S1[x ^ a] == b1 and S2[x ^ k] ^ S2[x ^ k ^ a] == b2)
            ok &= Fraction(cnt, 256) == Fraction(D1[a][b1] * D2[a][b2], 256)
print("\nkey nibble on the S2 input, averaged over the key: one-round differential probability "
      "= product of the two S-box probabilities for all (a, b1, b2):", ok)
