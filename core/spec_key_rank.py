"""Key rank per part if the official round function (rotations (32,16)) were keyed by the
key schedule of the published specification instead of the code's.

Spec (Li et al., IEEE IoT-J 2025, Sec. II-D, Eq. (2)-(4)), LLBC-128-128, K = K_0 || K_1:
    K_{j+2} = G^2(K_{j+1}) xor K_j xor RC_j,      RK_j = K_j || K_{j+1}       (j = 0..N-1)
    G(x) = (x[0~15]<<<3 xor x[32~47]) || (x[16~31]<<<1 xor x[48~63]) || x[0~15]<<<3 || x[16~31]<<<1
with [0~15] the most significant 16 bits.  Round j uses K_j and K_{j+1}, so 20 rounds use
K_0 .. K_20.  Constants do not change ranks and are dropped.

Part k of a 64-bit word = bits i (weight 2^i) with (i mod 16) // 4 == k, as in key_rank.py.
Compare: the code's key schedule gives rank 80 per part (key_rank.py).

Run: python spec_key_rank.py
"""
from key_rank import rank

ROUNDS = 20


def word(vars_):  # 64-bit word as a list of 64 masks, index = bit weight (LSB first)
    return list(vars_)


def rol16(v, a):
    return [v[(b - a) % 16] for b in range(16)]


def G(x):
    # 16-bit words from most significant: a = x[48..63], b = x[32..47], c = x[16..31], d = x[0..15]
    a, b, c, d = x[48:64], x[32:48], x[16:32], x[0:16]
    a3, b1 = rol16(a, 3), rol16(b, 1)
    hi = [p ^ q for p, q in zip(a3, c)]
    mid = [p ^ q for p, q in zip(b1, d)]
    return b1 + a3 + mid + hi  # LSB-first concatenation of (hi, mid, a3, b1)


def main():
    K0 = [1 << (64 + i) for i in range(64)]  # K_0 = high half of the master key
    K1 = [1 << i for i in range(64)]
    Ks = [K0, K1]
    for j in range(ROUNDS - 1):
        g2 = G(G(Ks[-1]))
        Ks.append([p ^ q for p, q in zip(g2, Ks[-2])])
    assert len(Ks) == ROUNDS + 1
    for k in range(4):
        bits = [w[i] for w in Ks for i in range(64) if (i % 16) // 4 == k]
        print(f"part {k}: {len(bits)} round-key bits (K_0..K_{ROUNDS}), rank over the master key = {rank(bits)}")
    allbits = [w[i] for w in Ks for i in range(64)]
    print("all parts together:", rank(allbits))


if __name__ == "__main__":
    main()
