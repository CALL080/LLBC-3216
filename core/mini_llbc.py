"""Each part of the official LLBC code is itself a 32-bit LLBC-type cipher.

Part k of a 64-bit half = cells k, k+4, k+8, k+12 (bits with (i mod 16)//4 == k).
Pack them as a 16-bit word w_k with nibble j = cell k + 4j.  Then rotating the 64-bit half by
16 (resp. 32) bits moves cell c to c+4 (resp. c+8), i.e. nibble j to j+1 (resp. j+2): a 4-bit
(resp. 8-bit) rotation of w_k.  Claim checked here: part k of the full 20-round cipher equals the
cipher with 16-bit halves, the same S-boxes and MixWord with rotations (8, 4), keyed with part k
of every round key and constant.

Run: python mini_llbc.py
"""
import random

from split_test import S1, S2, RC_ENC, key_schedule, make_cipher

M16 = 0xFFFF


def pack(x64, k):
    return sum(((x64 >> (4 * (k + 4 * j))) & 0xF) << (4 * j) for j in range(4))


def rol16(x, r):
    return ((x << r) | (x >> (16 - r))) & M16


def sbox16(x, s):
    return sum(s[(x >> (4 * j)) & 0xF] << (4 * j) for j in range(4))


def mini_encrypt(L, R, A, B, rc):
    """A, B, rc: per-round 16-bit keys/constant for this part."""
    for j in range(len(A)):
        s1 = sbox16(rol16(L, 8), S1)
        s2 = sbox16(rol16(L, 4), S2)
        t = s1 ^ s2
        m0 = rol16(t, 8) ^ s1
        m1 = rol16(m0, 4) ^ t
        m2 = rol16(m1, 8) ^ m0
        m3 = rol16(m2, 4) ^ m1
        L, R = m2 ^ m3 ^ R ^ B[j], L ^ A[j] ^ rc
    return L, R


def main():
    rng = random.Random(20261026)
    enc = make_cipher(32, 16)
    ok = 0
    trials = 1000
    for _ in range(trials):
        K = [rng.getrandbits(32) for _ in range(4)]
        P = [rng.getrandbits(32) for _ in range(4)]
        rk = key_schedule(K)
        C = enc(P, rk)
        L, R = (P[0] << 32) | P[1], (P[2] << 32) | P[3]
        CL, CR = (C[0] << 32) | C[1], (C[2] << 32) | C[3]
        good = True
        for k in range(4):
            A = [pack((r[0] << 32) | r[1], k) for r in rk]
            B = [pack((r[2] << 32) | r[3], k) for r in rk]
            got = mini_encrypt(pack(L, k), pack(R, k), A, B, pack(RC_ENC << 32, k))
            good &= got == (pack(CL, k), pack(CR, k))
        ok += good
    print(f"part k of LLBC(32,16) == 32-bit LLBC-type cipher with rotations (8,4): {ok}/{trials}")
    print("constant RC64 = RC_ENC||0 per part:", [f"{pack(RC_ENC << 32, k):04X}" for k in range(4)])


if __name__ == "__main__":
    main()
