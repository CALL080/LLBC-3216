"""Where the official LLBC code (LLBC_cipher.rar) departs from the published specification,
checked by running code rather than by reading.

Spec: Li et al., IEEE IoT-J 12(21), 2025, Sec. II (Alg. 1, Fig. 2, Eq. (1)-(4), Tables III-V).
Code: LLBC_cipher/C/LLBC.cpp (transcribed in split_test.py) and LLBC_cipher/Hardware/LLBCRound.v
(transcribed below with Verilog's MSB-first [0:n-1] bit numbering).

Run: python spec_vs_impl.py
"""
import random

from split_test import make_cipher, key_schedule, RC_ENC, RC_KS, M32, M16

M64 = (1 << 64) - 1
S1 = [0x4, 0x7, 0x5, 0x1, 0xc, 0xb, 0xd, 0x8, 0xe, 0xf, 0x6, 0x3, 0xa, 0x9, 0x2, 0x0]
S2 = [0xe, 0xc, 0x8, 0x4, 0xb, 0xd, 0x1, 0x5, 0xf, 0x3, 0xa, 0x2, 0x9, 0x7, 0x0, 0x6]


# ---------------- Verilog transcription (LLBCRound.v, LLBC.v) ----------------
def sl(x, n, i, j):
    """x[i:j] of an [0:n-1] vector (bit 0 = MSB)."""
    return (x >> (n - 1 - j)) & ((1 << (j - i + 1)) - 1)


def cat(*parts):
    v = 0
    for val, w in parts:
        v = (v << w) | val
    return v


def sub(x, s):
    return cat(*[(s[sl(x, 64, 4 * i, 4 * i + 3)], 4) for i in range(16)])


def v_round(cl, cr, ck):
    perm_left = cat((sl(cl, 64, 32, 63), 32), (sl(cl, 64, 0, 31), 32))
    perm_right = cat((sl(cl, 64, 16, 63), 48), (sl(cl, 64, 0, 15), 16))
    sL, sR = sub(perm_left, S1), sub(perm_right, S2)
    a = sL ^ sR
    b = cat((sl(a, 64, 32, 63), 32), (sl(a, 64, 0, 31), 32)) ^ sL
    c = cat((sl(b, 64, 16, 63), 48), (sl(b, 64, 0, 15), 16)) ^ a
    d = cat((sl(c, 64, 32, 63), 32), (sl(c, 64, 0, 31), 32)) ^ b
    e = cat((sl(d, 64, 16, 63), 48), (sl(d, 64, 0, 15), 16)) ^ c
    f = e ^ d
    kl, kr = sl(ck, 128, 0, 63), sl(ck, 128, 64, 127)
    p1, p2, p3, p4 = (sl(kr, 64, 16 * i, 16 * i + 15) for i in range(4))
    r3 = lambda x: cat((sl(x, 16, 3, 15), 13), (sl(x, 16, 0, 2), 3))
    r1 = lambda x: cat((sl(x, 16, 1, 15), 15), (sl(x, 16, 0, 0), 1))
    rk11, rk21 = r3(p1), r1(p2)
    tk11, tk21, tk31, tk41 = rk11 ^ p3, rk21 ^ p4, rk11, rk21
    rk12, rk22 = r3(tk11), r1(tk21)
    tk12, tk22, tk32, tk42 = rk12 ^ tk31, rk22 ^ tk41, rk12, rk22
    kr1, kr2 = cat((tk12, 16), (tk22, 16)), cat((tk32, 16), (tk42, 16))
    kl1, kl2 = sl(kl, 64, 0, 31), sl(kl, 64, 32, 63)
    nk = cat((kl1 ^ kr1, 32), (kl2 ^ kr2 ^ 0x243F6A88, 32), (kr1 ^ 0x85A308D3, 32), (kr2, 32))
    nl = f ^ cr ^ sl(nk, 128, 64, 127)
    nr = cl ^ sl(nk, 128, 0, 63) ^ cat((0x243F6A88, 32), (0, 32))
    return nl, nr, nk


def v_encrypt(p, k, rounds=20):
    L, R = p >> 64, p & M64
    for _ in range(rounds):
        L, R, k = v_round(L, R, k)
    return (L << 64) | R


def words(x):
    return [(x >> (96 - 32 * i)) & M32 for i in range(4)]


def unwords(w):
    return (w[0] << 96) | (w[1] << 64) | (w[2] << 32) | w[3]


# ---------------- spec pieces (Eq. (3), MSB-first [0~15] = top 16 bits) ----------------
def rol16(x, a):
    return ((x << a) | (x >> (16 - a))) & M16


def G(k):
    a, b, c, d = (k >> 48) & M16, (k >> 32) & M16, (k >> 16) & M16, k & M16
    return (((rol16(a, 3) ^ c) << 48) | ((rol16(b, 1) ^ d) << 32) | (rol16(a, 3) << 16) | rol16(b, 1))


def main():
    rng = random.Random(2026)
    enc = make_cipher(32, 16)
    print("== 1. Verilog = C ? (transcriptions, 2000 random (P, K), 20 rounds)")
    tv = v_encrypt(0, (1 << 128) - 1)
    print("   Verilog, P=0, K=FF..FF ->", " ".join(f"{w:08X}" for w in words(tv)), "(C: 0C40B9A6 69053EC4 F1DD3EEA 85366E1E)")
    same = 0
    for _ in range(2000):
        p, k = rng.getrandbits(128), rng.getrandbits(128)
        same += words(v_encrypt(p, k)) == enc(words(p), key_schedule(words(k)))
    print(f"   agree: {same}/2000")
    tb = 0x0123456789ABCDEFFEDCBA9876543210
    print("   testbench vector P=K=0123..3210 ->", " ".join(f"{w:08X}" for w in words(v_encrypt(tb, tb))),
          "(LLBC_tb.v prints nothing and checks nothing)")

    print("== 2. Key schedule of the code vs Eq. (2)-(4)")
    K = [rng.getrandbits(32) for _ in range(4)]
    rk = key_schedule(K)
    A = lambda j: (rk[j][0] << 32) | rk[j][1]   # XORed into L (-> T, -> R_{i+1})
    B = lambda j: (rk[j][2] << 32) | rk[j][3]   # XORed after MixWord
    K0, K1 = (K[0] << 32) | K[1], (K[2] << 32) | K[3]
    print("   round-0 key equals the master key (Eq. (4): RK_0 = K_0 || K_1)?", (A(0), B(0)) == (K0, K1))
    ok = all(B(j + 1) == G(G(B(j))) ^ (RC_KS << 32) for j in range(19)) and B(0) == G(G(K1)) ^ (RC_KS << 32)
    print("   right half: B_{j+1} = G^2(B_j) xor (RC_KS || 0) for all rounds?", ok)
    ok = all(A(j + 1) == A(j) ^ G(G(B(j))) ^ RC_ENC for j in range(19)) and A(0) == K0 ^ G(G(K1)) ^ RC_ENC
    print("   left half:  A_{j+1} = A_j xor G^2(B_j) xor (0 || RC_ENC) for all rounds?", ok)
    flips, dep = 200, 0
    for _ in range(flips):
        K2 = list(K); K2[rng.randrange(2)] ^= 1 << rng.randrange(32)
        rk2 = key_schedule(K2)
        dep += any(rk2[j][2:] != rk[j][2:] for j in range(20))
    print(f"   flipping a bit of K_0 changes the post-MixWord keys B_j: {dep}/{flips} times "
          "(B_0..B_19 depend only on the 64 bits K_1)")
    spec_K2 = G(G(K1)) ^ K0  # Eq. (2) without RC_0 (where RC_0 goes in 64 bits is not stated)
    print("   code's A_0 = Eq. (2)'s K_2 up to the constant?", (A(0) ^ RC_ENC) == spec_K2)
    print("   code's B_0 = Eq. (2)'s K_3 = G^2(K_2) xor K_1 up to constants?",
          (B(0) ^ (RC_KS << 32)) == (G(G(spec_K2)) ^ K1))

    print("== 3. Round constants")
    print(f"   code: RC_ENC = {RC_ENC:08X} in every round (spec: RC_i = RC_0 <<< i); "
          f"RC_KS = {RC_KS:08X} = Table III's RC_0 for LLBC-128-256")

    print("== 4. Round function of the code: S2 input")
    L = rng.getrandbits(64); R = rng.getrandbits(64)
    st = [L >> 32, L & M32, R >> 32, R & M32]
    a = enc(st, [[0, 0, 0, 0]])
    b = enc(st, [[rng.getrandbits(32), rng.getrandbits(32), 0, 0]])
    print("   changing the round key on the T path changes L_{i+1}?", a[:2] != b[:2],
          "(spec Fig. 2 / Alg. 1 line 6: S2 takes T_i = L_i xor RK xor RC, so it should)")


if __name__ == "__main__":
    main()
