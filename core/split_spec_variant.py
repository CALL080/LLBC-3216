"""Does the four-way split need anything besides the rotation amounts (32, 16)?

Variant = the published round function (Alg. 1 / Fig. 2: S2 takes the keyed T_i, RC_i = RC_0 <<< i,
M_L = S2 output, R_{i+1} = T_i) with rotations (32, 16) and independent random round keys.
If no flip leaks into another part, the split comes from the rotations alone, whatever the key
schedule, the constants or the M_L/M_R order.
"""
import random
from split_test import S1, S2, PART, sbox_word, M32

M64 = (1 << 64) - 1
RC0 = 0x243F6A88


def rol(x, r):
    return ((x << r) | (x >> (64 - r))) & M64


def sub(x, s):
    return (sbox_word(x >> 32, s) << 32) | sbox_word(x & M32, s)


def mixword(ml, mr, b1, b2):  # Eq. (1)
    t0 = rol(ml ^ mr, b1) ^ ml
    t1 = rol(t0, b2) ^ (ml ^ mr)
    u = rol(t1, b1) ^ t0
    return u ^ (rol(u, b2) ^ t1)


def enc(L, R, rks, b1, b2):
    for i, (ka, kb) in enumerate(rks):
        rc = ((RC0 << (i % 32)) | (RC0 >> (32 - i % 32))) & M32 if i % 32 else RC0
        T = L ^ ka ^ (rc << 32)
        L, R = mixword(sub(rol(T, b2), S2), sub(rol(L, b1), S1), b1, b2) ^ R ^ kb, T
    return L, R


def part64(x, k):
    m = (PART[k] << 32) | PART[k]
    return x & m


def run(b1, b2, trials=2000, seed=3):
    rng = random.Random(seed)
    leaks = [[0] * 4 for _ in range(4)]
    for t in range(trials):
        rks = [(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(20)]
        L, R, k = rng.getrandbits(64), rng.getrandbits(64), t % 4
        dL, dR = part64(rng.getrandbits(64), k) or part64(M64, k), part64(rng.getrandbits(64), k)
        a, b = enc(L, R, rks, b1, b2), enc(L ^ dL, R ^ dR, rks, b1, b2)
        for j in range(4):
            if part64(a[0] ^ b[0], j) or part64(a[1] ^ b[1], j):
                leaks[k][j] += 1
    return leaks


if __name__ == "__main__":
    for b1, b2 in [(32, 16), (3, 19)]:
        print(f"published round function, rotations ({b1},{b2}), random round keys; rows = part flipped (500 each)")
        for k, row in enumerate(run(b1, b2)):
            print("  part", k, row)
