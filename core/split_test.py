"""Does LLBC as implemented (rotations 32/16) split into four independent 32-bit ciphers?

Transcription of LLBC_cipher/C/LLBC.cpp (sub_cell, m_function, encrypt, key_schedule)
on 32-bit words; ROL64 acts on the word pair (x[0] = high word).
Part k (k = 0..3) of each 32-bit word = bits b with (b mod 16) // 4 == k.
Test: flip random bits of one part of the plaintext and check which parts of the
ciphertext change. Control: the same cipher with S-box/MixWord rotations (3, 19).
"""
import random

S1 = [0x4, 0x7, 0x5, 0x1, 0xc, 0xb, 0xd, 0x8, 0xe, 0xf, 0x6, 0x3, 0xa, 0x9, 0x2, 0x0]
S2 = [0xe, 0xc, 0x8, 0x4, 0xb, 0xd, 0x1, 0x5, 0xf, 0x3, 0xa, 0x2, 0x9, 0x7, 0x0, 0x6]
RC_ENC, RC_KS = 0x243F6A88, 0x85A308D3
M32, M16 = 0xFFFFFFFF, 0xFFFF
PART = [0x000F000F, 0x00F000F0, 0x0F000F00, 0xF000F000]


def rol64(x, a):
    v = (x[0] << 32) | x[1]
    v = ((v << a) | (v >> (64 - a))) & ((1 << 64) - 1)
    return [v >> 32, v & M32]


def sbox_word(w, s):
    return sum(s[(w >> (4 * i)) & 0xF] << (4 * i) for i in range(8))


def make_cipher(r1, r2):
    def sub_cell(st):
        t1, t2 = rol64(st[0:2], r1), rol64(st[0:2], r2)
        return [sbox_word(t1[0], S1), sbox_word(t1[1], S1), sbox_word(t2[0], S2), sbox_word(t2[1], S2)]

    def m_function(left, right):
        tl, tr = list(left), list(right)
        tr = [tr[0] ^ tl[0], tr[1] ^ tl[1]]
        t2 = list(tr)
        tr = rol64(tr, r1)
        tl = [tl[0] ^ tr[0], tl[1] ^ tr[1]]
        t3 = list(tl)
        tl = rol64(tl, r2)
        t2 = [t2[0] ^ tl[0], t2[1] ^ tl[1]]
        tr = list(t2)
        t2 = rol64(t2, r1)
        t3 = [t3[0] ^ t2[0], t3[1] ^ t2[1]]
        t4 = list(t3)
        t3 = rol64(t3, r2)
        tr = [tr[0] ^ t3[0], tr[1] ^ t3[1]]
        return [t4[0] ^ tr[0], t4[1] ^ tr[1]]

    def encrypt(state, rk):
        st = list(state)
        for j in range(len(rk)):
            left = [st[0] ^ rk[j][0] ^ RC_ENC, st[1] ^ rk[j][1]]
            right_start = [st[2], st[3]]
            s = sub_cell(st)
            m = m_function(s[0:2], s[2:4])
            m = [m[0] ^ rk[j][2], m[1] ^ rk[j][3]]
            st = [m[0] ^ right_start[0], m[1] ^ right_start[1], left[0], left[1]]
        return st

    return encrypt


def rol16(x, a):
    return ((x << a) | (x >> (16 - a))) & M16


def key_schedule(key, rounds=20):
    key = list(key)
    out = []
    for _ in range(rounds):
        k2, k3 = key[2], key[3]
        x3, x2, x1, x0 = k2 >> 16, k2 & M16, k3 >> 16, k3 & M16
        for _ in range(2):
            tmp = rol16(x3, 3); x3 = tmp ^ x1; x1 = tmp
            tmp = rol16(x2, 1); x2 = tmp ^ x0; x0 = tmp
        key[2] = (x3 << 16) | x2
        key[3] = (x1 << 16) | x0
        key[0] ^= key[2]
        key[1] ^= key[3]
        key[1] ^= RC_ENC
        key[2] ^= RC_KS
        out.append(list(key))
    return out


def run(r1, r2, trials=2000, seed=1):
    rng = random.Random(seed)
    enc = make_cipher(r1, r2)
    leaks = [[0] * 4 for _ in range(4)]  # leaks[k][j]: flips in part k changed part j
    for t in range(trials):
        rk = key_schedule([rng.getrandbits(32) for _ in range(4)])
        p = [rng.getrandbits(32) for _ in range(4)]
        k = t % 4
        d = [rng.getrandbits(32) & PART[k] for _ in range(4)]
        if not any(d):
            d[0] = PART[k] & 0x1 or PART[k]
        q = [p[i] ^ d[i] for i in range(4)]
        c, c2 = enc(p, rk), enc(q, rk)
        for j in range(4):
            if any((c[i] ^ c2[i]) & PART[j] for i in range(4)):
                leaks[k][j] += 1
    return leaks


if __name__ == "__main__":
    # sanity: the official test vector (P = 0, K = all ones), printed for comparison with the C program
    enc = make_cipher(32, 16)
    c = enc([0, 0, 0, 0], key_schedule([M32] * 4))
    print("official-style 20 rounds, P=0, K=FF..FF ->", " ".join(f"{w:08X}" for w in c))
    for r1, r2 in [(32, 16), (3, 19)]:
        L = run(r1, r2)
        print(f"rotations ({r1},{r2}): rows = part flipped, cols = parts of C that changed (out of 500 each)")
        for k in range(4):
            print("  part", k, L[k])
