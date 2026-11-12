"""How many master-key bits does each of the four 32-bit parts depend on?

The official key schedule (LLBC_cipher/C/LLBC.cpp, key_schedule) is affine over GF(2):
16-bit rotations by 3 and 1, XORs and constants. We track every state bit as a 128-bit
mask over the master key and compute, for each part k, the GF(2) rank of the round-key
bits that part k uses (bits b of each 32-bit round-key word with (b mod 16) // 4 == k).

Rounds: the round keys of rounds 0..R-1. The first round's (rk0, rk1) only goes into the
new right half and the last round's (rk2, rk3) is a final whitening, so the rank is an upper
bound on the effective key length of a part, not the exact value.
"""
import sys

PART = [0x000F000F, 0x00F000F0, 0x0F000F00, 0xF000F000]


def sym_word(w):
    # master key word w (0..3), bit b -> variable index 32*w + b
    return [1 << (32 * w + b) for b in range(32)]


def rol16(v, a):
    return [v[(b - a) % 16] for b in range(16)]


def xor(a, b):
    return [x ^ y for x, y in zip(a, b)]


def key_schedule(rounds):
    key = [sym_word(w) for w in range(4)]  # constants are ignored: they do not change the rank
    out = []
    for _ in range(rounds):
        k2, k3 = key[2], key[3]
        x3, x2, x1, x0 = k2[16:], k2[:16], k3[16:], k3[:16]
        for _ in range(2):
            tmp = rol16(x3, 3); x3 = xor(tmp, x1); x1 = tmp
            tmp = rol16(x2, 1); x2 = xor(tmp, x0); x0 = tmp
        key[2] = x2 + x3
        key[3] = x0 + x1
        key[0] = xor(key[0], key[2])
        key[1] = xor(key[1], key[3])
        out.append([list(w) for w in key])
    return out


def rank(vectors):
    basis = {}
    for v in vectors:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return len(basis)


def main(rounds=20):
    rks = key_schedule(rounds)
    print(f"rounds = {rounds}")
    all_vecs = []
    for k in range(4):
        vecs = [rk[w][b] for rk in rks for w in range(4) for b in range(32) if (PART[k] >> b) & 1]
        all_vecs += vecs
        print(f"part {k}: {len(vecs)} round-key bits, rank over the master key = {rank(vecs)}")
    print(f"all parts together: rank = {rank(all_vecs)}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
