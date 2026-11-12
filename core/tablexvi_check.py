"""Check the five test vectors of Table XVI of the original LLBC paper.

Prints one line per vector: "<match for rotations (32,16)> <match for rotations (3,19)>".
The official implementation uses (32,16); the paper's own parameters are (3,19), so the
vectors are reproduced by (32,16) only.  Output equals tableXVI_check_log.txt.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from split_test import make_cipher, key_schedule, M32  # noqa: E402

F = M32
VECTORS = [
    ([0, 0, 0, 0], [0, 0, 0, 0], "6A7BE51D 43DC7B8E A328EE23 63B764BF"),
    ([0, 0, 0, 0], [F, F, F, F], "0C40B9A6 69053EC4 F1DD3EEA 85366E1E"),
    ([F, F, F, F], [0, 0, 0, 0], "601DDECB 1B7346B7 9A6BB8FF 32F6FF58"),
    ([F, F, F, F], [F, F, F, F], "A4FF76B7 91CE1FF2 57613730 E8EABC4F"),
    ([0x01234567, 0x89ABCDEF, 0xFEDCBA98, 0x76543210],
     [0x01234567, 0x89ABCDEF, 0xFEDCBA98, 0x76543210], "7A03B217 F2F19162 A32A1105 DB69CD28"),
]

if __name__ == "__main__":
    encs = {rot: make_cipher(*rot) for rot in [(32, 16), (3, 19)]}
    for p, k, expect in VECTORS:
        want = [int(x, 16) for x in expect.split()]
        rk = key_schedule(k)
        print(*[encs[rot](p, rk) == want for rot in [(32, 16), (3, 19)]])
