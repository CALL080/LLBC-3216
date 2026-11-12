"""Decryption with known plaintexts (Section 6.1): estimate and small-scale experiment.

With N random known plaintexts, a row of the full table of one part (2^n rows) is filled with
probability 1 - exp(-N / 2^n).  A ciphertext can be decrypted when the values of all four parts
are in the tables: probability (1 - exp(-N / 2^n))^4.

The experiment replaces the four 32-bit part ciphers by four random permutations of 2^12 values
(the attack uses nothing but the split) and counts how many of 2000 new ciphertexts can be
decrypted from N = 4 * 2^12 known plaintexts.  Seed: random.Random(20261112 + trial).
"""
import math
import random

for e in (32, 33, 34, 35, 36):
    print(f"n = 32, N = 2^{e}: (1 - exp(-N/2^32))^4 = {(1 - math.exp(-2.0 ** (e - 32))) ** 4:.4f}")

n, tests = 12, 2000
size = 1 << n
N = 4 * size
print(f"\nsmall-scale experiment: parts of {n} bits, N = 4 * 2^{n} = {N} known plaintexts, {tests} test ciphertexts")
print(f"formula: {(1 - math.exp(-N / size)) ** 4:.4f}")
rates = []
for trial in range(10):
    rng = random.Random(20261112 + trial)
    perms = []
    for _ in range(4):
        p = list(range(size))
        rng.shuffle(p)
        perms.append(p)
    known = [set() for _ in range(4)]            # ciphertext part values seen among the known pairs
    for _ in range(N):
        for k in range(4):
            known[k].add(perms[k][rng.randrange(size)])
    ok = sum(all(perms[k][rng.randrange(size)] in known[k] for k in range(4)) for _ in range(tests))
    rates.append(ok / tests)
    print(f"  trial {trial}: {ok}/{tests} = {ok / tests:.4f}")
print(f"mean over 10 trials: {sum(rates) / len(rates):.4f}  (min {min(rates):.4f}, max {max(rates):.4f})")
