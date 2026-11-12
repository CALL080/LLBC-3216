"""Driver for measure.c: part-k round keys from random master keys (official key schedule),
self-test against mini_llbc.py, then the measurements C (differential counts, linear
correlations) and D (empirical integral sums) on one 32-bit part.

Run: python measure.py test | diff | lin | integ      (results appended to results/measure_*.txt)
"""
import os, random, subprocess, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "core"))
from split_test import key_schedule, RC_ENC
from mini_llbc import pack, mini_encrypt

PART = 0
from run_min import wsl_path, USE_WSL  # noqa: E402
# compiled measure.c (gcc -O3 -march=native -fopenmp -o measure measure.c); override with LLBC_MEASURE
EXE = os.environ.get("LLBC_MEASURE", wsl_path(os.path.join(HERE, "measure")))


def part_keys(K, rounds):
    rk = key_schedule(K, rounds)
    A = [pack((r[0] << 32) | r[1], PART) for r in rk]
    B = [pack((r[2] << 32) | r[3], PART) for r in rk]
    return A, B, pack(RC_ENC << 32, PART)


def keyline(K, rounds):
    A, B, rc = part_keys(K, rounds)
    return f"{rounds} {rc:x} " + " ".join(f"{a:x} {b:x}" for a, b in zip(A, B))


def run(args, stdin):
    p = subprocess.run((["wsl", "-e"] if USE_WSL else []) + ["env", f"OMP_NUM_THREADS={os.environ.get('MEASURE_THREADS', '20')}", EXE] + args,
                       input=stdin, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    return p.stdout.strip(), p.stderr.strip()


def keys(n, seed):
    rng = random.Random(seed)
    return [[rng.getrandbits(32) for _ in range(4)] for _ in range(n)]


def selftest():
    rng = random.Random(3)
    ok = 0
    for K in keys(5, 11):
        for rounds in (1, 7, 20):
            pts = [(rng.getrandbits(16), rng.getrandbits(16)) for _ in range(20)]
            out, _ = run(["test"], keyline(K, rounds) + "\n" + "\n".join(f"{l:x} {r:x}" for l, r in pts) + "\n")
            A, B, rc = part_keys(K, rounds)
            for (l, r), line in zip(pts, out.split()):
                L, R = mini_encrypt(l, r, A, B, rc)
                ok += int(line, 16) == (L << 16) | R
    print(f"self-test: measure.c == mini_llbc.py on {ok}/300")
    return ok == 300


# best characteristics from results/*_n16_r*.log (input / output of the r rounds)
DIFF = {6: (0x40, 0x5155, 0x1555, 0x4, 24), 7: (0x0, 0x4000, 0x40, 0x0, 24),
        8: (0x8000, 0x9991, 0x80, 0x0, 28), 9: (0x0, 0x4, 0xc8c4, 0x4000, 32)}
LIN = {9: (0x303, 0x1010, 0x1010, 0x303, 12), 10: (0x808, 0x0, 0x0, 0x2020, 12),
       11: (0x9090, 0x0, 0x0, 0xe0e0, 14), 12: (0x202, 0x8080, 0x202, 0x9090, 16),
       13: (0x4040, 0x0, 0x0, 0x4040, 16), 14: (0x101, 0x1010, 0x0, 0x1010, 18),
       20: (0x808, 0x0, 0x808, 0x4020, 26)}


def main():
    mode = sys.argv[1]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    if not selftest():
        return
    if mode == "test":
        return
    log = open(os.path.join(HERE, "results", f"measure_{mode}.txt"), "a")
    say = lambda s: (print(s, flush=True), log.write(s + "\n"), log.flush())
    if mode == "diff":
        for r, (a, b, c, d, w) in DIFF.items():
            for i, K in enumerate(keys(4, 100 + r)):
                out, t = run(["diff", f"{a:x}", f"{b:x}", f"{c:x}", f"{d:x}"], keyline(K, r))
                say(f"diff r={r} w={w} key{i}: count={out} (ordered x; trail predicts 2^(32-{w}) = {2**(32-w):g}) [{t}]")
    elif mode == "lin":
        for r, (a, b, c, d, w) in LIN.items():
            for i, K in enumerate(keys(4, 200 + r)):
                out, t = run(["lin", f"{a:x}", f"{b:x}", f"{c:x}", f"{d:x}"], keyline(K, r))
                s = int(out)
                cor = s / 2**32
                lg = f"{math.log2(abs(cor)):.2f}" if s else "-inf"
                say(f"lin r={r} w={w} key{i}: correlation={cor:+.3e} (log2|c|={lg}; trail 2^-{w}; "
                    f"random 32-bit permutation ~2^-16) [{t}]")
    elif mode == "integ":
        # python measure.py integ [rounds ...]      (default: 7 8 9; number of keys: INTEG_KEYS, default 16)
        # For every fixed input bit b, measure.c (integall) returns the XOR of E(x) over the 2^31 inputs
        # with x[b] = 0.  (E is a permutation, so the half with x[b] = 1 has the same XOR.)
        # Statistic: for each b, the output bits whose XOR is zero for ALL keys; then the maximum over b.
        # A random permutation gives a zero XOR with probability 1/2 per output bit and key.
        nk = int(os.environ.get("INTEG_KEYS", "16"))
        raw = open(os.path.join(HERE, "results", "measure_integ_raw.txt"), "a")
        pc = lambda v: bin(v).count("1")
        for r in (list(map(int, sys.argv[2:])) or [7, 8, 9]):
            common = [0xFFFFFFFF] * 32
            for i, K in enumerate(keys(nk, 300 + r)):
                out, t = run(["integall"], keyline(K, r))
                sums = [int(v, 16) for v in out.split()]
                assert len(sums) == 32
                raw.write(f"r={r} key{i} K={' '.join(f'{k:08x}' for k in K)} xor sums for b=0..31: "
                          + " ".join(f"{s:08x}" for s in sums) + "\n")
                raw.flush()
                zeros = [32 - pc(s) for s in sums]
                common = [c & ~s & 0xFFFFFFFF for c, s in zip(common, sums)]
                say(f"integ r={r} key{i}: zero-sum output bits per fixed bit: min {min(zeros)}, "
                    f"mean {sum(zeros) / 32:.1f}, max {max(zeros)} [{t}]")
            counts = [pc(c) for c in common]
            best = max(range(32), key=lambda b: counts[b])
            say(f"integ r={r} ({nk} keys): number of output bits with zero sum for all keys, for b=0..31: "
                + " ".join(map(str, counts)))
            say(f"integ r={r} ({nk} keys): maximum = {counts[best]} at fixed bit {best}, mask {common[best]:08x}; "
                f"distinct masks: {sorted({f'{c:08x}' for c in common})}")


if __name__ == "__main__":
    main()
