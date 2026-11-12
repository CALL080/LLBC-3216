"""Exhaustive counts (differential) / correlations (linear) for the input and output of the best
characteristics in results/<mode>_n16_r<r>.log, 4 keys each (key schedule of the reference
implementation, part 0; master keys from random.Random(700 + r)).

mode: jdiff / jlin = joint S-box model, diff / lin = independent S-box model.
The count is the number of inputs x (ordered) with E(x) xor E(x xor din) = dout, so it includes every
characteristic with the same input and output differences; the correlation is that of the linear hull.

Run: MEASURE_THREADS=20 python measure_best.py jdiff 6 7 8 9 10      (results/measure_best_jdiff.txt)
Table 5 of the paper: diff 6; jdiff 3 5; jlin 5 6.   Table 8: jdiff 6..10; jlin 7..11.
"""
import math, os, re, sys
import measure as m

HERE = os.path.dirname(os.path.abspath(__file__))
mode, rounds = sys.argv[1], list(map(int, sys.argv[2:]))
kind = "diff" if mode in ("diff", "jdiff") else "lin"
log = open(os.path.join(HERE, "results", f"measure_best_{mode}.txt"), "a")
for r in rounds:
    sat = [l for l in open(os.path.join(HERE, "results", f"{mode}_n16_r{r}.log")) if ": sat" in l]
    line = min(sat, key=lambda l: int(re.search(r"checked weight=(\d+)", l).group(1)))
    w = int(re.search(r"checked weight=(\d+)", line).group(1))
    v = dict((k, int(x, 16)) for k, x in re.findall(r"([LR]\d+)=([0-9a-f]+)", line))
    args = [f"{v['L0']:x}", f"{v['R0']:x}", f"{v[f'L{r}']:x}", f"{v[f'R{r}']:x}"]
    for i, K in enumerate(m.keys(4, 700 + r)):
        out, t = m.run([kind] + args, m.keyline(K, r))
        if kind == "diff":
            s = f"{mode} r={r} w={w} key{i}: count={out} (ordered x; predicted 2^(32-{w}) = {2**(32-w):g}) [{t}]"
        else:
            c = int(out) / 2**32
            s = f"{mode} r={r} w={w} key{i}: correlation={c:+.3e} (log2|c|={math.log2(abs(c)) if c else float('-inf'):.2f}; trail 2^-{w}) [{t}]"
        print(s, flush=True); log.write(s + "\n"); log.flush()
