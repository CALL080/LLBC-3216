# Artifact: full-round attacks on the reference implementation of LLBC

Code and logs for the paper "Four Ciphers in One: Full-Round Attacks on the Reference Implementation
of LLBC". The reference implementation of the LLBC block cipher uses the rotation amounts (32,16)
and splits into four independent 32-bit ciphers.

Section, table, and figure numbers below are those of the paper.
In code comments and logs, "official code" or "official implementation" means the reference
implementation, that is, the code released by the designers.

## Requirements

- Python 3 (standard library only).
- gcc with OpenMP (for `attack/` and `trail/measure.c`).
- STP 2.3.3, official static build (only for re-running the searches in `trail/`).
- Windows: WSL. The Linux tools (STP, `measure`) are called through `wsl -e bash -lc`; on Linux they are called directly.
- Environment variables: `LLBC_STP` (STP binary, default `~/stp-official/stp`), `LLBC_MEASURE`
  (compiled `trail/measure`, default `trail/measure`), `MEASURE_THREADS` (default 20).

## Designers' code (not included)

The designers' code has no redistribution licence, so it is not part of this artifact.
`attack/exp_a.cpp` and `attack/exp_b.cpp` include `../bench/llbc_lib.cpp`. Create it as follows:

1. Download `LLBC_cipher.rar` from GitHub `13est-Wishes/LLBC-cipher`
   (SHA-1 `a9f914fadd3027b4c1232789ab81e1a194ae02f6`; Wayback copy:
   https://web.archive.org/web/20261026114738/https://raw.githubusercontent.com/13est-Wishes/LLBC-cipher/main/LLBC_cipher.rar).
2. Extract `LLBC_cipher/C/LLBC.cpp`.
3. Rename its `int main()` to `int main_orig()` and save as `bench/llbc_lib.cpp`:
   `mkdir bench; sed 's/^int main()/int main_orig()/' LLBC.cpp > bench/llbc_lib.cpp`
   (this is the only edit).
4. Build and run, e.g. `cd attack && g++ -O3 -march=native -fopenmp -w -o exp_a exp_a.cpp && ./exp_a 1000 20261026`
   (`exp_b` likewise, `./exp_b 32 20261026`). Experiment A needs about 3 minutes on 20 threads.

All other scripts use Python transcriptions of the designers' code (`core/split_test.py`) and do not need it.

## Paper result -> script -> log

Run the `core/` scripts from the artifact folder, e.g. `python core/split_test.py`.

| Paper result | Script | Log |
|---|---|---|
| Test vectors of Table XVI of the proposal (Section 2.3, Table 2) | `core/tablexvi_check.py` | `core/tableXVI_check_log.txt` |
| C and Verilog logic agree, 2000/2000 (Section 2.1); differences from the specification (Table 2) | `core/spec_vs_impl.py` | `core/spec_vs_impl_log.txt` |
| Proposition 1: 3969 rotation pairs, exact dependencies and cell model (Section 3.2) | `core/split_condition.py` (about 1 minute) | `core/split_condition_log.txt` |
| The split on a transcription of the designers' code (Section 3) | `core/split_test.py` | stdout (the check on the C code itself is Experiment A) |
| Corollary 1: each part is a 32-bit LLBC-type cipher with rotations (8,4), 1000/1000 | `core/mini_llbc.py` | `core/mini_llbc_log.txt` |
| Remark 1: the split also holds with the round function of the specification | `core/split_spec_variant.py` | `core/split_spec_variant_log.txt` |
| Experiment A: 2^32 chosen plaintexts, 175.5 s, 1000/1000 (Section 4.2) | `attack/exp_a.cpp` | `attack/exp_a_log.txt`, `attack/time_a.txt` |
| Table 3: key ranks 32..80 and 80/96/112/128 (Section 4.3) | `core/key_rank.py` | stdout |
| Experiment B: reduced key recovery (Section 4.3) | `attack/exp_b.cpp` | `attack/exp_b_log.txt`, `attack/time_b.txt` |
| Remark 3: rank 118 with the key schedule of the specification | `core/spec_key_rank.py` | `core/spec_key_rank_log.txt` |
| Table 4 and Section 5.1: joint S-box model (8 output pairs per input difference, 7 masks with correlation +-1, key-averaged identity) | `core/joint_sbox.py` | `core/joint_sbox_log.txt` |
| Section 6.1: decryption with known plaintexts (formula and small-scale experiment) | `core/known_plaintext.py` | `core/known_plaintext_log.txt` |

Experiment A compares 1000 target ciphertexts during the 2^32 queries; it does not store the 64-GiB tables.

### Section 5 (`trail/`, Tables 5-8)

The model generator is `trail/gen_model.py` (`python gen_model.py {diff|lin|jdiff|jlin} n b1 b2 rounds W > model.cvc`;
`j*` = joint S-box model, `diff`/`lin` = independent S-box model; n=16 with (8,4) is one 32-bit part,
n=64 with (32,16) is the full cipher).
`check_trail.py` independently re-checks every characteristic printed by STP.
The `.cvc`/`.out` model files are not included; they are regenerated into `trail/models/`.

Build the measurement program once: `cd trail && gcc -O3 -march=native -fopenmp -o measure measure.c`.
`python measure.py test` compares it with `core/mini_llbc.py` (300/300).
All measurements run over all 2^32 inputs of part 0, with round keys from the key schedule of the
reference implementation and master keys from seeded `random.Random` (seeds are in the scripts).
A differential count is the number of inputs x with E(x) xor E(x xor din) = dout, so it includes every
characteristic with the same input and output differences; a correlation is that of the linear hull.

| Result | Command (in `trail/`) | Results file(s) in `trail/results/` |
|---|---|---|
| Table 7: minimum weights (joint model `jdiff`/`jlin`; independent model `diff`/`lin`) | `run_min.py`, `par_min.py` (parallel over bounds W), `seq_par.py`, `cap_run.py` (time-capped), `chain.py` (round chain) | `jdiff_n16_r*.log`, `jlin_n16_r*.log`, `diff_n16_r*.log`, `lin_n16_r*.log`, `lin_n64_r*.log`, `par_jdiff7.txt`, `par_jlin7.txt`, `chain_*.txt`, `run_jdiff8_10.txt`, `run_jlin*.txt` |
| Table 5: check of the S-box models, 4 keys | `python measure_best.py diff 6`, `python measure_best.py jdiff 3 5`, `python measure_best.py jlin 5 6` | `measure_best_diff.txt`, `measure_best_jdiff.txt` (r=3, 5), `measure_best_jlin.txt` (r=5, 6) |
| Table 8: differential and linear rows, 4 keys | `python measure_best.py jdiff 6 7 8 9 10`, `python measure_best.py jlin 7 8 9 10 11` | `measure_best_jdiff.txt`, `measure_best_jlin.txt`, `run_measure_*.txt` |
| Table 8: output-difference distribution (1 key; 20 rounds as reference) | `python ddist.py` | `ddist.txt` |
| Table 8: integral rows, 16 keys | `python measure.py integ 7 8 9` | `measure_integ.txt` (summary), `measure_integ_raw.txt` (keys and all XOR sums) |
| Table 6: impossible differential / zero-correlation searches on 20 rounds of the full cipher | `impossible.py`, `imp_one.py` | `impossible_r3.txt`, `impossible_r5.txt`, `impossible_r20.txt`, `run_impossible20.txt`, `run_imp_*.txt` |

Notes on Table 8.

- `measure_best_jlin.txt` contains two blocks for 11 rounds. The block with `w=14` is the minimum-weight
  characteristic of Table 7 and is the one used in Table 8.
- Integral rows: for every fixed input bit b, `measure integall` returns the XOR of E(x) over the 2^31
  inputs with x[b] = 0 (the other half has the same XOR, because E is a permutation). The statistic is
  the number of output bits whose XOR is zero for all 16 keys, and the table gives the maximum over
  the 32 choices of b. For a random permutation each output bit has a zero XOR with probability 1/2
  per key, so a bit is zero for all 16 keys with probability 2^-16.
- Random reference for the output-difference distribution: the number of pairs per output difference is
  Poisson with mean 1/2; over 2^32 output differences the maximum is 9 or 10 pairs (18 or 20 inputs)
  with probability about 0.97.
- The 20-round lower bounds of Table 7 (69 and 26) are sums of the minimum weights for 8+6+6 rounds
  (differential) and 12+8 rounds (linear). The minimum numbers of active S-boxes in Section 5.3 are
  2w/3 for the differential weights w of Table 7.

Example (generation only, no solver): `python trail/gen_model.py jdiff 16 8 4 2 3 > model.cvc`.

## Licence

MIT (see `LICENSE`). This covers the code in this repository only. The designers' code is not included and keeps its own terms.
