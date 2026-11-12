// Experiment A: decrypt arbitrary ciphertexts of full-round LLBC (official implementation)
// with 2^32 chosen plaintexts, using that the cipher splits into four independent 32-bit ciphers.
//
// Part k of a 32-bit word = bits b with (b mod 16) / 4 == k (mask PART[k], 8 bits per word),
// so part k of the 128-bit state is 32 bits.  Query the plaintexts P(v) whose four parts all
// equal v, for v = 0 .. 2^32-1.  Part k of E(P(v)) is E_k(v), a permutation of v.  For each
// target ciphertext C, part k of C equals E_k(v) for exactly one v, and that v is part k of the
// plaintext.  The attack only sees C; the plaintexts are kept to check the answer.
//
// The oracle is the official encrypt() with roundKeys from the official key_schedule().
// Build (WSL): g++ -O3 -march=native -fopenmp -w -o exp_a exp_a.cpp
// Run:         ./exp_a [targets=1000] [seed]
#include "../bench/llbc_lib.cpp"
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <random>
#include <vector>
#include <algorithm>

static const u32 PART[4] = {0x000F000F, 0x00F000F0, 0x0F000F00, 0xF000F000};

static inline u32 get_part(const u32 s[4], int k) {
    u32 v = 0;
    for (int w = 0; w < 4; w++) v |= _pext_u32(s[w], PART[k]) << (8 * w);
    return v;
}

static inline void set_part(u32 s[4], int k, u32 v) {
    for (int w = 0; w < 4; w++) s[w] = (s[w] & ~PART[k]) | _pdep_u32((v >> (8 * w)) & 0xFF, PART[k]);
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 1000;
    unsigned long long seed = argc > 2 ? strtoull(argv[2], 0, 0) : std::random_device{}();
    std::mt19937_64 rng(seed);
    printf("experiment A: targets=%d seed=%llu threads=%d\n", n, seed, omp_get_max_threads());

    // secret key and targets (the attack below uses only the ciphertexts)
    u32 key[4];
    for (int i = 0; i < 4; i++) key[i] = (u32)rng();
    printf("secret key: %08X %08X %08X %08X\n", key[0], key[1], key[2], key[3]);
    u32 ks[4] = {key[0], key[1], key[2], key[3]};
    key_schedule(ks, ROUNDS);

    // check the harness on the official test vector would overwrite roundKeys, so it is done
    // in exp_b; here only check encrypt/decrypt round trip of the official code
    std::vector<u32> P(4 * n), C(4 * n);
    for (int i = 0; i < n; i++) {
        for (int w = 0; w < 4; w++) P[4 * i + w] = (u32)rng();
        memcpy(&C[4 * i], &P[4 * i], 16);
        encrypt(&C[4 * i], ROUNDS);
        u32 d[4]; memcpy(d, &C[4 * i], 16); decrypt(d, ROUNDS);
        if (memcmp(d, &P[4 * i], 16)) { printf("official decrypt does not invert encrypt\n"); return 1; }
    }

    // per part: sorted target values, and a 2^20-bit prefilter on the top 20 bits
    std::vector<u32> tv[4];
    std::vector<unsigned long long> bm[4];
    for (int k = 0; k < 4; k++) {
        bm[k].assign(1 << 14, 0);
        for (int i = 0; i < n; i++) {
            u32 t = get_part(&C[4 * i], k);
            tv[k].push_back(t);
            bm[k][t >> 26] |= 1ULL << ((t >> 20) & 63);
        }
        std::sort(tv[k].begin(), tv[k].end());
        tv[k].erase(std::unique(tv[k].begin(), tv[k].end()), tv[k].end());
    }

    // the 2^32 chosen-plaintext queries
    std::vector<std::pair<u32, u32>> hit[4];  // (E_k(v), v)
    unsigned long long n_query = 0;
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:n_query)
    for (long long hi = 0; hi < 256; hi++) {
        std::vector<std::pair<u32, u32>> loc[4];
        for (u32 lo = 0; lo < (1u << 24); lo++) {
            u32 v = ((u32)hi << 24) | lo;
            u32 s[4] = {0, 0, 0, 0};
            for (int k = 0; k < 4; k++) set_part(s, k, v);
            encrypt(s, ROUNDS);
            n_query++;
            for (int k = 0; k < 4; k++) {
                u32 e = get_part(s, k);
                if (!((bm[k][e >> 26] >> ((e >> 20) & 63)) & 1)) continue;
                if (std::binary_search(tv[k].begin(), tv[k].end(), e)) loc[k].push_back({e, v});
            }
        }
        #pragma omp critical
        {
            for (int k = 0; k < 4; k++) hit[k].insert(hit[k].end(), loc[k].begin(), loc[k].end());
            if (hi % 32 == 31) { printf("  ... %lld/256 blocks, %.0f s\n", hi + 1, omp_get_wtime() - t0); fflush(stdout); }
        }
    }
    double dt = omp_get_wtime() - t0;
    printf("queries: %llu (= 2^%.2f), time %.1f s\n", n_query, log2((double)n_query), dt);

    // E_k is a permutation, so every distinct target value must be hit exactly once
    for (int k = 0; k < 4; k++) {
        std::sort(hit[k].begin(), hit[k].end());
        bool dup = std::adjacent_find(hit[k].begin(), hit[k].end(),
            [](auto &a, auto &b) { return a.first == b.first; }) != hit[k].end();
        printf("part %d: %zu distinct target values, %zu hits%s\n", k, tv[k].size(), hit[k].size(),
               dup ? "  (DUPLICATE: not a permutation?)" : "");
        if (hit[k].size() != tv[k].size() || dup) return 1;
    }

    // reconstruct each plaintext from the ciphertext alone, then check
    int ok_p = 0, ok_c = 0;
    for (int i = 0; i < n; i++) {
        u32 r[4] = {0, 0, 0, 0};
        for (int k = 0; k < 4; k++) {
            u32 t = get_part(&C[4 * i], k);
            auto it = std::lower_bound(hit[k].begin(), hit[k].end(), std::make_pair(t, 0u));
            set_part(r, k, it->second);
        }
        ok_p += !memcmp(r, &P[4 * i], 16);
        encrypt(r, ROUNDS);
        ok_c += !memcmp(r, &C[4 * i], 16);
    }
    printf("recovered plaintext == original: %d / %d\n", ok_p, n);
    printf("re-encryption == target:          %d / %d\n", ok_c, n);
    printf("example: C = %08X %08X %08X %08X -> P = %08X %08X %08X %08X\n",
           C[0], C[1], C[2], C[3], P[0], P[1], P[2], P[3]);
    return ok_p == n && ok_c == n ? 0 : 1;
}
