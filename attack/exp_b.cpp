// Experiment B: master-key recovery of full-round LLBC (official implementation), part by part.
//
// The key schedule is affine over GF(2): rk(K) = L(K) xor rk(0).  Part k of the cipher only
// sees the round-key bits in PART[k], i.e. L_k(K) = mask_k(L(K)), of rank 80 over the master key.
//   phase 1: guess L_0(K) (2^80 in the real attack) and filter with part 0 of 3 known pairs
//   phase 2..4: inside ker L_0 (then ker L_0 cap ker L_1, ...) guess the new directions of
//               part 1, 2, 3 (2^16 each) and filter with that part of the same 3 known pairs
// Reduced run: in phase 1 the attacker is told 80-d of the 80 coordinates of
// L_0(K*) in a fixed basis, and searches the other d (2^d candidates, Gray code: one round-key
// delta XOR per step).  Phases 2..4 run at full size.  A 4th known pair confirms the key.
//
// encrypt_rk / key_schedule_rk are our own implementation from the paper's equations that take
// the round keys as an argument; they are checked against the designers' functions and the
// official test vector before the attack.
// Build (WSL): g++ -O3 -march=native -fopenmp -w -o exp_b exp_b.cpp
// Run:         ./exp_b [d=32] [seed]
#include "../bench/llbc_lib.cpp"
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <random>
#include <vector>

static const u32 PART[4] = {0x000F000F, 0x00F000F0, 0x0F000F00, 0xF000F000};
static const int NW = ROUNDS * 4;  // round-key words
static const int NKP = 4;          // known pairs: 0..2 filter, 3 confirms

typedef u32 RK[ROUNDS][4];

// Independent implementation of the official code's encryption and key schedule, written from
// the equations in the paper (Sec. 2: Eq. (1)-(4)), with the round keys passed as an argument.
// S-boxes and constants are those of the LLBC paper (Tables III and IV).  self_check() compares
// them with the designers' functions (official test vector and 1000 random keys) before use.
static const uint8_t SB1[16] = {0x4,0x7,0x5,0x1,0xc,0xb,0xd,0x8,0xe,0xf,0x6,0x3,0xa,0x9,0x2,0x0};
static const uint8_t SB2[16] = {0xe,0xc,0x8,0x4,0xb,0xd,0x1,0x5,0xf,0x3,0xa,0x2,0x9,0x7,0x0,0x6};
static const uint64_t K_RC = 0x243F6A88ULL;   // round constant
static const uint64_t K_KS = 0x85A308D3ULL;   // key-schedule constant

static inline uint64_t rotl64(uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }
static inline uint16_t rotl16(uint16_t x, int r) { return (uint16_t)((x << r) | (x >> (16 - r))); }
static inline uint64_t sbox64(uint64_t x, const uint8_t *t) {
    uint64_t y = 0;
    for (int c = 0; c < 16; c++) y |= (uint64_t)t[(x >> (4 * c)) & 0xF] << (4 * c);
    return y;
}
static inline uint64_t mixword(uint64_t y1, uint64_t y2) {      // Eq. (3) with (b1, b2) = (32, 16)
    uint64_t t = y1 ^ y2, m0 = rotl64(t, 32) ^ y1, m1 = rotl64(m0, 16) ^ t;
    uint64_t m2 = rotl64(m1, 32) ^ m0, m3 = rotl64(m2, 16) ^ m1;
    return m2 ^ m3;
}
static inline uint64_t G(uint64_t x) {                          // 16-bit words (a, b, c, d), a on top
    uint16_t a = x >> 48, b = x >> 32, c = x >> 16, d = x;
    uint16_t a3 = rotl16(a, 3), b1 = rotl16(b, 1);
    return ((uint64_t)(uint16_t)(a3 ^ c) << 48) | ((uint64_t)(uint16_t)(b1 ^ d) << 32) | ((uint64_t)a3 << 16) | b1;
}

static void encrypt_rk(u32 state[4], const RK rk) {             // Eq. (1)-(2)
    uint64_t L = ((uint64_t)state[0] << 32) | state[1], R = ((uint64_t)state[2] << 32) | state[3];
    for (int j = 0; j < ROUNDS; j++) {
        uint64_t A = ((uint64_t)rk[j][0] << 32) | rk[j][1], B = ((uint64_t)rk[j][2] << 32) | rk[j][3];
        uint64_t f = mixword(sbox64(rotl64(L, 32), SB1), sbox64(rotl64(L, 16), SB2));
        uint64_t nL = f ^ R ^ B;
        R = L ^ A ^ (K_RC << 32);
        L = nL;
    }
    state[0] = L >> 32; state[1] = (u32)L; state[2] = R >> 32; state[3] = (u32)R;
}

static void key_schedule_rk(const u32 master[4], RK rk) {      // Eq. (4)
    uint64_t A = ((uint64_t)master[0] << 32) | master[1], B = ((uint64_t)master[2] << 32) | master[3];
    for (int j = 0; j < ROUNDS; j++) {
        B = G(G(B));
        A ^= B ^ K_RC;
        B ^= K_KS << 32;
        rk[j][0] = A >> 32; rk[j][1] = (u32)A; rk[j][2] = B >> 32; rk[j][3] = (u32)B;
    }
}

// ---- GF(2) linear algebra: key vectors (128 bits) and their masked round-key images ----
struct Vec { u32 img[NW]; u32 key[4]; int piv; };
static RK RK0;  // rk(0)

static void lin_image(const u32 key[4], int k, u32 img[NW]) {
    RK rk; key_schedule_rk(key, rk);
    for (int i = 0; i < NW; i++) img[i] = (rk[i / 4][i % 4] ^ RK0[i / 4][i % 4]) & PART[k];
}

static int first_bit(const u32 img[NW]) {
    for (int i = 0; i < NW; i++) if (img[i]) return 32 * i + __builtin_ctz(img[i]);
    return -1;
}

static void xorv(Vec &a, const Vec &b) {
    for (int i = 0; i < NW; i++) a.img[i] ^= b.img[i];
    for (int i = 0; i < 4; i++) a.key[i] ^= b.key[i];
}

static bool bit(const u32 img[NW], int p) { return (img[p / 32] >> (p % 32)) & 1; }

// reduce v against rows (insertion order); returns the set of rows used
static std::vector<int> reduce(Vec &v, const std::vector<Vec> &rows) {
    std::vector<int> used;
    for (size_t j = 0; j < rows.size(); j++)
        if (bit(v.img, rows[j].piv)) { xorv(v, rows[j]); used.push_back(j); }
    return used;
}

// split span(V) for part k: rows = vectors with independent L_k images, ker = basis of the kernel
static void split(const std::vector<Vec> &V, int k, std::vector<Vec> &rows, std::vector<Vec> &ker) {
    rows.clear(); ker.clear();
    for (const Vec &x : V) {
        Vec v; memcpy(v.key, x.key, 16); lin_image(v.key, k, v.img);
        reduce(v, rows);
        v.piv = first_bit(v.img);
        (v.piv < 0 ? ker : rows).push_back(v);
    }
}

static inline bool part_eq(const u32 a[4], const u32 b[4], int k) {
    return !(((a[0] ^ b[0]) | (a[1] ^ b[1]) | (a[2] ^ b[2]) | (a[3] ^ b[3])) & PART[k]);
}

static u32 P[NKP][4], C[NKP][4];

static bool pass_part(const RK rk, int k) {
    for (int i = 0; i < 3; i++) {
        u32 s[4] = {P[i][0], P[i][1], P[i][2], P[i][3]};
        encrypt_rk(s, rk);
        if (!part_eq(s, C[i], k)) return false;
    }
    return true;
}

static bool self_check(std::mt19937_64 &rng) {
    u32 key[4] = {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}, s[4] = {0, 0, 0, 0};
    const u32 tv[4] = {0x0C40B9A6, 0x69053EC4, 0xF1DD3EEA, 0x85366E1E};
    RK rk; key_schedule_rk(key, rk); encrypt_rk(s, rk);
    if (memcmp(s, tv, 16)) { printf("test vector mismatch\n"); return false; }
    for (int t = 0; t < 1000; t++) {
        u32 k1[4], k2[4], a[4], b[4];
        for (int i = 0; i < 4; i++) { k1[i] = k2[i] = (u32)rng(); a[i] = b[i] = (u32)rng(); }
        key_schedule(k1, ROUNDS); encrypt(a, ROUNDS);         // official (global roundKeys)
        key_schedule_rk(k2, rk); encrypt_rk(b, rk);           // own implementation
        if (memcmp(rk, roundKeys, sizeof rk) || memcmp(a, b, 16)) { printf("copy differs from official\n"); return false; }
    }
    return true;
}

int main(int argc, char **argv) {
    int d = argc > 1 ? atoi(argv[1]) : 32;
    unsigned long long seed = argc > 2 ? strtoull(argv[2], 0, 0) : std::random_device{}();
    std::mt19937_64 rng(seed);
    printf("experiment B: d=%d seed=%llu threads=%d\n", d, seed, omp_get_max_threads());
    if (d < 8 || d > 48) { printf("d must be in 8..48\n"); return 1; }
    if (!self_check(rng)) return 1;
    printf("self-check: test vector and 1000 random (key, plaintext) agree with the official code\n");

    u32 zero[4] = {0, 0, 0, 0};
    key_schedule_rk(zero, RK0);

    // secret key and known pairs (official code as the oracle)
    u32 Ks[4];
    for (int i = 0; i < 4; i++) Ks[i] = (u32)rng();
    { u32 t[4] = {Ks[0], Ks[1], Ks[2], Ks[3]}; key_schedule(t, ROUNDS); }
    for (int i = 0; i < NKP; i++) {
        for (int w = 0; w < 4; w++) P[i][w] = C[i][w] = (u32)rng();
        encrypt(C[i], ROUNDS);
    }
    printf("secret key: %08X %08X %08X %08X\n", Ks[0], Ks[1], Ks[2], Ks[3]);

    // ---- phase 1: part 0 ----
    std::vector<Vec> V(128), rows, ker;
    for (int b = 0; b < 128; b++) { memset(V[b].key, 0, 16); V[b].key[b / 32] = 1u << (b % 32); }
    split(V, 0, rows, ker);
    printf("phase 1: part 0 depends on %zu directions of the master key (kernel %zu)\n", rows.size(), ker.size());
    if ((int)rows.size() < d) return 1;

    Vec t; memcpy(t.key, Ks, 16); lin_image(t.key, 0, t.img);
    std::vector<int> used = reduce(t, rows);  // L_0(K*) = sum of rows[used].img
    unsigned long long a_true = 0;            // coordinates 0..d-1: unknown to the attacker
    u32 Kstart[4] = {0, 0, 0, 0};             // known part: coordinates d..79
    for (int j : used) {
        if (j < d) a_true |= 1ULL << j;
        else for (int i = 0; i < 4; i++) Kstart[i] ^= rows[j].key[i];
    }
    RK base; key_schedule_rk(Kstart, base);

    std::vector<unsigned long long> surv;
    unsigned long long n_enc = 0;
    int cb = d - 8;  // 256 chunks of 2^cb Gray-code steps
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:n_enc)
    for (long long c = 0; c < 256; c++) {
        unsigned long long i0 = (unsigned long long)c << cb, i1 = i0 + (1ULL << cb);
        unsigned long long g = i0 ^ (i0 >> 1);
        RK rk; memcpy(rk, base, sizeof rk);
        u32 *w = &rk[0][0];
        for (int j = 0; j < d; j++) if ((g >> j) & 1) for (int x = 0; x < NW; x++) w[x] ^= rows[j].img[x];
        for (unsigned long long i = i0; i < i1; i++) {
            u32 s[4] = {P[0][0], P[0][1], P[0][2], P[0][3]};
            encrypt_rk(s, rk);
            n_enc++;
            if (part_eq(s, C[0], 0) && pass_part(rk, 0)) {
                n_enc += 2;
                #pragma omp critical
                surv.push_back(i ^ (i >> 1));
            }
            if (i + 1 < i1) {
                int j = __builtin_ctzll(i + 1);
                for (int x = 0; x < NW; x++) w[x] ^= rows[j].img[x];
            }
        }
        if (c % 32 == 31) {
            #pragma omp critical
            { printf("  ... %lld/256 chunks done, %.0f s\n", c + 1, omp_get_wtime() - t0); fflush(stdout); }
        }
    }
    double dt = omp_get_wtime() - t0;
    printf("phase 1: %llu encryptions (2^%.2f), %.1f s, %zu survivor(s)\n", n_enc, log2((double)n_enc), dt, surv.size());
    for (auto a : surv) printf("  candidate %s\n", a == a_true ? "= true coordinates" : "!= true coordinates");

    // ---- phases 2..4: parts 1, 2, 3 inside the remaining kernel ----
    int found = 0;
    for (auto a : surv) {
        u32 K[4] = {Kstart[0], Kstart[1], Kstart[2], Kstart[3]};
        for (int j = 0; j < d; j++) if ((a >> j) & 1) for (int i = 0; i < 4; i++) K[i] ^= rows[j].key[i];
        std::vector<Vec> Vk = ker;
        bool alive = true;
        for (int k = 1; k <= 3 && alive; k++) {
            std::vector<Vec> rk_rows, rk_ker;
            split(Vk, k, rk_rows, rk_ker);
            int r = rk_rows.size();
            std::vector<unsigned long long> sv;
            for (unsigned long long b = 0; b < (1ULL << r); b++) {
                u32 Kc[4] = {K[0], K[1], K[2], K[3]};
                for (int j = 0; j < r; j++) if ((b >> j) & 1) for (int i = 0; i < 4; i++) Kc[i] ^= rk_rows[j].key[i];
                RK rkc; key_schedule_rk(Kc, rkc);
                if (pass_part(rkc, k)) sv.push_back(b);
            }
            printf("phase %d: part %d adds %d new directions (kernel left %zu), 2^%d candidates, %zu survivor(s)\n",
                   k + 1, k, r, rk_ker.size(), r, sv.size());
            if (sv.size() != 1) { alive = false; break; }  // a reduced run expects exactly one
            for (int j = 0; j < r; j++) if ((sv[0] >> j) & 1) for (int i = 0; i < 4; i++) K[i] ^= rk_rows[j].key[i];
            Vk = rk_ker;
        }
        if (!alive) continue;
        RK rkf; key_schedule_rk(K, rkf);
        u32 s[4] = {P[3][0], P[3][1], P[3][2], P[3][3]};
        encrypt_rk(s, rkf);
        bool conf = !memcmp(s, C[3], 16), same = !memcmp(K, Ks, 16);
        printf("recovered key: %08X %08X %08X %08X  (kernel left %zu)\n", K[0], K[1], K[2], K[3], Vk.size());
        printf("  4th known pair encrypts correctly: %s;  equals the secret key: %s\n", conf ? "yes" : "no", same ? "yes" : "no");
        found += conf && same;
    }
    return found == 1 ? 0 : 1;
}
