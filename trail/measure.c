// Exhaustive measurements on one 32-bit part of the reference implementation of LLBC
// (mini_llbc.py: halves of 16 bits, nibble j = cell k+4j, S1 on L<<<8, S2 on L<<<4,
//  MixWord with rotations (8,4); L' = MW xor R xor B_j, R' = L xor A_j xor rc).
//
// Round keys come on stdin: "r rc A_0 B_0 A_1 B_1 ..." (hex), computed in measure.py from a
// master key with the official key schedule (part k bits).
//
//   ./measure test                         read keys, then lines "L R" -> print E(L,R)
//   ./measure diff  dL dR oL oR            #{x in 2^32 : E(x) xor E(x xor din) = dout}
//   ./measure lin   gL gR hL hR            sum_x (-1)^(<g,x> xor <h,E(x)>)
//   ./measure integ bit val                xor of E(x) over the 2^31 x with x[bit] = val
//   ./measure integall                     32 lines: for b = 0..31, xor of E(x) over the x with x[b] = 0
//   ./measure ddist dL dR                  distribution of E(x) xor E(x xor din) over all x
// All numeric arguments are hexadecimal (bit 10 is written "a").
// x = (L << 16) | R.  Build: gcc -O3 -march=native -fopenmp -o measure measure.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

static const int S1[16] = {0x4,0x7,0x5,0x1,0xc,0xb,0xd,0x8,0xe,0xf,0x6,0x3,0xa,0x9,0x2,0x0};
static const int S2[16] = {0xe,0xc,0x8,0x4,0xb,0xd,0x1,0x5,0xf,0x3,0xa,0x2,0x9,0x7,0x0,0x6};
static uint16_t T1[65536], T2[65536];
static int R;
static uint16_t RC, KA[64], KB[64];

static inline uint16_t rol(uint16_t x, int r) { return (uint16_t)((x << r) | (x >> (16 - r))); }

static inline uint32_t enc(uint32_t x) {
    uint16_t L = x >> 16, Rr = x & 0xFFFF;
    for (int j = 0; j < R; j++) {
        uint16_t s1 = T1[rol(L, 8)], s2 = T2[rol(L, 4)];
        uint16_t t = s1 ^ s2;
        uint16_t m0 = rol(t, 8) ^ s1;
        uint16_t m1 = rol(m0, 4) ^ t;
        uint16_t m2 = rol(m1, 8) ^ m0;
        uint16_t m3 = rol(m2, 4) ^ m1;
        uint16_t nL = m2 ^ m3 ^ Rr ^ KB[j];
        Rr = L ^ KA[j] ^ RC;
        L = nL;
    }
    return ((uint32_t)L << 16) | Rr;
}

int main(int argc, char **argv) {
    for (int x = 0; x < 65536; x++) {
        int a = 0, b = 0;
        for (int c = 0; c < 4; c++) {
            a |= S1[(x >> (4 * c)) & 0xF] << (4 * c);
            b |= S2[(x >> (4 * c)) & 0xF] << (4 * c);
        }
        T1[x] = a; T2[x] = b;
    }
    unsigned v;
    if (scanf("%d %x", &R, &v) != 2) return 1;
    RC = v;
    for (int j = 0; j < R; j++) { unsigned a, b; if (scanf("%x %x", &a, &b) != 2) return 1; KA[j] = a; KB[j] = b; }
    const char *mode = argv[1];
    if (!strcmp(mode, "test")) {
        unsigned l, r;
        while (scanf("%x %x", &l, &r) == 2) printf("%08x\n", enc((l << 16) | r));
        return 0;
    }
    unsigned p[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4 && i + 2 < argc; i++) p[i] = strtoul(argv[i + 2], 0, 16);
    double t0 = omp_get_wtime();
    if (!strcmp(mode, "diff")) {
        uint32_t din = (p[0] << 16) | p[1], dout = (p[2] << 16) | p[3];
        unsigned long long cnt = 0;
        #pragma omp parallel for schedule(static) reduction(+:cnt)
        for (long long hi = 0; hi < 65536; hi++)
            for (uint32_t lo = 0; lo < 65536; lo++) {
                uint32_t x = ((uint32_t)hi << 16) | lo;
                cnt += (enc(x) ^ enc(x ^ din)) == dout;
            }
        printf("%llu\n", cnt);
    } else if (!strcmp(mode, "lin")) {
        uint32_t g = (p[0] << 16) | p[1], h = (p[2] << 16) | p[3];
        long long s = 0;
        #pragma omp parallel for schedule(static) reduction(+:s)
        for (long long hi = 0; hi < 65536; hi++)
            for (uint32_t lo = 0; lo < 65536; lo++) {
                uint32_t x = ((uint32_t)hi << 16) | lo;
                s += (__builtin_parity(x & g) ^ __builtin_parity(enc(x) & h)) ? -1 : 1;
            }
        printf("%lld\n", s);
    } else if (!strcmp(mode, "ddist")) {
        // full output-difference distribution for one input difference (ordered x, 2^32 inputs)
        uint32_t din = (p[0] << 16) | p[1];
        uint16_t *h = (uint16_t *)calloc(1ull << 32, sizeof(uint16_t));
        if (!h) { fprintf(stderr, "no memory\n"); return 1; }
        #pragma omp parallel for schedule(static)
        for (long long hi = 0; hi < 65536; hi++)
            for (uint32_t lo = 0; lo < 65536; lo++) {
                uint32_t x = ((uint32_t)hi << 16) | lo;
                uint32_t d = enc(x) ^ enc(x ^ din);
                __atomic_fetch_add(&h[d], 1, __ATOMIC_RELAXED);
            }
        unsigned long long tail[8] = {0};  // #outputs with count >= 8, 16, 24, 32, 48, 64, 128, 256
        const int th[8] = {8, 16, 24, 32, 48, 64, 128, 256};
        uint32_t top[5] = {0}; unsigned topc[5] = {0};
        for (unsigned long long d = 0; d < (1ull << 32); d++) {
            unsigned c = h[d];
            for (int t = 0; t < 8; t++) tail[t] += c >= th[t];
            if (c > topc[4]) {
                int k = 4;
                while (k > 0 && c > topc[k - 1]) { topc[k] = topc[k - 1]; top[k] = top[k - 1]; k--; }
                topc[k] = c; top[k] = (uint32_t)d;
            }
        }
        printf("top:");
        for (int k = 0; k < 5; k++) printf(" %08x:%u", top[k], topc[k]);
        printf("\ntail(>=8,16,24,32,48,64,128,256):");
        for (int t = 0; t < 8; t++) printf(" %llu", tail[t]);
        printf("\n");
        free(h);
    } else if (!strcmp(mode, "integall")) {
        // one pass: for every bit b, xor of E(x) over the x with x[b] = 0 (2^31 texts each)
        uint32_t acc[32] = {0};
        #pragma omp parallel
        {
            uint32_t loc[32] = {0};
            #pragma omp for schedule(static)
            for (long long hi = 0; hi < 65536; hi++)
                for (uint32_t lo = 0; lo < 65536; lo++) {
                    uint32_t x = ((uint32_t)hi << 16) | lo, e = enc(x), z = ~x;
                    for (int b = 0; b < 32; b++) loc[b] ^= e & -((z >> b) & 1);
                }
            #pragma omp critical
            for (int b = 0; b < 32; b++) acc[b] ^= loc[b];
        }
        for (int b = 0; b < 32; b++) printf("%08x\n", acc[b]);
    } else if (!strcmp(mode, "integ")) {
        int bit = p[0] & 31, val = p[1] & 1;
        uint32_t acc = 0;
        #pragma omp parallel for schedule(static) reduction(^:acc)
        for (long long hi = 0; hi < 32768; hi++)
            for (uint32_t lo = 0; lo < 65536; lo++) {
                uint32_t y = ((uint32_t)hi << 16) | lo;           // 31 free bits
                uint32_t low = y & ((1u << bit) - 1), high = bit == 31 ? 0 : (y >> bit) << (bit + 1);
                uint32_t x = high | ((uint32_t)val << bit) | low;
                acc ^= enc(x);
            }
        printf("%08x\n", acc);
    }
    fprintf(stderr, "%.1f s\n", omp_get_wtime() - t0);
    return 0;
}
