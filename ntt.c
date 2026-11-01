#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"

const int32_t dilithium_zetas[N] = {
         0,    25847, -2608894,  -518909,   237124,  -777960,  -876248,   466468,
   1826347,  2353451,  -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
   2725464,  1024112, -1079900,  3585928,  -549488, -1119584,  2619752, -2108549,
  -2118186, -3859737, -1399561, -3277672,  1757237,   -19422,  4010497,   280005,
   2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
  -3861115, -3043716,  3574422, -2867647,  3539968,  -300467,  2348700,  -539299,
  -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
    811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
  -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
  -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
   3412210,  -983419,  2147896,  2715295, -2967645, -3693493,  -411027, -2477047,
   -671102, -1228525,   -22981, -1308169,  -381987,  1349076,  1852771, -1430430,
  -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
  -3724342,    -8578,  1653064, -3249728,  2389356,  -210977,   759969, -1316856,
    189548, -3553272,  3159746, -1851402, -2409325,  -177440,  1315589,  1341330,
   1285669, -1584928,  -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
   2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
    266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
    900702,  1859098,   909542,   819034,   495491, -1613174,   -43260,  -522500,
   -655327, -3122442,  2031748,  3207046, -3556995,  -525098,  -768622, -3595838,
    342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
   2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
  -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
  -1333058,  1237275, -3318210, -1430225,  -451100,  1312455,  3306115, -1962642,
  -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
   -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
  -2013608,  2432395,  2454455,  -164721,  1957272,  3369112,   185531, -1207385,
  -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
  -3038916,  3523897,  3866901,   269760,  2213111,  -975884,  1717735,   472078,
   -426683,  1723600, -1803090,  1910376, -1667432, -1104333,  -260646, -3833893,
  -2939036, -2235985,  -420899, -2286327,   183443,  -976891,  1612842, -3545687,
   -554416,  3919660,   -48306, -1362209,  3937738,  1400424,  -846154,  1976782
};

/*
 * Twiddle factors reordered for the forward 2+2+2+2 merged NTT.
 * Entry 0 remains unused.  For each merged parent group the following
 * three entries are consecutive: parent zeta, left-child zeta,
 * right-child zeta.  The original dilithium_zetas table is retained for
 * the inverse NTT and the C reference implementations.
 */
const int32_t dilithium_merge_zetas[N] = {
          0,     25847,  -2608894,   -518909,    237124,   1826347,   2353451,   -777960,
    -359251,  -2091905,   -876248,   3119733,  -2884855,    466468,   3111497,   2680103,
    2725464,   2706023,     95776,   1024112,   3077325,   3530437,  -1079900,  -1661693,
   -3592148,   3585928,  -2537516,   3915439,   -549488,  -3861115,  -3043716,  -1119584,
    3574422,  -2867647,   2619752,   3539968,   -300467,  -2108549,   2348700,   -539299,
   -2118186,  -1699267,  -1643818,  -3859737,   3505694,  -3821735,  -1399561,   3507263,
   -2140649,  -3277672,  -1600420,   3699596,   1757237,    811944,    531354,    -19422,
     954230,   3881043,   4010497,   3900724,  -2556880,    280005,   2071892,  -2797779,
   -3930395,   2091667,   3407706,  -1528703,   2316500,   3817976,  -3677745,  -3342478,
    2244091,  -3041255,  -2446433,  -3562462,  -1452451,    266997,   2434439,   3475950,
   -1235728,   3513181,   2176455,  -3520352,  -3759364,  -1585221,  -1197226,  -3193378,
   -1257611,    900702,   1859098,   1939314,    909542,    819034,  -4083598,    495491,
   -1613174,  -1000202,    -43260,   -522500,  -3190144,   -655327,  -3122442,  -3157330,
    2031748,   3207046,  -3632928,  -3556995,   -525098,    126922,   -768622,  -3595838,
    3412210,    342297,    286988,   -983419,  -2437823,   4108315,   2147896,   3437287,
   -3342277,   2715295,   1735879,    203044,  -2967645,   2842341,   2691481,  -3693493,
   -2590150,   1265009,   -411027,   4055324,   1247620,  -2477047,   2486353,   1595974,
    -671102,  -3767016,   1250494,  -1228525,   2635921,  -3548272,    -22981,  -2994039,
    1869119,  -1308169,   1903435,  -1050970,   -381987,  -1333058,   1237275,   1349076,
   -3318210,  -1430225,   1852771,   -451100,   1312455,  -1430430,   3306115,  -1962642,
   -3343383,  -1279661,   1917081,    264944,  -2546312,  -1374803,    508951,   1500165,
     777191,   3097992,   2235880,   3406031,     44288,   -542412,  -2831860,  -1100098,
   -1671176,  -1846953,    904516,  -2584293,  -3724270,   3958618,    594136,  -3776993,
   -3724342,  -2013608,   2432395,     -8578,   2454455,   -164721,   1653064,   1957272,
    3369112,  -3249728,    185531,  -1207385,   2389356,  -3183426,    162844,   -210977,
    1616392,   3014001,    759969,    810149,   1652634,  -1316856,  -3694233,  -1799107,
     189548,  -3038916,   3523897,  -3553272,   3866901,    269760,   3159746,   2213111,
    -975884,  -1851402,   1717735,    472078,  -2409325,   -426683,   1723600,   -177440,
   -1803090,   1910376,   1315589,  -1667432,  -1104333,   1341330,   -260646,  -3833893,
    1285669,  -2939036,  -2235985,  -1584928,   -420899,  -2286327,   -812732,    183443,
    -976891,  -1439742,   1612842,  -3545687,  -3019102,   -554416,   3919660,  -3881060,
     -48306,  -1362209,  -3628969,   3937738,   1400424,   3839961,   -846154,   1976782
};

/*
 * Twiddle factors for the development-only 2+2+2+2 merged inverse NTT.
 *
 * Every three consecutive entries are already-negated Montgomery-domain
 * twiddles in [left child, right child, parent] order.  The four merged
 * stages contain 64, 16, 4 and 1 triples respectively.  This file is
 * generated from dilithium_zetas by tools/gen_invntt_merge_zetas.py; do not
 * edit the constants by hand.
 */
const int32_t dilithium_inv_merge_zetas[255] = {
#include "invntt_merge_zetas.inc"
};

/*************************************************
* Name:        ntt
*
* Description: Forward NTT, in-place. No modular reduction is performed after
*              additions or subtractions. Output vector is in bitreversed order.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
void ntt_c_ref(int32_t a[N]) {
  unsigned int len, start, j, k;
  int32_t zeta, t;

  k = 0;
  for(len = 128; len > 0; len >>= 1) {
    for(start = 0; start < N; start = j + len) {
      zeta = dilithium_zetas[++k];
      ntt_montgomery_set_precompute_b(zeta);
      for(j = start; j < start + len; ++j) {
        t = ntt_montgomery_modmul_precompute_b_new0715_sched(a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont_c_ref
*
* Description: Inverse NTT and multiplication by Montgomery factor 2^32.
*              In-place. No modular reductions after additions or
*              subtractions; input coefficients need to be smaller than
*              Q in absolute value. Output coefficient are smaller than Q in
*              absolute value.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
void invntt_tomont_c_ref(int32_t a[N]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
  const int32_t f = 41978; // mont^2/256

  k = 256;
  for(len = 1; len < N; len <<= 1) {
    for(start = 0; start < N; start = j + len) {
      zeta = -dilithium_zetas[--k];
      ntt_montgomery_set_precompute_b(zeta);
      for(j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = ntt_montgomery_modmul_precompute_b_new0715_sched(a[j + len]);
      }
    }
  }

  ntt_montgomery_set_precompute_b(f);
  for(j = 0; j < N; ++j) {
    a[j] = ntt_montgomery_modmul_precompute_b_new0715_sched(a[j]);
  }
}

/*
 * Development oracle for a two-layer merged inverse NTT.
 *
 * This deliberately retains the same Montgomery multiplication primitive and
 * the separate final scale as invntt_tomont_c_ref.  It is not wired to the
 * production invntt_tomont/invntt_asm entry points.
 */
#ifdef INVNTT_2LAYER_MODEL
static void invntt_tomont_2layer_c_core(
    int32_t a[N]
#ifdef INVNTT_2LAYER_DEBUG
    , int32_t stages[4][N]
#endif
) {
  unsigned int block, j, len, zi;
#ifdef INVNTT_2LAYER_DEBUG
  unsigned int stage;
#endif
  int32_t x0, x1, x2, x3;
  int32_t t0, t1, t2, t3;
  int32_t zeta_left, zeta_right, zeta_parent;
  const int32_t f = 41978; /* mont^2/256 */

  zi = 0;
#ifdef INVNTT_2LAYER_DEBUG
  stage = 0;
#endif
  for(len = 1; len < N; len <<= 2) {
    for(block = 0; block < N; block += 4*len) {
      zeta_left   = dilithium_inv_merge_zetas[zi++];
      zeta_right  = dilithium_inv_merge_zetas[zi++];
      zeta_parent = dilithium_inv_merge_zetas[zi++];

      for(j = block; j < block + len; ++j) {
        x0 = a[j];
        x1 = a[j + len];
        x2 = a[j + 2*len];
        x3 = a[j + 3*len];

        t0 = x0 + x1;
        ntt_montgomery_set_precompute_b(zeta_left);
        t1 = ntt_montgomery_modmul_precompute_b_new0715_sched(x0 - x1);

        t2 = x2 + x3;
        ntt_montgomery_set_precompute_b(zeta_right);
        t3 = ntt_montgomery_modmul_precompute_b_new0715_sched(x2 - x3);

        a[j] = t0 + t2;
        ntt_montgomery_set_precompute_b(zeta_parent);
        a[j + 2*len] =
            ntt_montgomery_modmul_precompute_b_new0715_sched(t0 - t2);

        a[j + len] = t1 + t3;
        a[j + 3*len] =
            ntt_montgomery_modmul_precompute_b_new0715_sched(t1 - t3);
      }
    }

#ifdef INVNTT_2LAYER_DEBUG
    if(stages != 0) {
      for(j = 0; j < N; ++j)
        stages[stage][j] = a[j];
    }
    ++stage;
#endif
  }

  /*
   * Keep scaling independent so this oracle tests only two-layer butterfly
   * merging, not a second scale-fusion transformation.
   */
  ntt_montgomery_set_precompute_b(f);
  for(j = 0; j < N; ++j)
    a[j] = ntt_montgomery_modmul_precompute_b_new0715_sched(a[j]);

  /* 85 triples = 255 entries; retained as an auditable loop invariant. */
  (void)zi;
}

void invntt_tomont_2layer_c_ref(int32_t a[N]) {
  invntt_tomont_2layer_c_core(
      a
#ifdef INVNTT_2LAYER_DEBUG
      , 0
#endif
  );
}

#ifdef INVNTT_2LAYER_DEBUG
void invntt_tomont_2layer_c_ref_debug(int32_t a[N],
                                     int32_t stages[4][N]) {
  invntt_tomont_2layer_c_core(a, stages);
}
#endif
#endif
