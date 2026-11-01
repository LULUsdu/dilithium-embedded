#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT -4186625 // 2^32 % Q
#define QINV 58728449 // q^(-1) mod 2^32

#define montgomery_reduce DILITHIUM_NAMESPACE(_montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

void ntt_montgomery_set_precompute_b(int32_t b);
int32_t ntt_montgomery_modmul_precompute_b(int32_t a);
int32_t ntt_montgomery_modmul_precompute_b_new0715_sched(int32_t a);
int32_t ntt_montgomery_modmul_direct_new0715_sched(int32_t a, int32_t b);
int32_t ntt_signed_montgomery_modmul_precompute_b(int32_t a);

/* Exact signed 64-bit product; Cortex-M0 signed-halfword implementation. */
int64_t mul_s32_s32_direct_cortex_m0(int32_t a, int32_t b);
/* Same result bits as set_precompute_b(b) + new0715_sched(a), no global state. */
int32_t ntt_montgomery_modmul_direct_signed16_cortex_m0(int32_t a, int32_t b);

#define reduce32 DILITHIUM_NAMESPACE(_reduce32)
int32_t reduce32(int32_t a);

#define caddq DILITHIUM_NAMESPACE(_caddq)
int32_t caddq(int32_t a);

#define freeze DILITHIUM_NAMESPACE(_freeze)
int32_t freeze(int32_t a);

#endif
