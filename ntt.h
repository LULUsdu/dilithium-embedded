#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define ntt DILITHIUM_NAMESPACE(_ntt)
void ntt(int32_t a[N]);

/*
 * Forward-NTT entry points used by the dedicated correctness self-test.
 * ntt_asm_abi_probe returns zero on success and a nonzero ABI failure mask.
 */
void ntt_asm(int32_t a[N]);
void ntt_c_ref(int32_t a[N]);
int ntt_asm_abi_probe(int32_t a[N]);

#define invntt_tomont DILITHIUM_NAMESPACE(_invntt_tomont)
void invntt_tomont(int32_t a[N]);
void invntt_asm(int32_t a[N]);
void invntt_tomont_c_ref(int32_t a[N]);
#ifdef INVNTT_2LAYER_MODEL
void invntt_tomont_2layer_c_ref(int32_t a[N]);
#ifdef INVNTT_2LAYER_DEBUG
void invntt_tomont_2layer_c_ref_debug(int32_t a[N],
                                     int32_t stages[4][N]);
#endif
#endif
int invntt_asm_abi_probe(int32_t a[N]);

#endif
