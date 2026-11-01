#include <stdint.h>
#include <string.h>

#include "poly_mul_12289.h"

#define Q12289 12289

void NTT_forward12289new2(const uint16_t *src, int32_t *workspace,
                          uint16_t *dst);
void PairMul12289new2(const uint16_t *a_ntt, const uint16_t *b_ntt,
                      uint16_t *acc_ntt);
void NTT_backward12289new3(int32_t *workspace, uint16_t *srcdst);

/*
 * The assembly routines use packed 16-bit inputs/outputs and a separate
 * 32-bit workspace.  Static storage avoids placing 3584 bytes on the small
 * Cortex-M0 stack.  Consequently poly_mul_12289 is not reentrant and must
 * not be called concurrently or from an interrupt that can preempt it.
 */
static uint16_t input_a[N];
static uint16_t input_b[N];
static uint16_t ntt_a[N];
static uint16_t ntt_b[N];
static uint16_t ntt_product[N];
static int32_t ntt_workspace[N];
static uint16_t prepared_s1_12289[L][N];
static uint16_t prepared_s2_12289[K][N];

volatile unsigned int poly_mul_12289_forward_calls;

static uint16_t encode_mod_12289(int32_t value)
{
  int32_t reduced = value % Q12289;

  if(reduced < 0)
    reduced += Q12289;
  return (uint16_t)reduced;
}

static int32_t decode_centered_12289(uint16_t value)
{
  int32_t reduced = (int32_t)(int16_t)value % Q12289;

  if(reduced < 0)
    reduced += Q12289;
  if(reduced > Q12289/2)
    reduced -= Q12289;
  return reduced;
}

static void encode_poly_12289(uint16_t dst[N], const poly *src)
{
  unsigned int i;

  for(i = 0; i < N; ++i)
    dst[i] = encode_mod_12289(src->coeffs[i]);
}

static void decode_poly_12289(poly *dst, const uint16_t src[N])
{
  unsigned int i;

  for(i = 0; i < N; ++i)
    dst->coeffs[i] = decode_centered_12289(src[i]);
}

static void mul_with_precomputed_a_12289(poly *c, const poly *b)
{
  // encode_poly_12289(input_b, b);
  NTT_forward12289new2(input_b, ntt_workspace, ntt_b);
  ++poly_mul_12289_forward_calls;

  /* PairMul12289new2 adds into its third argument. */
  memset(ntt_product, 0, sizeof(ntt_product));
  PairMul12289new2(ntt_a, ntt_b, ntt_product);
  NTT_backward12289new3(ntt_workspace, ntt_product);
  // decode_poly_12289(c, ntt_product);
}

static void mul_prepared_12289(poly *c, const uint16_t b_ntt[N])
{
  /* PairMul12289new2 adds into its third argument. */
  memset(ntt_product, 0, sizeof(ntt_product));
  PairMul12289new2(ntt_a, b_ntt, ntt_product);
  NTT_backward12289new3(ntt_workspace, ntt_product);
  decode_poly_12289(c, ntt_product);
}

void poly_mul_12289(poly *c, const poly *a, const poly *b)
{
//  encode_poly_12289(input_a, a);
  NTT_forward12289new2(input_a, ntt_workspace, ntt_a);
  ++poly_mul_12289_forward_calls;
  mul_with_precomputed_a_12289(c, b);
}

void polyvecl_mul_challenge_12289(polyvecl *r,
                                  const poly *challenge,
                                  const polyvecl *s1)
{
  unsigned int i;

  // encode_poly_12289(input_a, challenge);
  NTT_forward12289new2(input_a, ntt_workspace, ntt_a);
  ++poly_mul_12289_forward_calls;

  for(i = 0; i < L; ++i)
    mul_with_precomputed_a_12289(&r->vec[i], &s1->vec[i]);
}

void polyveck_mul_challenge_12289(polyveck *r,
                                  const poly *challenge,
                                  const polyveck *s2)
{
  unsigned int i;

  // encode_poly_12289(input_a, challenge);
  NTT_forward12289new2(input_a, ntt_workspace, ntt_a);
  ++poly_mul_12289_forward_calls;

  for(i = 0; i < K; ++i)
    mul_with_precomputed_a_12289(&r->vec[i], &s2->vec[i]);
}

void poly_mul_12289_prepare_secrets(const polyvecl *s1,
                                    const polyveck *s2)
{
  unsigned int i;

  for(i = 0; i < L; ++i) {
    encode_poly_12289(input_b, &s1->vec[i]);
    NTT_forward12289new2(input_b, ntt_workspace, prepared_s1_12289[i]);
    ++poly_mul_12289_forward_calls;
  }

  for(i = 0; i < K; ++i) {
    encode_poly_12289(input_b, &s2->vec[i]);
    NTT_forward12289new2(input_b, ntt_workspace, prepared_s2_12289[i]);
    ++poly_mul_12289_forward_calls;
  }
}

void poly_mul_12289_prepare_challenge(const poly *challenge)
{
  encode_poly_12289(input_a, challenge);
  NTT_forward12289new2(input_a, ntt_workspace, ntt_a);
  ++poly_mul_12289_forward_calls;
}

void polyvecl_mul_prepared_12289(polyvecl *r)
{
  unsigned int i;

  for(i = 0; i < L; ++i)
    mul_prepared_12289(&r->vec[i], prepared_s1_12289[i]);
}

void polyveck_mul_prepared_12289(polyveck *r)
{
  unsigned int i;

  for(i = 0; i < K; ++i)
    mul_prepared_12289(&r->vec[i], prepared_s2_12289[i]);
}

void poly_mul_12289_clear_prepared_secrets(void)
{
  volatile uint16_t *p;
  volatile int32_t *w;
  unsigned int i;

  p = &prepared_s1_12289[0][0];
  for(i = 0; i < L*N; ++i)
    p[i] = 0;

  p = &prepared_s2_12289[0][0];
  for(i = 0; i < K*N; ++i)
    p[i] = 0;

  p = input_b;
  for(i = 0; i < N; ++i)
    p[i] = 0;

  p = ntt_product;
  for(i = 0; i < N; ++i)
    p[i] = 0;

  w = ntt_workspace;
  for(i = 0; i < N; ++i)
    w[i] = 0;
}
