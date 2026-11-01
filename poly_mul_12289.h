#ifndef POLY_MUL_12289_H
#define POLY_MUL_12289_H

#include "polyvec.h"

#define poly_mul_12289 DILITHIUM_NAMESPACE(_poly_mul_12289)
/*
 * Uses shared static NTT buffers to keep large arrays off the Cortex-M0
 * stack. These multiplication entry points are not reentrant and must not
 * be called concurrently or from a preempting interrupt.
 */
void poly_mul_12289(poly *c, const poly *a, const poly *b);

#define polyvecl_mul_challenge_12289 \
        DILITHIUM_NAMESPACE(_polyvecl_mul_challenge_12289)
void polyvecl_mul_challenge_12289(polyvecl *r,
                                  const poly *challenge,
                                  const polyvecl *s1);

#define polyveck_mul_challenge_12289 \
        DILITHIUM_NAMESPACE(_polyveck_mul_challenge_12289)
void polyveck_mul_challenge_12289(polyveck *r,
                                  const poly *challenge,
                                  const polyveck *s2);

#define poly_mul_12289_prepare_secrets \
        DILITHIUM_NAMESPACE(_poly_mul_12289_prepare_secrets)
void poly_mul_12289_prepare_secrets(const polyvecl *s1,
                                    const polyveck *s2);

#define poly_mul_12289_prepare_challenge \
        DILITHIUM_NAMESPACE(_poly_mul_12289_prepare_challenge)
void poly_mul_12289_prepare_challenge(const poly *challenge);

#define polyvecl_mul_prepared_12289 \
        DILITHIUM_NAMESPACE(_polyvecl_mul_prepared_12289)
void polyvecl_mul_prepared_12289(polyvecl *r);

#define polyveck_mul_prepared_12289 \
        DILITHIUM_NAMESPACE(_polyveck_mul_prepared_12289)
void polyveck_mul_prepared_12289(polyveck *r);

#define poly_mul_12289_clear_prepared_secrets \
        DILITHIUM_NAMESPACE(_poly_mul_12289_clear_prepared_secrets)
void poly_mul_12289_clear_prepared_secrets(void);

extern volatile unsigned int poly_mul_12289_forward_calls;

#endif
