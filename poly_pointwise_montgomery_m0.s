    .syntax unified
    .cpu cortex-m0
    .thumb

    .equ PPM_M0_Q, 8380417
    .equ PPM_M0_QINV, 58728449
    .equ PPM_M0_POLY_BYTES, 1024

    .text
    .align 2

/*
 * Exact signedHigh(t * Q), using Q = 2^23 - 2^13 + 1.
 * This is the arithmetic used by ntt_dilithium_m0.s.
 */
    .macro PPM_M0_SPARSE_Q_SIGNED_HIGH T, OUTLO, OUTHI, LO13, HI13, SIGN
    lsls    \OUTLO, \T, #23
    asrs    \OUTHI, \T, #9
    lsls    \LO13, \T, #13
    asrs    \HI13, \T, #19
    asrs    \SIGN, \T, #31
    subs    \OUTLO, \OUTLO, \LO13
    sbcs    \OUTHI, \OUTHI, \HI13
    adds    \OUTLO, \OUTLO, \T
    adcs    \OUTHI, \OUTHI, \SIGN
    .endm

/*
 * void poly_pointwise_montgomery_asm(poly *c, const poly *a, const poly *b)
 *
 * For all 256 coefficients, compute exactly:
 *   c[i] = montgomery_reduce((int64_t)a[i] * b[i]).
 *
 * r8  = QINV, r9 = Q
 * r6 = signedHigh(a*b) across sparse-Q reduction; r10/r11 unused
 * [sp,#0/#4/#8/#12] = c/a/b/end pointers
 */
    .section .text.poly_pointwise_montgomery_asm,"ax",%progbits
    .align 2
    .global poly_pointwise_montgomery_asm
    .type poly_pointwise_montgomery_asm, %function
    .thumb_func
poly_pointwise_montgomery_asm:
    push    {r4-r7, lr}
    mov     r4, r8
    mov     r5, r9
    mov     r6, r10
    mov     r7, r11
    push    {r4-r7}
    sub     sp, sp, #20

    str     r0, [sp, #0]
    str     r1, [sp, #4]
    str     r2, [sp, #8]
    movs    r3, #128
    lsls    r3, r3, #3
    adds    r3, r3, r1
    str     r3, [sp, #12]

    ldr     r3, =PPM_M0_QINV
    mov     r8, r3
    ldr     r3, =PPM_M0_Q
    mov     r9, r3

ppm_m0_loop:
    ldr     r4, [sp, #4]
    ldr     r0, [r4]
    adds    r4, r4, #4
    str     r4, [sp, #4]

    ldr     r5, [sp, #8]
    ldr     r1, [r5]
    adds    r5, r5, #4
    str     r5, [sp, #8]

    /* Signed halfwords: a=a0+2^16*a1, b=b0+2^16*b1. */
    uxth    r2, r0
    asrs    r3, r0, #16
    uxth    r4, r1
    asrs    r5, r1, #16

    mov     r6, r2
    muls    r6, r4              /* p00 */
    muls    r2, r5              /* p01; a0 is dead */
    muls    r4, r3              /* p10; b0 is dead */
    muls    r3, r5              /* p11 */

    lsrs    r0, r6, #16
    uxth    r1, r2
    adds    r0, r0, r1
    uxth    r1, r4
    adds    r0, r0, r1          /* S=hi16(p00)+lo16(p01)+lo16(p10) */
    asrs    r2, r2, #16
    adds    r3, r3, r2
    asrs    r4, r4, #16
    adds    r3, r3, r4
    lsrs    r1, r0, #16
    adds    r3, r3, r1          /* signed C.hi including S carry */

    uxth    r1, r6
    lsls    r0, r0, #16
    orrs    r0, r1              /* C.lo from the same four products */
    mov     r6, r3              /* C.hi survives sparse-Q scratch r1-r5 */
    mov     r2, r8
    muls    r0, r2              /* t=signed32(low32(C.lo*QINV)) */

    /* signedHigh(t*Q), then REDC = signedHigh(a*b)-signedHigh(t*Q). */
    PPM_M0_SPARSE_Q_SIGNED_HIGH r0, r1, r2, r3, r4, r5
    subs    r0, r6, r2

    ldr     r1, [sp, #0]
    str     r0, [r1]
    adds    r1, r1, #4
    str     r1, [sp, #0]

    ldr     r0, [sp, #4]
    ldr     r1, [sp, #12]
    cmp     r0, r1
    bne     ppm_m0_loop

    add     sp, sp, #20
    pop     {r4-r7}
    mov     r8, r4
    mov     r9, r5
    mov     r10, r6
    mov     r11, r7
    pop     {r4-r7, pc}
    .ltorg
    .size poly_pointwise_montgomery_asm, .-poly_pointwise_montgomery_asm
