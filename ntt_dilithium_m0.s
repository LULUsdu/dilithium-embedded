    .syntax unified
    .cpu cortex-m0
    .thumb

    .equ MONT_QINV, 58728449
    .equ DILITHIUM_Q, 8380417

    /* Independent stage3 experiments; zero preserves the stage2 binary. */
    .ifndef NTT_QINV_IN_LR
    .equ NTT_QINV_IN_LR, 0
    .endif
    .ifndef NTT_SPLIT_PARENT_ZETA
    .equ NTT_SPLIT_PARENT_ZETA, 0
    .endif
    .ifndef INVNTT_PRECOMPUTE_SCALE
    .equ INVNTT_PRECOMPUTE_SCALE, 0
    .endif
    /* sp+16 holds this only in the four-specialized-stage inverse path. */
    .equ SCALE_QINV_PRODUCT, 0xff7fe3fa /* low32(41978 * 58728449) */

    .text
    .align 2

    .extern dilithium_zetas
    .extern dilithium_merge_zetas
    .extern dilithium_inv_merge_zetas

    /*
     * Exact signedHigh(t * Q), Q = 2^23 - 2^13 + 1.
     *
     * OUTLO/OUTHI form t<<23, LO13/HI13 form t<<13 and SIGN is the
     * sign-extension word of t.  The flag-setting word-pair subtraction and
     * addition propagate the low-word borrow/carry into OUTHI.  Three-operand
     * Thumb-1 immediate shifts avoid a separate MOV at every destination.
     */
    .macro NTT_SPARSE_Q_SIGNED_HIGH T, OUTLO, OUTHI, LO13, HI13, SIGN
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
     * Final-scale Montgomery reduction with CONST=41978 and FACTOR resident.
     * FACTOR is QINV normally, or low32(41978*QINV) in the scale experiment.
     * Since CONST is positive and fits 16 bits, for signed a:
     *   high(a*CONST) = ((a ASR 16)*CONST
     *                    + high16((a & 0xffff)*CONST)) ASR 16.
     * r6/r7 are deliberately outside the clobber set and remain live across
     * consecutive lanes. CORR is the current lane register whose input dies.
     */
    .macro NTT_FINAL_SCALE_MONT CONST, FACTOR, CORR
    .if INVNTT_SCALE_CALL_SIGNED16
    mov     r1, \CONST
    bl      ntt_montgomery_modmul_direct_signed16_cortex_m0
    .else
    mov     r2, r0
    .if !INVNTT_PRECOMPUTE_SCALE
    muls    r2, \CONST
    .endif
    muls    r2, \FACTOR
    mov     r12, r2

    uxth    r2, r0
    asrs    r3, r0, #16
    mov     r4, r2
    muls    r4, \CONST
    lsrs    r4, r4, #16
    muls    r3, \CONST
    adds    r3, r3, r4
    asrs    r3, r3, #16
    mov     \CORR, r3

    mov     r0, r12
    NTT_SPARSE_Q_SIGNED_HIGH r0, r1, r2, r3, r4, r5
    mov     r0, \CORR
    subs    r0, r0, r2
    .endif
    .endm

    .ifndef INVNTT_FINAL_SCALE_PATH
    .equ INVNTT_FINAL_SCALE_PATH, 0
    .endif

    .ifndef INVNTT_UNROLL_H16
    .equ INVNTT_UNROLL_H16, 4
    .endif
    .ifndef INVNTT_UNROLL_H64
    .equ INVNTT_UNROLL_H64, 1
    .endif
    .ifndef INVNTT_UNROLL_H256
    .equ INVNTT_UNROLL_H256, 4
    .endif

    /* One lane with an assembler-time offset inside a four-lane batch. */
    .macro INVNTT_FIXED_4WAY_UNIT HBYTES, OFF, SCALE
    ldr     r0, [sp, #0]
    .if \OFF
    adds    r0, r0, #\OFF
    .endif
    ldr     r2, [r0]
    mov     r8, r2
    .if \HBYTES == 256
    movs    r1, #128
    lsls    r1, r1, #1
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r9, r2
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r10, r2
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r11, r2

    INVNTT_MERGED_2LAYER_UNIT

    .if \SCALE
    .if INVNTT_FINAL_SCALE_PATH == 0
    mov     r0, r8
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r8
    mov     r8, r0
    mov     r0, r9
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r9
    mov     r9, r0
    mov     r0, r10
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r10
    mov     r10, r0
    mov     r0, r11
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r11
    mov     r11, r0
    .elseif INVNTT_FINAL_SCALE_PATH == 1
    ldr     r6, [sp, #20]
    .if INVNTT_PRECOMPUTE_SCALE
    ldr     r7, [sp, #16]
    .else
    ldr     r7, [sp, #40]
    .endif
    mov     r0, r8
    NTT_FINAL_SCALE_MONT r6, r7, r8
    mov     r8, r0
    mov     r0, r9
    NTT_FINAL_SCALE_MONT r6, r7, r9
    mov     r9, r0
    mov     r0, r10
    NTT_FINAL_SCALE_MONT r6, r7, r10
    mov     r10, r0
    mov     r0, r11
    NTT_FINAL_SCALE_MONT r6, r7, r11
    mov     r11, r0
    .elseif INVNTT_FINAL_SCALE_PATH == 2
    ldr     r6, [sp, #20]
    .if INVNTT_PRECOMPUTE_SCALE
    ldr     r7, [sp, #16]
    .else
    ldr     r7, [sp, #40]
    .endif
    mov     r0, r8
    NTT_FINAL_SCALE_MONT r6, r7, r8
    mov     r8, r0
    mov     r0, r9
    NTT_FINAL_SCALE_MONT r6, r7, r9
    mov     r9, r0
    ldr     r6, [sp, #20]
    .if INVNTT_PRECOMPUTE_SCALE
    ldr     r7, [sp, #16]
    .else
    ldr     r7, [sp, #40]
    .endif
    mov     r0, r10
    NTT_FINAL_SCALE_MONT r6, r7, r10
    mov     r10, r0
    mov     r0, r11
    NTT_FINAL_SCALE_MONT r6, r7, r11
    mov     r11, r0
    .else
    .error "INVNTT_FINAL_SCALE_PATH must be 0, 1 or 2"
    .endif
    .endif

    ldr     r0, [sp, #0]
    .if \OFF
    adds    r0, r0, #\OFF
    .endif
    mov     r2, r8
    str     r2, [r0]
    .if \HBYTES == 256
    movs    r1, #128
    lsls    r1, r1, #1
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r9
    str     r2, [r0]
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r10
    str     r2, [r0]
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r11
    str     r2, [r0]
    adds    r0, r0, #4
    mov     r3, r0
    .endm

    /*
     * Set to zero to build the retained stage-by-stage forward NTT.
     * The inverse NTT is independent of this switch.
     */
    .ifndef NTT_USE_2LAYER_MERGE
    .equ NTT_USE_2LAYER_MERGE, 1
    .endif

    /* Independent experiment switches. All zero preserves the default path.
     * BL obeys AAPCS: r0-r3/r12/LR/flags are dead at these boundaries;
     * live butterfly lanes are in r8-r11. The NTT prologues save LR and
     * keep SP aligned to 8 bytes. No global precompute slot is accessed.
     */
    .ifndef NTT_CALL_SIGNED16
    .equ NTT_CALL_SIGNED16, 0
    .endif
    .ifndef INVNTT_CALL_SIGNED16
    .equ INVNTT_CALL_SIGNED16, 0
    .endif
    .ifndef INVNTT_SCALE_CALL_SIGNED16
    .equ INVNTT_SCALE_CALL_SIGNED16, 0
    .endif
    /* Explicitly separate inlining from CALL: a true CALL flag always BLs. */
    .ifndef NTT_INLINE_SIGNED16
    .equ NTT_INLINE_SIGNED16, 0
    .endif
    .ifndef INVNTT_INLINE_SIGNED16
    .equ INVNTT_INLINE_SIGNED16, 0
    .endif
    .if (NTT_INLINE_SIGNED16 != 0) && (NTT_INLINE_SIGNED16 != 1)
    .error "NTT_INLINE_SIGNED16 must be 0 or 1"
    .endif
    .if (INVNTT_INLINE_SIGNED16 != 0) && (INVNTT_INLINE_SIGNED16 != 1)
    .error "INVNTT_INLINE_SIGNED16 must be 0 or 1"
    .endif
    .if (NTT_INLINE_SIGNED16 && NTT_CALL_SIGNED16) || (INVNTT_INLINE_SIGNED16 && INVNTT_CALL_SIGNED16)
    .error "signed16 INLINE and CALL are mutually exclusive per direction"
    .endif
    .if (NTT_CALL_SIGNED16 != 0) && (NTT_CALL_SIGNED16 != 1)
    .error "NTT_CALL_SIGNED16 must be 0 or 1"
    .endif
    .if (INVNTT_CALL_SIGNED16 != 0) && (INVNTT_CALL_SIGNED16 != 1)
    .error "INVNTT_CALL_SIGNED16 must be 0 or 1"
    .endif
    .if (INVNTT_SCALE_CALL_SIGNED16 != 0) && (INVNTT_SCALE_CALL_SIGNED16 != 1)
    .error "INVNTT_SCALE_CALL_SIGNED16 must be 0 or 1"
    .endif
    .if (NTT_CALL_SIGNED16 || NTT_INLINE_SIGNED16) && !NTT_USE_2LAYER_MERGE
    .error "signed16 experiment requires merged forward NTT"
    .endif

    /*
     * Number of leading inverse layer-pairs implemented as genuine two-layer
     * units.  Valid values are 0..4:
     *   0: retained eight-layer baseline
     *   1: merge len=1/2
     *   2: additionally merge len=4/8
     *   3: additionally merge len=16/32
     *   4: additionally merge len=64/128
     */
    .ifndef INVNTT_MERGED_PAIRS
    .equ INVNTT_MERGED_PAIRS, 0
    .endif
    .if (INVNTT_MERGED_PAIRS < 0) || (INVNTT_MERGED_PAIRS > 4)
    .error "INVNTT_MERGED_PAIRS must be in the range 0..4"
    .endif

    .ifndef INVNTT_FUSE_FINAL_SCALE
    .equ INVNTT_FUSE_FINAL_SCALE, 0
    .endif
    .if INVNTT_FUSE_FINAL_SCALE && (INVNTT_MERGED_PAIRS != 4)
    .error "INVNTT_FUSE_FINAL_SCALE requires INVNTT_MERGED_PAIRS=4"
    .endif
    .ifndef INVNTT_SPECIALIZED_STAGES
    .equ INVNTT_SPECIALIZED_STAGES, 4
    .endif
    .if (INVNTT_SPECIALIZED_STAGES < 0) || (INVNTT_SPECIALIZED_STAGES > 4)
    .error "INVNTT_SPECIALIZED_STAGES must be in the range 0..4"
    .endif
    .if (INVNTT_SPECIALIZED_STAGES > 0) && (INVNTT_MERGED_PAIRS != 4)
    .error "specialized invNTT stages require INVNTT_MERGED_PAIRS=4"
    .endif
    .ifndef INVNTT_FIRST_MERGE_STAGE_STOP
    .equ INVNTT_FIRST_MERGE_STAGE_STOP, 0
    .endif

    .if INVNTT_CALL_SIGNED16 || INVNTT_INLINE_SIGNED16 || INVNTT_SCALE_CALL_SIGNED16
    .if (INVNTT_SPECIALIZED_STAGES != 4) || !INVNTT_FUSE_FINAL_SCALE
    .error "signed16 experiment requires four specialized inverse pairs and fused scale"
    .endif
    .endif

    .if NTT_QINV_IN_LR || NTT_SPLIT_PARENT_ZETA || INVNTT_PRECOMPUTE_SCALE
    .if !NTT_INLINE_SIGNED16 || !INVNTT_INLINE_SIGNED16 || INVNTT_SCALE_CALL_SIGNED16
    .error "stage3 experiments require inline butterflies and dedicated scale"
    .endif
    .if (INVNTT_SPECIALIZED_STAGES != 4) || (INVNTT_FINAL_SCALE_PATH != 1) || !INVNTT_FUSE_FINAL_SCALE
    .error "stage3 experiments require four specialized stages and scale path 1"
    .endif
    .endif

    /*
     * Fixed NTT frame, relative to SP after the prologue:
     *   [sp, #0]  current low pointer
     *   [sp, #4]  current high pointer
     *   [sp, #8]  end of current low half
     *   [sp, #12] end of coefficient array
     *   [sp, #16] byte offset between group halves
     *   [sp, #20] baseline Montgomery temporary, or fused scale constant
     *   [sp, #24] immutable coefficient-array base
     *
     * Persistent high registers:
     *   r8  = MONT_QINV
     *   r9  = DILITHIUM_Q
     *   r10 = current zeta
     *   r11 = next zeta pointer
     */

    /*
     * Fused signed Montgomery multiplication and butterfly.
     *
     * The signed-high products are reconstructed from four unsigned 16x16
     * products and the sign correction.  This is the same arithmetic as
     * ntt_montgomery_modmul_precompute_b_new0715_sched, but zeta, QINV and Q
     * arrive through persistent registers and there is no per-butterfly call,
     * register save, stack allocation or global precompute-b access.
     */
    .macro NTT_FUSED_BUTTERFLY OFF
    ldr     r6, [sp, #4]
    ldr     r0, [r6, #\OFF]
    mov     r1, r10

    /* t1 = low32(a * zeta) * MONT_QINV. */
    mov     r4, r8
    mov     r2, r0
    muls    r2, r1
    muls    r2, r4
    str     r2, [sp, #20]

    /* Signed-high correction for a * zeta. */
    mov     r2, r0
    asrs    r2, r2, #31
    mov     r3, r1
    ands    r3, r2
    mov     r2, r1
    asrs    r2, r2, #31
    mov     r4, r0
    ands    r4, r2
    adds    r3, r3, r4
    mov     r12, r3

    /* unsignedHigh(a * zeta), using four 16x16 partial products. */
    uxth    r2, r0
    lsrs    r3, r0, #16
    uxth    r4, r1
    lsrs    r5, r1, #16

    mov     r6, r2
    muls    r6, r4
    mov     r7, r2
    muls    r7, r5
    mov     r0, r3
    muls    r0, r4
    mov     r1, r3
    muls    r1, r5

    mov     r4, r6
    lsrs    r4, r4, #16
    uxth    r2, r7
    adds    r4, r4, r2
    uxth    r2, r0
    adds    r4, r4, r2

    lsrs    r7, r7, #16
    adds    r1, r1, r7
    lsrs    r0, r0, #16
    adds    r1, r1, r0
    lsrs    r4, r4, #16
    adds    r1, r1, r4

    mov     r2, r12
    subs    r1, r1, r2
    mov     r7, r1

    /* signedHigh(t1 * Q), with Q resident in r9. */
    ldr     r0, [sp, #20]
    mov     r5, r9
    mov     r2, r0
    asrs    r2, r2, #31
    ands    r2, r5
    mov     r12, r2

    uxth    r1, r0
    lsrs    r3, r0, #16
    uxth    r4, r5
    lsrs    r5, r5, #16

    mov     r0, r1
    muls    r0, r4
    mov     r2, r1
    muls    r2, r5
    mov     r1, r3
    muls    r1, r4
    muls    r3, r5

    mov     r4, r0
    lsrs    r4, r4, #16
    uxth    r5, r2
    adds    r4, r4, r5
    uxth    r5, r1
    adds    r4, r4, r5

    lsrs    r2, r2, #16
    adds    r3, r3, r2
    lsrs    r1, r1, #16
    adds    r3, r3, r1
    lsrs    r4, r4, #16
    adds    r3, r3, r4

    mov     r2, r12
    subs    r3, r3, r2
    mov     r0, r7
    subs    r0, r0, r3

    /* Fuse the Montgomery result with the butterfly and stores. */
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    ldr     r1, [r6, #\OFF]
    mov     r2, r1
    subs    r2, r2, r0
    adds    r1, r1, r0
    str     r2, [r7, #\OFF]
    str     r1, [r6, #\OFF]
    .endm

    /*
     * Inlined signed Montgomery multiplication for a Cortex-M0 merged unit.
     *
     * Input:
     *   r0 = coefficient, r1 = zeta
     * Output:
     *   r0 = MontgomeryReduce(coefficient * zeta)
     * Clobbers:
     *   r1-r7, r12 and CORR.
     *
     * TEMP is deliberately a macro parameter so each call site can share the
     * same arithmetic while respecting the live butterfly values around it.
     */
    /* Signed16 core for butterfly boundaries only.
     * In: r0=a, r1=b; [sp,#40]=MONT_QINV in the existing outer NTT frame.
     * QINV_IN_LR uses LR instead, initialized after caller LR is saved.
     * CACHED=1 takes b0 from r12 and signed b1 from r7; r1 is unused.
     * The shared parent cache is live for exactly two consecutive cores.
     * Out: r0=MontgomeryReduce(a*b), bit-identical to the public interface.
     * Clobbers: r1-r6 and APSR. Preserves r7-r12, LR and SP; no stores/BL.
     * No incoming flags are consumed. r6 holds C.hi across the sparse-Q
     * expansion; its scratch r1-r5 cannot overwrite C.hi. In particular,
     * no CORR/high lane is borrowed, even when its current value is dead.
     */
    .macro NTT_SIGNED16_CORE_INLINE CACHED=0
    uxth    r2, r0
    asrs    r3, r0, #16
    .if \CACHED
    mov     r4, r12             /* cached unsigned b0; r7 holds signed b1 */
    .else
    uxth    r4, r1
    asrs    r5, r1, #16
    .endif
    mov     r6, r2
    muls    r6, r4              /* p00 = a0*b0 */
    .if \CACHED
    muls    r2, r7
    .else
    muls    r2, r5              /* p01 = a0*b1 */
    .endif
    muls    r4, r3              /* p10 = b0*a1 */
    .if \CACHED
    muls    r3, r7
    .else
    muls    r3, r5              /* p11 = a1*b1 */
    .endif
    lsrs    r0, r6, #16
    uxth    r1, r2
    adds    r0, r0, r1
    uxth    r1, r4
    adds    r0, r0, r1
    asrs    r2, r2, #16
    adds    r3, r3, r2
    asrs    r4, r4, #16
    adds    r3, r3, r4
    lsrs    r1, r0, #16
    adds    r3, r3, r1          /* C.hi includes the middle carry */
    uxth    r1, r6
    lsls    r0, r0, #16
    orrs    r0, r1              /* C.lo from these same four products */
    mov     r6, r3
    .if NTT_QINV_IN_LR
    mov     r2, lr
    .else
    ldr     r2, [sp, #40]
    .endif
    muls    r0, r2
    NTT_SPARSE_Q_SIGNED_HIGH r0, r1, r2, r3, r4, r5
    subs    r0, r6, r2
    .endm

    .macro NTT_MERGED_MONTGOMERY_INLINE TEMP, CORR, CALL=INVNTT_SCALE_CALL_SIGNED16, INLINE=0, CACHED=0
    .if \INLINE
    NTT_SIGNED16_CORE_INLINE \CACHED
    .elseif \CALL
    bl      ntt_montgomery_modmul_direct_signed16_cortex_m0
    .else
    /* t1 = low32(a * zeta) * MONT_QINV. */
    ldr     r4, [sp, #40]
    mov     r2, r0
    muls    r2, r1
    muls    r2, r4
    mov     r12, r2

    /* correction = ((a >> 31) & zeta) + ((zeta >> 31) & a). */
    mov     r2, r0
    asrs    r2, r2, #31
    mov     r3, r1
    ands    r3, r2
    mov     r2, r1
    asrs    r2, r2, #31
    mov     r4, r0
    ands    r4, r2
    adds    r3, r3, r4
    mov     \CORR, r3

    /* unsignedHigh(a * zeta), from four aligned 16x16 products. */
    uxth    r2, r0
    lsrs    r3, r0, #16
    uxth    r4, r1
    lsrs    r5, r1, #16

    mov     r6, r2
    muls    r6, r4
    mov     r7, r2
    muls    r7, r5
    mov     r0, r3
    muls    r0, r4
    mov     r1, r3
    muls    r1, r5

    mov     r4, r6
    lsrs    r4, r4, #16
    uxth    r2, r7
    adds    r4, r4, r2
    uxth    r2, r0
    adds    r4, r4, r2

    lsrs    r7, r7, #16
    adds    r1, r1, r7
    lsrs    r0, r0, #16
    adds    r1, r1, r0
    lsrs    r4, r4, #16
    adds    r1, r1, r4

    mov     r2, \CORR
    subs    r1, r1, r2
    mov     r7, r1

    /*
     * Exact signedHigh(t1 * Q), exploiting:
     *   Q = 2^23 - 2^13 + 1.
     *
     * Build the signed 64-bit product as high/low word pairs. SUBS/SBCS
     * propagate the low-word borrow; ADDS/ADCS propagate the final carry.
     */
    mov     r0, r12
    NTT_SPARSE_Q_SIGNED_HIGH r0, r1, r2, r3, r4, r5

    mov     r0, r7
    subs    r0, r0, r2
    .endif
    .endm

    /*
     * One genuine two-layer, four-coefficient forward unit.
     *
     * Persistent input:
     *   r8-r11 = x0, x1, x2, x3, then the four live butterfly values
     *   [sp,#8/#12/#16] = parent, left-child, right-child zetas
     *
     * The high registers survive each inlined Montgomery core, so no
     * butterfly intermediate is written to the stack. [sp,#24] is the only
     * per-core temporary. NTT_CALL_SIGNED16 selects direct BL calls;
     * NTT_INLINE_SIGNED16 selects the signed16 core without calls. Neither
     * uses CORR or accesses the global precompute-b slot.
     */
    .macro NTT_MERGED_2LAYER_UNIT
    /* Parent butterflies: (x0,x2) and (x1,x3). */
    ldr     r1, [sp, #8]
    .if NTT_SPLIT_PARENT_ZETA
    uxth    r2, r1
    mov     r12, r2
    asrs    r7, r1, #16
    .endif
    mov     r0, r10
    NTT_MERGED_MONTGOMERY_INLINE 24, r10, NTT_CALL_SIGNED16, NTT_INLINE_SIGNED16, NTT_SPLIT_PARENT_ZETA
    mov     r1, r8
    mov     r2, r1
    adds    r1, r1, r0
    subs    r2, r2, r0
    mov     r8, r1
    mov     r10, r2

    .if !NTT_SPLIT_PARENT_ZETA
    ldr     r1, [sp, #8]
    .endif
    mov     r0, r11
    NTT_MERGED_MONTGOMERY_INLINE 24, r11, NTT_CALL_SIGNED16, NTT_INLINE_SIGNED16, NTT_SPLIT_PARENT_ZETA
    mov     r1, r9
    mov     r2, r1
    adds    r1, r1, r0
    subs    r2, r2, r0
    mov     r9, r1
    mov     r11, r2

    /* Child butterfly (y0,y1), producing final x0/x1 positions. */
    ldr     r1, [sp, #12]
    mov     r0, r9
    NTT_MERGED_MONTGOMERY_INLINE 24, r9, NTT_CALL_SIGNED16, NTT_INLINE_SIGNED16
    mov     r1, r8
    mov     r2, r1
    adds    r1, r1, r0
    subs    r2, r2, r0
    mov     r8, r1
    mov     r9, r2

    /* Child butterfly (y2,y3), producing final x2/x3 positions. */
    ldr     r1, [sp, #16]
    mov     r0, r11
    NTT_MERGED_MONTGOMERY_INLINE 24, r11, NTT_CALL_SIGNED16, NTT_INLINE_SIGNED16
    mov     r1, r10
    mov     r2, r1
    adds    r1, r1, r0
    subs    r2, r2, r0
    mov     r10, r1
    mov     r11, r2
    .endm

    .macro NTT_STAGE_4WAY NAME, OFFSET
\NAME\()_start:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    movs    r0, #(\OFFSET / 4)
    lsls    r0, r0, #2
    str     r0, [sp, #16]
    mov     r7, r6
    adds    r7, r7, r0
    str     r7, [sp, #4]
    str     r7, [sp, #8]

\NAME\()_group:
    mov     r5, r11
    ldr     r0, [r5]
    adds    r5, r5, #4
    mov     r11, r5
    mov     r10, r0

\NAME\()_butterfly:
    NTT_FUSED_BUTTERFLY 0
    NTT_FUSED_BUTTERFLY 4
    NTT_FUSED_BUTTERFLY 8
    NTT_FUSED_BUTTERFLY 12

    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #16
    adds    r7, r7, #16
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r4, [sp, #8]
    cmp     r6, r4
    beq     \NAME\()_group_done
    b       \NAME\()_butterfly

\NAME\()_group_done:
    ldr     r0, [sp, #16]
    adds    r6, r6, r0
    adds    r7, r7, r0
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    str     r7, [sp, #8]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     \NAME\()_done
    b       \NAME\()_group

\NAME\()_done:
    .endm

    .macro NTT_MERGED_STAGE_IMM NAME, HBYTES, NEXT
\NAME\()_start:
    ldr     r0, [sp, #28]
    str     r0, [sp, #20]
\NAME\()_group:
    ldr     r1, [sp, #32]
    ldmia   r1!, {r2-r4}
    str     r2, [sp, #8]
    str     r3, [sp, #12]
    str     r4, [sp, #16]
    str     r1, [sp, #32]
    ldr     r0, [sp, #20]
    adds    r0, r0, #\HBYTES
    str     r0, [sp, #4]
\NAME\()_unit:
    ldr     r0, [sp, #20]
    ldr     r2, [r0]
    mov     r8, r2
    adds    r0, r0, #\HBYTES
    ldr     r2, [r0]
    mov     r9, r2
    adds    r0, r0, #\HBYTES
    ldr     r2, [r0]
    mov     r10, r2
    adds    r0, r0, #\HBYTES
    ldr     r2, [r0]
    mov     r11, r2
    NTT_MERGED_2LAYER_UNIT
    ldr     r0, [sp, #20]
    mov     r2, r8
    str     r2, [r0]
    adds    r0, r0, #\HBYTES
    mov     r2, r9
    str     r2, [r0]
    adds    r0, r0, #\HBYTES
    mov     r2, r10
    str     r2, [r0]
    adds    r0, r0, #\HBYTES
    mov     r2, r11
    str     r2, [r0]
    adds    r0, r0, #4
    mov     r3, r0
    ldr     r0, [sp, #20]
    adds    r0, r0, #4
    str     r0, [sp, #20]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     81f
    b       \NAME\()_unit
81:
    mov     r0, r3
    str     r0, [sp, #20]
    ldr     r1, [sp, #0]
    cmp     r0, r1
    beq     82f
    b       \NAME\()_group
82:
    b       \NEXT
    .endm

    .global pqcrystals_dilithium3_ref_ntt
    .global ntt_asm
    .type pqcrystals_dilithium3_ref_ntt, %function
    .type ntt_asm, %function
    .thumb_func
pqcrystals_dilithium3_ref_ntt:
ntt_asm:
    .if NTT_USE_2LAYER_MERGE
    /*
     * Merged frame after saving r4-r11 and lr (36 bytes):
     *   [sp,#0]  coefficient-array end
     *   [sp,#4]  end of current four-coefficient-unit lane
     *   [sp,#8]  parent zeta
     *   [sp,#12] left-child zeta
     *   [sp,#16] right-child zeta
     *   [sp,#20] parent-group base, then current x0 unit pointer
     *   [sp,#24] unused/reserved
     *   [sp,#28] immutable coefficient-array base
     *   [sp,#32] next merged-zeta triple
     *   [sp,#36] current half-group size in bytes (256,64,16,4)
     *   [sp,#40] MONT_QINV constant
     *
     * The 44-byte frame gives 36+44=80 bytes total stack use and preserves
     * 8-byte AAPCS alignment. With NTT_CALL_SIGNED16, each BL temporarily
     * adds the target's 16-byte frame (96 bytes total); otherwise no BL.
     */
    push    {r4-r7, lr}
    mov     r4, r8
    mov     r5, r9
    mov     r6, r10
    mov     r7, r11
    push    {r4-r7}
    sub     sp, sp, #44

    ldr     r1, =MONT_QINV
    str     r1, [sp, #40]
    .if NTT_QINV_IN_LR
    mov     lr, r1
    .endif
    str     r0, [sp, #28]
    movs    r1, #128
    lsls    r1, r1, #3
    adds    r1, r1, r0
    str     r1, [sp, #0]
    ldr     r1, =dilithium_merge_zetas + 4
    str     r1, [sp, #32]

    movs    r1, #128
    lsls    r1, r1, #1
    str     r1, [sp, #36]

    /*
     * All four merged stage pairs share this loop and the same inlined unit.
     * HBYTES advances 256 -> 64 -> 16 -> 4 without duplicating the core.
     */
ntt_merged_stage_start:
    ldr     r0, [sp, #28]
    str     r0, [sp, #20]

ntt_merged_group:
    ldr     r1, [sp, #32]
    /* Aligned parent/left/right triple; writeback keeps zeta monotonic. */
    ldmia   r1!, {r2-r4}
    str     r2, [sp, #8]
    str     r3, [sp, #12]
    str     r4, [sp, #16]
    str     r1, [sp, #32]

    ldr     r0, [sp, #20]
    ldr     r1, [sp, #36]
    adds    r0, r0, r1
    str     r0, [sp, #4]

ntt_merged_unit:
    /* Load each input coefficient once; r8-r11 then hold all live values. */
    ldr     r0, [sp, #20]
    ldr     r1, [sp, #36]
    ldr     r2, [r0]
    mov     r8, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r9, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r10, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r11, r2

    NTT_MERGED_2LAYER_UNIT

    /* Write each final output once, reconstructing addresses monotonically. */
    ldr     r0, [sp, #20]
    ldr     r1, [sp, #36]
    mov     r2, r8
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r9
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r10
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r11
    str     r2, [r0]
    adds    r0, r0, #4
    mov     r3, r0

    ldr     r0, [sp, #20]
    adds    r0, r0, #4
    str     r0, [sp, #20]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     ntt_merged_group_advance
    b       ntt_merged_unit

ntt_merged_group_advance:
    /* r3 is the address immediately after the just-written x3 lane. */
    mov     r0, r3
    str     r0, [sp, #20]
    ldr     r1, [sp, #0]
    cmp     r0, r1
    beq     ntt_merged_stage_advance
    b       ntt_merged_group

ntt_merged_stage_advance:
    ldr     r1, [sp, #36]
    lsrs    r1, r1, #2
    cmp     r1, #16
    bne     94f
    b       ntt_merged_h16_start
94:
    cmp     r1, #4
    bne     95f
    b       ntt_merged_contiguous_start
95:
    cmp     r1, #1
    bne     93f
    b       ntt_merged_complete
93:
    str     r1, [sp, #36]
    b       ntt_merged_stage_start

    .ltorg

    NTT_MERGED_STAGE_IMM ntt_merged_h16, 16, ntt_merged_contiguous_start

    .ltorg

    /*
     * The final len=2/1 merged stage has HBYTES=4, so each radix-4 unit is
     * contiguous. One LDM/STM pair replaces four scalar address walks.
     */
ntt_merged_contiguous_start:
    ldr     r0, [sp, #28]
    str     r0, [sp, #20]

ntt_merged_contiguous_group:
    ldr     r1, [sp, #32]
    ldmia   r1!, {r2-r4}
    str     r2, [sp, #8]
    str     r3, [sp, #12]
    str     r4, [sp, #16]
    str     r1, [sp, #32]

    /* Four aligned inputs, each loaded exactly once. */
    ldr     r0, [sp, #20]
    ldmia   r0!, {r2-r5}
    mov     r8, r2
    mov     r9, r3
    mov     r10, r4
    mov     r11, r5

    NTT_MERGED_2LAYER_UNIT

    /* Four final outputs, each stored exactly once; r4 becomes next group. */
    ldr     r4, [sp, #20]
    mov     r0, r8
    mov     r1, r9
    mov     r2, r10
    mov     r3, r11
    stmia   r4!, {r0-r3}
    str     r4, [sp, #20]
    ldr     r0, [sp, #0]
    cmp     r4, r0
    bne     92f
    b       ntt_merged_complete
92:
    b       ntt_merged_contiguous_group

    .ltorg

ntt_merged_complete:
    add     sp, sp, #44
    pop     {r4-r7}
    mov     r8, r4
    mov     r9, r5
    mov     r10, r6
    mov     r11, r7
    pop     {r4-r7, pc}
    .ltorg
    .else
    /*
     * Save r4-r11 once and allocate one 28-byte working frame.  Total stack
     * use is 20 + 16 + 28 = 64 bytes and remains fixed for the whole NTT.
     * The fused implementation contains no BL instructions.
     */
    push    {r4-r7, lr}
    mov     r4, r8
    mov     r5, r9
    mov     r6, r10
    mov     r7, r11
    push    {r4-r7}
    sub     sp, sp, #28

    str     r0, [sp, #0]
    str     r0, [sp, #24]
    movs    r1, #128
    lsls    r1, r1, #3
    adds    r1, r1, r0
    str     r1, [sp, #12]

    ldr     r1, =MONT_QINV
    mov     r8, r1
    ldr     r1, =DILITHIUM_Q
    mov     r9, r1
    ldr     r1, =dilithium_zetas + 4
    mov     r11, r1
    b       ntt_stage_len128_start
    .ltorg

    NTT_STAGE_4WAY ntt_stage_len128, 512
    NTT_STAGE_4WAY ntt_stage_len64, 256
    NTT_STAGE_4WAY ntt_stage_len32, 128
    NTT_STAGE_4WAY ntt_stage_len16, 64
    NTT_STAGE_4WAY ntt_stage_len8, 32

ntt_stage_len4_start:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    movs    r0, #16
    mov     r7, r6
    adds    r7, r7, r0
    str     r7, [sp, #4]
ntt_stage_len4_group:
    mov     r5, r11
    ldr     r0, [r5]
    adds    r5, r5, #4
    mov     r11, r5
    mov     r10, r0
    NTT_FUSED_BUTTERFLY 0
    NTT_FUSED_BUTTERFLY 4
    NTT_FUSED_BUTTERFLY 8
    NTT_FUSED_BUTTERFLY 12
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #32
    adds    r7, r7, #32
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     ntt_stage_len2_start
    b       ntt_stage_len4_group

ntt_stage_len2_start:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    movs    r0, #8
    mov     r7, r6
    adds    r7, r7, r0
    str     r7, [sp, #4]
ntt_stage_len2_group:
    mov     r5, r11
    ldr     r0, [r5]
    adds    r5, r5, #4
    mov     r11, r5
    mov     r10, r0
    NTT_FUSED_BUTTERFLY 0
    NTT_FUSED_BUTTERFLY 4
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #16
    adds    r7, r7, #16
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     ntt_stage_len1_start
    b       ntt_stage_len2_group

ntt_stage_len1_start:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    mov     r7, r6
    adds    r7, r7, #4
    str     r7, [sp, #4]
ntt_stage_len1_group:
    mov     r5, r11
    ldr     r0, [r5]
    adds    r5, r5, #4
    mov     r11, r5
    mov     r10, r0
    NTT_FUSED_BUTTERFLY 0
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #8
    adds    r7, r7, #8
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     ntt_done
    b       ntt_stage_len1_group

ntt_done:
    add     sp, sp, #28
    pop     {r4-r7}
    mov     r8, r4
    mov     r9, r5
    mov     r10, r6
    mov     r11, r7
    pop     {r4-r7, pc}
    .endif

    .size pqcrystals_dilithium3_ref_ntt, .-pqcrystals_dilithium3_ref_ntt
    .size ntt_asm, .-ntt_asm

    /*
     * Signed Montgomery core used only inside the inverse NTT.
     *
     * Input:
     *   r0 = signed coefficient
     *   r1 = signed multiplier
     *   r8 = MONT_QINV
     *   r9 = DILITHIUM_Q
     * Output:
     *   r0 = MontgomeryReduce(r0 * r1)
     *
     * The signed high words are reconstructed from four unsigned 16x16
     * products plus the two's-complement sign correction.  The arithmetic is
     * identical to ntt_montgomery_modmul_precompute_b_new0715_sched, but all
     * constants and the multiplier are already in registers.
     */
    .macro INVNTT_SIGNED_MONTGOMERY_CORE
    /* t1 = low32(a * b) * MONT_QINV. */
    mov     r4, r8
    mov     r2, r0
    muls    r2, r1
    muls    r2, r4
    str     r2, [sp, #20]

    /* Signed-high correction for a * b. */
    mov     r2, r0
    asrs    r2, r2, #31
    mov     r3, r1
    ands    r3, r2
    mov     r2, r1
    asrs    r2, r2, #31
    mov     r4, r0
    ands    r4, r2
    adds    r3, r3, r4
    mov     r12, r3

    /* unsignedHigh(a * b), using four 16x16 partial products. */
    uxth    r2, r0
    lsrs    r3, r0, #16
    uxth    r4, r1
    lsrs    r5, r1, #16

    mov     r6, r2
    muls    r6, r4
    mov     r7, r2
    muls    r7, r5
    mov     r0, r3
    muls    r0, r4
    mov     r1, r3
    muls    r1, r5

    mov     r4, r6
    lsrs    r4, r4, #16
    uxth    r2, r7
    adds    r4, r4, r2
    uxth    r2, r0
    adds    r4, r4, r2

    lsrs    r7, r7, #16
    adds    r1, r1, r7
    lsrs    r0, r0, #16
    adds    r1, r1, r0
    lsrs    r4, r4, #16
    adds    r1, r1, r4

    mov     r2, r12
    subs    r1, r1, r2
    mov     r7, r1

    /* signedHigh(t1 * Q), with Q resident in r9. */
    ldr     r0, [sp, #20]
    mov     r5, r9
    mov     r2, r0
    asrs    r2, r2, #31
    ands    r2, r5
    mov     r12, r2

    uxth    r1, r0
    lsrs    r3, r0, #16
    uxth    r4, r5
    lsrs    r5, r5, #16

    mov     r0, r1
    muls    r0, r4
    mov     r2, r1
    muls    r2, r5
    mov     r1, r3
    muls    r1, r4
    muls    r3, r5

    mov     r4, r0
    lsrs    r4, r4, #16
    uxth    r5, r2
    adds    r4, r4, r5
    uxth    r5, r1
    adds    r4, r4, r5

    lsrs    r2, r2, #16
    adds    r3, r3, r2
    lsrs    r1, r1, #16
    adds    r3, r3, r1
    lsrs    r4, r4, #16
    adds    r3, r3, r4

    mov     r2, r12
    subs    r3, r3, r2
    mov     r0, r7
    subs    r0, r0, r3
    .endm

    /*
     * Inverse butterfly.  Unlike the forward macro, this first stores the
     * low+high sum, then Montgomery-multiplies low-high by the negative zeta.
     */
    .macro INVNTT_FUSED_BUTTERFLY OFF
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    ldr     r0, [r6, #\OFF]
    ldr     r1, [r7, #\OFF]
    mov     r2, r0
    adds    r2, r2, r1
    subs    r0, r0, r1
    str     r2, [r6, #\OFF]
    mov     r1, r10
    INVNTT_SIGNED_MONTGOMERY_CORE
    ldr     r7, [sp, #4]
    str     r0, [r7, #\OFF]
    .endm

    .macro INVNTT_FUSED_SCALE OFF
    ldr     r6, [sp, #0]
    ldr     r0, [r6, #\OFF]
    mov     r1, r10
    INVNTT_SIGNED_MONTGOMERY_CORE
    ldr     r6, [sp, #0]
    str     r0, [r6, #\OFF]
    .endm

    /*
     * Genuine cross-layer inverse two-layer unit for every merged (L,2L)
     * stage.  This is not four independent butterflies from one layer.
     *
     * Input and output:
     *   r8-r11 = x0,x1,x2,x3, then y0,y1,y2,y3
     *   output positions = a[j],a[j+L],a[j+2L],a[j+3L]
     *   [sp,#28] = already-negated left-child zeta
     *   [sp,#32] = already-negated right-child zeta
     *   [sp,#36] = already-negated parent zeta
     *   [sp,#40] = MONT_QINV
     *
     * NTT_MERGED_MONTGOMERY_INLINE is direction-independent: it computes the
     * same signed Montgomery reduction as INVNTT_SIGNED_MONTGOMERY_CORE, but
     * leaves the other high registers live.  The CORR register at each call
     * is the butterfly lane whose old value is already dead.
     */
    .macro INVNTT_MERGED_2LAYER_UNIT
    /* t0=x0+x1; t1=Mont((x0-x1)*zeta_left). */
    mov     r0, r8
    mov     r2, r0
    mov     r3, r9
    subs    r0, r0, r3
    adds    r2, r2, r3
    mov     r8, r2
    ldr     r1, [sp, #28]
    NTT_MERGED_MONTGOMERY_INLINE 24, r9, INVNTT_CALL_SIGNED16, INVNTT_INLINE_SIGNED16
    mov     r9, r0

    /* t2=x2+x3; t3=Mont((x2-x3)*zeta_right). */
    mov     r0, r10
    mov     r2, r0
    mov     r3, r11
    subs    r0, r0, r3
    adds    r2, r2, r3
    mov     r10, r2
    ldr     r1, [sp, #32]
    NTT_MERGED_MONTGOMERY_INLINE 24, r11, INVNTT_CALL_SIGNED16, INVNTT_INLINE_SIGNED16
    mov     r11, r0

    /* y0=t0+t2; y2=Mont((t0-t2)*zeta_parent). */
    mov     r0, r8
    mov     r2, r0
    mov     r3, r10
    subs    r0, r0, r3
    adds    r2, r2, r3
    mov     r8, r2
    ldr     r1, [sp, #36]
    .if NTT_SPLIT_PARENT_ZETA
    uxth    r2, r1
    mov     r12, r2
    asrs    r7, r1, #16
    .endif
    NTT_MERGED_MONTGOMERY_INLINE 24, r10, INVNTT_CALL_SIGNED16, INVNTT_INLINE_SIGNED16, NTT_SPLIT_PARENT_ZETA
    mov     r10, r0

    /* y1=t1+t3; y3=Mont((t1-t3)*zeta_parent). */
    mov     r0, r9
    mov     r2, r0
    mov     r3, r11
    subs    r0, r0, r3
    adds    r2, r2, r3
    mov     r9, r2
    .if !NTT_SPLIT_PARENT_ZETA
    ldr     r1, [sp, #36]
    .endif
    NTT_MERGED_MONTGOMERY_INLINE 24, r11, INVNTT_CALL_SIGNED16, INVNTT_INLINE_SIGNED16, NTT_SPLIT_PARENT_ZETA
    mov     r11, r0
    .endm

    /*
     * Generate one stage with a compile-time byte offset.  Each invocation
     * emits independent group and butterfly loops; there is no run-time len
     * variable or outer stage dispatcher.
     */
    .macro INVNTT_STAGE NAME, OFFSET
\NAME:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    movs    r0, #(\OFFSET / 4)
    lsls    r0, r0, #2
    str     r0, [sp, #16]
    mov     r7, r6
    adds    r7, r7, r0
    str     r7, [sp, #4]
    str     r7, [sp, #8]

\NAME\()_group:
    mov     r5, r11
    subs    r5, r5, #4
    mov     r11, r5
    ldr     r0, [r5]
    negs    r0, r0
    mov     r10, r0

\NAME\()_butterfly:
    INVNTT_FUSED_BUTTERFLY 0
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #4
    adds    r7, r7, #4
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r4, [sp, #8]
    cmp     r6, r4
    bne     \NAME\()_butterfly

    ldr     r0, [sp, #16]
    adds    r6, r6, r0
    adds    r7, r7, r0
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    str     r7, [sp, #8]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    bne     \NAME\()_group
    .endm

    /*
     * Four-butterfly stage for len >= 8.  The current group keeps one zeta
     * while offsets 0, 4, 8 and 12 are processed.  Pointer maintenance and
     * the inner-loop branch occur once per four butterflies.
     */
    .macro INVNTT_STAGE_4WAY NAME, OFFSET
\NAME:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    movs    r0, #(\OFFSET / 4)
    lsls    r0, r0, #2
    str     r0, [sp, #16]
    mov     r7, r6
    adds    r7, r7, r0
    str     r7, [sp, #4]
    str     r7, [sp, #8]

\NAME\()_group:
    mov     r5, r11
    subs    r5, r5, #4
    mov     r11, r5
    ldr     r0, [r5]
    negs    r0, r0
    mov     r10, r0

\NAME\()_butterfly4:
    INVNTT_FUSED_BUTTERFLY 0
    INVNTT_FUSED_BUTTERFLY 4
    INVNTT_FUSED_BUTTERFLY 8
    INVNTT_FUSED_BUTTERFLY 12
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #16
    adds    r7, r7, #16
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r4, [sp, #8]
    cmp     r6, r4
    beq     1f
    b       \NAME\()_butterfly4
1:

    ldr     r0, [sp, #16]
    adds    r6, r6, r0
    adds    r7, r7, r0
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    str     r7, [sp, #8]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     2f
    b       \NAME\()_group
2:
    .endm

    .macro INVNTT_STAGE_LEN1
invntt_stage_len1:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    mov     r7, r6
    adds    r7, r7, #4
    str     r7, [sp, #4]
invntt_stage_len1_group:
    mov     r5, r11
    subs    r5, r5, #4
    mov     r11, r5
    ldr     r0, [r5]
    negs    r0, r0
    mov     r10, r0
    INVNTT_FUSED_BUTTERFLY 0
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #8
    adds    r7, r7, #8
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     1f
    b       invntt_stage_len1_group
1:
    .endm

    .macro INVNTT_STAGE_LEN2
invntt_stage_len2:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    mov     r7, r6
    adds    r7, r7, #8
    str     r7, [sp, #4]
invntt_stage_len2_group:
    mov     r5, r11
    subs    r5, r5, #4
    mov     r11, r5
    ldr     r0, [r5]
    negs    r0, r0
    mov     r10, r0
    INVNTT_FUSED_BUTTERFLY 0
    INVNTT_FUSED_BUTTERFLY 4
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #16
    adds    r7, r7, #16
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     1f
    b       invntt_stage_len2_group
1:
    .endm

    .macro INVNTT_STAGE_LEN4
invntt_stage_len4:
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
    mov     r7, r6
    adds    r7, r7, #16
    str     r7, [sp, #4]
invntt_stage_len4_group:
    mov     r5, r11
    subs    r5, r5, #4
    mov     r11, r5
    ldr     r0, [r5]
    negs    r0, r0
    mov     r10, r0
    INVNTT_FUSED_BUTTERFLY 0
    INVNTT_FUSED_BUTTERFLY 4
    INVNTT_FUSED_BUTTERFLY 8
    INVNTT_FUSED_BUTTERFLY 12
    ldr     r6, [sp, #0]
    ldr     r7, [sp, #4]
    adds    r6, r6, #32
    adds    r7, r7, #32
    str     r6, [sp, #0]
    str     r7, [sp, #4]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     1f
    b       invntt_stage_len4_group
1:
    .endm

    /*
     * One compile-time-specialized inverse two-layer stage.  HBYTES is an
     * assembler constant: no hot-path stage-width stack load is emitted.
     * SCALE is nonzero only for the public final HBYTES=256 stage.
     */
    .macro INVNTT_MERGED_STAGE_FIXED NAME, HBYTES, NEXT, SCALE
\NAME:
    ldr     r0, [sp, #24]
    str     r0, [sp, #0]

\NAME\()_group:
    ldr     r0, [sp, #8]
    ldmia   r0!, {r2-r4}
    str     r0, [sp, #8]
    str     r2, [sp, #28]
    str     r3, [sp, #32]
    str     r4, [sp, #36]

    .if (\HBYTES != 4) && !((\HBYTES == 16) && (INVNTT_UNROLL_H16 == 4))
    ldr     r0, [sp, #0]
    .if \HBYTES == 256
    movs    r1, #128
    lsls    r1, r1, #1
    adds    r2, r0, r1
    .else
    mov     r2, r0
    adds    r2, r2, #\HBYTES
    .endif
    str     r2, [sp, #4]
    .endif

\NAME\()_unit:
    .if (\HBYTES == 16) && (INVNTT_UNROLL_H16 == 4)
    INVNTT_FIXED_4WAY_UNIT 16, 0, 0
    INVNTT_FIXED_4WAY_UNIT 16, 4, 0
    INVNTT_FIXED_4WAY_UNIT 16, 8, 0
    INVNTT_FIXED_4WAY_UNIT 16, 12, 0
    b       \NAME\()_group_advance
    .elseif (\HBYTES == 16) && (INVNTT_UNROLL_H16 == 2)
    INVNTT_FIXED_4WAY_UNIT 16, 0, 0
    INVNTT_FIXED_4WAY_UNIT 16, 4, 0
    ldr     r0, [sp, #0]
    adds    r0, r0, #8
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    b       \NAME\()_unit
    .elseif (\HBYTES == 64) && (INVNTT_UNROLL_H64 == 4)
    INVNTT_FIXED_4WAY_UNIT 64, 0, 0
    INVNTT_FIXED_4WAY_UNIT 64, 4, 0
    INVNTT_FIXED_4WAY_UNIT 64, 8, 0
    INVNTT_FIXED_4WAY_UNIT 64, 12, 0
    ldr     r0, [sp, #0]
    adds    r0, r0, #16
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    ldr     r0, =\NAME\()_unit
    adds    r0, r0, #1
    bx      r0
    .ltorg
    .elseif (\HBYTES == 64) && (INVNTT_UNROLL_H64 == 2)
    INVNTT_FIXED_4WAY_UNIT 64, 0, 0
    INVNTT_FIXED_4WAY_UNIT 64, 4, 0
    ldr     r0, [sp, #0]
    adds    r0, r0, #8
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    b       \NAME\()_unit
    .elseif (\HBYTES == 256) && (INVNTT_UNROLL_H256 == 4)
    INVNTT_FIXED_4WAY_UNIT 256, 0, 1
    INVNTT_FIXED_4WAY_UNIT 256, 4, 1
    INVNTT_FIXED_4WAY_UNIT 256, 8, 1
    INVNTT_FIXED_4WAY_UNIT 256, 12, 1
    ldr     r0, [sp, #0]
    adds    r0, r0, #16
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    ldr     r0, =\NAME\()_unit
    adds    r0, r0, #1
    bx      r0
    .ltorg
    .elseif (\HBYTES == 256) && (INVNTT_UNROLL_H256 == 2)
    INVNTT_FIXED_4WAY_UNIT 256, 0, 1
    INVNTT_FIXED_4WAY_UNIT 256, 4, 1
    ldr     r0, [sp, #0]
    adds    r0, r0, #8
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    b       \NAME\()_unit
    .elseif (\HBYTES == 256) && (INVNTT_UNROLL_H256 == 1)
    INVNTT_FIXED_4WAY_UNIT 256, 0, 1
    ldr     r0, [sp, #0]
    adds    r0, r0, #4
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    b       \NAME\()_unit
    .else
    ldr     r0, [sp, #0]
    ldr     r2, [r0]
    mov     r8, r2
    .if \HBYTES == 256
    movs    r1, #128
    lsls    r1, r1, #1
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r9, r2
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r10, r2
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    ldr     r2, [r0]
    mov     r11, r2

    INVNTT_MERGED_2LAYER_UNIT

    .if \SCALE
    mov     r0, r8
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r8
    mov     r8, r0
    mov     r0, r9
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r9
    mov     r9, r0
    mov     r0, r10
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r10
    mov     r10, r0
    mov     r0, r11
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r11
    mov     r11, r0
    .endif

    ldr     r0, [sp, #0]
    mov     r2, r8
    str     r2, [r0]
    .if \HBYTES == 256
    movs    r1, #128
    lsls    r1, r1, #1
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r9
    str     r2, [r0]
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r10
    str     r2, [r0]
    .if \HBYTES == 256
    adds    r0, r0, r1
    .else
    adds    r0, r0, #\HBYTES
    .endif
    mov     r2, r11
    str     r2, [r0]
    adds    r0, r0, #4
    mov     r3, r0

    .if \HBYTES == 4
    b       \NAME\()_group_advance
    .else
    ldr     r0, [sp, #0]
    adds    r0, r0, #4
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     \NAME\()_group_advance
    b       \NAME\()_unit
    .endif
    .endif

\NAME\()_group_advance:
    mov     r0, r3
    str     r0, [sp, #0]
    ldr     r1, [sp, #12]
    cmp     r0, r1
    beq     \NAME\()_complete
    .if (\HBYTES == 16) || (\HBYTES == 256) || ((\HBYTES == 64) && (INVNTT_UNROLL_H64 == 4))
    ldr     r0, =\NAME\()_group
    adds    r0, r0, #1
    bx      r0
    .ltorg
    .else
    b       \NAME\()_group
    .endif

\NAME\()_complete:
    b       \NEXT
    .endm

    /*
     * Inverse NTT fixed frame, relative to SP after the prologue:
     *   [sp, #0]  current low pointer
     *   [sp, #4]  current high pointer
     *   [sp, #8]  end of current low half, or merged-zeta pointer
     *   [sp, #12] end of coefficient array
     *   [sp, #16] byte offset between group halves
     *   [sp, #20] Montgomery temporary t1
     *   [sp, #24] immutable coefficient-array base
     *   [sp, #28] merged left-child zeta
     *   [sp, #32] merged right-child zeta
     *   [sp, #36] merged parent zeta
     *   [sp, #40] MONT_QINV for merged Montgomery inline
     *   [sp, #44] HBYTES value after the last enabled merged stage
     *   [sp, #48] exact expected end of the enabled merged-zeta prefix
     *
     * Persistent high registers in the retained stage-by-stage path:
     *   r8  = MONT_QINV
     *   r9  = DILITHIUM_Q
     *   r10 = current negative zeta, then 41978
     *   r11 = reverse zeta pointer
     *
     * The 52-byte local frame plus 36 saved-register bytes uses 88 bytes and
     * preserves 8-byte AAPCS alignment.  The merged unit temporarily uses
     * r8-r11 for four live coefficients and restores the baseline persistent
     * constants before entering any retained radix-2 suffix.
     */
    .section .text.invntt_asm,"ax",%progbits
    .align 2
    .global pqcrystals_dilithium3_ref_invntt_tomont
    .global invntt_asm
    .type pqcrystals_dilithium3_ref_invntt_tomont, %function
    .type invntt_asm, %function
    .thumb_func
pqcrystals_dilithium3_ref_invntt_tomont:
invntt_asm:
    push    {r4-r7, lr}
    mov     r4, r8
    mov     r5, r9
    mov     r6, r10
    mov     r7, r11
    push    {r4-r7}
    .if INVNTT_SPECIALIZED_STAGES == 4
    sub     sp, sp, #44
    .else
    sub     sp, sp, #52
    .endif

    str     r0, [sp, #0]
    str     r0, [sp, #24]
    movs    r1, #128
    lsls    r1, r1, #3
    adds    r1, r1, r0
    str     r1, [sp, #12]

    ldr     r1, =MONT_QINV
    mov     r8, r1
    str     r1, [sp, #40]
    .if NTT_QINV_IN_LR
    mov     lr, r1
    .endif
    ldr     r1, =DILITHIUM_Q
    mov     r9, r1
    .if INVNTT_MERGED_PAIRS > 0
    .if INVNTT_FUSE_FINAL_SCALE
    ldr     r1, =41978
    str     r1, [sp, #20]
    .if INVNTT_PRECOMPUTE_SCALE
    ldr     r1, =SCALE_QINV_PRODUCT
    str     r1, [sp, #16]
    .endif
    .endif
    ldr     r1, =dilithium_inv_merge_zetas
    str     r1, [sp, #8]
    .if INVNTT_SPECIALIZED_STAGES > 0
    .if INVNTT_SPECIALIZED_STAGES < 4
    movs    r1, #1
    lsls    r1, r1, #10
    ldr     r2, =dilithium_inv_merge_zetas + 1020
    str     r1, [sp, #44]
    str     r2, [sp, #48]
    .endif
    b       invntt_fixed_h4
    .elseif INVNTT_MERGED_PAIRS == 1
    movs    r1, #4
    str     r1, [sp, #16]
    movs    r1, #16
    ldr     r2, =dilithium_inv_merge_zetas + 768
    .elseif INVNTT_MERGED_PAIRS == 2
    movs    r1, #64
    ldr     r2, =dilithium_inv_merge_zetas + 960
    .elseif INVNTT_MERGED_PAIRS == 3
    movs    r1, #1
    lsls    r1, r1, #8
    ldr     r2, =dilithium_inv_merge_zetas + 1008
    .elseif INVNTT_MERGED_PAIRS > 0
    movs    r1, #4
    str     r1, [sp, #16]
    movs    r1, #1
    lsls    r1, r1, #10
    ldr     r2, =dilithium_inv_merge_zetas + 1020
    .endif
    .if INVNTT_SPECIALIZED_STAGES == 0
    str     r1, [sp, #44]
    str     r2, [sp, #48]
    b       invntt_merged_stage_start
    .endif
    .else
    ldr     r1, =dilithium_zetas + 1024
    mov     r11, r1
    b       invntt_stage_len1
    .endif
    .ltorg

    .if INVNTT_SPECIALIZED_STAGES > 0
    INVNTT_MERGED_STAGE_FIXED invntt_fixed_h4, 4, invntt_after_fixed_h4, 0
invntt_after_fixed_h4:
    .if INVNTT_SPECIALIZED_STAGES > 1
    b       invntt_fixed_h16
    .else
    movs    r1, #16
    str     r1, [sp, #16]
    b       invntt_merged_stage_start
    .endif
    .endif

    .if INVNTT_SPECIALIZED_STAGES > 1
    INVNTT_MERGED_STAGE_FIXED invntt_fixed_h16, 16, invntt_after_fixed_h16, 0
invntt_after_fixed_h16:
    .if INVNTT_SPECIALIZED_STAGES > 2
    b       invntt_fixed_h64
    .else
    movs    r1, #64
    str     r1, [sp, #16]
    b       invntt_merged_stage_start
    .endif
    .endif

    .if INVNTT_SPECIALIZED_STAGES > 2
    INVNTT_MERGED_STAGE_FIXED invntt_fixed_h64, 64, invntt_after_fixed_h64, 0
invntt_after_fixed_h64:
    .if INVNTT_SPECIALIZED_STAGES > 3
    b       invntt_fixed_h256
    .else
    movs    r1, #128
    lsls    r1, r1, #1
    str     r1, [sp, #16]
    b       invntt_merged_stage_start
    .endif
    .endif

    .if INVNTT_SPECIALIZED_STAGES > 3
    INVNTT_MERGED_STAGE_FIXED invntt_fixed_h256, 256, invntt_fixed_complete, 1
invntt_fixed_complete:
    ldr     r0, [sp, #8]
    ldr     r1, =dilithium_inv_merge_zetas + 1020
    cmp     r0, r1
    beq     98f
    bkpt    #0
98:
    b       invntt_done
    .ltorg
    .endif

    .if (INVNTT_MERGED_PAIRS > 0) && (INVNTT_SPECIALIZED_STAGES < 4)
/*
 * One shared cross-layer stage loop handles HBYTES=4,16,64,256.  Each
 * logical 4L
 * coefficient block loads one [left,right,parent] triple, then advances its
 * lane j monotonically from base through base+L-1.
 */
invntt_merged_stage_start:
    ldr     r0, [sp, #24]
    str     r0, [sp, #0]

invntt_merged_group:
    /* Load one already-negated [left,right,parent] triple. */
    ldr     r0, [sp, #8]
    ldmia   r0!, {r2-r4}
    str     r0, [sp, #8]
    str     r2, [sp, #28]
    str     r3, [sp, #32]
    str     r4, [sp, #36]

    ldr     r0, [sp, #0]
    ldr     r1, [sp, #16]
    adds    r2, r0, r1
    str     r2, [sp, #4]

invntt_merged_unit:
    /* Load x0,x1,x2,x3 once at offsets 0,H,2H,3H. */
    ldr     r0, [sp, #0]
    ldr     r1, [sp, #16]
    ldr     r2, [r0]
    mov     r8, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r9, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r10, r2
    adds    r0, r0, r1
    ldr     r2, [r0]
    mov     r11, r2

    INVNTT_MERGED_2LAYER_UNIT

    /*
     * Experimental final-stage/scale fusion.  Only HBYTES=256 takes this
     * path.  It preserves the exact arithmetic order:
     *   y = Mont(butterfly_value * 41978)
     * and deliberately does not pre-combine a zeta with 41978.
     */
    .if INVNTT_FUSE_FINAL_SCALE
    ldr     r0, [sp, #16]
    lsrs    r0, r0, #8
    cmp     r0, #1
    beq     95f
    b       96f
95:

    mov     r0, r8
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r8
    mov     r8, r0

    mov     r0, r9
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r9
    mov     r9, r0

    mov     r0, r10
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r10
    mov     r10, r0

    mov     r0, r11
    ldr     r1, [sp, #20]
    NTT_MERGED_MONTGOMERY_INLINE 24, r11
    mov     r11, r0
96:
    .endif

    /* Store y0,y1,y2,y3 once at the same four monotonic addresses. */
    ldr     r0, [sp, #0]
    ldr     r1, [sp, #16]
    mov     r2, r8
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r9
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r10
    str     r2, [r0]
    adds    r0, r0, r1
    mov     r2, r11
    str     r2, [r0]
    adds    r0, r0, #4
    mov     r3, r0

    ldr     r0, [sp, #0]
    adds    r0, r0, #4
    str     r0, [sp, #0]
    ldr     r1, [sp, #4]
    cmp     r0, r1
    beq     invntt_merged_group_advance
    b       invntt_merged_unit

invntt_merged_group_advance:
    mov     r0, r3
    str     r0, [sp, #0]
    ldr     r1, [sp, #12]
    cmp     r0, r1
    beq     invntt_merged_stage_advance
    b       invntt_merged_group

invntt_merged_stage_advance:
    ldr     r0, [sp, #16]
    lsls    r0, r0, #2
    str     r0, [sp, #16]
    .if INVNTT_FIRST_MERGE_STAGE_STOP
    add     sp, sp, #52
    pop     {r4-r7}
    mov     r8, r4
    mov     r9, r5
    mov     r10, r6
    mov     r11, r7
    pop     {r4-r7, pc}
    .endif
    ldr     r1, [sp, #44]
    cmp     r0, r1
    beq     invntt_merged_complete
    b       invntt_merged_stage_start

invntt_merged_complete:
    /* Public-data bounds assertion: the enabled zeta prefix ends exactly. */
    ldr     r0, [sp, #8]
    ldr     r1, [sp, #48]
    cmp     r0, r1
    beq     97f
    bkpt    #0
97:
    ldr     r1, =MONT_QINV
    mov     r8, r1
    ldr     r1, =DILITHIUM_Q
    mov     r9, r1
    .if INVNTT_MERGED_PAIRS == 1
    ldr     r1, =dilithium_zetas + 256
    mov     r11, r1
    b       invntt_stage_len4
    .elseif INVNTT_MERGED_PAIRS == 2
    ldr     r1, =dilithium_zetas + 64
    mov     r11, r1
    b       invntt_stage_len16
    .elseif INVNTT_MERGED_PAIRS == 3
    ldr     r1, =dilithium_zetas + 16
    mov     r11, r1
    b       invntt_stage_len64
    .else
    .if INVNTT_FUSE_FINAL_SCALE
    b       invntt_done
    .else
    b       invntt_stage_scale
    .endif
    .endif
    .ltorg
    .endif

    .if INVNTT_MERGED_PAIRS == 0
    INVNTT_STAGE_LEN1
    INVNTT_STAGE_LEN2
    INVNTT_STAGE_LEN4
    INVNTT_STAGE_4WAY invntt_stage_len8, 32
    INVNTT_STAGE_4WAY invntt_stage_len16, 64
    INVNTT_STAGE_4WAY invntt_stage_len32, 128
    INVNTT_STAGE_4WAY invntt_stage_len64, 256
    INVNTT_STAGE_4WAY invntt_stage_len128, 512
    .elseif INVNTT_MERGED_PAIRS == 1
    INVNTT_STAGE_LEN4
    INVNTT_STAGE_4WAY invntt_stage_len8, 32
    INVNTT_STAGE_4WAY invntt_stage_len16, 64
    INVNTT_STAGE_4WAY invntt_stage_len32, 128
    INVNTT_STAGE_4WAY invntt_stage_len64, 256
    INVNTT_STAGE_4WAY invntt_stage_len128, 512
    .elseif INVNTT_MERGED_PAIRS == 2
    INVNTT_STAGE_4WAY invntt_stage_len16, 64
    INVNTT_STAGE_4WAY invntt_stage_len32, 128
    INVNTT_STAGE_4WAY invntt_stage_len64, 256
    INVNTT_STAGE_4WAY invntt_stage_len128, 512
    .elseif INVNTT_MERGED_PAIRS == 3
    INVNTT_STAGE_4WAY invntt_stage_len64, 256
    INVNTT_STAGE_4WAY invntt_stage_len128, 512
    .endif

invntt_stage_scale:
    ldr     r0, =41978
    mov     r10, r0
    ldr     r6, [sp, #24]
    str     r6, [sp, #0]
invntt_stage_scale_coefficient:
    INVNTT_FUSED_SCALE 0
    INVNTT_FUSED_SCALE 4
    INVNTT_FUSED_SCALE 8
    INVNTT_FUSED_SCALE 12
    ldr     r6, [sp, #0]
    adds    r6, r6, #16
    str     r6, [sp, #0]
    ldr     r0, [sp, #12]
    cmp     r6, r0
    beq     invntt_done
    b       invntt_stage_scale_coefficient

invntt_done:
    .if INVNTT_SPECIALIZED_STAGES == 4
    add     sp, sp, #44
    .else
    add     sp, sp, #52
    .endif
    pop     {r4-r7}
    mov     r8, r4
    mov     r9, r5
    mov     r10, r6
    mov     r11, r7
    pop     {r4-r7, pc}
    .ltorg

    .size pqcrystals_dilithium3_ref_invntt_tomont, .-pqcrystals_dilithium3_ref_invntt_tomont
    .size invntt_asm, .-invntt_asm

    .section .note.GNU-stack,"",%progbits
