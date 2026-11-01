#include <stdint.h>
#include <string.h>

#include "api.h"

#define MAX_MLEN 1024

#define FAIL_KEYPAIR     1
#define FAIL_SIGN        2
#define FAIL_SIGN_LEN    3
#define FAIL_OPEN        4
#define FAIL_OPEN_LEN    5
#define FAIL_OPEN_MSG    6
#define SIGN_ATTEMPT_COUNT 6

extern volatile uint32_t crypto_sign_last_attempts;

#if defined(BLINKY_MEMORY_MEASURE) || defined(BLINKY_MEMORY_BASELINE) || \
    defined(BLINKY_MEMORY_KEYPAIR) || defined(BLINKY_MEMORY_SIGN) || \
    defined(BLINKY_MEMORY_VERIFY)
/* Kept in every size-difference build so the measurement object has the same
 * link environment and its test-only control-stack reference is resolvable. */
__attribute__((aligned(8)))
uint8_t memory_measure_control_stack[2048];
#endif

#ifdef BLINKY_MEMORY_MEASURE

#define STACK_PATTERN       0xA5A5A5A5u
#define STACK_GUARD_BYTES   64u
#define CONTROL_STACK_SIZE  2048u

typedef int (*stack_measure_fn_t)(void *argument);

extern uint8_t *stack_measure_stack_base(void);
extern uint8_t *stack_measure_stack_limit(void);
extern int stack_measure_call(stack_measure_fn_t function, void *argument);
extern void stack_measure_start(void (*controller)(void));

volatile uint32_t keypair_stack_peak;
volatile uint32_t sign_stack_peak;
volatile uint32_t verify_stack_peak;
volatile uint32_t keypair_stack_overflow;
volatile uint32_t sign_stack_overflow;
volatile uint32_t verify_stack_overflow;
volatile int32_t keypair_measure_result;
volatile int32_t sign_measure_result;
volatile int32_t verify_measure_result;
volatile uint32_t sign_measure_attempts;
volatile uint32_t sign_measure_max_attempts;
volatile uint32_t measured_message_length;
volatile uint32_t memory_measurement_complete;

#endif

static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
static uint8_t sk[CRYPTO_SECRETKEYBYTES];
static uint8_t msg[MAX_MLEN];
static uint8_t sm[MAX_MLEN + CRYPTO_BYTES];
static uint8_t opened[MAX_MLEN + CRYPTO_BYTES];

volatile int test_status;
volatile int test_case_index;
volatile int test_fail_code;
volatile uint32_t sign_attempt_counts[SIGN_ATTEMPT_COUNT];

static void fill_message(size_t mlen)
{
    size_t i;

    for (i = 0; i < mlen; ++i) {
        msg[i] = (uint8_t)((i * 31u + mlen * 17u + 7u) & 0xffu);
    }
}

static int run_one_test(size_t mlen)
{
    size_t smlen = 0;
    size_t opened_len = 0;
    int ret;

    fill_message(mlen);

    ret = crypto_sign_keypair(pk, sk);
    if (ret != 0) {
        return FAIL_KEYPAIR;
    }

    ret = crypto_sign(sm, &smlen, msg, mlen, sk);
    if (ret != 0) {
        return FAIL_SIGN;
    }

    if (smlen != mlen + CRYPTO_BYTES) {
        return FAIL_SIGN_LEN;
    }

    ret = crypto_sign_open(opened, &opened_len, sm, smlen, pk);
    if (ret != 0) {
        return FAIL_OPEN;
    }

    if (opened_len != mlen) {
        return FAIL_OPEN_LEN;
    }

    if (memcmp(msg, opened, mlen) != 0) {
        return FAIL_OPEN_MSG;
    }

    return 0;
}

#ifdef BLINKY_MEMORY_MEASURE

typedef struct {
    size_t message_length;
    size_t signed_length;
    size_t opened_length;
} memory_measure_arguments_t;

static memory_measure_arguments_t memory_measure_arguments;

__attribute__((noinline))
static void memory_measurement_done_breakpoint(void)
{
    __asm volatile ("nop");
}

static void stack_fill(void)
{
    uint32_t *p = (uint32_t *)stack_measure_stack_base();
    uint32_t *end = (uint32_t *)stack_measure_stack_limit();

    while (p < end) {
        *p++ = STACK_PATTERN;
    }
}

static uint32_t stack_scan(uint32_t *overflow)
{
    const uint8_t *base = stack_measure_stack_base();
    const uint8_t *limit = stack_measure_stack_limit();
    const uint32_t *p = (const uint32_t *)base;
    const uint32_t *end = (const uint32_t *)limit;
    const uint32_t *guard_end = (const uint32_t *)(base + STACK_GUARD_BYTES);

    *overflow = 0;
    while (p < guard_end) {
        if (*p++ != STACK_PATTERN) {
            *overflow = 1;
            break;
        }
    }

    p = (const uint32_t *)base;
    while (p < end && *p == STACK_PATTERN) {
        ++p;
    }
    return (uint32_t)(limit - (const uint8_t *)p);
}

static int measure_keypair(void *argument)
{
    (void)argument;
    return crypto_sign_keypair(pk, sk);
}

static int measure_sign(void *argument)
{
    memory_measure_arguments_t *args = (memory_measure_arguments_t *)argument;
    return crypto_sign(sm, &args->signed_length, msg, args->message_length, sk);
}

static int measure_verify(void *argument)
{
    memory_measure_arguments_t *args = (memory_measure_arguments_t *)argument;
    return crypto_sign_open(opened, &args->opened_length,
                            sm, args->signed_length, pk);
}

static uint32_t run_stack_measurement(stack_measure_fn_t function,
                                      void *argument,
                                      volatile int32_t *result,
                                      volatile uint32_t *overflow)
{
    uint32_t used;
    uint32_t overflow_local;

    stack_fill();
    *result = stack_measure_call(function, argument);
    used = stack_scan(&overflow_local);
    *overflow = overflow_local;
    return used;
}

/* Runs permanently on memory_measure_control_stack with interrupts disabled. */
static void memory_measure_controller(void)
{
    static const size_t test_lengths[] = {
        0, 1, 33, 128, 135, 136, 137, 256, 1024
    };
    uint32_t used;
    size_t i;

    keypair_stack_peak = 0;
    sign_stack_peak = 0;
    verify_stack_peak = 0;
    keypair_stack_overflow = 0;
    sign_stack_overflow = 0;
    verify_stack_overflow = 0;
    sign_measure_max_attempts = 0;
    measured_message_length = 0;
    memory_measurement_complete = 0;

    used = run_stack_measurement(measure_keypair, 0,
                                 &keypair_measure_result,
                                 &keypair_stack_overflow);
    keypair_stack_peak = used;
    if (keypair_measure_result != 0 || keypair_stack_overflow != 0) {
        memory_measurement_complete = 0x80000001u;
        memory_measurement_done_breakpoint();
        for (;;) { }
    }

    for (i = 0; i < sizeof(test_lengths) / sizeof(test_lengths[0]); ++i) {
        memory_measure_arguments.message_length = test_lengths[i];
        memory_measure_arguments.signed_length = 0;
        fill_message(test_lengths[i]);

        used = run_stack_measurement(measure_sign, &memory_measure_arguments,
                                     &sign_measure_result,
                                     &sign_stack_overflow);
        sign_measure_attempts = crypto_sign_last_attempts;
        if (used > sign_stack_peak) {
            sign_stack_peak = used;
            measured_message_length = (uint32_t)test_lengths[i];
        }
        if (sign_measure_attempts > sign_measure_max_attempts) {
            sign_measure_max_attempts = sign_measure_attempts;
        }
        if (sign_measure_result != 0 || sign_stack_overflow != 0) {
            memory_measurement_complete = 0x80000002u;
            memory_measurement_done_breakpoint();
            for (;;) { }
        }

        memory_measure_arguments.opened_length = 0;
        used = run_stack_measurement(measure_verify, &memory_measure_arguments,
                                     &verify_measure_result,
                                     &verify_stack_overflow);
        if (used > verify_stack_peak) {
            verify_stack_peak = used;
        }
        if (verify_measure_result != 0 || verify_stack_overflow != 0 ||
            memory_measure_arguments.opened_length != test_lengths[i] ||
            memcmp(msg, opened, test_lengths[i]) != 0) {
            memory_measurement_complete = 0x80000003u;
            memory_measurement_done_breakpoint();
            for (;;) { }
        }
    }

    memory_measurement_complete = 1;
    memory_measurement_done_breakpoint();
    for (;;) { }
}

#endif

int main(void)
{
#ifdef BLINKY_MEMORY_MEASURE
    stack_measure_start(memory_measure_controller);
    for (;;) { }
#elif defined(BLINKY_MEMORY_BASELINE)
    test_status = 0;
    return test_status;
#elif defined(BLINKY_MEMORY_KEYPAIR)
    return crypto_sign_keypair(pk, sk);
#elif defined(BLINKY_MEMORY_SIGN)
    {
        size_t measured_smlen = 0;
        fill_message(MAX_MLEN);
        return crypto_sign(sm, &measured_smlen, msg, MAX_MLEN, sk);
    }
#elif defined(BLINKY_MEMORY_VERIFY)
    {
        size_t measured_opened_len = 0;
        return crypto_sign_open(opened, &measured_opened_len,
                                sm, CRYPTO_BYTES, pk);
    }
#else
    static const size_t test_lengths[] = {
        0, 1, 33, 128, 256, 1024
    };
    size_t i;
    int ret;

    test_status = 1;
    test_case_index = -1;
    test_fail_code = 0;

    for (i = 0; i < SIGN_ATTEMPT_COUNT; ++i) {
        sign_attempt_counts[i] = 0;
    }

    for (i = 0; i < sizeof(test_lengths) / sizeof(test_lengths[0]); ++i) {
        test_case_index = (int)i;
        ret = run_one_test(test_lengths[i]);
        if (ret != 0) {
            test_fail_code = ret;
            test_status = -1;
            return ret;
        }
        sign_attempt_counts[i] = crypto_sign_last_attempts;
    }

    test_status = 0;
    test_case_index = (int)(sizeof(test_lengths) / sizeof(test_lengths[0]));
    test_fail_code = 0;
    return 0;
#endif
}
