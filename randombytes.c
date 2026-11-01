#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "randombytes.h"
#include "fips202.h"

/*
 * Deterministic pseudo-random source for Keil/Cortex-M functional tests only.
 * This is not a cryptographically secure RNG and must not be used for real
 * key generation, signatures, or production security code.
 */
void randombytes(uint8_t *out, size_t outlen)
{
    static const uint8_t seed[32] = {
        0x44, 0x69, 0x6c, 0x69, 0x74, 0x68, 0x69, 0x75,
        0x6d, 0x33, 0x2d, 0x4d, 0x30, 0x2d, 0x74, 0x65,
        0x73, 0x74, 0x2d, 0x73, 0x65, 0x65, 0x64, 0x2d,
        0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x31
    };
    static uint64_t counter = 0;
    uint8_t input[40];
    uint8_t block[64];
    size_t take;
    unsigned int i;

    memcpy(input, seed, sizeof(seed));

    while (outlen > 0) {
        for (i = 0; i < 8; ++i) {
            input[32 + i] = (uint8_t)(counter >> (8 * i));
        }
        ++counter;

        shake256(block, sizeof(block), input, sizeof(input));

        take = outlen < sizeof(block) ? outlen : sizeof(block);
        memcpy(out, block, take);
        out += take;
        outlen -= take;
    }
}
