#pragma once
#include <cstdint>

struct Xoshiro256pp {
    uint64_t s[4];

    explicit Xoshiro256pp(uint64_t seed) {
        // SplitMix64 to initialize all four state words from one seed
        uint64_t x = seed;
        for (int i = 0; i < 4; i++) {
            x += 0x9e3779b97f4a7c15ULL;
            uint64_t z = x;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            s[i] = z ^ (z >> 31);
        }
    }

    uint64_t next() {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    // Uniform double in [0, 1)
    double next_double() {
        return (next() >> 11) * (1.0 / (1ULL << 53));
    }

private:
    static uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};
