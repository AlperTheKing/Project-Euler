#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 TARGET = 1'000'000'000'000'000'000;

void require(const bool condition, const char* description) {
    if (!condition) throw std::runtime_error(std::string("Check failed: ") + description);
}

u64 prefix_xor(const u64 n) {
    switch (n % 4) {
        case 0: return n;
        case 1: return 1;
        case 2: return n + 1;
        default: return 0;
    }
}

u64 block_xor(const u64 first, const u64 count, const u64 offset) {
    const u64 low = (offset - 1) % 4;
    const u64 begin = first + (offset - 1) / 4;
    const u64 high = prefix_xor(begin + count - 1) ^ (begin == 0 ? 0 : prefix_xor(begin - 1));
    return (high << 2) ^ (count % 2 == 0 ? 0 : low);
}

u64 solve(const u64 limit) {
    u64 result = 0;
    u64 first = 0;
    u64 length = 1;
    bool high = true;

    const auto consume = [&](const u64 size, const u64 offset) {
        if (limit < offset || first > (limit - offset) / 4) return false;
        const u64 count = std::min(size, (limit - offset) / 4 - first + 1);
        result ^= block_xor(first, count, offset);
        first += count;
        return count == size;
    };

    // The offsets b[n] - 4*n alternate between 6/4 blocks of length 3^k,
    // separated by 5 blocks of length (3^(k+1) + (-1)^k) / 4.
    while (consume(length, high ? 6 : 4)) {
        const u64 middle = high ? (3 * length + 1) / 4 : (3 * length - 1) / 4;
        if (!consume(middle, 5)) break;
        length *= 3;
        high = !high;
    }
    return result;
}

u64 boundary_offset(const u64 index) {
    u64 first = 1;
    u64 second = 2;
    bool high = true;
    for (;;) {
        if (index < first) return high ? 6 : 4;
        if (index < second) return 5;
        first = 3 * first + (high ? 2 : 1);
        second = 3 * second + (high ? 1 : 2);
        high = !high;
    }
}

u64 count_b(const u64 limit) {
    u64 count = limit / 4;
    if (count != 0 && limit % 4 < boundary_offset(count - 1) - 4) --count;
    return count;
}

u64 complement_term(const u64 index) {
    const u64 upper = index + 1 + index / 3;
    const u64 lower = upper - 1;
    return lower - count_b(lower) > index ? lower : upper;
}

void check_limit(const u64 limit) {
    if (limit == 0) {
        require(solve(0) == 0, "empty range");
        return;
    }
    const bool included = count_b(limit) != count_b(limit - 1);
    require((solve(limit) ^ solve(limit - 1)) == (included ? limit - 1 : 0),
        "XOR change at a large limit");
}

void check_index(const u64 index) {
    const u64 value = 4 * index + boundary_offset(index);
    require(value == complement_term(index) + complement_term(index + 1)
        + complement_term(index + 2), "three-term sum in the complementary sequence");
    check_limit(value - 1);
    check_limit(value);
    check_limit(value + 1);
}

void run_tests() {
    require(solve(0) == 0, "empty range");
    require(solve(20) == 26, "given example");

    constexpr u64 DIRECT_LIMIT = 1'000'000;
    std::vector<u64> a;
    std::vector<u64> b;
    a.reserve(DIRECT_LIMIT);
    b.reserve(DIRECT_LIMIT);
    std::size_t next = 0;
    u64 expected = 0;
    for (u64 value = 1; value <= DIRECT_LIMIT; ++value) {
        if (next < b.size() && b[next] == value) {
            expected ^= value - 1;
            ++next;
        } else {
            a.push_back(value);
            if (a.size() >= 3) {
                const std::size_t size = a.size();
                const u64 term = a[size - 3] + a[size - 2] + a[size - 1];
                const u64 index = static_cast<u64>(b.size());
                require(term == 4 * index + boundary_offset(index), "directly generated b term");
                b.push_back(term);
            }
        }
        require(solve(value) == expected, "every directly generated XOR prefix");
    }
    constexpr std::array<u64, 16> GIVEN_A{1, 2, 3, 4, 5, 7, 8, 10, 11, 13, 14, 15, 17, 18, 19, 21};
    constexpr std::array<u64, 5> GIVEN_B{6, 9, 12, 16, 20};
    require(std::equal(GIVEN_A.begin(), GIVEN_A.end(), a.begin()), "given a prefix");
    require(std::equal(GIVEN_B.begin(), GIVEN_B.end(), b.begin()), "given b prefix");

    constexpr u64 MAXIMUM = std::numeric_limits<u64>::max();
    constexpr u64 MAX_INDEX = (MAXIMUM - 6) / 4;
    u64 first = 1;
    u64 second = 2;
    bool high = true;
    while (first <= MAX_INDEX) {
        for (const u64 boundary : {first, second}) {
            const u64 left = boundary > 4 ? boundary - 4 : 0;
            const u64 right = std::min(MAX_INDEX, boundary + 4);
            for (u64 index = left; index <= right; ++index) check_index(index);
        }
        first = 3 * first + (high ? 2 : 1);
        second = 3 * second + (high ? 1 : 2);
        high = !high;
    }
    u64 state = 1017;
    for (unsigned i = 0; i < 4096; ++i) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        check_index(state % (MAX_INDEX + 1));
        check_limit(state);
    }
    for (unsigned bit = 0; bit < 64; ++bit) {
        const u64 power = u64{1} << bit;
        check_limit(power - 1);
        check_limit(power);
        check_limit(power + 1);
    }
    check_index(MAX_INDEX);
    check_limit(MAXIMUM);
    check_limit(TARGET);
    std::cout << "All checks passed; every limit through 1000000, block boundaries, "
                 "4096 large samples, and 64-bit limits.\n";
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
        } else {
            if (argc > 2) throw std::invalid_argument("Usage: Euler1017 [limit | --self-test]");
            u64 limit = TARGET;
            if (argc == 2) {
                const std::string argument(argv[1]);
                const auto parsed = std::from_chars(argument.data(), argument.data() + argument.size(), limit);
                if (parsed.ec != std::errc{} || parsed.ptr != argument.data() + argument.size()) {
                    throw std::invalid_argument("The limit must be an unsigned 64-bit integer");
                }
            }
            std::cout << solve(limit) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
