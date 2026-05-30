#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 MOD = 1'234'567'891ULL;
constexpr u64 INV2 = (MOD + 1) / 2;
constexpr u64 TARGET_N = 1'000'000'000'000'000'003ULL;

u64 normalize(const i64 value) {
    i64 result = value % static_cast<i64>(MOD);
    if (result < 0) {
        result += static_cast<i64>(MOD);
    }
    return static_cast<u64>(result);
}

u64 sub_mod(const u64 lhs, const u64 rhs) {
    return lhs >= rhs ? lhs - rhs : lhs + MOD - rhs;
}

u64 mul_mod(const u64 lhs, const u64 rhs) {
    return static_cast<u64>((static_cast<u128>(lhs) * rhs) % MOD);
}

u64 pow_mod(u64 base, u64 exponent) {
    u64 result = 1;
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exponent >>= 1;
    }
    return result;
}

struct Field {
    std::array<u64, 4> c{};

    Field() = default;

    explicit Field(const i64 c0) : c{normalize(c0), 0, 0, 0} {}

    Field(const i64 c0, const i64 c1, const i64 c2, const i64 c3)
        : c{normalize(c0), normalize(c1), normalize(c2), normalize(c3)} {}

    bool is_base() const {
        return c[1] == 0 && c[2] == 0 && c[3] == 0;
    }

    bool is_rho_multiple() const {
        return c[0] == 0 && c[2] == 0 && c[3] == 0;
    }
};

Field operator-(const Field& lhs, const Field& rhs) {
    Field result;
    for (int i = 0; i < 4; ++i) {
        result.c[static_cast<std::size_t>(i)] =
            sub_mod(lhs.c[static_cast<std::size_t>(i)], rhs.c[static_cast<std::size_t>(i)]);
    }
    return result;
}

Field operator-(const Field& value) {
    Field result;
    for (int i = 0; i < 4; ++i) {
        const u64 current = value.c[static_cast<std::size_t>(i)];
        result.c[static_cast<std::size_t>(i)] = current == 0 ? 0 : MOD - current;
    }
    return result;
}

Field operator*(const Field& lhs, const Field& rhs) {
    std::array<u128, 7> product{};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            product[static_cast<std::size_t>(i + j)] +=
                static_cast<u128>(lhs.c[static_cast<std::size_t>(i)]) *
                rhs.c[static_cast<std::size_t>(j)];
        }
    }

    for (int i = 6; i >= 4; --i) {
        product[static_cast<std::size_t>(i - 4)] += 2 * product[static_cast<std::size_t>(i)];
    }

    Field result;
    for (int i = 0; i < 4; ++i) {
        result.c[static_cast<std::size_t>(i)] =
            static_cast<u64>(product[static_cast<std::size_t>(i)] % MOD);
    }
    return result;
}

bool operator==(const Field& lhs, const Field& rhs) {
    return lhs.c == rhs.c;
}

Field square(const Field& value) {
    return value * value;
}

Field cube(const Field& value) {
    return value * value * value;
}

using Block = std::array<Field, 9>;

const Field RHO(0, 1, 0, 0);
const Field INV_RHO(0, 0, 0, static_cast<i64>(INV2));

Field get_value(const Block& block, const i64 center, const i64 index) {
    const i64 offset = index - center;
    assert(-4 <= offset && offset <= 4);
    return block[static_cast<std::size_t>(offset + 4)];
}

Field odd_value(const Block& block, const i64 center, const i64 k) {
    return get_value(block, center, k + 1) * cube(get_value(block, center, k - 1)) -
           get_value(block, center, k - 2) * cube(get_value(block, center, k));
}

Field even_value(const Block& block, const i64 center, const i64 k) {
    return get_value(block, center, k) * INV_RHO *
           (get_value(block, center, k + 2) * square(get_value(block, center, k - 1)) -
            get_value(block, center, k - 2) * square(get_value(block, center, k + 1)));
}

void advance(Block& block, i64& center, const int bit) {
    const Block previous = block;
    const i64 previous_center = center;
    const i64 next_center = 2 * center + bit;
    Block next{};

    for (int offset = -4; offset <= 4; ++offset) {
        const i64 index = next_center + offset;
        if (index % 2 == 0) {
            next[static_cast<std::size_t>(offset + 4)] = even_value(previous, previous_center, index / 2);
        } else {
            next[static_cast<std::size_t>(offset + 4)] = odd_value(previous, previous_center, (index + 1) / 2);
        }
    }

    block = next;
    center = next_center;
}

Field eds_value(const u64 n) {
    assert(n > 0);

    Block block{Field(1), -RHO, Field(-1), Field(0), Field(1), RHO, Field(-1), Field(0, -2, 0, 0), Field(-3)};
    i64 center = 1;

    int highest = 63;
    while (((n >> highest) & 1ULL) == 0) {
        --highest;
    }

    for (int bit_index = highest - 1; bit_index >= 0; --bit_index) {
        advance(block, center, static_cast<int>((n >> bit_index) & 1ULL));
    }

    return block[4];
}

u64 sequence_value(const u64 n) {
    const Field value = eds_value(n);
    if ((n & 1ULL) != 0) {
        assert(value.is_base());
        if ((((n - 1) / 2) & 1ULL) == 0) {
            return value.c[0];
        }
        return value.c[0] == 0 ? 0 : MOD - value.c[0];
    }

    assert(value.is_rho_multiple());
    if ((((n - 2) / 2) & 1ULL) == 0) {
        return value.c[1];
    }
    return value.c[1] == 0 ? 0 : MOD - value.c[1];
}

u64 direct_value(const int n) {
    std::vector<u64> a(static_cast<std::size_t>(n + 5), 0);
    a[1] = 1;
    a[2] = 1;
    a[3] = 1;
    a[4] = 2;

    for (int k = 3; k + 2 <= n; ++k) {
        const u64 u = (k % 2 == 0) ? 1 : 2;
        const u64 left = mul_mod(a[static_cast<std::size_t>(k)], a[static_cast<std::size_t>(k)]);
        const u64 right = mul_mod(u, mul_mod(a[static_cast<std::size_t>(k + 1)], a[static_cast<std::size_t>(k - 1)]));
        const u64 denominator = a[static_cast<std::size_t>(k - 2)];
        assert(denominator != 0);
        a[static_cast<std::size_t>(k + 2)] = mul_mod(sub_mod(left, right), pow_mod(denominator, MOD - 2));
    }

    return a[static_cast<std::size_t>(n)];
}

void run_checkpoints() {
    assert(square(square(RHO)) == Field(2));
    for (int n = 1; n <= 200; ++n) {
        assert(sequence_value(static_cast<u64>(n)) == direct_value(n));
    }
    assert(sequence_value(13) == 23'321);
    assert(direct_value(1003) == 231'906'014);
    assert(sequence_value(1003) == 231'906'014);
}

}  // namespace

int main() {
    run_checkpoints();
    std::cout << sequence_value(TARGET_N) << '\n';
    return 0;
}
