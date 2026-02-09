#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr int MOD = 100000007;
constexpr int CELLS = 9;
constexpr int FULL_MASK = (1 << CELLS) - 1;
constexpr int STATE_COUNT = 1 << CELLS;

int add_mod(const int a, const int b) {
    int x = a + b;
    if (x >= MOD) {
        x -= MOD;
    }
    return x;
}

int sub_mod(const int a, const int b) {
    int x = a - b;
    if (x < 0) {
        x += MOD;
    }
    return x;
}

int mul_mod(const int a, const int b) {
    return static_cast<int>((static_cast<i64>(a) * static_cast<i64>(b)) % MOD);
}

int pow_mod(int base, u64 exp) {
    int result = 1;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return result;
}

int inv_mod(const int a) {
    return pow_mod(a, static_cast<u64>(MOD - 2));
}

void dfs_layer_fill(const int filled_mask, const int next_mask, std::array<int, STATE_COUNT>& ways) {
    if (filled_mask == FULL_MASK) {
        ++ways[static_cast<std::size_t>(next_mask)];
        return;
    }

    const int free_bits = (~filled_mask) & FULL_MASK;
    const int bit = __builtin_ctz(static_cast<unsigned int>(free_bits));

    // Place a domino across layers: marks this cell in next_mask.
    dfs_layer_fill(filled_mask | (1 << bit), next_mask | (1 << bit), ways);

    const int r = bit / 3;
    const int c = bit % 3;

    // Place inside the current layer along columns.
    if (c < 2 && ((filled_mask & (1 << (bit + 1))) == 0)) {
        dfs_layer_fill(filled_mask | (1 << bit) | (1 << (bit + 1)), next_mask, ways);
    }

    // Place inside the current layer along rows.
    if (r < 2 && ((filled_mask & (1 << (bit + 3))) == 0)) {
        dfs_layer_fill(filled_mask | (1 << bit) | (1 << (bit + 3)), next_mask, ways);
    }
}

std::array<std::vector<std::pair<int, int>>, STATE_COUNT> build_transitions() {
    std::array<std::vector<std::pair<int, int>>, STATE_COUNT> transitions;

    for (int mask = 0; mask < STATE_COUNT; ++mask) {
        std::array<int, STATE_COUNT> ways{};
        dfs_layer_fill(mask, 0, ways);

        auto& edges = transitions[static_cast<std::size_t>(mask)];
        edges.reserve(64);
        for (int nxt = 0; nxt < STATE_COUNT; ++nxt) {
            const int cnt = ways[static_cast<std::size_t>(nxt)];
            if (cnt > 0) {
                edges.emplace_back(nxt, cnt % MOD);
            }
        }
    }

    return transitions;
}

std::vector<int> generate_sequence(const int terms,
                                   const std::array<std::vector<std::pair<int, int>>, STATE_COUNT>& transitions) {
    std::vector<int> seq(static_cast<std::size_t>(terms + 1), 0);

    std::array<int, STATE_COUNT> current{};
    std::array<int, STATE_COUNT> next{};
    current[0] = 1;
    seq[0] = 1;

    for (int step = 1; step <= terms; ++step) {
        next.fill(0);
        for (int mask = 0; mask < STATE_COUNT; ++mask) {
            const int value = current[static_cast<std::size_t>(mask)];
            if (value == 0) {
                continue;
            }
            for (const auto& [nxt, cnt] : transitions[static_cast<std::size_t>(mask)]) {
                const int add = mul_mod(value, cnt);
                next[static_cast<std::size_t>(nxt)] = add_mod(next[static_cast<std::size_t>(nxt)], add);
            }
        }
        current = next;
        seq[static_cast<std::size_t>(step)] = current[0];
    }

    return seq;
}

std::vector<int> berlekamp_massey(const std::vector<int>& sequence) {
    std::vector<int> C{1};
    std::vector<int> B{1};

    int L = 0;
    int m = 1;
    int b = 1;

    for (int n = 0; n < static_cast<int>(sequence.size()); ++n) {
        int d = sequence[static_cast<std::size_t>(n)];
        for (int i = 1; i <= L; ++i) {
            d = add_mod(d, mul_mod(C[static_cast<std::size_t>(i)], sequence[static_cast<std::size_t>(n - i)]));
        }

        if (d == 0) {
            ++m;
            continue;
        }

        const std::vector<int> T = C;
        const int coef = mul_mod(d, inv_mod(b));

        if (static_cast<int>(C.size()) < static_cast<int>(B.size()) + m) {
            C.resize(static_cast<std::size_t>(static_cast<int>(B.size()) + m), 0);
        }

        for (int i = 0; i < static_cast<int>(B.size()); ++i) {
            const int idx = i + m;
            C[static_cast<std::size_t>(idx)] = sub_mod(C[static_cast<std::size_t>(idx)],
                                                       mul_mod(coef, B[static_cast<std::size_t>(i)]));
        }

        if (2 * L <= n) {
            L = n + 1 - L;
            B = T;
            b = d;
            m = 1;
        } else {
            ++m;
        }
    }

    std::vector<int> recurrence(static_cast<std::size_t>(L), 0);
    for (int i = 1; i <= L; ++i) {
        recurrence[static_cast<std::size_t>(i - 1)] = (MOD - C[static_cast<std::size_t>(i)]) % MOD;
    }
    return recurrence;
}

std::vector<int> combine_polynomials(const std::vector<int>& a, const std::vector<int>& b,
                                     const std::vector<int>& recurrence) {
    const int k = static_cast<int>(recurrence.size());
    std::vector<int> prod(static_cast<std::size_t>(2 * k), 0);

    for (int i = 0; i < k; ++i) {
        if (a[static_cast<std::size_t>(i)] == 0) {
            continue;
        }
        for (int j = 0; j < k; ++j) {
            if (b[static_cast<std::size_t>(j)] == 0) {
                continue;
            }
            const int idx = i + j;
            prod[static_cast<std::size_t>(idx)] = add_mod(
                prod[static_cast<std::size_t>(idx)],
                mul_mod(a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(j)]));
        }
    }

    for (int i = 2 * k - 2; i >= k; --i) {
        const int coef = prod[static_cast<std::size_t>(i)];
        if (coef == 0) {
            continue;
        }
        for (int j = 1; j <= k; ++j) {
            const int idx = i - j;
            prod[static_cast<std::size_t>(idx)] = add_mod(
                prod[static_cast<std::size_t>(idx)], mul_mod(coef, recurrence[static_cast<std::size_t>(j - 1)]));
        }
    }

    prod.resize(static_cast<std::size_t>(k));
    return prod;
}

int linear_recurrence_nth(const std::vector<int>& initial, const std::vector<int>& recurrence,
                          const std::vector<int>& bits_lsb_first) {
    const int k = static_cast<int>(recurrence.size());

    std::vector<int> result(static_cast<std::size_t>(k), 0);
    result[0] = 1;

    std::vector<int> x_poly(static_cast<std::size_t>(k), 0);
    if (k == 1) {
        x_poly[0] = recurrence[0];
    } else {
        x_poly[1] = 1;
    }

    for (const int bit : bits_lsb_first) {
        if (bit != 0) {
            result = combine_polynomials(result, x_poly, recurrence);
        }
        x_poly = combine_polynomials(x_poly, x_poly, recurrence);
    }

    int answer = 0;
    for (int i = 0; i < k; ++i) {
        answer = add_mod(answer, mul_mod(result[static_cast<std::size_t>(i)], initial[static_cast<std::size_t>(i)]));
    }
    return answer;
}

std::vector<int> bits_from_u64(u64 n) {
    std::vector<int> bits;
    bits.reserve(64);
    while (n > 0ULL) {
        bits.push_back(static_cast<int>(n & 1ULL));
        n >>= 1ULL;
    }
    if (bits.empty()) {
        bits.push_back(0);
    }
    return bits;
}

std::vector<int> bits_from_decimal(std::string value) {
    std::vector<int> bits;
    bits.reserve(value.size() * 4);

    while (!(value.size() == 1 && value[0] == '0')) {
        int carry = 0;
        std::string next;
        next.reserve(value.size());

        for (char ch : value) {
            const int cur = carry * 10 + static_cast<int>(ch - '0');
            const int q = cur / 2;
            carry = cur % 2;
            if (!next.empty() || q != 0) {
                next.push_back(static_cast<char>('0' + q));
            }
        }

        bits.push_back(carry);
        value = next.empty() ? "0" : next;
    }

    if (bits.empty()) {
        bits.push_back(0);
    }
    return bits;
}

int nth_term(const u64 n, const std::vector<int>& initial, const std::vector<int>& recurrence) {
    if (n < initial.size()) {
        return initial[static_cast<std::size_t>(n)];
    }
    return linear_recurrence_nth(initial, recurrence, bits_from_u64(n));
}

int nth_term_decimal(const std::string& n_decimal, const std::vector<int>& initial,
                     const std::vector<int>& recurrence) {
    return linear_recurrence_nth(initial, recurrence, bits_from_decimal(n_decimal));
}

bool validate_recurrence(const std::vector<int>& sequence, const std::vector<int>& recurrence) {
    const int k = static_cast<int>(recurrence.size());
    for (int n = k; n < static_cast<int>(sequence.size()); ++n) {
        int expected = 0;
        for (int i = 1; i <= k; ++i) {
            expected = add_mod(expected, mul_mod(recurrence[static_cast<std::size_t>(i - 1)],
                                                 sequence[static_cast<std::size_t>(n - i)]));
        }
        if (expected != sequence[static_cast<std::size_t>(n)]) {
            return false;
        }
    }
    return true;
}

bool run_checkpoints(const std::vector<int>& initial, const std::vector<int>& recurrence) {
    if (nth_term(2ULL, initial, recurrence) != 229) {
        std::cerr << "Checkpoint failed: f(2)\n";
        return false;
    }
    if (nth_term(4ULL, initial, recurrence) != 117805) {
        std::cerr << "Checkpoint failed: f(4)\n";
        return false;
    }
    if (nth_term(10ULL, initial, recurrence) != 96149360) {
        std::cerr << "Checkpoint failed: f(10)\n";
        return false;
    }
    if (nth_term(1000ULL, initial, recurrence) != 24806056) {
        std::cerr << "Checkpoint failed: f(10^3)\n";
        return false;
    }
    if (nth_term(1000000ULL, initial, recurrence) != 30808124) {
        std::cerr << "Checkpoint failed: f(10^6)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    const auto transitions = build_transitions();
    const std::vector<int> sequence = generate_sequence(1300, transitions);
    const std::vector<int> recurrence = berlekamp_massey(sequence);

    if (recurrence.empty()) {
        std::cerr << "Failed to derive recurrence\n";
        return 2;
    }

    if (!validate_recurrence(sequence, recurrence)) {
        std::cerr << "Recurrence validation failed\n";
        return 3;
    }

    std::vector<int> initial(recurrence.size(), 0);
    for (std::size_t i = 0; i < recurrence.size(); ++i) {
        initial[i] = sequence[i];
    }

    if (!skip_checkpoints && !run_checkpoints(initial, recurrence)) {
        return 4;
    }

    std::string exponent = "1";
    exponent.append(10000, '0');  // 10^10000

    const int answer = nth_term_decimal(exponent, initial, recurrence);
    std::cout << answer << '\n';
    return 0;
}
