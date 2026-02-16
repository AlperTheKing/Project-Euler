#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using mpf = boost::multiprecision::cpp_dec_float_100;

constexpr int kDefaultM = 5000000;
constexpr int kMinE = -180;
constexpr int kMaxE = 90;
constexpr i64 kM = 100000000000000000LL;

struct Options {
    int m = kDefaultM;
    bool run_checkpoints = true;
};

bool parse_nonnegative_int(const std::string& text, int& out) {
    if (text.empty()) return false;
    std::uint64_t value = 0;
    for (char ch : text) {
        if (ch < '0' || ch > '9') return false;
        value = value * 10ULL + static_cast<std::uint64_t>(ch - '0');
        if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            return false;
        }
    }
    out = static_cast<int>(value);
    return true;
}

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    return parse_nonnegative_int(arg.substr(prefix.size()), value);
}

bool parse_arguments(int argc, char** argv, Options& options) {
    bool seen_positional_m = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        int parsed = 0;
        if (parse_int_after_prefix(arg, "--m=", parsed)) {
            options.m = parsed;
            continue;
        }

        if (!seen_positional_m && parse_nonnegative_int(arg, parsed)) {
            options.m = parsed;
            seen_positional_m = true;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.m <= 0) {
        std::cerr << "m must be positive.\n";
        return false;
    }

    return true;
}

mpf cbrt_mpf(const mpf& x) {
    if (x == 0) return mpf(0);
    mpf y = exp(log(x) / 3);
    for (int it = 0; it < 30; ++it) {
        y = (2 * y + x / (y * y)) / 3;
    }
    return y;
}

struct RootTable {
    std::vector<i64> u1;
    std::vector<i64> u2;
    std::vector<i64> u3;
};

RootTable build_root_table() {
    RootTable table;
    const int len = kMaxE - kMinE + 1;
    table.u1.resize(static_cast<std::size_t>(len));
    table.u2.resize(static_cast<std::size_t>(len));
    table.u3.resize(static_cast<std::size_t>(len));

    const mpf p = (mpf(29) - 3 * sqrt(mpf(93))) / 2;
    const mpf q = mpf(29) - p;
    const mpf r = (mpf(1) + cbrt_mpf(p) + cbrt_mpf(q)) / 3;

    mpf z = 1;
    for (int i = 0; i < -kMinE; ++i) z /= r;

    for (int i = 0; i < len; ++i) {
        mpf x = z;
        const i64 y1 = x.convert_to<i64>();
        table.u1[static_cast<std::size_t>(i)] = y1;

        x = mpf(kM) * (x - mpf(y1));
        const i64 y2 = x.convert_to<i64>();
        table.u2[static_cast<std::size_t>(i)] = y2;

        x = mpf(kM) * (x - mpf(y2));
        const i64 y3 = x.convert_to<i64>();
        table.u3[static_cast<std::size_t>(i)] = y3;

        z *= r;
    }

    return table;
}

inline bool leq_lex(const i64 a1, const i64 a2, const i64 a3,
                    const i64 b1, const i64 b2, const i64 b3) {
    if (a1 != b1) return a1 < b1;
    if (a2 != b2) return a2 < b2;
    return a3 <= b3;
}

struct SolveResult {
    u64 s_m = 0;
    u64 s_10 = 0;
    u64 s_1000 = 0;
};

SolveResult solve_with_table(const int m, const RootTable& table) {
    SolveResult out;

    const int len = static_cast<int>(table.u1.size());
    int i0 = -kMinE;

    for (int j = 1; j <= m; ++j) {
        const i64 n = static_cast<i64>(j) * static_cast<i64>(j);

        if (n == 1) {
            out.s_m += 1;
            if (j == 10) out.s_10 = out.s_m;
            if (j == 1000) out.s_1000 = out.s_m;
            continue;
        }

        while (i0 < len && table.u1[static_cast<std::size_t>(i0)] < n) ++i0;
        if (i0 >= len) {
            throw std::runtime_error("Exponent window exceeded; increase [minE,maxE].");
        }

        int i = i0 - 1;
        i64 rep = 1;

        i64 n1 = n - table.u1[static_cast<std::size_t>(i)] - 1;
        i64 n2 = kM - table.u2[static_cast<std::size_t>(i)] - 1;
        i64 n3 = kM - table.u3[static_cast<std::size_t>(i)];

        while (true) {
            i -= 3;
            while (i >= 0 &&
                   !leq_lex(table.u1[static_cast<std::size_t>(i)],
                            table.u2[static_cast<std::size_t>(i)],
                            table.u3[static_cast<std::size_t>(i)],
                            n1, n2, n3)) {
                --i;
            }
            if (i < 0) {
                throw std::runtime_error("Exponent window underflow; decrease minE.");
            }

            ++rep;

            n1 -= table.u1[static_cast<std::size_t>(i)];
            n2 -= table.u2[static_cast<std::size_t>(i)];
            n3 -= table.u3[static_cast<std::size_t>(i)];

            if (n3 < 0) {
                n3 += kM;
                --n2;
            }
            if (n2 < 0) {
                n2 += kM;
                --n1;
            }

            if (n1 == 0 && n2 == 0 && n3 < 1000) {
                out.s_m += static_cast<u64>(rep);
                break;
            }
        }

        if (j == 10) out.s_10 = out.s_m;
        if (j == 1000) out.s_1000 = out.s_m;
    }

    return out;
}

bool run_validations(const RootTable& table) {
    const SolveResult sample = solve_with_table(1000, table);
    if (sample.s_10 != 61ULL) {
        std::cerr << "Validation failed: S(10) should be 61, got " << sample.s_10 << '\n';
        return false;
    }
    if (sample.s_1000 != 19403ULL) {
        std::cerr << "Validation failed: S(1000) should be 19403, got " << sample.s_1000 << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) return 1;

    try {
        const RootTable table = build_root_table();

        if (options.run_checkpoints && !run_validations(table)) {
            return 1;
        }

        const SolveResult result = solve_with_table(options.m, table);
        std::cout << result.s_m << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
