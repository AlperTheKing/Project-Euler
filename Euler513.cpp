#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>

namespace {

using i64 = long long;

constexpr bool kClosed = true;
constexpr bool kOpen = false;

struct Options {
    int n = 100000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;

    unsigned long long parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') return false;
        const unsigned long long d = static_cast<unsigned long long>(ch - '0');
        if (parsed > (std::numeric_limits<unsigned long long>::max() - d) / 10ULL) return false;
        parsed = parsed * 10ULL + d;
    }
    if (parsed > static_cast<unsigned long long>(std::numeric_limits<int>::max())) return false;
    value = static_cast<int>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) continue;
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

inline i64 floor_div(i64 a, i64 b) {
    assert(b > 0);
    if (a >= 0) return a / b;
    return -((-a + b - 1) / b);
}

i64 isqrt_floor(i64 n) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

i64 points_in_trapezoid(i64 slope,
                        i64 intercept,
                        i64 denominator,
                        i64 lower_domain,
                        i64 upper_domain,
                        bool boundary) {
    i64 result = 0;
    while (true) {
        assert(denominator > 0);
        if (std::llabs(upper_domain - lower_domain) <= 8) {
            i64 s = 0;
            const i64 adjustment = boundary ? 0 : 1;
            if (upper_domain > lower_domain) {
                for (i64 x = lower_domain + 1; x <= upper_domain; ++x) {
                    s += floor_div(slope * x + intercept - adjustment, denominator);
                }
            } else {
                for (i64 x = upper_domain + 1; x <= lower_domain; ++x) {
                    s += floor_div(slope * x + intercept - adjustment, denominator);
                }
                s = -s;
            }
            result += s;
            break;
        }

        result += (upper_domain - lower_domain) * floor_div(intercept, denominator);
        intercept = ((intercept % denominator) + denominator) % denominator;
        result += ((upper_domain - lower_domain) * (upper_domain + lower_domain + 1) / 2) *
                  floor_div(slope, denominator);
        slope = ((slope % denominator) + denominator) % denominator;

        if (slope != 0) {
            const i64 upper_value = floor_div(slope * upper_domain + intercept, denominator);
            const i64 lower_value = floor_div(slope * lower_domain + intercept, denominator);
            result += upper_domain * upper_value - lower_domain * lower_value;
            lower_domain = upper_value;
            upper_domain = lower_value;
            std::swap(slope, denominator);
            intercept = -intercept;
            boundary = !boundary;
        } else {
            if (intercept == 0 && !boundary) result -= upper_domain - lower_domain;
            break;
        }
    }
    return result;
}

i64 points_in_trapezoid_mod2(i64 slope,
                             i64 intercept,
                             i64 denominator,
                             i64 lower_domain,
                             i64 upper_domain,
                             bool boundary,
                             i64 x_residue,
                             i64 y_residue) {
    if ((y_residue & 1LL) != 0) intercept += denominator;
    if ((x_residue & 1LL) != 0) {
        intercept -= slope;
        ++lower_domain;
        ++upper_domain;
    }
    return points_in_trapezoid(2 * slope, intercept, 2 * denominator,
                               floor_div(lower_domain, 2), floor_div(upper_domain, 2), boundary);
}

i64 f_value(i64 n) {
    i64 result = 0;
    const i64 three_halves_n = n + n / 2;
    const i64 root = isqrt_floor(three_halves_n);
    const int ij[3][2] = {{0, 1}, {1, 0}, {1, 1}};

    for (int tc = 0; tc < 3; ++tc) {
        const i64 i = ij[tc][0];
        const i64 j = ij[tc][1];

        i64 max_t = 1;
        for (i64 s = 2; s < root; ++s) {
            if (3 * (max_t + 1) * (max_t + 1) <= s * s) ++max_t;
            if (i == j || (s & 1LL) == 0) {
                for (i64 t = ((s - 1) & 1LL) + 1; t <= max_t; t += 2) {
                    const i64 v_mid = t * n / ((s - t) * (s + t));
                    const i64 v_max = n * (s + t) / (s * s + 2 * s * t - t * t);
                    result += points_in_trapezoid_mod2(s, 0, t, 0, v_mid, kClosed, j, i);
                    result += points_in_trapezoid_mod2(t, n, s, v_mid, v_max, kClosed, j, i);
                    result -= points_in_trapezoid_mod2(s + 3 * t, 0, s + t, 0, v_max, kOpen, j, i);
                }
            }
        }

        const i64 u_max = three_halves_n / root;
        for (i64 u = 1 + ((i + 1) & 1LL); u <= u_max; u += 2) {
            const i64 v_max_outer = std::min(n / 2, u - 1);
            for (i64 v = 1 + ((j + 1) & 1LL); v <= v_max_outer; v += 2) {
                const i64 mid_s = (v + n) / u;
                const int residue_count = (i == j ? 2 : 1);
                for (int sr = 0; sr < residue_count; ++sr) {
                    const i64 s_residue = (i == j ? sr : 0);

                    i64 min_s = 0, max_s = -1;
                    i64 slope0 = 0, slope1 = 1;
                    if (u * u < 3 * v * v) {
                        min_s = root;
                        max_s = n * (3 * v - u) / (2 * u * v + v * v - u * u);
                        slope0 = u - v;
                        slope1 = 3 * v - u;
                    } else {
                        min_s = std::max(root, floor_div(u + v - 1, v));
                        max_s = n * u / ((u - v) * (u + v));
                        slope0 = v;
                        slope1 = u;
                    }

                    if (max_s >= min_s) {
                        result += points_in_trapezoid_mod2(slope0, 0, slope1,
                                                           min_s - 1, max_s, kClosed,
                                                           s_residue, s_residue);
                        if (mid_s < max_s) {
                            result -= points_in_trapezoid_mod2(u, -n, v,
                                                               std::max(mid_s, min_s - 1), max_s,
                                                               kOpen, s_residue, s_residue);
                        }
                    }
                }
            }
        }
    }
    return result;
}

std::unordered_map<i64, i64> memo_F;

i64 F(i64 n) {
    if (n <= 0) return 0;
    const auto it = memo_F.find(n);
    if (it != memo_F.end()) return it->second;

    i64 result = f_value(n);
    i64 k = 3;
    i64 n_over_k = n / k;
    while (k <= n_over_k) {
        result -= F(n_over_k);
        k += 2;
        n_over_k = n / k;
    }

    i64 min_k = n / (n_over_k + 1);
    while (n_over_k) {
        const i64 max_k = n / n_over_k;
        const i64 left = (min_k + 1) + (min_k & 1LL);
        const i64 right = max_k - ((max_k + 1) & 1LL);
        i64 count = 0;
        if (right >= left) count = (right - left) / 2 + 1;
        result -= F(n_over_k) * count;
        --n_over_k;
        min_k = max_k;
    }

    memo_F.emplace(n, result);
    return result;
}

bool run_checkpoints() {
    if (F(10) != 3) {
        std::cerr << "Validation failed: F(10)\n";
        return false;
    }
    if (F(50) != 165) {
        std::cerr << "Validation failed: F(50)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    memo_F.reserve(1 << 15);
    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    std::cout << F(options.n) << '\n';
    return 0;
}
