#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

struct Fraction {
    std::uint64_t num = 0;
    std::uint64_t den = 1;
};

Fraction reduce(Fraction f) {
    std::uint64_t g = std::gcd(f.num, f.den);
    f.num /= g;
    f.den /= g;
    return f;
}

bool less_frac(const Fraction& a, const Fraction& b) {
    return static_cast<__int128>(a.num) * b.den < static_cast<__int128>(b.num) * a.den;
}

bool greater_frac(const Fraction& a, const Fraction& b) {
    return static_cast<__int128>(a.num) * b.den > static_cast<__int128>(b.num) * a.den;
}

Fraction next_transition(const Fraction& r, int base_m, int margin, int& used_m) {
    auto compute_best = [&](int M, int& end_i, bool& overflow, int& best_i) {
        std::vector<std::uint64_t> p(static_cast<std::size_t>(M) + 1, 0);
        p[1] = 1;
        p[2] = 2;
        int d = 1;
        Fraction best{0, 1};
        bool best_set = false;
        best_i = -1;
        end_i = 2;
        overflow = false;

        for (int i = 3; i <= M; ++i) {
            while (d + 1 <= i - 1) {
                const std::uint64_t left = p[static_cast<std::size_t>(i - 1)];
                const std::uint64_t right = p[static_cast<std::size_t>(i - (d + 1))];
                if (static_cast<__int128>(left) * r.den <=
                    static_cast<__int128>(r.num) * right) {
                    ++d;
                } else {
                    break;
                }
            }

            if (i - (d + 1) >= 1) {
                Fraction cand{p[static_cast<std::size_t>(i - 1)],
                              p[static_cast<std::size_t>(i - (d + 1))]};
                if (greater_frac(cand, r) && (!best_set || less_frac(cand, best))) {
                    best = cand;
                    best_set = true;
                    best_i = i;
                }
            }

            const std::uint64_t a = p[static_cast<std::size_t>(i - 1)];
            const std::uint64_t b = p[static_cast<std::size_t>(i - d)];
            if (a > std::numeric_limits<std::uint64_t>::max() - b) {
                overflow = true;
                end_i = i - 1;
                break;
            }
            p[static_cast<std::size_t>(i)] = a + b;
            end_i = i;
        }
        if (!best_set) {
            best = Fraction{0, 1};
        }
        return best;
    };

    const int max_m = 1 << 20;
    Fraction prev_best{0, 1};
    bool have_prev = false;
    int M = base_m;
    while (true) {
        int end_i = 2;
        bool overflow = false;
        int best_i = -1;
        Fraction best = compute_best(M, end_i, overflow, best_i);

        if (overflow || M >= max_m) {
            used_m = end_i;
            return reduce(best);
        }
        if (have_prev && best.num == prev_best.num && best.den == prev_best.den &&
            best_i < end_i - margin) {
            used_m = end_i;
            return reduce(best);
        }
        prev_best = best;
        have_prev = true;
        M *= 2;
    }
}

std::vector<std::uint64_t> compute_p_brut(const Fraction& r, int nmax) {
    std::vector<int> g(nmax + 1, 0);
    g[0] = 1'000'000'000;
    g[1] = 1;
    for (int n = 2; n <= nmax; ++n) {
        for (int k = 1; k <= n; ++k) {
            if (static_cast<__int128>(k) * r.num <
                static_cast<__int128>(g[n - k]) * r.den) {
                g[n] = k;
                break;
            }
        }
        if (g[n] == 0) g[n] = n;
    }
    std::vector<std::uint64_t> p;
    for (int n = 1; n <= nmax; ++n) {
        if (g[n] == n) p.push_back(static_cast<std::uint64_t>(n));
    }
    return p;
}

std::vector<std::uint64_t> compute_p_recur(const Fraction& r, int nmax) {
    std::vector<std::uint64_t> p;
    p.reserve(64);
    p.push_back(1);
    p.push_back(2);
    int d = 1;
    std::vector<std::uint64_t> work;
    work.reserve(256);
    work.push_back(0);
    work.push_back(1);
    work.push_back(2);

    while (true) {
        int i = static_cast<int>(work.size());
        while (d + 1 <= i - 1) {
            std::uint64_t left = work[static_cast<std::size_t>(i - 1)];
            std::uint64_t right = work[static_cast<std::size_t>(i - (d + 1))];
            if (static_cast<__int128>(left) * r.den <=
                static_cast<__int128>(r.num) * right) {
                ++d;
            } else {
                break;
            }
        }
        std::uint64_t next = work[static_cast<std::size_t>(i - 1)] +
                             work[static_cast<std::size_t>(i - d)];
        if (next > static_cast<std::uint64_t>(nmax)) break;
        work.push_back(next);
        p.push_back(next);
    }
    return p;
}

void validate() {
    const int base_m = 512;
    const int margin = 32;
    int used_m = 0;
    Fraction r{1, 1};
    Fraction t2{0, 1};
    Fraction t22{0, 1};
    for (int i = 2; i <= 22; ++i) {
        r = next_transition(r, base_m, margin, used_m);
        if (i == 2) t2 = r;
        if (i == 22) t22 = r;
    }
    if (!(t2.num == 2 && t2.den == 1)) {
        std::cerr << "Validation failed for T(2).\n";
        std::exit(1);
    }
    if (!(t22.num == 145 && t22.den == 23)) {
        std::cerr << "Validation failed for T(22).\n";
        std::exit(1);
    }

    const Fraction checks[] = {
        Fraction{1, 1}, Fraction{2, 1}, Fraction{5, 2},
        Fraction{3, 1}, Fraction{7, 2}, Fraction{11, 3},
    };
    for (const auto& rr : checks) {
        auto p_brut = compute_p_brut(rr, 200);
        auto p_rec = compute_p_recur(rr, 200);
        if (p_brut != p_rec) {
            std::cerr << "Validation failed for P-positions at r="
                      << rr.num << "/" << rr.den << ".\n";
            std::exit(1);
        }
    }
}

}

int main(int argc, char** argv) {
    int target = 123456;
    int base_m = 6000;
    bool do_validate = true;
    bool target_set = false;
    bool base_set = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-validate") {
            do_validate = false;
            continue;
        }
        char* end = nullptr;
        long val = std::strtol(arg.c_str(), &end, 10);
        if (end && *end == '\0') {
            if (!target_set) {
                target = std::max(1L, val);
                target_set = true;
            } else if (!base_set) {
                base_m = std::max(128L, val);
                base_set = true;
            }
        }
    }

    if (do_validate) {
        validate();
    }

    Fraction ans{1, 1};
    if (target == 1) {
        ans = Fraction{1, 1};
    } else {
        Fraction r{1, 1};
        const int margin = 64;
        int used_m = base_m;
        for (int i = 2; i <= target; ++i) {
            ans = next_transition(r, base_m, margin, used_m);
            r = ans;
        }
    }

    long double value = static_cast<long double>(ans.num) / static_cast<long double>(ans.den);
    std::cout << std::fixed << std::setprecision(10) << value << '\n';
    return 0;
}
