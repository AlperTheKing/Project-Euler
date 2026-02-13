#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using boost::multiprecision::cpp_int;

template <std::size_t P>
struct PairInfo {
    long double weight;
    long double p_pos;
    long double p_zero;
    long double p_neg;
    std::array<std::int16_t, P> exp;
};

template <std::size_t P>
struct State {
    long double score;
    long double prob;
    std::array<std::int16_t, P> exp;
};

std::vector<int> collect_primes_for_pairs() {
    std::vector<int> primes;
    auto add_factor = [&](int x) {
        int t = x;
        for (int p = 2; p * p <= t; ++p) {
            if (t % p != 0) {
                continue;
            }
            primes.push_back(p);
            while (t % p == 0) {
                t /= p;
            }
        }
        if (t > 1) {
            primes.push_back(t);
        }
    };

    for (int k = 25; k <= 49; ++k) {
        add_factor(k);
        add_factor(100 - k);
    }

    std::sort(primes.begin(), primes.end());
    primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
    return primes;
}

template <std::size_t P>
std::array<std::int16_t, P> exponent_vector(int num, int den, const std::array<int, P>& primes) {
    std::array<std::int16_t, P> e{};

    int a = num;
    int b = den;
    for (std::size_t i = 0; i < P; ++i) {
        const int p = primes[i];
        int cn = 0;
        int cd = 0;
        while (a % p == 0) {
            a /= p;
            ++cn;
        }
        while (b % p == 0) {
            b /= p;
            ++cd;
        }
        e[i] = static_cast<std::int16_t>(2 * (cn - cd));
    }

    return e;
}

template <std::size_t P>
std::vector<PairInfo<P>> build_pairs(const std::array<int, P>& primes) {
    std::vector<PairInfo<P>> pairs;
    pairs.reserve(25);

    for (int k = 25; k <= 49; ++k) {
        const long double q = static_cast<long double>(k) / 100.0L;
        const long double w = 2.0L * std::log((1.0L - q) / q);

        PairInfo<P> info{};
        info.weight = w;
        info.p_pos = (1.0L - q) * (1.0L - q);
        info.p_zero = 2.0L * q * (1.0L - q);
        info.p_neg = q * q;
        info.exp = exponent_vector<P>(100 - k, k, primes);
        pairs.push_back(info);
    }

    return pairs;
}

template <std::size_t P>
std::vector<State<P>> enumerate_states(const std::vector<PairInfo<P>>& pairs, int begin, int end) {
    std::vector<State<P>> states;
    states.push_back(State<P>{0.0L, 1.0L, {}});

    for (int i = begin; i < end; ++i) {
        const auto& pair = pairs[i];
        std::vector<State<P>> next;
        next.reserve(states.size() * 3);

        for (const auto& st : states) {
            State<P> pos = st;
            pos.score += pair.weight;
            pos.prob *= pair.p_pos;
            for (std::size_t j = 0; j < P; ++j) {
                pos.exp[j] = static_cast<std::int16_t>(pos.exp[j] + pair.exp[j]);
            }
            next.push_back(pos);

            State<P> zero = st;
            zero.prob *= pair.p_zero;
            next.push_back(zero);

            State<P> neg = st;
            neg.score -= pair.weight;
            neg.prob *= pair.p_neg;
            for (std::size_t j = 0; j < P; ++j) {
                neg.exp[j] = static_cast<std::int16_t>(neg.exp[j] - pair.exp[j]);
            }
            next.push_back(neg);
        }

        states.swap(next);
    }

    return states;
}

template <std::size_t P>
int compare_exact_product(const std::array<std::int16_t, P>& a,
                         const std::array<std::int16_t, P>& b,
                         const std::array<int, P>& primes) {
    std::array<int, P> c{};
    bool all_zero = true;
    for (std::size_t i = 0; i < P; ++i) {
        c[i] = static_cast<int>(a[i]) + static_cast<int>(b[i]);
        if (c[i] != 0) {
            all_zero = false;
        }
    }
    if (all_zero) {
        return 0;
    }

    cpp_int num = 1;
    cpp_int den = 1;

    for (std::size_t i = 0; i < P; ++i) {
        const int e = c[i];
        if (e > 0) {
            for (int t = 0; t < e; ++t) {
                num *= primes[i];
            }
        } else if (e < 0) {
            for (int t = 0; t < -e; ++t) {
                den *= primes[i];
            }
        }
    }

    if (num > den) {
        return 1;
    }
    if (num < den) {
        return -1;
    }
    return 0;
}

template <std::size_t P>
long double solve_pairs(const std::vector<PairInfo<P>>& pairs, const std::array<int, P>& primes) {
    const int split = static_cast<int>(pairs.size()) / 2;

    auto left = enumerate_states(pairs, 0, split);
    auto right = enumerate_states(pairs, split, static_cast<int>(pairs.size()));

    std::sort(right.begin(), right.end(), [](const State<P>& a, const State<P>& b) {
        return a.score < b.score;
    });

    std::vector<long double> right_scores(right.size());
    std::vector<long double> right_prefix(right.size() + 1, 0.0L);
    for (std::size_t i = 0; i < right.size(); ++i) {
        right_scores[i] = right[i].score;
        right_prefix[i + 1] = right_prefix[i] + right[i].prob;
    }
    const long double total_right = right_prefix.back();

    constexpr long double eps = 1e-14L;
    long double answer = 0.0L;

    for (const auto& ls : left) {
        const long double threshold = -ls.score;

        const auto it_lo = std::upper_bound(right_scores.begin(), right_scores.end(), threshold - eps);
        const auto it_hi = std::upper_bound(right_scores.begin(), right_scores.end(), threshold + eps);
        const std::size_t lo = static_cast<std::size_t>(it_lo - right_scores.begin());
        const std::size_t hi = static_cast<std::size_t>(it_hi - right_scores.begin());

        long double p_correct = total_right - right_prefix[hi];

        for (std::size_t i = lo; i < hi; ++i) {
            const int cmp = compare_exact_product(ls.exp, right[i].exp, primes);
            if (cmp > 0) {
                p_correct += right[i].prob;
            } else if (cmp == 0) {
                p_correct += 0.5L * right[i].prob;
            }
        }

        answer += ls.prob * p_correct;
    }

    return answer;
}

long double brute_accuracy(const std::vector<long double>& p) {
    const int n = static_cast<int>(p.size());
    assert(n <= 20);

    long double acc = 0.0L;
    for (int mask = 0; mask < (1 << n); ++mask) {
        long double prob = 1.0L;
        long double score = 0.0L;
        for (int i = 0; i < n; ++i) {
            const bool say_head = ((mask >> i) & 1) != 0;
            if (say_head) {
                prob *= (1.0L - p[i]);
                score += std::log((1.0L - p[i]) / p[i]);
            } else {
                prob *= p[i];
                score -= std::log((1.0L - p[i]) / p[i]);
            }
        }
        if (score > 0) {
            acc += prob;
        } else if (score == 0) {
            acc += 0.5L * prob;
        }
    }
    return acc;
}

long double solve_problem() {
    const std::vector<int> prime_vec = collect_primes_for_pairs();
    std::array<int, 21> primes{};
    assert(prime_vec.size() == primes.size());
    for (std::size_t i = 0; i < primes.size(); ++i) {
        primes[i] = prime_vec[i];
    }

    const auto pairs = build_pairs(primes);
    return solve_pairs(pairs, primes);
}

void validate() {
    {
        const std::vector<long double> p = {0.2L, 0.4L, 0.6L, 0.8L};
        const long double x = brute_accuracy(p);
        assert(std::fabsl(x - 0.832L) < 1e-15L);
    }

    {
        const std::vector<long double> p = {0.25L, 0.26L, 0.50L, 0.74L, 0.75L};
        const long double brute = brute_accuracy(p);

        const std::array<int, 21> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
                                            37, 41, 43, 47, 53, 59, 61, 67, 71, 73};

        std::vector<PairInfo<21>> pairs;
        for (int k : {25, 26}) {
            const long double q = static_cast<long double>(k) / 100.0L;
            PairInfo<21> info{};
            info.weight = 2.0L * std::log((1.0L - q) / q);
            info.p_pos = (1.0L - q) * (1.0L - q);
            info.p_zero = 2.0L * q * (1.0L - q);
            info.p_neg = q * q;
            info.exp = exponent_vector<21>(100 - k, k, primes);
            pairs.push_back(info);
        }

        const long double mitm = solve_pairs(pairs, primes);
        assert(std::fabsl(brute - mitm) < 1e-15L);
    }
}

}  // namespace

int main() {
    validate();
    std::cout << std::fixed << std::setprecision(10) << solve_problem() << '\n';
    return 0;
}
