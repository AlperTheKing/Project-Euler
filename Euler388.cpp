#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;

struct Options {
    i64 n = 10000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    i64 parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

std::string to_string_i128(i128 value) {
    if (value == 0) {
        return "0";
    }
    bool neg = false;
    if (value < 0) {
        neg = true;
        value = -value;
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (neg) {
        s.push_back('-');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::string first9_last9_token(const std::string& s) {
    if (s.size() <= 18U) {
        return s;
    }
    return s.substr(0, 9) + s.substr(s.size() - 9);
}

std::vector<int> mu_sieve(const int limit) {
    std::vector<int> primes;
    primes.reserve(limit / 10);
    std::vector<int> lp(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> mu(static_cast<std::size_t>(limit + 1), 0);
    mu[1] = 1;

    for (int i = 2; i <= limit; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0) {
            lp[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const i64 x = 1LL * i * p;
            if (x > limit) {
                break;
            }
            lp[static_cast<std::size_t>(x)] = p;
            if (p == lp[static_cast<std::size_t>(i)]) {
                mu[static_cast<std::size_t>(x)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(x)] = -mu[static_cast<std::size_t>(i)];
        }
    }
    return mu;
}

class MertensPrefix {
   public:
    explicit MertensPrefix(const i64 n) {
        const i64 estimated = static_cast<i64>(std::pow(static_cast<long double>(n), 2.0L / 3.0L));
        sieve_limit_ = static_cast<int>(std::max<i64>(1000000LL, estimated + 1000));

        const std::vector<int> mu = mu_sieve(sieve_limit_);
        prefix_.assign(static_cast<std::size_t>(sieve_limit_ + 1), 0LL);
        for (int i = 1; i <= sieve_limit_; ++i) {
            prefix_[static_cast<std::size_t>(i)] =
                prefix_[static_cast<std::size_t>(i - 1)] + static_cast<i64>(mu[static_cast<std::size_t>(i)]);
        }
        cache_.reserve(1 << 20);
    }

    i64 M(const i64 n) {
        if (n <= sieve_limit_) {
            return prefix_[static_cast<std::size_t>(n)];
        }
        const auto it = cache_.find(n);
        if (it != cache_.end()) {
            return it->second;
        }

        i64 ans = 1;
        for (i64 l = 2; l <= n;) {
            const i64 q = n / l;
            const i64 r = n / q;
            ans -= (r - l + 1) * M(q);
            l = r + 1;
        }

        cache_[n] = ans;
        return ans;
    }

   private:
    int sieve_limit_{0};
    std::vector<i64> prefix_;
    std::unordered_map<i64, i64> cache_;
};

i128 distinct_lines(const i64 n) {
    MertensPrefix mertens(n);
    i128 answer = 0;

    for (i64 l = 1; l <= n;) {
        const i64 q = n / l;
        const i64 r = n / q;
        const i64 mu_block = mertens.M(r) - mertens.M(l - 1);

        const i128 qq = static_cast<i128>(q);
        const i128 f = (qq + 1) * (qq + 1) * (qq + 1) - 1;
        answer += f * static_cast<i128>(mu_block);

        l = r + 1;
    }

    return answer;
}

i128 brute_distinct_lines(const int n) {
    i128 count = 0;
    for (int a = 0; a <= n; ++a) {
        for (int b = 0; b <= n; ++b) {
            for (int c = 0; c <= n; ++c) {
                if (a == 0 && b == 0 && c == 0) {
                    continue;
                }
                const int g = std::gcd(a, std::gcd(b, c));
                if (g == 1) {
                    ++count;
                }
            }
        }
    }
    return count;
}

bool run_checkpoints() {
    if (distinct_lines(1) != 7) {
        std::cerr << "Checkpoint failed: D(1)\n";
        return false;
    }

    const int brute_n = 40;
    const i128 brute = brute_distinct_lines(brute_n);
    const i128 fast = distinct_lines(brute_n);
    if (fast != brute) {
        std::cerr << "Checkpoint failed against brute force at N=" << brute_n << '\n';
        return false;
    }

    const i128 known = distinct_lines(1000000);
    if (to_string_i128(known) != "831909254469114121") {
        std::cerr << "Checkpoint failed: D(10^6), got " << to_string_i128(known) << '\n';
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

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const i128 answer = distinct_lines(options.n);
    const std::string full = to_string_i128(answer);
    std::cout << first9_last9_token(full) << '\n';
    return 0;
}
