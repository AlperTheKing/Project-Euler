#include <boost/rational.hpp>

#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
#include <thread>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

using i64 = std::int64_t;
using Fraction = boost::rational<i64>;
using cpp_int = boost::multiprecision::cpp_int;

namespace {

const std::string kDigits = "123456789";
constexpr int kTotalSplits = 8;
constexpr int kOperatorCount = 5;
constexpr int kOperators = kOperatorCount;
constexpr int kTotalExpressions = 390625;  // 5^8

char op_from_code(int code) {
    if (code == 0) return '+';
    if (code == 1) return '-';
    if (code == 2) return '*';
    if (code == 3) return '/';
    return '#'; // concat marker
}

Fraction apply(Fraction lhs, const Fraction& rhs, char op) {
    switch (op) {
        case '+':
            return lhs + rhs;
        case '-':
            return lhs - rhs;
        case '*':
            return lhs * rhs;
        case '/':
            if (rhs == 0) {
                throw std::overflow_error("division by zero");
            }
            return lhs / rhs;
    }
    return lhs;
}

std::vector<Fraction> all_results(const std::vector<i64>& values, const std::vector<char>& ops) {
    const int n = static_cast<int>(values.size());
    std::vector<std::vector<std::vector<Fraction>>> dp(static_cast<size_t>(n),
                                                      std::vector<std::vector<Fraction>>(static_cast<size_t>(n)));

    for (int i = 0; i < n; ++i) {
        dp[static_cast<size_t>(i)][static_cast<size_t>(i)] = {Fraction(values[static_cast<size_t>(i)], 1)};
    }

    for (int len = 2; len <= n; ++len) {
        for (int l = 0; l + len <= n; ++l) {
            const int r = l + len - 1;
            std::set<Fraction> seen;
            std::vector<Fraction> combined;

            for (int k = l; k < r; ++k) {
                const auto& left = dp[static_cast<size_t>(l)][static_cast<size_t>(k)];
                const auto& right = dp[static_cast<size_t>(k + 1)][static_cast<size_t>(r)];
                const char op = ops[static_cast<size_t>(k)];

                for (const Fraction& a : left) {
                    for (const Fraction& b : right) {
                        if (op == '/' && b == 0) {
                            continue;
                        }
                        Fraction v = apply(a, b, op);
                        if (seen.insert(v).second) {
                            combined.push_back(v);
                        }
                    }
                }
            }
            dp[static_cast<size_t>(l)][static_cast<size_t>(r)] = std::move(combined);
        }
    }

    return dp[0][static_cast<size_t>(n - 1)];
}

void collect_for_range(int start, int end, std::vector<i64>& out) {
    for (int code = start; code < end; ++code) {
        int x = code;
        std::vector<i64> values;
        values.reserve(9);
        std::vector<char> ops;
        ops.reserve(8);

        i64 current = kDigits[0] - '0';
        for (std::size_t i = 0; i < kDigits.size() - 1; ++i) {
            const int choice = x % kOperators;
            x /= kOperators;
            const i64 next_digit = static_cast<i64>(kDigits[i + 1] - '0');

            const char op = op_from_code(choice);
            if (op == '#') {
                current = current * 10 + next_digit;
            } else {
                values.push_back(current);
                ops.push_back(op);
                current = next_digit;
            }
        }
        values.push_back(current);

        const auto results = all_results(values, ops);
        for (const Fraction& v : results) {
            if (v.denominator() == 1 && v.numerator() > 0) {
                out.push_back(v.numerator());
            }
        }
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

}  // namespace

int main() {
    const unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    const int base = kTotalExpressions / static_cast<int>(threads);
    const int rem = kTotalExpressions % static_cast<int>(threads);

    std::vector<std::vector<i64>> thread_results(threads);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    int start = 0;
    for (unsigned t = 0; t < threads; ++t) {
        const int extra = t < static_cast<unsigned>(rem) ? 1 : 0;
        const int end = start + base + extra;
        workers.emplace_back(collect_for_range, start, end, std::ref(thread_results[static_cast<size_t>(t)]));
        start = end;
    }

    for (auto& worker : workers) {
        worker.join();
    }

    std::vector<i64> all;
    for (const auto& vec : thread_results) {
        all.insert(all.end(), vec.begin(), vec.end());
    }
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());

    cpp_int sum = 0;
    for (const i64 v : all) {
        sum += v;
    }

    std::cout << "Total number of unique reachable numbers: " << all.size() << '\n';
    std::cout << "Sum of all unique reachable integers: " << sum << '\n';
    return 0;
}
