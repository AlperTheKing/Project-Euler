#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;

class DuodigitSolver {
  public:
    DuodigitSolver() {
        for (int a = 0; a <= 9; ++a) {
            for (int b = a; b <= 9; ++b) {
                digit_pairs_.emplace_back(a, b);
            }
        }
    }

    std::string d_of_n(const int n) const {
        std::vector<int> visited(static_cast<std::size_t>(n), 0);
        std::vector<int> parent(static_cast<std::size_t>(n), -2);
        std::vector<int> depth(static_cast<std::size_t>(n), 0);
        std::vector<unsigned char> parent_digit(static_cast<std::size_t>(n), 0);
        std::vector<int> queue;
        queue.reserve(static_cast<std::size_t>(n));

        int token = 0;
        std::string best;
        int best_len = 0;

        for (const auto [a, b] : digit_pairs_) {
            if (a == 0 && b == 0) {
                continue;
            }

            int edges[2] = {a, b};
            int edge_count = (a == b) ? 1 : 2;
            if (edge_count == 2 && edges[0] > edges[1]) {
                std::swap(edges[0], edges[1]);
            }

            int starts[2];
            int start_count = 0;
            for (int i = 0; i < edge_count; ++i) {
                if (edges[i] != 0 && (start_count == 0 || starts[start_count - 1] != edges[i])) {
                    starts[start_count++] = edges[i];
                }
            }
            if (start_count == 0) {
                continue;
            }

            ++token;
            queue.clear();
            int head = 0;
            int found = -1;

            for (int i = 0; i < start_count; ++i) {
                const int d = starts[i];
                const int r = d % n;
                if (visited[static_cast<std::size_t>(r)] == token) {
                    continue;
                }
                visited[static_cast<std::size_t>(r)] = token;
                parent[static_cast<std::size_t>(r)] = -1;
                parent_digit[static_cast<std::size_t>(r)] = static_cast<unsigned char>(d);
                depth[static_cast<std::size_t>(r)] = 1;
                queue.push_back(r);
            }

            if (visited[0] == token) {
                found = 0;
            }

            while (found == -1 && head < static_cast<int>(queue.size())) {
                const int r = queue[static_cast<std::size_t>(head++)];
                const int next_depth = depth[static_cast<std::size_t>(r)] + 1;

                if (best_len > 0 && next_depth > best_len) {
                    continue;
                }

                for (int i = 0; i < edge_count; ++i) {
                    const int d = edges[i];
                    const int nr = static_cast<int>((static_cast<long long>(r) * 10 + d) % n);
                    if (visited[static_cast<std::size_t>(nr)] == token) {
                        continue;
                    }

                    visited[static_cast<std::size_t>(nr)] = token;
                    parent[static_cast<std::size_t>(nr)] = r;
                    parent_digit[static_cast<std::size_t>(nr)] = static_cast<unsigned char>(d);
                    depth[static_cast<std::size_t>(nr)] = next_depth;

                    if (nr == 0) {
                        found = 0;
                        break;
                    }

                    queue.push_back(nr);
                }
            }

            if (found == -1) {
                continue;
            }

            std::string candidate(static_cast<std::size_t>(depth[0]), '0');
            int cur = 0;
            for (int i = depth[0] - 1; i >= 0; --i) {
                candidate[static_cast<std::size_t>(i)] =
                    static_cast<char>('0' + parent_digit[static_cast<std::size_t>(cur)]);
                cur = parent[static_cast<std::size_t>(cur)];
                if (cur == -1) {
                    break;
                }
            }

            if (best.empty() || candidate.size() < best.size() ||
                (candidate.size() == best.size() && candidate < best)) {
                best = std::move(candidate);
                best_len = static_cast<int>(best.size());
            }
        }

        return best;
    }

    cpp_int d_of_n_int(const int n) const {
        const std::string s = d_of_n(n);
        cpp_int value = 0;
        for (const char c : s) {
            value = value * 10 + (c - '0');
        }
        return value;
    }

    cpp_int D(const int k) const {
        cpp_int total = 0;
        for (int n = 1; n <= k; ++n) {
            total += d_of_n_int(n);
        }
        return total;
    }

  private:
    std::vector<std::pair<int, int>> digit_pairs_;
};

std::string scientific_13sig(cpp_int value) {
    std::string digits = value.convert_to<std::string>();
    int exponent = static_cast<int>(digits.size()) - 1;

    std::string sig;
    if (digits.size() >= 13) {
        sig = digits.substr(0, 13);
        const int round_digit = (digits.size() > 13) ? (digits[13] - '0') : 0;
        if (round_digit >= 5) {
            int i = 12;
            while (i >= 0 && sig[static_cast<std::size_t>(i)] == '9') {
                sig[static_cast<std::size_t>(i)] = '0';
                --i;
            }
            if (i >= 0) {
                ++sig[static_cast<std::size_t>(i)];
            } else {
                sig = "1" + sig;
                ++exponent;
            }
        }
        if (sig.size() > 13) {
            sig = sig.substr(0, 13);
        }
    } else {
        sig = digits;
        sig.append(13 - sig.size(), '0');
    }

    return std::string(1, sig[0]) + "." + sig.substr(1) + "e" + std::to_string(exponent);
}

}  // namespace

int main() {
    const DuodigitSolver solver;

    assert(solver.d_of_n(12) == "12");
    assert(solver.d_of_n(102) == "1122");
    assert(solver.d_of_n(103) == "515");
    assert(solver.d_of_n(290) == "11011010");
    assert(solver.d_of_n(317) == "211122");

    assert(solver.D(110) == cpp_int("11047"));
    assert(solver.D(150) == cpp_int("53312"));
    assert(solver.D(500) == cpp_int("29570988"));

    const cpp_int answer = solver.D(50'000);
    std::cout << scientific_13sig(answer) << '\n';
    return 0;
}
