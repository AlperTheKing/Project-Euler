#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr int ALPHA = 4;

int char_id(char ch) {
    switch (ch) {
        case 'A':
            return 0;
        case 'E':
            return 1;
        case 'F':
            return 2;
        default:
            return 3;
    }
}

struct Node {
    std::array<int, ALPHA> next{};
    int link = 0;
    std::array<int, 4> out{};

    Node() {
        next.fill(-1);
        out.fill(0);
    }
};

struct Automaton {
    std::vector<Node> nodes;

    Automaton() { nodes.emplace_back(); }

    void add_pattern(const std::string& s, int id) {
        int v = 0;
        for (char ch : s) {
            const int c = char_id(ch);
            if (nodes[static_cast<std::size_t>(v)].next[static_cast<std::size_t>(c)] == -1) {
                nodes[static_cast<std::size_t>(v)].next[static_cast<std::size_t>(c)] =
                    static_cast<int>(nodes.size());
                nodes.emplace_back();
            }
            v = nodes[static_cast<std::size_t>(v)].next[static_cast<std::size_t>(c)];
        }
        nodes[static_cast<std::size_t>(v)].out[static_cast<std::size_t>(id)] += 1;
    }

    void build() {
        std::queue<int> q;

        for (int c = 0; c < ALPHA; ++c) {
            int& to = nodes[0].next[static_cast<std::size_t>(c)];
            if (to == -1) {
                to = 0;
            } else {
                nodes[static_cast<std::size_t>(to)].link = 0;
                q.push(to);
            }
        }

        while (!q.empty()) {
            const int v = q.front();
            q.pop();

            const int lk = nodes[static_cast<std::size_t>(v)].link;
            for (int k = 0; k < 4; ++k) {
                nodes[static_cast<std::size_t>(v)].out[static_cast<std::size_t>(k)] +=
                    nodes[static_cast<std::size_t>(lk)].out[static_cast<std::size_t>(k)];
            }

            for (int c = 0; c < ALPHA; ++c) {
                int& to = nodes[static_cast<std::size_t>(v)].next[static_cast<std::size_t>(c)];
                if (to == -1) {
                    to = nodes[static_cast<std::size_t>(lk)].next[static_cast<std::size_t>(c)];
                } else {
                    nodes[static_cast<std::size_t>(to)].link =
                        nodes[static_cast<std::size_t>(lk)].next[static_cast<std::size_t>(c)];
                    q.push(to);
                }
            }
        }
    }
};

int encode_counts(int a, int b, int c, int d) {
    return ((a * 3 + b) * 3 + c) * 3 + d;
}

u64 f(int n) {
    Automaton ac;
    ac.add_pattern("FREE", 0);
    ac.add_pattern("FARE", 1);
    ac.add_pattern("AREA", 2);
    ac.add_pattern("REEF", 3);
    ac.build();

    const int states = static_cast<int>(ac.nodes.size());
    constexpr int CNT = 81;

    std::vector<u64> dp(static_cast<std::size_t>(states * CNT), 0ULL);
    std::vector<u64> next(static_cast<std::size_t>(states * CNT), 0ULL);

    dp[0] = 1ULL;

    for (int pos = 0; pos < n; ++pos) {
        std::fill(next.begin(), next.end(), 0ULL);

        for (int st = 0; st < states; ++st) {
            const auto& out = ac.nodes[static_cast<std::size_t>(st)].out;
            for (int code = 0; code < CNT; ++code) {
                const u64 ways = dp[static_cast<std::size_t>(st * CNT + code)];
                if (ways == 0ULL) {
                    continue;
                }

                int t = code;
                int c3 = t % 3;
                t /= 3;
                int c2 = t % 3;
                t /= 3;
                int c1 = t % 3;
                t /= 3;
                int c0 = t;

                for (int ch = 0; ch < ALPHA; ++ch) {
                    const int ns = ac.nodes[static_cast<std::size_t>(st)].next[static_cast<std::size_t>(ch)];
                    const auto& add = ac.nodes[static_cast<std::size_t>(ns)].out;

                    const int n0 = (c0 + add[0] >= 2) ? 2 : (c0 + add[0]);
                    const int n1 = (c1 + add[1] >= 2) ? 2 : (c1 + add[1]);
                    const int n2 = (c2 + add[2] >= 2) ? 2 : (c2 + add[2]);
                    const int n3 = (c3 + add[3] >= 2) ? 2 : (c3 + add[3]);

                    if (n0 == 2 || n1 == 2 || n2 == 2 || n3 == 2) {
                        continue;
                    }

                    const int ncode = encode_counts(n0, n1, n2, n3);
                    next[static_cast<std::size_t>(ns * CNT + ncode)] += ways;
                }
            }
        }

        dp.swap(next);
    }

    const int target = encode_counts(1, 1, 1, 1);
    u64 ans = 0ULL;
    for (int st = 0; st < states; ++st) {
        ans += dp[static_cast<std::size_t>(st * CNT + target)];
    }
    return ans;
}

}  // namespace

int main() {
    assert(f(9) == 1ULL);
    assert(f(15) == 72'863ULL);

    std::cout << f(30) << "\n";
    return 0;
}
