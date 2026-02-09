#include <array>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 n = 1000000000000000000ULL;
    bool run_checkpoints = true;
};

struct Node {
    std::array<int, 10> next{};
    int fail = 0;
    bool out = false;
    Node() {
        next.fill(-1);
    }
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL;
}

std::string mul_by_11(const std::string& s) {
    std::string out;
    out.resize(s.size());
    int carry = 0;
    for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) {
        const int v = 11 * (s[static_cast<std::size_t>(i)] - '0') + carry;
        out[static_cast<std::size_t>(i)] = static_cast<char>('0' + (v % 10));
        carry = v / 10;
    }
    while (carry > 0) {
        out.insert(out.begin(), static_cast<char>('0' + (carry % 10)));
        carry /= 10;
    }
    std::size_t pos = 0;
    while (pos + 1 < out.size() && out[pos] == '0') {
        ++pos;
    }
    if (pos > 0) {
        out.erase(0, pos);
    }
    return out;
}

std::vector<std::string> generate_forbidden(int max_len) {
    std::vector<std::string> patterns;
    std::string cur = "11";  // powers 11^k, k>=1
    while (static_cast<int>(cur.size()) <= max_len) {
        patterns.push_back(cur);
        cur = mul_by_11(cur);
    }
    return patterns;
}

std::vector<Node> build_aho(const std::vector<std::string>& patterns) {
    std::vector<Node> nodes(1);

    for (const std::string& p : patterns) {
        int state = 0;
        for (const char ch : p) {
            const int d = ch - '0';
            int nx = nodes[static_cast<std::size_t>(state)].next[static_cast<std::size_t>(d)];
            if (nx == -1) {
                nx = static_cast<int>(nodes.size());
                nodes[static_cast<std::size_t>(state)].next[static_cast<std::size_t>(d)] = nx;
                nodes.push_back(Node{});
            }
            state = nx;
        }
        nodes[static_cast<std::size_t>(state)].out = true;
    }

    std::queue<int> q;
    for (int d = 0; d <= 9; ++d) {
        int& nx = nodes[0].next[static_cast<std::size_t>(d)];
        if (nx == -1) {
            nx = 0;
        } else {
            nodes[static_cast<std::size_t>(nx)].fail = 0;
            q.push(nx);
        }
    }

    while (!q.empty()) {
        const int v = q.front();
        q.pop();

        const int f = nodes[static_cast<std::size_t>(v)].fail;
        if (nodes[static_cast<std::size_t>(f)].out) {
            nodes[static_cast<std::size_t>(v)].out = true;
        }

        for (int d = 0; d <= 9; ++d) {
            int& nx = nodes[static_cast<std::size_t>(v)].next[static_cast<std::size_t>(d)];
            if (nx == -1) {
                nx = nodes[static_cast<std::size_t>(f)].next[static_cast<std::size_t>(d)];
            } else {
                nodes[static_cast<std::size_t>(nx)].fail =
                    nodes[static_cast<std::size_t>(f)].next[static_cast<std::size_t>(d)];
                q.push(nx);
            }
        }
    }

    return nodes;
}

u128 cap_add(u128 a, u128 b, u128 cap) {
    const u128 s = a + b;
    if (s >= cap || s < a) {
        return cap;
    }
    return s;
}

std::string nth_eleven_free(u64 n_target) {
    const int max_digits = 64;
    const std::vector<std::string> patterns = generate_forbidden(max_digits);
    const std::vector<Node> automaton = build_aho(patterns);
    const int states = static_cast<int>(automaton.size());

    const u128 cap = static_cast<u128>(n_target) + 1ULL;

    std::vector<std::vector<u128>> ways(
        static_cast<std::size_t>(max_digits + 1),
        std::vector<u128>(static_cast<std::size_t>(states), 0));

    for (int s = 0; s < states; ++s) {
        ways[0][static_cast<std::size_t>(s)] = 1;
    }

    for (int rem = 1; rem <= max_digits; ++rem) {
        for (int s = 0; s < states; ++s) {
            if (automaton[static_cast<std::size_t>(s)].out) {
                ways[static_cast<std::size_t>(rem)][static_cast<std::size_t>(s)] = 0;
                continue;
            }
            u128 sum = 0;
            for (int d = 0; d <= 9; ++d) {
                const int ns = automaton[static_cast<std::size_t>(s)].next[static_cast<std::size_t>(d)];
                if (automaton[static_cast<std::size_t>(ns)].out) {
                    continue;
                }
                sum = cap_add(sum, ways[static_cast<std::size_t>(rem - 1)][static_cast<std::size_t>(ns)], cap);
            }
            ways[static_cast<std::size_t>(rem)][static_cast<std::size_t>(s)] = sum;
        }
    }

    u128 remaining = static_cast<u128>(n_target);
    int len = 0;
    for (int L = 1; L <= max_digits; ++L) {
        u128 cnt = 0;
        for (int d = 1; d <= 9; ++d) {
            const int ns = automaton[0].next[static_cast<std::size_t>(d)];
            if (automaton[static_cast<std::size_t>(ns)].out) {
                continue;
            }
            cnt = cap_add(cnt, ways[static_cast<std::size_t>(L - 1)][static_cast<std::size_t>(ns)], cap);
        }
        if (remaining > cnt) {
            remaining -= cnt;
        } else {
            len = L;
            break;
        }
    }

    if (len == 0) {
        return {};
    }

    std::string answer;
    answer.reserve(static_cast<std::size_t>(len));
    int state = 0;

    for (int pos = 0; pos < len; ++pos) {
        const int start_digit = (pos == 0) ? 1 : 0;
        for (int d = start_digit; d <= 9; ++d) {
            const int ns = automaton[static_cast<std::size_t>(state)].next[static_cast<std::size_t>(d)];
            if (automaton[static_cast<std::size_t>(ns)].out) {
                continue;
            }
            const u128 cnt =
                ways[static_cast<std::size_t>(len - pos - 1)][static_cast<std::size_t>(ns)];
            if (remaining > cnt) {
                remaining -= cnt;
            } else {
                answer.push_back(static_cast<char>('0' + d));
                state = ns;
                break;
            }
        }
    }

    return answer;
}

bool run_checkpoints() {
    if (nth_eleven_free(3) != "3") {
        std::cerr << "Checkpoint failed: E(3)\n";
        return false;
    }
    if (nth_eleven_free(200) != "213") {
        std::cerr << "Checkpoint failed: E(200)\n";
        return false;
    }
    if (nth_eleven_free(500000) != "531563") {
        std::cerr << "Checkpoint failed: E(500000)\n";
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

    std::cout << nth_eleven_free(options.n) << '\n';
    return 0;
}
