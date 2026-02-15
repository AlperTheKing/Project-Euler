#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kDefaultN = 100'000'000ULL;

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        u64 parsed = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed)) {
            options.n = parsed;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

inline u64 next_value(const u64 x) {
    return static_cast<u64>((static_cast<u128>(x) * static_cast<u128>(x) + 45ULL) % kMod);
}

i64 compute_score(const std::vector<int>& seq) {
    std::vector<i64> stack;
    stack.reserve(seq.size());

    i64 sum = 0;
    for (int v : seq) {
        sum += v;
        stack.push_back(static_cast<i64>(v));
        while (stack.size() >= 3) {
            const std::size_t n = stack.size();
            const i64 a = stack[n - 3];
            const i64 b = stack[n - 2];
            const i64 c = stack[n - 1];
            if (b < a || b < c) {
                break;
            }
            stack[n - 3] = a - b + c;
            stack.resize(n - 2);
        }
    }

    i64 delta = 0;
    std::size_t i = 0;
    std::size_t j = (stack.empty() ? 0 : stack.size() - 1);
    int sign = 1;
    while (!stack.empty() && i <= j) {
        if (stack[i] > stack[j]) {
            delta += sign * stack[i];
            ++i;
        } else {
            delta += sign * stack[j];
            if (j == 0) {
                break;
            }
            --j;
        }
        sign = -sign;
    }

    return (sum + delta) / 2;
}

i64 solve_n(const u64 n) {
    if (n == 0) {
        return 0;
    }

    std::unordered_map<int, int> seen;
    seen.reserve(100000);
    std::vector<int> prefix;
    prefix.reserve(100000);

    int s = 0;
    int cycle_start = -1;
    int cycle_len = 0;

    for (u64 i = 0; i < n; ++i) {
        auto it = seen.find(s);
        if (it != seen.end()) {
            cycle_start = it->second;
            cycle_len = static_cast<int>(prefix.size()) - cycle_start;
            break;
        }
        seen.emplace(s, static_cast<int>(prefix.size()));
        prefix.push_back(s);
        s = static_cast<int>(next_value(static_cast<u64>(s)));
    }

    if (cycle_start < 0) {
        return compute_score(prefix);
    }

    const u64 p = static_cast<u64>(cycle_start);
    const u64 cycle = static_cast<u64>(cycle_len);

    const u64 a = (n - p) % cycle;
    const u64 b = cycle - a;

    auto append_range = [&](std::vector<int>& out, u64 start, u64 len) {
        if (len == 0) {
            return;
        }
        out.insert(out.end(),
                   prefix.begin() + static_cast<std::ptrdiff_t>(start),
                   prefix.begin() + static_cast<std::ptrdiff_t>(start + len));
    };

    std::vector<int> U;
    U.reserve(static_cast<std::size_t>(p + a));
    append_range(U, 0, p);
    append_range(U, p, a);

    std::vector<int> V;
    V.reserve(static_cast<std::size_t>(b + a));
    append_range(V, p + a, b);
    append_range(V, p, a);

    std::vector<int> T;
    T.reserve(U.size() + V.size());
    T.insert(T.end(), U.begin(), U.end());
    T.insert(T.end(), V.begin(), V.end());

    const i64 s1 = compute_score(U);
    const i64 s2 = compute_score(T);

    const u64 u = U.size();
    const u64 v = V.size();
    const u64 blocks = (n - u) / v;

    const __int128 result =
        static_cast<__int128>(s1) +
        static_cast<__int128>(s2 - s1) * static_cast<__int128>(blocks);
    return static_cast<i64>(result);
}

bool run_checkpoints() {
    struct Checkpoint {
        u64 n;
        i64 expected;
    };
    const Checkpoint checkpoints[] = {
        {2, 45},
        {4, 4'284'990},
        {100, 26'365'463'243},
        {10'000, 2'495'838'522'951},
    };

    for (const auto& checkpoint : checkpoints) {
        const i64 got = solve_n(checkpoint.n);
        if (got != checkpoint.expected) {
            std::cerr << "Checkpoint failed at n=" << checkpoint.n
                      << ": got " << got << ", expected " << checkpoint.expected << '\n';
            return false;
        }
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
        return 1;
    }

    const i64 answer = solve_n(options.n);
    std::cout << answer << '\n';
    return 0;
}
