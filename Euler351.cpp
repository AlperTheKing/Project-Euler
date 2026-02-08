#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 order = 100'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0ULL;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--order=", options.order)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.order >= 1ULL;
}

class TotientSummatory {
public:
    TotientSummatory() {
        memo_[0ULL] = 0ULL;
        memo_[1ULL] = 1ULL;
    }

    u64 sum(const u64 n) {
        const auto it = memo_.find(n);
        if (it != memo_.end()) {
            return it->second;
        }
        u128 total = static_cast<u128>(n) * static_cast<u128>(n + 1ULL) / 2ULL;
        for (u64 left = 2ULL; left <= n;) {
            const u64 quotient = n / left;
            const u64 right = n / quotient;
            total -= static_cast<u128>(right - left + 1ULL) * static_cast<u128>(sum(quotient));
            left = right + 1ULL;
        }
        const u64 result = static_cast<u64>(total);
        memo_[n] = result;
        return result;
    }

private:
    std::unordered_map<u64, u64> memo_;
};

u64 hidden_points(const u64 order) {
    TotientSummatory summatory;
    const u64 triangle = order * (order + 1ULL) / 2ULL;
    const u64 visible_per_sector = summatory.sum(order);
    const u64 hidden_per_sector = triangle - visible_per_sector;
    return hidden_per_sector * 6ULL;
}

u64 brute_hidden_points(const int order) {
    u64 hidden_per_sector = 0ULL;
    for (int r = 1; r <= order; ++r) {
        int visible = 0;
        for (int x = 1; x <= r; ++x) {
            if (std::gcd(x, r) == 1) {
                ++visible;
            }
        }
        hidden_per_sector += static_cast<u64>(r - visible);
    }
    return hidden_per_sector * 6ULL;
}

bool run_checkpoints() {
    if (hidden_points(5ULL) != 30ULL) {
        std::cerr << "Checkpoint failed for H(5)" << '\n';
        return false;
    }
    if (hidden_points(10ULL) != 138ULL) {
        std::cerr << "Checkpoint failed for H(10)" << '\n';
        return false;
    }
    if (hidden_points(1'000ULL) != 1'177'848ULL) {
        std::cerr << "Checkpoint failed for H(1000)" << '\n';
        return false;
    }
    if (hidden_points(200ULL) != brute_hidden_points(200)) {
        std::cerr << "Checkpoint failed for brute-force cross-check at order=200" << '\n';
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
    std::cout << hidden_points(options.order) << '\n';
    return 0;
}
