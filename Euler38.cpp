#include <algorithm>
#include <iostream>
#include <string>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

bool is_1_to_9_pandigital(const std::string& s) {
    if (s.size() != 9U) {
        return false;
    }

    bool used[10] = {false};
    for (const char c : s) {
        const int d = c - '0';
        if (d <= 0 || d >= 10 || used[d]) {
            return false;
        }
        used[d] = true;
    }
    for (int d = 1; d <= 9; ++d) {
        if (!used[d]) {
            return false;
        }
    }
    return true;
}

std::string concatenated_product(const int base, const int max_multiplier) {
    std::string out;
    out.reserve(16U);

    for (int n = 1; n <= max_multiplier; ++n) {
        out += std::to_string(base * n);
    }

    return out;
}

int solve() {
    int best = 0;
    for (int base = 1; base <= 9999; ++base) {
        std::string concat;
        for (int n = 1; concat.size() < 9U; ++n) {
            concat += std::to_string(base * n);
            if (concat.size() == 9U && n >= 2 && is_1_to_9_pandigital(concat)) {
                best = std::max(best, std::stoi(concat));
            }
        }
    }
    return best;
}

bool run_checkpoints() {
    if (concatenated_product(192, 3) != "192384576") {
        std::cerr << "Checkpoint failed for base=192" << '\n';
        return false;
    }
    if (!is_1_to_9_pandigital("123456789")) {
        std::cerr << "Checkpoint failed for pandigital helper" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
