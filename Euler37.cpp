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

bool is_prime(const int n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (int p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

bool is_truncatable_prime(const int n) {
    if (n < 10 || !is_prime(n)) {
        return false;
    }

    const std::string s = std::to_string(n);
    for (std::size_t cut = 1; cut < s.size(); ++cut) {
        const int left = std::stoi(s.substr(cut));
        const int right = std::stoi(s.substr(0, s.size() - cut));
        if (!is_prime(left) || !is_prime(right)) {
            return false;
        }
    }

    return true;
}

int solve() {
    int count = 0;
    int total = 0;

    for (int n = 11; count < 11; n += 2) {
        if (is_truncatable_prime(n)) {
            total += n;
            ++count;
        }
    }

    return total;
}

bool run_checkpoints() {
    if (!is_truncatable_prime(3797)) {
        std::cerr << "Checkpoint failed for 3797" << '\n';
        return false;
    }
    if (is_truncatable_prime(47)) {
        std::cerr << "Checkpoint failed for 47" << '\n';
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
