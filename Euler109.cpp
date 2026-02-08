#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 100;
    bool run_checkpoints = true;
};

struct Throw {
    std::string label;
    int score = 0;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 2;
}

std::vector<Throw> all_nonfinal_throws() {
    std::vector<Throw> out;
    out.reserve(62);

    for (int n = 1; n <= 20; ++n) {
        out.push_back({"S" + std::to_string(n), n});
        out.push_back({"D" + std::to_string(n), 2 * n});
        out.push_back({"T" + std::to_string(n), 3 * n});
    }

    out.push_back({"S25", 25});
    out.push_back({"D25", 50});
    return out;
}

std::vector<Throw> final_doubles() {
    std::vector<Throw> out;
    out.reserve(21);

    for (int n = 1; n <= 20; ++n) {
        out.push_back({"D" + std::to_string(n), 2 * n});
    }
    out.push_back({"D25", 50});
    return out;
}

int solve(const int limit) {
    const std::vector<Throw> first_throws = all_nonfinal_throws();
    const std::vector<Throw> doubles = final_doubles();

    int count = 0;

    for (const Throw& d : doubles) {
        if (d.score < limit) {
            ++count;
        }
    }

    for (const Throw& a : first_throws) {
        for (const Throw& d : doubles) {
            if (a.score + d.score < limit) {
                ++count;
            }
        }
    }

    for (std::size_t i = 0; i < first_throws.size(); ++i) {
        for (std::size_t j = i; j < first_throws.size(); ++j) {
            const int first_two = first_throws[i].score + first_throws[j].score;
            for (const Throw& d : doubles) {
                if (first_two + d.score < limit) {
                    ++count;
                }
            }
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(6) != 11) {
        std::cerr << "Checkpoint failed for limit=6" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
