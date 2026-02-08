#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    std::string names_file = "resources/documents/0022_names.txt";
    bool run_checkpoints = true;
};

bool parse_string_after_prefix(const std::string& arg,
                               const std::string& prefix,
                               std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    value = tail;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_string_after_prefix(arg, "--names-file=", options.names_file)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

std::vector<std::string> parse_csv_names(const std::string& payload) {
    std::vector<std::string> names;
    std::string current;
    bool inside_quotes = false;

    for (char ch : payload) {
        if (ch == '"') {
            inside_quotes = !inside_quotes;
            if (!inside_quotes) {
                names.push_back(current);
                current.clear();
            }
            continue;
        }
        if (inside_quotes) {
            current.push_back(ch);
        }
    }

    return names;
}

i64 name_value(const std::string& name) {
    i64 value = 0;
    for (char c : name) {
        if (c >= 'A' && c <= 'Z') {
            value += static_cast<i64>(c - 'A' + 1);
        }
    }
    return value;
}

i64 score_sum(std::vector<std::string> names) {
    std::sort(names.begin(), names.end());

    i64 total = 0;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const i64 position = static_cast<i64>(i + 1);
        total += position * name_value(names[i]);
    }
    return total;
}

i64 solve(const std::string& names_file) {
    std::ifstream input(names_file);
    if (!input) {
        throw std::runtime_error("Could not open names file: " + names_file);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    return score_sum(parse_csv_names(buffer.str()));
}

bool run_checkpoints() {
    if (name_value("COLIN") != 53) {
        std::cerr << "Checkpoint failed for COLIN value" << '\n';
        return false;
    }

    {
        const std::vector<std::string> sample = {"B", "A"};
        if (score_sum(sample) != 5) {
            std::cerr << "Checkpoint failed for small score set" << '\n';
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
        return 2;
    }

    try {
        std::cout << solve(options.names_file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
