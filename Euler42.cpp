#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    std::string words_file = "resources/documents/0042_words.txt";
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
        if (parse_string_after_prefix(arg, "--words-file=", options.words_file)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

std::vector<std::string> parse_csv_words(const std::string& payload) {
    std::vector<std::string> words;
    std::string current;
    bool inside_quotes = false;

    for (char c : payload) {
        if (c == '"') {
            inside_quotes = !inside_quotes;
            if (!inside_quotes) {
                words.push_back(current);
                current.clear();
            }
            continue;
        }
        if (inside_quotes) {
            current.push_back(c);
        }
    }

    return words;
}

int word_value(const std::string& word) {
    int value = 0;
    for (char c : word) {
        if (c >= 'A' && c <= 'Z') {
            value += c - 'A' + 1;
        }
    }
    return value;
}

std::unordered_set<int> triangle_set_up_to(const int limit) {
    std::unordered_set<int> triangles;
    for (int n = 1;; ++n) {
        const int t = n * (n + 1) / 2;
        if (t > limit) {
            break;
        }
        triangles.insert(t);
    }
    return triangles;
}

int solve(const std::string& words_file) {
    std::ifstream input(words_file);
    if (!input) {
        throw std::runtime_error("Could not open words file: " + words_file);
    }

    std::ostringstream buf;
    buf << input.rdbuf();
    const std::vector<std::string> words = parse_csv_words(buf.str());

    int max_word_value = 0;
    for (const std::string& w : words) {
        max_word_value = std::max(max_word_value, word_value(w));
    }

    const auto triangles = triangle_set_up_to(max_word_value);

    int count = 0;
    for (const std::string& w : words) {
        if (triangles.count(word_value(w)) > 0U) {
            ++count;
        }
    }

    return count;
}

bool run_checkpoints() {
    if (word_value("SKY") != 55) {
        std::cerr << "Checkpoint failed for SKY value" << '\n';
        return false;
    }
    const auto triangles = triangle_set_up_to(100);
    if (triangles.count(55) == 0U || triangles.count(54) != 0U) {
        std::cerr << "Checkpoint failed for triangle-number set" << '\n';
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

    try {
        std::cout << solve(options.words_file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
