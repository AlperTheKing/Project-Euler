#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0089_roman.txt";
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
        if (parse_string_after_prefix(arg, "--file=", options.file)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

int roman_value(const char c) {
    switch (c) {
        case 'I':
            return 1;
        case 'V':
            return 5;
        case 'X':
            return 10;
        case 'L':
            return 50;
        case 'C':
            return 100;
        case 'D':
            return 500;
        case 'M':
            return 1000;
        default:
            return 0;
    }
}

int roman_to_int(const std::string& s) {
    int total = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int v = roman_value(s[i]);
        if (i + 1 < s.size() && roman_value(s[i + 1]) > v) {
            total -= v;
        } else {
            total += v;
        }
    }
    return total;
}

std::string int_to_minimal_roman(int value) {
    static const std::vector<std::pair<int, std::string>> table = {
        {1000, "M"}, {900, "CM"}, {500, "D"}, {400, "CD"}, {100, "C"}, {90, "XC"},
        {50, "L"},   {40, "XL"},  {10, "X"},  {9, "IX"},   {5, "V"},   {4, "IV"},
        {1, "I"},
    };

    std::string out;
    for (const auto& [v, tok] : table) {
        while (value >= v) {
            out += tok;
            value -= v;
        }
    }
    return out;
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open roman file: " + file_path);
    }

    int saved = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        const int value = roman_to_int(line);
        const std::string minimal = int_to_minimal_roman(value);
        saved += static_cast<int>(line.size() - minimal.size());
    }

    return saved;
}

bool run_checkpoints() {
    if (int_to_minimal_roman(roman_to_int("VIIII")) != "IX") {
        std::cerr << "Checkpoint failed for VIIII" << '\n';
        return false;
    }
    if (int_to_minimal_roman(roman_to_int("LXXXX")) != "XC") {
        std::cerr << "Checkpoint failed for LXXXX" << '\n';
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
        std::cout << solve(options.file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
