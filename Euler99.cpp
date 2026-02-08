#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    std::string file = "resources/documents/0099_base_exp.txt";
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

int solve_from_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;

    int best_line = -1;
    long double best_score = -1.0L;
    int line_number = 0;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        ++line_number;
        std::size_t comma = line.find(',');
        if (comma == std::string::npos) {
            throw std::runtime_error("Malformed line in base/exp list");
        }

        const long double base = std::stold(line.substr(0, comma));
        const long double exponent = std::stold(line.substr(comma + 1));
        const long double score = exponent * std::log(base);

        if (score > best_score) {
            best_score = score;
            best_line = line_number;
        }
    }

    if (best_line < 0) {
        throw std::runtime_error("No data rows found in base/exp list");
    }

    return best_line;
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open base/exp file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return solve_from_text(buffer.str());
}

bool run_checkpoints() {
    const std::string sample = "2,11\n3,7\n";
    if (solve_from_text(sample) != 2) {
        std::cerr << "Checkpoint failed for small base/exp sample" << '\n';
        return false;
    }

    const std::string statement_pair = "632382,518061\n519432,525806\n";
    if (solve_from_text(statement_pair) != 1) {
        std::cerr << "Checkpoint failed for statement comparison pair" << '\n';
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
