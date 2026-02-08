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
    std::string file = "resources/documents/0102_triangles.txt";
    bool run_checkpoints = true;
};

struct Point {
    i64 x = 0;
    i64 y = 0;
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

i64 cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

bool contains_origin(const Point& a, const Point& b, const Point& c) {
    const i64 c1 = cross(a, b);
    const i64 c2 = cross(b, c);
    const i64 c3 = cross(c, a);

    const bool all_positive = (c1 > 0 && c2 > 0 && c3 > 0);
    const bool all_negative = (c1 < 0 && c2 < 0 && c3 < 0);
    return all_positive || all_negative;
}

std::vector<i64> parse_csv_i64(const std::string& line) {
    std::vector<i64> values;
    std::istringstream input(line);
    std::string token;
    while (std::getline(input, token, ',')) {
        values.push_back(std::stoll(token));
    }
    return values;
}

int solve_from_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    int count = 0;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const std::vector<i64> v = parse_csv_i64(line);
        if (v.size() != 6) {
            throw std::runtime_error("Malformed triangle row");
        }

        const Point a{v[0], v[1]};
        const Point b{v[2], v[3]};
        const Point c{v[4], v[5]};
        if (contains_origin(a, b, c)) {
            ++count;
        }
    }

    return count;
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open triangles file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return solve_from_text(buffer.str());
}

bool run_checkpoints() {
    const std::string sample =
        "-340,495,-153,-910,835,-947\n"
        "-175,41,-421,-714,574,-645\n";

    if (solve_from_text(sample) != 1) {
        std::cerr << "Checkpoint failed for statement sample" << '\n';
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
