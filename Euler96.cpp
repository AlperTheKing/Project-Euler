#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0096_sudoku.txt";
    bool run_checkpoints = true;
};

struct Sudoku {
    std::array<int, 81> cell{};
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

int box_index(const int r, const int c) { return (r / 3) * 3 + (c / 3); }

bool solve_sudoku(Sudoku& s,
                  std::array<int, 9>& row_mask,
                  std::array<int, 9>& col_mask,
                  std::array<int, 9>& box_mask) {
    int best_pos = -1;
    int best_count = 10;
    int best_candidates = 0;

    for (int pos = 0; pos < 81; ++pos) {
        if (s.cell[static_cast<std::size_t>(pos)] != 0) {
            continue;
        }

        const int r = pos / 9;
        const int c = pos % 9;
        const int b = box_index(r, c);

        const int used = row_mask[static_cast<std::size_t>(r)] | col_mask[static_cast<std::size_t>(c)] |
                         box_mask[static_cast<std::size_t>(b)];
        const int candidates = (~used) & 0x3FE;  // bits 1..9

        const int count = __builtin_popcount(static_cast<unsigned>(candidates));
        if (count == 0) {
            return false;
        }
        if (count < best_count) {
            best_count = count;
            best_pos = pos;
            best_candidates = candidates;
            if (count == 1) {
                break;
            }
        }
    }

    if (best_pos == -1) {
        return true;
    }

    const int r = best_pos / 9;
    const int c = best_pos % 9;
    const int b = box_index(r, c);

    for (int d = 1; d <= 9; ++d) {
        if (((best_candidates >> d) & 1) == 0) {
            continue;
        }

        s.cell[static_cast<std::size_t>(best_pos)] = d;
        row_mask[static_cast<std::size_t>(r)] |= (1 << d);
        col_mask[static_cast<std::size_t>(c)] |= (1 << d);
        box_mask[static_cast<std::size_t>(b)] |= (1 << d);

        if (solve_sudoku(s, row_mask, col_mask, box_mask)) {
            return true;
        }

        row_mask[static_cast<std::size_t>(r)] &= ~(1 << d);
        col_mask[static_cast<std::size_t>(c)] &= ~(1 << d);
        box_mask[static_cast<std::size_t>(b)] &= ~(1 << d);
        s.cell[static_cast<std::size_t>(best_pos)] = 0;
    }

    return false;
}

bool initialize_masks(const Sudoku& s,
                      std::array<int, 9>& row_mask,
                      std::array<int, 9>& col_mask,
                      std::array<int, 9>& box_mask) {
    row_mask.fill(0);
    col_mask.fill(0);
    box_mask.fill(0);

    for (int pos = 0; pos < 81; ++pos) {
        const int value = s.cell[static_cast<std::size_t>(pos)];
        if (value == 0) {
            continue;
        }

        const int r = pos / 9;
        const int c = pos % 9;
        const int b = box_index(r, c);
        const int bit = 1 << value;

        if ((row_mask[static_cast<std::size_t>(r)] & bit) ||
            (col_mask[static_cast<std::size_t>(c)] & bit) ||
            (box_mask[static_cast<std::size_t>(b)] & bit)) {
            return false;
        }

        row_mask[static_cast<std::size_t>(r)] |= bit;
        col_mask[static_cast<std::size_t>(c)] |= bit;
        box_mask[static_cast<std::size_t>(b)] |= bit;
    }

    return true;
}

int solve_single(Sudoku s) {
    std::array<int, 9> row_mask;
    std::array<int, 9> col_mask;
    std::array<int, 9> box_mask;

    if (!initialize_masks(s, row_mask, col_mask, box_mask)) {
        throw std::runtime_error("Invalid sudoku grid");
    }
    if (!solve_sudoku(s, row_mask, col_mask, box_mask)) {
        throw std::runtime_error("Sudoku has no solution");
    }

    return 100 * s.cell[0] + 10 * s.cell[1] + s.cell[2];
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open sudoku file: " + file_path);
    }

    int total = 0;
    std::string line;
    Sudoku current;
    int row = 0;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.rfind("Grid", 0U) == 0U) {
            row = 0;
            current = Sudoku{};
            continue;
        }

        for (int c = 0; c < 9; ++c) {
            current.cell[static_cast<std::size_t>(row * 9 + c)] = line[static_cast<std::size_t>(c)] - '0';
        }
        ++row;
        if (row == 9) {
            total += solve_single(current);
        }
    }

    return total;
}

bool run_checkpoints() {
    Sudoku easy{};
    const std::string rows[9] = {
        "034678912", "672195348", "198342567", "859761423", "426853791",
        "713924856", "961537284", "287419635", "345286179",
    };
    for (int r = 0; r < 9; ++r) {
        for (int c = 0; c < 9; ++c) {
            easy.cell[static_cast<std::size_t>(r * 9 + c)] = rows[r][static_cast<std::size_t>(c)] - '0';
        }
    }

    if (solve_single(easy) != 534) {
        std::cerr << "Checkpoint failed for easy sudoku" << '\n';
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
