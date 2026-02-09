#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u16 = std::uint16_t;
using u64 = std::uint64_t;

constexpr int LETTERS = 10;
constexpr int MAX_LEN = 9;
constexpr int MAX_SUM = 45;
constexpr u16 LETTER_MASK = (1u << LETTERS) - 1u;     // digits 0..9
constexpr u16 CELL_MASK = LETTER_MASK & ~1u;           // digits 1..9

struct Options {
    std::string file = "p424_kakuro200.txt";
    bool run_checkpoints = true;
};

struct Run {
    std::vector<int> cells;
    std::array<int, 2> clue_letters{};
    int clue_len = 0;  // 1 or 2
};

struct Puzzle {
    int n = 0;
    std::vector<int> white_index;   // n*n -> white id or -1
    std::vector<int> cell_letter;   // white id -> letter id or -1
    std::vector<Run> runs;
};

struct State {
    std::array<u16, LETTERS> letter_dom{};
    std::vector<u16> cell_dom;
};

std::array<std::array<std::vector<u64>, MAX_SUM + 1>, MAX_LEN + 1> tuple_codes;

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg.rfind("--file=", 0U) == 0U) {
            options.file = arg.substr(7);
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return !options.file.empty();
}

inline int popcount(u16 x) {
    return __builtin_popcount(static_cast<unsigned>(x));
}

inline bool is_single(u16 x) {
    return x != 0 && (x & (x - 1u)) == 0;
}

inline int single_value(u16 x) {
    return __builtin_ctz(static_cast<unsigned>(x));
}

void gen_tuples(int len, int pos, int used_mask, int sum, u64 code) {
    if (pos == len) {
        tuple_codes[len][sum].push_back(code);
        return;
    }
    for (int d = 1; d <= 9; ++d) {
        if (used_mask & (1 << d)) {
            continue;
        }
        gen_tuples(len, pos + 1, used_mask | (1 << d), sum + d, code | (static_cast<u64>(d) << (4 * pos)));
    }
}

void prepare_tuples() {
    for (int len = 1; len <= MAX_LEN; ++len) {
        gen_tuples(len, 0, 0, 0, 0);
    }
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(16);
    int depth = 0;
    for (char ch : line) {
        if (ch == ',' && depth == 0) {
            out.push_back(cur);
            cur.clear();
        } else if (ch != ' ' && ch != '\r' && ch != '\n' && ch != '\t') {
            cur.push_back(ch);
            if (ch == '(') {
                ++depth;
            } else if (ch == ')') {
                --depth;
            }
        }
    }
    out.push_back(cur);
    return out;
}

int letter_id(char ch) {
    if (ch < 'A' || ch > 'J') {
        return -1;
    }
    return ch - 'A';
}

bool is_white_token(const std::string& tok) {
    if (tok == "O") {
        return true;
    }
    return tok.size() == 1 && letter_id(tok[0]) != -1;
}

std::vector<std::string> split_by_comma_inside(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == ',') {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    if (!cur.empty()) {
        out.push_back(cur);
    }
    return out;
}

bool parse_puzzle_line(const std::string& line, Puzzle& puzzle) {
    const std::vector<std::string> tokens = split_csv(line);
    if (tokens.empty()) {
        return false;
    }
    int n = 0;
    try {
        n = std::stoi(tokens[0]);
    } catch (...) {
        return false;
    }
    if (n <= 0 || static_cast<int>(tokens.size()) != 1 + n * n) {
        return false;
    }

    std::vector<std::string> grid(tokens.begin() + 1, tokens.end());
    puzzle.n = n;
    puzzle.white_index.assign(static_cast<std::size_t>(n * n), -1);
    puzzle.cell_letter.clear();
    puzzle.runs.clear();

    int white_count = 0;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            const std::string& tok = grid[static_cast<std::size_t>(r * n + c)];
            if (!is_white_token(tok)) {
                continue;
            }
            puzzle.white_index[static_cast<std::size_t>(r * n + c)] = white_count++;
            if (tok == "O") {
                puzzle.cell_letter.push_back(-1);
            } else {
                puzzle.cell_letter.push_back(letter_id(tok[0]));
            }
        }
    }

    auto add_run = [&](int r, int c, bool horizontal, const std::string& clue_str) {
        Run run;
        run.clue_len = static_cast<int>(clue_str.size());
        if (run.clue_len < 1 || run.clue_len > 2) {
            return;
        }
        run.clue_letters[0] = letter_id(clue_str[0]);
        run.clue_letters[1] = (run.clue_len == 2 ? letter_id(clue_str[1]) : -1);
        if (run.clue_letters[0] < 0 || (run.clue_len == 2 && run.clue_letters[1] < 0)) {
            return;
        }

        if (horizontal) {
            for (int cc = c + 1; cc < n; ++cc) {
                const int idx = puzzle.white_index[static_cast<std::size_t>(r * n + cc)];
                if (idx < 0) {
                    break;
                }
                run.cells.push_back(idx);
            }
        } else {
            for (int rr = r + 1; rr < n; ++rr) {
                const int idx = puzzle.white_index[static_cast<std::size_t>(rr * n + c)];
                if (idx < 0) {
                    break;
                }
                run.cells.push_back(idx);
            }
        }
        if (!run.cells.empty()) {
            puzzle.runs.push_back(std::move(run));
        }
    };

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            const std::string& tok = grid[static_cast<std::size_t>(r * n + c)];
            if (tok.size() < 2 || tok.front() != '(' || tok.back() != ')') {
                continue;
            }
            const std::string inside = tok.substr(1, tok.size() - 2);
            const std::vector<std::string> parts = split_by_comma_inside(inside);
            for (const std::string& part : parts) {
                if (part.size() < 2) {
                    return false;
                }
                if (part[0] == 'h') {
                    add_run(r, c, true, part.substr(1));
                } else if (part[0] == 'v') {
                    add_run(r, c, false, part.substr(1));
                } else {
                    return false;
                }
            }
        }
    }

    return true;
}

bool load_puzzles(const std::string& file, std::vector<Puzzle>& puzzles) {
    std::ifstream in(file);
    if (!in) {
        return false;
    }

    puzzles.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        Puzzle p;
        if (!parse_puzzle_line(line, p)) {
            return false;
        }
        puzzles.push_back(std::move(p));
    }
    return !puzzles.empty();
}

bool propagate_links(const Puzzle& puzzle, State& st, bool& changed) {
    for (std::size_t i = 0; i < puzzle.cell_letter.size(); ++i) {
        const int letter = puzzle.cell_letter[i];
        if (letter < 0) {
            continue;
        }
        const u16 allowed = static_cast<u16>(st.letter_dom[letter] & CELL_MASK);
        u16 cdom = st.cell_dom[i];
        u16 ndom = static_cast<u16>(cdom & allowed);
        if (ndom == 0) {
            return false;
        }
        if (ndom != cdom) {
            st.cell_dom[i] = ndom;
            changed = true;
        }

        u16 ldom = st.letter_dom[letter];
        u16 nldom = static_cast<u16>(ldom & ndom);
        if (nldom == 0) {
            return false;
        }
        if (nldom != ldom) {
            st.letter_dom[letter] = nldom;
            changed = true;
        }
    }
    return true;
}

bool propagate_all_diff_letters(State& st, bool& changed) {
    bool local_changed = true;
    while (local_changed) {
        local_changed = false;

        for (int i = 0; i < LETTERS; ++i) {
            if (st.letter_dom[i] == 0) {
                return false;
            }
        }

        u16 assigned_mask = 0;
        for (int i = 0; i < LETTERS; ++i) {
            if (is_single(st.letter_dom[i])) {
                assigned_mask |= st.letter_dom[i];
            }
        }

        for (int i = 0; i < LETTERS; ++i) {
            if (is_single(st.letter_dom[i])) {
                continue;
            }
            u16 ndom = static_cast<u16>(st.letter_dom[i] & ~assigned_mask);
            if (ndom == 0) {
                return false;
            }
            if (ndom != st.letter_dom[i]) {
                st.letter_dom[i] = ndom;
                local_changed = true;
                changed = true;
            }
        }

        for (int d = 0; d <= 9; ++d) {
            const u16 bit = static_cast<u16>(1u << d);
            int count = 0;
            int where = -1;
            for (int i = 0; i < LETTERS; ++i) {
                if (st.letter_dom[i] & bit) {
                    ++count;
                    where = i;
                }
            }
            if (count == 0) {
                return false;
            }
            if (count == 1 && st.letter_dom[where] != bit) {
                st.letter_dom[where] = bit;
                local_changed = true;
                changed = true;
            }
        }
    }
    return true;
}

bool propagate_run(const Run& run, State& st, bool& changed) {
    const int len = static_cast<int>(run.cells.size());
    if (len <= 0 || len > MAX_LEN) {
        return false;
    }

    const int a = run.clue_letters[0];
    const int b = run.clue_letters[1];

    if (run.clue_len == 2) {
        u16 dom_a = st.letter_dom[a];
        u16 ndom_a = static_cast<u16>(dom_a & ~1u);  // no leading zero in two-digit sum.
        if (ndom_a == 0) {
            return false;
        }
        if (ndom_a != dom_a) {
            st.letter_dom[a] = ndom_a;
            changed = true;
        }
    }

    std::array<bool, MAX_SUM + 1> sum_allowed{};
    if (run.clue_len == 1) {
        const u16 dom = st.letter_dom[a];
        for (int d = 0; d <= 9; ++d) {
            if ((dom >> d) & 1u) {
                sum_allowed[d] = true;
            }
        }
    } else {
        const u16 dom_a = st.letter_dom[a];
        const u16 dom_b = st.letter_dom[b];
        for (int da = 0; da <= 9; ++da) {
            if (!((dom_a >> da) & 1u)) {
                continue;
            }
            for (int db = 0; db <= 9; ++db) {
                if (!((dom_b >> db) & 1u)) {
                    continue;
                }
                const int s = 10 * da + db;
                if (s <= MAX_SUM) {
                    sum_allowed[s] = true;
                }
            }
        }
    }

    std::array<u16, MAX_LEN> union_masks{};
    std::array<bool, MAX_SUM + 1> sum_valid{};
    bool any_tuple = false;

    for (int s = 0; s <= MAX_SUM; ++s) {
        if (!sum_allowed[s]) {
            continue;
        }
        const auto& tuples = tuple_codes[len][s];
        if (tuples.empty()) {
            continue;
        }
        bool any_for_sum = false;
        for (u64 code : tuples) {
            bool ok = true;
            for (int pos = 0; pos < len; ++pos) {
                const int d = static_cast<int>((code >> (4 * pos)) & 0xFULL);
                const u16 bit = static_cast<u16>(1u << d);
                if ((st.cell_dom[static_cast<std::size_t>(run.cells[pos])] & bit) == 0) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                continue;
            }
            any_for_sum = true;
            any_tuple = true;
            for (int pos = 0; pos < len; ++pos) {
                const int d = static_cast<int>((code >> (4 * pos)) & 0xFULL);
                union_masks[static_cast<std::size_t>(pos)] |= static_cast<u16>(1u << d);
            }
        }
        if (any_for_sum) {
            sum_valid[s] = true;
        }
    }

    if (!any_tuple) {
        return false;
    }

    for (int pos = 0; pos < len; ++pos) {
        const int cid = run.cells[pos];
        const u16 ndom = static_cast<u16>(st.cell_dom[static_cast<std::size_t>(cid)] & union_masks[static_cast<std::size_t>(pos)]);
        if (ndom == 0) {
            return false;
        }
        if (ndom != st.cell_dom[static_cast<std::size_t>(cid)]) {
            st.cell_dom[static_cast<std::size_t>(cid)] = ndom;
            changed = true;
        }
    }

    if (run.clue_len == 1) {
        u16 allowed = 0;
        for (int s = 0; s <= 9; ++s) {
            if (sum_valid[s]) {
                allowed |= static_cast<u16>(1u << s);
            }
        }
        const u16 ndom = static_cast<u16>(st.letter_dom[a] & allowed);
        if (ndom == 0) {
            return false;
        }
        if (ndom != st.letter_dom[a]) {
            st.letter_dom[a] = ndom;
            changed = true;
        }
    } else if (a == b) {
        u16 allowed = 0;
        for (int d = 0; d <= 9; ++d) {
            const int s = 11 * d;
            if (s <= MAX_SUM && sum_valid[s]) {
                allowed |= static_cast<u16>(1u << d);
            }
        }
        allowed = static_cast<u16>(allowed & ~1u);
        const u16 ndom = static_cast<u16>(st.letter_dom[a] & allowed);
        if (ndom == 0) {
            return false;
        }
        if (ndom != st.letter_dom[a]) {
            st.letter_dom[a] = ndom;
            changed = true;
        }
    } else {
        u16 allowed_a = 0;
        u16 allowed_b = 0;
        const u16 dom_a = st.letter_dom[a];
        const u16 dom_b = st.letter_dom[b];
        for (int da = 0; da <= 9; ++da) {
            if (!((dom_a >> da) & 1u)) {
                continue;
            }
            for (int db = 0; db <= 9; ++db) {
                if (!((dom_b >> db) & 1u)) {
                    continue;
                }
                const int s = 10 * da + db;
                if (s <= MAX_SUM && sum_valid[s]) {
                    allowed_a |= static_cast<u16>(1u << da);
                    allowed_b |= static_cast<u16>(1u << db);
                }
            }
        }
        allowed_a = static_cast<u16>(allowed_a & ~1u);
        const u16 ndom_a = static_cast<u16>(st.letter_dom[a] & allowed_a);
        const u16 ndom_b = static_cast<u16>(st.letter_dom[b] & allowed_b);
        if (ndom_a == 0 || ndom_b == 0) {
            return false;
        }
        if (ndom_a != st.letter_dom[a]) {
            st.letter_dom[a] = ndom_a;
            changed = true;
        }
        if (ndom_b != st.letter_dom[b]) {
            st.letter_dom[b] = ndom_b;
            changed = true;
        }
    }

    return true;
}

bool propagate(const Puzzle& puzzle, State& st) {
    while (true) {
        bool changed = false;

        if (!propagate_links(puzzle, st, changed)) {
            return false;
        }
        if (!propagate_all_diff_letters(st, changed)) {
            return false;
        }
        if (!propagate_links(puzzle, st, changed)) {
            return false;
        }

        for (const Run& run : puzzle.runs) {
            if (!propagate_run(run, st, changed)) {
                return false;
            }
        }

        if (!changed) {
            break;
        }
    }
    return true;
}

bool solved(const State& st) {
    for (int i = 0; i < LETTERS; ++i) {
        if (!is_single(st.letter_dom[i])) {
            return false;
        }
    }
    for (u16 dom : st.cell_dom) {
        if (!is_single(dom)) {
            return false;
        }
    }
    return true;
}

struct Choice {
    bool is_letter = true;
    int idx = -1;
    u16 dom = 0;
    int size = 1000;
};

Choice choose_variable(const State& st) {
    Choice best;
    for (int i = 0; i < LETTERS; ++i) {
        const int sz = popcount(st.letter_dom[i]);
        if (sz > 1 && sz < best.size) {
            best = {true, i, st.letter_dom[i], sz};
        }
    }
    for (int i = 0; i < static_cast<int>(st.cell_dom.size()); ++i) {
        const int sz = popcount(st.cell_dom[static_cast<std::size_t>(i)]);
        if (sz > 1 && sz < best.size) {
            best = {false, i, st.cell_dom[static_cast<std::size_t>(i)], sz};
        }
    }
    return best;
}

bool dfs_solve(const Puzzle& puzzle, State& st) {
    if (!propagate(puzzle, st)) {
        return false;
    }
    if (solved(st)) {
        return true;
    }

    const Choice choice = choose_variable(st);
    if (choice.idx < 0) {
        return false;
    }

    for (int d = 0; d <= 9; ++d) {
        const u16 bit = static_cast<u16>(1u << d);
        if ((choice.dom & bit) == 0) {
            continue;
        }

        State next = st;
        if (choice.is_letter) {
            next.letter_dom[static_cast<std::size_t>(choice.idx)] = bit;
        } else {
            next.cell_dom[static_cast<std::size_t>(choice.idx)] = bit;
        }

        if (dfs_solve(puzzle, next)) {
            st = std::move(next);
            return true;
        }
    }
    return false;
}

bool solve_puzzle(const Puzzle& puzzle, u64& encrypted_value) {
    State st;
    st.letter_dom.fill(LETTER_MASK);
    st.cell_dom.assign(puzzle.cell_letter.size(), CELL_MASK);

    if (!dfs_solve(puzzle, st)) {
        return false;
    }

    encrypted_value = 0ULL;
    for (int i = 0; i < LETTERS; ++i) {
        encrypted_value = encrypted_value * 10ULL + static_cast<u64>(single_value(st.letter_dom[i]));
    }
    return true;
}

bool run_checkpoints(const std::vector<Puzzle>& puzzles) {
    if (puzzles.empty()) {
        std::cerr << "Checkpoint failed: no puzzles loaded\n";
        return false;
    }

    u64 first_value = 0ULL;
    if (!solve_puzzle(puzzles[0], first_value)) {
        std::cerr << "Checkpoint failed: first puzzle unsolved\n";
        return false;
    }
    if (first_value != 8426039571ULL) {
        std::cerr << "Checkpoint failed: first puzzle encrypted value\n";
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

    prepare_tuples();

    std::vector<Puzzle> puzzles;
    if (!load_puzzles(options.file, puzzles)) {
        std::cerr << "Failed to load puzzle file: " << options.file << '\n';
        return 2;
    }

    if (options.run_checkpoints && !run_checkpoints(puzzles)) {
        return 3;
    }

    u64 total = 0ULL;
    for (const Puzzle& puzzle : puzzles) {
        u64 value = 0ULL;
        if (!solve_puzzle(puzzle, value)) {
            std::cerr << "Failed to solve a puzzle\n";
            return 4;
        }
        total += value;
    }

    std::cout << total << '\n';
    return 0;
}
