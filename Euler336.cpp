#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

void reverse_suffix(std::string& s, const int idx) {
    std::reverse(s.begin() + idx, s.end());
}

int simon_moves(std::string s) {
    const std::string sorted = [&]() {
        std::string t = s;
        std::sort(t.begin(), t.end());
        return t;
    }();

    int moves = 0;
    const int n = static_cast<int>(s.size());

    for (int i = 0; i <= n - 2; ++i) {
        const char target = sorted[static_cast<std::size_t>(i)];
        int pos = -1;
        for (int j = i; j < n; ++j) {
            if (s[static_cast<std::size_t>(j)] == target) {
                pos = j;
                break;
            }
        }

        if (pos == i) {
            continue;
        }

        if (pos != n - 1) {
            reverse_suffix(s, pos);
            ++moves;
        }

        reverse_suffix(s, i);
        ++moves;
    }

    return moves;
}

std::vector<std::string> generate_maximix_arrangements(const int n) {
    std::string start;
    start.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        start.push_back(static_cast<char>('A' + i));
    }

    std::vector<std::string> states{start};

    // Inverse of Simon's mandatory last move (for i = n-2).
    for (std::string& s : states) {
        reverse_suffix(s, n - 2);
    }

    // For i = n-3 ... 0, inverse of two moves: R(p) after R(i), with p in [i+1, n-2].
    for (int i = n - 3; i >= 0; --i) {
        std::vector<std::string> next;
        const std::size_t branch = static_cast<std::size_t>(n - i - 2);
        next.reserve(states.size() * branch);

        for (const std::string& after : states) {
            std::string after_ri = after;
            reverse_suffix(after_ri, i);

            for (int p = i + 1; p <= n - 2; ++p) {
                std::string before = after_ri;
                reverse_suffix(before, p);
                next.push_back(std::move(before));
            }
        }

        states.swap(next);
    }

    std::sort(states.begin(), states.end());
    return states;
}

bool run_checkpoints() {
    {
        const std::vector<std::string> mm4 = generate_maximix_arrangements(4);
        if (mm4.size() != 2 || mm4[0] != "DACB" || mm4[1] != "DBAC") {
            std::cerr << "Checkpoint failed: n=4 maximix set mismatch\n";
            return false;
        }
    }

    {
        const std::vector<std::string> mm6 = generate_maximix_arrangements(6);
        if (mm6.size() != 24 || mm6[9] != "DFAECB") {
            std::cerr << "Checkpoint failed: n=6 count or 10th lexicographic mismatch\n";
            return false;
        }

        const int max_moves = 2 * 6 - 3;
        for (const std::string& s : mm6) {
            if (simon_moves(s) != max_moves) {
                std::cerr << "Checkpoint failed: generated arrangement is not maximix for n=6\n";
                return false;
            }
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const std::vector<std::string> mm11 = generate_maximix_arrangements(11);
    std::cout << mm11[2010] << '\n';
    return 0;
}
