#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    std::string file = "resources/documents/0098_words.txt";
    bool run_checkpoints = true;
};

struct SquareCache {
    std::unordered_map<std::string, std::vector<std::string>> pattern_to_squares;
    std::unordered_set<u64> square_values;
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

std::vector<std::string> parse_words_csv(const std::string& text) {
    std::vector<std::string> words;
    std::string token;
    std::istringstream input(text);

    while (std::getline(input, token, ',')) {
        std::string word;
        for (char c : token) {
            if (c >= 'A' && c <= 'Z') {
                word.push_back(c);
            }
        }
        if (!word.empty()) {
            words.push_back(word);
        }
    }

    return words;
}

std::string pattern_key(const std::string& s) {
    std::unordered_map<char, int> index;
    int next_index = 0;
    std::string key;

    for (char c : s) {
        auto it = index.find(c);
        if (it == index.end()) {
            it = index.emplace(c, next_index++).first;
        }
        key += std::to_string(it->second);
        key.push_back(',');
    }

    return key;
}

u64 pow10_u64(int n) {
    u64 value = 1;
    for (int i = 0; i < n; ++i) {
        value *= 10ULL;
    }
    return value;
}

SquareCache build_square_cache_for_length(const int length) {
    const u64 low = pow10_u64(length - 1);
    const u64 high = pow10_u64(length) - 1ULL;

    const u64 start = static_cast<u64>(std::ceil(std::sqrt(static_cast<long double>(low))));
    const u64 end = static_cast<u64>(std::floor(std::sqrt(static_cast<long double>(high))));

    SquareCache cache;
    for (u64 x = start; x <= end; ++x) {
        const u64 sq = x * x;
        const std::string digits = std::to_string(sq);
        cache.square_values.insert(sq);
        cache.pattern_to_squares[pattern_key(digits)].push_back(digits);
    }

    return cache;
}

bool map_word_to_digits(const std::string& word,
                        const std::string& digits,
                        std::array<int, 26>& letter_to_digit,
                        std::array<int, 10>& digit_to_letter) {
    letter_to_digit.fill(-1);
    digit_to_letter.fill(-1);

    for (std::size_t i = 0; i < word.size(); ++i) {
        const int letter = word[i] - 'A';
        const int digit = digits[i] - '0';

        const int mapped_digit = letter_to_digit[static_cast<std::size_t>(letter)];
        if (mapped_digit != -1 && mapped_digit != digit) {
            return false;
        }
        const int mapped_letter = digit_to_letter[static_cast<std::size_t>(digit)];
        if (mapped_letter != -1 && mapped_letter != letter) {
            return false;
        }

        letter_to_digit[static_cast<std::size_t>(letter)] = digit;
        digit_to_letter[static_cast<std::size_t>(digit)] = letter;
    }

    return true;
}

u64 mapped_value(const std::string& word, const std::array<int, 26>& letter_to_digit) {
    if (word.empty()) {
        return 0;
    }

    if (letter_to_digit[static_cast<std::size_t>(word[0] - 'A')] == 0) {
        return 0;
    }

    u64 value = 0;
    for (char c : word) {
        const int digit = letter_to_digit[static_cast<std::size_t>(c - 'A')];
        if (digit < 0) {
            return 0;
        }
        value = value * 10ULL + static_cast<u64>(digit);
    }
    return value;
}

u64 solve_from_words(const std::vector<std::string>& words) {
    std::unordered_map<std::string, std::vector<std::string>> anagram_groups;
    for (const std::string& word : words) {
        std::string sorted = word;
        std::sort(sorted.begin(), sorted.end());
        anagram_groups[sorted].push_back(word);
    }

    std::unordered_map<int, SquareCache> cache_by_length;
    for (const auto& [_, group] : anagram_groups) {
        if (group.size() < 2) {
            continue;
        }
        const int length = static_cast<int>(group[0].size());
        if (cache_by_length.find(length) == cache_by_length.end()) {
            cache_by_length.emplace(length, build_square_cache_for_length(length));
        }
    }

    u64 best = 0;
    std::array<int, 26> letter_to_digit{};
    std::array<int, 10> digit_to_letter{};

    for (const auto& [_, group] : anagram_groups) {
        if (group.size() < 2) {
            continue;
        }

        const int length = static_cast<int>(group[0].size());
        const SquareCache& cache = cache_by_length.at(length);

        for (std::size_t i = 0; i < group.size(); ++i) {
            for (std::size_t j = i + 1; j < group.size(); ++j) {
                const std::string& a = group[i];
                const std::string& b = group[j];

                const std::string key = pattern_key(a);
                auto it = cache.pattern_to_squares.find(key);
                if (it == cache.pattern_to_squares.end()) {
                    continue;
                }

                for (const std::string& sq_digits : it->second) {
                    if (!map_word_to_digits(a, sq_digits, letter_to_digit, digit_to_letter)) {
                        continue;
                    }
                    if (letter_to_digit[static_cast<std::size_t>(a[0] - 'A')] == 0) {
                        continue;
                    }

                    const u64 va = std::stoull(sq_digits);
                    const u64 vb = mapped_value(b, letter_to_digit);
                    if (vb == 0) {
                        continue;
                    }
                    if (cache.square_values.find(vb) == cache.square_values.end()) {
                        continue;
                    }

                    best = std::max(best, std::max(va, vb));
                }
            }
        }
    }

    return best;
}

u64 solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open words file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return solve_from_words(parse_words_csv(buffer.str()));
}

bool run_checkpoints() {
    const std::vector<std::string> words = {"CARE", "RACE", "STOP"};
    if (solve_from_words(words) != 9216ULL) {
        std::cerr << "Checkpoint failed for CARE/RACE" << '\n';
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
