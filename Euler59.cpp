#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    std::string file = "resources/documents/0059_cipher.txt";
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

std::vector<int> parse_cipher_csv(const std::string& payload) {
    std::vector<int> values;
    std::string token;

    for (char c : payload) {
        if ((c >= '0' && c <= '9')) {
            token.push_back(c);
        } else if (!token.empty()) {
            values.push_back(std::stoi(token));
            token.clear();
        }
    }
    if (!token.empty()) {
        values.push_back(std::stoi(token));
    }

    return values;
}

std::string decrypt_with_key(const std::vector<int>& cipher, const std::string& key) {
    std::string out;
    out.resize(cipher.size());

    for (std::size_t i = 0; i < cipher.size(); ++i) {
        out[i] = static_cast<char>(cipher[i] ^ key[i % key.size()]);
    }

    return out;
}

int text_score(const std::string& text) {
    int score = 0;

    for (char c : text) {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (uc < 9 || uc > 126) {
            return -1000000;
        }
        if (c == ' ') {
            score += 3;
        }
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) {
            score += 2;
        }
    }

    auto count_occurrences = [&](const std::string& needle) {
        int cnt = 0;
        std::size_t pos = 0;
        while (true) {
            pos = text.find(needle, pos);
            if (pos == std::string::npos) {
                break;
            }
            ++cnt;
            pos += needle.size();
        }
        return cnt;
    };

    score += 40 * count_occurrences(" the ");
    score += 25 * count_occurrences(" and ");
    score += 20 * count_occurrences(" of ");
    score += 20 * count_occurrences(" to ");
    return score;
}

i64 ascii_sum(const std::string& text) {
    i64 total = 0;
    for (char c : text) {
        total += static_cast<unsigned char>(c);
    }
    return total;
}

i64 solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open cipher file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::vector<int> cipher = parse_cipher_csv(buffer.str());

    int best_score = std::numeric_limits<int>::min();
    std::string best_text;

    for (char a = 'a'; a <= 'z'; ++a) {
        for (char b = 'a'; b <= 'z'; ++b) {
            for (char c = 'a'; c <= 'z'; ++c) {
                const std::string key = {a, b, c};
                const std::string text = decrypt_with_key(cipher, key);
                const int score = text_score(text);
                if (score > best_score) {
                    best_score = score;
                    best_text = text;
                }
            }
        }
    }

    return ascii_sum(best_text);
}

bool run_checkpoints() {
    const std::string key = "abc";
    const std::string plain = "hello world";
    std::vector<int> cipher;
    cipher.reserve(plain.size());
    for (std::size_t i = 0; i < plain.size(); ++i) {
        cipher.push_back(static_cast<unsigned char>(plain[i]) ^ key[i % key.size()]);
    }

    if (decrypt_with_key(cipher, key) != plain) {
        std::cerr << "Checkpoint failed for XOR roundtrip" << '\n';
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
