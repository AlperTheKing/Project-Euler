#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0054_poker.txt";
    bool run_checkpoints = true;
};

struct Card {
    int rank = 0;
    char suit = 'X';
};

struct HandRank {
    int category = 0;
    std::vector<int> tiebreak;
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

int rank_from_char(const char c) {
    if (c >= '2' && c <= '9') {
        return c - '0';
    }
    switch (c) {
        case 'T':
            return 10;
        case 'J':
            return 11;
        case 'Q':
            return 12;
        case 'K':
            return 13;
        case 'A':
            return 14;
        default:
            return -1;
    }
}

Card parse_card(const std::string& s) {
    if (s.size() != 2U) {
        throw std::runtime_error("Invalid card token: " + s);
    }
    Card c;
    c.rank = rank_from_char(s[0]);
    c.suit = s[1];
    if (c.rank < 2 || c.rank > 14) {
        throw std::runtime_error("Invalid card rank: " + s);
    }
    return c;
}

HandRank evaluate_hand(const std::array<Card, 5>& hand) {
    std::array<int, 15> freq{};
    std::vector<int> ranks;
    ranks.reserve(5U);

    bool flush = true;
    const char suit0 = hand[0].suit;

    for (const Card& c : hand) {
        ++freq[static_cast<std::size_t>(c.rank)];
        ranks.push_back(c.rank);
        if (c.suit != suit0) {
            flush = false;
        }
    }

    std::sort(ranks.begin(), ranks.end());
    bool straight = true;
    for (int i = 1; i < 5; ++i) {
        if (ranks[static_cast<std::size_t>(i)] != ranks[static_cast<std::size_t>(i - 1)] + 1) {
            straight = false;
            break;
        }
    }

    std::vector<std::pair<int, int>> groups;
    for (int r = 2; r <= 14; ++r) {
        if (freq[static_cast<std::size_t>(r)] > 0) {
            groups.emplace_back(freq[static_cast<std::size_t>(r)], r);
        }
    }
    std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
            return a.first > b.first;
        }
        return a.second > b.second;
    });

    HandRank result;
    if (straight && flush) {
        result.category = 8;
        result.tiebreak = {ranks.back()};
        return result;
    }
    if (groups[0].first == 4) {
        result.category = 7;
        result.tiebreak = {groups[0].second, groups[1].second};
        return result;
    }
    if (groups[0].first == 3 && groups[1].first == 2) {
        result.category = 6;
        result.tiebreak = {groups[0].second, groups[1].second};
        return result;
    }
    if (flush) {
        result.category = 5;
        result.tiebreak.assign(ranks.rbegin(), ranks.rend());
        return result;
    }
    if (straight) {
        result.category = 4;
        result.tiebreak = {ranks.back()};
        return result;
    }
    if (groups[0].first == 3) {
        result.category = 3;
        result.tiebreak = {groups[0].second};
        for (std::size_t i = 1; i < groups.size(); ++i) {
            result.tiebreak.push_back(groups[i].second);
        }
        return result;
    }
    if (groups[0].first == 2 && groups[1].first == 2) {
        result.category = 2;
        result.tiebreak = {groups[0].second, groups[1].second, groups[2].second};
        return result;
    }
    if (groups[0].first == 2) {
        result.category = 1;
        result.tiebreak = {groups[0].second};
        for (std::size_t i = 1; i < groups.size(); ++i) {
            result.tiebreak.push_back(groups[i].second);
        }
        return result;
    }

    result.category = 0;
    result.tiebreak.assign(ranks.rbegin(), ranks.rend());
    return result;
}

bool player1_wins(const std::array<Card, 5>& p1, const std::array<Card, 5>& p2) {
    const HandRank h1 = evaluate_hand(p1);
    const HandRank h2 = evaluate_hand(p2);

    if (h1.category != h2.category) {
        return h1.category > h2.category;
    }

    return h1.tiebreak > h2.tiebreak;
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open poker file: " + file_path);
    }

    int wins = 0;
    std::string line;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::array<Card, 5> p1;
        std::array<Card, 5> p2;

        for (int i = 0; i < 5; ++i) {
            std::string token;
            iss >> token;
            p1[static_cast<std::size_t>(i)] = parse_card(token);
        }
        for (int i = 0; i < 5; ++i) {
            std::string token;
            iss >> token;
            p2[static_cast<std::size_t>(i)] = parse_card(token);
        }

        if (player1_wins(p1, p2)) {
            ++wins;
        }
    }

    return wins;
}

bool run_checkpoints() {
    {
        const std::array<Card, 5> p1 = {
            parse_card("5H"), parse_card("5C"), parse_card("6S"), parse_card("7S"), parse_card("KD")};
        const std::array<Card, 5> p2 = {
            parse_card("2C"), parse_card("3S"), parse_card("8S"), parse_card("8D"), parse_card("TD")};
        if (player1_wins(p1, p2)) {
            std::cerr << "Checkpoint failed for example hand 1" << '\n';
            return false;
        }
    }
    {
        const std::array<Card, 5> p1 = {
            parse_card("4D"), parse_card("6S"), parse_card("9H"), parse_card("QH"), parse_card("QC")};
        const std::array<Card, 5> p2 = {
            parse_card("3D"), parse_card("6D"), parse_card("7H"), parse_card("QD"), parse_card("QS")};
        if (!player1_wins(p1, p2)) {
            std::cerr << "Checkpoint failed for example hand 4" << '\n';
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
        std::cout << solve(options.file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
