#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

using LL = long long;
using VI = std::vector<int>;

std::string source = "thereisasyetinsufficientdataforameaningfulanswer";

std::map<std::pair<VI, int>, LL> dp;

LL count_subtree(VI counts, const int remaining) {
    if (remaining == 0) {
        return 1;
    }

    std::sort(counts.begin(), counts.end(), std::greater<int>());
    while (!counts.empty() && counts.back() == 0) {
        counts.pop_back();
    }

    auto key = std::make_pair(std::move(counts), remaining);
    auto it = dp.find(key);
    if (it != dp.end()) {
        return it->second;
    }

    LL total = 1;
    for (std::size_t i = 0; i < key.first.size(); ++i) {
        if (key.first[i] == 0) {
            continue;
        }
        VI next = key.first;
        --next[i];
        total += count_subtree(std::move(next), remaining - 1);
    }

    dp[std::move(key)] = total;
    return total;
}

std::string build_word(VI counts, const int remaining, LL pos) {
    if (pos == 1) {
        return "";
    }

    --pos;
    for (int c = 0; c < 26; ++c) {
        if (counts[c] == 0) {
            continue;
        }
        --counts[c];
        const LL cnt = count_subtree(counts, remaining - 1);
        if (cnt < pos) {
            ++counts[c];
            pos -= cnt;
            continue;
        }
        return std::string(1, static_cast<char>('a' + c)) + build_word(std::move(counts), remaining - 1, pos);
    }

    return "";
}

LL position_of(const VI& counts, const int remaining, const std::string_view target) {
    if (target.empty()) {
        return 1;
    }

    const int first = target.front() - 'a';
    LL pos = 1;

    for (int c = 0; c < first; ++c) {
        if (counts[c] == 0) {
            continue;
        }
        VI next = counts;
        --next[c];
        pos += count_subtree(next, remaining - 1);
    }

    VI next = counts;
    --next[first];
    pos += position_of(next, remaining - 1, target.substr(1));
    return pos;
}

int main() {
    VI counts(26, 0);
    for (const char ch : source) {
        ++counts[static_cast<std::size_t>(ch - 'a')];
    }

    const auto position_of_word = [&](const std::string_view word) {
        return position_of(counts, 15, word) - 1;
    };

    const auto word_at = [&](const LL rank) {
        return build_word(counts, 15, rank + 1);
    };

    const LL target =
        position_of_word("legionary") +
        position_of_word("calorimeters") -
        position_of_word("annihilate") +
        position_of_word("orchestrated") -
        position_of_word("fluttering");

    std::cout << word_at(target) << '\n';
    return 0;
}
