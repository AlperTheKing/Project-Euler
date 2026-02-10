#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Project Euler 610: build a Markov chain over minimal Roman numerals for 1..999 (thousands handled by
// a one-line equation for the leading M-run) and solve expected value by backward substitution.

static int roman_value(char c) {
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

static int roman_to_int(const std::string& s) {
    int total = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int v = roman_value(s[i]);
        if (i + 1 < s.size() && roman_value(s[i + 1]) > v) total -= v;
        else total += v;
    }
    return total;
}

static std::string int_to_minimal_roman(int value) {
    if (value <= 0) return "";
    static const std::pair<int, const char*> table[] = {
        {1000, "M"}, {900, "CM"}, {500, "D"}, {400, "CD"}, {100, "C"}, {90, "XC"}, {50, "L"},
        {40, "XL"},  {10, "X"},   {9, "IX"},  {5, "V"},   {4, "IV"},  {1, "I"},
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

int main() {
    constexpr long double p_stop = 0.02L;
    constexpr long double p_letter = 0.14L;

    const char letters[7] = {'I', 'V', 'X', 'L', 'C', 'D', 'M'};

    std::vector<std::string> minr(1000);
    for (int v = 0; v <= 999; ++v) minr[v] = int_to_minimal_roman(v);

    // transitions r -> r' for appending a letter, or r itself if rejected
    std::vector<std::array<int, 7>> nxt(1000);
    for (int r = 1; r <= 999; ++r) {
        const std::string& s = minr[r];
        for (int i = 0; i < 7; ++i) {
            std::string cand = s;
            cand.push_back(letters[i]);
            const int v = roman_to_int(cand);
            if (0 < v && v <= 999 && minr[v] == cand) nxt[r][i] = v;
            else nxt[r][i] = r;
        }
    }

    std::vector<long double> E(1000, 0.0L);
    for (int r = 999; r >= 1; --r) {
        int rejected = 0;
        long double sum_next = 0.0L;
        for (int i = 0; i < 7; ++i) {
            const int v = nxt[r][i];
            if (v == r) ++rejected;
            else sum_next += E[v];
        }
        const long double denom = 1.0L - p_letter * (long double)rejected;
        E[r] = (p_stop * (long double)r + p_letter * sum_next) / denom;
    }

    // Leading M-run: E0 satisfies 0.86*E0 = 140 + 0.14*sum(E[value(letter)]), letter in {I,V,X,L,C,D}.
    const long double S_other = E[1] + E[5] + E[10] + E[50] + E[100] + E[500];
    const long double E0 = (140.0L + 0.14L * S_other) / 0.86L;

    std::cout << std::fixed << std::setprecision(8) << (double)E0 << "\n";
    return 0;
}
