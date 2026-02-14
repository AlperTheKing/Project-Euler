#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>
#include <functional>

using i64 = std::int64_t;

static i64 binom(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    i64 res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - k + i) / i;
    }
    return res;
}

static long double expected_formula(int decks,
                                    int suits,
                                    int ranks,
                                    int jokers_per_deck,
                                    bool need_suits,
                                    bool need_ranks,
                                    bool need_decks) {
    const i64 total_cards = static_cast<i64>(decks) * (static_cast<i64>(suits) * ranks + jokers_per_deck);

    const int max_s = need_suits ? suits : 0;
    const int max_r = need_ranks ? ranks : 0;
    const int max_d = need_decks ? decks : 0;

    long double ans = 0.0L;

    for (int s = 0; s <= max_s; ++s) {
        const i64 cs = need_suits ? binom(suits, s) : 1;
        for (int r = 0; r <= max_r; ++r) {
            const i64 cr = need_ranks ? binom(ranks, r) : 1;
            for (int d = 0; d <= max_d; ++d) {
                const i64 cd = need_decks ? binom(decks, d) : 1;

                const int missed = s + r + d;
                if (missed == 0) {
                    continue;
                }

                const i64 safe_per_deck = static_cast<i64>(suits - s) * (ranks - r) + jokers_per_deck;
                const i64 safe_cards = static_cast<i64>(decks - d) * safe_per_deck;
                const i64 bad_cards = total_cards - safe_cards;

                const long double weight = static_cast<long double>(total_cards + 1) /
                                           static_cast<long double>(bad_cards + 1);
                const long double mult = static_cast<long double>(cs) *
                                         static_cast<long double>(cr) *
                                         static_cast<long double>(cd);

                if (missed & 1) {
                    ans += mult * weight;
                } else {
                    ans -= mult * weight;
                }
            }
        }
    }

    return ans;
}

struct TinyCard {
    int deck_bit;
    int suit_bit;
    int rank_bit;
};

static long double tiny_bruteforce_expectation() {
    std::vector<TinyCard> cards;
    for (int d = 0; d < 2; ++d) {
        for (int s = 0; s < 2; ++s) {
            for (int r = 0; r < 2; ++r) {
                cards.push_back({1 << d, 1 << s, 1 << r});
            }
        }
    }

    const int n = static_cast<int>(cards.size());
    const int full_deck = (1 << 2) - 1;
    const int full_suit = (1 << 2) - 1;
    const int full_rank = (1 << 2) - 1;

    std::vector<long double> memo(1 << n, -1.0L);

    std::function<long double(int)> dfs = [&](int mask) -> long double {
        long double& ref = memo[mask];
        if (ref >= 0.0L) {
            return ref;
        }

        int dmask = 0;
        int smask = 0;
        int rmask = 0;
        int drawn = 0;

        for (int i = 0; i < n; ++i) {
            if ((mask >> i) & 1) {
                ++drawn;
                dmask |= cards[i].deck_bit;
                smask |= cards[i].suit_bit;
                rmask |= cards[i].rank_bit;
            }
        }

        if (dmask == full_deck && smask == full_suit && rmask == full_rank) {
            ref = 0.0L;
            return ref;
        }

        const int rem = n - drawn;
        long double sum = 0.0L;
        for (int i = 0; i < n; ++i) {
            if (((mask >> i) & 1) == 0) {
                sum += dfs(mask | (1 << i));
            }
        }

        ref = 1.0L + sum / static_cast<long double>(rem);
        return ref;
    };

    return dfs(0);
}

int main() {
    const long double check_rank = expected_formula(1, 4, 13, 2, false, true, false);
    assert(std::fabsl(check_rank - 29.05361725L) < 1e-8L);

    const long double tiny_formula = expected_formula(2, 2, 2, 0, true, true, true);
    const long double tiny_exact = tiny_bruteforce_expectation();
    assert(std::fabsl(tiny_formula - tiny_exact) < 1e-12L);

    const long double ans = expected_formula(10, 4, 13, 2, true, true, true);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(ans) << '\n';
    return 0;
}
