#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;
constexpr unsigned LIMIT = 1000;
constexpr u64 MOD = 1'000'000'000;

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

std::vector<unsigned> grundy_values(const unsigned limit) {
    std::vector<unsigned> grundy(limit + 1), seen(limit + 2);
    for (unsigned n = 2; n <= limit; ++n) {
        for (unsigned left = 0; left <= (n - 2) / 2; ++left) {
            const unsigned value = grundy[left] ^ grundy[n - 2 - left];
            if (value < seen.size()) seen[value] = n;
        }
        while (seen[grundy[n]] == n) ++grundy[n];
    }
    return grundy;
}

int character_sign(unsigned value) {
    int sign = 1;
    for (; value; value &= value - 1) sign = -sign;
    return sign;
}

u64 residue(i128 value, const u64 modulus) {
    value %= static_cast<i128>(modulus);
    if (value < 0) value += modulus;
    return static_cast<u64>(value);
}

std::vector<u64> balanced_character(const std::vector<unsigned>& grundy,
                                  const unsigned character, const u64 modulus) {
    const unsigned degree = static_cast<unsigned>(grundy.size() - 1) / 2;
    std::vector<int> odd(degree), even(degree + 1), numerator(degree + 1), denominator(degree + 1);
    for (unsigned n = 0; n < degree; ++n) odd[n] = character_sign(character & grundy[2 * n + 1]);
    for (unsigned n = 1; n <= degree; ++n) even[n] = character_sign(character & grundy[2 * n]);
    numerator[0] = denominator[0] = 1;
    for (unsigned n = 1; n <= degree; ++n) {
        int common = 0;
        for (unsigned i = 1; i < n; ++i) common += even[i] * even[n - i];
        for (unsigned i = 0; i < n; ++i) common -= odd[i] * odd[n - 1 - i];
        numerator[n] = common + 2 * even[n];
        denominator[n] = common - 2 * even[n];
    }

    // The balanced-row series is sqrt(((1+E)^2-A^2)/((1-E)^2-A^2)).
    const u64 double_modulus = 2 * modulus;
    std::vector<u64> square(degree + 1), root(degree + 1);
    square[0] = root[0] = 1;
    for (unsigned n = 1; n <= degree; ++n) {
        i128 coefficient = numerator[n];
        for (unsigned i = 1; i <= n; ++i) {
            coefficient -= static_cast<i128>(denominator[i]) * square[n - i];
        }
        square[n] = residue(coefficient, double_modulus);
        u128 products = 0;
        for (unsigned i = 1; 2 * i < n; ++i) products += u128{2} * root[i] * root[n - i];
        if (n % 2 == 0) products += static_cast<u128>(root[n / 2]) * root[n / 2];
        // For even M, this symmetric sum modulo 2M needs earlier roots only modulo M.
        const u64 twice_root = (square[n] + double_modulus
            - static_cast<u64>(products % double_modulus)) % double_modulus;
        require(twice_root % 2 == 0, "integral square-root coefficient");
        root[n] = twice_root / 2;
    }
    return root;
}

struct CharacterWork {
    const std::vector<unsigned>* grundy = nullptr;
    std::vector<std::vector<u64>>* coefficients = nullptr;
    u64 modulus = 0;
    unsigned first = 0, stride = 1;
    std::exception_ptr error;
};

void* calculate_characters(void* argument) {
    auto& work = *static_cast<CharacterWork*>(argument);
    try {
        for (unsigned h = work.first; h < work.coefficients->size(); h += work.stride) {
            (*work.coefficients)[h] = balanced_character(*work.grundy, h, work.modulus);
        }
    } catch (...) {
        work.error = std::current_exception();
    }
    return nullptr;
}

std::vector<u64> fair_counts(const unsigned limit, unsigned thread_count = 0) {
    if (limit > LIMIT) throw std::invalid_argument("Supported length is at most 1000");
    const auto grundy = grundy_values(limit);
    unsigned characters = 1;
    while (characters <= *std::max_element(grundy.begin(), grundy.end())) characters *= 2;
    const u64 modulus = MOD * characters;
    if (thread_count == 0) {
        thread_count = limit < 128 ? 1 : std::min(8U, std::max(1U, std::thread::hardware_concurrency()));
    }
    thread_count = std::min(thread_count, characters);
    std::vector<std::vector<u64>> coefficients(characters);
    std::vector<CharacterWork> work(thread_count);
    std::vector<pthread_t> threads(thread_count);
    std::vector<bool> started(thread_count);
    for (unsigned i = 0; i < thread_count; ++i) {
        work[i] = {&grundy, &coefficients, modulus, i, thread_count, nullptr};
    }
    for (unsigned i = 1; i < thread_count; ++i) {
        started[i] = pthread_create(&threads[i], nullptr, calculate_characters, &work[i]) == 0;
        if (!started[i]) calculate_characters(&work[i]);
    }
    calculate_characters(&work[0]);
    for (unsigned i = 1; i < thread_count; ++i) {
        if (started[i] && pthread_join(threads[i], nullptr) != 0) std::terminate();
    }
    for (const auto& job : work) {
        if (job.error) std::rethrow_exception(job.error);
    }
    std::vector<u64> counts(limit / 2 + 1);
    for (unsigned character = 0; character < characters; ++character) {
        for (unsigned n = 0; n < counts.size(); ++n) {
            counts[n] = (counts[n] + coefficients[character][n]) % modulus;
        }
    }
    for (u64& count : counts) {
        require(count % characters == 0, "integral Walsh average");
        count /= characters;
    }
    return counts;
}

std::vector<unsigned char> brute_outcomes(const unsigned length, std::vector<unsigned>& powers) {
    powers.assign(length + 1, 1);
    for (unsigned i = 1; i <= length; ++i) powers[i] = 3 * powers[i - 1];
    std::vector<unsigned char> outcomes(powers[length]);
    for (unsigned state = 1; state < outcomes.size(); ++state) {
        unsigned remaining = state;
        unsigned char result = 0;
        for (unsigned i = 0; i < length && result != 3; ++i) {
            const unsigned pin = remaining % 3;
            remaining /= 3;
            if (pin == 0) continue;
            const unsigned single = state - pin * powers[i];
            if (pin == 1 && !(outcomes[single] & 2)) result |= 1;
            if (pin == 2 && !(outcomes[single] & 1)) result |= 2;
            const unsigned next = remaining % 3;
            if (next == 0) continue;
            const unsigned pair = single - next * powers[i + 1];
            if ((pin == 1 || next == 1) && !(outcomes[pair] & 2)) result |= 1;
            if ((pin == 2 || next == 2) && !(outcomes[pair] & 1)) result |= 2;
        }
        outcomes[state] = result;
    }
    return outcomes;
}

unsigned predicted_outcome(unsigned state, const unsigned length, const std::vector<unsigned>& grundy) {
    int balance = 0;
    unsigned previous = 0, run = 0, nim_sum = 0;
    for (unsigned i = 0; i < length; ++i) {
        const unsigned pin = state % 3;
        state /= 3;
        if (pin == 1) ++balance;
        if (pin == 2) --balance;
        if (pin != 0 && previous != 0 && pin != previous) ++run;
        else {
            nim_sum ^= grundy[run];
            run = pin == 0 ? 0 : 1;
        }
        previous = pin;
    }
    nim_sum ^= grundy[run];
    if (balance > 0) return 1;
    if (balance < 0) return 2;
    return nim_sum == 0 ? 0 : 3;
}

void run_tests() {
    constexpr unsigned SMALL = 12;
    const auto counts = fair_counts(LIMIT);
    require(counts[0] == 1, "empty arrangement");
    require(counts[3] == 8, "given f(6)");
    require(counts[10] == 49'844, "given f(20)");
    require(fair_counts(LIMIT, 1) == counts, "serial and pthread agreement");

    std::vector<unsigned> powers;
    const auto outcomes = brute_outcomes(SMALL, powers);
    const auto grundy = grundy_values(LIMIT);
    for (unsigned state = 0; state < outcomes.size(); ++state) {
        require(outcomes[state] == predicted_outcome(state, SMALL, grundy),
                "original-game minimax agrees with balance and alternating blocks");
    }
    for (unsigned n = 0; n <= SMALL; ++n) {
        u64 expected = 0;
        for (unsigned word = 0; word < (1U << n); ++word) {
            unsigned state = 0;
            for (unsigned i = 0; i < n; ++i) state += (1 + ((word >> i) & 1U)) * powers[i];
            expected += outcomes[state] == 0;
        }
        require(expected == (n % 2 == 0 ? counts[n / 2] : 0), "exhaustive initial arrangements");
        require(fair_counts(n).back() == counts[n / 2], "short-prefix consistency");
    }

    std::vector<u64> binomial(LIMIT + 1);
    binomial[0] = 1;
    const auto unsigned_counts = balanced_character(grundy, 0, MOD);
    for (unsigned n = 1; n <= LIMIT; ++n) {
        for (unsigned k = n; k > 0; --k) binomial[k] = (binomial[k] + binomial[k - 1]) % MOD;
        if (n % 2 == 0) require(unsigned_counts[n / 2] == binomial[n / 2], "central binomial coefficients");
    }
    std::cout << "All checks passed; " << outcomes.size()
              << " original-game positions, given examples, and central binomial coefficients through " << LIMIT << ".\n";
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") run_tests();
        else {
            if (argc != 1) throw std::invalid_argument("Usage: Euler1015 [--self-test]");
            std::cout << std::setfill('0') << std::setw(9) << fair_counts(LIMIT).back() << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
