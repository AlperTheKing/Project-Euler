#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;
using Polynomial = std::vector<i64>;
constexpr unsigned DEFAULT_LENGTH = 1000;
constexpr unsigned MAX_LENGTH = 100'000'000;
constexpr unsigned PREPERIOD = 53;
constexpr unsigned PERIOD = 34;
constexpr unsigned CHARACTERS = 16;
constexpr unsigned TWO_EXPONENT = 13;
constexpr unsigned FIVE_EXPONENT = 9;
constexpr u64 MOD = 1'000'000'000;
constexpr u64 SCALED_MOD = CHARACTERS * MOD;

constexpr u64 power(const u64 base, const unsigned exponent) {
    u64 result = 1;
    for (unsigned i = 0; i < exponent; ++i) result *= base;
    return result;
}

constexpr u64 TWO_MOD = power(2, TWO_EXPONENT);
constexpr u64 FIVE_MOD = power(5, FIVE_EXPONENT);
static_assert(TWO_MOD * FIVE_MOD == SCALED_MOD);

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

std::vector<u64> quadratic_character(const std::vector<unsigned>& grundy,
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
        const u64 twice_root = (square[n] + double_modulus
            - static_cast<u64>(products % double_modulus)) % double_modulus;
        require(twice_root % 2 == 0, "integral square-root coefficient");
        root[n] = twice_root / 2;
    }
    return root;
}

void trim(Polynomial& polynomial) {
    while (polynomial.size() > 1 && polynomial.back() == 0) polynomial.pop_back();
}

Polynomial product(const Polynomial& a, const Polynomial& b) {
    Polynomial result(a.size() + b.size() - 1);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) result[i + j] += a[i] * b[j];
    }
    trim(result);
    return result;
}

Polynomial difference(Polynomial a, const Polynomial& b) {
    a.resize(std::max(a.size(), b.size()));
    for (std::size_t i = 0; i < b.size(); ++i) a[i] -= b[i];
    trim(a);
    return a;
}

Polynomial derivative(const Polynomial& polynomial) {
    Polynomial result(std::max<std::size_t>(1, polynomial.size() - 1));
    for (std::size_t i = 1; i < polynomial.size(); ++i) result[i - 1] = static_cast<i64>(i) * polynomial[i];
    trim(result);
    return result;
}

std::vector<unsigned> certify_period() {
    constexpr unsigned END = 2 * PREPERIOD + PERIOD + 2;
    auto grundy = grundy_values(END + PERIOD - 1);
    // Guy-Smith: this finite window certifies the period for every n >= PREPERIOD.
    for (unsigned n = PREPERIOD; n < END; ++n) {
        require(grundy[n] == grundy[n + PERIOD], "Grundy periodicity certificate");
    }
    require(*std::max_element(grundy.begin(), grundy.end()) < CHARACTERS, "Walsh dimension");
    return grundy;
}

struct Term {
    unsigned lag;
    i64 constant, linear;
};

struct Recurrence {
    unsigned width;
    std::vector<Term> terms;
};

Recurrence make_recurrence(const unsigned character, const std::vector<unsigned>& grundy) {
    const auto sign = [&](unsigned n) {
        if (n >= PREPERIOD) n = PREPERIOD + (n - PREPERIOD) % PERIOD;
        return character_sign(character & grundy[n]);
    };
    constexpr unsigned HALF_PERIOD = PERIOD / 2;
    constexpr unsigned ODD_LENGTH = (PREPERIOD - 1) / 2 + HALF_PERIOD;
    constexpr unsigned EVEN_LENGTH = (PREPERIOD + 1) / 2 + HALF_PERIOD;
    Polynomial odd(ODD_LENGTH), even(EVEN_LENGTH);
    for (unsigned i = 0; i < ODD_LENGTH; ++i) {
        odd[i] = sign(2 * i + 1) - (i >= HALF_PERIOD ? sign(2 * (i - HALF_PERIOD) + 1) : 0);
    }
    for (unsigned i = 1; i < EVEN_LENGTH; ++i) {
        even[i] = sign(2 * i) - (i > HALF_PERIOD ? sign(2 * (i - HALF_PERIOD)) : 0);
    }
    Polynomial plus = even, minus = even;
    for (i64& value : minus) value = -value;
    ++plus[0];
    ++minus[0];
    --plus[HALF_PERIOD];
    --minus[HALF_PERIOD];
    auto odd_square = product(odd, odd);
    odd_square.insert(odd_square.begin(), 0);
    const auto p = difference(product(plus, plus), odd_square);
    const auto q = difference(product(minus, minus), odd_square);
    const auto t = product(p, q);
    auto u = difference(product(derivative(p), q), product(p, derivative(q)));
    // H^2 = P/Q gives (PQ)H' = (P'Q-PQ')H/2.
    for (i64& value : u) {
        require(value % 2 == 0, "integral differential equation");
        value /= 2;
    }
    require(t[0] == 1, "unit leading coefficient");
    const unsigned width = static_cast<unsigned>(std::max(t.size() - 1, u.size()));
    Recurrence recurrence{width, {}};
    for (unsigned lag = 1; lag <= width; ++lag) {
        const i64 tj = lag < t.size() ? t[lag] : 0;
        const i64 uj = lag <= u.size() ? u[lag - 1] : 0;
        if (tj != 0 || uj != 0) recurrence.terms.push_back({lag, uj + static_cast<i64>(lag) * tj, -tj});
    }
    return recurrence;
}

const std::vector<Recurrence>& recurrences() {
    static const auto result = [] {
        const auto grundy = certify_period();
        std::vector<Recurrence> values;
        for (unsigned h = 0; h < CHARACTERS; ++h) values.push_back(make_recurrence(h, grundy));
        return values;
    }();
    return result;
}

u64 inverse(const u64 value, const u64 modulus) {
    u64 remainder = modulus, next = value % modulus;
    i128 coefficient = 0, next_coefficient = 1;
    while (next != 0) {
        const u64 quotient = remainder / next;
        const u64 temporary = remainder % next;
        remainder = next;
        next = temporary;
        const i128 next_value = coefficient - static_cast<i128>(quotient) * next_coefficient;
        coefficient = next_coefficient;
        next_coefficient = next_value;
    }
    require(remainder == 1, "invertible unit");
    return residue(coefficient, modulus);
}

u64 working_modulus(const unsigned prime, const unsigned exponent, const unsigned limit) {
    unsigned logarithm = 0;
    for (unsigned n = limit; n >= prime; n /= prime) ++logarithm;
    // Integrating a residual modulo p^E loses at most L=floor(log_p n) digits.
    // E >= e+L preserves the target; E >= 2L permits every division.
    const unsigned digits = std::max(exponent + logarithm, 2 * logarithm);
    u64 modulus = 1;
    for (unsigned i = 0; i < digits; ++i) {
        require(modulus <= std::numeric_limits<u64>::max() / prime, "working modulus fits");
        modulus *= prime;
    }
    return modulus;
}

u64 divide_step(const i128 numerator, unsigned n, const unsigned prime, const u64 modulus) {
    u64 divisor = 1;
    while (n % prime == 0) {
        n /= prime;
        divisor *= prime;
    }
    const u64 value = residue(numerator, modulus);
    if (value % divisor != 0) throw std::runtime_error("Nonintegral prime-power division");
    const u64 reduced = modulus / divisor;
    return static_cast<u64>(static_cast<u128>(value / divisor) * inverse(n, reduced) % reduced);
}

u64 combine(const u64 two, const u64 five) {
    static const u64 multiplier = inverse(FIVE_MOD % TWO_MOD, TWO_MOD);
    const u64 b = five % FIVE_MOD;
    const u64 multiple = ((two % TWO_MOD + TWO_MOD - b % TWO_MOD) * multiplier) % TWO_MOD;
    return b + FIVE_MOD * multiple;
}

std::vector<u64> character_coefficients(const Recurrence& recurrence, const std::vector<unsigned>& indices) {
    const unsigned last = indices.back();
    const u64 mod2 = working_modulus(2, TWO_EXPONENT, last);
    const u64 mod5 = working_modulus(5, FIVE_EXPONENT, last);
    struct Residues { u64 two = 0, five = 0; };
    std::vector<Residues> history(2 * recurrence.width);
    history[0] = history[recurrence.width] = {1, 1};
    std::vector<u64> values;
    values.reserve(indices.size());
    std::size_t next = 0;
    if (indices[0] == 0) {
        values.push_back(1);
        ++next;
    }
    for (unsigned n = 1; n <= last; ++n) {
        const unsigned position = n % recurrence.width;
        const unsigned end = recurrence.width + position;
        i128 numerator2 = 0, numerator5 = 0;
        for (const auto& term : recurrence.terms) {
            if (term.lag > n) break;
            const i64 factor = term.constant + static_cast<i64>(n) * term.linear;
            const auto& previous = history[end - term.lag];
            numerator2 += static_cast<i128>(factor) * previous.two;
            numerator5 += static_cast<i128>(factor) * previous.five;
        }
        const Residues current{divide_step(numerator2, n, 2, mod2), divide_step(numerator5, n, 5, mod5)};
        history[position] = history[end] = current;
        if (next < indices.size() && n == indices[next]) {
            values.push_back(combine(current.two, current.five));
            ++next;
        }
    }
    return values;
}

struct CharacterWork {
    const std::vector<Recurrence>* recurrences;
    const std::vector<unsigned>* indices;
    std::vector<std::vector<u64>>* coefficients;
    unsigned first, stride;
    std::exception_ptr error;
};

void* calculate_characters(void* argument) {
    auto& work = *static_cast<CharacterWork*>(argument);
    try {
        for (unsigned h = work.first; h < CHARACTERS; h += work.stride) {
            (*work.coefficients)[h] = character_coefficients((*work.recurrences)[h], *work.indices);
        }
    } catch (...) {
        work.error = std::current_exception();
    }
    return nullptr;
}

std::vector<u64> fair_counts_at(const std::vector<unsigned>& indices, unsigned thread_count = 0) {
    if (indices.empty()) return {};
    const auto& equations = recurrences();
    if (thread_count == 0) {
        thread_count = indices.back() < 64 ? 1 : std::min(8U, std::max(1U, std::thread::hardware_concurrency()));
    }
    thread_count = std::min(thread_count, CHARACTERS);
    std::vector<std::vector<u64>> coefficients(CHARACTERS);
    std::vector<CharacterWork> work(thread_count);
    std::vector<pthread_t> threads(thread_count);
    std::vector<bool> started(thread_count);
    for (unsigned i = 0; i < thread_count; ++i) {
        work[i] = {&equations, &indices, &coefficients, i, thread_count, nullptr};
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
    std::vector<u64> counts(indices.size());
    for (const auto& character : coefficients) {
        for (std::size_t i = 0; i < counts.size(); ++i) counts[i] = (counts[i] + character[i]) % SCALED_MOD;
    }
    for (u64& count : counts) {
        require(count % CHARACTERS == 0, "integral Walsh average");
        count /= CHARACTERS;
    }
    return counts;
}

u64 solve(const unsigned length, const unsigned thread_count = 0) {
    if (length > MAX_LENGTH) throw std::invalid_argument("Supported length is at most 100000000");
    if (length == 0) return 1;
    if (length % 2 != 0) return 0;
    return fair_counts_at({length / 2}, thread_count)[0];
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

u64 central_binomial_residue(const unsigned n, const unsigned prime, const unsigned exponent) {
    const u64 modulus = power(prime, exponent);
    u64 numerator = 1, denominator = 1;
    int valuation = 0;
    for (unsigned i = 1; i <= n; ++i) {
        unsigned top = n + i, bottom = i;
        while (top % prime == 0) {
            top /= prime;
            ++valuation;
        }
        while (bottom % prime == 0) {
            bottom /= prime;
            --valuation;
        }
        numerator = numerator * top % modulus;
        denominator = denominator * bottom % modulus;
    }
    require(valuation >= 0, "integral binomial coefficient");
    if (static_cast<unsigned>(valuation) >= exponent) return 0;
    return numerator * inverse(denominator, modulus) % modulus
        * power(prime, static_cast<unsigned>(valuation)) % modulus;
}

void run_tests() {
    constexpr unsigned SMALL = 12;
    constexpr unsigned CHECK_LENGTH = 1000;
    std::vector<unsigned> indices(CHECK_LENGTH / 2 + 1);
    for (unsigned i = 0; i < indices.size(); ++i) indices[i] = i;
    const auto counts = fair_counts_at(indices);
    require(counts[0] == 1, "empty arrangement");
    require(counts[3] == 8, "given f(6)");
    require(counts[10] == 49'844, "given f(20)");
    require(fair_counts_at(indices, 1) == counts, "serial and pthread agreement");

    const auto grundy = grundy_values(CHECK_LENGTH);
    std::vector<u64> expected_counts(indices.size()), unsigned_counts;
    for (unsigned h = 0; h < CHARACTERS; ++h) {
        const auto expected = quadratic_character(grundy, h, SCALED_MOD);
        require(character_coefficients(recurrences()[h], indices) == expected,
                "every character coefficient agrees with quadratic extraction");
        for (std::size_t i = 0; i < indices.size(); ++i) {
            expected_counts[i] = (expected_counts[i] + expected[i]) % SCALED_MOD;
        }
        if (h == 0) unsigned_counts = expected;
    }
    for (u64& count : expected_counts) count /= CHARACTERS;
    require(counts == expected_counts, "all fair counts agree with quadratic extraction");
    for (unsigned n = PREPERIOD; n <= CHECK_LENGTH; ++n) {
        require(grundy[n] == grundy[PREPERIOD + (n - PREPERIOD) % PERIOD], "periodic Grundy lookup");
    }

    std::vector<unsigned> powers;
    const auto outcomes = brute_outcomes(SMALL, powers);
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
        require(solve(n, 1) == expected, "short-prefix consistency");
    }

    for (const unsigned n : {14U, 16U, 18U, 30U, 32U, 48U, 50U, 62U, 64U, 126U, 128U,
                            248U, 250U, 254U, 256U, 510U, 512U, 998U, 1000U}) {
        require(solve(n, 1) == counts[n / 2], "agreement across working-precision boundaries");
    }
    require(solve(MAX_LENGTH - 1) == 0, "large odd length");

    std::vector<u64> binomial(CHECK_LENGTH + 1);
    binomial[0] = 1;
    for (unsigned n = 1; n <= CHECK_LENGTH; ++n) {
        for (unsigned k = n; k > 0; --k) binomial[k] = (binomial[k] + binomial[k - 1]) % MOD;
        if (n % 2 == 0) require(unsigned_counts[n / 2] % MOD == binomial[n / 2], "central binomial coefficients");
    }
    const std::vector<unsigned> large_indices{
        8191, 8192, 8193, 15624, 15625, 15626, 16383, 16384, 16385
    };
    const auto large_binomials = character_coefficients(recurrences()[0], large_indices);
    for (std::size_t i = 0; i < large_indices.size(); ++i) {
        const unsigned n = large_indices[i];
        require(large_binomials[i] % TWO_MOD == central_binomial_residue(n, 2, TWO_EXPONENT),
                "binomial product across powers of two");
        require(large_binomials[i] % FIVE_MOD == central_binomial_residue(n, 5, FIVE_EXPONENT),
                "binomial product across powers of five");
    }
    std::cout << "All checks passed; " << outcomes.size()
              << " original-game positions, " << CHARACTERS * indices.size()
              << " character coefficients, periodicity, and prime-power precision.\n";
}

unsigned parse_length(const std::string& text) {
    unsigned length = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), length);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || length > MAX_LENGTH) {
        throw std::invalid_argument("Length must be an integer from 0 to 100000000");
    }
    return length;
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") run_tests();
        else {
            if (argc > 2) throw std::invalid_argument("Usage: Euler1015 [n | --self-test]");
            const unsigned length = argc == 2 ? parse_length(argv[1]) : DEFAULT_LENGTH;
            std::cout << std::setfill('0') << std::setw(9) << solve(length) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
