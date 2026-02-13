#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

class PrimeStream {
public:
    PrimeStream() : primes_{2}, next_candidate_(3), first_(true) {}

    u32 next() {
        if (first_) {
            first_ = false;
            return 2;
        }
        while (true) {
            const u32 candidate = next_candidate_;
            next_candidate_ += 2;
            bool is_prime = true;
            for (u32 p : primes_) {
                if (static_cast<u64>(p) * p > candidate) {
                    break;
                }
                if (candidate % p == 0) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime) {
                primes_.push_back(candidate);
                return candidate;
            }
        }
    }

private:
    std::vector<u32> primes_;
    u32 next_candidate_;
    bool first_;
};

class AlphaGenerator {
public:
    AlphaGenerator() : emitted_a0_(false), ones_left_(0), need_two_(false), current_prime_(0) {}

    u32 next() {
        if (!emitted_a0_) {
            emitted_a0_ = true;
            current_prime_ = primes_.next();
            ones_left_ = current_prime_;
            need_two_ = false;
            return 2;
        }

        if (ones_left_ > 0) {
            --ones_left_;
            return 1;
        }

        if (!need_two_) {
            need_two_ = true;
            return 2;
        }

        current_prime_ = primes_.next();
        ones_left_ = current_prime_ - 1;
        need_two_ = false;
        return 1;
    }

private:
    PrimeStream primes_;
    bool emitted_a0_;
    u32 ones_left_;
    bool need_two_;
    u32 current_prime_;
};

struct HomographicState {
    u64 p = 2;
    u64 q = 3;
    u64 r = 3;
    u64 s = 2;

    bool can_emit() const {
        return r != 0 && s != 0 && (p / r) == (q / s);
    }

    u32 emit() {
        const u64 a = p / r;
        const u64 np = r;
        const u64 nq = s;
        const u64 nr = p - a * r;
        const u64 ns = q - a * s;
        p = np;
        q = nq;
        r = nr;
        s = ns;
        return static_cast<u32>(a);
    }

    void consume(u32 n) {
        const u64 np = p * n + q;
        const u64 nq = p;
        const u64 nr = r * n + s;
        const u64 ns = r;
        p = np;
        q = nq;
        r = nr;
        s = ns;
    }
};

std::vector<u32> first_coefficients(std::size_t count) {
    std::vector<u32> out;
    out.reserve(count);

    AlphaGenerator alpha;
    HomographicState st;

    while (out.size() < count) {
        if (st.can_emit()) {
            out.push_back(st.emit());
        } else {
            st.consume(alpha.next());
        }
    }

    return out;
}

u64 sum_first_coefficients(std::size_t count) {
    u64 sum = 0;
    std::size_t emitted = 0;

    AlphaGenerator alpha;
    HomographicState st;

    while (emitted < count) {
        if (st.can_emit()) {
            sum += st.emit();
            ++emitted;
        } else {
            st.consume(alpha.next());
        }
    }

    return sum;
}

void run_validations() {
    const std::array<u32, 10> expected = {0, 1, 5, 6, 16, 9, 1, 10, 16, 11};
    const auto got = first_coefficients(10);
    assert(std::equal(got.begin(), got.end(), expected.begin(), expected.end()));
    assert(std::accumulate(got.begin(), got.end(), 0ULL) == 75ULL);
    assert(sum_first_coefficients(10) == 75ULL);
}

}  // namespace

int main() {
    run_validations();
    constexpr std::size_t kCount = 100'000'000ULL;
    std::cout << sum_first_coefficients(kCount) << '\n';
    return 0;
}
