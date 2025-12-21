// Commit message: Euler905: reduce by sum-difference chain and unwind turn order.
// Approach: repeatedly replace the known-sum value with the absolute difference,
// record which role held the sum each step, then unwind turns to get F().
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

static uint64_t pow_int(uint64_t base, int exp) {
    uint64_t result = 1;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

static uint64_t next_turn_after(uint64_t t, int mod) {
    uint64_t candidate = t + 1;
    uint64_t r = candidate % 3;
    if (r == static_cast<uint64_t>(mod)) {
        return candidate;
    }
    uint64_t add = (static_cast<uint64_t>(mod) + 3 - r) % 3;
    return candidate + add;
}

static uint64_t compute_F(uint64_t A, uint64_t B, uint64_t C) {
    static constexpr int role_mod[3] = {1, 2, 0}; // A,B,C turns.
    static constexpr int role_base[3] = {1, 2, 3};
    std::vector<uint8_t> roles;
    int base_role = -1;

    // Reduce by replacing the sum with the difference, then unwind turns.
    while (true) {
        if (A == B + C) {
            roles.push_back(0);
            if (B == C) {
                base_role = 0;
                break;
            }
            A = (B >= C) ? (B - C) : (C - B);
        } else if (B == A + C) {
            roles.push_back(1);
            if (A == C) {
                base_role = 1;
                break;
            }
            B = (A >= C) ? (A - C) : (C - A);
        } else {
            roles.push_back(2);
            if (A == B) {
                base_role = 2;
                break;
            }
            C = (A >= B) ? (A - B) : (B - A);
        }
    }

    uint64_t turns = role_base[base_role];
    roles.pop_back();
    for (auto it = roles.rbegin(); it != roles.rend(); ++it) {
        turns = next_turn_after(turns, role_mod[*it]);
    }
    return turns;
}

static bool run_validation() {
    bool ok = true;
    if (compute_F(2, 1, 1) != 1) {
        std::cerr << "Validation failed: F(2,1,1) != 1\n";
        ok = false;
    }
    if (compute_F(2, 7, 5) != 5) {
        std::cerr << "Validation failed: F(2,7,5) != 5\n";
        ok = false;
    }
    if (compute_F(1, 1, 2) != 3) {
        std::cerr << "Validation failed: F(1,1,2) != 3\n";
        ok = false;
    }
    if (compute_F(1, 2, 3) != 3) {
        std::cerr << "Validation failed: F(1,2,3) != 3\n";
        ok = false;
    }
    if (ok) {
        std::cerr << "Validation checkpoints passed.\n";
    }
    return ok;
}

int main() {
    if (!run_validation()) {
        return 1;
    }

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 1;
    }

    std::vector<uint64_t> partial(threads, 0);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned int t = 0; t < threads; ++t) {
        pool.emplace_back([t, threads, &partial]() {
            uint64_t local = 0;
            for (int a = 1 + static_cast<int>(t); a <= 7; a += static_cast<int>(threads)) {
                for (int b = 1; b <= 19; ++b) {
                    uint64_t A = pow_int(static_cast<uint64_t>(a), b);
                    uint64_t B = pow_int(static_cast<uint64_t>(b), a);
                    uint64_t C = A + B;
                    local += compute_F(A, B, C);
                }
            }
            partial[t] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    uint64_t total = 0;
    for (uint64_t v : partial) {
        total += v;
    }

    std::cout << total << "\n";
    return 0;
}
