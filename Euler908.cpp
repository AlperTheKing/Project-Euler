#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kMod = 1'111'211'113LL;

i64 mod_pow(i64 a, i64 e) {
    i64 r = 1 % kMod;
    i64 x = a % kMod;
    i64 p = e;
    while (p > 0) {
        if (p & 1) {
            r = static_cast<i64>((__int128)r * x % kMod);
        }
        x = static_cast<i64>((__int128)x * x % kMod);
        p >>= 1;
    }
    return r;
}

std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; (i64)p * p <= n; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[q] = false;
        }
    }
    std::vector<int> primes;
    for (int x = 2; x <= n; ++x) {
        if (is_prime[x]) {
            primes.push_back(x);
        }
    }
    return primes;
}

u64 triangular_residue_count_bruteforce(int s) {
    std::vector<char> seen(s, 0);
    int cnt = 0;
    for (int n = 1; n <= 2 * s; ++n) {
        const int r = (int)((1LL * n * (n + 1) / 2) % s);
        if (!seen[r]) {
            seen[r] = 1;
            ++cnt;
        }
    }
    return (u64)cnt;
}

u64 u_prime_power(int p, int e) {
    if (p == 2) {
        return u64{1} << e;
    }
    u64 u = (u64)(p + 1) / 2;
    for (int k = 2; k <= e; ++k) {
        const u64 d = (k % 2 == 0) ? (u64)(p - 1) : (u64)(p - 1) / 2;
        u = (u64)p * u - d;
    }
    return u;
}

u64 u_from_factorization(int s, const std::vector<int>& primes) {
    int x = s;
    u64 u = 1;
    for (int p : primes) {
        if ((i64)p * p > x) {
            break;
        }
        if (x % p != 0) {
            continue;
        }
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        u *= u_prime_power(p, e);
    }
    if (x > 1) {
        u *= u_prime_power(x, 1);
    }
    return u;
}

std::vector<int> mobius_array(int n) {
    std::vector<int> mu(n + 1, 0);
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);
    mu[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            const i64 v = 1LL * i * p;
            if (v > n || p > spf[i]) {
                break;
            }
            spf[(int)v] = p;
            if (i % p == 0) {
                mu[(int)v] = 0;
                break;
            }
            mu[(int)v] = -mu[i];
        }
    }
    return mu;
}

struct State908 {
    u64 integer = 0;
    int count = 0;
    int prime_index = 0;
};

std::vector<State908> enumerate_states(int n_limit) {
    const std::vector<int> primes = primes_up_to(2 * n_limit);
    std::vector<State908> states;
    states.reserve((std::size_t)n_limit * 20);

    for (u64 x = 1; x <= (u64)2 * n_limit; x <<= 1) {
        states.push_back(State908{x, (int)x, 1});
    }

    for (std::size_t head = 0; head < states.size(); ++head) {
        const State908 st = states[head];
        if (st.count > n_limit) {
            continue;
        }
        const int max_new_value_factor = n_limit / st.count;
        if (max_new_value_factor < 2) {
            continue;
        }
        const int max_prime = 2 * max_new_value_factor - 1;
        for (int new_prime_index = st.prime_index; new_prime_index < (int)primes.size(); ++new_prime_index) {
            const int prime = primes[new_prime_index];
            if (prime > max_prime) {
                break;
            }
            u64 new_integer = st.integer * (u64)prime;
            u64 prime_power = (u64)prime;
            int new_value_factor = (prime + 1) / 2;
            while (new_value_factor <= max_new_value_factor) {
                const int new_value = st.count * new_value_factor;
                states.push_back(State908{new_integer, new_value, new_prime_index + 1});

                const u128 next_prime_power = (u128)prime_power * (u64)prime;
                if (next_prime_power > (u128)std::numeric_limits<u64>::max()) {
                    break;
                }
                prime_power = (u64)next_prime_power;
                new_value_factor = (int)(((u64)prime * prime_power) / (2ULL * (u64)prime + 2ULL) + 1ULL);

                const u128 next_integer = (u128)new_integer * (u64)prime;
                if (next_integer > (u128)std::numeric_limits<u64>::max()) {
                    break;
                }
                new_integer = (u64)next_integer;
            }
        }
    }

    return states;
}

i64 count_clock_sequences(int n_limit) {
    const auto states = enumerate_states(n_limit);

    std::vector<i64> inv(n_limit + 2, 0);
    inv[1] = 1;
    for (int i = 2; i <= n_limit + 1; ++i) {
        inv[i] = (kMod - (kMod / i) * inv[kMod % i] % kMod) % kMod;
    }

    std::vector<i64> rep(n_limit + 1, 0);

    struct Task {
        const std::vector<State908>* states = nullptr;
        const std::vector<i64>* inv = nullptr;
        int n_limit = 0;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::vector<i64> local_rep;
    };

    auto worker = [](void* arg) -> void* {
        auto* task = static_cast<Task*>(arg);
        task->local_rep.assign((std::size_t)task->n_limit + 1, 0);
        const auto& states_ref = *task->states;
        const auto& inv_ref = *task->inv;
        const int n_limit_local = task->n_limit;

        for (std::size_t idx = task->begin; idx < task->end; ++idx) {
            const auto& st = states_ref[idx];
            if (st.count > n_limit_local) {
                continue;
            }
            const u64 free = st.integer - (u64)st.count;
            const int max_k = (int)std::min<u64>(free, (u64)(n_limit_local - st.count));

            i64 comb = 1;
            int p = st.count;
            for (int k = 0; k <= max_k; ++k) {
                i64& cell = task->local_rep[p];
                cell += comb;
                if (cell >= kMod) {
                    cell -= kMod;
                }
                if (k == max_k) {
                    break;
                }
                comb = static_cast<i64>((__int128)comb * (i64)(free - (u64)k) % kMod);
                comb = static_cast<i64>((__int128)comb * inv_ref[k + 1] % kMod);
                ++p;
            }
        }
        return nullptr;
    };

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    std::size_t thread_count = 1;
    if (cpu_count > 1) {
        thread_count = static_cast<std::size_t>(cpu_count);
        if (thread_count > 16) {
            thread_count = 16;
        }
        if (thread_count > states.size()) {
            thread_count = states.size();
        }
        if (thread_count * 256 > states.size()) {
            thread_count = std::max<std::size_t>(1, states.size() / 256);
        }
        if (thread_count == 0) {
            thread_count = 1;
        }
    }

    if (thread_count == 1) {
        Task task;
        task.states = &states;
        task.inv = &inv;
        task.n_limit = n_limit;
        task.begin = 0;
        task.end = states.size();
        worker(&task);
        rep = std::move(task.local_rep);
    } else {
        std::vector<pthread_t> threads(thread_count);
        std::vector<Task> tasks(thread_count);
        const std::size_t base = states.size() / thread_count;
        const std::size_t rem = states.size() % thread_count;
        std::size_t cur = 0;
        for (std::size_t i = 0; i < thread_count; ++i) {
            const std::size_t len = base + (i < rem ? 1 : 0);
            tasks[i].states = &states;
            tasks[i].inv = &inv;
            tasks[i].n_limit = n_limit;
            tasks[i].begin = cur;
            tasks[i].end = cur + len;
            cur += len;
            const int rc = ::pthread_create(&threads[i], nullptr, worker, &tasks[i]);
            assert(rc == 0);
        }
        for (std::size_t i = 0; i < thread_count; ++i) {
            const int rc = ::pthread_join(threads[i], nullptr);
            assert(rc == 0);
            const auto& local = tasks[i].local_rep;
            for (int p = 1; p <= n_limit; ++p) {
                rep[p] += local[p];
                if (rep[p] >= kMod) {
                    rep[p] -= kMod;
                }
            }
        }
    }

    for (int i = 2; i <= n_limit; ++i) {
        rep[i] += rep[i - 1];
        if (rep[i] >= kMod) {
            rep[i] -= kMod;
        }
    }

    const auto mu = mobius_array(n_limit);
    i64 result = 0;
    for (int i = 1; i <= n_limit; ++i) {
        result += (i64)mu[i] * rep[n_limit / i];
        result %= kMod;
    }
    if (result < 0) {
        result += kMod;
    }
    return result;
}

void validate() {
    const auto primes = primes_up_to(500);
    for (int s = 1; s <= 200; ++s) {
        assert(u_from_factorization(s, primes) == triangular_residue_count_bruteforce(s));
    }

    assert(count_clock_sequences(3) == 3);
    assert(count_clock_sequences(4) == 7);
    assert(count_clock_sequences(10) == 561);
}

}  // namespace

int main() {
    validate();
    std::cout << count_clock_sequences(10'000) << '\n';
    return 0;
}
