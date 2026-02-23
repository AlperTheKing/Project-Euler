#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <unordered_set>
#include <vector>

namespace {

constexpr int kMaxN = 12;
constexpr int kLenBits = 4;
constexpr int kDigitBits = 5;
constexpr int kDigitMask = (1 << kDigitBits) - 1;

struct CycleData {
    int n = 0;
    std::array<int, kMaxN> d{};
};

struct ExpandWorker {
    const std::vector<std::uint64_t>* keys = nullptr;
    std::atomic<std::size_t>* next_idx = nullptr;
    int max_n = 0;
    std::array<std::vector<std::uint64_t>, kMaxN + 1> produced;
};

std::uint64_t encode_rotated(const std::array<int, kMaxN>& seq, int n, int start) {
    std::uint64_t key = static_cast<std::uint64_t>(n);
    for (int i = 0; i < n; ++i) {
        const int v = seq[static_cast<std::size_t>((start + i) % n)];
        key |= (static_cast<std::uint64_t>(v) << (kLenBits + kDigitBits * i));
    }
    return key;
}

std::uint64_t canonicalize_cycle(const std::array<int, kMaxN>& seq, int n) {
    int best = 0;
    for (int s = 1; s < n; ++s) {
        for (int k = 0; k < n; ++k) {
            const int a = seq[static_cast<std::size_t>((s + k) % n)];
            const int b = seq[static_cast<std::size_t>((best + k) % n)];
            if (a < b) {
                best = s;
                break;
            }
            if (a > b) {
                break;
            }
        }
    }
    return encode_rotated(seq, n, best);
}

CycleData decode_cycle(std::uint64_t key) {
    CycleData out;
    out.n = static_cast<int>(key & ((1u << kLenBits) - 1));
    for (int i = 0; i < out.n; ++i) {
        out.d[static_cast<std::size_t>(i)] =
            static_cast<int>((key >> (kLenBits + kDigitBits * i)) & kDigitMask);
    }
    return out;
}

bool is_trace_complex(const CycleData& c) {
    std::int64_t A = 1;
    std::int64_t B = 0;
    std::int64_t C = 0;
    std::int64_t D = 1;
    for (int i = 0; i < c.n; ++i) {
        const int a = c.d[static_cast<std::size_t>(i)];
        const __int128 nextA = static_cast<__int128>(A) * a + B;
        const __int128 nextB = -static_cast<__int128>(A);
        const __int128 nextC = static_cast<__int128>(C) * a + D;
        const __int128 nextD = -static_cast<__int128>(C);
        A = static_cast<std::int64_t>(nextA);
        B = static_cast<std::int64_t>(nextB);
        C = static_cast<std::int64_t>(nextC);
        D = static_cast<std::int64_t>(nextD);
    }
    const std::int64_t trace = A + D;
    return trace >= -1 && trace <= 1;
}

bool is_primitive_cycle(const CycleData& c) {
    for (int d = 1; d < c.n; ++d) {
        if (c.n % d != 0) {
            continue;
        }
        bool periodic = true;
        for (int i = 0; i < c.n; ++i) {
            if (c.d[static_cast<std::size_t>(i)] != c.d[static_cast<std::size_t>((i + d) % c.n)]) {
                periodic = false;
                break;
            }
        }
        if (periodic) {
            return false;
        }
    }
    return true;
}

void emit_expand1(const CycleData& c, std::vector<std::uint64_t>& out) {
    if (c.n + 1 > kMaxN) {
        return;
    }
    for (int i = 0; i < c.n; ++i) {
        const int j = (i + 1) % c.n;
        CycleData t;
        t.n = c.n + 1;
        int pos = 0;
        t.d[static_cast<std::size_t>(pos++)] = c.d[static_cast<std::size_t>(i)] + 1;
        t.d[static_cast<std::size_t>(pos++)] = 1;
        t.d[static_cast<std::size_t>(pos++)] = c.d[static_cast<std::size_t>(j)] + 1;
        int k = (j + 1) % c.n;
        while (k != i) {
            t.d[static_cast<std::size_t>(pos++)] = c.d[static_cast<std::size_t>(k)];
            k = (k + 1) % c.n;
        }
        if (is_trace_complex(t)) {
            out.push_back(canonicalize_cycle(t.d, t.n));
        }
    }
}

void emit_expand0(const CycleData& c, std::vector<std::uint64_t>& out) {
    if (c.n + 2 > kMaxN) {
        return;
    }
    for (int i = 0; i < c.n; ++i) {
        const int v = c.d[static_cast<std::size_t>(i)];
        for (int a = 0; a <= v; ++a) {
            const int b = v - a;
            CycleData t;
            t.n = c.n + 2;
            int pos = 0;
            t.d[static_cast<std::size_t>(pos++)] = a;
            t.d[static_cast<std::size_t>(pos++)] = 0;
            t.d[static_cast<std::size_t>(pos++)] = b;
            int k = (i + 1) % c.n;
            while (k != i) {
                t.d[static_cast<std::size_t>(pos++)] = c.d[static_cast<std::size_t>(k)];
                k = (k + 1) % c.n;
            }
            if (is_trace_complex(t)) {
                out.push_back(canonicalize_cycle(t.d, t.n));
            }
        }
    }
}

void* expand_worker_entry(void* raw) {
    auto* worker = static_cast<ExpandWorker*>(raw);
    while (true) {
        const std::size_t idx = worker->next_idx->fetch_add(1, std::memory_order_relaxed);
        if (idx >= worker->keys->size()) {
            break;
        }
        const CycleData c = decode_cycle((*worker->keys)[idx]);
        emit_expand1(c, worker->produced[static_cast<std::size_t>(c.n + 1)]);
        emit_expand0(c, worker->produced[static_cast<std::size_t>(c.n + 2)]);
    }
    return nullptr;
}

std::array<std::unordered_set<std::uint64_t>, kMaxN + 1> generate_classes(int max_n) {
    std::array<std::unordered_set<std::uint64_t>, kMaxN + 1> classes;
    const std::vector<std::vector<int>> seeds = {
        {0}, {1}, {1, 1}, {1, 2}, {2, 1}, {1, 3}, {3, 1},
    };

    for (const auto& s : seeds) {
        CycleData c;
        c.n = static_cast<int>(s.size());
        for (int i = 0; i < c.n; ++i) {
            c.d[static_cast<std::size_t>(i)] = s[static_cast<std::size_t>(i)];
        }
        if (c.n <= max_n && is_trace_complex(c)) {
            classes[static_cast<std::size_t>(c.n)].insert(canonicalize_cycle(c.d, c.n));
        }
    }

    for (int n = 1; n <= max_n; ++n) {
        if (classes[static_cast<std::size_t>(n)].empty()) {
            continue;
        }

        std::vector<std::uint64_t> keys;
        keys.reserve(classes[static_cast<std::size_t>(n)].size());
        for (const std::uint64_t key : classes[static_cast<std::size_t>(n)]) {
            keys.push_back(key);
        }

        long cpu = ::sysconf(_SC_NPROCESSORS_ONLN);
        std::size_t thread_count = 1;
        if (cpu > 0) {
            thread_count = static_cast<std::size_t>(cpu);
        }
        thread_count = std::max<std::size_t>(1, std::min<std::size_t>(thread_count, keys.size()));

        std::atomic<std::size_t> next_idx{0};
        std::vector<pthread_t> threads(thread_count);
        std::vector<ExpandWorker> workers(thread_count);
        for (std::size_t i = 0; i < thread_count; ++i) {
            workers[i].keys = &keys;
            workers[i].next_idx = &next_idx;
            workers[i].max_n = max_n;
            pthread_create(&threads[i], nullptr, expand_worker_entry, &workers[i]);
        }
        for (std::size_t i = 0; i < thread_count; ++i) {
            pthread_join(threads[i], nullptr);
        }

        for (std::size_t t = 0; t < thread_count; ++t) {
            for (int len = n + 1; len <= std::min(max_n, n + 2); ++len) {
                auto& dst = classes[static_cast<std::size_t>(len)];
                auto& src = workers[t].produced[static_cast<std::size_t>(len)];
                if (!src.empty()) {
                    dst.insert(src.begin(), src.end());
                }
            }
        }
    }

    return classes;
}

std::uint64_t brute_exact_period_small(int m) {
    const int cap = m + 1;
    std::array<int, kMaxN> seq{};
    std::uint64_t cnt = 0;

    auto is_primitive = [&](int n) {
        for (int d = 1; d < n; ++d) {
            if (n % d != 0) {
                continue;
            }
            bool periodic = true;
            for (int i = 0; i < n; ++i) {
                if (seq[static_cast<std::size_t>(i)] != seq[static_cast<std::size_t>(i % d)]) {
                    periodic = false;
                    break;
                }
            }
            if (periodic) {
                return false;
            }
        }
        return true;
    };

    auto dfs = [&](auto&& self, int depth, std::int64_t A, std::int64_t B, std::int64_t C,
                   std::int64_t D) -> void {
        if (depth == m) {
            const std::int64_t trace = A + D;
            if (trace >= -1 && trace <= 1 && is_primitive(m)) {
                ++cnt;
            }
            return;
        }
        for (int a = 0; a <= cap; ++a) {
            seq[static_cast<std::size_t>(depth)] = a;
            const __int128 nextA = static_cast<__int128>(A) * a + B;
            const __int128 nextB = -static_cast<__int128>(A);
            const __int128 nextC = static_cast<__int128>(C) * a + D;
            const __int128 nextD = -static_cast<__int128>(C);
            self(self, depth + 1, static_cast<std::int64_t>(nextA), static_cast<std::int64_t>(nextB),
                 static_cast<std::int64_t>(nextC), static_cast<std::int64_t>(nextD));
        }
    };

    dfs(dfs, 0, 1, 0, 0, 1);
    return cnt;
}

}  // namespace

int main() {
    const auto classes = generate_classes(kMaxN);

    std::array<std::uint64_t, kMaxN + 1> exact{};
    for (int n = 1; n <= kMaxN; ++n) {
        std::uint64_t class_cnt = 0;
        for (const std::uint64_t key : classes[static_cast<std::size_t>(n)]) {
            const CycleData c = decode_cycle(key);
            if (is_primitive_cycle(c)) {
                ++class_cnt;
            }
        }
        exact[static_cast<std::size_t>(n)] = class_cnt * static_cast<std::uint64_t>(n);
    }

    const std::uint64_t q1 = exact[1];
    const std::uint64_t q2 = exact[1] + exact[2];
    if (q1 != 2) {
        std::cerr << "Validation failed: Q(1) = " << q1 << ", expected 2\n";
        return 1;
    }
    if (q2 != 6) {
        std::cerr << "Validation failed: Q(2) = " << q2 << ", expected 6\n";
        return 1;
    }
    if (exact[2] != 4) {
        std::cerr << "Validation failed: exact period 2 = " << exact[2] << ", expected 4\n";
        return 1;
    }

    for (int n = 1; n <= 7; ++n) {
        const std::uint64_t brute = brute_exact_period_small(n);
        if (brute != exact[static_cast<std::size_t>(n)]) {
            std::cerr << "Validation failed: exact period " << n << " = "
                      << exact[static_cast<std::size_t>(n)] << ", brute force gives " << brute << '\n';
            return 1;
        }
    }

    std::uint64_t q12 = 0;
    for (int n = 1; n <= kMaxN; ++n) {
        q12 += exact[static_cast<std::size_t>(n)];
    }
    std::cout << q12 << '\n';
    return 0;
}
