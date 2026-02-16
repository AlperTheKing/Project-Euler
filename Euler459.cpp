#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

using namespace std;

namespace {

struct SeqResult {
    vector<long long> counts;
    int total_xor = 0;
    bool overflow = false;
};

static vector<int> build_triangular(int n) {
    vector<int> res;
    for (long long k = 1;; ++k) {
        long long t = k * (k + 1) / 2;
        if (t > n) {
            break;
        }
        res.push_back(static_cast<int>(t));
    }
    return res;
}

static vector<int> build_squares(int n) {
    vector<int> res;
    for (long long k = 1; k * k <= n; ++k) {
        res.push_back(static_cast<int>(k * k));
    }
    return res;
}

// Nimber multiplication using Fermat 2-power decomposition.
struct NimMul {
    vector<int> fermat;
    unordered_map<uint64_t, int> memo;

    NimMul() {
        fermat = {1, 2, 4, 16, 256, 65536};
        memo.reserve(1 << 20);
    }

    int largest_fermat(int x) const {
        for (int i = static_cast<int>(fermat.size()) - 1; i >= 0; --i) {
            if (fermat[i] <= x) {
                return fermat[i];
            }
        }
        return 1;
    }

    int mul(int a, int b) {
        if (a < b) {
            std::swap(a, b);
        }
        uint64_t key = (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
        auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        int res = 0;
        if (a == 0 || b == 0) {
            res = 0;
        } else if (a == 1) {
            res = b;
        } else if (b == 1) {
            res = a;
        } else {
            int n = largest_fermat(a);
            int a1 = a / n;
            int a0 = a % n;
            int b1 = b / n;
            int b0 = b % n;

            int p = mul(a0, b0);
            int q = mul(a0, b1) ^ mul(a1, b0);
            int r = mul(a1, b1);

            int shift = __builtin_ctz(n);
            res = p ^ ((q ^ r) << shift) ^ mul(r, n >> 1);
        }

        memo[key] = res;
        return res;
    }
};

static int nim_pow(int base, int exp, NimMul& nim) {
    int res = 1;
    int cur = base;
    int e = exp;
    while (e > 0) {
        if (e & 1) {
            res = nim.mul(res, cur);
        }
        e >>= 1;
        if (e > 0) {
            cur = nim.mul(cur, cur);
        }
    }
    return res;
}

static SeqResult compute_sequence(int n, const vector<int>& lengths, int limit, bool validate) {
    vector<int> prefix(n + 1, 0);
    vector<long long> counts(limit, 0);
    vector<uint16_t> freq(limit, 0);
    vector<uint16_t> used(limit, 0);

    int* const prefix_data = prefix.data();
    long long* const counts_data = counts.data();
    uint16_t* const freq_data = freq.data();
    uint16_t* const used_data = used.data();
    const int* const lengths_data = lengths.data();
    const int lengths_size = static_cast<int>(lengths.size());

    int len_count = 0;
    for (int i = 1; i <= n; ++i) {
        while (len_count < lengths_size && lengths_data[len_count] <= i) {
            ++len_count;
        }

        const int p_prev = prefix_data[i - 1];
        int used_len = 0;
        for (int idx = 0; idx < len_count; ++idx) {
            const int l = lengths_data[idx];
            const int val = p_prev ^ prefix_data[i - l];
            if (val >= limit) {
                SeqResult overflow_res;
                overflow_res.overflow = true;
                return overflow_res;
            }
            if (freq_data[val] == 0) {
                used_data[used_len++] = static_cast<uint16_t>(val);
            }
            ++freq_data[val];
        }

        int g = 0;
        while (g < limit && freq_data[g] > 0) {
            ++g;
        }
        if (g >= limit) {
            SeqResult overflow_res;
            overflow_res.overflow = true;
            return overflow_res;
        }

        prefix_data[i] = p_prev ^ g;
        for (int k = 0; k < used_len; ++k) {
            const int val = static_cast<int>(used_data[k]);
            counts_data[val ^ g] += freq_data[val];
            freq_data[val] = 0;
        }
    }

    if (validate) {
        long long total_segments = 0;
        for (int l : lengths) {
            total_segments += static_cast<long long>(n - l + 1);
        }
        long long sum_counts = 0;
        for (long long c : counts) {
            sum_counts += c;
        }
        assert(sum_counts == total_segments);
    }

    SeqResult res;
    res.counts = std::move(counts);
    res.total_xor = prefix_data[n];
    return res;
}

struct SequenceTask {
    int n = 0;
    const vector<int>* lengths = nullptr;
    int limit = 0;
    bool validate = false;
    SeqResult* out = nullptr;
};

static void* compute_sequence_worker(void* raw) {
    auto* task = static_cast<SequenceTask*>(raw);
    *task->out = compute_sequence(task->n, *task->lengths, task->limit, task->validate);
    return nullptr;
}

static long long solve(int n, bool validate) {
    vector<int> heights = build_triangular(n);
    vector<int> widths = build_squares(n);
    int limit = 512;

    SeqResult rows;
    SeqResult cols;
    while (true) {
        if (n >= 50000) {
            SequenceTask row_task{n, &heights, limit, validate, &rows};
            pthread_t thread_id{};
            const int create_rc =
                pthread_create(&thread_id, nullptr, compute_sequence_worker, &row_task);
            if (create_rc == 0) {
                cols = compute_sequence(n, widths, limit, validate);
                pthread_join(thread_id, nullptr);
            } else {
                rows = compute_sequence(n, heights, limit, validate);
                cols = compute_sequence(n, widths, limit, validate);
            }
        } else {
            rows = compute_sequence(n, heights, limit, validate);
            cols = compute_sequence(n, widths, limit, validate);
        }

        if (!rows.overflow && !cols.overflow) {
            break;
        }
        if (limit >= (1 << 20)) {
            cerr << "Exceeded maximum nimber limit while computing sequences\n";
            std::exit(1);
        }
        limit <<= 1;
    }

    NimMul nim;

    int total_xor = nim.mul(rows.total_xor, cols.total_xor);
    int max_value = total_xor;
    for (int i = 0; i < limit; ++i) {
        if (rows.counts[i] != 0 || cols.counts[i] != 0) {
            max_value = max(max_value, i);
        }
    }

    int field_size = 2;
    while (field_size <= max_value) {
        long long next = 1LL * field_size * field_size;
        field_size = static_cast<int>(next);
    }
    int inv_exp = field_size - 2;

    vector<int> inv(limit, -1);
    inv[0] = 0;
    if (limit > 1) {
        inv[1] = 1;
    }

    long long sum_rows = 0;
    long long sum_cols = 0;
    for (long long c : rows.counts) {
        sum_rows += c;
    }
    for (long long c : cols.counts) {
        sum_cols += c;
    }

    long long answer = 0;
    if (total_xor == 0) {
        answer = rows.counts[0] * sum_cols + cols.counts[0] * sum_rows
            - rows.counts[0] * cols.counts[0];
    } else {
        for (int a = 1; a < limit; ++a) {
            long long count_a = rows.counts[a];
            if (count_a == 0) {
                continue;
            }
            if (inv[a] == -1) {
                inv[a] = nim_pow(a, inv_exp, nim);
            }
            int b = nim.mul(inv[a], total_xor);
            if (b >= 0 && b < limit) {
                answer += count_a * cols.counts[b];
            }
        }
    }

    return answer;
}

static void run_validations() {
    const vector<pair<int, long long>> tests = {
        {1, 1},
        {2, 0},
        {5, 8},
        {100, 31395},
    };
    for (const auto& test : tests) {
        long long got = solve(test.first, false);
        if (got != test.second) {
            cerr << "Validation failed for N=" << test.first << ": got "
                 << got << ", expected " << test.second << "\n";
            std::exit(1);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    int n = 1'000'000;
    bool run_validation = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--validate") {
            run_validation = true;
        } else if (arg == "--no-validate") {
            run_validation = false;
        } else {
            n = stoi(arg);
        }
    }

    if (run_validation) {
        run_validations();
    }

    cout << solve(n, false) << "\n";
    return 0;
}
