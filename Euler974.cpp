#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>
#include <map>
#include <array>
#include <mutex>
#include <future>
#include <iomanip>

using namespace std;

// Represents the counts of digits 1, 3, 5, 7, 9
struct Counts {
    int c[5]; // Indices: 0->1, 1->3, 2->5, 3->7, 4->9

    bool operator<(const Counts& other) const {
        for (int i = 0; i < 5; ++i) {
            if (c[i] != other.c[i]) return c[i] < other.c[i];
        }
        return false;
    }

    bool operator==(const Counts& other) const {
        for (int i = 0; i < 5; ++i) {
            if (c[i] != other.c[i]) return false;
        }
        return true;
    }
};

const int DIGITS[] = {1, 3, 5, 7, 9};

// Global memoization cache
// Map from Counts -> array of 7 longs (count of permutations for each remainder)
map<Counts, array<long long, 7>> memo;
mutex memo_mutex;

// Modular powers of 10 for efficient remainder calculation
long long pow10mod7[50];

void init_pow() {
    pow10mod7[0] = 1;
    for (int i = 1; i < 50; ++i) {
        pow10mod7[i] = (pow10mod7[i - 1] * 10) % 7;
    }
}

// Factorials for multinomial check (optional, but good for debugging total counts)
long long factorial[30];
void init_fact() {
    factorial[0] = 1;
    for (int i = 1; i < 30; ++i) factorial[i] = factorial[i - 1] * i;
}

// DP function to count permutations
// Returns array[7] where index is remainder mod 7
array<long long, 7> get_rem_counts(Counts counts) {
    int total_len = 0;
    for (int x : counts.c) total_len += x;

    if (total_len == 0) {
        array<long long, 7> res = {0};
        res[0] = 1;
        return res;
    }

    {
        lock_guard<mutex> lock(memo_mutex);
        if (memo.count(counts)) {
            return memo[counts];
        }
    }

    array<long long, 7> total_res = {0};

    // Try placing each available digit at the CURRENT position (most significant of the remaining block? No, let's build from right to left or left to right)
    // Actually, order matters for mod 7. 
    // Is it simpler to place the *first* digit (most significant)?
    // Remainder = (d * 10^(L-1) + rest) % 7.
    // Yes.
    
    // We strictly need to be consistent. Let's assume we are placing the most significant digit.
    // Power needed is 10^(total_len - 1).
    int power = pow10mod7[total_len - 1];

    for (int i = 0; i < 5; ++i) {
        if (counts.c[i] > 0) {
            Counts next_counts = counts;
            next_counts.c[i]--;
            
            // Recurse
            // No need to lock for recursion, but cache access needs lock
            // To avoid deadlock or high contention, we might need thread-local caches or better locking strategy.
            // For now, simple lock is safest.
            
            // NOTE: To avoid holding lock during recursion, we release it. 
            // We already released it after checking cache.
            
            array<long long, 7> sub_res = get_rem_counts(next_counts);

            int digit = DIGITS[i];
            int current_rem_contribution = (digit * power) % 7;

            for (int r = 0; r < 7; ++r) {
                if (sub_res[r] > 0) {
                    int new_rem = (current_rem_contribution + r) % 7;
                    total_res[new_rem] += sub_res[r];
                }
            }
        }
    }

    {
        lock_guard<mutex> lock(memo_mutex);
        memo[counts] = total_res;
    }
    return total_res;
}

// Optimized DP without locking for single-threaded reconstruction
array<long long, 7> get_rem_counts_nolock(Counts counts, map<Counts, array<long long, 7>>& local_memo) {
    int total_len = 0;
    for (int x : counts.c) total_len += x;

    if (total_len == 0) {
        array<long long, 7> res = {0};
        res[0] = 1;
        return res;
    }

    if (local_memo.count(counts)) return local_memo[counts];

    array<long long, 7> total_res = {0};
    int power = pow10mod7[total_len - 1];

    for (int i = 0; i < 5; ++i) {
        if (counts.c[i] > 0) {
            Counts next_counts = counts;
            next_counts.c[i]--;
            
            array<long long, 7> sub_res = get_rem_counts_nolock(next_counts, local_memo);
            int digit = DIGITS[i];
            int current_rem_contribution = (digit * power) % 7;

            for (int r = 0; r < 7; ++r) {
                 total_res[(current_rem_contribution + r) % 7] += sub_res[r];
            }
        }
    }
    return local_memo[counts] = total_res;
}


// Worker function for a batch of count tuples
long long process_tuples(const vector<Counts>& tuples, int length) {
    long long count = 0;
    for (const auto& c : tuples) {
        // Condition: Sum of Digits % 3 == 0
        long long sum_digits = 0;
        for (int i = 0; i < 5; ++i) sum_digits += (long long)c.c[i] * DIGITS[i];
        if (sum_digits % 3 != 0) continue;

        // Condition 2: Divisible by 5 (Last digit is 5)
        // We need permutations of (Counts - one 5) ending in 5.
        // So we look at permutations of the remaining set.
        // Remaining set R has counts c.
        // But the full number is P * 10 + 5.
        // We want (P * 10 + 5) % 7 == 0 => 3*P + 5 == 0 mod 7 => 3*P == 2 => P == 3 mod 7.
        // So we need P (permutation of remaining digits) to have remainder 3.
        
        if (c.c[2] < 1) continue; // Must have at least one 5 to be the last digit
        
        Counts remaining = c;
        remaining.c[2]--; // Remove the last 5
        
        array<long long, 7> rems = get_rem_counts(remaining);
        count += rems[3]; // We need remainder 3
    }
    return count;
}

// Generate partitions
void generate_partitions(int index, int remaining_len, Counts current, vector<Counts>& result) {
    if (index == 5) {
        if (remaining_len == 0) {
            result.push_back(current);
        }
        return;
    }
    
    // Each count must be odd
    // Max count is remaining_len.
    // Try 1, 3, 5...
    for (int k = 1; k <= remaining_len; k += 2) {
        current.c[index] = k;
        // Optimization: remaining digits must sum to remaining_len - k
        // Minimum for future digits is 1 each.
        // So k cannot be such that remaining_len - k < (5 - 1 - index) * 1
        if (remaining_len - k >= (4 - index)) {
            generate_partitions(index + 1, remaining_len - k, current, result);
        }
    }
}

// Find logic for N-th number
// Returns string representation
string find_nth(long long N) {
    long long current_count = 0;
    
    // 1. Find Length
    for (int L = 5; ; L += 2) {
        // Generate tuples
        vector<Counts> all_tuples;
        generate_partitions(0, L, {0}, all_tuples);
        
        // Parallel processing to count
        int num_threads = thread::hardware_concurrency();
        vector<future<long long>> futures;
        int chunk_size = (all_tuples.size() + num_threads - 1) / num_threads;
        
        // Reset/Clear memo for new length to save memory? 
        // L=19 states ~3000. Not needed to clear.
        
        for (int i = 0; i < num_threads; ++i) {
            int start = i * chunk_size;
            if (start >= all_tuples.size()) break;
            int end = min((int)all_tuples.size(), start + chunk_size);
            
            // Copy sub-vector
            vector<Counts> chunk(all_tuples.begin() + start, all_tuples.begin() + end);
            futures.push_back(async(launch::async, process_tuples, chunk, L));
        }
        
        long long count_for_L = 0;
        for (auto& f : futures) count_for_L += f.get();
        
        if (current_count + count_for_L >= N) {
            // Target is in this length.
            // Find which Tuple
            
            // We need to iterate tuples in a deterministic order (lexicographical usually implied by problem?)
            // "Define Theta(n) be the nth very odd number".
            // Numbers are ordered by value.
            // Value order: Smallest length comes first.
            // Within same length, numbers are sorted naturally.
            // So we cannot just pick a tuple. We have to mix all tuples.
            // Actually, we construct digit by digit.
            // We know the length is L.
            // We know the set of valid tuples.
            
            map<Counts, array<long long, 7>> local_memo; // For fast reconstruction

            string res = "";
            
            // We have 'current_count' numbers before this Length.
            long long target_in_L = N - current_count;
            
            // We need to construct L digits.
            // Wait, the last digit is FIXED to 5.
            // So we really construct L-1 digits.
            // And we have a constraint that the counts of the FULL number (including the last 5) must match ONE of the valid tuples.
            // AND the digits used so far affect the remaining needed counts.
            
            // Let's refine the reconstruction state.
            // We are building a prefix.
            // At each step, we try placing a digit d (1, 3, 5, 7, 9).
            // Length remaining: rem_len.
            // Counts used so far.
            // We need to count how many ways to complete this prefix to form a VALID number (valid tuple + divisible by 7 + divisible by 3 + ends in 5).
            
            // State for reconstruction:
            // - Current index (0 to L-1).
            // - Current remainder mod 7.
            // - Current counts of digits used.
            // - Remainder condition for Div 3?
            //   Div 3 condition is on the SUM of digits. 
            //   Sum so far + Sum of future = 0 mod 3.
            
            // But actually, we only care if the FINAL counts match ANY valid tuple.
            // Valid tuples are those where counts are all odd and sum % 3 == 0.
            
            // So, for a trial digit d:
            // count_valid_completions = Sum over all Valid Tuples T of (ways to form T given we have used counts so far + d).
            
            Counts current_used = {{0}};
            int current_rem = 0;
            // The last digit is 5.
            // The first L-1 digits determine divisibility if we consider the last 5.
            
            // Let's iterate positions 0 to L-2 (since L-1 is fixed 5).
            for (int pos = 0; pos < L - 1; ++pos) {
                int power = pow10mod7[L - 1 - pos];
                
                for (int d_idx = 0; d_idx < 5; ++d_idx) {
                    int d = DIGITS[d_idx];
                    
                    // Try placing 'd'
                    // Calculate how many valid numbers start with 'res + d'
                    
                    long long ways = 0;
                    
                    // Temporarily update state
                    current_used.c[d_idx]++;
                    int next_rem = (current_rem + d * power) % 7;
                    
                    // Sum ways over valid target tuples
                    // A target tuple T is valid if T.c[i] >= current_used.c[i] for all i
                    // AND T is a valid tuple (odd counts, sum%3==0)
                    // AND the remaining counts form a suffix that satisfies the modulo 7 requirement.
                    // Modulo requirement:
                    // N = Prefix * 10^(suffix_len + 1) + Suffix * 10 + 5
                    // We maintain next_rem which is Prefix * 10^(suffix_len + 1)
                    // So we need: next_rem + 10 * Suffix + 5 == 0 (mod 7)
                    // 3 * Suffix == -(next_rem + 5) (mod 7)
                    // Suffix == -(next_rem + 5) * 5 (mod 7)
                    
                    long long term = (next_rem + 5) % 7;
                    int needed_suffix_rem = (int)(((-term * 5) % 7 + 7) % 7);
                    
                    // Check if SuffixCounts are non-negative
                    
                    for (const auto& T : all_tuples) {
                        // Check div 3
                        long long s = 0;
                        for(int k=0; k<5; ++k) s+=T.c[k]*DIGITS[k];
                        if (s%3 != 0) continue;
                        
                        // Check if T is reachable
                        bool possible = true;
                        Counts suffix_counts = {{0}};
                        for(int k=0; k<5; ++k) {
                            int needed = T.c[k] - current_used.c[k];
                            if (k == 2) needed--; // Reserved for last digit
                            if (needed < 0) { possible = false; break; }
                            suffix_counts.c[k] = needed;
                        }
                        if (!possible) continue;
                        
                        // Count permutations of suffix_counts matching needed_rem
                        array<long long, 7> c_rems = get_rem_counts_nolock(suffix_counts, local_memo);
                        ways += c_rems[needed_suffix_rem];
                    }
                    
                    if (target_in_L <= ways) {
                        // Accepted this digit
                        res += to_string(d);
                        current_rem = next_rem;
                        // current_used is already incremented
                        goto next_position; // Break from digit loop, continue to next pos
                    } else {
                        target_in_L -= ways;
                        current_used.c[d_idx]--; // Backtrack
                    }
                }
                // If we reach here, we failed to find a digit for this position
                cerr << "Error: Failed to reconstructed digit at pos " << pos << " L=" << L << " TargetRemaining=" << target_in_L << endl;
                return "ERROR";

                next_position:;
            }
            
            res += "5"; // Append last digit
            return res;
            
        } else {
            current_count += count_for_L;
        }
    }
}

int main() {
    init_pow();
    init_fact();
    
    // Validation
    cout << "Validating Theta(1)..." << endl;
    string t1 = find_nth(1);
    cout << "Theta(1) = " << t1 << " (Expected: 1117935)" << endl;
    if (t1 != "1117935") {
        cerr << "Validation Failed for Theta(1)!" << endl;
        // return 1;
    }

    cout << "Validating Theta(1000)..." << endl;
    string t1000 = find_nth(1000);
    cout << "Theta(1000) = " << t1000 << " (Expected: 11137955115)" << endl;
    if (t1000 != "11137955115") {
         cerr << "Validation Failed for Theta(1000)!" << endl;
         // return 1;
    }

    cout << "Calculating Theta(10^16)..." << endl;
    string ans = find_nth(10000000000000000LL);
    cout << "Theta(10^16) = " << ans << endl;

    return 0;
}
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using u64 = std::uint64_t;

namespace {

constexpr std::array<int, 5> kDigits{1, 3, 5, 7, 9};
constexpr int kAllOddMask = (1 << 5) - 1;

int digit_bit(int digit) {
    switch (digit) {
        case 1: return 1 << 0;
        case 3: return 1 << 1;
        case 5: return 1 << 2;
        case 7: return 1 << 3;
        case 9: return 1 << 4;
        default: return 0;
    }
}

u64 sat_add(u64 a, u64 b, u64 limit) {
    if (a >= limit) return limit;
    const u64 remaining = limit - a;
    return (b >= remaining) ? limit : (a + b);
}

struct DpTable {
    int length = 0;
    u64 limit = 0;
    std::vector<u64> dp;

    static constexpr int kR3 = 3;
    static constexpr int kR7 = 7;
    static constexpr int kMask = 32;
    static constexpr int kStates = kR3 * kR7 * kMask;  // 672

    std::size_t index(int pos, int r3, int r7, int mask) const {
        return ((static_cast<std::size_t>(pos) * kR3 + static_cast<std::size_t>(r3)) * kR7 +
                static_cast<std::size_t>(r7)) *
                   kMask +
               static_cast<std::size_t>(mask);
    }

    u64 at(int pos, int r3, int r7, int mask) const { return dp[index(pos, r3, r7, mask)]; }
};

DpTable build_dp(int length, u64 limit) {
    DpTable table;
    table.length = length;
    table.limit = limit;
    table.dp.assign(static_cast<std::size_t>(length + 1) * DpTable::kStates, 0);

    for (int r3 = 0; r3 < 3; ++r3) {
        for (int r7 = 0; r7 < 7; ++r7) {
            for (int mask = 0; mask < 32; ++mask) {
                table.dp[table.index(length, r3, r7, mask)] =
                    (r3 == 0 && r7 == 0 && mask == kAllOddMask) ? 1ULL : 0ULL;
            }
        }
    }

    for (int pos = length - 1; pos >= 0; --pos) {
        const bool is_last = (pos == length - 1);
        for (int r3 = 0; r3 < 3; ++r3) {
            for (int r7 = 0; r7 < 7; ++r7) {
                for (int mask = 0; mask < 32; ++mask) {
                    u64 sum = 0;
                    if (is_last) {
                        const int d = 5;
                        const int nr3 = (r3 + d) % 3;
                        const int nr7 = (r7 * 10 + d) % 7;
                        const int nmask = mask ^ digit_bit(d);
                        sum = table.at(pos + 1, nr3, nr7, nmask);
                    } else {
                        for (int d : kDigits) {
                            const int nr3 = (r3 + d) % 3;
                            const int nr7 = (r7 * 10 + d) % 7;
                            const int nmask = mask ^ digit_bit(d);
                            sum = sat_add(sum, table.at(pos + 1, nr3, nr7, nmask), limit);
                            if (sum == limit) break;
                        }
                    }
                    table.dp[table.index(pos, r3, r7, mask)] = sum;
                }
            }
        }
    }

    return table;
}

u64 count_len(int length, u64 limit) {
    const DpTable dp = build_dp(length, limit);
    return dp.at(0, 0, 0, 0);
}

std::vector<u64> count_lengths_parallel(const std::vector<int>& lengths, u64 limit) {
    std::vector<u64> counts(lengths.size(), 0);
    unsigned num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    num_threads = std::min<unsigned>(num_threads, static_cast<unsigned>(lengths.size()));

    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    for (unsigned t = 0; t < num_threads; ++t) {
        workers.emplace_back([&, t]() {
            for (std::size_t i = t; i < lengths.size(); i += num_threads) {
                counts[i] = count_len(lengths[i], limit);
            }
        });
    }
    for (auto& th : workers) th.join();
    return counts;
}

std::string unrank_in_length(const DpTable& dp, u64 k) {
    const int L = dp.length;
    int r3 = 0, r7 = 0, mask = 0;
    std::string out;
    out.reserve(static_cast<std::size_t>(L));

    for (int pos = 0; pos < L; ++pos) {
        const bool is_last = (pos == L - 1);
        if (is_last) {
            const int d = 5;
            out.push_back(static_cast<char>('0' + d));
            r3 = (r3 + d) % 3;
            r7 = (r7 * 10 + d) % 7;
            mask ^= digit_bit(d);
            continue;
        }

        bool chosen = false;
        for (int d : kDigits) {
            const int nr3 = (r3 + d) % 3;
            const int nr7 = (r7 * 10 + d) % 7;
            const int nmask = mask ^ digit_bit(d);
            const u64 cnt = dp.at(pos + 1, nr3, nr7, nmask);
            if (k > cnt) {
                k -= cnt;
                continue;
            }

            out.push_back(static_cast<char>('0' + d));
            r3 = nr3;
            r7 = nr7;
            mask = nmask;
            chosen = true;
            break;
        }
        if (!chosen) throw std::runtime_error("Unranking failed: k out of range");
    }

    if (!(r3 == 0 && r7 == 0 && mask == kAllOddMask)) {
        throw std::runtime_error("Unranking produced invalid terminal state");
    }
    return out;
}

std::string theta(u64 n) {
    if (n == 0) throw std::invalid_argument("n must be >= 1");

    u64 cumulative = 0;
    int next_length = 5;
    constexpr int kBatch = 16;  // odd lengths: 32 digits per batch

    while (true) {
        std::vector<int> lengths;
        lengths.reserve(kBatch);
        for (int i = 0; i < kBatch; ++i) {
            lengths.push_back(next_length);
            next_length += 2;
        }

        const std::vector<u64> counts = count_lengths_parallel(lengths, n);
        for (std::size_t i = 0; i < lengths.size(); ++i) {
            const u64 cnt = counts[i];
            if (sat_add(cumulative, cnt, n) >= n) {
                const int L = lengths[i];
                const u64 k = n - cumulative;
                const DpTable dp = build_dp(L, n);
                return unrank_in_length(dp, k);
            }
            cumulative = sat_add(cumulative, cnt, n);
        }
    }
}

bool is_very_odd(const std::string& s) {
    std::array<int, 10> counts{};
    int mod105 = 0;
    for (char c : s) {
        if (c < '0' || c > '9') return false;
        const int d = c - '0';
        if (d % 2 == 0) return false;
        mod105 = (mod105 * 10 + d) % 105;
        counts[d]++;
    }
    if (s.empty() || s.back() != '5') return false;
    if (mod105 != 0) return false;
    for (int d : kDigits) {
        if (counts[d] % 2 == 0) return false;
    }
    return true;
}

void require_equal(const std::string& got, const std::string& expected, int exit_code) {
    if (got != expected) {
        std::cerr << "Validation failed: expected " << expected << " but got " << got << "\n";
        std::exit(exit_code);
    }
}

}  // namespace

int main() {
    {
        const std::vector<std::string> expected_first10 = {
            "1117935", "1131795", "1197315", "1311975", "1371195",
            "1711395", "1713915", "1739115", "1773975", "1911735",
        };
        for (std::size_t i = 0; i < expected_first10.size(); ++i) {
            const std::string got = theta(static_cast<u64>(i + 1));
            require_equal(got, expected_first10[i], 10 + static_cast<int>(i));
        }
    }

    require_equal(theta(1), "1117935", 1);
    require_equal(theta(1000), "11137955115", 2);

    const u64 target = 10'000'000'000'000'000ULL;
    const auto t0 = std::chrono::steady_clock::now();
    const std::string ans = theta(target);
    const auto t1 = std::chrono::steady_clock::now();

    if (!is_very_odd(ans)) {
        std::cerr << "Validation failed: Theta(10^16) does not satisfy constraints\n";
        return 3;
    }

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << ans << "\n";
    std::cerr << "[OK] Theta(10^16) computed in " << ms << " ms\n";
    return 0;
}
