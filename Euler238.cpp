#include <atomic>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

// Type aliases for easier reading
using u64 = unsigned long long;

// BBS Parameters
const u64 MODULUS = 20300713;
const u64 SEED = 14025256;
const u64 TARGET_K_SMALL = 1000;
const u64 EXPECTED_SMALL_SUM = 4742;

// Function to generate the sequence of digits and find the period
struct BBSSequence {
  std::string w_str;
  std::vector<int> w_digits;
  std::vector<u64> prefix_sum;
  u64 total_sum_per_period;

  void generate() {
    u64 s = SEED;
    std::unordered_map<u64, int> seen;
    int index = 0;

    // We know it will cycle eventually. Since modulus is ~2e7, cycle is shorter
    // than that. Actually, we need to generate until we hit a repeated state
    // s_n. The problem statement says s_0 = 14025256, s_{n+1} = s_n^2 mod
    // 20300713. Digits are concatenated.

    while (seen.find(s) == seen.end()) {
      seen[s] = index;

      // Append digits of s to w
      std::string s_str = std::to_string(s);
      w_str += s_str;
      for (char c : s_str) {
        w_digits.push_back(c - '0');
      }

      s = (s * s) % MODULUS;
      index++;
    }

    // Note: The problem implies w starts from s0.
    // "Concatenate these numbers s0s1s2... to create a string w"
    // Cycle detected at state 's'.
    // However, for the purpose of the sum of digits, we treat w as the sequence
    // of digits. The "infinite length" string is formed by repeating the cycle.
    // Wait, does the string w start repeating exactly?
    // s0 -> s1 -> ... -> sk -> ... -> sk+L -> ... where sk = sk+L
    // The digits from s0 to sk-1 are the pre-period.
    // But usually BBS cycles are pure if s0 is a quadratic residue?
    // Or if we enter a cycle, the digits corresponding to the cycle repeat.
    // Let's verify if s0 is part of the cycle or pre-period.

    std::cout << "Cycle detected at index " << index << " (value: " << s << ")"
              << std::endl;
    std::cout << "First seen at index " << seen[s] << std::endl;

    // For this specific problem, let's assume valid range logic later.
    // Just construct the full digit sequence for one period if it's pure,
    // or enough to cover the non-periodic part + one period.

    // Prefix sums
    prefix_sum.push_back(0);
    u64 current_sum = 0;
    for (int d : w_digits) {
      current_sum += d;
      prefix_sum.push_back(current_sum);
    }
    total_sum_per_period =
        current_sum; // This might be slightly wrong if pre-period exists.
    // If pre-period exists, the "infinite string" isn't just repeating
    // 'total_sum_per_period' from the start. But let's check if pre-period is
    // empty (0).
    if (seen[s] != 0) {
      std::cerr
          << "Warning: Pre-period detected. Simplification assumption failed."
          << std::endl;
    }
  }
};

// Full solver
void solve_full() {
  BBSSequence bbs;
  std::cout << "Generating sequence..." << std::endl;
  bbs.generate();

  u64 S = bbs.total_sum_per_period;
  size_t L = bbs.w_digits.size();

  std::cout << "S (Sum): " << S << std::endl;
  std::cout << "L (Length): " << L << std::endl;

  // Build set R of residues present in prefix sums
  // std::vector<bool> is space efficient (1 bit per bool)
  std::cout << "Building residue table..." << std::endl;
  std::vector<bool> R(S + 1, false);
  std::vector<u64> P_mod_S(L); // Cache P[i] % S for fast access

  // P[0] = 0
  R[0] = true;
  // Fill R and P_mod_S
  // prefix_sum has L+1 entries: 0, w1, w1+w2...
  // prefix_sum[i] corresponds to sum of first i digits.
  // P_mod_S[i] should store prefix_sum[i] % S.
  // NOTE: In the search logic, we use P[z-1].
  // z=1 -> P[0]. z=2 -> P[1].
  // So we need P_mod_S to store the sequence of prefix sums modulo S.

  // Reconstruct just the mods to save memory if needed (though bbs has them)
  // We already have bbs.prefix_sum.
  for (size_t i = 0; i < L; ++i) {
    u64 val = bbs.prefix_sum[i] % S;
    P_mod_S[i] = val;
    // The residue is present.
    // Wait, R represents the set of prefix sums that exist.
    // We iterate through period.
    // prefix_sum[0]...prefix_sum[L-1].
    // prefix_sum[L] = S == 0 mod S.
    // So we just mark all prefix_sum[i] % S.
  }

  for (size_t i = 0; i <= L; ++i) { // Include L to ensure 0/S is marked
    R[bbs.prefix_sum[i] % S] = true;
  }

  std::cout << "Computing p(k) for all k..." << std::endl;

  // We need to sum p(k) for k = 1 to S.
  // Store them to handle the full range multiplication.

  std::vector<u64> p_vals(S + 1);
  std::atomic<u64> processed_count(0);

  unsigned int num_threads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads;
  u64 chunk_size = S / num_threads;

  auto worker = [&](u64 start_k, u64 end_k) {
    for (u64 k = start_k; k <= end_k; ++k) {
      // Find minimal z such that (P[z-1] + k) % S is in R.
      u64 z = 1;
      while (true) {
        // accessing P_mod_S[(z-1) % L]
        // and comparing with k
        u64 start_residue = P_mod_S[(z - 1) % L];
        u64 target = (start_residue + k);
        if (target >= S)
          target -= S; // manually mod S

        if (R[target]) {
          p_vals[k] = z;
          break;
        }
        z++;
      }
    }
    processed_count += (end_k - start_k + 1);
  };

  for (unsigned int i = 0; i < num_threads; ++i) {
    u64 start = 1 + i * chunk_size;
    u64 end = (i == num_threads - 1) ? S : (start + chunk_size - 1);
    if (start <= end) {
      threads.emplace_back(worker, start, end);
    }
  }

  // Monitor progress
  // Move to join
  for (auto &t : threads) {
    if (t.joinable())
      t.join();
  }

  std::cout << "All p(k) computed." << std::endl;

  // Calculate total answer
  // Limit: 2 * 10^15
  u64 LIMIT = 2000000000000000ULL;
  u64 num_periods = LIMIT / S;
  u64 remainder = LIMIT % S;

  u64 sum_p_k_one_period = 0;
  u64 sum_p_k_remainder = 0;

  // Parallel reduction for sum? Or just linear is fast enough (10^8 adds ~
  // 0.1s)
  for (u64 k = 1; k <= S; ++k) {
    sum_p_k_one_period += p_vals[k];
  }

  for (u64 k = 1; k <= remainder; ++k) {
    sum_p_k_remainder += p_vals[k];
  }

  u64 total_ans = num_periods * sum_p_k_one_period + sum_p_k_remainder;

  std::cout << "Sum for one period: " << sum_p_k_one_period << std::endl;
  std::cout << "Remainder sum (" << remainder << "): " << sum_p_k_remainder
            << std::endl;
  std::cout << "Total Answer: " << total_ans << std::endl;
}

int main() {
  // Run validation if argument provided or just default
  // For now, let's run full solver which includes generation
  // Maybe verify small case again inside logic?
  // The previous main verified small case.
  // Let's keep a quick check.

  // Re-instantiate for small check? No, waste of time.
  // Just run the fast solver.
  solve_full();
  return 0;
}
