#include <algorithm>
#include <iostream>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

const int NUM_PRODUCTS = 5;
const int A[] = {5248, 1312, 2624, 5760, 3936};
const int B[] = {640, 1888, 3776, 3776, 5664};

const int GROUP_ID[] = {0, 1, 1, 2, 1};

struct Fraction {
  long long u, v;
  bool operator<(const Fraction &other) const {
    if (u != other.u)
      return u < other.u;
    return v < other.v;
  }
  bool operator==(const Fraction &other) const {
    return u == other.u && v == other.v;
  }
};

long long gcd(long long a, long long b) { return b == 0 ? a : gcd(b, a % b); }

Fraction reduce(long long u, long long v) {
  long long g = gcd(u, v);
  return {u / g, v / g};
}

long long SA = 0, SB = 0;

std::mutex result_mutex;
std::set<Fraction> solutions;

void check_m(Fraction m) {
  if (m.u <= m.v)
    return;

  long long u = m.u;
  long long v = m.v;
  long long u2 = u * u;
  long long v2 = v * v;

  long long nk[] = {5, 59, 59};
  long long dk[] = {41, 41, 90};

  long long C[3];
  long long min_S[3] = {1, 3, 1};
  long long max_S[3] = {0, 0, 0};

  for (int g = 0; g < 3; ++g) {
    long long num = u * nk[g];
    long long den = v * dk[g];
    Fraction ratio = reduce(num, den);
    long long N_gm = ratio.u;
    long long D_gm = ratio.v;

    long long current_group_max_S = 0;
    for (int i = 0; i < 5; ++i) {
      if (GROUP_ID[i] == g) {
        if (D_gm > A[i]) {
          current_group_max_S = -1;
          break;
        }
        if (N_gm > B[i]) {
          current_group_max_S = -1;
          break;
        }
        long long max_s_i = std::min(A[i] / D_gm, B[i] / N_gm);
        if (max_s_i < 1) {
          current_group_max_S = -1;
          break;
        }
        current_group_max_S += max_s_i;
      }
    }

    if (current_group_max_S == -1)
      return;
    max_S[g] = current_group_max_S;

    long long LCM_DK = 3690;
    long long term1 = v2 * SB * LCM_DK;
    long long term2 = u2 * SA * nk[g] * (LCM_DK / dk[g]);

    C[g] = D_gm * (term1 - term2);
  }

  for (long long s0 = min_S[0]; s0 <= max_S[0]; ++s0) {
    for (long long s2 = min_S[2]; s2 <= max_S[2]; ++s2) {
      long long partial = C[0] * s0 + C[2] * s2;
      if (C[1] == 0) {
        if (partial == 0) {
          {
            std::lock_guard<std::mutex> lock(result_mutex);
            solutions.insert(m);
            return;
          }
        }
      } else {
        if (-partial % C[1] == 0) {
          long long s1 = -partial / C[1];
          if (s1 >= min_S[1] && s1 <= max_S[1]) {
            std::lock_guard<std::mutex> lock(result_mutex);
            solutions.insert(m);
            return;
          }
        }
      }
    }
  }
}

void solve_chunk(const std::vector<Fraction> &candidates, int start, int end) {
  for (int i = start; i < end; ++i) {
    check_m(candidates[i]);
  }
}

int main() {
  for (int i : A)
    SA += i;
  for (int i : B)
    SB += i;

  std::set<Fraction> candidates;
  int a_max = 5248;
  int b_max = 640;

  for (int a = 1; a <= a_max; ++a) {
    for (int b = 1; b <= b_max; ++b) {
      long long num = (long long)b * 41;
      long long den = (long long)a * 5;
      Fraction m = reduce(num, den);
      if (m.u > m.v) {
        candidates.insert(m);
      }
    }
  }

  std::vector<Fraction> cand_vec(candidates.begin(), candidates.end());
  std::cout << "Generated " << cand_vec.size() << " candidates." << std::endl;

  int num_threads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads;
  int chunk = cand_vec.size() / num_threads;

  for (int i = 0; i < num_threads; ++i) {
    int start = i * chunk;
    int end = (i == num_threads - 1) ? cand_vec.size() : start + chunk;
    threads.emplace_back(solve_chunk, std::cref(cand_vec), start, end);
  }

  for (auto &t : threads)
    t.join();

  std::cout << "Found " << solutions.size() << " solutions." << std::endl;

  if (!solutions.empty()) {
    std::vector<Fraction> sorted_sols(solutions.begin(), solutions.end());
    std::sort(sorted_sols.begin(), sorted_sols.end(),
              [](const Fraction &a, const Fraction &b) {
                return (double)a.u / a.v < (double)b.u / b.v;
              });

    Fraction best = sorted_sols.back();
    std::cout << "Largest m: " << best.u << "/" << best.v << std::endl;
  }

  return 0;
}
