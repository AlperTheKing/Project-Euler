#include <algorithm>
#include <cmath>
#include <future>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

using int64 = long long;

const int64 MOD = 100000;
const int64 M1 = 32;
const int64 M2 = 3125;
const int64 PHI_M2 = 2500;

int64 P_table[M2 + 1];

int64 power(int64 base, int64 exp, int64 mod) {
  int64 res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1)
      res = (res * base) % mod;
    base = (base * base) % mod;
    exp /= 2;
  }
  return res;
}

int64 modInverse(int64 n, int64 mod) {
  return power(n, mod > 2000 ? PHI_M2 - 1 : mod / 2, mod);
}

void init_P_table() {
  P_table[0] = 1;
  for (int64 i = 1; i <= M2; ++i) {
    if (i % 5 != 0) {
      P_table[i] = (P_table[i - 1] * i) % M2;
    } else {
      P_table[i] = P_table[i - 1];
    }
  }
}

int64 get_P_mod_3125(int64 n) {
  int64 q = n / M2;
  int64 r = n % M2;
  int64 res = (q % 2 == 0) ? 1 : (M2 - 1);
  res = (res * P_table[r]) % M2;
  return res;
}

int64 factorial_no_5_mod_3125(int64 n) {
  if (n == 0)
    return 1;
  int64 part1 = get_P_mod_3125(n);
  int64 part2 = factorial_no_5_mod_3125(n / 5);
  return (part1 * part2) % M2;
}

int64 count_factors_5(int64 n) {
  int64 count = 0;
  while (n > 0) {
    count += n / 5;
    n /= 5;
  }
  return count;
}

int64 solve_crt(int64 res3125) {
  int64 inv21 = 29;
  int64 target = (32 - (res3125 % 32)) % 32;
  int64 k = (target * inv21) % 32;
  return res3125 + k * M2;
}

int64 f(int64 n) {
  int64 fact_no_5 = factorial_no_5_mod_3125(n);
  int64 num_5s = count_factors_5(n);

  int64 two_exp = num_5s % PHI_M2;
  int64 inv_2 = power(2, PHI_M2 - two_exp, M2);

  int64 res3125 = (fact_no_5 * inv_2) % M2;

  if (n >= 10) {
    return solve_crt(res3125);
  } else {
    int64 res = 1;
    for (int i = 1; i <= n; ++i)
      res *= i;
    while (res % 10 == 0)
      res /= 10;
    return res % 100000;
  }
}

int main() {
  init_P_table();

  auto handle9 = std::async(std::launch::async, f, 9);
  auto handle10 = std::async(std::launch::async, f, 10);
  auto handle20 = std::async(std::launch::async, f, 20);

  int64 r9 = handle9.get();
  int64 r10 = handle10.get();
  int64 r20 = handle20.get();

  std::cout << "Validation:" << std::endl;
  std::cout << "f(9)  = " << r9 << " (Expected: 36288)" << std::endl;
  std::cout << "f(10) = " << r10 << " (Expected: 36288)" << std::endl;
  std::cout << "f(20) = " << r20 << " (Expected: 17664)" << std::endl;

  bool ok = (r9 == 36288 && r10 == 36288 && r20 == 17664);
  if (!ok) {
    std::cerr << "Validation Failed!" << std::endl;
    return 1;
  }

  std::cout << "Validation Passed. Computing f(10^12)..." << std::endl;

  int64 N = 1000000000000LL;
  int64 result = f(N);

  std::cout << "f(10^12) = " << std::setfill('0') << std::setw(5) << result
            << std::endl;

  return 0;
}
