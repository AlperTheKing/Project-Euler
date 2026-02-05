#include <cmath>
#include <future>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

const long long MOD_MAIN = 1475789056;

using int128 = __int128;
using ll = long long;

ll power(ll base, ll exp, ll mod) {
  ll res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1)
      res = (int128)res * base % mod;
    base = (int128)base * base % mod;
    exp /= 2;
  }
  return res;
}

ll get_phi(ll n) {
  ll result = n;
  for (ll i = 2; i * i <= n; i++) {
    if (n % i == 0) {
      while (n % i == 0)
        n /= i;
      result -= result / i;
    }
  }
  if (n > 1)
    result -= result / n;
  return result;
}

ll hyper_tower_2(ll height, ll m) {
  if (m == 1)
    return 0;
  if (height == 0)
    return 1 % m;
  if (height == 1)
    return 2 % m;
  if (height == 2)
    return 4 % m;
  if (height == 3)
    return 16 % m;
  if (height == 4)
    return 65536 % m;

  ll phi = get_phi(m);
  ll exponent_val = hyper_tower_2(height - 1, phi);

  ll eff_exponent = exponent_val + phi;
  return power(2, eff_exponent, m);
}

ll solve_n0() { return 1; }
ll solve_n1() { return 3; }
ll solve_n2() { return 7; }
ll solve_n3() { return 61; }

ll solve_n4() {
  ll term = hyper_tower_2(7, MOD_MAIN);
  return (term - 3 + MOD_MAIN) % MOD_MAIN;
}

ll solve_large_n() {
  ll term = hyper_tower_2(200, MOD_MAIN);
  return (term - 3 + MOD_MAIN) % MOD_MAIN;
}

int main() {
  auto f0 = std::async(std::launch::async, solve_n0);
  auto f1 = std::async(std::launch::async, solve_n1);
  auto f2 = std::async(std::launch::async, solve_n2);
  auto f3 = std::async(std::launch::async, solve_n3);
  auto f4 = std::async(std::launch::async, solve_n4);
  auto f5 = std::async(std::launch::async, solve_large_n);
  auto f6 = std::async(std::launch::async, solve_large_n);

  ll s0 = f0.get();
  ll s1 = f1.get();
  ll s2 = f2.get();
  ll s3 = f3.get();
  ll s4 = f4.get();
  ll s5 = f5.get();
  ll s6 = f6.get();

  ll total = (s0 + s1 + s2 + s3 + s4 + s5 + s6) % MOD_MAIN;

  std::cout << total << std::endl;

  return 0;
}
