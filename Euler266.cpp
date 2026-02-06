#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

struct BigInt {
  static const int WORDS = 6;
  uint64_t data[WORDS];

  BigInt() { std::memset(data, 0, sizeof(data)); }

  BigInt(uint64_t val) {
    std::memset(data, 0, sizeof(data));
    data[0] = val;
  }

  bool operator<(const BigInt &other) const {
    for (int i = WORDS - 1; i >= 0; --i) {
      if (data[i] != other.data[i])
        return data[i] < other.data[i];
    }
    return false;
  }

  bool operator<=(const BigInt &other) const { return !(other < *this); }

  bool operator>(const BigInt &other) const { return other < *this; }

  bool operator==(const BigInt &other) const {
    for (int i = 0; i < WORDS; ++i)
      if (data[i] != other.data[i])
        return false;
    return true;
  }

  void multiplyBy(uint64_t p) {
    unsigned __int128 carry = 0;
    for (int i = 0; i < WORDS; ++i) {
      unsigned __int128 prod = (unsigned __int128)data[i] * p + carry;
      data[i] = (uint64_t)prod;
      carry = prod >> 64;
    }
  }

  BigInt operator*(const BigInt &other) const {
    BigInt res;
    for (int i = 0; i < WORDS; ++i) {
      unsigned __int128 carry = 0;
      for (int j = 0; j < WORDS - i; ++j) {
        unsigned __int128 prod = (unsigned __int128)data[i] * other.data[j] +
                                 res.data[i + j] + carry;
        res.data[i + j] = (uint64_t)prod;
        carry = prod >> 64;
      }
    }
    return res;
  }

  std::string toString() const {
    BigInt temp = *this;
    std::string s = "";
    if (temp.isZero())
      return "0";

    while (!temp.isZero()) {
      uint64_t rem = temp.divMod(10);
      s += std::to_string(rem);
    }
    std::reverse(s.begin(), s.end());
    return s;
  }

  uint64_t divMod(uint64_t divisor) {
    unsigned __int128 rem = 0;
    for (int i = WORDS - 1; i >= 0; --i) {
      unsigned __int128 cur = data[i] + (rem << 64);
      data[i] = (uint64_t)(cur / divisor);
      rem = cur % divisor;
    }
    return (uint64_t)rem;
  }

  bool isZero() const {
    for (int i = 0; i < WORDS; ++i)
      if (data[i] != 0)
        return false;
    return true;
  }

  double toDouble() const {
    double res = 0;
    for (int i = WORDS - 1; i >= 0; --i) {
      res = res * 18446744073709551616.0 + data[i];
    }
    return res;
  }

  uint64_t mod10_16() const {
    BigInt temp = *this;
    return temp.divMod(10000000000000000ULL);
  }
};

struct Entry {
  BigInt val;
  double logVal;

  bool operator<(const Entry &other) const { return logVal < other.logVal; }
};

std::vector<uint64_t> primes;
void generatePrimes(int limit) {
  std::vector<bool> is_prime(limit, true);
  is_prime[0] = is_prime[1] = false;
  for (int p = 2; p * p < limit; p++) {
    if (is_prime[p]) {
      for (int i = p * p; i < limit; i += p)
        is_prime[i] = false;
    }
  }
  for (int p = 2; p < limit; p++) {
    if (is_prime[p])
      primes.push_back(p);
  }
}

BigInt P;

void getSubsets(const std::vector<uint64_t> &p_list,
                std::vector<Entry> &entries) {
  entries.push_back({BigInt(1), 0.0});

  for (uint64_t p : p_list) {
    int n = entries.size();
    double logP = std::log((double)p);
    for (int i = 0; i < n; ++i) {
      Entry newEntry = entries[i];
      newEntry.val.multiplyBy(p);
      newEntry.logVal += logP;
      entries.push_back(newEntry);
    }
  }
}

bool squareLeq(const BigInt &c, const BigInt &p) {
  BigInt sq = c * c;
  return sq <= p;
}

int main() {
  auto start_time = std::chrono::high_resolution_clock::now();

  generatePrimes(190);

  P = BigInt(1);
  for (uint64_t p : primes) {
    P.multiplyBy(p);
  }

  std::cout << "Primes count: " << primes.size() << std::endl;

  std::vector<uint64_t> set1, set2;
  for (size_t i = 0; i < primes.size(); ++i) {
    if (i % 2 == 0)
      set1.push_back(primes[i]);
    else
      set2.push_back(primes[i]);
  }

  std::vector<Entry> left, right;
  getSubsets(set1, left);
  getSubsets(set2, right);

  std::sort(left.begin(), left.end());
  std::sort(right.begin(), right.end());

  std::cout << "Left size: " << left.size() << ", Right size: " << right.size()
            << std::endl;

  double totalLogP = 0;
  for (auto p : primes)
    totalLogP += std::log((double)p);
  double targetLog = totalLogP / 2.0;

  std::cout << "Target Log: " << targetLog << std::endl;

  BigInt bestVal(0);
  double bestLog = -1.0;

  for (const auto &l_entry : left) {
    double rem = targetLog - l_entry.logVal;
    auto it =
        std::lower_bound(right.begin(), right.end(), Entry{BigInt(0), rem});

    int idx = std::distance(right.begin(), it);
    int start = std::min((int)right.size() - 1, idx + 50);
    int end = std::max(0, idx - 50);

    for (int k = start; k >= end; --k) {
      if (l_entry.logVal + right[k].logVal < bestLog - 1e-9) {
        break;
      }

      BigInt prod = l_entry.val * right[k].val;

      if (squareLeq(prod, P)) {
        if (bestLog < 0 || prod > bestVal) {
          bestVal = prod;
          bestLog = l_entry.logVal + right[k].logVal;
        }
        break;
      }
    }
  }

  uint64_t last16 = bestVal.mod10_16();
  std::cout << "PSR mod 10^16: " << last16 << std::endl;

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << "Time: " << elapsed.count() << "s" << std::endl;

  return 0;
}
