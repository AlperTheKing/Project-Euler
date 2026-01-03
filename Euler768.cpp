#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

namespace {

struct BigInt {
    static constexpr uint32_t kBase = 1000000000U;
    vector<uint32_t> digits;

    BigInt(uint64_t value = 0) {
        if (value == 0) {
            digits.push_back(0);
            return;
        }
        while (value > 0) {
            digits.push_back(static_cast<uint32_t>(value % kBase));
            value /= kBase;
        }
    }

    bool is_zero() const {
        return digits.size() == 1 && digits[0] == 0;
    }

    void trim() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
    }

    void add(const BigInt& other) {
        uint64_t carry = 0;
        size_t n = max(digits.size(), other.digits.size());
        if (digits.size() < n) digits.resize(n, 0);
        for (size_t i = 0; i < n; ++i) {
            uint64_t sum = carry + digits[i];
            if (i < other.digits.size()) sum += other.digits[i];
            digits[i] = static_cast<uint32_t>(sum % kBase);
            carry = sum / kBase;
        }
        if (carry) digits.push_back(static_cast<uint32_t>(carry));
    }

    BigInt operator*(const BigInt& other) const {
        if (is_zero() || other.is_zero()) return BigInt(0);
        BigInt res;
        res.digits.assign(digits.size() + other.digits.size(), 0);
        for (size_t i = 0; i < digits.size(); ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < other.digits.size() || carry; ++j) {
                uint64_t cur = res.digits[i + j] + carry;
                if (j < other.digits.size()) {
                    cur += static_cast<uint64_t>(digits[i]) * other.digits[j];
                }
                res.digits[i + j] = static_cast<uint32_t>(cur % kBase);
                carry = cur / kBase;
            }
        }
        res.trim();
        return res;
    }

    bool equals_uint64(uint64_t value) const {
        BigInt other(value);
        return digits == other.digits;
    }
};

ostream& operator<<(ostream& os, const BigInt& value) {
    if (value.digits.empty()) return os << 0;
    os << value.digits.back();
    for (int i = static_cast<int>(value.digits.size()) - 2; i >= 0; --i) {
        os << setw(9) << setfill('0') << value.digits[i];
    }
    return os;
}

struct Poly {
    vector<BigInt> coeffs;

    explicit Poly(int max_deg = 0) : coeffs(max_deg + 1, BigInt(0)) {}

    static Poly one() {
        Poly p(0);
        p.coeffs[0] = BigInt(1);
        return p;
    }

    Poly multiply(const Poly& other, int max_deg) const {
        int deg_a = static_cast<int>(coeffs.size()) - 1;
        int deg_b = static_cast<int>(other.coeffs.size()) - 1;
        int new_deg = min(max_deg, deg_a + deg_b);
        Poly res(new_deg);
        for (int i = 0; i <= deg_a; ++i) {
            if (coeffs[i].is_zero()) continue;
            for (int j = 0; j <= deg_b && i + j <= new_deg; ++j) {
                if (other.coeffs[j].is_zero()) continue;
                BigInt prod = coeffs[i] * other.coeffs[j];
                res.coeffs[i + j].add(prod);
            }
        }
        return res;
    }
};

vector<long long> trim_poly(vector<long long> poly) {
    while (poly.size() > 1 && poly.back() == 0) {
        poly.pop_back();
    }
    return poly;
}

vector<long long> divide_poly(vector<long long> a, const vector<long long>& b) {
    int da = static_cast<int>(a.size()) - 1;
    int db = static_cast<int>(b.size()) - 1;
    vector<long long> q(max(0, da - db + 1), 0);
    for (int i = da; i >= db; --i) {
        long long coef = a[i];
        if (coef == 0) continue;
        q[i - db] = coef;
        for (int j = 0; j <= db; ++j) {
            a[i - db + j] -= coef * b[j];
        }
    }
    return trim_poly(q);
}

vector<long long> cyclotomic(int n) {
    vector<vector<long long>> phi(n + 1);
    phi[1] = {-1, 1};
    for (int k = 2; k <= n; ++k) {
        vector<long long> poly(k + 1, 0);
        poly[0] = -1;
        poly[k] = 1;
        for (int d = 1; d < k; ++d) {
            if (k % d == 0) {
                poly = divide_poly(poly, phi[d]);
            }
        }
        phi[k] = trim_poly(poly);
    }
    return phi[n];
}

vector<int> multiply_by_x(const vector<int>& v, const vector<int>& rel) {
    int d = static_cast<int>(v.size());
    vector<int> res(d, 0);
    for (int i = 0; i + 1 < d; ++i) {
        res[i + 1] += v[i];
    }
    int h = v[d - 1];
    if (h != 0) {
        for (int i = 0; i < d; ++i) {
            res[i] -= h * rel[i];
        }
    }
    return res;
}

struct VectorHash {
    size_t operator()(const vector<int>& v) const noexcept {
        uint64_t h = 0xcbf29ce484222325ULL;
        for (int x : v) {
            uint64_t y = static_cast<uint64_t>(static_cast<int64_t>(x) ^ (static_cast<int64_t>(x) >> 32));
            h ^= y + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return static_cast<size_t>(h);
    }
};

vector<uint64_t> count_zero_sums(int S, int m) {
    vector<long long> phi = cyclotomic(S);
    int d = static_cast<int>(phi.size()) - 1;
    vector<int> rel(d, 0);
    for (int i = 0; i < d; ++i) rel[i] = static_cast<int>(phi[i]);

    vector<vector<int>> powers(S, vector<int>(d, 0));
    powers[0][0] = 1;
    for (int i = 1; i < S; ++i) {
        powers[i] = multiply_by_x(powers[i - 1], rel);
    }

    int mid = S / 2;
    int right_n = S - mid;
    uint64_t left_masks = 1ULL << mid;
    uint64_t right_masks = 1ULL << right_n;

    unordered_map<vector<int>, vector<uint32_t>, VectorHash> left_sums;
    left_sums.reserve(left_masks * 2);

    vector<int> sum(d, 0);
    for (uint64_t mask = 0; mask < left_masks; ++mask) {
        fill(sum.begin(), sum.end(), 0);
        int cnt = 0;
        for (int k = 0; k < mid; ++k) {
            if (mask & (1ULL << k)) {
                ++cnt;
                const auto& pv = powers[k];
                for (int i = 0; i < d; ++i) sum[i] += pv[i];
            }
        }
        if (cnt > m) continue;
        auto it = left_sums.find(sum);
        if (it == left_sums.end()) {
            it = left_sums.emplace(sum, vector<uint32_t>(m + 1, 0)).first;
        }
        it->second[cnt] += 1;
    }

    unsigned hw = thread::hardware_concurrency();
    int threads = hw == 0 ? 1 : static_cast<int>(hw);
    if (right_masks < 4096) threads = 1;
    if (threads > 8) threads = 8;
    if (right_masks < static_cast<uint64_t>(threads)) {
        threads = static_cast<int>(right_masks);
    }

    vector<vector<uint64_t>> partial(threads, vector<uint64_t>(m + 1, 0));
    vector<thread> workers;
    uint64_t chunk = (right_masks + threads - 1) / threads;

    auto worker = [&](int idx, uint64_t start, uint64_t end) {
        vector<int> local_sum(d, 0);
        auto& local = partial[idx];
        for (uint64_t mask = start; mask < end; ++mask) {
            fill(local_sum.begin(), local_sum.end(), 0);
            int cnt = 0;
            for (int k = 0; k < right_n; ++k) {
                if (mask & (1ULL << k)) {
                    ++cnt;
                    const auto& pv = powers[mid + k];
                    for (int i = 0; i < d; ++i) local_sum[i] += pv[i];
                }
            }
            if (cnt > m) continue;
            for (int i = 0; i < d; ++i) local_sum[i] = -local_sum[i];
            auto it = left_sums.find(local_sum);
            if (it == left_sums.end()) continue;
            const auto& counts = it->second;
            for (int lc = 0; lc + cnt <= m; ++lc) {
                uint32_t ways = counts[lc];
                if (ways) local[lc + cnt] += ways;
            }
        }
    };

    for (int t = 0; t < threads; ++t) {
        uint64_t start = static_cast<uint64_t>(t) * chunk;
        uint64_t end = min(right_masks, start + chunk);
        workers.emplace_back(worker, t, start, end);
    }
    for (auto& th : workers) th.join();

    vector<uint64_t> total(m + 1, 0);
    for (int t = 0; t < threads; ++t) {
        for (int k = 0; k <= m; ++k) {
            total[k] += partial[t][k];
        }
    }
    return total;
}

int radical(int n) {
    int result = 1;
    int x = n;
    for (int p = 2; p * p <= x; ++p) {
        if (x % p == 0) {
            result *= p;
            while (x % p == 0) x /= p;
        }
    }
    if (x > 1) result *= x;
    return result;
}

BigInt solve_chandelier(int n, int m) {
    if (m < 0 || m > n) return BigInt(0);
    if (m == 0) return BigInt(1);

    int S = radical(n);
    int G = n / S;

    vector<uint64_t> counts = count_zero_sums(S, m);
    Poly P(m);
    for (int k = 0; k <= m; ++k) P.coeffs[k] = BigInt(counts[k]);

    Poly Final = Poly::one();
    for (int i = 0; i < G; ++i) {
        Final = Final.multiply(P, m);
    }
    return m < static_cast<int>(Final.coeffs.size()) ? Final.coeffs[m] : BigInt(0);
}

}  // namespace

int main() {
    BigInt v12 = solve_chandelier(12, 4);
    cout << "f(12, 4)   = " << v12;
    cout << (v12.equals_uint64(15) ? " [PASS]" : " [FAIL]") << '\n';

    BigInt v36 = solve_chandelier(36, 6);
    cout << "f(36, 6)   = " << v36;
    cout << (v36.equals_uint64(876) ? " [PASS]" : " [FAIL]") << '\n';

    cout << "f(360, 20) = " << solve_chandelier(360, 20) << '\n';
    return 0;
}
