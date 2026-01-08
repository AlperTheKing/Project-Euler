#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

using Real = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<300>>;

namespace {

constexpr int N = 350;
constexpr int DIM = 7;
const std::array<int, DIM> kSmallPrimes = {2, 3, 5, 7, 11, 13, 17};

struct LengthInfo {
    int len;
    std::array<int, DIM> exp;
};

struct StateSpace {
    std::array<int, DIM> dims{};
    std::array<int, DIM> mult{};
    int size = 0;
    std::vector<std::array<int, DIM>> vecs;

    void init(const std::array<int, DIM>& dims_in) {
        dims = dims_in;
        mult[DIM - 1] = 1;
        for (int i = DIM - 2; i >= 0; --i) {
            mult[i] = mult[i + 1] * dims[i + 1];
        }
        size = mult[0] * dims[0];
        vecs.resize(size);
        for (int idx = 0; idx < size; ++idx) {
            int rem = idx;
            std::array<int, DIM> v{};
            for (int i = 0; i < DIM; ++i) {
                v[i] = rem / mult[i];
                rem %= mult[i];
            }
            vecs[idx] = v;
        }
    }

    int index_of(const std::array<int, DIM>& v) const {
        int idx = 0;
        for (int i = 0; i < DIM; ++i) {
            idx += v[i] * mult[i];
        }
        return idx;
    }
};

bool dominates(const std::array<int, DIM>& a, const std::array<int, DIM>& b) {
    for (int i = 0; i < DIM; ++i) {
        if (a[i] < b[i]) return false;
    }
    return true;
}

std::array<int, DIM> clamp_exp(const std::array<int, DIM>& v, const std::array<int, DIM>& cap) {
    std::array<int, DIM> out{};
    for (int i = 0; i < DIM; ++i) {
        out[i] = std::min(v[i], cap[i]);
    }
    return out;
}

std::vector<std::vector<int>> build_allowed(const std::vector<LengthInfo>& lengths,
                                            const StateSpace& space) {
    std::vector<std::vector<int>> allowed(space.size);
    for (int idx = 0; idx < space.size; ++idx) {
        const auto& v = space.vecs[idx];
        auto& list = allowed[idx];
        for (const auto& info : lengths) {
            if (dominates(v, info.exp)) {
                list.push_back(info.len);
            }
        }
    }
    return allowed;
}

void mobius_transform_inplace(std::vector<Real>& data, const StateSpace& space, int n_max) {
    const int stride = n_max + 1;
    for (int dim = 0; dim < DIM; ++dim) {
        if (space.dims[dim] <= 1) continue;
        int step = space.mult[dim];
        int block = step * space.dims[dim];
        for (int base = 0; base < space.size; base += block) {
            for (int offset = 0; offset < step; ++offset) {
                for (int e = space.dims[dim] - 1; e >= 1; --e) {
                    int idx = base + offset + e * step;
                    int prev = idx - step;
                    Real* row = &data[idx * stride];
                    Real* prow = &data[prev * stride];
                    for (int n = 0; n <= n_max; ++n) {
                        row[n] -= prow[n];
                    }
                }
            }
        }
    }
}

template <typename Func>
void parallel_for(int total, Func func, unsigned threads) {
    if (threads <= 1 || total < 64) {
        func(0, total);
        return;
    }
    threads = std::min<unsigned>(threads, static_cast<unsigned>(total));
    int chunk = (total + static_cast<int>(threads) - 1) / static_cast<int>(threads);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        int start = static_cast<int>(t) * chunk;
        int end = std::min(total, start + chunk);
        if (start >= end) continue;
        pool.emplace_back(func, start, end);
    }
    for (auto& th : pool) th.join();
}

void compute_F(const std::vector<std::vector<int>>& allowed,
               const std::vector<Real>& inv,
               int n_max,
               std::vector<Real>& out,
               unsigned threads) {
    const int stride = n_max + 1;
    auto worker = [&](int start, int end) {
        for (int idx = start; idx < end; ++idx) {
            Real* row = &out[idx * stride];
            std::fill(row, row + stride, 0.0L);
            row[0] = 1.0L;
            const auto& list = allowed[idx];
            for (int n = 0; n < n_max; ++n) {
                Real sum = 0.0L;
                for (int len : list) {
                    if (len > n + 1) break;
                    sum += row[n + 1 - len];
                }
                row[n + 1] = sum * inv[n + 1];
            }
        }
    };
    parallel_for(static_cast<int>(allowed.size()), worker, threads);
}

void convolve_series(const Real* A, const Real* B, Real* out, int n_max) {
    const int stride = n_max + 1;
    std::fill(out, out + stride, 0.0L);
    for (int i = 0; i <= n_max; ++i) {
        Real ai = A[i];
        if (ai == 0.0L) continue;
        for (int j = 0; j <= n_max - i; ++j) {
            Real bj = B[j];
            if (bj == 0.0L) continue;
            out[i + j] += ai * bj;
        }
    }
}

std::string format_scientific_no_plus(Real value, int sig_digits) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(sig_digits - 1) << value;
    std::string s = oss.str();
    std::size_t pos = s.find('e');
    if (pos != std::string::npos && pos + 1 < s.size() && s[pos + 1] == '+') {
        s.erase(pos + 1, 1);
    }
    return s;
}

} // namespace

int main() {
    std::array<int, DIM> max_exp{};
    for (int i = 0; i < DIM; ++i) {
        int p = kSmallPrimes[i];
        int e = 0;
        long long v = p;
        while (v <= N) {
            ++e;
            v *= p;
        }
        max_exp[i] = e;
    }

    std::vector<Real> inv(N + 1, 0.0L);
    for (int i = 1; i <= N; ++i) inv[i] = Real(1) / static_cast<Real>(i);

    std::map<int, std::vector<LengthInfo>> groups;
    std::array<int, DIM> max_local{};
    for (int k = 1; k <= N; ++k) {
        int tmp = k;
        std::array<int, DIM> exp{};
        for (int i = 0; i < DIM; ++i) {
            int p = kSmallPrimes[i];
            while (tmp % p == 0) {
                tmp /= p;
                ++exp[i];
            }
        }
        int p_large = tmp; // 1 or prime >= 19
        groups[p_large].push_back({k, exp});
        if (p_large != 1) {
            for (int i = 0; i < DIM; ++i) {
                max_local[i] = std::max(max_local[i], exp[i]);
            }
        }
    }
    for (auto& [key, vec] : groups) {
        std::sort(vec.begin(), vec.end(), [](const LengthInfo& a, const LengthInfo& b) {
            return a.len < b.len;
        });
    }

    StateSpace global_space;
    global_space.init({
        max_exp[0] + 1, max_exp[1] + 1, max_exp[2] + 1, max_exp[3] + 1,
        max_exp[4] + 1, max_exp[5] + 1, max_exp[6] + 1});

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    // F1(E): permutations using only small-prime lengths, with max exponents <= E.
    const auto& group1 = groups[1];
    auto allowed_global = build_allowed(group1, global_space);
    std::vector<Real> F_global(global_space.size * (N + 1), 0.0L);
    compute_F(allowed_global, inv, N, F_global, threads);

    // Build local state space for large-prime lengths (s <= 18).
    StateSpace local_space;
    local_space.init({
        max_local[0] + 1, max_local[1] + 1, max_local[2] + 1, max_local[3] + 1,
        max_local[4] + 1, max_local[5] + 1, max_local[6] + 1});

    std::vector<Real> F_large(local_space.size * (N + 1), 0.0L);
    for (int idx = 0; idx < local_space.size; ++idx) {
        F_large[idx * (N + 1)] = 1.0L;
    }

    std::vector<Real> Fp(local_space.size * (N + 1), 0.0L);
    std::vector<Real> new_large(local_space.size * (N + 1), 0.0L);

    for (const auto& [p, lengths] : groups) {
        if (p == 1) continue;
        auto allowed_local = build_allowed(lengths, local_space);
        std::fill(Fp.begin(), Fp.end(), 0.0L);
        compute_F(allowed_local, inv, N, Fp, 1);

        Real p2 = static_cast<Real>(p) * static_cast<Real>(p);
        for (int idx = 0; idx < local_space.size; ++idx) {
            Real* row = &Fp[idx * (N + 1)];
            row[0] = 1.0L;
            for (int n = 1; n <= N; ++n) {
                row[n] *= p2;
            }
        }

        std::fill(new_large.begin(), new_large.end(), 0.0L);
        auto worker = [&](int start, int end) {
            std::vector<Real> temp(N + 1);
            for (int idx = start; idx < end; ++idx) {
                const Real* A = &F_large[idx * (N + 1)];
                const Real* B = &Fp[idx * (N + 1)];
                convolve_series(A, B, temp.data(), N);
                Real* out = &new_large[idx * (N + 1)];
                std::copy(temp.begin(), temp.end(), out);
            }
        };
        parallel_for(local_space.size, worker, threads);
        F_large.swap(new_large);
    }

    // Convolve F1 with F_large (clamped) to get F_total in-place.
    std::vector<int> local_index_for_global(global_space.size, 0);
    for (int idx = 0; idx < global_space.size; ++idx) {
        auto v = global_space.vecs[idx];
        auto vloc = clamp_exp(v, max_local);
        local_index_for_global[idx] = local_space.index_of(vloc);
    }

    auto conv_worker = [&](int start, int end) {
        std::vector<Real> temp(N + 1);
        for (int idx = start; idx < end; ++idx) {
            Real* row = &F_global[idx * (N + 1)];
            const Real* B = &F_large[local_index_for_global[idx] * (N + 1)];
            convolve_series(row, B, temp.data(), N);
            std::copy(temp.begin(), temp.end(), row);
        }
    };
    parallel_for(global_space.size, conv_worker, threads);

    // Convert max<=E to exact max=E.
    mobius_transform_inplace(F_global, global_space, N);

    // Precompute L_small^2 for each global state.
    std::array<std::vector<Real>, DIM> prime_pows;
    for (int i = 0; i < DIM; ++i) {
        prime_pows[i].resize(max_exp[i] + 1, 1.0L);
        for (int e = 1; e <= max_exp[i]; ++e) {
            prime_pows[i][e] = prime_pows[i][e - 1] * static_cast<Real>(kSmallPrimes[i]);
        }
    }
    std::vector<Real> lcm_sq(global_space.size, 1.0L);
    for (int idx = 0; idx < global_space.size; ++idx) {
        Real val = 1.0L;
        const auto& v = global_space.vecs[idx];
        for (int i = 0; i < DIM; ++i) {
            val *= prime_pows[i][v[i]];
        }
        lcm_sq[idx] = val * val;
    }

    auto sum_for_n = [&](int target) {
        Real sum = 0.0L;
        for (int idx = 0; idx < global_space.size; ++idx) {
            sum += F_global[idx * (N + 1) + target] * lcm_sq[idx];
        }
        return sum;
    };

    Real g3 = sum_for_n(3);
    Real g5 = sum_for_n(5);
    Real g20 = sum_for_n(20);
    Real g350 = sum_for_n(350);

    auto check = [](const Real& got, const Real& expected) {
        using boost::multiprecision::abs;
        Real diff = abs(got - expected);
        Real tol = std::max(Real("1e-6"), abs(expected) * Real("1e-9"));
        return diff <= tol;
    };

    std::cout << "--- Validation Checkpoints ---\n";
    std::cout << std::fixed << std::setprecision(12);
    const Real expected3 = Real(31) / Real(6);
    const Real expected5 = Real(2081) / Real(120);
    const Real expected20 = Real("5106.136147");
    std::cout << "g(3) = " << g3 << (check(g3, expected3) ? " [PASS]" : " [FAIL]") << "\n";
    std::cout << "g(5) = " << g5 << (check(g5, expected5) ? " [PASS]" : " [FAIL]") << "\n";
    std::cout << "g(20) = " << g20 << (check(g20, expected20) ? " [PASS]" : " [FAIL]") << "\n";

    std::cout << "\n--- Final Solution ---\n";
    std::cout << format_scientific_no_plus(g350, 10) << "\n";

    return 0;
}
