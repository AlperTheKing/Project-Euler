#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <thread>
#include <vector>

using std::int64_t;

namespace {

struct Coeff {
    int c00;
    int c11;
    int c22;
    int cab;
    int cac;
    int cbc;
};

class Bitset {
public:
    explicit Bitset(size_t bits) : data_((bits + 63) / 64, 0) {}

    inline void set(size_t idx) {
        data_[idx >> 6] |= 1ULL << (idx & 63);
    }

    void or_with(const Bitset& other) {
        const size_t n = data_.size();
        for (size_t i = 0; i < n; ++i) {
            data_[i] |= other.data_[i];
        }
    }

    uint64_t count() const {
        uint64_t total = 0;
        for (uint64_t w : data_) {
            total += static_cast<uint64_t>(__builtin_popcountll(w));
        }
        return total;
    }

private:
    std::vector<uint64_t> data_;
};

std::vector<Coeff> build_coeffs() {
    std::set<std::array<int, 6>> uniq;
    for (int bits = 0; bits < (1 << 9); ++bits) {
        int B[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                B[i][j] = (bits >> (i * 3 + j)) & 1;
            }
        }
        int rows[3] = {
            B[0][0] | (B[0][1] << 1) | (B[0][2] << 2),
            B[1][0] | (B[1][1] << 1) | (B[1][2] << 2),
            B[2][0] | (B[2][1] << 1) | (B[2][2] << 2),
        };
        int cols[3] = {
            B[0][0] | (B[1][0] << 1) | (B[2][0] << 2),
            B[0][1] | (B[1][1] << 1) | (B[2][1] << 2),
            B[0][2] | (B[1][2] << 1) | (B[2][2] << 2),
        };
        if (!(rows[0] != rows[1] && rows[0] != rows[2] && rows[1] != rows[2])) {
            continue;
        }
        if (!(cols[0] != cols[1] && cols[0] != cols[2] && cols[1] != cols[2])) {
            continue;
        }
        std::array<int, 3> rsorted = {rows[0], rows[1], rows[2]};
        std::array<int, 3> csorted = {cols[0], cols[1], cols[2]};
        std::sort(rsorted.begin(), rsorted.end());
        std::sort(csorted.begin(), csorted.end());
        if (rsorted != csorted) {
            continue;
        }
        const int c00 = B[0][0];
        const int c11 = B[1][1];
        const int c22 = B[2][2];
        const int cab = B[0][1] + B[1][0];
        const int cac = B[0][2] + B[2][0];
        const int cbc = B[1][2] + B[2][1];
        uniq.insert({c00, c11, c22, cab, cac, cbc});
    }
    std::vector<Coeff> coeffs;
    coeffs.reserve(uniq.size());
    for (const auto& v : uniq) {
        coeffs.push_back({v[0], v[1], v[2], v[3], v[4], v[5]});
    }
    return coeffs;
}

int64_t compute_C(int n, unsigned threads, const std::vector<Coeff>& coeffs) {
    const int64_t n2 = static_cast<int64_t>(n) * n;
    const size_t bits = static_cast<size_t>(n2) + 1;

    Bitset le3(bits);
    for (int a = 0; a <= n; ++a) {
        for (int b = 0; b <= n; ++b) {
            const int64_t k = static_cast<int64_t>(a) * b;
            le3.set(static_cast<size_t>(k));
            le3.set(static_cast<size_t>(n2 - k));
        }
    }

    std::vector<int64_t> sq(n + 1, 0);
    for (int i = 0; i <= n; ++i) {
        sq[i] = static_cast<int64_t>(i) * i;
    }

    const unsigned tcount = std::max(1u, threads);
    if (tcount == 1) {
        Bitset local(bits);
        for (int a = 0; a <= n; ++a) {
            const int bmax = (n - a) / 2;
            if (bmax < a) continue;
            const int64_t aa = sq[a];
            for (int b = a; b <= bmax; ++b) {
                const int c = n - a - b;
                const int64_t bb = sq[b];
                const int64_t cc = sq[c];
                const int64_t ab = static_cast<int64_t>(a) * b;
                const int64_t ac = static_cast<int64_t>(a) * c;
                const int64_t bc = static_cast<int64_t>(b) * c;
                for (const auto& cf : coeffs) {
                    const int64_t k = cf.c00 * aa + cf.c11 * bb + cf.c22 * cc
                                      + cf.cab * ab + cf.cac * ac + cf.cbc * bc;
                    local.set(static_cast<size_t>(k));
                }
            }
        }
        le3.or_with(local);
    } else {
        const unsigned cap_threads = std::min(threads, 8u);
        std::vector<Bitset> locals;
        locals.reserve(cap_threads);
        for (unsigned t = 0; t < cap_threads; ++t) {
            locals.emplace_back(bits);
        }
        std::vector<std::thread> workers;
        workers.reserve(cap_threads);
        for (unsigned t = 0; t < cap_threads; ++t) {
            workers.emplace_back([&, t]() {
                Bitset& local = locals[t];
                for (int a = static_cast<int>(t); a <= n; a += static_cast<int>(cap_threads)) {
                    const int bmax = (n - a) / 2;
                    if (bmax < a) continue;
                    const int64_t aa = sq[a];
                    for (int b = a; b <= bmax; ++b) {
                        const int c = n - a - b;
                        const int64_t bb = sq[b];
                        const int64_t cc = sq[c];
                        const int64_t ab = static_cast<int64_t>(a) * b;
                        const int64_t ac = static_cast<int64_t>(a) * c;
                        const int64_t bc = static_cast<int64_t>(b) * c;
                        for (const auto& cf : coeffs) {
                            const int64_t k = cf.c00 * aa + cf.c11 * bb + cf.c22 * cc
                                              + cf.cab * ab + cf.cac * ac + cf.cbc * bc;
                            local.set(static_cast<size_t>(k));
                        }
                    }
                }
            });
        }
        for (auto& th : workers) th.join();
        for (const auto& local : locals) {
            le3.or_with(local);
        }
    }

    Bitset s2(bits);
    for (int a = 0; a <= n; ++a) {
        const int b = n - a;
        const int64_t w1 = sq[a];
        const int64_t w2 = sq[b];
        const int64_t w3 = 2LL * a * b;
        const int64_t vals[7] = {
            w1, w2, w3, w1 + w2, w1 + w3, w2 + w3, w1 + w2 + w3
        };
        for (int i = 0; i < 7; ++i) {
            const int64_t k = vals[i];
            if (k == 0 || k == n2) continue;
            s2.set(static_cast<size_t>(k));
        }
    }

    const uint64_t count_le3 = le3.count();
    const uint64_t count_s2 = s2.count();
    const uint64_t total = static_cast<uint64_t>(n2) + 1;
    const uint64_t n4 = total - count_le3;
    const uint64_t C = 3 * total - 4 - count_s2 + n4;
    return static_cast<int64_t>(C);
}

}  // namespace

int main() {
    const auto coeffs = build_coeffs();

    if (compute_C(5, 1, coeffs) != 64) {
        std::cerr << "Validation failed for C(5)." << '\n';
        return 1;
    }
    if (compute_C(10, 1, coeffs) != 274) {
        std::cerr << "Validation failed for C(10)." << '\n';
        return 1;
    }
    if (compute_C(20, 1, coeffs) != 1150) {
        std::cerr << "Validation failed for C(20)." << '\n';
        return 1;
    }

    const int n = 10000;
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = std::min(threads, 8u);
    std::cout << compute_C(n, threads, coeffs) << '\n';
    return 0;
}
