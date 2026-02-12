#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

constexpr i64 MOD = 1'000'000'009LL;
constexpr int D = 10;

using Matrix = std::array<std::array<i64, D>, D>;

Matrix mat_mul(const Matrix& A, const Matrix& B) {
    Matrix C{};
    for (int i = 0; i < D; ++i) {
        for (int k = 0; k < D; ++k) {
            const i64 aik = A[i][k];
            if (aik == 0) {
                continue;
            }
            for (int j = 0; j < D; ++j) {
                C[i][j] = static_cast<i64>((C[i][j] + static_cast<__int128>(aik) * B[k][j]) % MOD);
            }
        }
    }
    return C;
}

Matrix mat_pow(Matrix base, u64 exp) {
    Matrix result{};
    for (int i = 0; i < D; ++i) {
        result[i][i] = 1;
    }

    while (exp > 0) {
        if (exp & 1ULL) {
            result = mat_mul(base, result);
        }
        base = mat_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

std::array<i64, D> digit_counts(u64 M, u64 p10) {
    const u64 N = M + 1;
    const u64 block = 10ULL * p10;
    const u64 full = N / block;
    const u64 rem = N % block;

    std::array<i64, D> c{};
    for (int d = 0; d < D; ++d) {
        u64 cnt = full * p10;
        const u64 shift = static_cast<u64>(d) * p10;
        if (rem > shift) {
            const u64 extra = rem - shift;
            cnt += (extra < p10) ? extra : p10;
        }
        c[d] = static_cast<i64>(cnt % MOD);
    }
    return c;
}

i64 contribution_for_position(u64 R, u64 M, u64 p10) {
    const auto c = digit_counts(M, p10);

    Matrix T{};
    for (int from = 0; from < D; ++from) {
        for (int d = 0; d < D; ++d) {
            const int to = (from * d) % 10;
            T[to][from] += c[d];
            if (T[to][from] >= MOD) {
                T[to][from] -= MOD;
            }
        }
    }

    const Matrix P = mat_pow(T, R);

    i64 s = 0;
    for (int r = 0; r < D; ++r) {
        s = static_cast<i64>((s + static_cast<__int128>(r) * P[r][1]) % MOD);
    }
    return s;
}

i64 F_mod(u64 R, u64 M) {
    i64 ans = 0;
    i64 place_mod = 1;

    for (u64 p10 = 1; p10 <= M; p10 *= 10ULL) {
        const i64 s = contribution_for_position(R, M, p10);
        ans = static_cast<i64>((ans + static_cast<__int128>(place_mod) * s) % MOD);
        place_mod = (place_mod * 10LL) % MOD;
    }

    return ans;
}

}  // namespace

int main() {
    assert(F_mod(2, 7) == 204);
    assert(F_mod(23, 76) == 5'870'548);

    std::cout << F_mod(234'567, 765'432) << '\n';
    return 0;
}
