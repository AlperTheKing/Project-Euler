#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1000000000LL;
constexpr int kSize = 18;

using Matrix = std::array<std::array<i64, kSize>, kSize>;
using Vector = std::array<i64, kSize>;

Matrix zero_matrix() {
    Matrix m{};
    for (int i = 0; i < kSize; ++i) {
        for (int j = 0; j < kSize; ++j) {
            m[i][j] = 0;
        }
    }
    return m;
}

Matrix identity_matrix() {
    Matrix m = zero_matrix();
    for (int i = 0; i < kSize; ++i) {
        m[i][i] = 1;
    }
    return m;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
    Matrix c = zero_matrix();
    for (int i = 0; i < kSize; ++i) {
        for (int k = 0; k < kSize; ++k) {
            const i64 aik = a[i][k];
            if (aik == 0) {
                continue;
            }
            for (int j = 0; j < kSize; ++j) {
                if (b[k][j] == 0) {
                    continue;
                }
                c[i][j] = (c[i][j] + static_cast<__int128>(aik) * b[k][j]) % kMod;
            }
        }
    }
    return c;
}

Vector apply_matrix_vector(const Matrix& a, const Vector& v) {
    Vector out{};
    for (int i = 0; i < kSize; ++i) {
        __int128 sum = 0;
        for (int j = 0; j < kSize; ++j) {
            if (a[i][j] == 0 || v[j] == 0) {
                continue;
            }
            sum += static_cast<__int128>(a[i][j]) * v[j];
        }
        out[i] = static_cast<i64>(sum % kMod);
    }
    return out;
}

Matrix build_transition() {
    Matrix t = zero_matrix();

    // a_{n+1} = a_n + a_{n-1} + ... + a_{n-8}
    for (int i = 0; i < 9; ++i) {
        t[0][i] = 1;
    }
    for (int i = 1; i < 9; ++i) {
        t[i][i - 1] = 1;
    }

    // f_{n+1} = 10*(f_n+...+f_{n-8}) + 1*a_n + 2*a_{n-1} + ... + 9*a_{n-8}
    for (int i = 0; i < 9; ++i) {
        t[9][i] = i + 1;
        t[9][9 + i] = 10;
    }
    for (int i = 10; i < 18; ++i) {
        t[i][i - 1] = 1;
    }

    return t;
}

std::vector<Matrix> precompute_powers(const Matrix& base) {
    std::vector<Matrix> powers;
    powers.reserve(64);
    powers.push_back(base);
    for (int i = 1; i < 64; ++i) {
        powers.push_back(multiply(powers.back(), powers.back()));
    }
    return powers;
}

i64 f_of_n(const u64 n, const std::vector<Matrix>& powers) {
    // state at n=0:
    // [a_0, a_{-1},...,a_{-8}, f_0, f_{-1},...,f_{-8}] with a_0=1, others 0.
    Vector state{};
    for (int i = 0; i < kSize; ++i) {
        state[i] = 0;
    }
    state[0] = 1;

    u64 e = n;
    int bit = 0;
    while (e > 0) {
        if ((e & 1ULL) != 0ULL) {
            state = apply_matrix_vector(powers[static_cast<std::size_t>(bit)], state);
        }
        e >>= 1ULL;
        ++bit;
    }

    // f_n is at index 9.
    return state[9];
}

bool run_checkpoints(const std::vector<Matrix>& powers) {
    if (f_of_n(5ULL, powers) != 17891LL) {
        std::cerr << "Checkpoint failed: f(5)\n";
        return false;
    }

    // Cross-check with direct DP for small n.
    std::vector<i64> a(30, 0);
    std::vector<i64> f(30, 0);
    a[0] = 1;
    for (int n = 1; n < 30; ++n) {
        for (int d = 1; d <= 9; ++d) {
            if (n - d < 0) {
                continue;
            }
            a[n] = (a[n] + a[n - d]) % kMod;
            f[n] = (f[n] + 10LL * f[n - d] + static_cast<i64>(d) * a[n - d]) % kMod;
        }
    }
    for (int n = 1; n < 30; ++n) {
        if (f_of_n(static_cast<u64>(n), powers) != f[n]) {
            std::cerr << "Checkpoint failed: matrix/DP mismatch at n=" << n << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    const Matrix transition = build_transition();
    const std::vector<Matrix> powers = precompute_powers(transition);

    if (!skip_checkpoints && !run_checkpoints(powers)) {
        return 2;
    }

    i64 answer = 0;
    u64 p13 = 1ULL;
    for (int i = 1; i <= 17; ++i) {
        p13 *= 13ULL;
        answer += f_of_n(p13, powers);
        answer %= kMod;
    }

    std::cout << answer << '\n';
    return 0;
}
