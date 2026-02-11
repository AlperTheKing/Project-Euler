#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kPhi = kMod - 1ULL;
constexpr int kWords = 4;  // enough for width <= 199

struct Mat {
    int w = 0;
    int words = 0;
    std::vector<std::array<u64, kWords>> row;

    Mat() = default;
    explicit Mat(const int width) : w(width), words((width + 63) / 64), row(static_cast<std::size_t>(width)) {
        for (auto& r : row) {
            r.fill(0ULL);
        }
    }
};

inline void set_bit(std::array<u64, kWords>& r, const int c) {
    r[static_cast<std::size_t>(c >> 6)] |= (1ULL << (c & 63));
}

inline bool get_bit(const std::array<u64, kWords>& r, const int c) {
    return ((r[static_cast<std::size_t>(c >> 6)] >> (c & 63)) & 1ULL) != 0ULL;
}

Mat identity(const int w) {
    Mat I(w);
    for (int i = 0; i < w; ++i) {
        set_bit(I.row[static_cast<std::size_t>(i)], i);
    }
    return I;
}

Mat zero_mat(const int w) {
    return Mat(w);
}

Mat build_H(const int w) {
    Mat H(w);
    for (int i = 0; i < w; ++i) {
        set_bit(H.row[static_cast<std::size_t>(i)], i);
        if (i > 0) {
            set_bit(H.row[static_cast<std::size_t>(i)], i - 1);
        }
        if (i + 1 < w) {
            set_bit(H.row[static_cast<std::size_t>(i)], i + 1);
        }
    }
    return H;
}

Mat add(const Mat& A, const Mat& B) {
    Mat C(A.w);
    for (int i = 0; i < A.w; ++i) {
        for (int k = 0; k < A.words; ++k) {
            C.row[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                A.row[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] ^
                B.row[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
        }
    }
    return C;
}

Mat mul(const Mat& A, const Mat& B) {
    Mat C(A.w);
    for (int i = 0; i < A.w; ++i) {
        std::array<u64, kWords> acc{};
        acc.fill(0ULL);

        for (int wb = 0; wb < A.words; ++wb) {
            u64 bits = A.row[static_cast<std::size_t>(i)][static_cast<std::size_t>(wb)];
            while (bits != 0ULL) {
                const int b = __builtin_ctzll(bits);
                const int col = (wb << 6) + b;
                if (col < A.w) {
                    const auto& brow = B.row[static_cast<std::size_t>(col)];
                    for (int k = 0; k < A.words; ++k) {
                        acc[static_cast<std::size_t>(k)] ^= brow[static_cast<std::size_t>(k)];
                    }
                }
                bits &= bits - 1ULL;
            }
        }

        C.row[static_cast<std::size_t>(i)] = acc;
    }
    return C;
}

Mat left_mul_H(const Mat& P) {
    Mat R(P.w);
    for (int i = 0; i < P.w; ++i) {
        std::array<u64, kWords> acc = P.row[static_cast<std::size_t>(i)];
        if (i > 0) {
            for (int k = 0; k < P.words; ++k) {
                acc[static_cast<std::size_t>(k)] ^= P.row[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(k)];
            }
        }
        if (i + 1 < P.w) {
            for (int k = 0; k < P.words; ++k) {
                acc[static_cast<std::size_t>(k)] ^= P.row[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(k)];
            }
        }
        R.row[static_cast<std::size_t>(i)] = acc;
    }
    return R;
}

int rank_gf2(Mat A) {
    int r = 0;
    for (int c = 0; c < A.w && r < A.w; ++c) {
        int piv = -1;
        for (int i = r; i < A.w; ++i) {
            if (get_bit(A.row[static_cast<std::size_t>(i)], c)) {
                piv = i;
                break;
            }
        }
        if (piv == -1) {
            continue;
        }

        if (piv != r) {
            std::swap(A.row[static_cast<std::size_t>(piv)], A.row[static_cast<std::size_t>(r)]);
        }

        for (int i = r + 1; i < A.w; ++i) {
            if (!get_bit(A.row[static_cast<std::size_t>(i)], c)) {
                continue;
            }
            for (int k = 0; k < A.words; ++k) {
                A.row[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] ^=
                    A.row[static_cast<std::size_t>(r)][static_cast<std::size_t>(k)];
            }
        }
        ++r;
    }
    return r;
}

u64 mod_pow2(u64 e) {
    u64 base = 2ULL;
    u64 out = 1ULL;
    while (e > 0ULL) {
        if (e & 1ULL) {
            out = static_cast<u64>((static_cast<u128>(out) * static_cast<u128>(base)) % kMod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * static_cast<u128>(base)) % kMod);
        e >>= 1ULL;
    }
    return out;
}

u64 F_value(const int w, const u64 h) {
    const Mat H = build_H(w);
    Mat U_prev = zero_mat(w);   // U_0
    Mat U_cur = identity(w);    // U_1

    for (u64 i = 0; i < h; ++i) {
        Mat U_next = add(mul(H, U_cur), U_prev);  // U_{i+2}=H U_{i+1}+U_i
        U_prev = std::move(U_cur);
        U_cur = std::move(U_next);
    }

    const int rk = rank_gf2(U_cur);  // U_{h+1}
    const int nullity = w - rk;

    u64 exp = static_cast<u64>((static_cast<u128>(w % static_cast<int>(kPhi)) *
                                static_cast<u128>(h % kPhi)) %
                               static_cast<u128>(kPhi));
    exp = (exp + kPhi - static_cast<u64>(nullity % static_cast<int>(kPhi))) % kPhi;
    return mod_pow2(exp);
}

u64 S_value(const int w, const int n) {
    const Mat H = build_H(w);
    Mat A_prev = identity(w);  // U_{f1}
    Mat B_prev = H;            // U_{f1+1}
    Mat A_cur = identity(w);   // U_{f2}
    Mat B_cur = H;             // U_{f2+1}

    u64 fib_prev = 1ULL;  // f1
    u64 fib_cur = 1ULL;   // f2

    auto term_from_B = [&](const Mat& Bk, const u64 fib_k) -> u64 {
        const int rk = rank_gf2(Bk);
        const int nullity = w - rk;
        u64 exp = static_cast<u64>((static_cast<u128>(w % static_cast<int>(kPhi)) *
                                    static_cast<u128>(fib_k % kPhi)) %
                                   static_cast<u128>(kPhi));
        exp = (exp + kPhi - static_cast<u64>(nullity % static_cast<int>(kPhi))) % kPhi;
        return mod_pow2(exp);
    };

    u64 sum = 0ULL;
    if (n >= 1) {
        sum += term_from_B(B_prev, fib_prev);
        if (sum >= kMod) {
            sum -= kMod;
        }
    }
    if (n >= 2) {
        sum += term_from_B(B_cur, fib_cur);
        if (sum >= kMod) {
            sum -= kMod;
        }
    }

    for (int k = 3; k <= n; ++k) {
        const Mat P = mul(A_cur, A_prev);
        const Mat Q = mul(A_cur, B_prev);
        const Mat R = mul(B_cur, A_prev);
        const Mat S = mul(B_cur, B_prev);
        const Mat HP = left_mul_H(P);

        Mat A_next = add(add(Q, R), HP);
        Mat B_next = add(S, P);

        const u64 fib_next = (fib_cur + fib_prev) % kPhi;
        const u64 term = term_from_B(B_next, fib_next);
        sum += term;
        if (sum >= kMod) {
            sum -= kMod;
        }

        A_prev = std::move(A_cur);
        B_prev = std::move(B_cur);
        A_cur = std::move(A_next);
        B_cur = std::move(B_next);
        fib_prev = fib_cur;
        fib_cur = fib_next;
    }

    return sum;
}

}  // namespace

int main() {
    assert(F_value(1, 2) == 2ULL);
    assert(F_value(3, 3) == 512ULL);
    assert(F_value(4, 4) == 4096ULL);
    assert(F_value(7, 11) == 270016253ULL);

    assert(S_value(3, 3) == 32ULL);
    assert(S_value(4, 5) == 1052960ULL);
    assert(S_value(5, 7) == 346547294ULL);

    std::cout << S_value(199, 199) << '\n';
    return 0;
}

