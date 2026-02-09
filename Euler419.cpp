#include <array>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr int K = 92;
constexpr u64 MOD = 1ULL << 30;
constexpr u64 MASK = MOD - 1ULL;
constexpr u64 TARGET_N = 1000000000000ULL;

struct Options {
    u64 n = TARGET_N;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << "\n";
        return false;
    }
    return options.n >= 1;
}

const std::array<const char*, K> kElementStrings = {
    "22",  //  1 H
    "13112221133211322112211213322112",  //  2 He
    "312211322212221121123222112",  //  3 Li
    "111312211312113221133211322112211213322112",  //  4 Be
    "1321132122211322212221121123222112",  //  5 B
    "3113112211322112211213322112",  //  6 C
    "111312212221121123222112",  //  7 N
    "132112211213322112",  //  8 O
    "31121123222112",  //  9 F
    "111213322112",  // 10 Ne
    "123222112",  // 11 Na
    "3113322112",  // 12 Mg
    "1113222112",  // 13 Al
    "1322112",  // 14 Si
    "311311222112",  // 15 P
    "1113122112",  // 16 S
    "132112",  // 17 Cl
    "3112",  // 18 Ar
    "1112",  // 19 K
    "12",  // 20 Ca
    "3113112221133112",  // 21 Sc
    "11131221131112",  // 22 Ti
    "13211312",  // 23 V
    "31132",  // 24 Cr
    "111311222112",  // 25 Mn
    "13122112",  // 26 Fe
    "32112",  // 27 Co
    "11133112",  // 28 Ni
    "131112",  // 29 Cu
    "312",  // 30 Zn
    "13221133122211332",  // 31 Ga
    "31131122211311122113222",  // 32 Ge
    "11131221131211322113322112",  // 33 As
    "13211321222113222112",  // 34 Se
    "3113112211322112",  // 35 Br
    "11131221222112",  // 36 Kr
    "1321122112",  // 37 Rb
    "3112112",  // 38 Sr
    "1112133",  // 39 Y
    "12322211331222113112211",  // 40 Zr
    "1113122113322113111221131221",  // 41 Nb
    "13211322211312113211",  // 42 Mo
    "311322113212221",  // 43 Tc
    "132211331222113112211",  // 44 Ru
    "311311222113111221131221",  // 45 Rh
    "111312211312113211",  // 46 Pd
    "132113212221",  // 47 Ag
    "3113112211",  // 48 Cd
    "11131221",  // 49 In
    "13211",  // 50 Sn
    "3112221",  // 51 Sb
    "1322113312211",  // 52 Te
    "311311222113111221",  // 53 I
    "11131221131211",  // 54 Xe
    "13211321",  // 55 Cs
    "311311",  // 56 Ba
    "11131",  // 57 La
    "1321133112",  // 58 Ce
    "31131112",  // 59 Pr
    "111312",  // 60 Nd
    "132",  // 61 Pm
    "311332",  // 62 Sm
    "1113222",  // 63 Eu
    "13221133112",  // 64 Gd
    "3113112221131112",  // 65 Tb
    "111312211312",  // 66 Dy
    "1321132",  // 67 Ho
    "311311222",  // 68 Er
    "11131221133112",  // 69 Tm
    "1321131112",  // 70 Yb
    "311312",  // 71 Lu
    "11132",  // 72 Hf
    "13112221133211322112211213322113",  // 73 Ta
    "312211322212221121123222113",  // 74 W
    "111312211312113221133211322112211213322113",  // 75 Re
    "1321132122211322212221121123222113",  // 76 Os
    "3113112211322112211213322113",  // 77 Ir
    "111312212221121123222113",  // 78 Pt
    "132112211213322113",  // 79 Au
    "31121123222113",  // 80 Hg
    "111213322113",  // 81 Tl
    "123222113",  // 82 Pb
    "3113322113",  // 83 Bi
    "1113222113",  // 84 Po
    "1322113",  // 85 At
    "311311222113",  // 86 Rn
    "1113122113",  // 87 Fr
    "132113",  // 88 Ra
    "3113",  // 89 Ac
    "1113",  // 90 Th
    "13",  // 91 Pa
    "3"  // 92 U
};

const std::array<std::array<int, 6>, K> kDecay = {
    {
        { 0, -1, -1, -1, -1, -1},  //  1 H
        {71, 90,  0, 19,  2, -1},  //  2 He
        { 1, -1, -1, -1, -1, -1},  //  3 Li
        {31, 19,  2, -1, -1, -1},  //  4 Be
        { 3, -1, -1, -1, -1, -1},  //  5 B
        { 4, -1, -1, -1, -1, -1},  //  6 C
        { 5, -1, -1, -1, -1, -1},  //  7 N
        { 6, -1, -1, -1, -1, -1},  //  8 O
        { 7, -1, -1, -1, -1, -1},  //  9 F
        { 8, -1, -1, -1, -1, -1},  // 10 Ne
        { 9, -1, -1, -1, -1, -1},  // 11 Na
        {60, 10, -1, -1, -1, -1},  // 12 Mg
        {11, -1, -1, -1, -1, -1},  // 13 Al
        {12, -1, -1, -1, -1, -1},  // 14 Si
        {66, 13, -1, -1, -1, -1},  // 15 P
        {14, -1, -1, -1, -1, -1},  // 16 S
        {15, -1, -1, -1, -1, -1},  // 17 Cl
        {16, -1, -1, -1, -1, -1},  // 18 Ar
        {17, -1, -1, -1, -1, -1},  // 19 K
        {18, -1, -1, -1, -1, -1},  // 20 Ca
        {66, 90,  0, 19, 26, -1},  // 21 Sc
        {20, -1, -1, -1, -1, -1},  // 22 Ti
        {21, -1, -1, -1, -1, -1},  // 23 V
        {22, -1, -1, -1, -1, -1},  // 24 Cr
        {23, 13, -1, -1, -1, -1},  // 25 Mn
        {24, -1, -1, -1, -1, -1},  // 26 Fe
        {25, -1, -1, -1, -1, -1},  // 27 Co
        {29, 26, -1, -1, -1, -1},  // 28 Ni
        {27, -1, -1, -1, -1, -1},  // 29 Cu
        {28, -1, -1, -1, -1, -1},  // 30 Zn
        {62, 19, 88,  0, 19, 29},  // 31 Ga
        {66, 30, -1, -1, -1, -1},  // 32 Ge
        {31, 10, -1, -1, -1, -1},  // 33 As
        {32, -1, -1, -1, -1, -1},  // 34 Se
        {33, -1, -1, -1, -1, -1},  // 35 Br
        {34, -1, -1, -1, -1, -1},  // 36 Kr
        {35, -1, -1, -1, -1, -1},  // 37 Rb
        {36, -1, -1, -1, -1, -1},  // 38 Sr
        {37, 91, -1, -1, -1, -1},  // 39 Y
        {38,  0, 19, 42, -1, -1},  // 40 Zr
        {67, 39, -1, -1, -1, -1},  // 41 Nb
        {40, -1, -1, -1, -1, -1},  // 42 Mo
        {41, -1, -1, -1, -1, -1},  // 43 Tc
        {62, 19, 42, -1, -1, -1},  // 44 Ru
        {66, 43, -1, -1, -1, -1},  // 45 Rh
        {44, -1, -1, -1, -1, -1},  // 46 Pd
        {45, -1, -1, -1, -1, -1},  // 47 Ag
        {46, -1, -1, -1, -1, -1},  // 48 Cd
        {47, -1, -1, -1, -1, -1},  // 49 In
        {48, -1, -1, -1, -1, -1},  // 50 Sn
        {60, 49, -1, -1, -1, -1},  // 51 Sb
        {62, 19, 50, -1, -1, -1},  // 52 Te
        {66, 51, -1, -1, -1, -1},  // 53 I
        {52, -1, -1, -1, -1, -1},  // 54 Xe
        {53, -1, -1, -1, -1, -1},  // 55 Cs
        {54, -1, -1, -1, -1, -1},  // 56 Ba
        {55, -1, -1, -1, -1, -1},  // 57 La
        {56,  0, 19, 26, -1, -1},  // 58 Ce
        {57, -1, -1, -1, -1, -1},  // 59 Pr
        {58, -1, -1, -1, -1, -1},  // 60 Nd
        {59, -1, -1, -1, -1, -1},  // 61 Pm
        {60, 19, 29, -1, -1, -1},  // 62 Sm
        {61, -1, -1, -1, -1, -1},  // 63 Eu
        {62, 19, 26, -1, -1, -1},  // 64 Gd
        {66, 63, -1, -1, -1, -1},  // 65 Tb
        {64, -1, -1, -1, -1, -1},  // 66 Dy
        {65, -1, -1, -1, -1, -1},  // 67 Ho
        {66, 60, -1, -1, -1, -1},  // 68 Er
        {67, 19, 26, -1, -1, -1},  // 69 Tm
        {68, -1, -1, -1, -1, -1},  // 70 Yb
        {69, -1, -1, -1, -1, -1},  // 71 Lu
        {70, -1, -1, -1, -1, -1},  // 72 Hf
        {71, 90,  0, 19, 73, -1},  // 73 Ta
        {72, -1, -1, -1, -1, -1},  // 74 W
        {31, 19, 73, -1, -1, -1},  // 75 Re
        {74, -1, -1, -1, -1, -1},  // 76 Os
        {75, -1, -1, -1, -1, -1},  // 77 Ir
        {76, -1, -1, -1, -1, -1},  // 78 Pt
        {77, -1, -1, -1, -1, -1},  // 79 Au
        {78, -1, -1, -1, -1, -1},  // 80 Hg
        {79, -1, -1, -1, -1, -1},  // 81 Tl
        {80, -1, -1, -1, -1, -1},  // 82 Pb
        {60, 81, -1, -1, -1, -1},  // 83 Bi
        {82, -1, -1, -1, -1, -1},  // 84 Po
        {83, -1, -1, -1, -1, -1},  // 85 At
        {66, 84, -1, -1, -1, -1},  // 86 Rn
        {85, -1, -1, -1, -1, -1},  // 87 Fr
        {86, -1, -1, -1, -1, -1},  // 88 Ra
        {87, -1, -1, -1, -1, -1},  // 89 Ac
        {88, -1, -1, -1, -1, -1},  // 90 Th
        {89, -1, -1, -1, -1, -1},  // 91 Pa
        {90, -1, -1, -1, -1, -1}  // 92 U
    }
};

const std::array<int, K> kDecayCount = {
    1,  //  1 H
    5,  //  2 He
    1,  //  3 Li
    3,  //  4 Be
    1,  //  5 B
    1,  //  6 C
    1,  //  7 N
    1,  //  8 O
    1,  //  9 F
    1,  // 10 Ne
    1,  // 11 Na
    2,  // 12 Mg
    1,  // 13 Al
    1,  // 14 Si
    2,  // 15 P
    1,  // 16 S
    1,  // 17 Cl
    1,  // 18 Ar
    1,  // 19 K
    1,  // 20 Ca
    5,  // 21 Sc
    1,  // 22 Ti
    1,  // 23 V
    1,  // 24 Cr
    2,  // 25 Mn
    1,  // 26 Fe
    1,  // 27 Co
    2,  // 28 Ni
    1,  // 29 Cu
    1,  // 30 Zn
    6,  // 31 Ga
    2,  // 32 Ge
    2,  // 33 As
    1,  // 34 Se
    1,  // 35 Br
    1,  // 36 Kr
    1,  // 37 Rb
    1,  // 38 Sr
    2,  // 39 Y
    4,  // 40 Zr
    2,  // 41 Nb
    1,  // 42 Mo
    1,  // 43 Tc
    3,  // 44 Ru
    2,  // 45 Rh
    1,  // 46 Pd
    1,  // 47 Ag
    1,  // 48 Cd
    1,  // 49 In
    1,  // 50 Sn
    2,  // 51 Sb
    3,  // 52 Te
    2,  // 53 I
    1,  // 54 Xe
    1,  // 55 Cs
    1,  // 56 Ba
    1,  // 57 La
    4,  // 58 Ce
    1,  // 59 Pr
    1,  // 60 Nd
    1,  // 61 Pm
    3,  // 62 Sm
    1,  // 63 Eu
    3,  // 64 Gd
    2,  // 65 Tb
    1,  // 66 Dy
    1,  // 67 Ho
    2,  // 68 Er
    3,  // 69 Tm
    1,  // 70 Yb
    1,  // 71 Lu
    1,  // 72 Hf
    5,  // 73 Ta
    1,  // 74 W
    3,  // 75 Re
    1,  // 76 Os
    1,  // 77 Ir
    1,  // 78 Pt
    1,  // 79 Au
    1,  // 80 Hg
    1,  // 81 Tl
    1,  // 82 Pb
    2,  // 83 Bi
    1,  // 84 Po
    1,  // 85 At
    2,  // 86 Rn
    1,  // 87 Fr
    1,  // 88 Ra
    1,  // 89 Ac
    1,  // 90 Th
    1,  // 91 Pa
    1  // 92 U
};

using Vector = std::array<u64, K>;
using Matrix = std::array<std::array<u64, K>, K>;

Matrix build_transition_matrix() {
    Matrix m{};
    for (int src = 0; src < K; ++src) {
        for (int j = 0; j < kDecayCount[src]; ++j) {
            const int dst = kDecay[src][j];
            m[static_cast<std::size_t>(dst)][static_cast<std::size_t>(src)] += 1ULL;
        }
    }
    return m;
}

std::array<std::array<u64, 3>, K> build_digit_counts() {
    std::array<std::array<u64, 3>, K> counts{};
    for (int i = 0; i < K; ++i) {
        for (const char ch : std::string(kElementStrings[i])) {
            if (ch == '1') {
                counts[i][0] += 1ULL;
            } else if (ch == '2') {
                counts[i][1] += 1ULL;
            } else {
                counts[i][2] += 1ULL;
            }
        }
    }
    return counts;
}

std::string look_and_say_step(const std::string& s) {
    std::string out;
    out.reserve(s.size() * 2);
    for (std::size_t i = 0; i < s.size();) {
        std::size_t j = i + 1;
        while (j < s.size() && s[j] == s[i]) {
            ++j;
        }
        out.push_back(static_cast<char>('0' + (j - i)));
        out.push_back(s[i]);
        i = j;
    }
    return out;
}

std::array<u64, 3> direct_counts(u64 n) {
    std::string s = "1";
    for (u64 step = 1; step < n; ++step) {
        s = look_and_say_step(s);
    }
    std::array<u64, 3> cnt{0ULL, 0ULL, 0ULL};
    for (const char ch : s) {
        if (ch == '1') {
            ++cnt[0];
        } else if (ch == '2') {
            ++cnt[1];
        } else {
            ++cnt[2];
        }
    }
    return cnt;
}

Matrix multiply_matrix(const Matrix& a, const Matrix& b) {
    Matrix c{};
    for (int i = 0; i < K; ++i) {
        for (int k = 0; k < K; ++k) {
            if (a[i][k] == 0ULL) {
                continue;
            }
            const u64 aik = a[i][k];
            for (int j = 0; j < K; ++j) {
                if (b[k][j] == 0ULL) {
                    continue;
                }
                const u128 cur = static_cast<u128>(c[i][j]) +
                                 static_cast<u128>(aik) * static_cast<u128>(b[k][j]);
                c[i][j] = static_cast<u64>(cur & MASK);
            }
        }
    }
    return c;
}

Vector apply_matrix(const Matrix& m, const Vector& v) {
    Vector out{};
    for (int i = 0; i < K; ++i) {
        u128 acc = 0;
        for (int j = 0; j < K; ++j) {
            if (m[i][j] == 0ULL || v[j] == 0ULL) {
                continue;
            }
            acc += static_cast<u128>(m[i][j]) * static_cast<u128>(v[j]);
        }
        out[i] = static_cast<u64>(acc & MASK);
    }
    return out;
}

Vector evolve_vector(Vector v, u64 steps) {
    Matrix power = build_transition_matrix();
    u64 e = steps;
    while (e > 0ULL) {
        if (e & 1ULL) {
            v = apply_matrix(power, v);
        }
        e >>= 1ULL;
        if (e > 0ULL) {
            power = multiply_matrix(power, power);
        }
    }
    return v;
}

std::array<u64, 3> counts_from_vector(const Vector& v,
                                      const std::array<std::array<u64, 3>, K>& digit_counts) {
    std::array<u64, 3> out{0ULL, 0ULL, 0ULL};
    for (int i = 0; i < K; ++i) {
        if (v[i] == 0ULL) {
            continue;
        }
        for (int d = 0; d < 3; ++d) {
            const u128 add = static_cast<u128>(out[d]) +
                             static_cast<u128>(v[i]) * static_cast<u128>(digit_counts[i][d]);
            out[d] = static_cast<u64>(add & MASK);
        }
    }
    return out;
}

std::array<u64, 3> solve(u64 n) {
    if (n < 8ULL) {
        return direct_counts(n);
    }

    // The 8th term (1113213211) decomposes as Hf + Sn in Conway's element table.
    Vector v{};
    v[71] = 1ULL;  // Hf (72)
    v[49] = 1ULL;  // Sn (50)

    v = evolve_vector(v, n - 8ULL);
    const auto digit_counts = build_digit_counts();
    return counts_from_vector(v, digit_counts);
}

bool run_checkpoints() {
    const auto c8 = solve(8);
    if (!(c8[0] == 6ULL && c8[1] == 2ULL && c8[2] == 2ULL)) {
        std::cerr << "Checkpoint failed: term 8 counts\n";
        return false;
    }

    const auto c40 = solve(40);
    if (!(c40[0] == 31254ULL && c40[1] == 20259ULL && c40[2] == 11625ULL)) {
        std::cerr << "Checkpoint failed: term 40 counts\n";
        return false;
    }

    const auto c20_direct = direct_counts(20);
    const auto c20_fast = solve(20);
    if (c20_direct != c20_fast) {
        std::cerr << "Checkpoint failed: direct vs fast at term 20\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const auto ans = solve(options.n);
    std::cout << ans[0] << "," << ans[1] << "," << ans[2] << "\n";
    return 0;
}
