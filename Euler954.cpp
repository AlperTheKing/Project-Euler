#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u8 = std::uint8_t;

struct Option {
    u8 mask_all;
    u8 mask_nonzero;
    u8 sum_mod;
    u32 ways;
};

int norm7i(int x) {
    x %= 7;
    if (x < 0) {
        x += 7;
    }
    return x;
}

class HeptaphobiaCounter {
public:
    HeptaphobiaCounter() {
        build_rotations();
        build_first_ways();
        build_options();
    }

    u64 count_pow10(int exponent) const {
        u64 total = 0;
        for (int len = 1; len <= exponent; ++len) {
            total += count_length(len);
        }
        return total;
    }

private:
    static constexpr std::array<int, 6> kPow10Mod7 = {1, 3, 2, 6, 4, 5};
    static constexpr std::array<int, 7> kInvMod7 = {0, 1, 4, 5, 2, 3, 6};

    std::array<std::array<u8, 7>, 128> rotate_{};
    std::array<u32, 7> first_ways_{};
    std::array<std::vector<Option>, 4> options_by_count_;

    static int norm7(int x) {
        x %= 7;
        if (x < 0) {
            x += 7;
        }
        return x;
    }

    static int weight_pow10_mod7(int exp) {
        return kPow10Mod7[static_cast<std::size_t>(exp % 6)];
    }

    void build_rotations() {
        for (int mask = 0; mask < 128; ++mask) {
            for (int sh = 0; sh < 7; ++sh) {
                int out = 0;
                for (int r = 0; r < 7; ++r) {
                    if (((mask >> r) & 1) != 0) {
                        out |= (1 << ((r + sh) % 7));
                    }
                }
                rotate_[static_cast<std::size_t>(mask)][static_cast<std::size_t>(sh)] = static_cast<u8>(out);
            }
        }
    }

    void build_first_ways() {
        first_ways_.fill(0);
        for (int d = 1; d <= 9; ++d) {
            ++first_ways_[static_cast<std::size_t>(d % 7)];
        }
    }

    void add_option_count(std::unordered_map<u32, u32>& cnt, u8 mask_all, u8 mask_nonzero, u8 sum_mod) {
        const u32 key = static_cast<u32>(mask_all) |
                        (static_cast<u32>(mask_nonzero) << 7U) |
                        (static_cast<u32>(sum_mod) << 14U);
        ++cnt[key];
    }

    std::vector<Option> make_options_for_count(int cnt_digits) {
        std::unordered_map<u32, u32> cnt;
        if (cnt_digits == 0) {
            add_option_count(cnt, 0U, 0U, 0U);
        } else if (cnt_digits == 1) {
            for (int d0 = 0; d0 <= 9; ++d0) {
                const int r0 = d0 % 7;
                const u8 m_all = static_cast<u8>(1U << r0);
                const u8 m_nz = static_cast<u8>(d0 == 0 ? 0U : (1U << r0));
                add_option_count(cnt, m_all, m_nz, static_cast<u8>(r0));
            }
        } else if (cnt_digits == 2) {
            for (int d0 = 0; d0 <= 9; ++d0) {
                for (int d1 = 0; d1 <= 9; ++d1) {
                    const int r0 = d0 % 7;
                    const int r1 = d1 % 7;
                    const u8 m_all = static_cast<u8>((1U << r0) | (1U << r1));
                    u8 m_nz = 0U;
                    if (d0 != 0) {
                        m_nz = static_cast<u8>(m_nz | (1U << r0));
                    }
                    if (d1 != 0) {
                        m_nz = static_cast<u8>(m_nz | (1U << r1));
                    }
                    add_option_count(cnt, m_all, m_nz, static_cast<u8>((r0 + r1) % 7));
                }
            }
        } else {
            for (int d0 = 0; d0 <= 9; ++d0) {
                for (int d1 = 0; d1 <= 9; ++d1) {
                    for (int d2 = 0; d2 <= 9; ++d2) {
                        const int r0 = d0 % 7;
                        const int r1 = d1 % 7;
                        const int r2 = d2 % 7;
                        const u8 m_all = static_cast<u8>((1U << r0) | (1U << r1) | (1U << r2));
                        u8 m_nz = 0U;
                        if (d0 != 0) {
                            m_nz = static_cast<u8>(m_nz | (1U << r0));
                        }
                        if (d1 != 0) {
                            m_nz = static_cast<u8>(m_nz | (1U << r1));
                        }
                        if (d2 != 0) {
                            m_nz = static_cast<u8>(m_nz | (1U << r2));
                        }
                        add_option_count(cnt, m_all, m_nz, static_cast<u8>((r0 + r1 + r2) % 7));
                    }
                }
            }
        }

        std::vector<Option> options;
        options.reserve(cnt.size());
        for (const auto& [key, ways] : cnt) {
            const u8 m_all = static_cast<u8>(key & 127U);
            const u8 m_nz = static_cast<u8>((key >> 7U) & 127U);
            const u8 sum_mod = static_cast<u8>((key >> 14U) & 7U);
            options.push_back(Option{m_all, m_nz, sum_mod, ways});
        }
        return options;
    }

    void build_options() {
        for (int c = 0; c <= 3; ++c) {
            options_by_count_[static_cast<std::size_t>(c)] = make_options_for_count(c);
        }
    }

    u64 count_length(int len) const {
        u64 total = 0;
        for (int m = 1; m <= 6; ++m) {
            total += count_length_for_mod(len, m);
        }
        return total;
    }

    u64 count_length_for_mod(int len, int mod_target) const {
        std::vector<int> weights(static_cast<std::size_t>(len), 0);
        for (int i = 0; i < len; ++i) {
            const int exp = len - 1 - i;
            weights[static_cast<std::size_t>(i)] = weight_pow10_mod7(exp);
        }
        const int first_weight = weights[0];

        std::array<int, 7> counts{};
        counts.fill(0);
        for (int i = 1; i < len; ++i) {
            ++counts[static_cast<std::size_t>(weights[static_cast<std::size_t>(i)])];
        }

        struct VarData {
            int weight{0};
            const std::vector<Option>* options{nullptr};
            std::vector<int> contrib_mod;
            std::vector<u32> ways;
            u64 full_mask{0};
        };

        std::vector<VarData> vars;
        vars.reserve(6);
        for (int w : {1, 2, 3, 4, 5, 6}) {
            const int cnt = counts[static_cast<std::size_t>(w)];
            if (cnt == 0) {
                continue;
            }
            VarData v;
            v.weight = w;
            v.options = &options_by_count_[static_cast<std::size_t>(cnt)];
            v.contrib_mod.reserve(v.options->size());
            v.ways.reserve(v.options->size());
            for (const Option& opt : *v.options) {
                v.contrib_mod.push_back(norm7(w * static_cast<int>(opt.sum_mod)));
                v.ways.push_back(opt.ways);
            }
            v.full_mask = (v.options->size() == 64) ? ~0ULL : ((1ULL << v.options->size()) - 1ULL);
            vars.push_back(std::move(v));
        }

        const int k = static_cast<int>(vars.size());
        const int full_assigned = (1 << k) - 1;

        std::vector<std::vector<std::vector<u64>>> compat(
            static_cast<std::size_t>(k),
            std::vector<std::vector<u64>>(static_cast<std::size_t>(k)));

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                if (i == j) {
                    continue;
                }
                const int wi = vars[static_cast<std::size_t>(i)].weight;
                const int wj = vars[static_cast<std::size_t>(j)].weight;
                const int diff = norm7(wi - wj);
                const int c = norm7(-mod_target * kInvMod7[static_cast<std::size_t>(diff)]);

                const auto& oi = *vars[static_cast<std::size_t>(i)].options;
                const auto& oj = *vars[static_cast<std::size_t>(j)].options;
                compat[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)].assign(oi.size(), 0ULL);

                for (std::size_t a = 0; a < oi.size(); ++a) {
                    const u8 shifted = rotate_[static_cast<std::size_t>(oi[a].mask_all)][static_cast<std::size_t>(c)];
                    u64 mask = 0ULL;
                    for (std::size_t b = 0; b < oj.size(); ++b) {
                        if ((shifted & oj[b].mask_all) == 0U) {
                            mask |= (1ULL << b);
                        }
                    }
                    compat[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)][a] = mask;
                }
            }
        }

        u64 total = 0;
        for (int first_residue = 0; first_residue <= 6; ++first_residue) {
            const u32 ways_first = first_ways_[static_cast<std::size_t>(first_residue)];
            if (ways_first == 0U) {
                continue;
            }

            std::array<u64, 6> candidates{};
            candidates.fill(0ULL);
            bool ok = true;
            for (int i = 0; i < k; ++i) {
                candidates[static_cast<std::size_t>(i)] = vars[static_cast<std::size_t>(i)].full_mask;
            }

            for (int i = 0; i < k; ++i) {
                const int w = vars[static_cast<std::size_t>(i)].weight;
                if (w == first_weight) {
                    continue;
                }
                const int diff = norm7(first_weight - w);
                const int c = norm7(-mod_target * kInvMod7[static_cast<std::size_t>(diff)]);
                const int forbid = (first_residue + c) % 7;

                u64 keep = 0ULL;
                const auto& opts = *vars[static_cast<std::size_t>(i)].options;
                for (std::size_t oi = 0; oi < opts.size(); ++oi) {
                    if (((opts[oi].mask_nonzero >> forbid) & 1U) == 0U) {
                        keep |= (1ULL << oi);
                    }
                }
                candidates[static_cast<std::size_t>(i)] &= keep;
                if (candidates[static_cast<std::size_t>(i)] == 0ULL) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                continue;
            }

            const int target_rest = norm7(mod_target - first_residue * first_weight);

            auto dfs = [&](auto&& self,
                           int assigned_mask,
                           int cur_mod,
                           const std::array<u64, 6>& cands) -> u64 {
                if (assigned_mask == full_assigned) {
                    return static_cast<u64>(cur_mod == target_rest ? 1 : 0);
                }

                int pick = -1;
                int best = 1'000'000'000;
                for (int i = 0; i < k; ++i) {
                    if (((assigned_mask >> i) & 1) != 0) {
                        continue;
                    }
                    const int pc = std::popcount(cands[static_cast<std::size_t>(i)]);
                    if (pc == 0) {
                        return 0ULL;
                    }
                    if (pc < best) {
                        best = pc;
                        pick = i;
                    }
                }

                const auto& var = vars[static_cast<std::size_t>(pick)];
                u64 bits = cands[static_cast<std::size_t>(pick)];
                u64 subtotal = 0;

                while (bits != 0ULL) {
                    const u64 bit = bits & (~bits + 1ULL);
                    bits -= bit;
                    const int oi = std::countr_zero(bit);

                    std::array<u64, 6> next = cands;
                    const int next_assigned = assigned_mask | (1 << pick);
                    const int next_mod = (cur_mod + var.contrib_mod[static_cast<std::size_t>(oi)]) % 7;

                    bool good = true;
                    for (int j = 0; j < k; ++j) {
                        if (((next_assigned >> j) & 1) != 0) {
                            continue;
                        }
                        next[static_cast<std::size_t>(j)] &=
                            compat[static_cast<std::size_t>(pick)][static_cast<std::size_t>(j)][static_cast<std::size_t>(oi)];
                        if (next[static_cast<std::size_t>(j)] == 0ULL) {
                            good = false;
                            break;
                        }
                    }
                    if (!good) {
                        continue;
                    }

                    subtotal += static_cast<u64>(var.ways[static_cast<std::size_t>(oi)]) *
                                self(self, next_assigned, next_mod, next);
                }

                return subtotal;
            };

            total += static_cast<u64>(ways_first) * dfs(dfs, 0, 0, candidates);
        }

        return total;
    }
};

bool is_heptaphobic_bruteforce(u64 n) {
    if (n % 7ULL == 0ULL) {
        return false;
    }

    std::vector<int> digits;
    {
        u64 x = n;
        while (x > 0) {
            digits.push_back(static_cast<int>(x % 10ULL));
            x /= 10ULL;
        }
        std::reverse(digits.begin(), digits.end());
    }

    const int len = static_cast<int>(digits.size());
    std::vector<int> weights(static_cast<std::size_t>(len), 0);
    constexpr std::array<int, 6> kPow10Mod7 = {1, 3, 2, 6, 4, 5};
    for (int i = 0; i < len; ++i) {
        const int exp = len - 1 - i;
        weights[static_cast<std::size_t>(i)] = kPow10Mod7[static_cast<std::size_t>(exp % 6)];
    }

    const int mod = static_cast<int>(n % 7ULL);
    for (int i = 0; i < len; ++i) {
        for (int j = i + 1; j < len; ++j) {
            if (i == 0 && digits[static_cast<std::size_t>(j)] == 0) {
                continue;
            }
            const int delta_digit = digits[static_cast<std::size_t>(j)] - digits[static_cast<std::size_t>(i)];
            const int delta_weight = weights[static_cast<std::size_t>(i)] - weights[static_cast<std::size_t>(j)];
            const int delta = norm7i(delta_digit * delta_weight);
            if ((mod + delta) % 7 == 0) {
                return false;
            }
        }
    }

    return true;
}

u64 brute_count_pow10(int exponent) {
    u64 limit = 1ULL;
    for (int i = 0; i < exponent; ++i) {
        limit *= 10ULL;
    }
    u64 count = 0;
    for (u64 n = 1; n < limit; ++n) {
        if (is_heptaphobic_bruteforce(n)) {
            ++count;
        }
    }
    return count;
}

void run_validations() {
    HeptaphobiaCounter counter;
    assert(counter.count_pow10(2) == 74ULL);
    assert(counter.count_pow10(4) == 3'737ULL);
    assert(counter.count_pow10(3) == brute_count_pow10(3));
}

}  // namespace

int main() {
    run_validations();
    HeptaphobiaCounter counter;
    std::cout << counter.count_pow10(13) << '\n';
    return 0;
}
