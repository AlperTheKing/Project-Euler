#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

class Solver {
public:
    Solver() {
        const std::string phrase = "thereisasyetinsufficientdataforameaningfulanswer";
        std::array<int, 26> full_counts{};
        for (const char ch : phrase) {
            ++full_counts[static_cast<std::size_t>(ch - 'a')];
        }

        char_to_idx_.fill(-1);
        for (int ch = 0; ch < 26; ++ch) {
            if (full_counts[static_cast<std::size_t>(ch)] > 0) {
                char_to_idx_[static_cast<std::size_t>(ch)] = static_cast<int>(init_counts_.size());
                init_counts_.push_back(static_cast<std::uint8_t>(full_counts[static_cast<std::size_t>(ch)]));
                chars_.push_back(ch);
            }
        }

        radix_.assign(init_counts_.size(), 0ULL);
        if (!radix_.empty()) {
            radix_[0] = 1ULL;
            for (std::size_t i = 1; i < radix_.size(); ++i) {
                radix_[i] =
                    radix_[i - 1] * static_cast<u64>(init_counts_[i - 1] + 1U);
            }
        }
        for (auto& mp : memo_) {
            mp.reserve(600'000);
        }
    }

    u64 position_of(const std::string& word) {
        if (word.empty() || word.size() > 15U) {
            throw std::runtime_error("Invalid word length for ranking.");
        }

        std::vector<std::uint8_t> counts = init_counts_;
        u64 pos = 0ULL;

        for (std::size_t i = 0; i < word.size(); ++i) {
            const int target = word[i] - 'a';
            if (target < 0 || target >= 26) {
                throw std::runtime_error("Invalid character in word.");
            }

            for (int ch = 0; ch < target; ++ch) {
                const int idx = char_to_idx_[static_cast<std::size_t>(ch)];
                if (idx < 0 || counts[static_cast<std::size_t>(idx)] == 0U) {
                    continue;
                }
                --counts[static_cast<std::size_t>(idx)];
                pos += count_subtree(static_cast<int>(15 - (i + 1)), counts);
                ++counts[static_cast<std::size_t>(idx)];
            }

            const int idx = char_to_idx_[static_cast<std::size_t>(target)];
            if (idx < 0 || counts[static_cast<std::size_t>(idx)] == 0U) {
                throw std::runtime_error("Word cannot be formed from phrase letters.");
            }
            --counts[static_cast<std::size_t>(idx)];
            ++pos;  // current prefix itself
        }

        return pos;
    }

    std::string word_at(u64 p) {
        if (p == 0ULL) {
            throw std::runtime_error("Positions are 1-based.");
        }

        std::vector<std::uint8_t> counts = init_counts_;
        std::string out;
        out.reserve(15U);
        int rem = 14;

        while (true) {
            bool chosen = false;
            for (int ch = 0; ch < 26; ++ch) {
                const int idx = char_to_idx_[static_cast<std::size_t>(ch)];
                if (idx < 0 || counts[static_cast<std::size_t>(idx)] == 0U) {
                    continue;
                }
                --counts[static_cast<std::size_t>(idx)];
                const u64 subtree = count_subtree(rem, counts);
                if (p > subtree) {
                    p -= subtree;
                    ++counts[static_cast<std::size_t>(idx)];
                    continue;
                }

                out.push_back(static_cast<char>('a' + ch));
                chosen = true;
                if (p == 1ULL) {
                    return out;
                }
                --p;
                --rem;
                break;
            }
            if (!chosen) {
                throw std::runtime_error("Requested position is out of range.");
            }
        }
    }

private:
    std::array<int, 26> char_to_idx_{};
    std::vector<int> chars_;
    std::vector<std::uint8_t> init_counts_;
    std::vector<u64> radix_;
    std::array<std::unordered_map<u32, u64>, 15> memo_;

    u32 encode_id(const std::vector<std::uint8_t>& counts) const {
        u64 key = 0ULL;
        for (std::size_t i = 0; i < counts.size(); ++i) {
            key += static_cast<u64>(counts[i]) * radix_[i];
        }
        return static_cast<u32>(key);
    }

    u64 count_subtree(const int rem, std::vector<std::uint8_t>& counts) {
        if (rem == 0) {
            return 1ULL;
        }
        const u32 id = encode_id(counts);
        auto& bucket = memo_[static_cast<std::size_t>(rem)];
        const auto it = bucket.find(id);
        if (it != bucket.end()) {
            return it->second;
        }

        u64 total = 1ULL;  // current word (empty extension)
        for (std::size_t i = 0; i < counts.size(); ++i) {
            if (counts[i] == 0U) {
                continue;
            }
            --counts[i];
            total += count_subtree(rem - 1, counts);
            ++counts[i];
        }
        bucket.emplace(id, total);
        return total;
    }
};

bool run_checkpoints(Solver& solver) {
    if (solver.word_at(10ULL) != "aaaaaacdee") {
        std::cerr << "Checkpoint failed: W(10)\n";
        return false;
    }
    if (solver.position_of("aaaaaacdee") != 10ULL) {
        std::cerr << "Checkpoint failed: P(aaaaaacdee)\n";
        return false;
    }
    if (solver.word_at(115246685191495243ULL) != "euler") {
        std::cerr << "Checkpoint failed: W(115246685191495243)\n";
        return false;
    }
    if (solver.position_of("euler") != 115246685191495243ULL) {
        std::cerr << "Checkpoint failed: P(euler)\n";
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

    Solver solver;
    if (options.run_checkpoints && !run_checkpoints(solver)) {
        return 1;
    }

    const i128 p_legionary = static_cast<i128>(solver.position_of("legionary"));
    const i128 p_calorimeters = static_cast<i128>(solver.position_of("calorimeters"));
    const i128 p_annihilate = static_cast<i128>(solver.position_of("annihilate"));
    const i128 p_orchestrated = static_cast<i128>(solver.position_of("orchestrated"));
    const i128 p_fluttering = static_cast<i128>(solver.position_of("fluttering"));

    const i128 target =
        p_legionary + p_calorimeters - p_annihilate + p_orchestrated - p_fluttering;
    if (target <= 0) {
        std::cerr << "Computed position is non-positive.\n";
        return 1;
    }

    const std::string answer = solver.word_at(static_cast<u64>(target));
    std::cout << answer << '\n';
    return 0;
}
