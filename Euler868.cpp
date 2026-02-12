#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;

static u64 rank_sjt(const std::vector<int>& p) {
    int n = static_cast<int>(p.size());
    if (n <= 1) return 0;

    int mx = n - 1;
    int pos = 0;
    while (p[pos] != mx) ++pos;

    std::vector<int> q;
    q.reserve(n - 1);
    for (int i = 0; i < n; ++i) {
        if (i != pos) q.push_back(p[i]);
    }

    u64 b = rank_sjt(q);
    u64 off = (b & 1ULL) ? static_cast<u64>(pos) : static_cast<u64>(n - 1 - pos);
    return b * static_cast<u64>(n) + off;
}

static u64 swaps_to_reach(const std::string& word) {
    std::string sorted = word;
    std::sort(sorted.begin(), sorted.end());

    std::array<int, 256> idx{};
    idx.fill(-1);
    for (int i = 0; i < static_cast<int>(sorted.size()); ++i) {
        unsigned char c = static_cast<unsigned char>(sorted[i]);
        assert(idx[c] == -1);
        idx[c] = i;
    }

    std::vector<int> perm;
    perm.reserve(word.size());
    for (unsigned char c : word) {
        assert(idx[c] != -1);
        perm.push_back(idx[c]);
    }

    return rank_sjt(perm);
}

int main() {
    assert(swaps_to_reach("CBA") == 3ULL);
    assert(swaps_to_reach("BELFRY") == 59ULL);

    std::cout << swaps_to_reach("NOWPICKBELFRYMATHS") << '\n';
    return 0;
}
