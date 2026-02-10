#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string s;
    while (value > 0) {
        int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

// For a valid cut, with integer piece areas (a,b,c,d) and b<=c, one can solve for the
// (generally rational) areas X,Y of triangles BCP and CAP in terms of (a,b,c). The
// total area is:
//   T = a*(a+b)*(a+c) / (a^2 - b*c)
// and the cut exists with integer (a,b,c,d) iff (a^2 - b*c) > 0 and T is an integer
// (then d = T-a-b-c is automatically a positive integer).
static u128 S(int n) {
    u128 total = 0;

    for (int a = 1; a <= n; ++a) {
        const int a_room = n - a;
        if (a_room < 3) continue;  // need b,c,d >= 1

        // Since b<=c and d>=1: a + 2b + 1 <= n.
        int bmax = std::min(a - 1, (a_room - 1) / 2);
        if (bmax < 1) continue;

        // T is increasing in c. The minimum for given (a,b) happens at c=b:
        //   T_min = a*(a+b)/(a-b).
        // If T_min > n then no c works, which gives an additional b upper bound.
        const u64 bound_num = static_cast<u64>(a) * static_cast<u64>(n - a);
        const u64 bound_den = static_cast<u64>(a + n);
        const int bmax_by_Tmin = static_cast<int>(bound_num / bound_den);
        bmax = std::min(bmax, bmax_by_Tmin);
        if (bmax < 1) continue;

        const u64 a2 = static_cast<u64>(a) * static_cast<u64>(a);

        for (int b = 1; b <= bmax; ++b) {
            const int rem = n - a - b;
            if (rem <= 2) continue;  // need c>=b and d>=1

            // From T<=n:
            //   c * (a(a+b) + n b) <= a^2 (n-a-b)
            const u64 rhs = a2 * static_cast<u64>(rem);
            const u64 lhs_coeff = static_cast<u64>(a) * static_cast<u64>(a + b) + static_cast<u64>(n) * static_cast<u64>(b);
            const u64 cmax_T = rhs / lhs_coeff;

            // From bc < a^2 and c integer:
            const u64 cmax_bc = (a2 - 1) / static_cast<u64>(b);

            // From a+b+c+d<=n and d>=1:
            const u64 cmax_sum = static_cast<u64>(n - a - b - 1);

            u64 cmax = std::min({cmax_T, cmax_bc, cmax_sum});
            if (cmax < static_cast<u64>(b)) continue;

            const u64 base = static_cast<u64>(a) * static_cast<u64>(a + b);  // fits in u64 for n<=10000
            for (int c = b; c <= static_cast<int>(cmax); ++c) {
                const u64 bc = static_cast<u64>(b) * static_cast<u64>(c);
                const u64 denom = a2 - bc;  // >0 by construction

                const u64 numer = base * static_cast<u64>(a + c);
                if (numer % denom != 0) continue;

                const u64 T = numer / denom;
                if (T > static_cast<u64>(n)) continue;

                const int d = static_cast<int>(T) - a - b - c;
                if (d <= 0) continue;

                total += T;
            }
        }
    }

    return total;
}

static std::vector<std::array<int, 4>> quads_with_total_area(int target_total) {
    std::vector<std::array<int, 4>> quads;
    for (int a = 1; a <= target_total - 3; ++a) {
        for (int b = 1; b <= std::min(a - 1, target_total - a - 2); ++b) {
            for (int c = b; c <= target_total - a - b - 1; ++c) {
                if (1LL * b * c >= 1LL * a * a) break;
                const u64 a2 = static_cast<u64>(a) * static_cast<u64>(a);
                const u64 denom = a2 - static_cast<u64>(b) * static_cast<u64>(c);
                const u128 numer = static_cast<u128>(a) * static_cast<u128>(a + b) * static_cast<u128>(a + c);
                if (numer % denom != 0) continue;
                const u64 T = static_cast<u64>(numer / denom);
                if (static_cast<int>(T) != target_total) continue;
                const int d = target_total - a - b - c;
                if (d <= 0) continue;
                quads.push_back({a, b, c, d});
            }
        }
    }
    std::sort(quads.begin(), quads.end());
    return quads;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    assert(to_string_u128(S(20)) == "259");

    {
        const auto quads = quads_with_total_area(55);
        const std::vector<std::array<int, 4>> expected = {
            {20, 2, 24, 9},
            {22, 8, 11, 14},
        };
        assert(quads == expected);
    }

    std::cout << to_string_u128(S(10000)) << '\n';
    return 0;
}
