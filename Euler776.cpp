#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace {

using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::cpp_int;

struct Stats {
    cpp_int count = 0;
    cpp_int sum = 0;
};

cpp_dec_float_100 F_of_u64(std::uint64_t N) {
    const std::string s = std::to_string(N);
    const int L = static_cast<int>(s.size());
    const int MAX_SUM = 9 * L;

    std::vector<Stats> tight(MAX_SUM + 1), loose(MAX_SUM + 1);
    tight[0].count = 1;

    for (int pos = 0; pos < L; ++pos) {
        const int limit = s[static_cast<std::size_t>(pos)] - '0';
        std::vector<Stats> ntight(MAX_SUM + 1), nloose(MAX_SUM + 1);

        for (int sum = 0; sum <= MAX_SUM; ++sum) {
            if (tight[sum].count != 0) {
                for (int d = 0; d <= limit; ++d) {
                    const int ns = sum + d;
                    Stats& dst = (d == limit) ? ntight[ns] : nloose[ns];
                    dst.count += tight[sum].count;
                    dst.sum += tight[sum].sum * 10 + tight[sum].count * d;
                }
            }

            if (loose[sum].count != 0) {
                for (int d = 0; d <= 9; ++d) {
                    const int ns = sum + d;
                    Stats& dst = nloose[ns];
                    dst.count += loose[sum].count;
                    dst.sum += loose[sum].sum * 10 + loose[sum].count * d;
                }
            }
        }

        tight.swap(ntight);
        loose.swap(nloose);
    }

    cpp_dec_float_100 ans = 0;
    for (int sum = 1; sum <= MAX_SUM; ++sum) {
        const cpp_int total_sum = tight[sum].sum + loose[sum].sum;
        if (total_sum == 0) {
            continue;
        }
        ans += cpp_dec_float_100(total_sum) / cpp_dec_float_100(sum);
    }
    return ans;
}

void assert_close(const cpp_dec_float_100& got, const cpp_dec_float_100& expected,
                  const cpp_dec_float_100& tol) {
    cpp_dec_float_100 diff = got - expected;
    if (diff < 0) {
        diff = -diff;
    }
    assert(diff <= tol);
}

}  // namespace

int main() {
    assert(F_of_u64(10) == cpp_dec_float_100(19));
    assert_close(F_of_u64(123), cpp_dec_float_100("1.187764610390e3"),
                 cpp_dec_float_100("1e-9"));
    assert_close(F_of_u64(12345), cpp_dec_float_100("4.855801996238e6"),
                 cpp_dec_float_100("1e-6"));

    const cpp_dec_float_100 ans = F_of_u64(1'234'567'890'123'456'789ULL);
    std::cout << std::scientific << std::setprecision(12) << ans << '\n';
    return 0;
}
