#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace {
constexpr int GRID = 250;
constexpr int DIM = GRID + 1;
constexpr double kNeg = -1e100;
constexpr double kEps = 1e-12;

struct Vec {
  int u;
  int v;
  double len;
  double slope;
};

int idx(int x, int y) { return x * DIM + y; }

std::vector<Vec> build_vectors() {
  std::vector<Vec> vecs;
  vecs.reserve(40000);

  vecs.push_back({1, 0, 1.0, 0.0});
  for (int u = 1; u <= GRID; ++u) {
    for (int v = 1; v <= GRID; ++v) {
      if (std::gcd(u, v) == 1) {
        double len = std::sqrt(static_cast<double>(u) * u +
                               static_cast<double>(v) * v);
        double slope = static_cast<double>(v) / u;
        vecs.push_back({u, v, len, slope});
      }
    }
  }
  vecs.push_back({0, 1, 1.0, 1e100});

  std::sort(vecs.begin(), vecs.end(),
            [](const Vec &a, const Vec &b) { return a.slope < b.slope; });
  return vecs;
}

struct Result {
  double area;
  double len;
  double score;
};

Result solve_for_lambda(double lambda, const std::vector<Vec> &vecs) {
  static std::vector<double> score(DIM * DIM);
  static std::vector<double> area(DIM * DIM);
  static std::vector<double> length(DIM * DIM);

  std::fill(score.begin(), score.end(), kNeg);
  std::fill(area.begin(), area.end(), 0.0);
  std::fill(length.begin(), length.end(), 0.0);

  score[idx(0, 0)] = 0.0;

  for (const auto &vec : vecs) {
    int u = vec.u;
    int v = vec.v;
    double len = vec.len;
    double base_area = static_cast<double>(u) * (GRID - 0.5 * v);
    double base_score = base_area - lambda * len;

    for (int i = 0; i <= GRID - u; ++i) {
      double *score_row = &score[idx(i, 0)];
      double *area_row = &area[idx(i, 0)];
      double *len_row = &length[idx(i, 0)];
      double *score_row_to = &score[idx(i + u, 0)];
      double *area_row_to = &area[idx(i + u, 0)];
      double *len_row_to = &length[idx(i + u, 0)];

      double area_gain = base_area;
      double score_gain = base_score;

      for (int j = 0; j <= GRID - v; ++j) {
        double cur_score = score_row[j];
        if (cur_score > kNeg / 2) {
          int tj = j + v;
          double cand_score = cur_score + score_gain;
          if (cand_score > score_row_to[tj] + kEps) {
            score_row_to[tj] = cand_score;
            area_row_to[tj] = area_row[j] + area_gain;
            len_row_to[tj] = len_row[j] + len;
          }
        }
        area_gain -= u;
        score_gain -= u;
      }
    }
  }

  int end = idx(GRID, GRID);
  return {area[end], length[end], score[end]};
}
} // namespace

int main() {
  auto vecs = build_vectors();

  double lambda = 132.0;
  for (int iter = 0; iter < 40; ++iter) {
    Result res = solve_for_lambda(lambda, vecs);
    double new_lambda = res.area / res.len;
    if (std::abs(new_lambda - lambda) < 1e-12 ||
        std::abs(res.score) < 1e-11) {
      lambda = new_lambda;
      break;
    }
    lambda = new_lambda;
  }

  std::cout << std::fixed << std::setprecision(8) << lambda << '\n';
  return 0;
}
