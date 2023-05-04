#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// изменить на long double

double f1(double x) { return (double)1 + 4 / (x * x + 1); }

double f1der1(double x) { return (double)-8 * x / ((x * x + 1) * (x * x + 1)); }

double f1der2(double x) {
  return (double)8 * (3 * x * x - 1) /
         ((x * x + 1) * (x * x + 1) * (x * x + 1));
}

double f2(double x) { return (double)x * x * x; }

double f2der1(double x) { return (double)3 * x * x; }

double f2der2(double x) { return (double)6 * x; }

double f3(double x) { return (double)1 / std::pow(2, x); }

double f3der1(double x) { return (double)-1 / std::pow(2, x) * log(2); }

double f3der2(double x) { return (double)1 / std::pow(2, x) * log(2) * log(2); }

double root(double a, double b, double (*f)(double), double (*fder1)(double),
            double (*fder2)(double), double (*g)(double),
            double (*gder1)(double), double (*gder2)(double), double eps) {
  double f_a = f(a);
  double f_b = f(b);
  double g_a = g(a);
  double g_b = g(b);
  double c = (double)(b * fabs(g_a - f_a) + a * fabs(g_b - f_b)) /
             (fabs(g_b - f_b) + fabs(g_a - f_a));
  double f1 = fder1(c);
  double g1 = gder1(c);

  double g_c = g(c);
  double f_c = f(c);
  // проверить корректность
  if (g_c > f_c && g1 > 0 && f1 > 0 && f1 > g1) goto AC;
  if (g_c > f_c && g1 > 0 && f1 < 0) goto BC;
  if (g_c > f_c && g1 < 0 && f1 > 0) goto AC;
  if (g_c > f_c && g1 < 0 && f1 < 0 && g1 > f1) goto BC;
  if (g_c < f_c && g1 > 0 && f1 > 0 && g1 < f1) goto BC;
  if (g_c < f_c && g1 > 0 && f1 < 0) goto AC;
  if (g_c < f_c && g1 < 0 && f1 > 0) goto BC;
  if (g_c < f_c && g1 < 0 && f1 < 0 && f1 < g1) goto AC;

AC:
  if ((g_c - f_c) * (g(c + eps) - f(c + eps)) < 0)
    return c;
  else
    return root(c, b, f, fder1, fder2, g, gder1, gder2, eps);
BC:
  if ((g_c - f_c) * (g(c - eps) - f(c - eps)) < 0)
    return c;
  else
    return root(a, c, f, fder1, fder2, g, gder1, gder2, eps);
}

void getRoots(std::vector<std::pair<double, double>>& range,
              std::vector<double>& roots, double eps) {
  range[0].first = -2.0;
  range[0].second = -1.0;
  range[1].first = 0.0;
  range[1].second = 1.0;
  range[2].first = 1.0;
  range[2].second = 2.0;

  roots[0] = root(range[0].first, range[0].second, f1, f1der1, f1der2, f3,
                  f3der1, f3der2, eps);
  roots[1] = root(range[1].first, range[1].second, f2, f2der1, f2der2, f3,
                  f3der1, f3der2, eps);
  roots[2] = root(range[2].first, range[2].second, f1, f1der1, f1der2, f2,
                  f2der1, f2der2, eps);
  return;
}

double square(double a, double b, int n, double (*f)(double)) {
  double h = (double)(b - a) / n;
  double sum = 0.0;
  double c = a, d = a + h;
  while (d <= b) {
    sum += (f(c) + f(d)) / 2.0;
    c = d;
    d += h;
  }
  return h * sum;
}

double getSquare(double a, double b, double eps, double (*f)(double)) {
  double p = 3.0;
  double s = 0.0, s2 = 0.0;
  int n = 1000;
  do {
    s = square(a, b, n, f);
    n *= 2;
    s2 = square(a, b, n, f);
  } while (!(fabs(s - s2) / p <= eps));
  return s;
}

double getAnswer(std::vector<double>& roots, double eps) {
  double s = getSquare(roots[0], roots[2], eps, f1);
  double s1 = getSquare(roots[0], roots[1], eps, f3);
  double s2 = getSquare(roots[1], roots[2], eps, f2);
  return s - s1 - s2;
}

int main() {
  double eps;
  std::cin >> eps;
  std::vector<std::pair<double, double>> range(3);
  std::vector<double> roots(3);

  double eps1 = eps / 3.0;
  double eps2 = eps * 1E-3;
  getRoots(range, roots, eps2);

  // for (int i = 0; i < 3; ++i) {
  //   std::cout << roots[i] << " ";
  // }
  std::cout << std::fixed << std::setprecision(10) << getAnswer(roots, eps1);

  return 0;
}

// 6.5911043047