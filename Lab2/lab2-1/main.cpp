#include <iostream>
#include <cmath>

using namespace std;

int iter_count = 0;

double f(double x) {
    return x * x * x * x - 2.0 * x - 1.0;
}

double phi(double x) {
    return sqrt(sqrt(2 * x + 1));
}

double phi_s(double x) {
    return 1 / (2 * sqrt(sqrt((2 * x + 1) * (2 * x + 1) * (2 * x + 1))));
}

double f_s(double x) {
    return 4.0 * x * x * x - 2.0;
}

double f_ss(double x) {
    return 12.0 * x * x;
}

double iter_solve(double l, double r, double eps) {
    iter_count = 0;
    double x_k = r;
    double dx = 1.0;
    double q = max(abs(phi_s(l)), abs(phi_s(r)));
    double eps_coef = q / (1.0 - q);
    do {
        double x_k1 = phi(x_k);
        dx = eps_coef * abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

double newton_solve(double l, double r, double eps) {
    double x0 = l;
    if (!(f(x0) * f_ss(x0) > eps)) {
        x0 = r;
    }
    iter_count = 0;
    double x_k = x0;
    double dx = 1.0;
    do {
        double x_k1 = x_k - f(x_k) / f_s(x_k);
        dx = abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

int main() {
    cout.precision(9);
    double l, r, eps;
    cin >> l >> r >> eps;
    double root;
    root = iter_solve(l, r, eps);
    cout << "Метод простой итерации" << endl;
    cout << "x0 = " << root << endl;
    cout << "Решение получено за " << iter_count << " итераций" << "\n\n";
    root = newton_solve(l, r, eps);
    cout << "Метод Ньютона" << endl;
    cout << "x0 = " << root << endl;
    cout << "Решение получена за " << iter_count << " итераций" << "\n\n";
}
