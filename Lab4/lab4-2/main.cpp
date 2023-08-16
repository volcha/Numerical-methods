#include <iostream>
#include <cmath>
#include <exception>
#include <vector>
#include <tuple>
#include <functional>

using namespace std;

const double EPS = 0.000000001;
template<class T>
class tridiag_t {
private:
    int n;
    std::vector<T> a;
    std::vector<T> b;
    std::vector<T> c;
public:
    tridiag_t(const int & _n) : n(_n), a(n), b(n), c(n) {}
    tridiag_t(const std::vector<T> & _a, const std::vector<T> & _b, const std::vector<T> & _c) {
        if (!(_a.size() == _b.size() and _a.size() == _c.size())) {
            throw std::invalid_argument("Sizes of a, b, c are invalid");
        }
        n = _a.size();
        a = _a;
        b = _b;
        c = _c;
    }
    std::vector<T> solve(const std::vector<T> & d) {
        int m = d.size();
        if (n != m) {
            throw std::invalid_argument("Size of vector d is invalid");
        }
        std::vector<T> p(n);
        p[0] = -c[0] / b[0];
        std::vector<T> q(n);
        q[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
            q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
        }
        std::vector<T> x(n);
        x.back() = q.back();
        for (int i = n - 2; i >= 0; --i) {
            x[i] = p[i] * x[i + 1] + q[i];
        }
        return x;
    }
    friend std::istream & operator >> (std::istream & in, tridiag_t<T> & tridiag) {
        in >> tridiag.b[0] >> tridiag.c[0];
        for (int i = 1; i < tridiag.n - 1; ++i) {
            in >> tridiag.a[i] >> tridiag.b[i] >> tridiag.c[i];
        }
        in >> tridiag.a.back() >> tridiag.b.back();
        return in;
    }
    ~tridiag_t() = default;
};
void printData(const std::vector<std::tuple<double, double, double>> & v) {
    std::cout << "x = [";
    for (int i = 0; i < v.size(); ++i) {
        if (i) {
            std::cout << ", ";
        }
        std::cout << std::get<0>(v[i]);
    }
    std::cout << "]\n";
    std::cout << "y = [";
    for (int i = 0; i < v.size(); ++i) {
        if (i) {
            std::cout << ", ";
        }
        std::cout << std::get<1>(v[i]);
    }
    std::cout << "]\n";
}
bool border(double a, double b) {
    return (a < b) or (std::abs(b - a) < EPS);
}
class rungeKutta {
private:
    double l, r, y0, z0;
    std::function<double(double, double, double)> f, g;
public:
    rungeKutta(const double _l, const double _r,
        const std::function<double(double, double, double)> _f, const std::function<double(double, double, double)> _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}
    std::vector<std::tuple<double, double, double>> solve(double h) {
        double xi = l;
        double yi = y0;
        double zi = z0;
        std::vector<std::tuple<double, double, double>> res;
        res.push_back(std::make_tuple(xi, yi, zi));
        while (border(xi + h, r)) {
            double k1 = h * f(xi, yi, zi);
            double l1 = h * g(xi, yi, zi);
            double k2 = h * f(xi + 0.5 * h, yi + 0.5 * k1, zi + 0.5 * l1);
            double l2 = h * g(xi + 0.5 * h, yi + 0.5 * k1, zi + 0.5 * l1);
            double k3 = h * f(xi + 0.5 * h, yi + 0.5 * k2, zi + 0.5 * l2);
            double l3 = h * g(xi + 0.5 * h, yi + 0.5 * k2, zi + 0.5 * l2);
            double k4 = h * f(xi + h, yi + k3, zi + l3);
            double l4 = h * g(xi + h, yi + k3, zi + l3);
            double dy = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            double dz = (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;
            xi += h;
            yi += dy;
            zi += dz;
            res.push_back(std::make_tuple(xi, yi, zi));
        }
        return res;
    }
};
double rungeRomberg(const std::vector<std::tuple<double, double, double>> & y_2h, const std::vector<std::tuple<double, double, double>> & y_h, double p) {
    double coef = 1.0 / (std::pow(2, p) - 1.0);
    double res = 0.0;
    for (int i = 0; i < y_2h.size(); ++i) {
        res = std::max(res, coef * std::abs(std::get<1>(y_2h[i]) - std::get<1>(y_h[2 * i])));
    }
    return res;
}
class shooting {
private:
    double a, b;
    std::function<double(double, double, double)> f, g;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    shooting(const double _a, const double _b,
             const std::function<double(double, double, double)> _f,
             const std::function<double(double, double, double)> _g,
             const double _alpha, const double _beta, const double _y0,
             const double _delta, const double _gamma, const double _y1)
            : a(_a), b(_b), f(_f), g(_g),
              alpha(_alpha), beta(_beta), y0(_y0),
              delta(_delta), gamma(_gamma), y1(_y1) {}

    double get_start_cond(double eta) {
        return (y0 - alpha * eta) / beta;
    }

    double get_eta_next(double eta_prev, double eta, const std::vector<std::tuple<double, double, double>> sol_prev, const std::vector<std::tuple<double, double, double>> sol) {
        double yb_prev = std::get<1>(sol_prev.back());
        double zb_prev = std::get<2>(sol_prev.back());
        double phi_prev = delta * yb_prev + gamma * zb_prev - y1;
        double yb = std::get<1>(sol.back());
        double zb = std::get<2>(sol.back());
        double phi = delta * yb + gamma * zb - y1;
        return eta - (eta - eta_prev) / (phi - phi_prev) * phi;
    }
    std::vector<std::tuple<double, double, double>> solve(double h, double eps) {
        double eta_prev = 1.0;
        double eta = 0.8;
        while (1) {
            double rungeKutta_z0_prev = get_start_cond(eta_prev);
            rungeKutta de_solver_prev(a, b, f, g, eta_prev, rungeKutta_z0_prev);
            std::vector<std::tuple<double, double, double>> sol_prev = de_solver_prev.solve(h);

            double rungeKutta_z0 = get_start_cond(eta);
            rungeKutta de_solver(a, b, f, g, eta, rungeKutta_z0);
            std::vector<std::tuple<double, double, double>> sol = de_solver.solve(h);

            double eta_next = get_eta_next(eta_prev, eta, sol_prev, sol);
            if (std::abs(eta_next - eta) < eps) {
                return sol;
            } else {
                eta_prev = eta;
                eta = eta_next;
            }
        }
    }
};
class finDif {
private:
    using fx = std::function<double(double)>;
    using tridiag = tridiag_t<double>;

    double a, b;
    fx p, q, f;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    finDif(const double _a, const double _b,
        const fx _p, const fx _q, const fx _f,
        const double _alpha, const double _beta, const double _y0,
        const double _delta, const double _gamma, const double _y1)
        : a(_a), b(_b), p(_p), q(_q), f(_f),
        alpha(_alpha), beta(_beta), y0(_y0),
        delta(_delta), gamma(_gamma), y1(_y1) {}

    std::vector<std::tuple<double, double, double>> solve(double h) {
        int n = (b - a) / h;
        std::vector<double> xk(n + 1);
        for (int i = 0; i <= n; ++i) {
            xk[i] = a + h * i;
        }
        std::vector<double> a(n + 1);
        std::vector<double> b(n + 1);
        std::vector<double> c(n + 1);
        std::vector<double> d(n + 1);
        b[0] = h * alpha - beta;
        c[0] = beta;
        d[0] = h * y0;
        a.back() = -gamma;
        b.back() = h * delta + gamma;
        d.back() = h * y1;
        for (int i = 1; i < n; ++i) {
            a[i] = 1.0 - p(xk[i]) * h * 0.5;
            b[i] = -2.0 + h * h * q(xk[i]);
            c[i] = 1.0 + p(xk[i]) * h * 0.5;
            d[i] = h * h * f(xk[i]);
        }
        tridiag sys_eq(a, b, c);
        std::vector<double> yk = sys_eq.solve(d);
        std::vector<std::tuple<double, double, double>> res;
        for (int i = 0; i <= n; ++i) {
            res.push_back(std::make_tuple(xk[i], yk[i], NAN));
        }
        return res;
    }
};
double f(double x, double y, double z) {
    return z;
}
double fx(double x) {
    return 0.0;
}
double g(double x, double y, double z) {
    return ((-4 * x) * z - (4 * x * x + 2) * y);
}
double px(double x) {
    return 4 * x;
}
double qx(double x) {
    return 4 * x * x + 2;
}
int main() {
    cout.precision(5);
    cout << fixed;
    double h, eps;
    cin >> h >> eps;
    /*
    Краевые условия 3 рода
    alpha * y(a) + beta * y'(a) = y0
    delta * y(b) + gamma * y'(b) = y1
    */
    double a = 0, b = 2;
    double alpha = 0, beta = 1, y0 = 1;
    double delta = 4, gamma = -1, y1 = 23 * exp(-4);
    shooting mShooting(a, b, f, g, alpha, beta, y0, delta, gamma, y1);
    vector<tuple<double, double, double>> shootingSol = mShooting.solve(h, eps);
    cout << "Метод стрельбы:" << endl;
    printData(shootingSol);
    double shootingError = rungeRomberg(mShooting.solve(h, eps), mShooting.solve(h / 2, eps), 4);
    cout << "Погрешность вычислений = " << shootingError << "\n\n";
    finDif mFinDif(a, b, px, qx, fx, alpha, beta, y0, delta, gamma, y1);
    vector<tuple<double, double, double>> finDifSol = mFinDif.solve(h);
    cout << "Конечно-разностный метод:" << endl;
    printData(finDifSol);
    double finDifError = rungeRomberg(mFinDif.solve(h), mFinDif.solve(h / 2), 2);
    cout << "Погрешность вычислений = " << finDifError << "\n\n";
}
