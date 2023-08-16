#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class polynom {
private:
    vector<double> data;
    constexpr static double EPS = 0.000000001;
    size_t n;
public:
    polynom() : data(1), n(1) {}
    polynom(int _n) : data(_n), n(_n) {}
    polynom(const vector<double> & coef) : data(coef), n(data.size()) {}
    size_t size() const {
        return n;
    }
    double & operator [] (size_t id) {
        return data[id];
    }
    const double & operator [] (size_t id) const {
        return data[id];
    }
    friend polynom operator + (const polynom & lhs, const polynom & rhs) {
        polynom res(max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i) {
            res[i] += lhs[i];
        }
        for (size_t i = 0; i < rhs.size(); ++i) {
            res[i] += rhs[i];
        }
        return res;
    }
    friend polynom operator - (const polynom & lhs, const polynom & rhs) {
        polynom res(max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i) {
            res[i] += lhs[i];
        }
        for (size_t i = 0; i < rhs.size(); ++i) {
            res[i] -= rhs[i];
        }
        return res;
    }
    friend polynom operator * (double lambda, const polynom & p) {
        polynom res(p);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] *= lambda;
        }
        return res;
    }
    friend polynom operator / (const polynom & p, double lambda) {
        polynom res(p);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] /= lambda;
        }
        return res;
    }
    friend polynom operator * (const polynom & lhs, const polynom & rhs) {
        polynom res(lhs.size() + rhs.size());
        for (size_t i = 0; i < lhs.size(); ++i) {
            for (size_t j = 0; j < rhs.size(); ++j) {
                res[i + j] += lhs[i] * rhs[j];
            }
        }
        while (res.n > 1 and abs(res.data.back()) < EPS) {
            res.data.pop_back();
            --res.n;
        }
        return res;
    }
    polynom integrate() {
        polynom res(n + 1);
        for (size_t i = 1; i < n + 1; ++i) {
            res.data[i] = data[i - 1] / (double)i;
        }
        return res;
    }
    double integrate(double l, double r) {
        polynom F = integrate();
        return F(r) - F(l);
    }
    polynom derivative() {
        polynom res(n - 1);
        for (size_t i = 1; i < n; ++i) {
            res[i - 1] = data[i] * i;
        }
        return res;
    }
    double operator () (double x) {
        double res = 0.0;
        double xi = 1.0;
        for (double elem : data) {
            res += elem * xi;
            xi *= x;
        }
        return res;
    }
    friend ostream & operator << (ostream & out, const polynom & poly) {
        bool flag = false;
        int deg = 0;
        for (double elem : poly.data) {
            if (!(abs(elem) < EPS)) {
                if (flag and deg) {
                    out << (elem > EPS ? " + " : " - ");
                }
                out << abs(elem);
                flag = true;
                if (deg) {
                    out << " * x";
                    if (deg > 1) {
                        out << " ^ " << deg;
                    }
                }
            }
            ++deg;
        }
        if (!flag) {
            out << 0;
        }
        return out;
    }
    ~polynom() = default;
};
class interLagrange {
    vector<double> x;
    vector<double> y;
    size_t n;
public:
    interLagrange(const vector<double> & _x, const vector<double> & _y) : x(_x), y(_y), n(x.size()) {};
    polynom operator () () {
        polynom res(vector<double>({0}));
        for (size_t i = 0; i < n; ++i) {
            polynom li(vector<double>({1}));
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                polynom xij(vector<double>({-x[j], 1}));
                li = li * xij / (x[i] - x[j]);
            }
            res = res + y[i] * li;
        }
        return res;
    }
};
class interNewton {
private:
    vector<double> x;
    vector<double> y;
    size_t n;
    vector<vector<bool>> calc;
    vector<vector<double>> memo;
    double f(int l, int r) {
        if (calc[l][r]) {
            return memo[l][r];
        }
        calc[l][r] = true;
        double res;
        if (l + 1 == r) {
            res = (y[l] - y[r]) / (x[l] - x[r]);
        } else {
            res = (f(l, r - 1) - f(l + 1, r)) / (x[l] - x[r]);
        }
        return memo[l][r] = res;
    }
public:
    interNewton(const vector<double> & _x, const vector<double> & _y) : x(_x), y(_y), n(x.size()) {
        calc.resize(n, vector<bool>(n));
        memo.resize(n, vector<double>(n));
    };
    polynom operator () () {
        polynom li(vector<double>({-x[0], 1}));
        polynom res(vector<double>({y[0]}));
        int r = 0;
        for (size_t i = 1; i < n; ++i) {
            res = res + f(0, ++r) * li;
            li = li * polynom(vector<double>({-x[i], 1}));
        }
        return res;
    }
};
int main() {
    int n;
    cin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        cin >> x[i];
        y[i] = asin(x[i]) + x[i];
    }
    double x_error;
    cin >> x_error;
    interLagrange myLagrange(x, y);
    polynom lagrange = myLagrange();
    cout << "Интерполяционный многочлен Лагранжа:\n" << lagrange << endl;
    cout << "Погрешность в точке X' = " << abs(lagrange(x_error) - asin(x_error) - x_error) << "\n\n";
    interNewton myNewton(x, y);
    polynom newton = myNewton();
    cout << "Интерполяционный многочлен Ньютона:\n" << newton << endl;
    cout << "Погрешность в точке X' = " << abs(newton(x_error) - asin(x_error) - x_error) << "\n\n";
}
