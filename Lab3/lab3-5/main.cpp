#include <iostream>
#include <vector>
#include <cmath>

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
using func = double(double);
double integrRec(double l, double r, double h, func f) {
    double res = 0;
    double x0 = l;
    double x1 = l + h;
    while (x0 < r) {
        res += h * f((x0 + x1) / 2);
        x0 = x1;
        x1 += h;
    }
    return res;
}
double integrTrap(double l, double r, double h, func f) {
    double res = 0;
    double x0 = l;
    double x1 = l + h;
    while (x0 < r) {
        res += h * (f(x0) + f(x1));
        x0 = x1;
        x1 += h;
    }
    return res / 2;
}
double integrSimp(double l, double r, double h, func f) {
    double res = 0;
    double x0 = l;
    double x1 = l + h;
    while (x0 < r) {
        vector<double> x = {x0, (x0 + x1) * 0.5, x1};
        vector<double> y = {f(x[0]), f(x[1]), f(x[2])};
        interLagrange lagr(x, y);
        res += lagr().integrate(x0, x1);
        x0 = x1;
        x1 += h;
    }
    return res / 3;
}
double rungeRomberg(double fh, double fhk, double k, double d) {
    return (fh - fhk) / (pow(k, d) - 1.0);
}
double f(double x) {
    return (x * x) / (625.0 - pow(x, 4));
}
int main() {
    double l, r;
    cin >> l >> r;
    double h1, h2;
    cin >> h1 >> h2;
    double rec1 = integrRec(l, r, h1, f);
    double trap1 = integrTrap(l, r, h1, f);
    double simp1 = integrSimp(l, r, h1, f);
    cout.precision(5);
    cout << "Интеграл по методу прямоугольников с шагом " << h1 << " = " << rec1 << endl;
    cout << "Интеграл по методу трапеций с шагом " << h1 << " = " << trap1 << endl;
    cout << "Интеграл по методу Симпсона с шагом " << h1 << " = " << simp1 << "\n\n";
    double rec2 = integrRec(l, r, h2, f);
    double trap2 = integrTrap(l, r, h2, f);
    double simp2 = integrSimp(l, r, h2, f);
    cout << "Интеграл по метод прямоугольников с шагом " << h2 << " = " << rec2 << endl;
    cout << "Интеграл по метод трапеций с шагом " << h2 << " = " << trap2 << endl;
    cout << "Интеграл по метод Симпсона с шагом " << h2 << " = " << simp2 << "\n\n";
    double recError = abs(rungeRomberg(rec1, rec2, h2 / h1, 2));
    double trapError = abs(rungeRomberg(trap1, trap2, h2 / h1, 2));
    double simpError = abs(rungeRomberg(simp1, simp2, h2 / h1, 2));
    cout << "Погрешность вычислений методом прямоугольников = " << recError << endl;
    cout << "Погрешность вычислений методом трапеций = " << trapError << endl;
    cout << "Погрешность вычислений методом Симпсона = " << simpError << endl;
}
