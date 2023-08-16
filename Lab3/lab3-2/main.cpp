#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template<class T>
class tridiag_t {
private:
    const double EPS = 0.000001;
    int n;
    vector<T> a;
    vector<T> b;
    vector<T> c;
public:
    tridiag_t(const int & _n) : n(_n), a(n), b(n), c(n) {}
    tridiag_t(const vector<T> & _a, const vector<T> & _b, const vector<T> & _c) {
        if (!(_a.size() == _b.size() and _a.size() == _c.size())) {
            throw invalid_argument("Sizes of a, b, c are invalid");
        }
        n = _a.size();
        a = _a;
        b = _b;
        c = _c;
    }
    vector<T> solve(const vector<T> & d) {
        int m = d.size();
        if (n != m) {
            throw invalid_argument("Size of vector d is invalid");
        }
        vector<T> p(n);
        p[0] = -c[0] / b[0];
        vector<T> q(n);
        q[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
            q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
        }
        vector<T> x(n);
        x.back() = q.back();
        for (int i = n - 2; i >= 0; --i) {
            x[i] = p[i] * x[i + 1] + q[i];
        }
        return x;
    }
    friend istream & operator >> (istream & in, tridiag_t<T> & tridiag) {
        in >> tridiag.b[0] >> tridiag.c[0];
        for (int i = 1; i < tridiag.n - 1; ++i) {
            in >> tridiag.a[i] >> tridiag.b[i] >> tridiag.c[i];
        }
        in >> tridiag.a.back() >> tridiag.b.back();
        return in;
    }
    ~tridiag_t() = default;
};
class cub_spline_t {
    size_t n;
    vector<double> x;
    vector<double> y;
    vector<double> a, b, c, d;
    void buildSpline() {
        vector<double> h(n + 1);
        h[0] = NAN;
        for (size_t i = 1; i <= n; ++i) {
            h[i] = x[i] - x[i - 1];
        }
        vector<double> func1(n - 1);
        vector<double> func2(n - 1);
        vector<double> func3(n - 1);
        vector<double> func4(n - 1);
        for (size_t i = 2; i <= n; ++i) {
            func1[i - 2] = h[i - 1];
            func2[i - 2] = 2.0 * (h[i - 1] + h[i]);
            func3[i - 2] = h[i];
            func4[i - 2] = 3.0 * ((y[i] - y[i - 1]) / h[i] - (y[i - 1] - y[i - 2]) / h[i - 1]);
        }
        func1[0] = 0.0;
        func3.back() = 0.0;
        tridiag_t<double> systemOfFunc(func1, func2, func3);
        vector<double> cSolved = systemOfFunc.solve(func4);
        for (size_t i = 2; i <= n; ++i) {
            c[i] = cSolved[i - 2];
        }
        for (size_t i = 1; i <= n; ++i) {
            a[i] = y[i - 1];
        }
        for (size_t i = 1; i < n; ++i) {
            b[i] = (y[i] - y[i - 1]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
        c[1] = 0.0;
        b[n] = (y[n] - y[n - 1]) / h[n] - (2.0 / 3.0) * h[n] * c[n];
        d[n] = -c[n] / (3.0 * h[n]);
    }
public:
    cub_spline_t(const vector<double> & _x, const vector<double> & _y) {
        if (_x.size() != _y.size()) {
            throw invalid_argument("Sizes does not match");
        }
        x = _x;
        y = _y;
        n = x.size() - 1;
        a.resize(n + 1);
        b.resize(n + 1);
        c.resize(n + 1);
        d.resize(n + 1);
        buildSpline();

    }
    double operator () (double x0) {
        for (size_t i = 1; i <= n; ++i) {
            if (x[i - 1] <= x0 and x0 <= x[i]) {
                double x1 = x0 - x[i - 1];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                return a[i] + b[i] * x1 + c[i] * x2 + d[i] * x3;
            }
        }
        return NAN;
    }
    friend ostream & operator << (ostream & out, const cub_spline_t & spline) {
        for (size_t i = 1; i <= spline.n; ++i) {
            out << "i = " << i << ", a = " << spline.a[i] << ", b = " << spline.b[i] << ", c = " << spline.c[i] << ", d = " << spline.d[i] << '\n';
        }
        return out;
    }
};
int main() {
    int n;
    cin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        cin >> x[i];
    }
    for (int i = 0; i < n; ++i) {
        cin >> y[i];
    }
    double x0;
    cin >> x0;
    cout.precision(6);
    cub_spline_t func(x, y);
    cout << "Cплайн:\n" << func << endl;
    cout << "f(X') = " << func(x0) << endl;
}
