#include <iostream>
#include <vector>

using namespace std;

template<class T>
class tridiag_t {
private:
    const T EPS = 0.000001;
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

int main() {
    cout.precision(5);
    cout << fixed;
    int n;
    cin >> n;
    tridiag_t<double> tridiag_a(n);
    cin >> tridiag_a;
    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vector<double> x = tridiag_a.solve(b);
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
