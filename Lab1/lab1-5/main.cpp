#include <iostream>
#include <complex>
#include <vector>

using namespace std;

const double EPS = 1e-9;

template<class T>
vector<T> operator + (const vector<T> & a, const vector<T> & b) {
    size_t n = a.size();
    vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}

template<class T>
vector<T> operator - (const vector<T> & a, const vector<T> & b) {
    size_t n = a.size();
    vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}

template<class T>
class matrix_t {
private:
    size_t n, m;
    vector<vector<T>> data;
public:
    matrix_t() : n(1), m(1), data(1) {}

    matrix_t(size_t _n) : n(_n), m(_n) {
        data.resize(n, vector<T>(n));
    }

    matrix_t(size_t _n, size_t _m) : n(_n), m(_m) {
        data.resize(n, vector<T>(m));
    }

    matrix_t(const matrix_t<T> & other) {
        n = other.n;
        m = other.m;
        data = other.data;
    }

    matrix_t<T> & operator = (const matrix_t<T> & other) {
        if (this == &other) {
            return *this;
        }
        n = other.n;
        m = other.m;
        data = other.data;
        return *this;
    }

    static matrix_t<T> identity(size_t n) {
        matrix_t<T> res(n, n);
        for (size_t i = 0; i < n; ++i) {
            res[i][i] = T(1);
        }
        return res;
    }

    matrix_t t() const {
        matrix_t<T> res(m, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[j][i] = data[i][j];
            }
        }
        return res;
    }

    size_t rows() const {
        return n;
    }

    size_t cols() const {
        return m;
    }

    void swap_rows(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < m; ++k) {
            swap(data[i][k], data[j][k]);
        }
    }

    void swap_cols(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < n; ++k) {
            swap(data[k][i], data[k][j]);
        }
    }

    friend matrix_t<T> operator + (const matrix_t<T> & a, const matrix_t<T> & b) {
        if (a.rows() != b.rows() or a.cols() != b.cols()) {
            throw invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        matrix_t<T> res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = a[i][j] + b[i][j];
            }
        }
        return res;
    }

    friend matrix_t<T> operator - (const matrix_t<T> & a, const matrix_t<T> & b) {
        if (a.rows() != b.rows() or a.cols() != b.cols()) {
            throw invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        matrix_t<T> res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = a[i][j] - b[i][j];
            }
        }
        return res;
    }

    friend matrix_t<T> operator * (const matrix_t<T> & a, const matrix_t<T> & b) {
        if (a.cols() != b.rows()) {
            throw invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t k = a.cols();
        size_t m = b.cols();
        matrix_t<T> res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t ii = 0; ii < k; ++ii) {
                    res[i][j] += a[i][ii] * b[ii][j];
                }
            }
        }
        return res;
    }

    friend vector<T> operator * (const matrix_t<T> & a, const vector<T> & b) {
        if (a.cols() != b.size()) {
            throw invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        vector<T> c(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                c[i] += a[i][j] * b[j];
            }
        }
        return c;
    }

    friend matrix_t<T> operator * (T lambda, const matrix_t<T> & a) {
        size_t n = a.rows();
        size_t m = a.cols();
        matrix_t<T> res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = lambda * a[i][j];
            }
        }
        return res;
    }

    vector<T> & operator [] (size_t i) {
        return data[i];
    }

    const vector<T> & operator [] (size_t i) const {
        return data[i];
    }

    friend istream & operator >> (istream & in, matrix_t<T> & matr) {
        for (size_t i = 0; i < matr.rows(); ++i) {
            for (size_t j = 0; j < matr.cols(); ++j) {
                in >> matr[i][j];
            }
        }
        return in;
    }

    friend ostream & operator << (ostream & out, const matrix_t<T> & matr) {
        for (size_t i = 0; i < matr.rows(); ++i) {
            for (size_t j = 0; j < matr.cols(); ++j) {
                if (j) {
                    out << ", ";
                }
                out << matr[i][j];
            }
            out << '\n';
        }
        return out;
    }

    ~matrix_t() = default;
};

class qr {
private:
    static constexpr double INF = 1e18;
    static constexpr std::complex<double> COMPLEX_INF = std::complex<double>(INF, INF);
    size_t n;
    double eps;
    matrix_t<double> a;
    vector<std::complex<double>> eigen;

    double sign(double x) {
        if (x < eps) {
            return -1.0;
        }
        else if (x > eps) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }

    matrix_t<double> vvt(const vector<double> & b) {
        size_t n_b = b.size();
        matrix_t<double> res(n_b);
        for (size_t i = 0; i < n_b; ++i) {
            for (size_t j = 0; j < n_b; ++j) {
                res[i][j] = b[i] * b[j];
            }
        }
        return res;
    }

    double vtv(const vector<double> & v) {
        double res = 0;
        for (double elem : v) {
            res += elem * elem;
        }
        return res;
    }

    double norm(const vector<double> & v) {
        return sqrt(vtv(v));
    }

    matrix_t<double> householder(const vector<double> & b, int id) {
        vector<double> v(b);
        v[id] += sign(b[id]) * norm(b);
        return matrix_t<double>::identity(n) - (2.0 / vtv(v)) * vvt(v);
    }

    pair<std::complex<double>, std::complex<double>> solve_sq(double a11, double a12, double a21, double a22) {
        double a = 1.0;
        double b = -(a11 + a22);
        double c = a11 * a22 - a12 * a21;
        double d_sq = b * b - 4.0 * a * c;
        if (d_sq > eps) {
            std::complex<double> bad(NAN, NAN);
            return make_pair(bad, bad);
        }
        std::complex<double> d(0.0, sqrt(-d_sq));
        std::complex<double> x1 = (-b + d) / (2.0 * a);
        std::complex<double> x2 = (-b - d) / (2.0 * a);
        return make_pair(x1, x2);
    }

    void calc_eigen() {
        for (size_t i = 0; i < n; ++i) {
            if (i < n - 1 and !(abs(a[i + 1][i]) < eps)) {
                auto [l1, l2] = solve_sq(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1]);
                if (isnan(l1.real())) {
                    eigen[i] = COMPLEX_INF;
                    ++i;
                    eigen[i] = COMPLEX_INF;
                    continue;
                }
                eigen[i] = l1;
                eigen[++i] = l2;
            }
            else {
                eigen[i] = a[i][i];
            }
        }
    }

    bool check_diag() {
        for (size_t i = 0; i < n; ++i) {
            double col_sum = 0;
            for (size_t j = i + 2; j < n; ++j) {
                col_sum += a[j][i] * a[j][i];
            }
            double norm = sqrt(col_sum);
            if (!(norm < eps)) {
                return false;
            }
        }
        return true;
    }

    bool check_eps() {
        if (!check_diag()) {
            return false;
        }
        vector<std::complex<double>> prev_eigen(eigen);
        calc_eigen();
        for (size_t i = 0; i < n; ++i) {
            bool bad = (std::norm(eigen[i] - COMPLEX_INF) < eps);
            if (bad) {
                return false;
            }
            double delta = std::norm(eigen[i] - prev_eigen[i]);
            if (delta > eps) {
                return false;
            }
        }
        return true;
    }

    void build() {
        iter_count = 0;
        while (!check_eps()) {
            ++iter_count;
            matrix_t<double> q = matrix_t<double>::identity(n);
            matrix_t<double> r(a);
            for (size_t i = 0; i < n - 1; ++i) {
                vector<double> b(n);
                for (size_t j = i; j < n; ++j) {
                    b[j] = r[j][i];
                }
                matrix_t<double> h = householder(b, i);
                q = q * h;
                r = h * r;
            }
            a = r * q;
        }
    }
public:
    int iter_count;

    qr(const matrix_t<double> & _a, double _eps) {
        if (_a.rows() != _a.cols()) {
            throw invalid_argument("Matrix is not square");
        }
        n = _a.rows();
        a = matrix_t<double>(_a);
        eps = _eps;
        eigen.resize(n, COMPLEX_INF);
        build();
    };

    vector<std::complex<double>> get_values() {
        calc_eigen();
        return eigen;
    }

    ~qr() = default;
};

string format(complex<double> c) {
    if (abs(c.imag()) < EPS) {
        return to_string(c.real());
    }
    else {
        return to_string(c.real()) + " + i * (" + to_string(c.imag()) + ")";
    }
}

int main() {
    int n;
    double eps;
    cin >> n >> eps;
    matrix_t<double> a(n);
    cin >> a;
    qr my_qr(a, eps);
    vector<complex<double>> lambda = my_qr.get_values();
    cout << "Собственные значения:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "l_" << i + 1 << " = " << format(lambda[i]) << endl;
    }
    cout << "\nРешение получено за " << my_qr.iter_count << " итерации" << endl;
}
