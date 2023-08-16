#include <iostream>
#include <vector>

using namespace std;

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

class solver_si_ze {
private:
    matrix_t<double> a;
    size_t n;
    double eps;
    static constexpr double INF = 1e18;
public:
    int iter_count;

    solver_si_ze(const matrix_t<double> & _a, double _eps = 0.000001) {
        if (_a.rows() != _a.cols()) {
            throw invalid_argument("Matrix is not square");
        }
        a = matrix_t<double>(_a);
        n = a.rows();
        eps = _eps;
    }

    static double norm(const vector<double> & v) {
        double res = -INF;
        for (double elem : v) {
            res = max(res, abs(elem));
        }
        return res;
    }

    static double norm(const matrix_t<double> & m) {
        double res = -INF;
        for (size_t i = 0; i < m.rows(); ++i) {
            double s = 0;
            for (double elem : m[i]) {
                s += abs(elem);
            }
            res = max(res, s);
        }
        return res;
    }

    pair<matrix_t<double>, vector<double>> precalc_ab(const vector<double> & b, matrix_t<double> & alpha, vector<double> & beta) {
        for (size_t i = 0; i < n; ++i) {
            beta[i] = b[i] / a[i][i];
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    alpha[i][j] = -a[i][j] / a[i][i];
                }
            }
        }
        return make_pair(alpha, beta);
    }

    vector<double> solve_simple(const vector<double> & b) {
        matrix_t<double> alpha(n);
        vector<double> beta(n);
        precalc_ab(b, alpha, beta);
        double eps_coef = 1.0;
        if (norm(alpha) - 1.0 < eps) {
            eps_coef = norm(alpha) / (1.0 - norm(alpha));
        }
        double eps_k = 1.0;
        vector<double> x(beta);
        iter_count = 0;
        while (eps_k > eps) {
            vector<double> x_k = beta + alpha * x;
            eps_k = eps_coef * norm(x_k - x);
            x = x_k;
            ++iter_count;
        }
        return x;
    }

    vector<double> zeidel(const vector<double> & x, const matrix_t<double> & alpha, const vector<double> & beta) {
        vector<double> x_k(beta);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                x_k[i] += x_k[j] * alpha[i][j];
            }
            for (size_t j = i; j < n; ++j) {
                x_k[i] += x[j] * alpha[i][j];
            }
        }
        return x_k;
    }

    vector<double> solve_zeidel(const vector<double> & b) {
        matrix_t<double> alpha(n);
        vector<double> beta(n);
        precalc_ab(b, alpha, beta);
        matrix_t<double> c(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                c[i][j] = alpha[i][j];
            }
        }
        double eps_coef = 1.0;
        if (norm(alpha) - 1.0 < eps) {
            eps_coef = norm(c) / (1.0 - norm(alpha));
        }
        double eps_k = 1.0;
        vector<double> x(beta);
        iter_count = 0;
        while (eps_k > eps) {
            vector<double> x_k = zeidel(x, alpha, beta);
            eps_k = eps_coef * norm(x_k - x);
            x = x_k;
            ++iter_count;
        }
        return x;
    }

    ~solver_si_ze() = default;
};

int main() {
    cout.precision(5);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix_t<double> a(n);
    cin >> a;
    solver_si_ze solver(a, eps);
    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vector<double> x = solver.solve_simple(b);
    cout << "Метод простых итераций" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "Решени получено за " << solver.iter_count << " итераций" << "\n\n";
    vector<double> ze = solver.solve_zeidel(b);
    cout << "Метод Зейделя" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << ze[i] << endl;
    }
    cout << "Решени получено за " << solver.iter_count << " итераций" << endl;
}
