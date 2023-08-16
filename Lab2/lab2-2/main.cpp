#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
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

template<class T>
class lu_t {
private:
    const T EPS = 0.000001;
    matrix_t<T> l;
    matrix_t<T> u;
    T det;
    vector<pair<size_t, size_t>> swaps;

    void do_swaps(vector<T> & x) {
        for (pair<size_t, size_t> elem : swaps) {
            swap(x[elem.first], x[elem.second]);
        }
    }

    void decompose() {
        size_t n = u.rows();
        for (size_t i = 0; i < n; ++i) {
            size_t max_el_ind = i;
            for (size_t j = i + 1; j < n; ++j) {
                if (abs(u[j][i]) > abs(u[max_el_ind][i])) {
                    max_el_ind = j;
                }
            }
            if (max_el_ind != i) {
                pair<size_t, size_t> perm = make_pair(i, max_el_ind);
                swaps.push_back(perm);
                u.swap_rows(i, max_el_ind);
                l.swap_rows(i, max_el_ind);
                l.swap_cols(i, max_el_ind);
            }
            for (size_t j = i + 1; j < n; ++j) {
                if (abs(u[i][i]) < EPS) {
                    continue;
                }
                T mu = u[j][i] / u[i][i];
                l[j][i] = mu;
                for (size_t k = 0; k < n; ++k) {
                    u[j][k] -= mu * u[i][k];
                }
            }
        }
        det = (swaps.size() & 1 ? -1 : 1);
        for (size_t i = 0; i < n; ++i) {
            det *= u[i][i];
        }
    }
public:
    lu_t(const matrix_t<T> & matr) {
        if (matr.rows() != matr.cols()) {
            throw invalid_argument("Matrix is not square");
        }
        l = matrix_t<T>::identity(matr.rows());
        u = matrix_t<T>(matr);
        decompose();
    }

    vector<T> solve(vector<T> b) {
        int n = b.size();
        do_swaps(b);
        vector<T> z(n);
        for (int i = 0; i < n; ++i) {
            T summary = b[i];
            for (int j = 0; j < i; ++j) {
                summary -= z[j] * l[i][j];
            }
            z[i] = summary;
        }
        vector<T> x(n);
        for (int i = n - 1; i >= 0; --i) {
            if (abs(u[i][i]) < EPS) {
                continue;
            }
            T summary = z[i];
            for (int j = n - 1; j > i; --j) {
                summary -= x[j] * u[i][j];
            }
            x[i] = summary / u[i][i];
        }
        return x;
    }

    T get_det() {
        return det;
    }

    matrix_t<T> inv_matrix() {
        size_t n = l.rows();
        matrix_t<T> res(n);
        for (size_t i = 0; i < n; ++i) {
            vector<T> b(n);
            b[i] = T(1);
            vector<T> x = solve(b);
            for (size_t j = 0; j < n; ++j) {
                res[j][i] = x[j];
            }
        }
        return res;
    }

    friend ostream & operator << (ostream & out, const lu_t<T> & lu) {
        out << "Matrix L:\n" << lu.l << "\nMatrix U:\n" << lu.u << '\n';
        return out;
    }

    ~lu_t() = default;
};

int iter_count = 0;

const double a = 1;

double f1(double x1, double x2) {
    return x1 * x1 - 2 * log10(x2) - 1;
}

double f2(double x1, double x2) {
    return x1 * x1 - a * x1 * x2 + a;
}

double phi1(double x1, double x2) {
    return sqrt(2 * log10(x2) + 1);
}

double phi1_der(double x1, double x2) {
    return 1 / (x2 * sqrt(log(10)) * sqrt(2 * log(x2) + log(10)));
}

double phi2(double x1, double x2) {
    return (x1 * x1 + a)/(x1 * a);
}

double phi2_der(double x1, double x2) {
    return (a - 1 / (x1 * x1 * a));
}

double phi(double x1, double x2) {
    return phi1_der(x1, x2) * phi2_der(x1, x2);
}

pair<double, double> iter_solve(double l1, double r1, double l2, double r2, double eps) {
    iter_count = 0;
    double x1_k = r1;
    double x2_k = r2;
    double q = -1;
    q = max(q, abs(phi(l1, r1)));
    q = max(q, abs(phi(l1, r2)));
    q = max(q, abs(phi(l2, r1)));
    q = max(q, abs(phi(l2, r2)));
    double eps_coeff = q / (1 - q);
    double dx = 1;
    while (dx > eps) {
        double x1_k1 = phi1(x1_k, x2_k);
        double x2_k1 = phi2(x1_k, x2_k);
        dx = eps_coeff * (abs(x1_k1 - x1_k) + abs(x2_k1 - x2_k));
        ++iter_count;
        x1_k = x1_k1;
        x2_k = x2_k1;
    }
    return make_pair(x1_k, x2_k);
}

matrix_t<double> Jacobi(double x1, double x2) {
    matrix_t<double> result(2);
    result[0][0] = 2 * x1;
    result[0][1] = -2 / (x2 * log(10));
    result[1][0] = 2 * x1 - a * x2;
    result[1][1] = -a * x1;
    return result;
}

double norm(const vector<double> &vec) {
    double res = 0;
    for (auto elem: vec) {
        res = max(res, abs(elem));
    }
    return res;
}

pair<double, double> newton_solve(double x1_0, double x2_0, double eps) {
    iter_count = 0;
    vector<double> x_k = {x1_0, x2_0};
    double dx = 1;
    while (dx > eps) {
        double x1 = x_k[0];
        double x2 = x_k[1];
        lu_t<double> J(Jacobi(x1, x2));
        vector<double> f_k = {f1(x1,x2), f2(x1,x2)};
        vector<double> d_x = J.solve(f_k);
        vector<double> x_k1 = x_k - d_x;
        dx = norm(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    }
    return make_pair(x_k[0], x_k[1]);
}

int main() {
    cout.precision(8);
    cout << fixed;
    double l1, r1, l2, r2, eps;
    cin >> l1 >> r1 >> l2 >> r2 >> eps;
    auto [x0, y0] = iter_solve(l1, r1, l2, r2, eps);
    cout << "x1 = " << x0 << ", x2 = " << y0 << endl;
    cout << "Решение методом простой итерации получено за " << iter_count << " итераций" << "\n\n";
    auto [x0_n, y0_n] = newton_solve(r1, r2, eps);
    cout << "x1 = " << x0_n << ", x2 = " << y0_n << endl;
    cout << "Решение методом Ньютона получено за " << iter_count << " итераций" << "\n\n";
}
