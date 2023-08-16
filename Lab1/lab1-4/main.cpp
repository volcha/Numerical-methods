#include <iostream>
#include <cmath>
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

class rotation {
private:
    static constexpr double GLOBAL_EPS = 1e-9;
    size_t n;
    matrix_t<double> a;
    matrix_t<double> v;
    double eps;

    matrix_t<double> create_rotation(size_t i, size_t j, double phi) {
        matrix_t<double> u = matrix_t<double>::identity(n);
        u[i][i] = cos(phi);
        u[i][j] = -sin(phi);
        u[j][i] = sin(phi);
        u[j][j] = cos(phi);
        return u;
    }

    double calc_phi(size_t i, size_t j) {
        if (abs(a[i][i] - a[j][j]) < GLOBAL_EPS) {
            return atan2(1.0, 1.0);
        }
        else {
            return 0.5 * atan2(2 * a[i][j], a[i][i] - a[j][j]);
        }
    }

    static double norm(const matrix_t<double> & m) {
        double res = 0;
        for (size_t i = 0; i < m.rows(); ++i) {
            for (size_t j = 0; j < m.cols(); ++j) {
                if (i == j) {
                    continue;
                }
                res += m[i][j] * m[i][j];
            }
        }
        return sqrt(res);
    }

    void build() {
        iter_count = 0;
        while (norm(a) > eps) {
            ++iter_count;
            size_t i = 0, j = 1;
            for (size_t ii = 0; ii < n; ++ii) {
                for (size_t jj = 0; jj < n; ++jj) {
                    if (ii == jj) {
                        continue;
                    }
                    if (abs(a[ii][jj]) > abs(a[i][j])) {
                        i = ii;
                        j = jj;
                    }
                }
            }
            double phi = calc_phi(i, j);
            matrix_t<double> u = create_rotation(i, j, phi);
            v = v * u;
            a = u.t() * a * u;
        }
    }
public:
    int iter_count;

    rotation(const matrix_t<double> & _a, double _eps) {
        if (_a.rows() != _a.cols()) {
            throw invalid_argument("Matrix is not square");
        }
        a = matrix_t<double>(_a);
        n = a.rows();
        eps = _eps;
        v = matrix_t<double>::identity(n);
        build();
    };

    vector<double> get_values() {
        vector<double> res(n);
        for (size_t i = 0; i < n; ++i) {
            res[i] = a[i][i];
        }
        return res;
    }

    matrix_t<double> get_vectors() {
        return v;
    }

    ~rotation() = default;
};

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix_t<double> a(n);
    cin >> a;
    rotation rot(a, eps);
    vector<double> lambda = rot.get_values();
    cout << "Собственные значения:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "l_" << i + 1 << " = " << lambda[i] << endl;
    }
    matrix_t<double> v = rot.get_vectors();
    cout << "\nСобственные векторы:" << endl << v;
    cout << "\nРешение получено за " << rot.iter_count << " итераций" << endl;
}
