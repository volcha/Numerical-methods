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

int main() {
    cout.precision(5);
    cout << fixed;
    int n;
    cin >> n;
    matrix_t<double> a(n);
    cin >> a;
    lu_t<double> lu_a(a);
    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vector<double> x = lu_a.solve(b);
    cout << lu_a;
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "\ndet(A) = " << lu_a.get_det() << endl;
    cout << "\nA^(-1): " << endl;
    cout << lu_a.inv_matrix();
}
