#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

template<class T>
std::vector<T> operator + (const std::vector<T> & a, const std::vector<T> & b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}
template<class T>
std::vector<T> operator - (const std::vector<T> & a, const std::vector<T> & b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}
template<class T>
class matrix_t {
private:
    size_t n, m;
    std::vector<std::vector<T>> data;
public:
    matrix_t() : n(1), m(1), data(1) {}
    matrix_t(size_t _n) : n(_n), m(_n) {
        data.resize(n, std::vector<T>(n));
    }
    matrix_t(size_t _n, size_t _m) : n(_n), m(_m) {
        data.resize(n, std::vector<T>(m));
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
    matrix_t<T> t() const {
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
            std::swap(data[i][k], data[j][k]);
        }
    }
    void swap_cols(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < n; ++k) {
            std::swap(data[k][i], data[k][j]);
        }
    }
    friend matrix_t<T> operator + (const matrix_t<T> & a, const matrix_t<T> & b) {
        if (a.rows() != b.rows() or a.cols() != b.cols()) {
            throw std::invalid_argument("Sizes does not match");
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
            throw std::invalid_argument("Sizes does not match");
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
            throw std::invalid_argument("Sizes does not match");
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
    friend std::vector<T> operator * (const matrix_t<T> & a, const std::vector<T> & b) {
        if (a.cols() != b.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        std::vector<T> c(n);
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
    std::vector<T> & operator [] (size_t i) {
        return data[i];
    }
    const std::vector<T> & operator [] (size_t i) const {
        return data[i];
    }
    friend std::istream & operator >> (std::istream & in, matrix_t<T> & matr) {
        for (size_t i = 0; i < matr.rows(); ++i) {
            for (size_t j = 0; j < matr.cols(); ++j) {
                in >> matr[i][j];
            }
        }
        return in;
    }
    friend std::ostream & operator << (std::ostream & out, const matrix_t<T> & matr) {
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
    std::vector<std::pair<size_t, size_t>> swaps;
    void do_swaps(std::vector<T> & x) {
        for (std::pair<size_t, size_t> elem : swaps) {
            std::swap(x[elem.first], x[elem.second]);
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
                std::pair<size_t, size_t> perm = std::make_pair(i, max_el_ind);
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
            throw std::invalid_argument("Matrix is not square");
        }
        l = matrix_t<T>::identity(matr.rows());
        u = matrix_t<T>(matr);
        decompose();
    }
    std::vector<T> solve(std::vector<T> b) {
        int n = b.size();
        do_swaps(b);
        std::vector<T> z(n);
        for (int i = 0; i < n; ++i) {
            T summary = b[i];
            for (int j = 0; j < i; ++j) {
                summary -= z[j] * l[i][j];
            }
            z[i] = summary;
        }
        std::vector<T> x(n);
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
            std::vector<T> b(n);
            b[i] = T(1);
            std::vector<T> x = solve(b);
            for (size_t j = 0; j < n; ++j) {
                res[j][i] = x[j];
            }
        }
        return res;
    }
    friend std::ostream & operator << (std::ostream & out, const lu_t<T> & lu) {
        out << "Matrix L:\n" << lu.l << "Matrix U:\n" << lu.u;
        return out;
    }
    ~lu_t() = default;
};
class polynom {
private:
    std::vector<double> data;
    constexpr static double EPS = 0.000000001;
    size_t n;
public:
    polynom() : data(1), n(1) {}
    polynom(int _n) : data(_n), n(_n) {}
    polynom(const std::vector<double> & coef) : data(coef), n(data.size()) {}
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
        polynom res(std::max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i) {
            res[i] += lhs[i];
        }
        for (size_t i = 0; i < rhs.size(); ++i) {
            res[i] += rhs[i];
        }
        return res;
    }
    friend polynom operator - (const polynom & lhs, const polynom & rhs) {
        polynom res(std::max(lhs.size(), rhs.size()));
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
    friend std::ostream & operator << (std::ostream & out, const polynom & poly) {
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
class minimal_square_t {
    size_t n;
    size_t m;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> a;
    std::vector<std::function<double(double)>> phi;
    double get(double x0) {
        double res = 0.0;
        for (size_t i = 0; i < m; ++i) {
            res += a[i] * phi[i](x0);
        }
        return res;
    }
    void build() {
        matrix_t<double> lhs(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                lhs[i][j] = phi[j](x[i]);
            }
        }
        matrix_t<double> lhs_t = lhs.t();
        lu_t<double> lhs_lu(lhs_t * lhs);
        std::vector<double> rhs = lhs_t * y;
        a = lhs_lu.solve(rhs);
    }
public:
    minimal_square_t(const std::vector<double> & _x, const std::vector<double> & _y, const std::vector<std::function<double(double)>> & _phi) {
        if (_x.size() != _y.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        x = _x;
        y = _y;
        n = _x.size();
        m = _phi.size();
        a.resize(m);
        phi = _phi;
        build();
    }
    double operator () (double x0) {
        return get(x0);
    }
    double mmse() {
        double res = 0;
        for (size_t i = 0; i < n; ++i) {
            res += pow(get(x[i]) - y[i], 2.0);
        }
        return res;
    }
    friend std::ostream & operator << (std::ostream & out, const minimal_square_t & item) {
        for (size_t i = 0; i < item.m; ++i) {
            if (i) {
                out << ' ';
            }
            out << item.a[i];
        }
        return out;
    }
};
double f0(double x0) {
    return 1.0;
}
double f1(double x0) {
    return x0;
}
double f2(double x0) {
    return x0 * x0;
}
int main() {
    int n;
    std::cin >> n;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> x[i];
    }
    for (int i = 0; i < n; ++i) {
        std::cin >> y[i];
    }
    std::cout.precision(6);
    std::cout << std::fixed;
    std::vector<std::function<double(double)>> phi1 = {f0, f1};
    minimal_square_t ms1(x, y, phi1);
    std::cout << "Коэффициенты приближающего многочлена 1-ой степени a0, a1: " << std::endl;
    std::cout << ms1 << std::endl;
    std::cout << "Сумма квадратов ошибок многочлена 1-ой степени = " << ms1.mmse() << "\n\n";
    std::vector<std::function<double(double)>> phi2 = {f0, f1, f2};
    minimal_square_t ms2(x, y, phi2);
    std::cout << "Коэффициенты приближающего многочлена 2-ой степени a0, a1, a2: " << std::endl;
    std::cout << ms2 << std::endl;
    std::cout << "Сумма квадратов ошибок многочлена 2-ой степени = " << ms2.mmse() << "\n\n";
}
