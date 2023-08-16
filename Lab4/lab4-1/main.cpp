#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <tuple>

using namespace std;

const double EPS = 0.000000001;
void printData(const vector<tuple<double, double, double>> & v) {
    cout << "x = [";
    for (int i = 0; i < v.size(); ++i) {
        if (i) {
            cout << ", ";
        }
        cout << get<0>(v[i]);
    }
    cout << "]\n";
    cout << "y = [";
    for (int i = 0; i < v.size(); ++i) {
        if (i) {
            cout << ", ";
        }
        cout << get<1>(v[i]);
    }
    cout << "]\n";
}
bool border(double a, double b) {
    return (a < b) or (abs(b - a) < EPS);
}
class euler {
private:
    double l, r, y0, z0;
    function<double(double, double, double)> f, g;
public:
    euler(const double _l, const double _r,
        const function<double(double, double, double)> _f, const function<double(double, double, double)> _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}
    vector<tuple<double, double, double>> solve(double h) {
        double xi = l;
        double yi = y0;
        double zi = z0;
        vector<tuple<double, double, double>> res;
        res.push_back(make_tuple(xi, yi, zi));
        while (border(xi + h, r)) {
            double dy = h * f(xi, yi, zi);
            double dz = h * g(xi, yi, zi);
            xi += h;
            yi += dy;
            zi += dz;
            res.push_back(make_tuple(xi, yi, zi));
        }
        return res;
    }
};
class rungeKutta {
private:
    double l, r, y0, z0;
    function<double(double, double, double)> f, g;
public:
    rungeKutta(const double _l, const double _r,
        const function<double(double, double, double)> _f, const function<double(double, double, double)> _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}
    vector<tuple<double, double, double>> solve(double h) {
        double xi = l;
        double yi = y0;
        double zi = z0;
        vector<tuple<double, double, double>> res;
        res.push_back(make_tuple(xi, yi, zi));
        while (border(xi + h, r)) {
            double k1 = h * f(xi, yi, zi);
            double l1 = h * g(xi, yi, zi);
            double k2 = h * f(xi + 0.5 * h, yi + 0.5 * k1, zi + 0.5 * l1);
            double l2 = h * g(xi + 0.5 * h, yi + 0.5 * k1, zi + 0.5 * l1);
            double k3 = h * f(xi + 0.5 * h, yi + 0.5 * k2, zi + 0.5 * l2);
            double l3 = h * g(xi + 0.5 * h, yi + 0.5 * k2, zi + 0.5 * l2);
            double k4 = h * f(xi + h, yi + k3, zi + l3);
            double l4 = h * g(xi + h, yi + k3, zi + l3);
            double dy = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            double dz = (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;
            xi += h;
            yi += dy;
            zi += dz;
            res.push_back(make_tuple(xi, yi, zi));
        }
        return res;
    }
};
class adams {
private:
    double l, r, y0, z0;
    function<double(double, double, double)> f, g;
public:
    adams(const double _l, const double _r,
        const function<double(double, double, double)> _f, const function<double(double, double, double)> _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}
    double calc_tuple(function<double(double, double, double)> f, tuple<double, double, double> xyz) {
        return f(get<0>(xyz), get<1>(xyz), get<2>(xyz));
    }
    vector<tuple<double, double, double>> solve(double h) {
        if (l + 3.0 * h > r) {
            throw invalid_argument("Шаг слишком большой");
        }
        rungeKutta firstPoints(l, l + 3.0 * h, f, g, y0, z0);
        vector<tuple<double, double, double>> res = firstPoints.solve(h);
        size_t cnt = res.size();
        double xk = get<0>(res.back());
        double yk = get<1>(res.back());
        double zk = get<2>(res.back());
        while (border(xk + h, r)) {
            double dy = (h / 24.0) * (55.0 * calc_tuple(f, res[cnt - 1]) - 59.0 * calc_tuple(f, res[cnt - 2]) + 37.0 * calc_tuple(f, res[cnt - 3]) - 9.0 * calc_tuple(f, res[cnt - 4]));
            double dz = (h / 24.0) * (55.0 * calc_tuple(g, res[cnt - 1]) - 59.0 * calc_tuple(g, res[cnt - 2]) + 37.0 * calc_tuple(g, res[cnt - 3]) - 9.0 * calc_tuple(g, res[cnt - 4]));
            double xk1 = xk + h;
            double yk1 = yk + dy;
            double zk1 = zk + dz;
            res.push_back(make_tuple(xk1, yk1, zk1));
            ++cnt;
            dy = (h / 24.0) * (9.0 * calc_tuple(f, res[cnt - 1]) + 19.0 * calc_tuple(f, res[cnt - 2]) - 5.0 * calc_tuple(f, res[cnt - 3])
 + 1.0 * calc_tuple(f, res[cnt - 4]));
            dz = (h / 24.0) * (9.0 * calc_tuple(g, res[cnt - 1]) + 19.0 * calc_tuple(g, res[cnt - 2]) - 5.0 * calc_tuple(g, res[cnt - 3])
 + 1.0 * calc_tuple(g, res[cnt - 4]));
            xk += h;
            yk += dy;
            zk += dz;
            res.pop_back();
            res.push_back(make_tuple(xk, yk, zk));
        }
        return res;
    }
};
double rungeRomberg(const vector<tuple<double, double, double>> & y_2h, const vector<tuple<double, double, double>> & y_h, double p) {
    double coef = 1.0 / (pow(2, p) - 1.0);
    double res = 0.0;
    for (int i = 0; i < y_2h.size(); ++i) {
        res = max(res, coef * abs(get<1>(y_2h[i]) - get<1>(y_h[2 * i])));
    }
    return res;
}
double f(double x, double y, double z) {
    return z;
}
double g(double x, double y, double z) {
    return (12*y/(x*x));
}
int main() {
    cout.precision(5);
    cout << fixed;
    double l, r, y0, z0, h;
    cin >> l >> r;
    cin >> y0 >> z0 >> h;
    euler mEuler(l, r, f, g, y0, z0);
    vector<tuple<double, double, double>> solEuler = mEuler.solve(h);
    cout << "Метод Эйлера:" << endl;
    printData(solEuler);
    double eulerError = rungeRomberg(mEuler.solve(h), mEuler.solve(h / 2), 1);
    cout << "Погрешность  = " << eulerError << "\n\n";
    rungeKutta mRunge(l, r, f, g, y0, z0);
    vector<tuple<double, double, double>> solRunge = mRunge.solve(h);
    cout << "Метод Рунге-Кутты:" << endl;
    printData(solRunge);
    double rungeError = rungeRomberg(mRunge.solve(h), mRunge.solve(h / 2), 4);
    cout << "Погрешность = " << rungeError << "\n\n";
    adams mAdams(l, r, f, g, y0, z0);
    vector<tuple<double, double, double>> solAdams = mAdams.solve(h);
    cout << "Метод Адамса:" << endl;
    printData(solAdams);
    double adamsError = rungeRomberg(mAdams.solve(h), mAdams.solve(h / 2), 4);
    cout << "Погрешность = " << adamsError << endl;
}
