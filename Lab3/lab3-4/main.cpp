#include <iostream>
#include <vector>
#include <cmath>
#include <exception>

using namespace std;

const double EPS = 0.000000001;
bool border(double a, double b) {
    return (a < b) or (abs(b - a) < EPS);
}
class derivative_t {
    int n;
    vector<double> x;
    vector<double> y;
public:
    derivative_t(const vector<double> & xl, const vector<double> & yl) {
        if (xl.size() != yl.size()) {
            throw invalid_argument("Несоответствие размерностей");
        }
        x = xl;
        y = yl;
        n = x.size();
    }
    double derivative1(double x0) {
        for (int i = 0; i < n - 2; ++i) {
            if (x[i] < x0 && border(x0, x[i + 1])) {
                double dydx1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
                double res = dydx1 + (dydx2 - dydx1) * (2.0 * x0 - x[i] - x[i + 1]) / (x[i + 2] - x[i]);
                return res;
            }
        }
        return NAN;
    }
    double derivative2(double x0) {
        for (int i = 0; i < n - 2; ++i) {
            if (x[i] < x0 && border(x0, x[i + 1])) {
                double dydx1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
                double res = 2.0 * (dydx2 - dydx1) / (x[i + 2] - x[i]);
                return res;
            }
        }
        return NAN;
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
    derivative_t f(x, y);
    cout << "f'(" << x0 << ") = " << f.derivative1(x0) << endl;
    cout << "f''(" << x0 << ") = " << f.derivative2(x0) << endl;
}
