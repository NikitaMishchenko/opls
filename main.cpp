#include <iostream>

#include "opls.h"

int main()
{
    ///input
    const std::vector<double> x = {0, 3.36588, 3.63719, 0.56448, -3.02721, -3.8357, -1.11766, 2.62795, 3.95743, 1.64847};
    const std::vector<double> y = {3.9561, 74.84479, 89.44289, 6.46668, -14.53888, -34.55881, 1.70531, 43.80101, 109.1294, 18.81613};
    int polynom_degree = 3;

    ///oputput

    std::vector<double> b(polynom_degree + 1);
    std::vector<double> a(polynom_degree + 1);

    opls O(polynom_degree, x.size());
        O.print_c();
        O.print_QSP();
        ///O.manage_opls(polynom_degree, x, y, b, a);

    ///O.Test();
}
