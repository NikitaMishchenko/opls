#ifndef OPLS_H_INCLUDED
#define OPLS_H_INCLUDED

///REFERENCE https://ru.stackoverflow.com/questions/383716/%D0%9C%D0%9D%D0%9A-%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC

#include <vector>

class opls
{
public:

    std::vector<std::vector<double>> P; /// values of polynoms in mesh nodes
    std::vector<double> Q;              /// sum Pi^2
    std::vector<double> S;              /// sum xi*Pi^2
    int polynom_degree = 0;

    opls(int polynom_degree, int data_size); /// init size of P, Q, S

    ///METHODS
    void QS_count(const std::vector<double> &x);

    void manage_opls(int polynom_order, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a);

    void Test();

    ///EXTRA
    void print_QSP();

};

#endif // OPLS_H_INCLUDED
