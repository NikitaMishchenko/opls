#ifndef OPLS_H_INCLUDED
#define OPLS_H_INCLUDED

///REFERENCE https://ru.stackoverflow.com/questions/383716/%D0%9C%D0%9D%D0%9A-%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC

#include <vector>

class opls
{
public:

    std::vector<std::vector<double>> P; /// values of orthogonal polynom in mesh nodes
    std::vector<std::vector<double>> c; /// coefficients of polynom from meshed nodes /// bottom diagonal matrix
    std::vector<double> Q;              /// sum Pi^2
    std::vector<double> S;              /// sum xi*Pi^2
    int polynom_degree = 0;

    opls();
    opls(int n_polynom_degree, int data_size); /// init size of P, Q, S


    ///METHODS
    void QSP_count(const std::vector<double> &x); ///count orthogonal polynom values P[] in mesh nodes
    void orthogonal_polynom_coefficients_count(); ///count coefficient c[] of the orthogonal polynom
    void regression_polynom_coefficient_count(const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a);

    void manage_opls(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a);



    ///EXTRA
    void info();
    void print_QSP();
    void print_c();

    void Test();
};


#endif // OPLS_H_INCLUDED
