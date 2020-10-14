#ifndef OPLS_H_INCLUDED
#define OPLS_H_INCLUDED

///REFERENCE https://ru.stackoverflow.com/questions/383716/%D0%9C%D0%9D%D0%9A-%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC

#include <vector>

class opls
{
public:

    std::vector<std::vector<double>> P;

    ///METHODS
    void get_QS(int j, const std::vector<double> &x, std::vector<double> &Q, std::vector<double> &S);

    void manage_opls(int polynom_order, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a);

    void Test();


};

#endif // OPLS_H_INCLUDED
