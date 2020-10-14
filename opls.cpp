
#include "opls.h"

#include <iostream>
#include <vector>

opls::opls(int polynom_degree, int data_size)
{
    P.resize(polynom_degree+1);
        for(int i = 0; i < P.size(); ++i)
            P[i].resize(data_size);

    Q.resize(data_size);
    S.resize(data_size);
}


void opls::QS_count(const std::vector<double> &x)
{
    for(int i = 0; i < P[0].size(); i++) ///first row
        P[0][i] = 1.0;

};

void opls::manage_opls(int polynom_order, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a)
{
    QS_count(x);
        print_QSP();

};

void opls::Test()
{
        std::vector<double> x = {0, 3.36588, 3.63719, 0.56448, -3.02721, -3.8357, -1.11766, 2.62795, 3.95743, 1.64847};
        std::vector<double> y = {3.9561, 74.84479, 89.44289, 6.46668, -14.53888, -34.55881, 1.70531, 43.80101, 109.1294, 18.81613};

        std::vector<double> rQ = {10,69.17098425001,328.7768235086,962.75401849569};
        std::vector<double> rS = {7.82083,-27.512890311149,-1.7129366948718,250.39738237352};

        std::vector<double> rb = {29.906462,15.897655646369, 2.3791287087584, 1.000001267701};
        std::vector<double> ra = {3.9560877250835, 2.9999883433859, 2.0000071554385, 1.000001267701};
}

void opls::print_QSP()
{
    for(int j = 0; j < Q.size(); ++j)
        std::cout << "j = " << j << " Q = " << Q[j] << " S = " << S[j] << std::endl;

    for(int i = 0; i <  P.size(); ++i){
        for(int j = 0; j < P[i].size(); ++j){
            std::cout << P[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}
