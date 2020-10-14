
#include "opls.h"

#include <iostream>
#include <vector>

void opls::get_QS(int j, const std::vector<double> &x, std::vector<double> &Q, std::vector<double> &S)
{

};

void opls::manage_opls(int polynom_order, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a)
{
    std::vector<double> Q(x.size(), 0);
    std::vector<double> S(x.size(), 0);

    for(int j = 0; j < x.size(); ++j){
        get_QS(j, x, Q, S);
            std::cout << "j = " << j << " Q = " << Q[j] << " S = " << S[j] << std::endl;
    }
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
