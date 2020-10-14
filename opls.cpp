
#include "opls.h"

#include <iostream>
#include <vector>
#include <cmath>

opls::opls()
{
    polynom_degree = 0.0;
    P.resize(0);
    c.resize(0);
    Q.resize(0);
    S.resize(0);
}

opls::opls(int n_polynom_degree, int data_size){
    polynom_degree = n_polynom_degree;

    P.resize(polynom_degree+1);
        for(size_t i = 0; i < P.size(); ++i)
            P[i].resize(data_size);

    c.resize(polynom_degree+2);
    for(size_t i = 0; i < c.size(); ++i)
        for(size_t j = 0; j <= i; ++j)
            c[i].resize(j+1);

    Q.resize(polynom_degree+1);
    S.resize(polynom_degree+1);
}


void opls::QSP_count( const std::vector<double> &x){ ///j - height, i - width
    int j = 0;
    ///first row
    if(j == 0){
        for(size_t i = 0; i < P[0].size(); i++)
            P[j][i] = 1.0;
        for(auto _x : x)
            Q[j] += _x;
        S[j] = P[j].size(); ///summ of Pj^2
        j++;
    }

    while(j <= polynom_degree){
        for(size_t i = 0; i < P[0].size(); i++)
            P[j][i] = (x[i] - Q[j-1]/S[j-1])*P[j-1][i] - (j > 1 ? S[j-1]/S[j-2]*P[j-2][i] : 0.0);
        for(size_t i = 0; i < x.size(); ++i)
            Q[j] += x[i]*P[j][i]*P[j][i];
        for(size_t i = 0; i < P[j].size(); ++i)
            S[j] += P[j][i]*P[j][i];
        j++;
    }
};

void opls::orthogonal_polynom_coefficients_count(){
    c[0][0] = 1.0;
    c[1][0] = -Q[0]/S[0];   c[1][1] = 1.0;

    for(size_t k = 2; k < c.size(); k++)
    {
        for(size_t j = 0; j <= k; j++)
        {
            if(k != j){
                c[k][j] = (j == 0 ? 0.0 : c[k-1][j-1])
                    - Q[k-1]/S[k-1]*c[k-1][j]
                        - (j == k ? 0.0 : S[k-1]/S[k-2]*c[k-2][j]);
            }else{
                c[k][j] = 1.0;
            }
        }
    }
}

void opls::regression_polynom_coefficient_count(const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a)
{
    ///g(x) = SUM(b0*p0(x) + b1*p1(x) + ... + bm*pm(x))
    for(int k = 0; k <= polynom_degree; ++k)
        for(int i = 0; i < y.size(); ++i)
            b[k] += y[i]*P[k][i]/S[k];

    ///g(x) = SUM(a0 + a1*x + a2*x^2 + ... + am*x^m)
    for(int k = 0; k <= polynom_degree; ++k)
        for(int i = k; i <= polynom_degree; ++i)
            a[k] += b[i]*c[i][k];
}

void opls::manage_opls(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a){

    QSP_count(x);
    orthogonal_polynom_coefficients_count();
    regression_polynom_coefficient_count(y,b,a);
};

void opls::Test(){
    ///true result(fixed_input)
    std::vector<double> rS = {10, 69.17098425001, 328.7768235086, 962.75401849569};
    std::vector<double> rQ = {7.82083, -27.512890311149, -1.7129366948718, 250.39738237352};

    std::vector<std::vector<double>> rP = {
                                            {1,1,1,1,1,1,1,1,1,1},
                                            {-0.782083, 2.583797,2.855107,-0.217603,-3.809293,-4.617783,-1.899743,1.845867,3.175347,0.866387},
                                            {-7.2281734230875,2.8073623836192,4.6030924341892,-7.1264829728247,3.0992779145833,8.9585998726813,-5.5494580486592,-1.3320551385999,6.9121153510662,-5.1442783729677},
                                            {3.6796621840027,-2.8171823354639,3.1956859234941,-3.0255967835539,8.7399448097287,-12.367028020194,15.20316889478,-12.281110609376,12.297473127759,-12.625017191178},
                                            };

    std::vector<std::vector<double>> rc = {
                                            {1},
                                            {-0.782083,1},
                                            {-7.2281734230875,-0.38433110143359,1},
                                            {3.6796621840027,-11.983278954678,-0.37912107270775,1},
                                            {20.209167857339,7.9217601883617,-14.812965853531,-0.63920555697172,1}
                                            };

    std::vector<double> rb = {29.906462, 15.897655646369, 2.3791287087584, 1.000001267701};
    std::vector<double> ra = {3.9560877250835, 2.9999883433859, 2.0000071554385, 1.000001267701};


    ///TEST START
    ///DATA
        ///fixed_input
        std::vector<double> x = {0, 3.36588, 3.63719, 0.56448, -3.02721, -3.8357, -1.11766, 2.62795, 3.95743, 1.64847};
        std::vector<double> y = {3.9561, 74.84479, 89.44289, 6.46668, -14.53888, -34.55881, 1.70531, 43.80101, 109.1294, 18.81613};
        int polynom_degree = 3;

        ///oputput empty
        std::vector<double> b(polynom_degree + 1, 0);
        std::vector<double> a(polynom_degree + 1, 0);

    ///MANAGING
    opls O(polynom_degree, x.size());
        O.manage_opls(x, y, b, a);

    ///REPROT
    std::cout << "Difference S[]:\n";
    for(int i = 0; i < O.S.size(); ++i)
        std::cout << "\t" << i << " " << fabs(O.S[i]-rS[i]) << "\n";
            std::cout << std::endl;

    std::cout << "Difference Q[]:\n";
    for(int i = 0; i < O.Q.size(); ++i)
        std::cout << "\t" << i << " " << fabs(O.Q[i]-rQ[i]) << "\n";
            std::cout << std::endl;

    std::cout << "Difference P[][]:\n";
    for(int i = 0; i < O.P.size(); ++i){
        for(int j = 0; j < O.P[i].size(); ++j)
            std::cout << "\t" << i << "x" << j << " " << fabs(O.P[i][j]-rP[i][j]) << "\t";
                std::cout << std::endl;
        }
    std::cout << std::endl;

    std::cout << "Difference c[][]:\n";
    for(int i = 0; i < O.c.size(); ++i){
        for(int j = 0; j < O.c[i].size(); ++j)
            std::cout << "\t" << i << "x" << j << " " << fabs(O.c[i][j]-rc[i][j]) << "\t";
                std::cout << std::endl;
        }
    std::cout << std::endl;


    std::cout << "Difference b[]:\n";
    for(int i = 0; i < b.size(); ++i)
        std::cout << "\t" << i << " " << fabs(b[i]-rb[i]) << "\n";
            std::cout << std::endl;

    std::cout << "Difference a[]:\n";
    for(int i = 0; i < a.size(); ++i)
        std::cout << "\t" << i << " " << fabs(a[i]-ra[i]) << "\n";
            std::cout << std::endl;

    std::cout << "Coefficients b[]:\n";
    for(auto _b : b)
        std::cout << _b << "\t";
            std::cout << "\n";

    std::cout << "Coefficients a[]:\n";
    for(auto _a : a)
        std::cout << _a << "\t";
            std::cout << "\n";
}

void opls::info(){
    std::cout << "opls object\n";
        std::cout << "polynom degree = " << polynom_degree << std::endl;
        std::cout << "P size = " << P.size() << "x" << P[0].size() << std::endl;
        std::cout << "c size = " << c.size() << "x" << c[0].size() << ".." << c[c.size()-1].size() << std::endl;
        std::cout << "Q size = " << Q.size() << std::endl;
        std::cout << "S size = " << S.size() << std::endl;
}

void opls::print_QSP(){
    for(size_t j = 0; j < Q.size(); ++j)
        std::cout << "j = " << j << " Q = " << Q[j] << " S = " << S[j] << std::endl;

    for(size_t i = 0; i <  P.size(); ++i){
        for(size_t j = 0; j < P[i].size(); ++j){
            std::cout << P[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void opls::print_c(){
    for(size_t i = 0; i < c.size(); ++i){
        for(size_t j = 0; j < c[i].size(); ++j){
            std::cout << c[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}
