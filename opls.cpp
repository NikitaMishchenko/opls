
#include "opls.h"

#include <iostream>
#include <vector>

opls::opls(int n_polynom_degree, int data_size){
    polynom_degree = n_polynom_degree;

    P.resize(polynom_degree+1);
        for(size_t i = 0; i < P.size(); ++i)
            P[i].resize(data_size);

    c.resize(polynom_degree+2);
    for(size_t i = 0; i < c.size(); ++i){
        for(size_t j = 0; j <= i; ++j){
            c[i].resize(j+1);
        }
    }

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
        //std::cout << j << std::endl; print_QSP();std::cout << std::endl;
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
    { //std::cout << "k = " << k;
        for(size_t j = 0; j <= k; j++)
        {
            //std::cout << " j = " << j <<" ";
            if(k != j){
                c[k][j] = (j == 0 ? 0.0 : c[k-1][j-1])
                    - Q[k-1]/S[k-1]*c[k-1][j]
                        - (j == k ? 0.0 : S[k-1]/S[k-2]*c[k-2][j]);
            }else{
                c[k][j] = 1.0;
                //std::cout << " c = " << c[k][j] << std::endl;
            }
        }
    }
}

void opls::manage_opls(int polynom_order, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &b, std::vector<double> &a){
    QSP_count(x);
    info();
        //print_QSP();

    orthogonal_polynom_coefficients_count();
        print_c();

};

void opls::Test(){
        std::vector<double> x = {0, 3.36588, 3.63719, 0.56448, -3.02721, -3.8357, -1.11766, 2.62795, 3.95743, 1.64847};
        std::vector<double> y = {3.9561, 74.84479, 89.44289, 6.46668, -14.53888, -34.55881, 1.70531, 43.80101, 109.1294, 18.81613};

        std::vector<double> rQ = {10,69.17098425001,328.7768235086,962.75401849569};
        std::vector<double> rS = {7.82083,-27.512890311149,-1.7129366948718,250.39738237352};

        std::vector<double> rb = {29.906462,15.897655646369, 2.3791287087584, 1.000001267701};
        std::vector<double> ra = {3.9560877250835, 2.9999883433859, 2.0000071554385, 1.000001267701};
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
