//
//  TrinomialTree.hpp
//  C++ IFI
//
//  Created by Itri-Raphaël Vacchiani on 21/05/2016.
//  Copyright © 2016 Itri-Raphaël Vacchiani. All rights reserved.
//

#ifndef TrinomialTree_hpp
#define TrinomialTree_hpp

#include <stdio.h>
#include<vector>
#include<cmath>
#include<fstream>
#include <string>
#include"BinomialTree.h"
using namespace std;

class TrinomTree {

protected:
    
    double T;
    double K;
    int N;
    double S0;
    double sig0;
    double sig;
    double r;


public:
    
    void setValeur(double T, int N, double S0, double K, double sig0, double r);
    
    double max(double a, double b);
    double min(double a, double b);
    
    
    double callPayoff(double St, double K, TrinomTree & TT);

    vector<vector<double>> calculActifTri(int N, double S0, double u, double d);
    
    vector<vector<double>> volatiliteSmile(int N, double S0, double u, double d, double sig0, double sig, TrinomTree & TT);
    
    double up(double T, int N, double sig, double r,TrinomTree & TT) const;
    
    double down(double T, int N, double sig, double r) const;
    
    vector<vector<double>> prob_p(double T, int N, vector<vector<double>> sigma, double r, TrinomTree & TT, double u, double d) const;
    
    vector<vector<double>> prob_q(double T, int N, vector<vector<double>> sigma, double r, TrinomTree & TT, double u, double d) const;
    
    double optCETri(double S0, double u, double d, double T, double r, vector<vector<double>> p,vector<vector<double>> q, double K, double N, TrinomTree & TT);
    
    
    void optEuropSVarTri(double M, double u, double d, double T, double r, vector<vector<double>> p,vector<vector<double>> q, double K, double N, TrinomTree & TT);

};



#endif /* TrinomialTree_hpp */
