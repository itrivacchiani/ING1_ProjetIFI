//
//  BinomialTree.h
//  C++ IFI
//
//  Created by Itri-Raphaël Vacchiani on 16/04/2016.
//  Copyright © 2016 Itri-Raphaël Vacchiani. All rights reserved.
//

#ifndef BINOMIALTREE_H_
#define BINOMIALTREE_H_
#include<iostream>
using namespace std;


class BinomTree {
    
protected:
    
    double T;
    double K;
    int N;
    double S0;
    double sig;
    double r;
    
public:
    
    void setValeur(double T, int N, double S0, double K, double sig, double r);
    
    double max(double a, double b);
    
    double CRR_up(double T, int N, double sig) const;
    
    double CRR_down(double T, int N, double sig) const;
    
    double CRR_p(double T, int N, double sig, double r, BinomTree & BT) const;
    
    vector<vector<double>> calculActif(int N, double S0, double u, double d);
    
    double callPayoff(double St, double K, BinomTree & BT);
    
    double putPayoff(double St, double K, BinomTree & BT);
    
    double optCE(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
 
    double optPE(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
    
    double optCA(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
    
    double optPA(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);

    void optEuropSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
    void optAmeSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
    
    void optCondor(double L, double K1, double K2, double K3, double K4, double u, double d, double T, double r, double p, double sig, double N, BinomTree & BT);
    
    void optBarrKO(double M, double u, double d, double T, double r, double p, double sig, double K, double N, double b_sup, BinomTree & BT);
    
    int boolBarrS(double S, double b_sup);
    
    void optBarrKOD(double M, double u, double d, double T, double r, double p, double sig, double K, double N, double b_sup, double b_inf, BinomTree & BT);
    
    int boolBarrD(double S, double b_sup, double b_inf);
    
    
    double bermudes(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
    
    void bermudeSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT);
};





#endif /* BinomialTree_h */
