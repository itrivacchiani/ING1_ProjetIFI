//
//  TrinomialTree.cpp
//  C++ IFI
//
//  Created by Itri-Raphaël Vacchiani on 21/05/2016.
//  Copyright © 2016 Itri-Raphaël Vacchiani. All rights reserved.
//

#include "TrinomialTree.hpp"

using namespace std;


void TrinomTree::setValeur(double _T, int _N, double _S0, double _K, double _sig, double _r){
    
    T=_T;
    N=_N;
    S0=_S0;
    K=_K;
    sig=_sig;
    r=_r;
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// FONCTION MAX/MIN /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double TrinomTree::max(double a, double b){
    
    return a>=b?a:b;
}

double TrinomTree::min(double a, double b){
    
    return a<=b?a:b;
}

double TrinomTree::callPayoff(double St, double K, TrinomTree & TT) {
    
    return TT.max((St-K),0);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// CALCUL ACTIF SS-JACENT //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<vector<double>> TrinomTree::calculActifTri(int N, double S0, double u, double d){
    
    vector<vector<double>> S(N+1,vector<double>(2*N));
    
    for(int n=0;n<=N;n++)
    {
        for (int i=0;i<=2*n;i++)
        {
            S[n][i]=(S0*pow(u,(i/2))*pow(d,n-(i/2)));
        }
    }
    
    return S;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// CALCUL VOLATILITE  //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<vector<double>> TrinomTree::volatiliteSmile(int N, double S0, double u, double d, double sig0,double sig, TrinomTree & TT){
    
    vector<vector<double>> volat_smile(N+1,vector<double>(2*N));
    
    vector<vector<double>> Sn(N+1,vector<double>(2*N));
    
    Sn = TT.calculActifTri(N, S0, u, d);
    
    for(int n=0;n<=N;n++)
        
    {
        for (int i=0;i<=2*n;i++)
        {
            volat_smile[n][i]= TT.min(sig,(sig0/sqrt(Sn[n][i])));
        }
    }
    
    return volat_smile;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// COEFFICIENTS OPT EURO //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double TrinomTree::up(double T, int N, double sig, double r,TrinomTree & TT) const {
    
    double delta_t=(T/N);
    double d = TT.down(T, N, sig, r);
    return (pow(1+r*delta_t,2))/d;
}

double TrinomTree::down(double T, int N, double sig, double r) const {
    
    double delta_t=(T/N);
    return 1+(r*delta_t)-(sig*sqrt(delta_t));
}

vector<vector<double>> TrinomTree::prob_p(double T, int N, vector<vector<double>> sigma, double r, TrinomTree & TT, double u, double d) const {
    
    double delta_t=(T/N);
    vector<vector<double>> p(N+1,vector<double>(2*N));
    
    for (int n=0;n<=N;n++)
    {
        for (int i=0;i<=2*n;i++)
        {
            p[n][i]=(pow(sigma[n][i],2)*delta_t)/((u-d)*(u-1-(r*delta_t)));
        }
    }
    
    return p;
}

vector<vector<double>> TrinomTree::prob_q(double T, int N, vector<vector<double>> sigma, double r, TrinomTree & TT, double u, double d) const {
    
    double delta_t=(T/N);
    vector<vector<double>> q(N+1,vector<double>(2*N));
    
    for (int n=0;n<=N;n++)
    {
        for (int i=0;i<=2*n;i++)
        {
            q[n][i]=(pow(sigma[n][i],2)*delta_t)/((u-d)*(1+(r*delta_t)-d));
        }
    }
    
    return q;
}


double TrinomTree::optCETri(double S0, double u, double d, double T, double r, vector<vector<double>> p, vector<vector<double>> q, double K, double N, TrinomTree & TT){
    
    vector<vector<double>>CE(N+1,vector<double>(2*N));
    vector<vector<double>>S(N+1,vector<double>(2*N));
    
    S=TT.calculActifTri(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (int i=0;i<=2*N;i++)
    {
        CE[N][i]=TT.callPayoff(S[N][i],K,TT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        for (int i=0;i<=2*n;i++)
        {
            CE[n][i]=((p[n][i]*CE[n+1][i+2])+(q[n][i]*CE[n+1][i])+((1-p[n][i]-q[n][i])*CE[n+1][i+1]))/R;
        }
    }
    
    return CE[0][0];
}



void TrinomTree::optEuropSVarTri(double M, double u, double d, double T, double r, vector<vector<double>> p,vector<vector<double>> q, double K, double N, TrinomTree & TT){
    
    vector<double>CE_Svar(M+1);
    
    ofstream fichier1("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/CallEurop_Svar_Trinomial.txt", ios::out);
    
    for (double i=0;i<=M;i++)
    {
        CE_Svar[i]=TT.optCETri(i, u, d, T, r, p, q, K, N, TT);
        fichier1 << CE_Svar[i] <<endl;
        
            }
    fichier1.close();
    
}




