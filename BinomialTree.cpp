//
//  BinomialTree.cpp
//  C++ IFI
//
//  Created by Itri-Raphaël Vacchiani on 16/04/2016.
//  Copyright © 2016 Itri-Raphaël Vacchiani. All rights reserved.
//

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include <string>
#include"BinomialTree.h"
using namespace std;

void BinomTree::setValeur(double _T, int _N, double _S0, double _K, double _sig, double _r){
    
    T=_T;
    N=_N;
    S0=_S0;
    K=_K;
    sig=_sig;
    r=_r;
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// FONCTION MAX /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BinomTree::max(double a, double b){
    
    return a>=b?a:b;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// COEFFICIENTS COX ROSS //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BinomTree::CRR_up(double T, int N, double sig) const {
    
    double delta_t=(T/N);
    return exp(sig*sqrt(delta_t));
}

double BinomTree::CRR_down(double T, int N, double sig) const {
    
    double delta_t=(T/N);
    return exp(-sig*sqrt(delta_t));
}

double BinomTree::CRR_p(double T, int N, double sig, double r, BinomTree & BT) const {
    
    double delta_t=(T/N);
    return ((exp(r*delta_t) - BT.CRR_down(T,N,sig))/(BT.CRR_up(T,N,sig) - BT.CRR_down(T,N,sig)));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// CALCUL ACTIF SS-JACENT //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<vector<double>> BinomTree::calculActif(int N, double S0, double u, double d){
    
    vector<vector<double>> S(N+1,vector<double>(N+1));
    
    for(int n=0;n<=N;n++)
    {
        for (int i=0;i<=n;i++)
        {
            S[n][i]=(S0*pow(u,i)*pow(d,n-i));
        }
    }
    
    return S;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// FONCTIONS PAYOFF //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BinomTree::callPayoff(double St, double K, BinomTree & BT) {
    
    return BT.max((St-K),0);
}

double BinomTree::putPayoff(double St, double K, BinomTree & BT) {
    
    return BT.max((K-St),0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// OPTIONS EUROPEENES //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BinomTree::optCE(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<vector<double>>CE(N+1,vector<double>(N+1));
    vector<vector<double>>S(N+1,vector<double>(N+1));

    S=BT.calculActif(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);

    for (int i=0;i<=N;i++)
    {
        CE[N][i]=BT.callPayoff(S[N][i],K,BT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        for (int i=0;i<=n;i++)
        {
            CE[n][i]=((p*CE[n+1][i+1])+((1-p)*CE[n+1][i]))/R;
        }
    }
    
    return CE[0][0];
}

double BinomTree::optPE(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<vector<double>>PE(N+1,vector<double>(N+1));
    vector<vector<double>>S(N+1,vector<double>(N+1));
    
    S=BT.calculActif(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (int i=0;i<=N;i++)
    {
        PE[N][i]=BT.putPayoff(S[N][i],K,BT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        for (int i=0;i<=n;i++)
        {
            PE[n][i]=((p*PE[n+1][i+1])+((1-p)*PE[n+1][i]))/R;
        }
    }
    
    return PE[0][0];
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////OPTIONS AMERICAINES //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BinomTree::optCA(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<vector<double>>CA(N+1,vector<double>(N+1));
    vector<vector<double>>S(N+1,vector<double>(N+1));
    
    S=BT.calculActif(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (int i=0;i<=N;i++)
    {
        CA[N][i]=BT.callPayoff(S[N][i],K,BT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        for (int i=0;i<=n;i++)
        {
            CA[n][i]=BT.max(callPayoff(S[n][i], K, BT),((p*CA[n+1][i+1])+((1-p)*CA[n+1][i]))/R);
        }
    }
    
    return CA[0][0];
}

double BinomTree::optPA(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<vector<double>>PA(N+1,vector<double>(N+1));
    vector<vector<double>>S(N+1,vector<double>(N+1));
    
    S=BT.calculActif(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (int i=0;i<=N;i++)
    {
        PA[N][i]=BT.putPayoff(S[N][i],K,BT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        for (int i=0;i<=n;i++)
        {
            PA[n][i]=BT.max(putPayoff(S[n][i], K, BT),((p*PA[n+1][i+1])+((1-p)*PA[n+1][i]))/R);
        }
    }
    
    return PA[0][0];
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////OPTIONS EUROP S VAR + FILES //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void BinomTree::optEuropSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<double>CE_Svar(M+1);
    vector<double>PE_Svar(M+1);
    
    ofstream fichier1("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/CallEurop_Svar.txt", ios::out);
    ofstream fichier2("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/PutEurop_Svar.txt", ios::out);

    
    for (double i=0;i<=M;i++)
    {
        CE_Svar[i]=BT.optCE(i, u, d, T, r, p, sig, K, N, BT);
        fichier1 << CE_Svar[i] <<endl;
        
        PE_Svar[i]=BT.optPE(i, u, d, T, r, p, sig, K, N, BT);
        fichier2 << PE_Svar[i] <<endl;
    }
    fichier1.close();
    fichier2.close();
    
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////OPTIONS AMERI S VAR + FILES //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   

void BinomTree::optAmeSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<double>CA_Svar(M+1);
    vector<double>PA_Svar(M+1);
    
    ofstream fichier1("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/CallAme_Svar.txt", ios::out);
    ofstream fichier2("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/PutAme_Svar.txt", ios::out);
    
    
    for (double i=0;i<=M;i++)
    {
        CA_Svar[i]=BT.optCA(i, u, d, T, r, p, sig, K, N, BT);
        fichier1 << CA_Svar[i] <<endl;
        
        PA_Svar[i]=BT.optPA(i, u, d, T, r, p, sig, K, N, BT);
        fichier2 << PA_Svar[i] <<endl;
    }
    fichier1.close();
    fichier2.close();
    
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// CONDOR PAYOFF + FILES //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void BinomTree::optCondor(double L, double K1, double K2, double K3, double K4, double u, double d, double T, double r, double p, double sig, double N, BinomTree & BT)
{
    vector<double>CoP(L+1);
    vector<double>CoOption(L+1);
    
    ofstream fichier1("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/CondorPayoff.txt", ios::out);
    ofstream fichier2("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/CondorOption.txt", ios::out);
    
    for (double i=0;i<=L;i++)
    {
        CoP[i]=BT.max((i-K1),0)+BT.max((i-K4),0)-BT.max((i-K2),0)-BT.max((i-K3),0);
        fichier1 << CoP[i] <<endl;
        
        CoOption[i]=BT.optCE(i, u, d, T, r, p, sig, K1, N, BT)+BT.optCE(i, u, d, T, r, p, sig, K4, N, BT)-BT.optCE(i, u, d, T, r, p, sig, K2, N, BT)-BT.optCE(i, u, d, T, r, p, sig, K3, N, BT);
        fichier2 << CoOption[i] <<endl;
        
    }
    fichier1.close();
    fichier2.close();
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// SIMPLE BARRIERE ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int BinomTree::boolBarrS(double S, double b_sup)
{
    if(S>=b_sup)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

void BinomTree::optBarrKO(double M, double u, double d, double T, double r, double p, double sig, double K, double N, double b_sup, BinomTree & BT)
{
    
    vector<double>Barr_S(M+1);
    vector<vector<double>> S(N+1,vector<double>(N+1));
    
    ofstream fichier("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/SimpleKnockOut.txt");
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (double St=0;St<=M;St++)
    {
        S = BT.calculActif(N,St,u,d);
        
        vector<vector<double>>CEBarr(N+1,vector<double>(N+1));
        
        for (double i=0;i<=N;i++)
        {
            CEBarr[N][i]=BT.boolBarrS(S[N][i], b_sup)*BT.callPayoff(S[N][i], K, BT);
        }

        for (double n=N-1;n>=0;n--)
        {
            for (double i=0;i<=n;i++)
            {
                CEBarr[n][i]=BT.boolBarrS(S[n][i], b_sup)*(((p*CEBarr[n+1][i+1])+((1-p)*CEBarr[n+1][i]))/R);
            }
        }
        
        fichier << CEBarr[0][0] <<endl;
    }
    
    fichier.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// DOUBLE KNOCK-OUT ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int BinomTree::boolBarrD(double S, double b_sup, double b_inf)
{
    if((S>=b_sup)or(S<=b_inf))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

void BinomTree::optBarrKOD(double M, double u, double d, double T, double r, double p, double sig, double K, double N, double b_sup, double b_inf, BinomTree & BT)
{
    
    vector<double>Barr_D(M+1);
    vector<vector<double>> S(N+1,vector<double>(N+1));
    
    ofstream fichier("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/DoubleKnockOut.txt");
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (double St=0;St<=M;St++)
    {
        S = BT.calculActif(N,St,u,d);
        
        vector<vector<double>>CEBarr(N+1,vector<double>(N+1));
        
        for (double i=0;i<=N;i++)
        {
            CEBarr[N][i]=BT.boolBarrD(S[N][i], b_sup, b_inf)*BT.callPayoff(S[N][i], K, BT);
        }
        
        for (double n=N-1;n>=0;n--)
        {
            for (double i=0;i<=n;i++)
            {
                CEBarr[n][i]=BT.boolBarrD(S[n][i], b_sup,b_inf)*(((p*CEBarr[n+1][i+1])+((1-p)*CEBarr[n+1][i]))/R);
            }
        }
        
        Barr_D[St]=CEBarr[0][0];
        fichier << Barr_D[St] <<endl;
    }
    
    fichier.close();
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// OPT BERMUDES ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double BinomTree::bermudes(double S0, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<vector<double>>berm(N+1,vector<double>(N+1));
    vector<vector<double>>S(N+1,vector<double>(N+1));
    
    S=BT.calculActif(N,S0,u,d);
    
    double delta_t=T/N;
    double R=exp(r*delta_t);
    
    for (int i=0;i<=N;i++)
    {
        berm[N][i]=BT.putPayoff(S[N][i],K,BT);
    }
    
    for (int n=N-1;n>=0;n--)
    {
        
        if (n%10 != 0){
            
            for (int i=0;i<=n;i++)
            {
                berm[n][i]=((p*berm[n+1][i+1])+((1-p)*berm[n+1][i]))/R;
            }
        }
        else{
            for (int i=0;i<=n;i++)
            {
                berm[n][i]=BT.max(BT.putPayoff(S[n][i],K,BT),((p*berm[n+1][i+1])+((1-p)*berm[n+1][i]))/R);
            }
        }
            }
    
    return berm[0][0];

}


void BinomTree::bermudeSVar(double M, double u, double d, double T, double r, double p, double sig, double K, double N, BinomTree & BT){
    
    vector<double>berm_Svar(M+1);
    
    ofstream fichier1("/Users/itrivacchiani/Library/Mobile Documents/com~apple~CloudDocs/GM1/Projet Semestre 2/C++ IFI/C++/Bermudes_Svar.txt", ios::out);
        for (double i=0;i<=M;i++)
    {
        berm_Svar[i]=BT.bermudes(i, u, d, T, r, p, sig, K, N, BT);
        fichier1 << berm_Svar[i] <<endl;
    }
    fichier1.close();
 }



