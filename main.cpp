//
//  main.cpp
//  C++ IFI
//
//  Created by Itri-Raphaël Vacchiani on 16/04/2016.
//  Copyright © 2016 Itri-Raphaël Vacchiani. All rights reserved.
//

#include<iostream>
#include<vector>
#include<cmath>
#include"BinomialTree.h"
#include "TrinomialTree.hpp"
using namespace std;

int main(){
    
    
    int choix;
    int _N;
    double _T;
    double _sig, _sig0;
    double _r;
    double _S0;
    double _K;
    double M;
    double L;
    double K1, K2, K3, K4;
    double u,d,p;
    double b_sup, b_inf;
    
    BinomTree BT;
    TrinomTree TT;
    
    BT.setValeur(_T, _N, _S0, _K, _sig, _r);
    
    TT.setValeur(_T, _N, _S0, _K, _sig0, _r);
    
    vector<vector<double>> sigma;
    vector<vector<double>>  p_tri;
    vector<vector<double>>  q_tri;
    
    do{
        cout<<"Les choix: "<<endl;
        cout<<"1: Options Européennes."<<endl;
        cout<<"2: Options Américaines."<<endl;
        cout<<"3: Option Condor."<<endl;
        cout<<"4: Simple Barrière Knock-Out."<<endl;
        cout<<"5: Double Knock-Out."<<endl;
        cout<<"6: Call Euro - Arbre Trinomial"<<endl;
        cout<<"7: Option Bermudes."<<endl;
        cout<<"0: Quitter"<<endl;
        cout<<"Entrer votre choix:"<<endl;
        cin>>choix;
        
        

        switch(choix)
        {
            case 1:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"S0 (valeur de l'option renvoyée dans le terminal) : ";
                cin>>_S0;
                cout<<endl;
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt) : ";
                cin>>M;
                cout<<endl;
                
               
                u=BT.CRR_up(_T, _N, _sig);
               
                d=1/u;
        
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.optEuropSVar(M, u, d, _T, _r, p, _sig,_K, _N, BT);

                double prix_ce;
                prix_ce=BT.optCE(_S0, u, d, _T, _r, p, _sig, _K, _N, BT);
                double prix_pe;
                prix_pe=BT.optPE(_S0, u, d, _T, _r, p, _sig, _K, _N, BT);
                cout<<"Prix du call"<<prix_ce<<"pour S0 = "<<_S0<<endl;
                cout<<"Prix du put"<<prix_pe<<"pour S0 = "<<_S0<<endl;
                
                
                
                break;
                
            case 2:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"S0 (valeur de l'option renvoyée dans le terminal) : ";
                cin>>_S0;
                cout<<endl;
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt): ";
                cin>>M;
                cout<<endl;
                
                u=BT.CRR_up(_T, _N, _sig);
                
                d=1/u;
             
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.optAmeSVar(M, u, d, _T, _r, p, _sig,_K, _N, BT);

                double prix_ca;
                prix_ca=BT.optCA(_S0, u, d, _T, _r, p, _sig, _K, _N, BT);
                double prix_pa;
                prix_pa=BT.optPA(_S0, u, d, _T, _r, p, _sig, _K, _N, BT);
                cout<<"Prix du call = "<<prix_ca<<"pour S0 = "<<_S0<<endl;
                cout<<"Prix du put = "<<prix_pa<<"pour S0 = "<<_S0<<endl;
                break;
                
                
                
            case 3:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"L : ";
                cin>>L;
                cout<<endl;
                
                cout<<"K1 : ";
                cin>>K1;
                cout<<endl;
                
                cout<<"K2  (>K1) :";
                cin>>K2;
                cout<<endl;
                
                cout<<"K3  (>K2) :";
                cin>>K3;
                cout<<endl;
                
                cout<<"K4  (>K3) :";
                cin>>K4;
                cout<<endl;
                
             
                u=BT.CRR_up(_T, _N, _sig);
               
                d=1/u;
               
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.optCondor(L, K1, K2, K3, K4, u, d, _T, _r, p, _sig, _N, BT);
                break;
                
            case 4:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt) : ";
                cin>>M;
                cout<<endl;
                
                cout<<"Barrière supérieure :";
                cin>>b_sup;
                cout<<endl;
                
                u=BT.CRR_up(_T, _N, _sig);
                
                d=1/u;
                
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.optBarrKO(M, u, d, _T, _r, p, _sig, _K, _N, b_sup, BT);
                break;
                
                
                
            case 5:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt) : ";
                cin>>M;
                cout<<endl;
                
                cout<<"Barrière supérieure :";
                cin>>b_sup;
                cout<<endl;
                
                cout<<"Barrière inférieure :";
                cin>>b_inf;
                cout<<endl;
                
                u=BT.CRR_up(_T, _N, _sig);
                
                d=1/u;
                
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.optBarrKOD(M, u, d, _T, _r, p, _sig, _K, _N, b_sup, b_inf, BT);
                
                break;

                
            case 6:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"S0 (valeur de l'option renvoyée dans le terminal) : ";
                cin>>_S0;
                cout<<endl;
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma 0 : ";
                cin>>_sig0;
                cout<<endl;
                
                cout<<"Sigma Chapeau : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt) : ";
                cin>>M;
                cout<<endl;
                
                
                d=TT.down(_T, _N, _sig, _r);
                
                u=TT.up(_T, _N, _sig, _r, TT);
                
                
                sigma=TT.volatiliteSmile(_N, _S0, u, d, _sig0, _sig, TT);
                
                p_tri = TT.prob_p(_T, _N, sigma, _r, TT, u, d);
                
                q_tri = TT.prob_q(_T, _N, sigma, _r, TT, u, d);
                
                
                TT.optEuropSVarTri(M, u, d, _T, _r, p_tri, q_tri, _K, _N, TT);
                
                
                double prix_cetri;
                prix_cetri=TT.optCETri(_S0, u, d, _T, _r, p_tri, q_tri, _K, _N, TT);
                
                cout<<"Prix du call"<<prix_cetri<<"pour S0 = "<<_S0<<endl;
                
                
                break;

                
            case 7:
                
                cout<<"T : ";
                cin>>_T;
                cout<<endl;
                
                cout<<"N : ";
                cin>>_N;
                cout<<endl;
                
                cout<<"S0 (valeur de l'option renvoyée dans le terminal) : ";
                cin>>_S0;
                cout<<endl;
                
                cout<<"Strike K : ";
                cin>>_K;
                cout<<endl;
                
                cout<<"Sigma : ";
                cin>>_sig;
                cout<<endl;
                
                cout<<"r : ";
                cin>>_r;
                cout<<endl;
                
                cout<<"Valeur max de S0 (valeurs de l'option enregistrées dans un fichier .txt) : ";
                cin>>M;
                cout<<endl;
                
                
                u=BT.CRR_up(_T, _N, _sig);
                
                d=1/u;
                
                p=BT.CRR_p(_T, _N, _sig, _r, BT);
                
                BT.bermudeSVar(M, u, d, _T, _r, p, _sig,_K, _N, BT);
                
                double prix_bermudes;
                prix_bermudes=BT.bermudes(_S0, u, d, _T, _r, p, _sig, _K, _N, BT);
                cout<<"Prix de l'option bermudes"<<prix_bermudes<<"pour S0 = "<<_S0<<endl;
                
                
                break;

                
            case 0:
                break;
                
            default:
                cout<<"Votre choix est incorrect"<<endl;
                break;
        }
        
                
    } while(choix!=0);


}

