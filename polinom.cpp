#ifndef POLINOM_CPP
#define POLINOM_CPP
#include <iostream>
#include "polinom.h"
#include <cmath>
#include <algorithm>

//----------------------------------------- KONSTRUKTORI -----------------------------------------

Polinom::Polinom(){
    koeficijenti = {0};
}

Polinom::Polinom(vector<double> k){
    int i = k.size() - 1;
    int pravaVelicina(k.size());
    if(k.size() == 0)
        k.push_back(0);

    while(k[i] == 0 && i >= 0){
        pravaVelicina--;
        i--;
    }
    if(pravaVelicina == 0)
        pravaVelicina++;
    for(i = 0 ; i < pravaVelicina ; i++)
        koeficijenti.push_back(k[i]);
}

//----------------------------------------- OPERATORI -----------------------------------------

istream& operator>>(istream &fin, Polinom &p){
    double k;
    char znak(' ');
    vector<pair<double,int>> monomi;
    int stepen(-1);

    while(znak != '\n'){
        znak = fin.peek();
        if(znak == 'x'){
            k = 1;
        }
        else if(znak >= '0' && znak <= '9'){
            fin>>k;
            if(fin.peek() == '*'){
                fin>>znak;
                if(fin.peek() != 'x')
                    throw("Neispravan unos!");
            }
        }
        else if(znak == '+'){
            fin>>znak;
            znak = fin.peek();
            if(znak == 'x')
                k = 1;
            else if(znak >= '0' && znak <= '9'){
                fin>>k;
                if(fin.peek() == '*'){
                    fin>>znak;
                    if(fin.peek() != 'x')
                        throw("Neispravan unos!");
                }
            }
            else
                throw("Neispravan unos!");
        }
        else if(znak == '-'){
            fin>>znak;
            znak = fin.peek();
            if(znak == 'x')
                k = -1;
            else if(znak >= '0' && znak <= '9'){
                fin>>k;
                k = -k;
                if(fin.peek() == '*'){
                    fin>>znak;
                    if(fin.peek() != 'x')
                        throw("Neispravan unos!");
                }
            }
            else
                throw("Neispravan unos!");
        }
        else
            throw("Neispravan unos!");

        znak = fin.peek();
        if(znak == 'x'){
            fin>>znak;
            znak = fin.peek();
            if(znak == '^'){
                fin>>znak;
                znak = fin.peek();
                if(!(znak >= '0' && znak <= '9'))
                    throw("Neispravan unos!");
                else
                    fin>>stepen;
            }
            else
                stepen = 1;
        }
        else
            stepen = 0;

        monomi.push_back({k,stepen});
        znak = fin.peek();
    }

    fin.ignore(1000,'\n');
    int max_stepen(-1);

    for(int i = 0 ; i < monomi.size() ; i++)
        if(monomi[i].second > max_stepen)
            max_stepen = monomi[i].second;

    vector<double> sviKoef(max_stepen+1,0);

    for(int i = 0 ; i < monomi.size() ; i++)
        sviKoef[monomi[i].second] += monomi[i].first;

    p.NapraviPolinom(sviKoef);
    return fin;
}

ostream& operator<<(ostream &fout, const Polinom &p){
    if(p.koeficijenti.size() == 1)
        fout<<p.koeficijenti[0];
    else{
        for(int i = p.koeficijenti.size()-1 ; i >=0 ; i--){
            if(p.koeficijenti[i] != 0){
                if(i != p.koeficijenti.size()-1 && p.koeficijenti[i] > 0)
                    fout<<"+";

                if(p.koeficijenti[i] != 1 && p.koeficijenti[i] != -1)
                    fout<<p.koeficijenti[i];
                else if(p.koeficijenti[i] == -1)
                    fout<<"-";
                if(i == 0 && (p.koeficijenti[i] == 1 || p.koeficijenti[i] == -1))
                   fout<<"1";

                if(i != 0 && i != 1)
                    fout<<"x^"<<i;
                else if(i == 1)
                    fout<<"x";
            }
        }
    }
    return fout;
}

double Polinom::operator()(double x){
    double Px(0);
    for(int i = koeficijenti.size() - 1; i >= 0; i--)
        Px = Px * x + koeficijenti[i];
    return Px;
}

Polinom operator-(const Polinom &p){
    vector<double> k;
    for(int i = 0 ; i < p.koeficijenti.size() ; i++)
        k.push_back(-p.koeficijenti[i]);
    return k;
}

Polinom operator+(const Polinom &p, const Polinom &q){
    vector<double> k;
    if(p.koeficijenti.size() > q.koeficijenti.size()){
        k = p.koeficijenti;
        for(int i = 0 ; i < q.koeficijenti.size() ; i++)
            k[i] += q.koeficijenti[i];
    }
    else{
        k = q.koeficijenti;
        for(int i = 0 ; i < p.koeficijenti.size() ; i++)
            k[i] += p.koeficijenti[i];
    }
    return k;
}

Polinom operator-(const Polinom &p, const Polinom &q){
    return p+(-q);
}

Polinom operator*(const Polinom &p, const Polinom &q){
    int n1(p.koeficijenti.size());
    int n2(q.koeficijenti.size());

    if(n1 == 1 || n2 == 1){
        vector<double> T;
        if(n1 == 1){
            T = q.koeficijenti;
            for(int i = 0 ; i < T.size() ; i++)
                T[i] *= p.koeficijenti[0];
        }
        else{
            T = p.koeficijenti;
            for(int i = 0 ; i < T.size() ; i++)
                T[i] *= q.koeficijenti[0];
        }
        return T;
    }
    int n = max(n1,n2);
    int m = ceil(double(n)/2);

    Polinom A,B,C,D;
    vector<double> tempKoef1;
    vector<double> tempKoef2;
    for(int i = 0 ; i < p.koeficijenti.size() ; i++){
        if(i < m)
            tempKoef1.push_back(p.koeficijenti[i]);
        else
            tempKoef2.push_back(p.koeficijenti[i]);
    }
    A.NapraviPolinom(tempKoef1);
    B.NapraviPolinom(tempKoef2);

    tempKoef1 = {};
    tempKoef2 ={};

    for(int i = 0 ; i < q.koeficijenti.size() ; i++){
        if(i < m)
            tempKoef1.push_back(q.koeficijenti[i]);
        else
            tempKoef2.push_back(q.koeficijenti[i]);
    }

    C.NapraviPolinom(tempKoef1);
    D.NapraviPolinom(tempKoef2);

    Polinom T1(A*C),T2(B*D),T3((A+B)*(C+D)-T1-T2);
    T3.PomnoziMonomom(m);
    Polinom T4(T2);
    T4.PomnoziMonomom(2*m);

    return T1 + T3 + T4;
}

Polinom operator^(const Polinom &p, int n){
    if(n < 0)
        throw("Stepen mora biti nenegativan!");
    if(n == 0)
        return vector<double>{1};
    if(n == 1)
        return p;
    Polinom p1(p^(n/2));
    if(n % 2 == 0)
        return p1*p1;
    else
        return p1*p1*p;
}

//----------------------------------------- DODATNE FUNKCIJE -----------------------------------------


void Polinom::NapraviPolinom(vector<double> k){
    koeficijenti = {};
    int i = k.size() - 1;
    int pravaVelicina(k.size());
    if(k.size() == 0)
        k.push_back(0);

    while(k[i] == 0 && i >= 0){
        pravaVelicina--;
        i--;
    }
    if(pravaVelicina == 0)
        pravaVelicina++;
    for(i = 0 ; i < pravaVelicina ; i++)
        koeficijenti.push_back(k[i]);
}

void Polinom::PomnoziMonomom(int m){
    koeficijenti.resize(koeficijenti.size() + m);
    for(int i = koeficijenti.size() - 1 ; i >= 0 ; i--)
        if(i >= m)
            koeficijenti[i] = koeficijenti[i-m];
        else
            koeficijenti[i] = 0;
}

double Polinom::NulaPolinoma(double a,double b,double e){
    if(this->operator()(a) == 0)
        return a;
    if(this->operator()(b) == 0)
        return b;

    if(this->operator()(a) * this->operator()(b) > 0)
        throw("Vrijednosti polinoma u tackama a i b moraju biti suprotnog znaka!");
    if(b < a)
        swap(a,b);

    double c((a+b)/2);
    while(abs(b-c) > e){
        if((this->operator()(a) > 0 && this->operator()(c) < 0) || (this->operator()(a) < 0 && this->operator()(c) > 0))
            b = c;
        else if((this->operator()(b) > 0 && this->operator()(c) < 0) || (this->operator()(b) < 0 && this->operator()(c) > 0))
            a = c;
        else
            return c;
        c = (a+b)/2;
    }
    return c;
}

double PresjekPolinoma(const Polinom &p, const Polinom &q, double a, double b, double e){
    try{
        return (p-q).NulaPolinoma(a,b,e);
    }
    catch(...){
        throw ("Vrijednost (P(a)-P(b))*(Q(a)-Q(b)) mora biti manja od 0!");
    }
}

void SortirajPolinome (vector<Polinom> &v){
    sort(v.begin(),v.end(),[](const Polinom &p, const Polinom &q){
         return (p.koeficijenti.size() > q.koeficijenti.size() ||
                (p.koeficijenti.size() == q.koeficijenti.size()
                 && abs(p.koeficijenti[p.koeficijenti.size()-1]) > abs(q.koeficijenti[q.koeficijenti.size()-1])));});
}

#endif // POLINOM_CPP
