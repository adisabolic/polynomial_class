#ifndef POLINOM_H
#define POLINOM_H

#include <vector>
#include <fstream>
using namespace std;

class Polinom{

    vector<double> koeficijenti;

public:

//----------------------------------------- KONSTRUKTORI -----------------------------------------
    Polinom();
    Polinom(vector<double> k);

//----------------------------------------- OPERATORI -----------------------------------------

    friend istream& operator>>(istream&, Polinom&);
    friend ostream& operator<<(ostream&, const Polinom&);
    double operator()(double);
    friend Polinom operator-(const Polinom&);
    friend Polinom operator+(const Polinom&, const Polinom&);
    friend Polinom operator-(const Polinom&, const Polinom&);
    friend Polinom operator*(const Polinom&, const Polinom&);

//----------------------------------------- DODATNE FUNKCIJE -----------------------------------------
    void NapraviPolinom(vector<double> k);
    void PomnoziMonomom(int m);
    double NulaPolinoma(double,double,double);
    friend void SortirajPolinome(vector<Polinom>&);

};

//----------------------------------------- DODATNE FUNKCIJE VAN KLASE-----------------------------------------
Polinom operator^(const Polinom&, int n);
void SortirajPolinome(vector<Polinom>&);
double PresjekPolinoma(const Polinom&, const Polinom&, double, double, double);

#endif // POLINOM_H
