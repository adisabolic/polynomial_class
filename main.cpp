#include <iostream>
#include "polinom.h"
#include <vector>

using namespace std;

int main(){

    Polinom p,q;
    try{

        cout<<"Unesi P(X): ";
        cin>>p;
        cout<<"Unesi Q(X): ";
        cin>>q;
        cout<<"P = "<<p<<endl;
        cout<<"Q = "<<q<<endl;
        cout<<"P + Q = "<<p+q<<endl;
        cout<<"P - Q = "<<p-q<<endl;
        cout<<"P * Q = "<<p*q<<endl;
        cout<<"P^2 = "<<(p^2)<<endl;
        cout<<"Q^3 = "<<(q^3)<<endl;

        vector<Polinom> v{p,q};
        SortirajPolinome(v);
        cout<<"Sortirani polinomi : "<<endl;
        for(int i = 0 ; i < v.size() ; i++)
            cout<<i+1<<".:"<<v[i]<<endl;

        cout<<"P(0) = "<<p(0)<<endl;
        cout<<"Q(-1) = "<<q(-1)<<endl;
        cout<<"Nula polinoma P: "<<p.NulaPolinoma(-10,10,0.000000001)<<endl;
        cout<<"Nula polinoma Q: "<<q.NulaPolinoma(-10,10,0.000000001)<<endl;
        cout<<"Presjek polinoma: "<<PresjekPolinoma(p,q,-10,10,0.000000001)<<endl;
    }
    catch(const char* poruka){
        cout<<poruka;
    }

}
