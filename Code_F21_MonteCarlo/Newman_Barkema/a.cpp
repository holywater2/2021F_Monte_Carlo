// Barkema Fig 1.1의 비밀을 알아냈다. Finite Lattice에서 계산을 한 결과를 이용하는 것이었다.


// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = const ant edge length of lattice
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;

#define L 3
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 25

int s[N];
int count;
int HH_c[N*N];
double Z[bin];
// double Z2[bin][N*N];
double ZU[bin];
// double ZU2[bin][N*N];

// Function for prints spin configuration
void print_s(){
    string res = "";
    for(int i = 0; i < N; i++){
        res += char(s[i]+'0');
        res += " ";
        if(i%L == L-1){
            res.pop_back();
            res += "\n";
        }
    }
    cout << res << '\n';
}

// Fucntion that makes all spins up configuration
void initialize(){
    int i;
    for(i = 0; i < N; i++)
        s[i] = 1;
}

// Function for calculate Hamiltonian
int sweep(){
    int i, k;
    int nn, sum = 0, sigma = 0;
    int HH = 0;
    for(i = 0; i < N; i++){
        sum = 0;
        if((nn = i + XNN)%L != 0)
            sum += s[nn];
        if(((nn = i - XNN)+1)%L != 0)
            sum += s[nn];
        if((nn = i + YNN) < N)
            sum += s[nn];
        if((nn = i - YNN) > 0)
            sum += s[nn];
        
        // cout << sum << " " << i << endl;
        sigma += s[i];
        HH -= J*sum*s[i];
    }
    // cout << endl;
    HH -= sigma*B;
    // print_s();
    // cout << sigma << " " << HH << endl << endl;
    for(k = 0; k < bin; k++){
        double Temp = 0.2*(k+1);
        ZU[k] += sigma*exp(-HH/Temp);
        Z[k] += exp(-HH/Temp);
        // if(count==93)
        //     cout << "&#& " << Z[k] << endl;
        // if(isnan(Z[k])|| isnan(ZU[k])) cout << endl << "error " << k << " " << count << " " << HH << " " << Temp <<endl;

        // cout << sigma*exp(-HH/Temp) << " " << exp(-HH/Temp) << endl;
    }
    return HH;
}

int calc(){
    // print_s();
    sweep();
    // HH_c[count] = sweep();
    count++;

}

void generator_sym(int i){
    if(i == N) return;
    if(i == 0) calc();
    s[i] = -1;
    calc();
    generator_sym(i+1);
    s[i] = 1;
    generator_sym(i+1);
}

void print_all(double* a, int len){
    string res = "";
    for(int i = 0; i < len; i++){
        res += to_string(a[i]) + " ";
    }
    res += "\n";
    cout << res << endl;
}

void print_all(int* a, int len){
    string res = "";
    for(int i = 0; i < len; i++){
        res += to_string(a[i]) + " ";
    }
    res += "\n";
    cout << res << endl;
}

int main(){
    initialize();
    generator_sym(0);
    // print_all(HH_c,N*N);
    print_all(ZU,25);
    print_all(Z,25);
    cout << count << endl;
    long double k = ZU[0]/Z[0];
    for(int i = 0; i < 25; i++){
        cout << (ZU[i]/Z[i])/k << " ";
    }
    cout << endl;
}