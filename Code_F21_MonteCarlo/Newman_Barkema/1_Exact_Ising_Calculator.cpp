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
int HH_c[16];
double Z[bin];
double ZU[bin];


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

void initialize(){
    int i;
    for(i = 0; i < N; i++)
        s[i] = 1;
}

int sweep(){
    int i, k;
    int nn, sum, sigma = 0;
    int HH = 0;
    for(i = 0; i < N; i++){
        // Periodical Boundary condition
        if((nn = i + XNN)%L == 0) nn -= L;
        sum = s[nn];
        if(((nn = i - XNN)-1)%L == 0) nn += L;
        sum += s[nn];
        if((nn = i + YNN) >= N) nn -= N;
        sum += s[nn];
        if((nn = i - YNN) < 0) nn += N;
        sum += s[nn];

        // Helical Boundary condition

        // if((nn = i + XNN) >= N) nn -= N;
        // sum = s[nn];
        // if((nn = i - XNN) <  0) nn += N;
        // sum += s[nn];
        // if((nn = i + YNN) >= N) nn -= N;
        // sum += s[nn];
        // if((nn = i - YNN) <  0) nn += N;
        // sum += s[nn];
        
        cout << sum << " " << i << endl;
        sigma += s[i];
        HH -= J*sum*s[i];
    }
    
    HH -= sigma*B;
    print_s();
    cout << sigma << " " << HH << endl << endl;
    for(k = 0; k < bin; k++){
        double Temp = 0.2*(k+1);
        ZU[k] += sigma*exp(-HH/Temp);
        Z[k] += exp(-HH/Temp);

        cout << sigma*exp(-HH/Temp) << " " << exp(-HH/Temp) << endl;
    }
    return HH;
}

int calc(){
    // print_s();
    HH_c[count] = sweep();
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
    print_all(HH_c,16);
    print_all(ZU,25);
    print_all(Z,25);
    cout << count << endl;
    long double k = ZU[0]/Z[0];
    for(int i = 0; i < 25; i++){
        cout << (ZU[i]/Z[i])/k << " ";
    }
    cout << endl;
}