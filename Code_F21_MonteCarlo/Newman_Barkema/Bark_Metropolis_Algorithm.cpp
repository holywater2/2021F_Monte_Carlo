// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = constant edge length of lattice
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
// drandom  = dobule precision random in the range 0 <= r < 1

using namespace std;

#define L 5
#define N (L*L)
#define XNN 1
#define YNN L

int s[N];
double prob[5];
double beta = 1;

void initialize(){
    int i;
    for(i=2; i<5; i+=2){
        prob[i] = exp(-2*beta*i);
    }
}

void sweep(){
    int i, k;
    int nn, sum, delta;

    for(k=0; k<N; k++){
        /* Choose a site */

        i = N*(rand()/(double)RAND_MAX);

        /* Calculate the sum of the neighbouring spins */
        if((nn = i + XNN) >= N) nn -= N;
        sum = s[nn];
        if((nn = i - XNN) < 0) nn += N;
        sum += s[nn];
        if((nn = i + YNN) >= N) nn -= N;
        sum += s[nn];
        if((nn = i - XNN) < 0) nn += N;
        sum += s[nn];

        /* Calculate the change in energy */

        delta = sum*s[i];

        /* Decide whether to flip spin */

        if(delta <= 0){
            s[i] = -s[i];
        } else if (rand()<prob[delta]) {
            s[i] = -s[i];
        }
    }
}

void print_s(int s[]){
    string res = "";
    for(int i = 0; i < N; i++){
        res += char(s[i]+'0');
        res += " ";
        if(i%L == L-1){
            res.pop_back();
            res += "\n";
        }
    }
    cout << res << endl;
}

int main(){
    initialize();
    print_s(s);
    sweep();
    print_s(s);
}