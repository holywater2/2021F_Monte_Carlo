// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = constant edge length of lattice
#include <math.h>
#include <stdlib.h>
#include <string.h>
// drandom  = dobule precision random in the range 0 <= r < 1

#define L 5
#define N (L*L)
#define XNN 1
#define YNN L

int s[N];
double prob[5];
double beta;

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

        i = N*rand();

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

int main(){
    initialize();
    sweep();   
}