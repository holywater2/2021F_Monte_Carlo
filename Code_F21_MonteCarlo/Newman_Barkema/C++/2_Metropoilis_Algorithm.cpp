// * beta   = inverse temperature
// * prob[] = array of acceptance probability
// * s[]    = lattice of spins with helical boundary conditions
// * L      = const ant edge length of lattice
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <tuple>
#include <random>
using namespace std;

#define L 100
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 25

static int s[N];
static double prob[5]; // 1 1 exp 1 exp
// double beta = 1;
static random_device rd;  // Will be used to obtain a seed for the random number engine
static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static uniform_real_distribution<> dis(0.0, 1.0);

void prob_calc(double b){
    int i;
    for(i = 2; i < 5; i += 2){
        prob[i] = exp(-2*b*i);
    }
    // cout << prob[2] << " " << prob[4] << "\n";
}

void initialize(double beta){
    int i;
    for(i = 0; i < N; i++){
        s[i] = 1;
    }
    prob_calc(beta);
}

int sweep_helical(int i){
    int nn, sum;

    if((nn = i + XNN) >= N) nn -= N;
    sum = s[nn];
    if((nn = i - XNN) < 0) nn += N;
    sum += s[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];
    if((nn = i - XNN) < 0) nn += N;
    sum += s[nn];

    return sum;
}

int helical(int i){
    int nn, sum;
    if((nn = i + XNN) >= N) nn -= N;
    sum = s[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];

    return sum;
}


tuple<int,int> measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = helical(i);
        res += J*sum*s[i];
        sigma += s[i];
    }
    
    HH = -res-B*sigma;

    return tuple<int,int>(HH,sigma);
}

int calculate(int n = N){
    int i, k, delta;
    for(i = 0; i < n; i++){
        k = int(N*dis(gen));
        
        delta = 2*J*sweep_helical(k)*s[k];
        // cout << delta << " "  <<  k <<  " "<< sweep_helical(k) << " " << s[k] << "\n";
        
        if(delta <= 0){
            s[k] *= -1;
        } else if(dis(gen) < prob[int(delta/2)]){
            // cout << "fliped" << "\n";
            s[k] *= -1;
        }
    }
    // cout << k << " " << delta << " " << s[k] << "\n";
}

int main(){
   
    // for(int i = 0; i < 10; i++){
    //     cout << int(10*dis(gen)) << " ";
    // }

    tuple<int,int> value;
    double m[25];
    double c[25];
    double res[4];

    int sigma, HH;

    for(int i = 0; i < 25; i++){
        for(int j = 0; j < 4; j++){
            res[j] = 0;
        }
        double beta = 1/(0.2*(i+1));

        initialize(beta);
        for(int j = 0; j < 2000; j++){
            calculate();
        }
        value = measure();
        sigma = get<1>(value);
        HH = get<0>(value);
        cout << sigma << " " <<HH << "\n";

        int epoch = 18000;
        for(int j = 0; j < epoch; j++){
            calculate();
            value = measure();
            sigma = get<1>(value);
            HH = get<0>(value);

            res[0] += abs(sigma)/double(epoch);
            res[1] += sigma*sigma/double(epoch);
            res[2] += HH/double(epoch);
            res[3] += HH*HH/double(epoch);
        }
        m[i] = res[0]/N;
        c[i] = beta*beta*(res[3]-res[2]*res[2])/N;
        cout << m[i] << " " << c[i] << "\n";
    }    
} 