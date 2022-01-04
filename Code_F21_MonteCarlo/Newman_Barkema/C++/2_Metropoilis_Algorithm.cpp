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
#include <fstream>
#include <ctime>
#include <vector>

using namespace std;

#define L 100
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 40

static int s[N];
static double prob[5]; // 1 1 exp 1 exp
static random_device rd;  // Will be used to obtain a seed for the random number engine
static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static uniform_real_distribution<> dis(0.0, 1.0);
static int Fliped_Step;
static int Total_Step;

void prob_calc(double beta){
    int i;
    for(i = 2; i < 5; i += 2){
        prob[i] = exp(-2*beta*i);
    }
}

void initialize(double beta){
    int i;
    for(i = 0; i < N; i++){
        s[i] = 1;
    }
    prob_calc(beta);
}

int sweep_helical(int i){
    int nn, sum = 0;

    if((nn = i - XNN) < 0) nn += N;
    sum += s[nn];
    if((nn = i + XNN) >= N) nn -= N;
    sum += s[nn];
    if((nn = i - YNN) < 0) nn += N;
    sum += s[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];

    return sum;
}

int helical(int i){
    int nn, sum;
    
    if((nn = i + XNN) == N) nn = 0;
    sum = s[nn];

    if((nn = i + YNN) >= N) nn -= N;
    sum += s[nn];

    return sum;
}

vector<int> measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = helical(i);
        res += J*sum*s[i];
        sigma += s[i];
    }
    
    HH = -res-B*sigma;

    return vector<int>({HH,sigma});
}

int calculate(int n = N){
    int i, k, delta;
    double a;
    for(i = 0; i < n; i++){
        k = N*dis(gen);

        delta = sweep_helical(k)*s[k];
        
        a = dis(gen);
        Total_Step++;

        if(delta <= 0){
            Fliped_Step++;
            s[k] = -s[k];
        } else if(a < prob[delta]){
            Fliped_Step++;
            s[k] = -s[k];
        }
    }
}

int main(){
    vector<int> value;
    vector<double> m(bin);
    vector<double> c(bin);
    vector<vector<double>> res(bin,vector<double>(4,0));

    int sigma, HH;

    for(int i = 0; i < bin; i++){
        Fliped_Step = 0;
        Total_Step = 0;
        double T = (5/double(bin))*(i+1);
        double beta = 1/T;

        initialize(beta);
        for(int j = 0; j < 2000; j++){
            calculate();
        }
        value = measure();
        HH = value[0];
        sigma = value[1];

        cout <<"idx: " << i << "||" << sigma << " " << HH << "\n";

        int epoch = 18000;
        for(int j = 0; j < epoch; j++){
            calculate();
            
            value = measure();
            HH = value[0];
            sigma = value[1];
           
            res[i][0] += abs(sigma)/double(epoch);
            res[i][1] += (sigma*sigma)/double(epoch);
            res[i][2] += HH/double(epoch);
            res[i][3] += (HH*HH)/double(epoch);
        }
        m[i] = res[i][0]/N;
        c[i] = (beta*beta)/N*(res[i][3]-(res[i][2]*res[i][2]));
        cout << m[i] << " " << c[i] << " " << Fliped_Step << " " << Total_Step << '\n';
    }
    
    bool save = true;

    if (save){
        ofstream myfile;
        myfile.open("Metropolis_c_"+to_string(L)+".csv");
        myfile << ",temperture,magnetization,specific heat,abs(sigma),sigma**2,HH,HH**2\n";
        for(int i = 0; i <bin; i++){
            string temp = to_string(i) + "," + to_string((5/double(bin))*(i+1)) + "," + to_string(m[i]) + "," + to_string(c[i]) + ",";
            temp = temp + to_string(res[i][0]) + "," + to_string(res[i][1]) + "," + to_string(res[i][2]) + "," + to_string(res[i][3]) + "\n";
            myfile << temp;
        }
        myfile.close();
    }
}