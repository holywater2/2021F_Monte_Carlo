#ifndef ____METROPOLIS____
#define ____METROPOLIS____ 

// File & IO System
// #include <iostream>

// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
// Mathmatics
#include <math.h>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>

using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

typedef tuple<int,int> duo;

class Model{
    public:
        const int L;
        const int N;
        const int Bin;
        const int B;
        const int J;

        // #ifdef _WIN32
        // static string Filename = ".\\Result\\Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif
        // #ifdef linux
        // static string Filename = "./Result/Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        short* sc; // Square lattice configuration of 2D Ising model
        double prob[5];

        long Fliped_Step = 0;
        long Total_Step  = 0;
        long Calc_call = 0;

        bool isTinf;
        Model(int L, int bin, int B, int J, int Tsrt, int Tfin, bool isTinf);
        Model(vector<double> args);
        void ProbCalc(double beta);
        void Initialize(double beta);
        // void Initialzie(int idx);
        int SweepHelical(int i);
        int BoundaryHelical(int i);
        duo Measure();
        void Calculate(int _n = 0);
        void IterateUntilEquilibrium(int equil_time);
};

Model::Model(int L, int bin, int B, int J, int Tsrt, int Tfin, bool isTinf) :L(L), N(L*L), Bin(bin), B(B), J(J){
    this-> isTinf = isTinf;
    this-> sc = new short[N];

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);

    for(int i = 0; i < bin; i++){ 
        if(!Tsrt){
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin))*(i+1);
        } else if(bin == 1){
            this->TV[i] = Tsrt;
        } else {
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin-1))*(i);
        }
        this->BetaV[i] = 1/TV[i];
    }
}

Model::Model(vector<double> args) :L(args[0]), N(L*L), Bin(args[1]), B(args[2]), J(args[3]){
    this-> isTinf = args[6];
    this-> sc = new short[N];

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);

    double Tsrt = args[4];
    double Tfin = args[5];
    double bin  = args[1];

    for(int i = 0; i < bin; i++){ 
        if(!Tsrt){
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin))*(i+1);
        } else if(bin == 1){
            this->TV[i] = Tsrt;
        } else {
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin-1))*(i);
        }
        this->BetaV[i] = 1/TV[i];
    }
}


void Model::ProbCalc(double beta){
    for(int i = 2; i < 5; i += 2){
        this->prob[i] = exp(-2*beta*i);
    }
}

void Model::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;
        // cout << i << " " << sc[i] << endl;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    this->ProbCalc(beta);
    // for(int i = 0; i < N; i++){
    //     cout << i << " " <<sc[i] << endl;
    // }
}

// void Model::Initialize(int idx){
//     this-> Initialzie(this-> BetaV[idx]);
// }

int Model::SweepHelical(int i){
    int nn, sum = 0;
    int XNN = 1, YNN = L;

    if((nn = i - XNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + XNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];
    if((nn = i - YNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + YNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];

    return sum;
}

int Model::BoundaryHelical(int i){
    int nn, sum = 0;
    int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo Model::Measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = this->BoundaryHelical(i);
        res += J*sum*sc[i];
        sigma += sc[i];
    }
    
    HH = -res-B*sigma;

    return make_tuple(HH,sigma);
}

void Model::Calculate(int _n){
    int i, k, delta, n;
    n = !_n ? (this->N) : _n;
    double a;
    for(i = 0; i < n; i++){
        // Sweep Randomly
        k = (this->N)*dis(gen);
        // Sweep Sequential
        // k = n

        delta = (this->SweepHelical(k))*(this->sc[k]);
        
        a = dis(gen);
        this->Total_Step++;

        if(delta <= 0){
            this->Fliped_Step++;
            this->sc[k] = -(this->sc[k]);
        } else if(a < prob[delta]){
            this->Fliped_Step++;
            this->sc[k] = -(this->sc[k]);
        }
    }
}

void Model::IterateUntilEquilibrium(int equil_time){
    for(int j = 0; j < equil_time; j++)
        Calculate();
}

#endif // ____METROPOLIS____ 