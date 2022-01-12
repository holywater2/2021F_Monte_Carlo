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
#include <sstream>
// #include <ifstream>


using namespace std;

#define L 5
#define N (L*L)
#define XNN 1
#define YNN L
#define B 0
#define J 1
#define bin 25

static int s[N];
static double prob[5]; // 1 1 exp 1 exp
static random_device rd;  // Will be used to obtain a seed for the random number engine
static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
static uniform_real_distribution<> dis(0.0, 1.0);
static int count;
static int count2;
vector<int> rand_int;
vector<double> rand_double;

void prob_calc(long double beta){
    int i;
    for(i = 2; i < 5; i += 2){
        prob[i] = expl(-2*beta*i);
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
    
    if((nn = i + XNN) == N) nn = 0;
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
    
    // HH = -res-B*sigma;
    HH = -res;

    return tuple<int,int>(HH,sigma);
}

int calculate(int n = N){
    int i, k, delta;
    double a;
    for(i = 0; i < n; i++){
        // k = N*dis(gen);
        k = rand_int[count2];

        delta = sweep_helical(k)*s[k];
        // cout << delta << " "  <<  k <<  " "<< sweep_helical(k) << " " << s[k] << "\n";
        
        // double a = dis(gen);
        a = rand_double[count2]; 
        count2++;

        if(delta <= 0){
            count++;
            s[k] = -s[k];
        } else if(a < prob[delta]){
            count++;
            s[k] = -s[k];
        }
    }
    // cout << k << " " << delta << " " << s[k] << "\n";
}



int main(){
    ifstream fin1;
    ifstream fin2;
    string line;
    string line2;
    string temp;
    stringstream ss;


    
    fin1.open("output.csv");
    fin2.open("output2.csv");

    fin1 >> line;
    fin2 >> line;

    while(!fin1.eof()){
        fin1 >> line;
        ss = stringstream(line);
        getline(ss, temp, ',');
        getline(ss, temp, ',');
        rand_int.push_back(stoi(temp));


        fin2 >> line;
        ss = stringstream(line);
        getline(ss, temp, ',');
        getline(ss, temp, ',');
        rand_double.push_back(stod(temp));
    }
    cout << rand_int.size() << '\n';
    cout << rand_double.size() << '\n';


    tuple<int,int> value;
    double m[bin];
    double c[bin];
    double res[bin][4];

    int sigma, HH;

    for(int i = 0; i < bin; i++){
        count = 0;
        count2 = 0;
        for(int j = 0; j < 4; j++) res[i][j] = 0;
        double T = (5/double(bin))*(i+1);
        long double beta = 1/T;
        // double beta = 1/(0.2*(i+1));

        initialize(beta);
        for(int j = 0; j < 2000; j++){
            calculate();
        }
        value = measure();
        HH = get<0>(value);
        sigma = get<1>(value);

        cout <<"idx: " << i << "||" << sigma << " " << HH << "\n";

        int epoch = 18000;
        for(int j = 0; j < epoch; j++){
            calculate();
            
            value = measure();
            HH = get<0>(value);
            sigma = get<1>(value);
           
            res[i][0] += abs(sigma)/double(epoch);
            res[i][1] += (sigma*sigma)/double(epoch);
            res[i][2] += HH/double(epoch);
            res[i][3] += (HH*HH)/double(epoch);
        }
        m[i] = res[i][0]/N;
        c[i] = (beta*beta)/N*(res[i][3]-(res[i][2]*res[i][2]));
        cout << m[i] << " " << c[i] << " " << count << " " << count2 << '\n';
        // cout << prob[2] << " " << prob[4] << "\n";
    }
    
    ofstream myfile;
    myfile.open("Metropolis_c_"+to_string(L)+".csv");
    myfile << "idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,HH,HH**2,exp(-4beta),exp(-8beta)\n";
    for(int i = 0; i <bin; i++){
        string temp = to_string(i) + "," + to_string(0.2*(i+1)) + "," + to_string(m[i]) + "," + to_string(c[i]) + ",";
        temp = temp + to_string(res[i][0]) + "," + to_string(res[i][1]) + "," + to_string(res[i][2]) + "," + to_string(res[i][3]) + "\n";
        myfile << temp;
    }
    myfile.close();
}