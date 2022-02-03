#include "Metropolis.hpp"
#include "Writer.hpp"

#include <iostream>
#include <iomanip>

const int kL = 100; /*Parameter: lattice size*/
const int kN = kL*kL;
const int kBin = 40; /*Parametr: Change binning of temperature*/
const int kB = 0;
const int kJ = 1;

const int Tsrt = 0;
const int Tfin = 5;

double isTinf = false;

#ifdef _WIN32
static string kFilename = ".\\Result\\Metropolis_c_Error_"+to_string(kL)+"_int"+to_string(kBin);
#endif
#ifdef linux
static string kFilename = "./Result/Metropolis_c_Error_"+to_string(kL)+"_int"+to_string(kBin);
#endif


vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << "Metropolis Algorithm\n";
    cout << "Radnomness test(seed): " << seed << endl;
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat << "\n";
    cout << "-------------------------------------------------------------------------------------------" << "\n";
    cout << "index---||-sig--HH-||magnetization-specific heat--Fliped--Total Step-----------------------" << "\n";
    cout << "-------------------------------------------------------------------------------------------" << "\n";
    cout << fixed << setprecision(6);

    __start__ = clock();
}

void Farewell(){
    __finish__ = clock();
    cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
}

int main(){
    Greetings();

    Model model = Model(args);
    Writer modelW = Writer(kFilename);
    Writer modelW2 = Writer(kFilename + "_auto");
    modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2,m_error\n");

    int sigma, HH, equil_time, epoch = 20000;
    double epoch_i = 1/double(epoch);

    for(int i = 0; i < kBin; i++){
        model.Initialize(model.BetaV[i]);
        model.res = vector<double>(5,0);

        // if(model.TV[i]<2.4 || model.TV[i]>2.0) equil_time =10000;

        // model.IterateUntilEquilibrium(2000);

        duo value = model.Measure();
        HH = get<0>(value);
        sigma =  get<1>(value);

        cout <<"idx: " << i << "\t|| " << sigma << "\t" << HH << setw(3) << "|| ";
        modelW2.WriteLine(to_string(model.TV[i]));

        epoch = 20000;
        string temp2 = "";
        for(int j = 0; j < epoch; j++){
            model.Calculate();

            value = model.Measure();
            HH = get<0>(value);
            sigma =  get<1>(value);
            
            model.res[0] += abs(sigma)*epoch_i;
            model.res[1] += (sigma*epoch_i*sigma);
            model.res[2] += (sigma*epoch_i*sigma)*(sigma*sigma);
            model.res[3] += HH*epoch_i;
            model.res[4] += HH*epoch_i*HH;
            
            temp2 += ","+ to_string(sigma);
        }
        modelW2.WriteLine(temp2 + "\n");
        delete temp2;

        model.MV[i] = model.res[0]/kN;
        model.CV[i] = (model.BetaV[i]*model.BetaV[i])/kN*(model.res[4]-model.res[3]*model.res[3]);        

        cout << model.MV[i] << "\t" << model.CV[i] << "\t" << model.Fliped_Step << setw(10) << "\t" << model.Total_Step << '\n';

        string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
        temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";
        modelW.WriteLine(temp);
    }
    modelW.CloseNewFile();
    modelW2.CloseNewFile();
    Farewell();
}