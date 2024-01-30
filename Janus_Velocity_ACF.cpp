#include <iostream>
#include <fstream>
#include <vector>

const int Nstep = 2000000;  // Steps
const int Ntau = 1000000;   // Lagtime

int main() {
    std::vector<double> tau(Ntau), C(Ntau);
    double Ct0 = 0.0;
    double Vx[Nstep], Vy[Nstep], Vz[Nstep], t[Nstep];

    std::ifstream inputFile("Input.dat");
    std::ofstream outputFile("VACF.dat");

    for (int i = 0; i < Nstep; ++i) {
        inputFile >> t[i] >> Vx[i] >> Vy[i] >> Vz[i];
        Ct0 += (Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i]);
    }

    Ct0 /= Nstep;

    std::cout << "Ct0: " << Ct0 << std::endl;

    for (int k = 0; k < Ntau; ++k) {
        tau[k] = k;
        C[k] = 0.0;

        for (int i = 0; i < (Nstep - k + 1); ++i) {
            double vprod = (Vx[i] * Vx[i + k - 1]) + (Vy[i] * Vy[i + k - 1]) + (Vz[i] * Vz[i + k - 1]);
            C[k] += vprod;
        }

        C[k] /= (Nstep - k + 1);
        C[k] /= Ct0;

        outputFile << tau[k] << " " << C[k] << std::endl;
    }

    inputFile.close();
    outputFile.close();

    return 0;
}
