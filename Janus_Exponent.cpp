#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

const int Nstep = 300000;
const int m = 1;

int main() {
    std::vector<double> Rx(Nstep), t(Nstep), dy(Nstep), dx(Nstep), A(Nstep);

    std::ifstream inputFile("file.dat");
    std::ofstream outputFile("output.dat");

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }

    for (int i = 0; i < Nstep; ++i) {
        inputFile >> t[i] >> Rx[i];
    }

    for (int i = 0; i < Nstep - 1; ++i) {
        dy[i] = Rx[i + m] - Rx[i];
        dx[i] = t[i + m] - t[i];
        A[i] = dy[i] / dx[i];
        outputFile << std::exp(t[i]) << " " << A[i] << std::endl;
    }

    inputFile.close();
    outputFile.close();

    return 0;
}
