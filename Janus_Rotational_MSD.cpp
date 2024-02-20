#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

const int Nstep = 5000001;
const int Ntau = 200000;

int main() {
    std::vector<float> tau(Ntau);
    std::vector<float> M(Ntau);
    float Rp[Nstep], Rq[Nstep], Rx[Nstep], Ry[Nstep];
    float cos_t[Nstep], sin_t[Nstep], theta[Nstep], cum[Nstep];
    float summ = 0.0f;
    float Pi = 3.1415927f;

    std::ifstream inputFile("input.dat");
    if (!inputFile) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    // Read data into arrays
    for (int i = 0; i < Nstep; ++i) {
        inputFile >> Rx[i] >> Ry[i] >> Rp[i] >> Rq[i];
    }
    inputFile.close();

    // Compute cos_theta and sin_theta
    for (int i = 0; i < Nstep; ++i) {
        float norm = std::sqrt(std::pow(Rp[i] - Rx[i], 2) + std::pow(Rq[i] - Ry[i], 2));
        cos_t[i] = (Rp[i] - Rx[i]) / norm;
        sin_t[i] = (Rq[i] - Ry[i]) / norm;
    }

    // Compute cumulative angle
    for (int i = 0; i < Nstep - 1; ++i) {
        float angle_diff = std::asin(0.5f * std::sqrt(std::pow(cos_t[i + 1] - cos_t[i], 2) +
                                                       std::pow(sin_t[i + 1] - sin_t[i], 2)));
        float sign_rot = ((cos_t[i] * sin_t[i + 1]) - (sin_t[i] * cos_t[i + 1])) > 0 ? 1.0f : -1.0f;
        theta[i] = 2 * angle_diff * sign_rot;
        summ += theta[i];
        cum[i] = summ;
    }

    std::ofstream outputFile("output.dat");
    if (!outputFile) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

    // Compute AMSD
    for (int j = 0; j < Ntau; ++j) {
        tau[j] = j;

        M[j] = 0.0f;
        for (int i = 0; i < Nstep - j; ++i) {
            float dR = std::pow(cum[i + j] - cum[i], 2);
            M[j] += dR;
        }
        M[j] /= Nstep - j;

        outputFile << tau[j] << " " << M[j] << std::endl;
    }

    outputFile.close();

    return 0;
}
