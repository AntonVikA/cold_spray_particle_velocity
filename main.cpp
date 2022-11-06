//
// Created by anton on 10/15/22.
//

#include <iostream>
#include <vector>
#include <functional>
#include "RKM.h"
#include <cmath>
#include <fstream>

const double rho_particle = 5000.0, d_particle = 5.0 * std::pow(10.0, -6.0), rho_air = 1.2, a = 330.0, mu = 18.27, C = 0.1;
double PI = 4.0 * std::atan(1.0);

double ExectSol(double vp) {
    double v = 600.0;
    return 1 / C * (v / (v - vp) + log(v - vp) - (1 + log(v)));
}

double V(double z) {
    return 600.0; //10.0 + 600.0 * std::tanh(0.1 * z);
}

double f1(double t, const std::vector<double>& y) {
    return y[1];
}

double f2(double t, const std::vector<double>& y) {
    return C * std::pow(std::abs(V(y[0]) - y[1]), 2.0);
}

int main() {

    const std::vector<std::function<double(double, const std::vector<double>&)>> functions = {f1, f2};
    std::vector<double> initial_condition{0., 0.};

    double t_begin = 0.0;
    double t_end = 5.00;
    double epsilon_computed;

    RKM rkm;
    std::ofstream file;
    file.open("speed.txt");

    rkm.init(functions);

    size_t iter = 0;
    for (auto ht : std::vector<double>{
                                       0.0001, 0.0002,
                                       0.0005,
                                       0.001, 0.002,
                                       0.005,
                                       0.01, 0.02}) {
        rkm.solve(t_begin, t_end, initial_condition, ht, epsilon_computed);
        std::cout << ht << " " << epsilon_computed << std::endl;
        file << std::log(ht) << " " << std::log(std::abs(epsilon_computed)) << std::endl;

    }

    file.close();


    return 0;
}
