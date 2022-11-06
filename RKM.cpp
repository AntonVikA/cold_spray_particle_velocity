//
// Created by anton on 10/15/22.
//

#include "RKM.h"

double ExectSol(double vp);
double V(double z);

std::vector<double> operator+ (const std::vector<double>& v_left, const std::vector<double>& v_right) {
    std::vector<double> res(v_left.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = v_left[i] + v_right[i];
    }

    return res;
}

std::vector<double> operator- (const std::vector<double>& v_left, const std::vector<double>& v_right) {
    std::vector<double> res(v_left.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = v_left[i] + v_right[i];
    }

    return res;
}

std::vector<double> operator* (double h, const std::vector<double>& v_right) {
    std::vector<double> res(v_right.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = h * v_right[i];
    }

    return res;
}

std::vector<double> operator/ (double h, const std::vector<double>& v_right) {
    std::vector<double> res(v_right.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = v_right[i] / h;
    }

    return res;
}

std::vector<double> operator* (const std::vector<double>& v_right, double h) {
    std::vector<double> res(v_right.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = h * v_right[i];
    }

    return res;
}

std::vector<double> operator/ (const std::vector<double>& v_right, double h) {
    std::vector<double> res(v_right.size(), 0);

    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = v_right[i] / h;
    }

    return res;
}

std::vector<double> RKM::operator() (double t, const std::vector<double>& y) {
    std::vector<double> res(number_of_equations, 0);
    for (size_t i = 0; i < number_of_equations; ++i) {
        res[i] = functions_[i](t, y);
    }

    return res;
}

RKM::RKM(const std::vector<std::function<function_type>>& functions) : functions_{functions}, number_of_equations{functions_.size()} {

}

void RKM::init(const std::vector<std::function<function_type>>& functions) {
    functions_ = functions;
    functions_[0](0.0, {0.0});
    number_of_equations = functions_.size();
    Y = std::vector<double>(number_of_equations, 0);
    dY = std::vector<double>(number_of_equations, 0);
}

double calc_epsilon(const std::vector<double>& epsilon_vec, const std::vector<double>& Y) {
    const double sum_of_epsilon_power2 = std::inner_product(epsilon_vec.begin(), epsilon_vec.end(), epsilon_vec.begin(), double(0));
    const double sum_of_y_power2 = std::inner_product(Y.begin(), Y.end(), Y.begin(), double(0));
    return sum_of_y_power2 < 1 ? 0.000001 : sqrt(sum_of_epsilon_power2 / sum_of_y_power2);
}

void RKM::solve(double t_begin, double t_end, const std::vector<double>& init_conditions, double h, double& epsilon_computed) {
    assert(init_conditions.size() == functions_.size());

    t_ = t_begin;

    Y = init_conditions;

    //file.open("speed.txt");
    size_t iter = 0;
    size_t n = t_end / h;
    std::cout << " n = " << n << std::endl;
    for (size_t i = 0; i < n; ++i) {
        /*std::cout << t_ << " " << Y[0] << std::endl;
        if (iter % 100 == 0) {
            file << std::fixed << t_;
            for (auto el: Y)
                file << " " << el;
            file << " " << V(Y[0]) << " "<< ExectSol(Y[1]) <<std::endl;
        }*/
        K1 = h * (*this)(t_, Y);
        K2 = h * (*this)(t_ + h / 2.0, Y + K1 / 2.0);
        K3 = h * (*this)(t_ + h / 2.0, Y + K2 / 2.0);
        K4 = h * (*this)(t_ + h, Y + K3);
        //K5 = h / 3.0 * (*this)(t_ + h, Y + 3.0 / 2.0 * K1 - 9.0 / 2.0 * K3 + 6.0 * K4);

        //epsilon_vec = 1.0 / 5.0 * (K1 + 4.0 * K4 - 9.0 / 2.0 * K3 - K5 / 2.0);
        //epsilon_t = calc_epsilon(epsilon_vec, Y);
        //h = h * std::pow(std::abs(epsilon_t / epsilon) + 0.001, -0.2);
        t_ += h;
        dY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
        Y = Y + dY;
        //file << Y[1] << " " << ExectSol(Y[1]) << " " << Y[0] << std::endl;
    }
    //std::cout << "solution: " << ExectSol(t_) << " " << Y[0] << " " << ExectSol(t_)  - Y[0] << std::endl;
    epsilon_computed = std::abs(ExectSol(Y[1]) - Y[0]);

    //file.close();
}