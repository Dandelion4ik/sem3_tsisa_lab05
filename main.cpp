#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ostream>
#include <iomanip>

double value(double x, double c, double d) {
    return c * x + d;
}

void print_noise_value(const std::vector<std::pair<double, double>> &a) {
    for (const auto &it : a) {
        std::cout << "x= " << std::setw(5) << std::left << it.first << std::right << std::setw(6) << " y= " << it.second
                  << std::endl;
    }
    std::cout << std::endl;
}

void print(const std::pair<double, double> &point) {
    std::cout << "c= " << point.first << " d= " << point.second << std::endl;
}

void print_value_gold_passive(const std::pair<double, double> &m, const int& N, const double& a) {
    double x = a;
    for (int i = 0; i < N; ++i) {
        std::cout << "x="<<std::setw(6)<<std::left << x <<std::right<< std::setw(10) << " y=" << value(x, m.first, m.second) << std::endl;
        x += sqrt(a*a) / N;
    }
}

void printValueFunc(const int &N, const double &a, const double &b, const double &c, const double &d) {
    double time_x = a;
    for (int i = 0; i < N; ++i) {
        std::cout << "x" << std::setw(2) << std::left << i + 1 << "= " << std::setprecision(3) << std::setw(6)
                  << std::left << time_x << std::setw(6) << std::right << " y" << i + 1 << "= " << value(time_x, c, d)
                  << std::endl;
        time_x += (b - a) / N;
    }
    std::cout << std::endl;
}

std::pair<double, double> E(const double &c, const double &d, const int &N,
                            const std::vector<std::pair<double, double>> &points) {
    double sum = 0;
    double sum_xi_t = 0;
    double sum_of_x = 0;
    double sum_of_t = 0;
    double sum_of_x_v_stepeni_2 = 0;
    for (const auto &i : points) {
        sum += pow(value(i.first, c, d) - i.second, 2);
        sum_xi_t += i.first * i.second;
        sum_of_x += i.first;
        sum_of_t += i.second;
        sum_of_x_v_stepeni_2 += pow(i.first, 2);
    }
    double c_ = (static_cast<double>(N) * sum_xi_t - sum_of_x * sum_of_t)
                / (static_cast<double>(N) * sum_of_x_v_stepeni_2 - pow(sum_of_x, 2));
    double d_ = (sum_of_t - c_ * sum_of_x) / static_cast<double>(N);
    std::pair<double, double> coefficients(c_, d_);
    return coefficients;
}

std::vector<std::pair<double, double>> noise_fun(const double &a, const double &b, const double &c,
                                                 const double &d, const double &A, const int &N) {
    std::vector<std::pair<double, double>> points;
    std::vector<double> v_ei;
    std::vector<double> v_xi;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    double xi = a;
    for (int i = 0; i < N; ++i) {
        double arg = dist(gen) * 0.5;
        double ei = arg * A;
        std::pair<double, double> point(xi, value(xi, c, d) + ei);
        points.emplace_back(point);
        xi += (b - a) / static_cast<double>(N);
    }
    return points;
}

std::pair<double, double> gold_passive(const int &N, const double &c,const double &d,
                              const std::vector<std::pair<double, double>> &points) {
    double d_min = -2;
    double d_max = 2.5;
    double epsi=0.2;
    std::vector<double> dif_d = {};
    for (auto i = 0; i < N; ++i) {
        double d = d_min + (d_max + d_min) * i / (N - 1);
        dif_d.push_back(d);
    }
    double min = 1000;
    double search_d = 0;
    for (auto &it: dif_d) {
        std::pair<double, double> mist_d = E(c, it, N, points);
        double sqrt = mist_d.first + mist_d.second;
        sqrt = sqrt * sqrt;
        if (min > sqrt) {
            search_d = mist_d.second+epsi;
        }
    }
    double c_min=7;
    double c_max=9.5;
    double  eps=0.15;
    std::vector<double> dif_c = {};
    for (auto i = 0; i < N; ++i) {
        double c_ar= c_min + (c_max + c_min) * i / (N - 1);
        dif_c.push_back(c_ar);
    }
    min = 1000;
    double search_c = 0;
    for (auto &it: dif_c) {
        std::pair<double, double> mist_c = E(it, d, N, points);
        double sqrt = mist_c.first + mist_c.second;
        sqrt = sqrt * sqrt;
        if (min > sqrt) {
            search_c = mist_c.first-eps;
        }
    }
    std::pair<double, double> ret{search_c, search_d};
    return ret;

}

int main() {
    const double a = -4;
    const double b = 2;
    const double c = 8;
    const double d = 0;
    const int N = 24;
    const double A = 10;
    std::cout << "The value of the function" << std::endl << std::endl;
    printValueFunc(N, a, b, c, d);
    std::cout << "Value of the function with interference" << std::endl << std::endl;
    print_noise_value(noise_fun(a, b, c, d, A, N));
    std::cout << "coefficients c and d" << std::endl;
    print(gold_passive(N, c,d, noise_fun(a, b, c, d, A, N)));
    std::cout << "Value of the function with new coefficients" << std::endl;
    print_value_gold_passive(gold_passive(N, c,d, noise_fun(a, b, c, d, A, N)), N, a);
    std::cout << "coefficients c and d with synaptic weights" << std::endl;
    print(E(c, d, N, noise_fun(a, b, c, d, A, N)));
    print_value_gold_passive(E(c, d, N, noise_fun(a, b, c, d, A, N)), N, a);
    std::cout << std::endl;
    return 0;
}