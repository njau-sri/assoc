#pragma once

#include <vector>

// ordinary least squares

class OLS
{
public:
    // X*b = y
    int fit(int p, const std::vector<double> &x, const std::vector<double> &y, double &dfe, double &sse);

    // X*b = y, subject to Z*b = 0
    int fit(int p, int q, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> zt, double &dfe, double &sse);

    const std::vector<double>& getb() const { return b_; }

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> b_;
};
