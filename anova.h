#pragma once

#include <memory>
#include <string>
#include <vector>

class ANOVA
{
public:

    struct Table
    {
        std::vector<std::string> names;
        std::vector<double> df;
        std::vector<double> ss;
        std::vector<double> ms;
        std::vector<double> f;
        std::vector<double> p;
        std::vector<double> error;
        std::vector<double> total;
    };

    struct Solution
    {
        std::vector<std::string> params;
        std::vector<double> coeffs;
    };

public:

    // regression effect, X
    void add_reg(const std::string &name, const std::vector<double> &x);

    // main effect, A
    void add_main(const std::string &name, const std::vector<std::string> &a);

    // crossed effect, A*B
    void add_crossed(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b);

    // nested effect, A(B)
    void add_nested(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b);

    Table solve1(const std::vector<double> &y) const;

    Table solve3(const std::vector<double> &y) const;

    Solution solution(const std::vector<double> &y) const;

private:

    struct Term
    {
        std::string name;
        std::vector<std::string> params;
        std::vector< std::vector<double> > data;
        std::vector< std::vector<double> > constr;
    };

private:

    std::vector< std::shared_ptr<Term> > tms_;
};
