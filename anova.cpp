#include <limits>
#include <functional>
#include "anova.h"
#include "util.h"
#include "ols.h"

namespace
{
    static const double kNaN = std::numeric_limits<double>::quiet_NaN();

} // namespace

void ANOVA::add_reg(const std::string &name, const std::vector<double> &x)
{
    auto tm = std::make_shared<Term>();
    tm->name = name;
    tm->params.push_back(name);
    tm->data.push_back(x);
    tms_.push_back(tm);
}

void ANOVA::add_main(const std::string &name, const std::vector<std::string> &a)
{
    std::vector<int> gi;
    std::vector<std::string> gn;
    factor(a, gn, gi);

    std::vector< std::vector<double> > x;
    design3(gi, x);

    auto tm = std::make_shared<Term>();
    tm->name = name;
    for (auto &e : gn)
        tm->params.push_back(name + "\t" + e);
    tm->constr.emplace_back(gn.size(), 1.0);
    tm->data.swap(x);
    tms_.push_back(tm);
}

void ANOVA::add_crossed(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    std::vector<int> gia;
    std::vector<std::string> gna;
    factor(a, gna, gia);

    std::vector<int> gib;
    std::vector<std::string> gnb;
    factor(b, gnb, gib);

    std::vector< std::vector<double> > xa;
    design3(gia, xa);

    std::vector< std::vector<double> > xb;
    design3(gib, xb);

    auto tm = std::make_shared<Term>();
    tm->name = name;

    int na = gna.size();
    int nb = gnb.size();

    // A1*B1 A1*B2 A1*B3, A2*B1 A2*B2 A2*B3
    for (int i = 0; i < na; ++i) {
        for (int j = 0; j < nb; ++j) {
            tm->params.push_back(name + "\t" + gna[i] + "*" + gnb[j]);
            auto v = xb[j];
            std::transform(v.begin(), v.end(), xa[i].begin(), v.begin(), std::multiplies<double>());
            tm->data.push_back(v);
        }
    }

    // zero-sum constraints for A
    for (int j = 0; j < nb; ++j) {
        std::vector<double> v(na*nb, 0);
        for (int i = 0; i < na; ++i)
            v[i*nb + j] = 1.0;
        tm->constr.push_back(v);
    }

    // zero-sum constraints for B
    for (int i = 0; i < na; ++i) {
        std::vector<double> v(na*nb, 0);
        std::fill_n(v.begin() + i*nb, nb, 1.0);
        tm->constr.push_back(v);
    }

    tms_.push_back(tm);
}

void ANOVA::add_nested(const std::string &name, const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    std::vector<int> gia;
    std::vector<std::string> gna;
    factor(a, gna, gia);

    std::vector<int> gib;
    std::vector<std::string> gnb;
    factor(b, gnb, gib);

    std::vector< std::vector<double> > xa;
    design3(gia, xa);

    std::vector< std::vector<double> > xb;
    design3(gib, xb);

    auto tm = std::make_shared<Term>();
    tm->name = name;

    int na = gna.size();
    int nb = gnb.size();

    // A1(B1) A2(B1) A3(B1), A1(B2) A2(B2) A3(B2)
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < na; ++j) {
            tm->params.push_back(name + "\t" + gna[j] + "(" + gnb[i] + ")");
            auto v = xa[j];
            std::transform(v.begin(), v.end(), xb[i].begin(), v.begin(), std::multiplies<double>());
            tm->data.push_back(v);
        }
        // zero-sum constraints for A
        std::vector<double> v(na*nb, 0);
        std::fill_n(v.begin() + i*na, na, 1.0);
        tm->constr.push_back(v);
    }

    tms_.push_back(tm);
}

ANOVA::Table ANOVA::solve1(const std::vector<double> &y) const
{
    int n = y.size();

    OLS ols;
    Table tab;
    std::vector<double> x;

    int q0 = 1;
    x.assign(n, 1);

    double dfe0 = n - 1;
    double sse0 = sumsqc(y);

    tab.total.push_back(dfe0);
    tab.total.push_back(sse0);

    for (auto &tm : tms_) {
        int q1 = q0 + tm->data.size();
        for (auto &e : tm->data)
            x.insert(x.end(), e.begin(), e.end());

        double dfe1, sse1;
        ols.fit(q1, x, y, dfe1, sse1);

        tab.names.push_back(tm->name);
        tab.df.push_back(dfe0 - dfe1);
        tab.ss.push_back(sse0 - sse1);

        q0 = q1;
        dfe0 = dfe1;
        sse0 = sse1;
    }

    tab.error.push_back(dfe0);
    tab.error.push_back(sse0);
    tab.error.push_back(sse0 / dfe0);

    int m = tms_.size();
    tab.ms.resize(m);
    tab.f.resize(m);
    tab.p.assign(m, kNaN);
    for (int j = 0; j < m; ++j) {
        tab.ms[j] = tab.ss[j] / tab.df[j];
        tab.f[j] = tab.ms[j] / tab.error[2];
        if (tab.df[j] > 0.0 && tab.ss[j] > 0.0)
            tab.p[j] = fpval(tab.f[j], tab.df[j], tab.error[0]);
    }

    return tab;
}

ANOVA::Table ANOVA::solve3(const std::vector<double> &y) const
{
    auto n = y.size();

    OLS ols;
    Table tab;
    std::vector<double> x;

    int q1 = 1;
    x.assign(n, 1);

    double dft = n - 1;
    double sst = sumsqc(y);

    tab.total.push_back(dft);
    tab.total.push_back(sst);

    for (auto &tm : tms_) {
        q1 += tm->data.size();
        for (auto &e : tm->data)
            x.insert(x.end(), e.begin(), e.end());
    }

    double dfe, sse;
    ols.fit(q1, x, y, dfe, sse);

    tab.error.push_back(dfe);
    tab.error.push_back(sse);
    tab.error.push_back(sse / dfe);

    for (auto &curr : tms_) {
        int q0 = 1;
        x.assign(n, 1);

        for (auto &tm : tms_) {
            if (&tm == &curr)
                continue;
            q0 += tm->data.size();
            for (auto &e : tm->data)
                x.insert(x.end(), e.begin(), e.end());
        }

        double dfe0, sse0;
        ols.fit(q0, x, y, dfe0, sse0);

        tab.names.push_back(curr->name);
        tab.df.push_back(dfe0 - dfe);
        tab.ss.push_back(sse0 - sse);
    }

    int m = tms_.size();
    tab.ms.resize(m);
    tab.f.resize(m);
    tab.p.assign(m, kNaN);
    for (int j = 0; j < m; ++j) {
        tab.ms[j] = tab.ss[j] / tab.df[j];
        tab.f[j] = tab.ms[j] / tab.error[2];
        if (tab.df[j] > 0.0 && tab.ss[j] > 0.0)
            tab.p[j] = fpval(tab.f[j], tab.df[j], tab.error[0]);
    }

    return tab;
}

ANOVA::Solution ANOVA::solution(const std::vector<double> &y) const
{
    int n = y.size();
    std::vector<double> x;
    std::vector<std::string> params;

    int p = 1;
    x.assign(n, 1);
    params.push_back("Constant");

    for (auto &tm : tms_) {
        params.insert(params.end(), tm->params.begin(), tm->params.end());
        p += tm->data.size();
        for (auto &e : tm->data)
            x.insert(x.end(), e.begin(), e.end());
    }

    std::vector<double> zt;
    int pos = 1, q = 0;
    for (auto &tm : tms_) {
        for (auto &e : tm->constr) {
            std::vector<double> v(p, 0);
            std::copy(e.begin(), e.end(), v.begin() + pos);
            zt.insert(zt.end(), v.begin(), v.end());
            q += 1;
        }
        pos += tm->data.size();
    }

    OLS ols;
    double dfe, sse;
    ols.fit(p, q, x, y, zt, dfe, sse);

    Solution sol;
    sol.params.swap(params);
    sol.coeffs = ols.getb();

    return sol;
}
