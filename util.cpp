#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions/fisher_f.hpp>
#include "util.h"

// 0/-1/1 coding, full rank, drop last level
void design1(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    int n = g.size();
    int m = g.empty() ? 0 : *std::max_element(g.begin(), g.end());

    x.assign(m, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        int k = g[i];
        if (k != m)
            x[k][i] = 1.0;
        else
            for (int j = 0; j < m; ++j)
                x[j][i] = -1.0;
    }
}

// 0/1 coding, full rank, drop last level
void design2(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    int n = g.size();
    int m = g.empty() ? 0 : *std::max_element(g.begin(), g.end());

    x.assign(m, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        if (g[i] != m)
            x[g[i]][i] = 1.0;
}

// 0/1 coding, overdetermined
void design3(const std::vector<int> &g, std::vector< std::vector<double> > &x)
{
    int n = g.size();
    int m = g.empty() ? 0 : *std::max_element(g.begin(), g.end());

    x.assign(m + 1, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        x[g[i]][i] = 1.0;
}

double fpval(double x, double df1, double df2)
{
    boost::math::fisher_f f(df1, df2);
    return boost::math::cdf(boost::math::complement(f, x));
}
