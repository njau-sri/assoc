#pragma once

#include <set>
#include <vector>
#include <numeric>
#include <iterator>
#include <algorithm>

using std::size_t;

void design1(const std::vector<int> &g, std::vector< std::vector<double> > &x);

void design2(const std::vector<int> &g, std::vector< std::vector<double> > &x);

void design3(const std::vector<int> &g, std::vector< std::vector<double> > &x);

double fpval(double x, double df1, double df2);

template<typename T1, typename T2>
size_t index(const std::vector<T1> &v, const T2 &a)
{
    return std::distance(v.begin(), std::find(v.begin(), v.end(), a));
}

template<typename T>
std::vector<size_t> order(const std::vector<T> &v)
{
    std::vector<size_t> z(v.size());
    std::iota(z.begin(), z.end(), size_t(0));
    std::sort(z.begin(), z.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });
    return z;
}

template<typename T>
std::vector<T> unique(std::vector<T> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

template<typename T>
std::vector<T> stable_unique(std::vector<T> v)
{
    std::set<T> seen;

    auto last = v.begin();
    for (auto itr = v.begin(); itr != v.end(); ++itr) {
        if (seen.insert(*itr).second) {
            if (last != itr)
                *last = *itr;
            ++last;
        }
    }

    v.erase(last, v.end());

    return v;
}

template<typename T1, typename T2>
std::vector<T1> subset(const std::vector<T1> &v, const std::vector<T2> &idx)
{
    std::vector<T1> z;
    z.reserve(idx.size());
    for (auto i : idx)
        z.push_back(v[i]);
    return z;
}

template<typename T>
std::vector<T> intersect(std::vector<T> a, std::vector<T> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::vector<T> c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));
    return c;
}

template<typename T>
std::vector<int> factor(const std::vector<T> &v)
{
    std::vector<int> gi;
    auto u = unique(v);
    gi.reserve(v.size());
    for (auto &e : v)
        gi.push_back(index(u, e));
    return gi;
}

template<typename T>
void factor(const std::vector<T> &v, std::vector<T> &gn, std::vector<int> &gi)
{
    auto u = unique(v);
    gi.clear();
    gi.reserve(v.size());
    for (auto &e : v)
        gi.push_back(index(u, e));
    gn.swap(u);
}

template<typename T>
double sumsqc(const std::vector<T> &v)
{
    auto sx = std::accumulate(v.begin(), v.end(), T(0));
    double mx = sx / v.size();
    double ssx = 0;
    for (auto x : v)
        ssx += (x - mx)*(x - mx);
    return ssx;
}
