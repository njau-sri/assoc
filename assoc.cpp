#include <cmath>
#include <cctype>
#include <limits>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>
#include "vcfio.h"
#include "cmdline.h"
#include "split.h"
#include "util.h"
#include "ols.h"
#include "anova.h"

using std::size_t;

namespace
{
    static const double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Par
    {
        std::string vcf;
        std::string pheno;
        std::string covar;
        std::string out;
        double rsq = 0.99;
        double alpha = 0.05;
        double preselect = 0.05;
        int mtc = 0;
        int sstype = 1;
        bool nogxe = false;
    };

    struct Phenotype
    {
        std::vector<std::string> ind;
        std::vector<std::string> phe;
        std::vector<std::string> env;
        std::vector<std::string> blk;
        std::vector< std::vector<double> > dat;
    };

    struct Covariate
    {
        std::vector<std::string> ind;
        std::vector<std::string> phe;
        std::vector< std::vector<double> > dat;
    };

    Par par;

    int read_pheno(const std::string &filename, Phenotype &pt)
    {
        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
            return 1;
        }

        size_t ln = 0;
        std::vector<std::string> colnames;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);
            split(line, " \t", colnames);
            if (!colnames.empty())
                break;
        }

        std::vector<size_t> jphe;
        size_t jenv = 0, jblk = 0;

        for (size_t j = 1; j < colnames.size(); ++j) {
            if (colnames[j] == "_ENV_")
                jenv = j;
            else if (colnames[j] == "_BLK_")
                jblk = j;
            else
                jphe.push_back(j);
        }

        for (auto j : jphe)
            pt.phe.push_back(colnames[j]);

        std::vector< std::vector<double> > dat;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);

            std::vector<std::string> vs;
            split(line, " \t", vs);
            if (vs.empty())
                continue;

            if (vs.size() != colnames.size()) {
                std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size() << "!=" << colnames.size() << "\n";
                return 1;
            }

            pt.ind.push_back(vs[0]);

            if (jenv > 0)
                pt.env.push_back(vs[jenv]);

            if (jblk > 0)
                pt.blk.push_back(vs[jblk]);

            std::vector<double> v;

            for (auto j : jphe) {
                if (vs[j] == "?" || vs[j] == "NA" || vs[j] == ".")
                    v.push_back(kNaN);
                else
                    v.push_back(std::stod(vs[j]));
            }

            dat.push_back(v);
        }

        int m = pt.phe.size();
        int n = pt.ind.size();

        for (int j = 0; j < m; ++j) {
            std::vector<double> v(n);
            for (int i = 0; i < n; ++i)
                v[i] = dat[i][j];
            pt.dat.push_back(v);
        }

        return 0;
    }

    int read_covar(const std::string &filename, Covariate &ct)
    {
        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
            return 1;
        }

        size_t ln = 0;
        std::vector<std::string> colnames;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);
            split(line, " \t", colnames);
            if (!colnames.empty())
                break;
        }

        if (colnames.size() < 2)
            return 0;

        ct.phe.assign(colnames.begin() + 1, colnames.end());

        std::vector< std::vector<double> > dat;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);

            std::vector<std::string> vs;
            split(line, " \t", vs);
            if (vs.empty())
                continue;

            if (vs.size() != colnames.size()) {
                std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size() << "!=" << colnames.size() << "\n";
                return 1;
            }

            ct.ind.push_back(vs[0]);
            vs.erase(vs.begin());

            int m = ct.phe.size();
            std::vector<double> v;

            for (auto &e : vs) {
                auto a = (e == "?" || e == "NA" || e == ".") ? kNaN : std::stod(e);
                if (!std::isfinite(a)) {
                    std::cerr << "ERROR: missing value is not allowed for covariates: " << e << "\n";
                    return 1;
                }
                v.push_back(a);
            }

            dat.push_back(v);
        }

        int m = ct.phe.size();
        int n = ct.ind.size();

        for (int j = 0; j < m; ++j) {
            std::vector<double> v(n);
            for (int i = 0; i < n; ++i)
                v[i] = dat[i][j];
            ct.dat.push_back(v);
        }

        return 0;
    }

    void filter_ind(const std::vector<int> &idx, Genotype &gt)
    {
        if (gt.ploidy == 1) {
            for (auto &v : gt.dat)
                subset(v, idx).swap(v);
        }
        else {
            std::vector<int> idx2;
            for (auto i : idx) {
                idx2.push_back(i * 2);
                idx2.push_back(i * 2 + 1);
            }
            for (auto &v : gt.dat)
                subset(v, idx2).swap(v);
        }

        subset(gt.ind, idx).swap(gt.ind);
    }

    void filter_ind(const std::vector<int> &idx, Phenotype &pt)
    {
        subset(pt.ind, idx).swap(pt.ind);

        if (!pt.env.empty())
            subset(pt.env, idx).swap(pt.env);

        if (!pt.blk.empty())
            subset(pt.blk, idx).swap(pt.blk);

        for (auto &v : pt.dat)
            subset(v, idx).swap(v);
    }

    void merge_data(Genotype &gt, Phenotype &pt, Covariate &ct, std::vector<int> &gi)
    {
        bool has_covar = !ct.phe.empty() && !ct.ind.empty();

        bool require = gt.ind != pt.ind;

        if (!require && has_covar)
            require = gt.ind != ct.ind;

        if (!require) {
            gi.resize(gt.ind.size());
            std::iota(gi.begin(), gi.end(), 0);
            return;
        }

        std::cerr << "INFO: merging data by individual...\n";

        auto ind = intersect(gt.ind, pt.ind);
        if (has_covar)
            ind = intersect(ind, ct.ind);

        std::vector<int> pi, ci;

        for (auto itr = pt.ind.begin(); itr != pt.ind.end(); ++itr) {
            if (std::binary_search(ind.begin(), ind.end(), *itr)) {
                pi.push_back(itr - pt.ind.begin());
                gi.push_back(index(gt.ind, *itr));
                if (has_covar)
                    ci.push_back(index(ct.ind, *itr));
            }
        }

        filter_ind(pi, pt);

        if (has_covar) {
            subset(ct.ind, ci).swap(ct.ind);
            for (auto &v : ct.dat)
                subset(v, ci).swap(v);
        }

        std::cerr << "INFO: after merging data, there are " << ind.size() << " individuals\n";
    }

    void parse_factor(Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic)
    {
        std::vector< std::vector<double> > xenv, xblk;

        if (!pt.env.empty()) {
            design1(factor(pt.env), xenv);
            ac.insert(ac.end(), xenv.begin(), xenv.end());
            ic.insert(ic.end(), xenv.begin(), xenv.end());
        }

        if (!pt.blk.empty()) {
            design1(factor(pt.blk), xblk);
            if (xenv.empty()) {
                ac.insert(ac.end(), xblk.begin(), xblk.end());
            }
            else {
                design3(factor(pt.env), xenv);
                for (auto &e : xenv) {
                    for (auto v : xblk) {
                        std::transform(e.begin(), e.end(), v.begin(), v.begin(), std::multiplies<double>());
                        ac.push_back(v);
                    }
                }
            }
        }
    }

    int glm(bool haploid, const std::vector< std::vector<allele_t> > &gt, const std::vector<int> &gi, const std::vector<double> &y,
        const std::vector< std::vector<double> > &ac, const std::vector< std::vector<double> > &ic, std::vector<double> &p1, std::vector<double> &p2)
    {
        int n = y.size();
        int m = gt.size();

        p1.assign(m, kNaN);
        if (!ic.empty())
            p2.assign(m, kNaN);

        OLS ols;

        std::vector<int> idx;
        std::vector<allele_t> g1;
        std::vector< std::pair<allele_t, allele_t> > g2;
        std::vector< std::vector<double> > xg;
        std::vector<double> xv, yv;

        int q0 = 1 + ac.size();
        std::vector<double> x0(n, 1.0);
        for (auto &e : ac)
            x0.insert(x0.end(), e.begin(), e.end());

        double dfe0, sse0;
        ols.fit(q0, x0, y, dfe0, sse0);

        for (int j = 0; j < m; ++j) {
            idx.clear();

            if (haploid) {
                g1.clear();
                for (int i = 0; i < n; ++i) {
                    auto a = gt[j][gi[i]];
                    if (a) {
                        idx.push_back(i);
                        g1.push_back(a);
                    }
                }
                if (idx.empty())
                    continue;
                design2(factor(g1), xg);
            }
            else {
                g2.clear();
                for (int i = 0; i < n; ++i) {
                    auto a = gt[j][gi[i] * 2];
                    auto b = gt[j][gi[i] * 2 + 1];
                    if (a && b) {
                        idx.push_back(i);
                        if (a > b)
                            std::swap(a, b);
                        g2.emplace_back(a, b);
                    }
                }
                if (idx.empty())
                    continue;
                design2(factor(g2), xg);
            }

            if (xg.empty())
                continue;

            int nv = idx.size();

            auto jdfe0 = dfe0;
            auto jsse0 = sse0;

            if (nv != n) {
                yv = subset(y, idx);

                xv.assign(nv, 1.0);
                for (auto &e : ac) {
                    auto v = subset(e, idx);
                    xv.insert(xv.end(), v.begin(), v.end());
                }
                ols.fit(q0, yv, xv, jdfe0, jsse0);
            }
            else {
                yv = y;
                xv = x0;
            }

            int q1 = q0 + xg.size();
            for (auto &e : xg)
                xv.insert(xv.end(), e.begin(), e.end());

            double jdfe1, jsse1;
            ols.fit(q1, xv, yv, jdfe1, jsse1);

            auto dfx = jdfe0 - jdfe1;
            auto ssx = jsse0 - jsse1;
            auto dfe = jdfe1;
            auto sse = jsse1;

            if (!ic.empty()) {
                auto q2 = q1 + xg.size()*ic.size();
                for (auto &e1 : xg) {
                    for (auto e2 : ic) {
                        if (nv != n) subset(e2, idx).swap(e2);
                        std::transform(e1.begin(), e1.end(), e2.begin(), e2.begin(), std::multiplies<double>());
                        xv.insert(xv.end(), e2.begin(), e2.end());
                    }
                }

                ols.fit(q2, xv, yv, dfe, sse);

                auto dfxi = jdfe1 - dfe;
                auto ssxi = jsse1 - sse;

                if (dfxi > 0.0 && ssxi > 0.0) {
                    auto f = (ssxi / dfxi) / (sse / dfe);
                    p2[j] = fpval(f, dfxi, dfe);
                }
            }

            if (dfx > 0.0 && ssx > 0.0) {
                auto f = (ssx / dfx) / (sse / dfe);
                p1[j] = fpval(f, dfx, dfe);
            }
        }

        return 0;
    }

    void forward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<int> &c1, const std::vector<double> &y,
        std::vector<double> &ps, std::vector<int> &in)
    {
        int n = y.size();
        int m = c1.size();

        ps.assign(m, kNaN);

        std::vector<const double*> xs;
        xs.push_back(x1.data());
        for (int i = 1; i < m; ++i)
            xs.push_back(xs[i - 1] + n*c1[i - 1]);

        double sst = sumsqc(y);

        OLS ols;
        std::vector<double> x;
        std::vector<bool> ignore(m, false);

        auto q0 = 1 + x0.size() / n;
        x.assign(n, 1.0);
        x.insert(x.end(), x0.begin(), x0.end());

        double dfe0, sse0;
        ols.fit(q0, x, y, dfe0, sse0);

        auto rsq = 1 - sse0 / sst;
        std::cerr << "INFO: forward selection 0: " << rsq << "\n";

        for (int step = 0; step < m; ++step) {
            int idx = 0;
            double dfe = 0.0, sse = 0.0, pval = 1.0;

            for (int j = 0; j < m; ++j) {
                if (ignore[j])
                    continue;

                int q1 = q0 + c1[j];
                x.resize(n*q0);
                x.insert(x.end(), xs[j], xs[j] + n*c1[j]);

                double dfe1, sse1;
                ols.fit(q1, x, y, dfe1, sse1);

                if (dfe1 > 0.0 && sse1 > 0.0 && dfe0 > dfe1 && sse0 > sse1) {
                    auto f = ((sse0 - sse1) / (dfe0 - dfe1)) / (sse1 / dfe1);
                    auto p = fpval(f, dfe0 - dfe1, dfe1);
                    if (p < pval) {
                        idx = j;
                        pval = p;
                        dfe = dfe1;
                        sse = sse1;
                    }
                    ps[j] = p;
                }
            }

            double alpha = par.alpha;
            if (par.mtc == 1)
                alpha /= m;
            else if (par.mtc == 2)
                alpha = alpha * (step + 1) / m;
            else if (par.mtc == 3)
                alpha /= (m - step);

            if (pval > alpha)
                break;

            rsq = 1 - sse / sst;
            if (rsq > par.rsq)
                break;

            in.push_back(idx);
            ignore[idx] = true;

            x.resize(n*q0);
            x.insert(x.end(), xs[idx], xs[idx] + n*c1[idx]);
            q0 += c1[idx];
            dfe0 = dfe;
            sse0 = sse;

            std::cerr << "INFO: forward selection " << step + 1 << ": " << rsq << " " << pval << "\n";
        }
    }

    void backward(const std::vector<double> &x0, const std::vector<double> &x1, const std::vector<int> &c1, const std::vector<double> &y,
        std::vector<double> &ps, std::vector<int> &in)
    {
        int n = y.size();
        int m = c1.size();

        std::vector<const double*> xs;
        xs.push_back(x1.data());
        for (int i = 1; i < m; ++i)
            xs.push_back(xs[i - 1] + n*c1[i - 1]);

        if (ps.empty())
            ps.assign(m, kNaN);

        OLS ols;
        std::vector<double> x, px;

        x.assign(n, 1.0);
        x.insert(x.end(), x0.begin(), x0.end());

        int step = 0;
        while (!in.empty()) {
            x.resize(n + x0.size());
            int q1 = 1 + x0.size() / n;
            for (auto i : in) {
                x.insert(x.end(), xs[i], xs[i] + n*c1[i]);
                q1 += c1[i];
            }

            double dfe1, sse1;
            ols.fit(q1, x, y, dfe1, sse1);

            if (dfe1 == 0.0 || sse1 == 0.0)
                return;

            px.clear();

            for (auto i : in) {
                x.resize(n + x0.size());
                int q0 = 1 + x0.size() / n;
                for (auto j : in) {
                    if (j != i) {
                        x.insert(x.end(), xs[j], xs[j] + n*c1[j]);
                        q0 += c1[j];
                    }
                }

                double dfe0, sse0;
                ols.fit(q0, x, y, dfe0, sse0);

                if (sse0 > sse1 && dfe0 > dfe1) {
                    auto f = ((sse0 - sse1) / (dfe0 - dfe1)) / (sse1 / dfe1);
                    auto p = fpval(f, dfe0 - dfe1, dfe1);
                    px.push_back(p);
                    ps[i] = p;
                }
                else
                    px.push_back(1.0);
            }

            auto itr = std::max_element(px.begin(), px.end());
            if (*itr <= par.alpha)
                break;

            in.erase(in.begin() + (itr - px.begin()));
            std::cerr << "INFO: backward elimination " << ++step << ": " << *itr << "\n";
        }
    }

    int stepwise(bool haploid, const std::vector< std::vector<allele_t> > &gt, const std::vector<int> &gi, const std::vector<double> &y,
        const std::vector< std::vector<double> > &ac, const std::vector< std::vector<double> > &ic, std::vector<double> &ps, std::vector<int> &in)
    {
        std::vector<double> x0;
        for (auto &e : ac)
            x0.insert(x0.end(), e.begin(), e.end());

        std::vector<allele_t> g1;
        std::vector<std::pair<allele_t, allele_t> > g2;
        std::vector<int> idx;
        std::vector<int> c1;
        std::vector<double> x1;

        int m = gt.size();
        for (int j = 0; j < m; ++j) {
            std::vector< std::vector<double> > xg;

            // !!! missing genotype is not allowed !!!

            if (haploid) {
                g1.clear();
                for (auto i : gi)
                    g1.push_back(gt[j][i]);
                design1(factor(g1), xg);
            }
            else {
                g2.clear();
                for (auto i : gi) {
                    auto a = gt[j][2 * i];
                    auto b = gt[j][2 * i + 1];
                    if (a > b)
                        std::swap(a, b);
                    g2.emplace_back(a, b);
                }
                design1(factor(g2), xg);
            }

            if (xg.empty())
                continue;

            idx.push_back(j);

            if (!ic.empty()) {
                int ng = xg.size();
                for (int k = 0; k < ng; ++k) {
                    for (auto v : ic) {
                        std::transform(v.begin(), v.end(), xg[k].begin(), v.begin(), std::multiplies<double>());
                        xg.push_back(v);
                    }
                }
            }

            for (auto &e : xg)
                x1.insert(x1.end(), e.begin(), e.end());
            c1.push_back(xg.size());
        }

        std::vector<double> px;
        forward(x0, x1, c1, y, px, in);
        backward(x0, x1, c1, y, px, in);

        subset(idx, in).swap(in);

        ps.assign(m, kNaN);

        int nx = idx.size();
        for (int i = 0; i < nx; ++i)
            ps[idx[i]] = px[i];

        return 0;
    }

    int fit_model(int t, const std::vector<int> &loci, const Genotype &gt, const std::vector<int> &gi, const Phenotype &pt, const Covariate &ct, std::string &stab, std::string &sest)
    {
        auto y = pt.dat[t];
        int n = y.size();

        std::vector<int> obs;
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(y[i]))
                obs.push_back(i);
        }

        int nv = obs.size();
        if (nv != n)
            subset(y, obs).swap(y);

        std::vector<std::string> env;
        if (!pt.env.empty())
            subset(pt.env, obs).swap(env);

        std::vector<std::string> blk;
        if (!pt.blk.empty())
            subset(pt.blk, obs).swap(blk);

        ANOVA aov;

        if (!env.empty())
            aov.add_main("_ENV_", env);

        if (!blk.empty()) {
            if (env.empty())
                aov.add_main("_BLK_", blk);
            else
                aov.add_nested("_BLK_(_ENV_)", blk, env);
        }

        if (!ct.phe.empty() && !ct.ind.empty()) {
            int m = ct.phe.size();
            for (int i = 0; i < m; ++i) {
                if (nv == n)
                    aov.add_reg(ct.phe[i], ct.dat[i]);
                else
                    aov.add_reg(ct.phe[i], subset(ct.dat[i], obs));
            }
        }

        subset(gi, obs).swap(obs);

        for (auto j : loci) {
            std::vector<std::string> gs;
            if (gt.ploidy == 1) {
                for (auto i : obs) {
                    auto a = gt.dat[j][i];
                    if (a)
                        gs.push_back(gt.allele[j][a - 1]);
                    else
                        gs.push_back("?");
                }
            }
            else {
                for (auto i : obs) {
                    auto a = gt.dat[j][i * 2];
                    auto b = gt.dat[j][i * 2 + 1];
                    if (a > b)
                        std::swap(a, b);
                    if (a && b)
                        gs.push_back(gt.allele[j][a - 1] + '/' + gt.allele[j][b - 1]);
                    else
                        gs.push_back("?/?");
                }
            }

            aov.add_main(gt.loc[j], gs);

            if (!env.empty() && !par.nogxe)
                aov.add_crossed(gt.loc[j] + "*_ENV_", gs, env);
        }

        auto tab = par.sstype == 1 ? aov.solve1(y) : aov.solve3(y);
        auto est = aov.solution(y);

        for (auto &e : tab.p) {
            if (e == 0.0)
                e = std::numeric_limits<double>::min();
        }

        std::ostringstream otab;
        otab << "Source\tDF\tSS\tMS\tF\tp\n";
        for (size_t m = tab.names.size(), i = 0; i < m; ++i)
            otab << tab.names[i] << "\t" << tab.df[i] << "\t" << tab.ss[i] << "\t"
            << tab.ms[i] << "\t" << tab.f[i] << "\t" << tab.p[i] << "\n";
        otab << "Error\t" << tab.error[0] << "\t" << tab.error[1] << "\t" << tab.error[2] << "\n";
        otab << "Total\t" << tab.total[0] << "\t" << tab.total[1] << "\n";
        stab = otab.str();

        std::ostringstream oest;
        oest << "Parameter\tEstimate\n";
        for (size_t m = est.params.size(), i = 0; i < m; ++i)
            oest << est.params[i] << "\t" << est.coeffs[i] << "\n";
        sest = oest.str();

        return 0;
    }

} // namespace

int assoc(int argc, char *argv[])
{
    std::cerr << "assoc (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd("assoc [options]");

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--pheno", "phenotype file", "");
    cmd.add("--covar", "covariate file", "");
    cmd.add("--out", "output file", "assoc.out");
    cmd.add("--alpha", "significance level", "0.05");
    cmd.add("--preselect", "pre-selection threshold", "0.05");
    cmd.add("--mtc", "multiple testing correction, BON/FDR", "");
    cmd.add("--rsq", "maximum model r-square", "0.99");
    cmd.add("--sstype", "sum of squares type", "1");

    cmd.add("--no-gxe", "ignore GxE interaction effect");

    if (argc < 2) {
        cmd.help();
        return 1;
    }

    cmd.parse(argc, argv);

    par.vcf = cmd.get("--vcf");
    par.pheno = cmd.get("--pheno");
    par.covar = cmd.get("--covar");
    par.out = cmd.get("--out");
    par.alpha = std::stod(cmd.get("--alpha"));
    par.preselect = std::stod(cmd.get("--preselect"));
    par.rsq = std::stod(cmd.get("--rsq"));
    par.sstype = std::stoi(cmd.get("--sstype"));

    par.nogxe = cmd.has("--no-gxe");

    if (!cmd.get("--mtc").empty()) {
        auto mtc = cmd.get("--mtc");
        std::transform(mtc.begin(), mtc.end(), mtc.begin(), ::toupper);
        if (mtc == "BON")
            par.mtc = 1;
        else if (mtc == "FDR")
            par.mtc = 2;
        else if (mtc == "HOLM")
            par.mtc = 3;
        else
            std::cerr << "WARNING: invalid command line argument value: " << mtc << "\n";
    }

    Genotype gt;
    Phenotype pt;
    Covariate ct;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::cerr << "INFO: reading phenotype file...\n";
    if (read_pheno(par.pheno, pt) != 0)
        return 1;
    std::cerr << "INFO: " << pt.ind.size() << " individuals, " << pt.phe.size() << " traits\n";

    if (!par.covar.empty()) {
        std::cerr << "INFO: reading covariate file...\n";
        if (read_covar(par.covar, ct) != 0)
            return 1;
        std::cerr << "INFO: " << ct.ind.size() << " individuals, " << ct.phe.size() << " covariates\n";
    }

    std::vector<int> gi;
    merge_data(gt, pt, ct, gi);

    if (gi.empty()) {
        std::cerr << "ERROR: no valid observations are found\n";
        return 1;
    }

    std::vector< std::vector<double> > ac, ic;
    parse_factor(pt, ac, ic);
    if (par.nogxe)
        ic.clear();

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    int m = gt.loc.size();
    int num_traits = pt.phe.size();
    bool haploid = gt.ploidy == 1;

    std::ofstream ofs_loc(par.out + ".loc");
    if (!ofs_loc) {
        std::cerr << "ERROR: can't open file: " << par.out << ".loc" << "\n";
        return 1;
    }

    std::ofstream ofs_aov(par.out + ".aov");
    if (!ofs_aov) {
        std::cerr << "ERROR: can't open file: " << par.out << ".aov" << "\n";
        return 1;
    }

    std::ofstream ofs_est(par.out + ".est");
    if (!ofs_est) {
        std::cerr << "ERROR: can't open file: " << par.out << ".est" << "\n";
        return 1;
    }

    std::vector< std::vector<double> > ps;

    for (int t = 0; t < num_traits; ++t) {
        std::vector<bool> keep;
        for (auto e : pt.dat[t])
            keep.push_back(std::isfinite(e));

        auto gi2 = gi;
        auto ac2 = ac;
        auto ic2 = ic;
        auto y = pt.dat[t];

        if (std::find(keep.begin(), keep.end(), false) != keep.end()) {
            gi2 = subset(gi, keep);
            for (auto &e : ac2)
                e = subset(e, keep);
            for (auto &e : ic2)
                e = subset(e, keep);
            y = subset(y, keep);
        }

        if (y.size() < 10) {
            std::cerr << "WARNING: not enough (<10) valid observations: " << y.size() << "\n";
            continue;
        }

        std::vector<int> loci;
        std::vector<double> p1;

        if (par.preselect > 0.0 && par.preselect < 1.0) {
            std::vector<double> p2;
            glm(haploid, gt.dat, gi2, y, ac, ic, p1, p2);

            std::vector<int> idx;
            for (int j = 0; j < m; ++j) {
                if ((std::isfinite(p1[j]) && p1[j] <= par.preselect) || (!ic.empty() && std::isfinite(p2[j]) && p2[j] <= par.preselect))
                    idx.push_back(j);
            }

            if (!idx.empty()) {
                std::cerr << "INFO: after pre-selection, there are " << idx.size() << " loci\n";
                auto geno = subset(gt.dat, idx);
                stepwise(haploid, geno, gi2, y, ac, ic, p2, loci);
                subset(idx, loci).swap(loci);
                for (size_t k = 0; k < idx.size(); ++k)
                    p1[idx[k]] = p2[k];
            }
        }
        else
            stepwise(haploid, gt.dat, gi2, y, ac, ic, p1, loci);

        ps.push_back(p1);

        ofs_loc << ">" << pt.phe[t] << "\n";
        for (auto j : loci)
            ofs_loc << gt.loc[j] << "\n";
        ofs_loc << "\n";

        std::string tab, est;
        auto s = fit_model(t, loci, gt, gi, pt, ct, tab, est);

        ofs_aov << ">" << pt.phe[t] << "\n" << tab << "\n";
        ofs_est << ">" << pt.phe[t] << "\n" << est << "\n";
    }

    std::ofstream ofs_ps(par.out + ".ps");
    if (!ofs_ps) {
        std::cerr << "ERROR: can't open file: " << par.out << ".ps" << "\n";
        return 1;
    }

    ofs_ps << "Locus\tChromosome\tPosition";
    for (auto &e : pt.phe)
        ofs_ps << "\t" << e;
    ofs_ps << "\n";

    for (int j = 0; j < m; ++j) {
        ofs_ps << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j];
        for (int t = 0; t < num_traits; ++t)
            ofs_ps << "\t" << ps[t][j];
        ofs_ps << "\n";
    }

    return 0;
}
