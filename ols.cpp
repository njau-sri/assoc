#include <cmath>
#include <algorithm>
#include "ols.h"

#ifdef BLAS_INT64
#include <cstddef>
using bint = std::ptrdiff_t;
static_assert(sizeof(bint) == 8, "BLAS integer size is not 8");
#else
using bint = int;
static_assert(sizeof(bint) == 4, "BLAS integer size is not 4");
#endif

extern "C"
{
    double dnrm2_(bint *n, const double *x, bint *incx);

    void dgemv_(char *trans, bint *m, bint *n, double *alpha, const double *a, bint *lda,
        const double *x, bint *incx, double *beta, double *y, bint *incy);

    void dgemm_(char *transa, char *transb, bint *m, bint *n, bint *k, double *alpha, const double *a, bint *lda,
        const double *b, bint *ldb, double *beta, double *c, bint *ldc);

    void dgels_(char *trans, bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb,
        double *work, bint *lwork, bint *info);

    void dgelsy_(bint *m, bint *n, bint *nrhs, double *a, bint *lda, double *b, bint *ldb, bint *jpvt,
        double *rcond, bint *rank, double *work, bint *lwork, bint *info);

    void dgeqp3_(bint *m, bint *n, double *a, bint *lda, bint *jpvt, double *tau, double *work, bint *lwork, bint *info);

    void dorgqr_(bint *m, bint *n, bint *k, double *a, bint *lda, const double *tau, double *work, bint *lwork, bint *info);

} // extern C

namespace
{
    double nrm2(bint n, const double *x, bint incx)
    {
        return dnrm2_(&n, x, &incx);
    }

    int gemv(char trans, bint m, bint n, double alpha, const double *a, bint lda, const double *x, bint incx, double beta, double *y, bint incy)
    {
        dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
        return 0;
    }

    int gemm(char transa, char transb, bint m, bint n, bint k, double alpha, const double *a, bint lda, const double *b, bint ldb, double beta, double *c, bint ldc)
    {
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
        return 0;
    }

    int gels(char trans, bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb)
    {
        bint info = 0;

        double wkopt;
        bint lwork = -1;
        dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info);

        if (info == 0) {
            lwork = static_cast<bint>(wkopt);
            double *work = new double[lwork];
            dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
            delete[] work;
        }

        return info;
    }

    int gelsy(bint m, bint n, bint nrhs, double *a, bint lda, double *b, bint ldb, bint *jpvt, double rcond, bint *rank)
    {
        bint info = 0;

        double wkopt;
        bint lwork = -1;
        dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &wkopt, &lwork, &info);

        if (info == 0) {
            lwork = static_cast<bint>(wkopt);
            double *work = new double[lwork];
            dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
            delete[] work;
        }

        return info;
    }

    int geqp3(bint m, bint n, double *a, bint lda, bint *jpvt, double *tau)
    {
        bint info = 0;

        double wkopt;
        bint lwork = -1;
        dgeqp3_(&m, &n, a, &lda, jpvt, tau, &wkopt, &lwork, &info);

        if (info == 0) {
            lwork = static_cast<bint>(wkopt);
            double *work = new double[lwork];
            dgeqp3_(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
            delete[] work;
        }

        return info;
    }

    int orgqr(bint m, bint n, bint k, double *a, bint lda, const double *tau)
    {
        bint info = 0;

        double wkopt;
        bint lwork = -1;
        dorgqr_(&m, &n, &k, a, &lda, tau, &wkopt, &lwork, &info);

        if (info == 0) {
            lwork = static_cast<bint>(wkopt);
            double *work = new double[lwork];
            dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
            delete[] work;
        }

        return info;
    }

} // namespace

int OLS::fit(int p, const std::vector<double> &x, const std::vector<double> &y, double &dfe, double &sse)
{
    static const double tol = 1e-8;

    int info = 0;
    int n = y.size();
    int ldy = std::max(n, p);
    int rank = std::min(n, p);

    dfe = sse = 0.0;

    x_ = x;
    y_.assign(ldy, 0.0);
    std::copy(y.begin(), y.end(), y_.begin());

    bool singular = true;

    if (n >= p) {
        singular = false;

        info = gels('N', n, p, 1, x_.data(), n, y_.data(), ldy);

        if (info == 0) {
            for (int i = 0; i < p; ++i) {
                if (std::fabs(x_[i*n + i]) < tol) {
                    singular = true;
                    break;
                }
            }
        }
        else
            singular = true;
    }

    if (singular) {
        x_ = x;
        y_.assign(ldy, 0.0);
        std::copy(y.begin(), y.end(), y_.begin());

        bint lrank = 0;
        std::vector<bint> jpvt(p, 0);
        double rcond = tol;

        info = gelsy(n, p, 1, x_.data(), n, y_.data(), ldy, jpvt.data(), rcond, &lrank);

        rank = lrank;
    }

    b_.assign(y_.begin(), y_.begin() + p);

    if (rank < n) {
        dfe = n - rank;
        y_ = y;
        gemv('N', n, p, -1.0, x.data(), n, b_.data(), 1, 1.0, y_.data(), 1);
        sse = nrm2(n, y_.data(), 1);
        sse *= sse;
    }

    return 0;
}

int OLS::fit(int p, int q, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> zt, double &dfe, double &sse)
{
    static const double tol = 1e-8;

    int n = y.size();
    int m = std::min(p, q);

    std::vector<bint> jpvt(q, 0);
    std::vector<double> tau(m);

    geqp3(p, q, zt.data(), p, jpvt.data(), tau.data());

    int k = 0;
    for (int i = 0; i < m; ++i) {
        if (std::fabs(zt[i*p + i]) > tol)
            ++k;
    }

    std::vector<double> Q(p*p);
    std::copy_n(zt.begin(), p*m, Q.begin());
    orgqr(p, p, m, Q.data(), p, tau.data());

    auto Q0 = Q.data() + p*k;

    std::vector<double> X(n*(p - k));
    gemm('N', 'N', n, p - k, p, 1.0, x.data(), n, Q0, p, 0.0, X.data(), n);

    fit(p - k, X, y, dfe, sse);

    std::vector<double> b(p);
    gemv('N', p, p - k, 1.0, Q0, p, b_.data(), 1, 0.0, b.data(), 1);

    b_ = b;

    return 0;
}
