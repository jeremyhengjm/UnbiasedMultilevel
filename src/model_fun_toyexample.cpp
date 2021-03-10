#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <time.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double inner_product_toyexample(NumericVector x0, NumericVector x1)
{
    double val=0.0;
    int k=0, n=x0.size();
    for(k = 0; k < n; k++)
    {
        val = val + x0[k] * x1[k];
    }
    return val;
}

// [[Rcpp::export]]
double mat_inner_product_toyexample(NumericMatrix mle, NumericVector x0, NumericVector x1)
{
    double val=0.0;
    int k=0, n=mle.nrow();
    for(k = 0; k < n; k++)
    {
        val = val + mle(k, 2) * x0[(int)mle(k, 0)] * x1[(int)mle(k, 1)];
    }
    return val;
}

// [[Rcpp::export]]
double compute_alpha_toyexample(NumericMatrix mle, NumericVector d, NumericVector r)
{
    int k0=0, k1=0, kml=0;
    int nml=mle.nrow(), n=r.size();
    double alpha=0.0, eps=1.0e-6;
    NumericVector y(n);
    for(k0 = 0; k0 < n; k0 ++)
    {
        y[k0] = 0.0;
    }
    for(kml = 0; kml < nml; kml ++)
    {
        y[(int)mle(kml, 0)] = y[(int)mle(kml, 0)] + mle(kml, 2) * d[(int)mle(kml, 1)];
    }
    if(inner_product_toyexample(y, y) > eps)
    {
        alpha = inner_product_toyexample(r, r) / inner_product_toyexample(y, y);
    }
    return alpha;
}

// [[Rcpp::export]]
double compute_beta_toyexample(NumericMatrix mle, NumericVector d, NumericVector rcg)
{
    int k0=0, k1=0, kml=0;
    int nml=mle.nrow(), n=rcg.size();
    double beta=0.0, eps=1.0e-6;
    NumericVector y0(n);
    NumericVector y1(n);
    for(k0 = 0; k0 < n; k0 ++)
    {
        y0[k0] = 0.0;
        y1[k0] = 0.0;
    }
    for(kml = 0; kml < nml; kml ++)
    {
        y0[(int)mle(kml, 0)] = y0[(int)mle(kml, 0)] + mle(kml, 2) * d[(int)mle(kml, 1)];
        y1[(int)mle(kml, 0)] = y1[(int)mle(kml, 0)] + mle(kml, 2) * rcg[(int)mle(kml, 1)];
    }
    if(inner_product_toyexample(y0, y0) > eps)
    {
        beta = inner_product_toyexample(y0, y1) / inner_product_toyexample(y0, y0);
    }
    return beta;
}

// [[Rcpp::export]]
NumericVector sparse_mat_vec_mult_toyexample(NumericMatrix mle, NumericVector x)
{
    int k=0, kml=0;
    int nml=mle.nrow(), n=x.size();
    NumericVector y(n);
    for(k = 0; k < n; k ++)
    {
        y[k] = 0.0;
    }
    for(kml = 0; kml < nml; kml ++)
    {
        y[(int)mle(kml, 0)] += mle(kml, 2) * x[(int)mle(kml, 1)];
    }
    return y;
}

// [[Rcpp::export]]
NumericVector inverse_u_vec_toyexample(NumericMatrix umat, NumericVector x)
{
    int k=0, kml=0;
    int n=x.size();
    NumericVector y(n);
    for(k = 0; k < n; k ++)
    {
        if(k > 0)
        {
            y[n - k - 1] = (x[n - k - 1] - umat(2 * n - 1 - 2 * k, 2) * y[n - k]) / umat(2 * n - 2 - 2 * k, 2);
        }
        else
        {
            y[n - k - 1] = x[n - k - 1] / umat(2 * n - 2 - 2 * k, 2);
        }
    }
    return y;
}

// [[Rcpp::export]]
NumericVector inverse_l_vec_toyexample(NumericMatrix lmat, NumericVector x)
{
    int k=0, kml=0;
    int n=x.size();
    NumericVector y(n);
    for(k = 0; k < n; k ++)
    {
        if(k > 0)
        {
            y[k] = (x[k] - lmat(2 * k - 1, 2) * y[k - 1]) / lmat(2 * k, 2);
        }
        else
        {
            y[k] = x[k] / lmat(2 * k, 2);
        }
    }
    return y;
}

// [[Rcpp::export]]
NumericMatrix sparse_u_factor_toyexample(NumericMatrix mle)
{   int nml=mle.nrow();
    int k=0, n=0, nlu=0;
    int indm = 0;
    n = (int)((nml + 2) / 3);
    nlu = 2 * n - 1;
    NumericMatrix lmat(nlu, 3);
    NumericMatrix umat(nlu, 3);
    for(k = 0; k < n; k ++)
    {
        if(k > 0)
        {
            umat(2 * k - 1, 0) = k - 1;
            umat(2 * k - 1, 1) = k;
            lmat(2 * k - 1, 0) = k;
            lmat(2 * k - 1, 1) = k - 1;
            if (k > 1)
            {
                indm = 2 + 3 * (k - 2) + 1;
            }
            else
            {
                indm = 1;
            }
            umat(2 * k - 1, 2) = mle(indm, 2) / lmat(2 * (k - 1), 2);

            if (k < (n - 1))
            {
                indm = 2 + 3 * (k - 1) + 2;
            }
            else
            {
                indm = 2 + 3 * (k - 1) + 1;
            }
            lmat(2 * k - 1, 2) = mle(indm, 2) / umat(2 * (k - 1), 2);


            umat(2 * k, 0) = k;
            umat(2 * k, 1) = k;
            lmat(2 * k, 0) = k;
            lmat(2 * k, 1) = k;
            umat(2 * k, 2) = sqrt(mle(2 + 3 * (k - 1), 2) - umat(2 * k - 1, 2) * lmat(2 * k - 1, 2));
            lmat(2 * k, 2) = umat(2 * k, 2);
        }
        else
        {
            lmat(0, 0) = 0;
            lmat(0, 1) = 0;
            lmat(0, 2) = sqrt(mle(0, 2));
            umat(0, 0) = 0;
            umat(0, 1) = 0;
            umat(0, 2) = sqrt(mle(0, 2));
        }
    }
    return umat;
}

// [[Rcpp::export]]
NumericMatrix sparse_l_factor_toyexample(NumericMatrix mle)
{   int nml=mle.nrow();
    int k=0, n=0, nlu=0;
    int indm = 0;
    n = (int)((nml + 2) / 3);
    nlu = 2 * n - 1;
    NumericMatrix lmat(nlu, 3);
    NumericMatrix umat(nlu, 3);
    for(k = 0; k < n; k ++)
    {
        if(k > 0)
        {
            umat(2 * k - 1, 0) = k - 1;
            umat(2 * k - 1, 1) = k;
            lmat(2 * k - 1, 0) = k;
            lmat(2 * k - 1, 1) = k - 1;
            if (k > 1)
            {
                indm = 2 + 3 * (k - 2) + 1;
            }
            else
            {
                indm = 1;
            }
            umat(2 * k - 1, 2) = mle(indm, 2) / lmat(2 * (k - 1), 2);

            if (k < (n - 1))
            {
                indm = 2 + 3 * (k - 1) + 2;
            }
            else
            {
                indm = 2 + 3 * (k - 1) + 1;
            }
            lmat(2 * k - 1, 2) = mle(indm, 2) / umat(2 * (k - 1), 2);


            umat(2 * k, 0) = k;
            umat(2 * k, 1) = k;
            lmat(2 * k, 0) = k;
            lmat(2 * k, 1) = k;
            umat(2 * k, 2) = sqrt(mle(2 + 3 * (k - 1), 2) - umat(2 * k - 1, 2) * lmat(2 * k - 1, 2));
            lmat(2 * k, 2) = umat(2 * k, 2);
        }
        else
        {
            lmat(0, 0) = 0;
            lmat(0, 1) = 0;
            lmat(0, 2) = sqrt(mle(0, 2));
            umat(0, 0) = 0;
            umat(0, 1) = 0;
            umat(0, 2) = sqrt(mle(0, 2));
        }
    }
    return lmat;
}

// [[Rcpp::export]]
NumericVector linear_solver_direct_toyexample(NumericMatrix mle, NumericVector d)
{
    int n = d.size();
    int nml = mle.nrow();
    int k0=0;
    NumericVector a(n - 1);  // elements lefto to diagonal
    NumericVector b(n);  // diagonal elements
    NumericVector c(n - 1);  // elements right to diagonal
    NumericVector dp(n);  // updated right part
    NumericVector cp(n - 1);  // updated values of matrix elements right to diagonal
    NumericVector x(n);  // solution vector
    // initialization of the input data for Tridiagonal Matrix Algorithm:
    // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    //
    // notations from the link above are preserved
    for(k0 = 0; k0 < n; k0++)
    {
        dp(k0) = 0.0;
        // elements that are left from the diagonal
        if(k0 < (n - 2))
        {
            a(k0) = mle(2 + 3 * k0 + 2, 2);
            cp(k0) = 0.0;
        }
        if(k0 == (n - 2))
        {
            a(k0) = mle(2 + 3 * k0 + 1, 2);
            cp(k0) = 0.0;

        }

        // elements that are right from the diagonal
        if(0 < k0 && k0 < (n - 1))
        {
            c(k0) = mle(2 + 3 * (k0 - 1) + 1, 2);
        }
        if(k0 == 0)
        {
            c(k0) = mle(1, 2);
        }
        // diagonal elements
        if(k0 > 0)
        {
            b(k0) = mle(2 + 3 * (k0 - 1), 2);
        }
        if(k0 == 0)
        {
            b(k0) = mle(0, 2);
        }
    }
    // cout << "a = " << a << endl;
    // cout << "b = " << b << endl;
    // cout << "c = " << c << endl;
    // calculation of the solution itself
    for(k0 = 0; k0 < (n - 1); k0 ++)
    {
        if(k0 == 0)
        {
            cp(k0) = c(k0) / b(k0);
        }
        if(k0 > 0)
        {
            cp(k0) = c(k0) / (b(k0) - a(k0 - 1) * cp(k0 - 1));
        }
    }
    for(k0 = 0; k0 < n; k0++)
    {
        if(k0 == 0)
        {
            dp(k0) = d(k0) / b(k0);
        }
        if(k0 > 0)
        {
            dp(k0) = (d(k0) - a(k0 - 1) * dp(k0 - 1)) / (b(k0) - a(k0 - 1) * cp(k0 - 1));
        }
    }
    // cout << "cp = " << cp << endl;
    // cout << "dp = " << dp << endl;
    for(k0 = (n - 1); k0 >= 0; k0 = (k0 - 1))
    {
        if(k0 == (n - 1))
        {
            x(k0) = dp(k0);
        }
        else
        {
            x(k0) = dp(k0) - cp(k0) * x(k0 + 1);
        }
    }
    // cout << "x = " << x << endl;
    return x;
}

// [[Rcpp::export]]
NumericVector linear_solver_toyexample(NumericMatrix mle, NumericVector b)
{
    int n = b.size();
    int nml = mle.nrow();
    int k0=0, k1=0, kml=0;
    int kiter=0, niter=50000;
    double factor=0.0, res=0.0, eps=1.0e-10;
    double alpha=0.0, beta=0.0;
    NumericVector x(n);
    NumericVector y(n);
    NumericVector z(n);
    NumericVector r(n);
    NumericVector d(n);
    NumericMatrix lmat(2 * n - 1, 3);
    NumericMatrix umat(2 * n - 1, 3);
    // initial guess
    for(k0 = 0; k0 < n; k0++)
    {
        x[k0] = 0.0;
    }
    for(k0 = 0; k0 < n; k0++)
    {
        r[k0] = 0.0;
        d[k0] = 0.0;
    }
    umat = sparse_u_factor_toyexample(mle);
    lmat = sparse_l_factor_toyexample(mle);
    z = inverse_l_vec_toyexample(lmat, b);
    for(k0 = 0; k0 < n; k0++)
    {
        r[k0] = z[k0];
        d[k0] = 0.0;
    }

    //cout << "z = " << z << endl;
    niter = 5;
    for(kiter = 0; kiter < niter; kiter++)
    {
        y = inverse_u_vec_toyexample(umat, x);
        y = sparse_mat_vec_mult_toyexample(mle, y);
        y = inverse_l_vec_toyexample(lmat, y);
        for(k0 = 0; k0 < n; k0++)
        {
            r[k0] = z[k0] - y[k0];
        }

        d = inverse_u_vec_toyexample(umat, r);
        d = sparse_mat_vec_mult_toyexample(mle, d);
        d = inverse_l_vec_toyexample(lmat, d);

        factor = 1.0 * inner_product_toyexample(d, r) / (inner_product_toyexample(d, d) + eps);
        for(k0 = 0; k0 < n; k0++)
        {
            x[k0] += factor * r[k0];
        }
        res = sqrt(inner_product_toyexample(r, r));
    }
    res = 0.0;
    for(k0 = 0; k0 < n; k0++)
    {
        if(res < abs(r[k0]))
        {
            res = abs(r[k0]);
        }
        // x[k0] += factor * r[k0];
    }

    //cout << "residual = " << res << endl;
    return inverse_u_vec_toyexample(umat, x);
}

// [[Rcpp::export]]
NumericVector cg_solver_toyexample(NumericMatrix mle, NumericVector b)
{
    int n = b.size();
    int nml = mle.nrow();
    int k0=0, k1=0, kml=0;
    int kiter=0, niter=50000;
    double factor=0.0, res=0.0, eps=1.0e-10;
    double alpha=0.0, beta=0.0;
    NumericVector x(n);
    NumericVector r(n);
    NumericVector rcg(n);
    NumericVector d(n);
    // initial guess
    for(k0 = 0; k0 < n; k0++)
    {
        x[k0] = 0.0;
    }
    for(k0 = 0; k0 < n; k0++)
    {
        r[k0] = 0.0;
        d[k0] = 0.0;
    }
    niter = 10;
    for(kiter = 0; kiter < niter; kiter++)
    {
        for(k0 = 0; k0 < n; k0++)
        {
            r[k0] = b[k0];
            d[k0] = 0.0;
        }
        for(kml = 0; kml < nml; kml++)
        {
            r[(int)mle(kml, 0)] = r[(int)mle(kml, 0)] - mle(kml, 2) * x[(int)mle(kml, 1)];
        }
        for(kml = 0; kml < nml; kml++)
        {
            d[(int)mle(kml, 0)] += mle(kml, 2) * r[(int)mle(kml, 1)];
        }
        factor = 1.0 * mat_inner_product_toyexample(mle, r, r) / (inner_product_toyexample(d, d) + eps);
        for(k0 = 0; k0 < n; k0++)
        {
            x[k0] += factor * r[k0];
        }
        res = sqrt(inner_product_toyexample(r, r));
    }
    niter = 0;
    for(k0 = 0; k0 < n; k0 ++)
    {
        r[k0] = b[k0];
        rcg[k0] = 0.0;
    }
    for(kml = 0; kml < nml; kml ++)
    {
        r[(int)mle(kml, 0)] = r[(int)mle(kml, 0)] - mle(kml, 2) * x[(int)mle(kml, 1)];
    }

    for(kml = 0; kml < nml; kml ++)
    {
        rcg[(int)mle(kml, 1)] = rcg[(int)mle(kml, 1)] + mle(kml, 2) * r[(int)mle(kml, 0)];
    }

    for(k0 = 0; k0 < n; k0 ++)
    {
        d[k0] = rcg[k0];
    }

    for(kiter = 0; kiter < niter; kiter++)
    {
        // calculation of alpha
        alpha = compute_alpha_toyexample(mle, d, rcg);
        // update of the solution vector
        for(k0 = 0; k0 < n; k0 ++)
        {
            x[k0] = x[k0] + alpha * d[k0];
            rcg[k0] = 0.0;
        }
        // update of the residual
        for(kml = 0; kml < nml; kml ++)
        {
            r[(int)mle(kml, 0)] = r[(int)mle(kml, 0)] - alpha * mle(kml, 2) * d[(int)mle(kml, 1)];
        }
        for(kml = 0; kml < nml; kml ++)
        {
            rcg[(int)mle(kml, 1)] = rcg[(int)mle(kml, 1)] + mle(kml, 2) * r[(int)mle(kml, 0)];
        }
        // update of the search direction
        beta = compute_beta_toyexample(mle, d, rcg);
        for(k0 = 0; k0 < n; k0++)
        {
            d[k0] = rcg[k0] - beta * d[k0];
        }
        res = sqrt(inner_product_toyexample(r, r));
        //cout << "residual = " << res << endl;
    }
    //cout << "residual = " << res << endl;
    return x;
}


// [[Rcpp::export]]
int compute_matrix_dimension_toyexample(int nx)
{
    return 3 * nx - 2;
}

// [[Rcpp::export]]
NumericMatrix init_system_matrix_toyexample(NumericVector param, int nx)
{
    int kx = 0, kml=0, i0=0, i1=0;
    int nml = 0;
    nml = compute_matrix_dimension_toyexample(nx);
    NumericMatrix mle(nml, 3);
    // boundary element
    kml = 0;
    i0 = 0;
    i1 = 0;
    mle(kml, 0) = i0;
    mle(kml, 1) = i1;
    mle(kml, 2) = 1.0;
    kml = kml + 1;

    i0 = 0;
    i1 = 1;
    mle(kml, 0) = i0;
    mle(kml, 1) = i1;
    mle(kml, 2) = 0.0;
    kml = kml + 1;

    // central elements
    for(kx = 1; kx < (nx - 1); kx++)
    {
        i0 = kx;
        i1 = kx;
        mle(kml, 0) = i0;
        mle(kml, 1) = i1;
        mle(kml, 2) = 2.0;
        kml = kml + 1;

        i0 = kx;
        i1 = kx + 1;
        mle(kml, 0) = i0;
        mle(kml, 1) = i1;
        mle(kml, 2) = -1.0;
        kml = kml + 1;

        i0 = kx;
        i1 = kx - 1;
        mle(kml, 0) = i0;
        mle(kml, 1) = i1;
        mle(kml, 2) = -1.0;
        kml = kml + 1;
    }

    i0 = nx - 1;
    i1 = nx - 1;
    mle(kml, 0) = i0;
    mle(kml, 1) = i1;
    mle(kml, 2) = 1.0;
    kml = kml + 1;

    i0 = nx - 1;
    i1 = nx - 2;
    mle(kml, 0) = i0;
    mle(kml, 1) = i1;
    mle(kml, 2) = 0.0;
    kml = kml + 1;

    return mle;
}

// [[Rcpp::export]]
NumericVector init_right_part_toyexample(NumericVector param, int nx)
{
    int k0=0, k1=0, kx=0;
    double dh = 2.0 * 3.1415926 / (nx - 1);
    double val0=0.0, val1=0.0;
    double scaling_factor = 1.0e0;
    NumericVector b(nx);
    for(kx = 0; kx < nx; kx++)
    {
        val0 = (2.0 * sin(2.0 * dh * kx) - sin(2.0 * dh * (kx + 1)) - sin(2.0 * dh * (kx - 1))) / 4.0;
        val1 = 2.0 * sin(dh * kx) - sin(dh * (kx + 1)) - sin(dh * (kx - 1));
        b[kx] = param[0] * val0 + param[1] * val1;
        // b[kx] = scaling_factor * b[kx];
    }
    b[0] = 0.0;
    b[nx - 1] = 0.0;
    return b;
}

// [[Rcpp::export]]
NumericVector numerical_solution_toyexample(NumericVector param, int nx)
{
    int nml = compute_matrix_dimension_toyexample(nx);
    NumericVector x(nx);
    NumericVector b(nx);
    NumericMatrix mle(nml, 3);

    mle = init_system_matrix_toyexample(param, nx);
    b = init_right_part_toyexample(param, nx);
    x = linear_solver_toyexample(mle, b);
    //x = linear_solver_direct_toyexample(mle, b);
    return x;
}

// [[Rcpp::export]]
NumericVector exact_solution_toyexample(NumericVector param, int nx)
{
    int kx=0;
    double dh = 2.0 * 3.1415926 / (nx - 1);
    double scaling_factor = 1.0e0;
    NumericVector x(nx);
    for(kx = 0; kx < nx; kx ++)
    {
        x[kx] = param[0] * sin(2.0 * dh * kx) / 4.0 + param[1] * sin(dh * kx);
        // x[kx] = scaling_factor * x[kx];
    }
    return x;
}

// [[Rcpp::export]]
NumericVector numerical_solution_values_toyexample(NumericVector x, int nval)
{
    int nx = x.size();
    int kx=0, kval=0;
    double dh = 2.0 * 3.1415926 / (nx - 1);
    double dhval = (2.0 * 3.1415926 - 0.0000000000001) / (nval - 1);
    double frac=0.0;
    NumericVector xval(nval);
    for(kval = 0; kval < nval; kval ++)
    {
        kx = (int) (kval * dhval / dh);
        frac = (kval * dhval / dh) - kx;
        xval[kval] = (1.0 - frac) * x[kx] + frac * x[kx + 1];
    }
    return xval;
}


// [[Rcpp::export]]
double l2_norm_toyexample(NumericVector param, int l)
{
    int kval = 0;
    int nval = 50000;
    int nx = 0;
    double res = 0.0;
    nx = 1 + (int)exp(1.0 * (1.0 * l + 3) * log(2.0));
    NumericVector x(nx);
    x = numerical_solution_toyexample(param, nx);
    // cout << "x = " << x << endl;
    NumericVector xnumer(nval);
    NumericVector xexact(nval);
    xnumer = numerical_solution_values_toyexample(x, nval);
    xexact = exact_solution_toyexample(param, nval);
    for(kval = 0; kval < nval; kval ++)
    {
        res += (xnumer[kval] - xexact[kval]) * (xnumer[kval] - xexact[kval]) / nval;
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix init_system_matrix_derivative_toyexample(NumericVector param, NumericVector pgrad, int nx)
{
    int kx = 0, kml=0, i0=0, i1=0;
    int nml = 0;
    nml = compute_matrix_dimension_toyexample(nx);
    NumericMatrix mle_grad(nml, 3);
    // boundary element
    kml = 0;
    i0 = 0;
    i1 = 0;
    mle_grad(kml, 0) = i0;
    mle_grad(kml, 1) = i1;
    mle_grad(kml, 2) = 0.0;
    kml = kml + 1;

    i0 = 0;
    i1 = 1;
    mle_grad(kml, 0) = i0;
    mle_grad(kml, 1) = i1;
    mle_grad(kml, 2) = 0.0;
    kml = kml + 1;

    // central elements
    for(kx = 1; kx < (nx - 1); kx++)
    {
        i0 = kx;
        i1 = kx;
        mle_grad(kml, 0) = i0;
        mle_grad(kml, 1) = i1;
        mle_grad(kml, 2) = 0.0;
        kml = kml + 1;

        i0 = kx;
        i1 = kx + 1;
        mle_grad(kml, 0) = i0;
        mle_grad(kml, 1) = i1;
        mle_grad(kml, 2) = 0.0;
        kml = kml + 1;

        i0 = kx;
        i1 = kx - 1;
        mle_grad(kml, 0) = i0;
        mle_grad(kml, 1) = i1;
        mle_grad(kml, 2) = 0.0;
        kml = kml + 1;
    }

    i0 = nx - 1;
    i1 = nx - 1;
    mle_grad(kml, 0) = i0;
    mle_grad(kml, 1) = i1;
    mle_grad(kml, 2) = 0.0;
    kml = kml + 1;

    i0 = nx - 1;
    i1 = nx - 2;
    mle_grad(kml, 0) = i0;
    mle_grad(kml, 1) = i1;
    mle_grad(kml, 2) = 0.0;
    kml = kml + 1;

    return mle_grad;
}

// [[Rcpp::export]]
NumericVector init_right_part_derivative_toyexample(NumericVector param, NumericVector pgrad, int nx)
{
    int k0=0, k1=0, kx=0;
    double dh = 2.0 * 3.1415926 / (nx - 1);
    double val0=0.0, val1=0.0;
    double scaling_factor = 1.0e0;
    NumericVector b(nx);
    for(kx = 0; kx < nx; kx++)
    {
        val0 = (2.0 * sin(2.0 * dh * kx) - sin(2.0 * dh * (kx + 1)) - sin(2.0 * dh * (kx - 1))) / 4.0;
        val1 = 2.0 * sin(dh * kx) - sin(dh * (kx + 1)) - sin(dh * (kx - 1));
        b[kx] = pgrad[0] * val0 + pgrad[1] * val1;
        // b[kx] = scaling_factor * b[kx];
    }
    b[0] = 0.0;
    b[nx - 1] = 0.0;
    return b;
}

// [[Rcpp::export]]
NumericVector numerical_solution_derivative_toyexample(NumericVector param, NumericVector pgrad, int nx)
{
    int kx=0, kml=0;
    int nml = 0;
    nml = compute_matrix_dimension_toyexample(nx);
    NumericVector x(nx);
    NumericVector xgrad(nx);
    NumericVector b(nx);
    NumericVector bgrad(nx);
    NumericMatrix mle(nml, 2);
    NumericMatrix mle_grad(nml, 2);
    b = init_right_part_toyexample(param, nx);
    bgrad = init_right_part_derivative_toyexample(param, pgrad, nx);
    mle = init_system_matrix_toyexample(param, nx);
    mle_grad = init_system_matrix_derivative_toyexample(param, pgrad, nx);
    for(kml = 0; kml < nml; kml ++)
    {
        bgrad[(int)mle_grad(kml, 0)] -= mle_grad(kml, 2) * x[(int)mle_grad(kml, 1)];
    }
    xgrad = linear_solver_toyexample(mle, bgrad);
    //xgrad = linear_solver_direct_toyexample(mle, bgrad);
    return xgrad;
}

// [[Rcpp::export]]
NumericVector observation_points_toyexample(int nvec)
{
    int kx=0, nx=nvec;
    double dh = 2.0 * 3.1415926 / nx;
    NumericVector x(nx);
    for(kx = 0; kx < nx; kx ++)
    {
        x[kx] = (0.5 + kx) * dh;
    }
    return x;
}

// [[Rcpp::export]]
NumericVector observation_toyexample(NumericVector param, int l)
{
    int kvec=0, nvec=50;
    int nx=0;
    int ind0=0;
    double w0=0.0;
    double eps=1.0e-3;
    double dh = 2.0 * 3.1415926;
    double scaling_factor = 1.0e0;
    nx = 1 + (int)(exp((l + 3.0) * log(2.0)) + eps);
    dh = dh / (nx - 1);
    NumericVector xvec(nvec);
    NumericVector xeval(nvec);
    NumericVector x(nx);
    x = numerical_solution_toyexample(param, nx);
    xeval = observation_points_toyexample(nvec);
    for(kvec = 0; kvec < nvec; kvec ++)
    {
        ind0 = (int)(xeval(kvec) / dh);
        w0 = xeval(kvec) / dh - ind0;
        xvec(kvec) = (1.0 - w0) * x[ind0] + w0 * x[ind0 + 1];
    }
    return xvec;
}

// [[Rcpp::export]]
NumericVector observation_grad_val_toyexample(NumericVector param, int l, int m)
{
    int kvec=0, nvec=50;
    int nx=0;
    int ind0=0;
    double w0=0.0;
    double eps=1.0e-3;
    double dh = 2.0 * 3.1415926;
    double scaling_factor = 1.0e0;
    nx = 1 + (int)(exp((l + 3.0) * log(2.0)) + eps);
    dh = dh / (nx - 1);
    NumericVector xvec(nvec);
    NumericVector xeval(nvec);
    NumericVector x(nx);
    NumericVector pgrad(param.size());
    for(kvec = 0; kvec < param.size(); kvec ++)
    {
        pgrad[kvec] = 0.0;
    }
    pgrad[m] = 1.0;
    x = numerical_solution_derivative_toyexample(param, pgrad, nx);
    xeval = observation_points_toyexample(nvec);
    for(kvec = 0; kvec < nvec; kvec ++)
    {
        ind0 = (int)(xeval(kvec) / dh);
        w0 = xeval(kvec) / dh - ind0;
        xvec(kvec) = (1.0 - w0) * x[ind0] +  w0 * x[ind0 + 1];
    }
    return xvec;
}


// [[Rcpp::export]]
NumericMatrix observation_grad_toyexample(NumericVector param, int l)
{
    int kvec=0, nvec=50;
    int k0=0, k1=0;
    int nx=0;
    int nparam = param.size();
    double eps=1.0e-3;
    NumericVector pgrad(nvec);
    NumericMatrix pgrad_mat(nvec, nparam);
    for(kvec = 0; kvec < nvec; kvec ++)
    {
        pgrad[kvec] = 0.0;
    }
    for(k0 = 0; k0 < nparam; k0 ++)
    {
        for (k1 = 0; k1 < nvec; k1 ++)
        {
            pgrad_mat(k1, k0) = 0.0;
        }
    }
    for(k0 = 0; k0 < nparam; k0 ++)
    {
        pgrad = observation_grad_val_toyexample(param, l, k0);
        for(k1 = 0; k1 < nvec; k1 ++)
        {
            pgrad_mat(k1, k0) = pgrad(k1);
        }
    }
    return pgrad_mat;
}
