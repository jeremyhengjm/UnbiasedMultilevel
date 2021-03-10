#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector flow_val_covid19(NumericVector x, NumericVector u)
{
    int nx = x.size();
    int nu = u.size();
    int kx = 0;
    double xsum = 0.0;
    double alpha = 0.775, beta = 0.125;
    // v - vector of flow direction
    NumericVector v(nx);
    for(kx = 0; kx < (nx - 1); kx ++)
    {
        xsum = xsum + x[kx];
    }

    // SIR-X - model
    v[0] = - alpha * x[0] * x[1] - u[0] * x[0];
    v[1] = alpha * x[0] * x[1] - beta * x[1] - u[0] * x[1] - u[1] * x[1];
    v[2] = (u[0] + u[1]) * x[1];
    v[3] = alpha * x[0] * x[1];
    return v;
}

// [[Rcpp::export]]
NumericMatrix flow_grad_x0_covid19(NumericVector x, NumericVector u)
{
    int nx = x.size();
    int nu = u.size();
    int kx = 0;
    NumericMatrix fmat(nx, nx);
    double xsum = 0.0;
    double alpha = 0.775, beta = 0.125;
    for(kx = 0; kx < (nx - 1); kx ++)
    {
        xsum = xsum + x[kx];
    }
    // derivatives of the right part
    // the first component is for the coordinate which is differentiated
    // the second component states for the coordinate with respect to which
    // differentiation is performed
    fmat(0, 0) = - alpha * x[1] - u[0];
    fmat(0, 1) = - alpha * x[0];
    fmat(0, 2) = 0.0;
    fmat(0, 3) = 0.0;

    fmat(1, 0) = alpha * x[1];
    fmat(1, 1) = alpha * x[0] - beta - u[0] - u[1];
    fmat(1, 2) = 0.0;
    fmat(1, 3) = 0.0;

    fmat(2, 0) = 0.0;
    fmat(2, 1) = u[0] + u[1];
    fmat(2, 2) = 0.0;
    fmat(2, 3) = 0.0;

    fmat(3, 0) = alpha * x[1];
    fmat(3, 1) = alpha * x[0];
    fmat(3, 2) = 0.0;
    fmat(3, 3) = 0.0;

    return fmat;
}

// [[Rcpp::export]]
NumericMatrix flow_grad_param_covid19(NumericVector x, NumericVector u)
{
    int nx = x.size();
    int nu = u.size();

    int kx = 0;
    NumericMatrix fmat(nx, nu);
    double xsum = 0.0;
    double alpha = 0.775, beta = 0.125;
    for(kx = 0; kx < (nx - 1); kx ++)
    {
        xsum = xsum + x[kx];
    }

    // derivatives of the right part
    // the first component is for the coordinate which is differentiated
    // the second component states for the parameter coordinate with respect to which
    // differentiation is performed
    fmat(0, 0) = -x[0];
    fmat(0, 1) = 0.0;

    fmat(1, 0) = -x[1];
    fmat(1, 1) = -x[1];

    fmat(2, 0) = x[1];
    fmat(2, 1) = x[1];

    fmat(3, 0) = 0.0;
    fmat(3, 1) = 0.0;
    return fmat;
}

// [[Rcpp::export]]
NumericVector flow_vec(NumericVector x, NumericVector u)
{
    int nxu = x.size();
    int nu = u.size();
    int nx = 4;
    int k0=0, k1=0, k2=0;
    NumericVector xval(nx);
    NumericMatrix fmat_x0(nx, nx);
    NumericMatrix fmat_u0(nx, nu);
    NumericVector fval(nx);
    NumericVector f(nxu);
    for(k0 = 0; k0 < nx; k0 ++)
    {
        xval[k0] = x[k0];
    }
    for(k0 = 0; k0 < nxu; k0 ++)
    {
        f[k0] = 0.0;
    }

    fmat_x0 = flow_grad_x0_covid19(xval, u);
    fmat_u0 = flow_grad_param_covid19(xval, u);
    fval = flow_val_covid19(xval, u);
    for(k0 = 0; k0 < nx; k0 ++)
    {
        f[k0] = fval[k0];
    }
    // Derivative with respect to model parameters
    for(k0 = 0; k0 < nx; k0 ++)
    {
        for(k1 = 0; k1 < nu; k1 ++)
        {
            f[nx + k0 * nu + k1] = fmat_u0(k0, k1);
            for(k2 = 0; k2 < nx; k2 ++)
            {
                f[nx + k0 * nu + k1] += fmat_x0(k0, k2) * x[nx + k2 * nu + k1];
            }
        }
    }

    // Derivative with respect to initial conditions
    for(k0 = 0; k0 < nx; k0 ++)
    {
        for(k1 = 0; k1 < nx; k1 ++)
        {
            for(k2 = 0; k2 < nx; k2 ++)
            {
                f[nx * (nu + 1) + k0 * nx + k1] += fmat_x0(k0, k2) * x[nx * (nu + 1) + k2 * nx + k1];
            }
        }
    }
    return f;
}

// [[Rcpp::export]]
NumericVector Single_Time_Step_Explicit_covid19(double dt, NumericVector x,
                                        NumericVector u,
                                        NumericVector wvec,
                                        NumericMatrix dmat)
{
    int nx = x.size();
    int nw = wvec.size();
    int k0=0, k1=0, k2=0;
    NumericMatrix fmat(nw, nx);
    NumericMatrix xmat(nw, nx);
    NumericVector f(nx);
    NumericVector dx(nx);
    NumericVector xup(nx);

    for(k0 = 0; k0 < nx; k0 ++)
    {
        xup[k0] = 0.0;
        dx[k0] = 0.0;
        for(k1 = 0; k1 < nw; k1 ++)
        {
            fmat(k1, k0) = 0.0;
            xmat(k1, k0) = x[k0];
        }
    }

    // calculation of the fmat;
    for(k0 = 0; k0 < nw; k0 ++)
    {
        for(k1 = 0; k1 < k0; k1 ++)
        {
            for(k2 = 0; k2 < nx; k2 ++)
            {
                xmat(k0, k2) = xmat(k0, k2) + dt * dmat(k0, k1) * fmat(k1, k2);
            }
        }
        for(k2 = 0; k2 < nx; k2 ++)
        {
            xup[k2] = xmat(k0, k2);
        }
        f = flow_vec(xup, u);
        for(k2 = 0; k2 < nx; k2 ++)
        {
            fmat(k0, k2) = f[k2];
        }
    }

    for(k1 = 0; k1 < nx; k1 ++)
    {
        xup[k1] = x[k1];
    }

    for(k0 = 0; k0 < nw; k0 ++)
    {
        for(k1 = 0; k1 < nx; k1 ++)
        {
            xup[k1] = xup[k1] + wvec[k0] * dt * fmat(k0, k1);
        }
    }
    // Demo for naive solver
    // f = flow_vec(x, u);
    // for(k0 = 0; k0 < nx; k0 ++)
    // {
    //     xup[k0] = x[k0] + dt * f[k0];
    // }
    return xup;
}

// [[Rcpp::export]]
NumericVector init_weight_vector_covid19()
{
    NumericVector w(4);
    w[0] = 1.0 / 6.0;
    w[1] = 1.0 / 3.0;
    w[2] = 1.0 / 3.0;
    w[3] = 1.0 / 6.0;
    return w;
}

// [[Rcpp::export]]
NumericMatrix init_d_matrix_covid19()
{
    int k0=0, k1=0;
    NumericMatrix dmat(4, 4);
    for(k0 = 0; k0 < dmat.nrow(); k0 ++)
    {
        for(k1 = 0; k1 < dmat.ncol(); k1 ++)
        {
            dmat(k0, k1) = 0.0;
        }
    }
    dmat(0, 0) = 0.0;
    dmat(1, 0) = 1.0 / 2.0;
    dmat(2, 1) = 1.0 / 2.0;
    dmat(3, 2) = 1.0;
    return dmat;
}

// [[Rcpp::export]]
NumericVector init_weight_vector_high_order_covid19()
{
    NumericVector w(6);
    w[0] = 16.0 / 135.0;
    w[1] = 0.0;
    w[2] = 6656.0 / 12825.0;
    w[3] = 28561.0 / 56430.0;
    w[4] = -9.0 / 50.0;
    w[5] = 2.0 / 55.0;
    return w;
}

// [[Rcpp::export]]
NumericMatrix init_d_matrix_high_order_covid19()
{
    int k0=0, k1=0;
    NumericMatrix dmat(6, 6);
    for(k0 = 0; k0 < dmat.nrow(); k0 ++)
    {
        for(k1 = 0; k1 < dmat.ncol(); k1 ++)
        {
            dmat(k0, k1) = 0.0;
        }
    }
    dmat(0, 0) = 0.0;
    dmat(1, 0) = 1.0 / 4.0;
    dmat(2, 0) = 3.0 / 32.0;
    dmat(2, 1) = 9.0 / 32.0;
    dmat(3, 0) = 1932.0 / 2197.0;
    dmat(3, 1) = -7200.0 / 2197.0;
    dmat(3, 2) = 7296.0 / 2197.0;
    dmat(4, 0) = 439.0 / 216.0;
    dmat(4, 1) = -8.0;
    dmat(4, 2) = 3680.0 / 513.0;
    dmat(3, 3) = -845.0 / 4104.0;
    dmat(5, 0) = -8.0 / 27.0;
    dmat(5, 1) = 2.0;
    dmat(5, 2) = -3544.0 / 2565.0;
    dmat(5, 3) = 1859.0 / 4104.0;
    dmat(5, 4) = -11.0 / 40.0;
    return dmat;
}

// [[Rcpp::export]]
NumericVector forward_simulation_covid19(NumericVector x0, NumericVector u, double tf,
                                  double tstep)
{
    int nx = x0.size();
    int k0=0, k1=0;
    int nt = int(tf / tstep);
    double dt;
    NumericVector wvec;
    NumericMatrix dmat;
    wvec = init_weight_vector_covid19();
    dmat = init_d_matrix_covid19();
    NumericVector x1(nx);

    for(k0 = 0 ; k0 < nx; k0 ++)
    {
        x1[k0] = x0[k0];
    }
    for(k0 = 0 ; k0 < nt; k0 ++)
    {
        dt = tstep;
        x1 = Single_Time_Step_Explicit_covid19(dt, x1, u, wvec, dmat);
    }
    return x1;
}

// [[Rcpp::export]]
NumericMatrix forward_simulation_output_covid19(NumericVector param, int l)
{
    double t0 = 0.0, tf = 52.0;
    double npp = 6.665e7;
    int neval = 52, nstep = 2;
    int nx = 4, nu=(param.size() - 1);
    int nxu = nx * (1 + nu + nx);
    int k0 = 0, k1 = 0;
    NumericVector x0(nx);
    NumericVector x1(nxu);
    NumericVector u(nu);
    for(k0 = 0; k0 < nu; k0 ++)
    {
        u[k0] = param[k0];
    }
    for(k0 = 0; k0 < l; k0 ++)
    {
        nstep = 2 * nstep;
    }
    NumericMatrix xmat(neval + 1, nxu);
    NumericMatrix xout(neval + 1, 2 + nu);
    for(k0 = 0; k0 < nx; k0 ++)
    {
        x0[k0] = 0.0;
    }
    x0[1] = 1.0 / npp;
    x0[0] = 1.0 - x0[1];
    x0[nx - 1] = x0[1];
    // x0[3] = x0[3] + 0.001 / npp;
    for(k0 = 0; k0 < xmat.nrow(); k0 ++)
    {
        for(k1 = 0; k1 < xmat.ncol(); k1 ++)
        {
            xmat(k0, k1) = 0.0;
        }
    }
    // initial condition
    for(k0 = 0; k0 < nx; k0 ++)
    {
        x1[k0] = x0[k0];
        xmat(0, k0) = x0[k0];
        xmat(0, nx * (nu + 1) + k0 * nx + k0) = 1.0;
        x1[nx * (nu + 1) + k0 * nx + k0] = 1.0;
    }
    for(k0 = 0; k0 < neval; k0 ++)
    {
        x1 = forward_simulation_covid19(x1, u, (t0 + param[nu]) / neval, (t0 + param[nu]) / neval / nstep);
    }

    for(k0 = 0; k0 < nxu; k0 ++)
    {
        xmat(0, k0) = x1[k0];
    }

    for(k0 = 0; k0 < neval; k0 ++)
    {
        x1 = forward_simulation_covid19(x1, u, tf / neval, tf / neval / nstep);
        // cout << "xmat = " << xmat << "\n";
        // cout << "\n";
        for(k1 = 0; k1 < nxu; k1 ++)
        {
            xmat(k0 + 1, k1) = x1[k1];
        }
    }
    for(k0 = 0; k0 < (neval + 1); k0 ++)
    {
        xout(k0, 0) = xmat(k0, 3);
        xout(k0, 1) = xmat(k0, 10);
        xout(k0, 2) = xmat(k0, 11);
        for(k1 = 0; k1 < nxu; k1 ++)
        {
            x1[k1] = xmat(k0, k1);
        }
        x1 = flow_vec(x1, u);
        xout(k0, 3) = x1[3];
    }
    // cout << "xmat = " << xmat << "\n";
    // return xmat;
    return xout;
}

// [[Rcpp::export]]
NumericVector observation_covid19(NumericVector u, int l)
{
    int neval = 52;
    int nu=u.size();
    int k0 = 0;
    NumericMatrix xmat(1 + neval, 2 + nu);
    NumericVector x(neval + 1);
    xmat = forward_simulation_output_covid19(u, l);
    for(k0 = 0; k0 < (neval + 1); k0 ++)
    {
        x(k0) = xmat(k0, 0);
    }
    return x;
}

// [[Rcpp::export]]
NumericMatrix observation_grad_covid19(NumericVector u, int l)
{
    int neval = 52;
    int nu=u.size();
    int k0 = 0, k1 = 0;
    NumericMatrix xmat(1 + neval, 2 + nu);
    NumericMatrix x(neval + 1, nu);
    xmat = forward_simulation_output_covid19(u, l);
    for(k0 = 0; k0 < (neval + 1); k0 ++)
    {
        for(k1 = 0; k1 < nu; k1 ++)
        {
            x(k0, k1) = xmat(k0, 1 + k1);
        }
    }
    long long int current_number=0;
    return x;
}




