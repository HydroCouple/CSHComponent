/*!
 *  \file    odesolver.cpp
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 *  Runge Kutta solves from
 *  This file and its associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2018
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#include "stdafx.h"
#include "odesolver.h"

#include <omp.h>
#include <math.h>

#define ODE_TINY 1.0e-30

ODESolver::ODESolver(int maxSize, SolverType solverType)
  :m_solverType(solverType)
{
  m_yscal.resize(maxSize, 0.0);
  m_yerr.resize(maxSize, 0.0);
  m_ytemp.resize(maxSize, 0.0);
  m_ak.resize(maxSize * 5, 0.0);

  switch (m_solverType)
  {
    case RKQS:
      m_solver = &ODESolver::rkqsDriver;
      break;
    default:
      m_solver = &ODESolver::rk4;
      break;
  }
}

ODESolver::~ODESolver()
{
}


int ODESolver::solve(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData)
{
  return (this->*m_solver)(y, n, t, dt, yout, derivs, userData);
}

int ODESolver::rk4(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData)
{
  double tdt, dtt, dt6, *dym, *dyt, *yt, *dydt;

  dym = new double[n]();
  dyt = new double[n]();
  yt = new double[n]();
  dydt = new double[n]();

  dtt = dt * 0.5;
  dt6 = dt / 6.0;

  tdt = t+dtt;

 derivs(t, y, dydt, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++)
  {
    yt[i] = y[i] + dtt * dydt[i]; //First step.
  }


  derivs(tdt, yt, dyt, userData); //Second step.

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++)
  {
    yt[i] = y[i] + dtt * dyt[i];
  }


  derivs(tdt, yt, dym, userData); //Third step.

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++)
  {
    yt[i] = y[i] + dt * dym[i];
    dym[i] += dyt[i];
  }

  derivs(t + dt, yt, dyt, userData); //Fourth step.

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++) //Accumulate increments with proper
  {
    yout[i] = y[i] + dt6 * (dydt[i] + dyt[i] + 2.0 * dym[i]); //weights.
  }

  delete[] yt;
  delete[] dyt;
  delete[] dym;
  delete[] dydt;

  return 0;
}


int ODESolver::rkqsDriver(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void* userData)
{

  double tDid, tNext;
  double t_est = t;
  double dt_est = dt;
  double t_end = t+dt;

  std::fill_n(m_yscal.begin(), n, 0.0);
  std::fill_n(m_yerr.begin(), n, 0.0);
  std::fill_n(m_ytemp.begin(), n, 0.0);
  std::fill_n(m_ak.begin(), 5 * n, 0.0);
  double *dydt = new double[n]();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i= 0; i < n; i++)
  {
    yout[i] = y[i];
  }

  for (int nstp = 1; nstp <= m_maxSteps; nstp++)
  {
    derivs(t_est, yout, dydt, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i= 0; i < n; i++)
    {
      m_yscal[i] = fabs(yout[i]) + fabs(dydt[i] * dt_est) + ODE_TINY;
    }

    if (((t_est + dt_est) - t_end) * (t_est + dt_est - t) > 0.0)
    {
      dt_est = t + dt - t_est;
    }

    if(rkqs(&t_est, dydt, yout, n, dt_est, &tDid, &tNext, derivs, userData))
    {
      break;
    }

    if( (t_est - t_end) * (t_end - t) >= 0.0)
    {
      delete[] dydt;
      return 0;
    }

    if (fabs(tNext) <= 0.0)
    {
      delete[] dydt;
      return 2;
    }

    t_est = tNext;
  }

  delete[] dydt;

  return 3;
}

int ODESolver::rkqs(double *t, double dydt[], double yout[], int n, double dtTry, double *dtDid, double *dtNext, ComputeDerivatives derivs, void* userData)
{
  double err, errmax, dt, dtTemp, tnew, told = *t;

  // --- set initial stepsize
  dt = dtTry;

  for (;;)
  {
      // --- take a Runge-Kutta-Cash-Karp step
      rkck(told, dydt, yout, n, dt, derivs, userData);

      // --- compute scaled maximum error
      errmax = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < n; i++)
      {
          err = fabs(m_yerr[i] / m_yscal[i]);

          if (err > errmax)
            errmax = err;
      }

      errmax /= m_eps;

      // --- error too large; reduce stepsize & repeat
      if (errmax > 1.0)
      {
          dtTemp = m_safety * dt * pow( errmax, m_pshrnk);

          if (dt >= 0)
          {
              if (dtTemp > 0.1 * dt)
                dt = dtTemp;
              else
                dt = 0.1 * dt;
          }
          else
          {
              if (dtTemp < 0.1 * dt)
                dt = dtTemp;
              else
                dt = 0.1 * dt;
          }

          tnew = told + dt;

          if (tnew == told)
            return 2;

          continue;
      }

      // --- step succeeded; compute size of next step
      else
      {
          if (errmax > m_errcon)
            *dtNext = m_safety * dt * pow(errmax,m_pgrow);
          else
            *dtNext = 5.0 * dt;

          *t += (*dtDid = dt);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for (int i = 0; i< n; i++)
          {
            yout[i] = m_ytemp[i];
          }

          return 0;
      }
  }
}

void ODESolver::rkck(double t, double dydt[], double yout[], int n, double dt, ComputeDerivatives derivs, void* userData)
{
  double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
         b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42= -0.9, b43=1.2,
         b51= -11.0/54.0, b52=2.5, b53= -70.0/27.0, b54=35.0/27.0,
         b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
         b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0,
         c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
         dc5= -277.0/14336.0;
  double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
         dc4=c4-13525.0/55296.0, dc6=c6-0.25;

  double *ak2 = &m_ak[0];
  double *ak3 = &m_ak[n];
  double *ak4 = &m_ak[2*n];
  double *ak5 = &m_ak[3*n];
  double *ak6 = &m_ak[4*n];

  for (int i = 0; i<n; i++)
      m_ytemp[i] = yout[i] + b21 * dt * dydt[i];

  derivs(t + a2 * dt, m_ytemp.data(), ak2, userData);

  for (int i = 0; i < n; i++)
      m_ytemp[i] = yout[i] + dt * (b31*dydt[i]+b32*ak2[i]);
  derivs(t + a3 * dt, m_ytemp.data(), ak3, userData);

  for (int i = 0; i<n; i++)
      m_ytemp[i] = yout[i] + dt *(b41*dydt[i]+b42*ak2[i] + b43*ak3[i]);
  derivs(t + a4 * dt, m_ytemp.data(), ak4, userData);

  for (int i = 0; i<n; i++)
      m_ytemp[i] = yout[i] + dt *(b51*dydt[i]+b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
  derivs(t + a5 * dt, m_ytemp.data(), ak5, userData);

  for (int i = 0; i<n; i++)
      m_ytemp[i] = yout[i] + dt * (b61 * dydt[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i]
                 + b65 * ak5[i]);
  derivs(t + a6 * dt, m_ytemp.data(), ak6, userData);

  for (int i = 0; i<n; i++)
      m_ytemp[i] = yout[i] + dt *(c1 * dydt[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);

  for (int i = 0; i<n; i++)
      m_yerr[i] = dt *(dc1 * dydt[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
}
