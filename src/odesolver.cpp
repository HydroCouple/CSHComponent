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

#include <cvode/cvode.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_serial.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <math.h>

#define ODE_TINY 1.0e-30

ODESolver::ODESolver(int size, SolverType solverType)
  : m_size(size),
    m_maxSteps(10000),
    m_order(6),
    m_safety(0.9),
    m_pgrow(-0.2),
    m_pshrnk(-0.25),
    m_errcon(1.89e-4),
    m_relTol(1e-4),
    m_absTol(1e-8),
    m_yscal(nullptr),
    m_yerr(nullptr),
    m_ytemp(nullptr),
    m_ak(nullptr),
    m_solverType(solverType),
    #ifdef USE_CVODE
    m_cvodeSolver(nullptr),
    m_cvy (nullptr)
  #endif
{
}

ODESolver::~ODESolver()
{
  clearMemory();
}

void ODESolver::initialize()
{
  clearMemory();

  switch (m_solverType)
  {
    case RKQS:
      {
        m_solver = &ODESolver::rkqsDriver;
        m_yscal = new double[m_size]();
        m_yerr = new double[m_size]();
        m_ytemp = new double[m_size]();
        m_ak = new double[m_size * 5]();
      }
      break;
#ifdef  USE_CVODE
    case CVODE_ADAMS:
      {

        m_cvodeSolver = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
        m_solver = &ODESolver::solveCVODE;

#ifdef USE_CVODE_OPENMP
        m_cvy = N_VNew_OpenMP(m_size, omp_get_max_threads());
#else
        m_cvy = N_VNew_Serial(m_size);
#endif

        CVodeInit(m_cvodeSolver, &ODESolver::ComputeDerivatives_CVODE, 0.0, m_cvy);

        CVodeSetMaxNumSteps(m_cvodeSolver, m_maxSteps);
        CVodeSetMaxOrd(m_cvodeSolver, std::min(m_order, 12));
        CVodeSStolerances(m_cvodeSolver, m_relTol, m_absTol);

      }
      break;
    case CVODE_BDF:
      {

        m_cvodeSolver = CVodeCreate(CV_BDF, CV_FUNCTIONAL);
        m_solver = &ODESolver::solveCVODE;

#ifdef USE_CVODE_OPENMP
        m_cvy = N_VNew_OpenMP(m_size, omp_get_max_threads());
#else
        m_cvy = N_VNew_Serial(m_size);
#endif

        CVodeInit(m_cvodeSolver, &ODESolver::ComputeDerivatives_CVODE, 0.0, m_cvy);

        CVodeSetMaxNumSteps(m_cvodeSolver, m_maxSteps);
        CVodeSetMaxOrd(m_cvodeSolver, std::min(m_order, 5));
        CVodeSStolerances(m_cvodeSolver, m_relTol, m_absTol);


      }
      break;
#endif
    default:
      {
        m_solver = &ODESolver::rk4;
      }
      break;
  }
}

int ODESolver::size() const
{
  return m_size;
}

void ODESolver::setSize(int size)
{
  m_size = size;
}

ODESolver::SolverType ODESolver::solverType() const
{
  return m_solverType;
}

void ODESolver::setSolverType(SolverType solverType)
{
  m_solverType = solverType;
}

int ODESolver::maxIterations() const
{
  return m_maxSteps;
}

void ODESolver::setMaxIterations(int iterations)
{
  if(iterations > 0)
    m_maxSteps = iterations;
}

int ODESolver::getIterations() const
{
  return m_currentIterations;
}

int ODESolver::order() const
{
  return m_order;
}

void ODESolver::setOrder(int order)
{
  if(order > 0)
    m_order = order;
}

double ODESolver::relativeTolerance() const
{
  return m_relTol;
}

void ODESolver::setRelativeTolerance(double tolerance)
{
  m_relTol = tolerance;
}

double ODESolver::absoluteTolerance() const
{
  return m_absTol;
}

void ODESolver::setAbsoluteTolerance(double tolerance)
{
  m_absTol = tolerance;
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

  std::fill(m_yscal, m_yscal + n, 0.0);
  std::fill(m_yerr, m_yerr + n, 0.0);
  std::fill(m_ytemp, m_ytemp + n, 0.0);
  std::fill(m_ak, m_ak + 5 * n, 0.0);
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
    m_currentIterations = nstp;

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

    errmax /= m_relTol;

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

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++)
    m_ytemp[i] = yout[i] + b21 * dt * dydt[i];

  derivs(t + a2 * dt, m_ytemp, ak2, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++)
    m_ytemp[i] = yout[i] + dt * (b31*dydt[i]+b32*ak2[i]);

  derivs(t + a3 * dt, m_ytemp, ak3, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i<n; i++)
    m_ytemp[i] = yout[i] + dt *(b41*dydt[i]+b42*ak2[i] + b43*ak3[i]);

  derivs(t + a4 * dt, m_ytemp, ak4, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i<n; i++)
    m_ytemp[i] = yout[i] + dt *(b51*dydt[i]+b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);

  derivs(t + a5 * dt, m_ytemp, ak5, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i<n; i++)
    m_ytemp[i] = yout[i] + dt * (b61 * dydt[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i]
                                 + b65 * ak5[i]);
  derivs(t + a6 * dt, m_ytemp, ak6, userData);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i<n; i++)
    m_ytemp[i] = yout[i] + dt *(c1 * dydt[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i<n; i++)
    m_yerr[i] = dt *(dc1 * dydt[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
}

#ifdef USE_CVODE

int ODESolver::solveCVODE(double y[], int n, double t, double dt, double yout[], ComputeDerivatives derivs, void *userData)
{
  RedirectionData redirectData; redirectData.deriv = derivs; redirectData.userData = userData;
  CVodeSetUserData(m_cvodeSolver, &redirectData);

#ifdef USE_CVODE_OPENMP
  double *y0 =  N_VGetArrayPointer_OpenMP(m_cvy);
#else
  double *y0 =  N_VGetArrayPointer_Serial(m_cvy);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < n; i++)
  {
    y0[i] = y[i];
  }

  CVodeReInit(m_cvodeSolver, t, m_cvy);

#ifdef USE_CVODE_OPENMP
  N_Vector ycvout = N_VMake_OpenMP(n,yout,omp_get_max_threads());
#else
  N_Vector ycvout = N_VMake_Serial(n,yout);
#endif

  double tNext = t+dt;
  double tOut = 0.0;

  int result = CVode(m_cvodeSolver, tNext, ycvout, &tOut, CV_NORMAL);

  long currentIterations = 0;
  CVodeGetNumSteps(m_cvodeSolver, &currentIterations);

  m_currentIterations = currentIterations;

  if(m_print)
    N_VPrint_Serial(ycvout);

#ifdef USE_CVODE_OPENMP
  N_VDestroy_OpenMP(ycvout);
#else
  N_VDestroy_Serial(ycvout);
#endif


  return result;
}

int ODESolver::ComputeDerivatives_CVODE(realtype t, N_Vector y, N_Vector dydt, void *user_data)
{
  RedirectionData *redirectDada = (RedirectionData*) user_data;

#ifdef USE_CVODE_OPENMP
  double *yData =  N_VGetArrayPointer_OpenMP(y);
  double *dydtData =  N_VGetArrayPointer_OpenMP(dydt);
#else
  double *yData = N_VGetArrayPointer(y);
  double *dydtData =  N_VGetArrayPointer(dydt);
#endif

  redirectDada->deriv(t, yData, dydtData, redirectDada->userData);

  return 0;
}

#endif

void ODESolver::clearMemory()
{
#ifdef USE_CVODE
  if(m_cvodeSolver)
  {
#ifdef USE_CVODE_OPENMP
    N_VDestroy_OpenMP(m_cvy);
    m_cvy = nullptr;
#else
    N_VDestroy_Serial(m_cvy);
    m_cvy = nullptr;
#endif

    CVodeFree(&m_cvodeSolver);
    m_cvodeSolver = nullptr;
  }
#endif

  if(m_yscal)
  {
    delete[] m_yscal; m_yscal = nullptr;
    delete[] m_yerr; m_yerr = nullptr;
    delete[] m_ytemp; m_ytemp = nullptr;
    delete[] m_ak; m_ak = nullptr;
  }
}
