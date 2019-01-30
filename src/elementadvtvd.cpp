#include "stdafx.h"
#include "elementadvtvd.h"
#include "element.h"
#include "elementjunction.h"
#include "elementadvupwind.h"
#include "cshmodel.h"

#include "math.h"

using namespace std;

void ElementAdvTVD::setAdvectionFunction(Element *element)
{
  if(element->flow >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      if(element->upstreamElement->upstreamElement &&!element->upstreamElement->upstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpJunction;
      }

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;
      }
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunction;

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunction;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;
      }
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        if(element->upstreamElement->upstreamElement &&!element->upstreamElement->upstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvTVD::fluxUpNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvTVD::fluxUpJunction;
        }

        if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
        }
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunction;

        if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunction;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
        }
      }
    }
  }
  else
  {
    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      if(element->downstreamElement->downstreamElement &&!element->downstreamElement->downstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownJunction;
      }

      if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;
      }
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxUpNeighbourDownstreamJunction;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunction;
      }
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
      {
        if(element->downstreamElement->downstreamElement &&!element->downstreamElement->downstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownJunction;
        }

        if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;
        }
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;

        if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxUpNeighbourDownstreamJunction;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunction;
        }
      }

    }
  }
}

double ElementAdvTVD::fluxUpNeighbour(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * T[element->upstreamElement->index] -
                        element->upstreamElement->upstreamElement->flow * T[element->upstreamElement->upstreamElement->index]) /
                        (element->upstreamElement->length + element->upstreamElement->upstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow * T[element->upstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * T[element->upstreamElement->index] -
                        element->upstreamElement->flow * element->upstreamElement->upstreamJunction->temperature.value) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow * T[element->upstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * T[element->index] - element->flow * element->upstreamJunction->temperature.value) /
                        element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]);

  return -outgoingFlux;
}


double ElementAdvTVD::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow * T[element->downstreamElement->index] -
                    element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow * T[element->downstreamElement->downstreamElement->index] -
                         element->downstreamElement->flow * T[element->downstreamElement->index]) /
                        (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow * T[element->downstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->downstreamElement->flow * T[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow * T[element->downstreamElement->index] -
                    element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow * element->downstreamElement->downstreamJunction->temperature.value -
                         element->downstreamElement->flow * T[element->downstreamElement->index]) /
                        (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow * T[element->downstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->downstreamElement->flow * T[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow * T[element->downstreamElement->index] -
                         element->flow * T[element->index]) /
                        (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow * T[element->upstreamElement->index] - element->flow * T[element->index]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow * element->downstreamJunction->temperature.value -
                         element->flow * T[element->index]) /
                         element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow * T[element->upstreamElement->index] - element->flow * T[element->index]);

  return outgoingFlux;
}


double ElementAdvTVD::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * S[element->upstreamElement->index] -
                        element->upstreamElement->upstreamElement->flow * S[element->upstreamElement->upstreamElement->index]) /
                        (element->upstreamElement->length + element->upstreamElement->upstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow * S[element->upstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * S[element->upstreamElement->index] -
                        element->upstreamElement->flow * element->upstreamElement->upstreamJunction->soluteConcs[soluteIndex].value) /
                       (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow * S[element->upstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * S[element->index] - element->flow * element->upstreamJunction->soluteConcs[soluteIndex].value) /
                        element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]);

  return -outgoingFlux;
}


double ElementAdvTVD::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow * S[element->downstreamElement->index] -
                    element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow * S[element->downstreamElement->downstreamElement->index] -
                         element->downstreamElement->flow * S[element->downstreamElement->index]) /
                        (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow * S[element->downstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->downstreamElement->flow * S[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow * S[element->downstreamElement->index] -
                    element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow * element->downstreamElement->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->downstreamElement->flow * S[element->downstreamElement->index]) /
                        (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow * S[element->downstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->downstreamElement->flow * S[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow * S[element->downstreamElement->index] -
                         element->flow * S[element->index]) /
                        (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow * S[element->upstreamElement->index] - element->flow * S[element->index]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow * element->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->flow * S[element->index]) /
                         element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow * S[element->upstreamElement->index] - element->flow * S[element->index]);

  return outgoingFlux;
}


double ElementAdvTVD::computeTVDLimiter(double r, int limiter)
{
  double r_func = 0;

  switch (limiter)
  {
    //Min-Mod
    case 0:
      {
        r_func = max(0.0 , min(1.0,r));
      }
      break;
      //Superbee
    case 1:
      {
        r_func = max({0.0, min(2.0 * r,1.0), min(r, 2.0)});
      }
      break;
      //Van Leer
    case 2:
      {
        r_func = (r + fabs(r))/(1.0 + fabs(r));
      }
      break;
      //MUSCL
    case 3:
      {
        r_func = max(0.0, min({2.0 * r, (r + 1.0)/2.0, 2.0}));
      }
      break;
      //Sweby
    case 4:
      {
        r_func = max({0.0, min(1.5 * r, 1.0), min(r, 1.5)});
      }
      break;
      //Van Albada
    case 5:
      {
        r_func = (r + r*r)/(1.0 + r*r);
      }
      break;
      //QUICK
    case 6:
      {
        r_func = max(0.0, min({2.0 * r, (3.0 + r)/4, 2.0}));
      }
      break;
      //UMIST
    case 7:
      {
        r_func = max(0.0, min({2.0 * r, (1.0 + 3.0 * r)/4, (3 + r)/4.0, 2.0}));
      }
      break;
  }

  return r_func;
}
