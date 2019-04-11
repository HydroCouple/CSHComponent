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
  if(element->flow.value >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      if(element->upstreamElement->upstreamElement &&!element->upstreamElement->upstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpNeighbour;
      }
      else if(element->upstreamElement->upstreamJunction->tIndex > -1)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpJunction;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvTVD::fluxUpJunctionBC;
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
      if(element->upstreamJunction->tIndex > -1)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunction;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunctionBC;
      }

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        if(element->upstreamJunction->tIndex > -1)
        {
          element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunction;
        }
        else
        {
          element->computeTempAdvDeriv[1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunctionBC;
        }
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
        else if(element->upstreamElement->upstreamJunction->sIndex[i] > -1)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvTVD::fluxUpJunction;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvTVD::fluxUpJunctionBC;
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

        if(element->upstreamJunction->sIndex[i] > -1)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunction;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunctionBC;
        }

        if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
        {
          if(element->upstreamJunction->sIndex[i] > -1)
          {
            element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunction;
          }
          else
          {
            element->computeSoluteAdvDeriv[i][1] = &ElementAdvTVD::fluxDownNeighbourUpstreamJunctionBC;
          }
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
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunctionBC;
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
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunctionBC;
        }
      }

    }
  }
}

double ElementAdvTVD::fluxUpNeighbour(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] -
                        element->upstreamElement->upstreamElement->flow.value * T[element->upstreamElement->upstreamElement->tIndex]) /
      (element->upstreamElement->length + element->upstreamElement->upstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] -
                        element->upstreamElement->flow.value * T[element->upstreamElement->upstreamJunction->tIndex]) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunctionBC(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] -
                        element->upstreamElement->flow.value * element->upstreamElement->upstreamJunction->temperature.value) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]);

  return incomingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * T[element->tIndex] - element->flow.value * T[element->upstreamJunction->tIndex]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunctionBC(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * T[element->tIndex] - element->flow.value * element->upstreamJunction->temperature.value) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return -outgoingFlux;
}


double ElementAdvTVD::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] -
                   element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow.value * T[element->downstreamElement->downstreamElement->tIndex] -
                        element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]) /
      (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] -
                   element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow.value * element->downstreamElement->downstreamJunction->temperature.value -
                         element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]) /
      (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunctionBC(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] -
                   element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow.value * element->downstreamElement->downstreamJunction->temperature.value -
                         element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]) /
      (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow.value * T[element->tIndex] - element->downstreamElement->flow.value * T[element->downstreamElement->tIndex]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] -
                        element->flow.value * T[element->tIndex]) /
      (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow.value * T[element->downstreamJunction->tIndex] -
                        element->flow.value * T[element->tIndex]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunctionBC(Element *element, double dt, double T[])
{
  double dedenom = (element->flow.value * T[element->tIndex] - element->upstreamElement->flow.value * T[element->upstreamElement->tIndex]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow.value * element->downstreamJunction->temperature.value -
                         element->flow.value * T[element->tIndex]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow.value * T[element->tIndex] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] - element->flow.value * T[element->tIndex]);

  return outgoingFlux;
}


double ElementAdvTVD::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] -
                        element->upstreamElement->upstreamElement->flow.value * S[element->upstreamElement->upstreamElement->sIndex[soluteIndex]]) /
      (element->upstreamElement->length + element->upstreamElement->upstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] +
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] -
                        element->upstreamElement->flow.value * S[element->upstreamElement->upstreamJunction->sIndex[soluteIndex]]) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] +
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]);

  return incomingFlux;
}

double ElementAdvTVD::fluxUpJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] -
                        element->upstreamElement->flow.value * element->upstreamElement->upstreamJunction->soluteConcs[soluteIndex].value) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] +
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]);

  return incomingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * S[element->sIndex[soluteIndex]] - element->flow.value * S[element->upstreamJunction->sIndex[soluteIndex]]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]);

  return -outgoingFlux;
}

double ElementAdvTVD::fluxDownNeighbourUpstreamJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow.value * S[element->sIndex[soluteIndex]] - element->flow.value * element->upstreamJunction->soluteConcs[soluteIndex].value) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]- element->flow.value * S[element->sIndex[soluteIndex]]);

  return -outgoingFlux;
}


double ElementAdvTVD::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]-
                   element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow.value * S[element->downstreamElement->downstreamElement->sIndex[soluteIndex]] -
                        element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]) /
      (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]+
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]-
                   element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow.value * S[element->downstreamElement->downstreamJunction->sIndex[soluteIndex]] -
                        element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]) /
      (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]+
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxDownJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]-
                   element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow.value * element->downstreamElement->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]) /
      (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = computeTVDLimiter(rw, element->model->m_TVDFluxLimiter, element, 1);

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]+
                        rw_func * interpFactor *
                        (element->flow.value * S[element->sIndex[soluteIndex]] - element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]);

  return -incomingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]]-
                        element->flow.value * S[element->sIndex[soluteIndex]]) /
      (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] - element->flow.value * S[element->sIndex[soluteIndex]]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow.value * S[element->downstreamJunction->sIndex[soluteIndex]] -
                        element->flow.value * S[element->sIndex[soluteIndex]]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] - element->flow.value * S[element->sIndex[soluteIndex]]);

  return outgoingFlux;
}

double ElementAdvTVD::fluxUpNeighbourDownstreamJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow.value * S[element->sIndex[soluteIndex]] - element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow.value * element->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->flow.value * S[element->sIndex[soluteIndex]]) /
      element->length / 2.0 / dedenom : 0.0;

  double re_func = computeTVDLimiter(re, element->model->m_TVDFluxLimiter, element);

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow.value * S[element->sIndex[soluteIndex]] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] - element->flow.value * S[element->sIndex[soluteIndex]]);

  return outgoingFlux;
}


double ElementAdvTVD::computeTVDLimiter(double r, TVDFluxLimiter limiter, Element *element, int upstream)
{
  double r_func = 0;

  switch (limiter)
  {
    //Min-Mod
    case MIN_MOD:
      {
        r_func = max(0.0 , min(1.0,r));
      }
      break;
      //Superbee
    case SUPERBEE:
      {
        r_func = max({0.0, min(2.0 * r,1.0), min(r, 2.0)});
      }
      break;
      //Van Leer
    case VAN_LEER:
      {
        r_func = (r + fabs(r))/(1.0 + fabs(r));
      }
      break;
      //MUSCL
    case MUSCL:
      {
        r_func = max(0.0, min({2.0 * r, (r + 1.0)/2.0, 2.0}));
      }
      break;
      //Sweby
    case SWEBY:
      {
        r_func = max({0.0, min(1.5 * r, 1.0), min(r, 1.5)});
      }
      break;
      //Van Albada
    case VAN_ALBADA:
      {
        r_func = (r + r*r)/(1.0 + r*r);
      }
      break;
      //QUICK
    case QUICK:
      {
        r_func = max(0.0, min({2.0 * r, (3.0 + r)/4, 2.0}));
      }
      break;
      //UMIST
    case UMIST:
      {
        r_func = max(0.0, min({2.0 * r, (1.0 + 3.0 * r)/4, (3 + r)/4.0, 2.0}));
      }
      break;
    case SOU:
      {
        r_func = r;
      }
      break;
    case FROMM:
      {
        r_func =  (1 + r) / 2.0;
      }
      break;
    case ULTIMATE_QUICKEST:
      {
        double cn = 0;

        switch (upstream)
        {
          case 0:
            {
              cn = element->upstreamCourantNumber;
            }
            break;
          default:
            {
              cn = element->downstreamCourantNumber;
            }
            break;
        }

        //        if(cn)
        {
          r_func = max(0.0,
                       min({0.5 * (1.0 * r) + (1.0 - r)*(1.0 - 2.0 * fabs(cn)) / 6.0,
                            2.0 / (1.0 - fabs(cn)),
                            2.0 * r / fabs(cn)})
                       );
        }
      }
      break;
    case SUPER_C:
      {
        double cn = 0;

        switch (upstream)
        {
          case 0:
            {
              cn = element->upstreamCourantNumber;
            }
            break;
          default:
            {
              cn = element->downstreamCourantNumber;
            }
            break;
        }

        if(r >= 0.0 && r <= 1.0)
        {
          r_func = min(1.0,  2.0 * r / fabs(cn));
        }
        else if(r > 1)
        {
          r_func = min(r , 2 / (1.0 - fabs(cn)));
        }
      }
      break;
    case HYPER_C:
      {
        double cn = 0;

        switch (upstream)
        {
          case 0:
            {
              cn = element->upstreamCourantNumber;
            }
            break;
          default:
            {
              cn = element->downstreamCourantNumber;
            }
            break;
        }

        if(r > 0)
        {
          r_func = min(2 * r / fabs(cn),  2 / (1 - fabs(cn)));
        }
      }
      break;
  }

  return r_func;
}
