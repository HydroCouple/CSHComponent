#include "stdafx.h"
#include "elementadvquick.h"
#include "element.h"
#include "elementjunction.h"
#include "elementadvupwind.h"


void ElementAdvQUICK::setAdvectionFunctions(Element *element)
{
  if(element->flow >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      if(element->upstreamElement->upstreamElement &&!element->upstreamElement->upstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvQUICK::fluxUpNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvQUICK::fluxUpJunction;
      }

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvQUICK::fluxDownNeighbourUpstreamNeighbour;
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
        element->computeTempAdvDeriv[1] = &ElementAdvQUICK::fluxDownNeighbourUpstreamJunction;
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
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvQUICK::fluxUpNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvQUICK::fluxUpJunction;
        }

        if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvQUICK::fluxDownNeighbourUpstreamNeighbour;
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
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvQUICK::fluxDownNeighbourUpstreamJunction;
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
        element->computeTempAdvDeriv[1] = &ElementAdvQUICK::fluxDownNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvQUICK::fluxDownJunction;
      }

      if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvQUICK::fluxUpNeighbourDownstreamNeighbour;
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
        element->computeTempAdvDeriv[1] = &ElementAdvQUICK::fluxUpNeighbourDownstreamJunction;
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
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvQUICK::fluxDownNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvQUICK::fluxDownJunction;
        }

        if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvQUICK::fluxUpNeighbourDownstreamNeighbour;
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
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvQUICK::fluxUpNeighbourDownstreamJunction;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunction;
        }
      }

    }
  }
}


double ElementAdvQUICK::fluxUpNeighbour(Element *element, double dt, double T[])
{


//  double fluxUpUp =
//  double fluxUp = (6.0/8.0) * element->upstreamElement->flow * T[element->upstreamElement->index];

//  double flux2

//  double incomingFlux = element->rho_cp * element->upstreamElement->flow * T[element->upstreamElement->index] +
//                        element->rho_cp * rw_func * interpFactor *
//                        (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]);

  return 0;
}

double ElementAdvQUICK::fluxUpJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * T[element->upstreamElement->index] -
                        element->upstreamElement->flow * element->upstreamElement->upstreamJunction->temperature.value) /
      (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->upstreamElement->flow * T[element->upstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvQUICK::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]);

  return -outgoingFlux;
}

double ElementAdvQUICK::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * T[element->index] - element->flow * element->upstreamJunction->temperature.value) /
                        element->length / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->downstreamElement->flow * T[element->downstreamElement->index] - element->flow * T[element->index]);

  return -outgoingFlux;
}



double ElementAdvQUICK::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow * T[element->downstreamElement->index] -
                    element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow * T[element->downstreamElement->downstreamElement->index] -
                         element->downstreamElement->flow * T[element->downstreamElement->index]) /
                        (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow * T[element->downstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->downstreamElement->flow * T[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvQUICK::fluxDownJunction(Element *element, double dt, double T[])
{
  double dwdenom = (element->downstreamElement->flow * T[element->downstreamElement->index] -
                    element->flow * T[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow * element->downstreamElement->downstreamJunction->temperature.value -
                         element->downstreamElement->flow * T[element->downstreamElement->index]) /
                        (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp * element->downstreamElement->flow * T[element->downstreamElement->index] +
                        element->rho_cp * rw_func * interpFactor *
                        (element->flow * T[element->index] - element->downstreamElement->flow * T[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvQUICK::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double T[])
{
  double dedenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow * T[element->downstreamElement->index] -
                         element->flow * T[element->index]) /
                        (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow * T[element->upstreamElement->index] - element->flow * T[element->index]);

  return outgoingFlux;
}

double ElementAdvQUICK::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double T[])
{
  double dedenom = (element->flow * T[element->index] - element->upstreamElement->flow * T[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow * element->downstreamJunction->temperature.value -
                         element->flow * T[element->index]) /
                         element->length / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->rho_cp * element->flow * T[element->index] +
                        element->rho_cp * re_func * interpFactor *
                        (element->upstreamElement->flow * T[element->upstreamElement->index] - element->flow * T[element->index]);

  return outgoingFlux;
}



double ElementAdvQUICK::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{

  double incomingFlux = 0.0;

  return incomingFlux;
}

double ElementAdvQUICK::fluxUpJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0;

  double rw = dwdenom ? (element->upstreamElement->flow * S[element->upstreamElement->index] -
                        element->upstreamElement->flow * element->upstreamElement->upstreamJunction->soluteConcs[soluteIndex].value) /
                       (element->upstreamElement->length / 2.0) / dwdenom : 0;


  double rw_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow * S[element->upstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]);

  return incomingFlux;
}

double ElementAdvQUICK::fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]) /
      (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
      (element->length + element->upstreamElement->length) / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]);

  return -outgoingFlux;
}

double ElementAdvQUICK::fluxDownNeighbourUpstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double re = dedenom ? (element->flow * S[element->index] - element->flow * element->upstreamJunction->soluteConcs[soluteIndex].value) /
                        element->length / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->downstreamElement->flow * S[element->downstreamElement->index] - element->flow * S[element->index]);

  return -outgoingFlux;
}



double ElementAdvQUICK::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow * S[element->downstreamElement->index] -
                    element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->downstreamElement->flow * S[element->downstreamElement->downstreamElement->index] -
                         element->downstreamElement->flow * S[element->downstreamElement->index]) /
                        (element->downstreamElement->downstreamElement->length + element->downstreamElement->length) / 2.0 / dwdenom : 0;

  double rw_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow * S[element->downstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->downstreamElement->flow * S[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvQUICK::fluxDownJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dwdenom = (element->downstreamElement->flow * S[element->downstreamElement->index] -
                    element->flow * S[element->index]) /
                   (element->downstreamElement->length + element->length) / 2.0;

  double rw = dwdenom ? (element->downstreamElement->flow * element->downstreamElement->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->downstreamElement->flow * S[element->downstreamElement->index]) /
                        (element->downstreamElement->length / 2.0) / dwdenom : 0;

  double rw_func = 0.0;

  double denom = (element->downstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->downstreamElement->length / 2.0 / denom;

  double incomingFlux = element->downstreamElement->flow * S[element->downstreamElement->index] +
                        rw_func * interpFactor *
                        (element->flow * S[element->index] - element->downstreamElement->flow * S[element->downstreamElement->index]);

  return -incomingFlux;
}

double ElementAdvQUICK::fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->downstreamElement->flow * S[element->downstreamElement->index] -
                         element->flow * S[element->index]) /
                        (element->downstreamElement->length + element->length) / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow * S[element->upstreamElement->index] - element->flow * S[element->index]);

  return outgoingFlux;
}

double ElementAdvQUICK::fluxUpNeighbourDownstreamJunction(Element *element, double dt, double S[], int soluteIndex)
{
  double dedenom = (element->flow * S[element->index] - element->upstreamElement->flow * S[element->upstreamElement->index]) /
                   (element->length + element->upstreamElement->length) / 2.0;

  double re = dedenom ? (element->flow * element->downstreamJunction->soluteConcs[soluteIndex].value -
                         element->flow * S[element->index]) /
                         element->length / 2.0 / dedenom : 0.0;

  double re_func = 0.0;

  double denom = (element->upstreamElement->length / 2.0) + (element->length / 2.0);
  double interpFactor = element->length / 2.0 / denom;

  double outgoingFlux = element->flow * S[element->index] +
                        re_func * interpFactor *
                        (element->upstreamElement->flow * S[element->upstreamElement->index] - element->flow * S[element->index]);

  return outgoingFlux;
}

