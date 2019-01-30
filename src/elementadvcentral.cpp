#include "stdafx.h"
#include "elementadvupwind.h"
#include "elementadvcentral.h"
#include "element.h"
#include "elementjunction.h"

void ElementAdvCentral::setAdvectionFunction(Element *element)
{
  if(element->flow >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvCentral::fluxUpNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunction;
    }

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[1] = &ElementAdvCentral::fluxDownNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvCentral::fluxUpNeighbour;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunction;
      }

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvCentral::fluxDownNeighbour;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
      }
    }
  }
  else
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvCentral::fluxUpNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;
    }

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[1] = &ElementAdvCentral::fluxDownNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunction;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvCentral::fluxUpNeighbour;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;
      }

      if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvCentral::fluxDownNeighbour;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunction;
      }
    }
  }
}

double ElementAdvCentral::fluxUpNeighbour(Element *element, double dt, double T[])
{
  double denom = (1.0 / element->upstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp  * (element->upstreamElement->flow * T[element->upstreamElement->index] * upstreamFactor +
                        element->flow  * T[element->index] * centerFactor);

  return incomingFlux;
}

double ElementAdvCentral::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double denom = (1.0 / element->downstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0 / denom;

  double outgoingFlux = -element->rho_cp * (element->downstreamElement->flow * T[element->downstreamElement->index] * downstreamFactor +
                 element->flow * T[element->index] * centerFactor);

  return outgoingFlux;
}

double ElementAdvCentral::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double denom = (1.0 / element->upstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow * S[element->upstreamElement->index] * upstreamFactor +
                        element->flow  * S[element->index] * centerFactor;

  return incomingFlux;
}

double ElementAdvCentral::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double denom = (1.0 / element->downstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0 / denom;

  double outgoingFlux = -(element->downstreamElement->flow * S[element->downstreamElement->index] * downstreamFactor +
                 element->flow * S[element->index] * centerFactor);

  return outgoingFlux;
}


