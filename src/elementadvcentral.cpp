#include "stdafx.h"
#include "elementadvupwind.h"
#include "elementadvcentral.h"
#include "element.h"
#include "elementjunction.h"

void ElementAdvCentral::setAdvectionFunction(Element *element)
{
  if(element->flow.value >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvCentral::fluxUpNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunctionBC;
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
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunctionBC;
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
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunctionBC;
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
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunctionBC;
      }
    }
  }
}

double ElementAdvCentral::fluxUpNeighbour(Element *element, double dt, double T[])
{
  double denom = (1.0 / element->upstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->rho_cp  * (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] * upstreamFactor +
                        element->flow.value  * T[element->tIndex] * centerFactor);

  return incomingFlux;
}

double ElementAdvCentral::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double denom = (1.0 / element->downstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0 / denom;

  double outgoingFlux = -element->rho_cp * (element->downstreamElement->flow.value * T[element->downstreamElement->tIndex] * downstreamFactor +
                 element->flow.value * T[element->tIndex] * centerFactor);

  return outgoingFlux;
}

double ElementAdvCentral::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double denom = (1.0 / element->upstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0 / denom;

  double incomingFlux = element->upstreamElement->flow.value * S[element->upstreamElement->tIndex] * upstreamFactor +
                        element->flow.value  * S[element->sIndex[soluteIndex]] * centerFactor;

  return incomingFlux;
}

double ElementAdvCentral::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double denom = (1.0 / element->downstreamElement->length / 2.0) + (1.0 / element->length/ 2.0);
  double centerFactor = 1.0 / element->length / 2.0 / denom;
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0 / denom;

  double outgoingFlux = -(element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]] * downstreamFactor +
                 element->flow.value * S[element->sIndex[soluteIndex]] * centerFactor);

  return outgoingFlux;
}


