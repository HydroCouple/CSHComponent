#include "stdafx.h"
#include "elementadvhybrid.h"
#include "elementadvcentral.h"
#include "elementadvupwind.h"
#include "element.h"
#include "elementjunction.h"

#include <math.h>

using namespace std;

void ElementAdvHybrid::setAdvectionFunction(Element *element)
{
  if(element->flow.value >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      if(fabs(element->upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
      {
        element->computeTempAdvDeriv[0] = &ElementAdvCentral::fluxUpNeighbour;
      }
      else if(element->upstreamPecletNumber < 2.0)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvHybrid::fluxUpNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpNeighbour;
      }
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunctionBC;
    }

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      if(fabs(element->downstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
      {
        element->computeTempAdvDeriv[1] = &ElementAdvCentral::fluxDownNeighbour;
      }
      else if(element->downstreamPecletNumber < 2.0)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvHybrid::fluxDownNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;
      }
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        if(fabs(element->upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvCentral::fluxUpNeighbour;
        }
        else if(element->upstreamPecletNumber < 2.0)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvHybrid::fluxUpNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpNeighbour;
        }
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunctionBC;
      }

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
      {
        if(fabs(element->downstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvCentral::fluxDownNeighbour;
        }
        else if(element->downstreamPecletNumber < 2.0)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvHybrid::fluxDownNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
        }
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
      if(fabs(element->upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
      {
        element->computeTempAdvDeriv[0] = &ElementAdvCentral::fluxUpNeighbour;
      }
      else if(element->upstreamPecletNumber > -2.0)
      {
        element->computeTempAdvDeriv[0] = &ElementAdvHybrid::fluxUpNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;
      }
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;
    }

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      if(fabs(element->downstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
      {
        element->computeTempAdvDeriv[1] = &ElementAdvCentral::fluxDownNeighbour;
      }
      else if(element->downstreamPecletNumber > -2.0)
      {
        element->computeTempAdvDeriv[1] = &ElementAdvHybrid::fluxDownNeighbour;
      }
      else
      {
        element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownNeighbor;
      }
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunctionBC;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        if(fabs(element->upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvCentral::fluxUpNeighbour;
        }
        else if(element->upstreamPecletNumber > -2.0)
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvHybrid::fluxUpNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;
        }
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;
      }

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
      {
        if(fabs(element->downstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvCentral::fluxDownNeighbour;
        }
        else if(element->downstreamPecletNumber > -2.0)
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvHybrid::fluxDownNeighbour;
        }
        else
        {
          element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownNeighbor;
        }
      }
      else
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunctionBC;
      }
    }
  }
}

double ElementAdvHybrid::fluxUpNeighbour(Element *element, double dt, double T[])
{
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0;
  double centerFactor = 1.0 / element->length / 2.0;

  double idwDenomFactor = upstreamFactor + centerFactor;

  upstreamFactor = upstreamFactor / idwDenomFactor;
  centerFactor   = centerFactor / idwDenomFactor;

  upstreamFactor = (1 + (1.0 / element->upstreamPecletNumber / upstreamFactor)) * upstreamFactor;
  centerFactor   = (1 - (1.0 / element->upstreamPecletNumber / centerFactor)) * centerFactor;

  double incomingFlux = element->rho_cp * (element->upstreamElement->flow.value * T[element->upstreamElement->tIndex] * upstreamFactor +
                        element->flow.value * T[element->tIndex] * centerFactor);

  return incomingFlux;
}

double ElementAdvHybrid::fluxDownNeighbour(Element *element, double dt, double T[])
{
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0;
  double centerFactor = 1.0 / element->length / 2.0;

  double idwDenomFactor = centerFactor + downstreamFactor;

  centerFactor /= idwDenomFactor;
  downstreamFactor   /= idwDenomFactor;

  centerFactor = (1 + (1.0 / element->downstreamPecletNumber/ centerFactor )) * centerFactor;
  downstreamFactor = (1 - (1.0 / element->downstreamPecletNumber / downstreamFactor)) * downstreamFactor;

  double outgoingFlux = element->rho_cp * (element->flow.value * T[element->tIndex] * centerFactor +
                          element-> downstreamElement->flow.value * T[element->downstreamElement->tIndex] * downstreamFactor);

  return -outgoingFlux;
}

double ElementAdvHybrid::fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double upstreamFactor = 1.0 / element->upstreamElement->length / 2.0;
  double centerFactor = 1.0 / element->length / 2.0;

  double idwDenomFactor = upstreamFactor + centerFactor;

  upstreamFactor = upstreamFactor / idwDenomFactor;
  centerFactor   = centerFactor / idwDenomFactor;

  upstreamFactor = (1 + (1.0 / element->upstreamPecletNumber / upstreamFactor)) * upstreamFactor;
  centerFactor   = (1 - (1.0 / element->upstreamPecletNumber / centerFactor)) * centerFactor;

  double incomingFlux = element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]] * upstreamFactor +
                        element->flow.value * S[element->sIndex[soluteIndex]] * centerFactor;

  return incomingFlux;
}

double ElementAdvHybrid::fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  double downstreamFactor = 1.0 / element->downstreamElement->length / 2.0;
  double centerFactor = 1.0 / element->length / 2.0;

  double idwDenomFactor = centerFactor + downstreamFactor;

  centerFactor /= idwDenomFactor;
  downstreamFactor   /= idwDenomFactor;

  centerFactor = (1 + (1.0 / element->downstreamPecletNumber/ centerFactor )) * centerFactor;
  downstreamFactor = (1 - (1.0 / element->downstreamPecletNumber / downstreamFactor)) * downstreamFactor;

  double outgoingFlux =   element->flow.value * S[element->sIndex[soluteIndex]] * centerFactor +
                          element-> downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]] * downstreamFactor;

  return -outgoingFlux;
}


