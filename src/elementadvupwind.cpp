#include "stdafx.h"
#include "elementadvupwind.h"
#include "element.h"
#include "cshmodel.h"
#include "elementjunction.h"

void ElementAdvUpwind::setAdvectionFunction(Element *element)
{
  if(element->flow >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpNeighbour;
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunction;
    }

    element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC)
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpNeighbour;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunction;
      }

      element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
    }
  }
  else
  {
    element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC)
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownNeighbor;
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunction;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC)
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownNeighbor;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunction;
      }
    }
  }
}

double ElementAdvUpwind::inFluxUpJunction(Element *element, double dt, double T[])
{
  return element->rho_cp * element->flow *  element->upstreamJunction->temperature.value;
}

double ElementAdvUpwind::inFluxUpNeighbour(Element *element, double dt, double T[])
{
  return element->rho_cp * element->upstreamElement->flow * T[element->upstreamElement->index];
}

double ElementAdvUpwind::inFluxSelf(Element *element, double dt, double T[])
{
  return element->rho_cp * element->flow * T[element->index];
}

double ElementAdvUpwind::outFluxSelf(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->flow * T[element->index];
}

double ElementAdvUpwind::outFluxDownJunction(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->flow * element->downstreamJunction->temperature.value;
}

double ElementAdvUpwind::outFluxDownNeighbor(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->downstreamElement->flow * T[element->downstreamElement->index];
}

double ElementAdvUpwind::inFluxUpJunction(Element *element, double dt, double S[], int soluteIndex)
{
  return element->flow * element->upstreamJunction->soluteConcs[soluteIndex].value;
}

double ElementAdvUpwind::inFluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  return element->upstreamElement->flow * S[element->upstreamElement->index];
}

double ElementAdvUpwind::inFluxSelf(Element *element, double dt, double S[], int soluteIndex)
{
  return element->flow * S[element->index];
}

double ElementAdvUpwind::outFluxSelf(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->flow * S[element->index];
}

double ElementAdvUpwind::outFluxDownJunction(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->flow * element->downstreamJunction->soluteConcs[soluteIndex].value;
}

double ElementAdvUpwind::outFluxDownNeighbor(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->downstreamElement->flow * S[element->downstreamElement->index];
}




