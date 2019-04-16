#include "stdafx.h"
#include "elementadvupwind.h"
#include "element.h"
#include "cshmodel.h"
#include "elementjunction.h"

void ElementAdvUpwind::setAdvectionFunction(Element *element)
{
  if(element->flow.value >= 0)
  {
    if(element->upstreamElement && !element->upstreamJunction->temperature.isBC && !element->upstreamJunction->tIndex >= 0)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpNeighbour;
    }
    else if(element->upstreamJunction->tIndex > -1)
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunction;
    }
    else
    {
      element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxUpJunctionBC;
    }

    element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxSelf;

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      if(element->upstreamElement && !element->upstreamJunction->soluteConcs[i].isBC && !element->upstreamJunction->sIndex[i] >= 0)
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpNeighbour;
      }
      else if(element->upstreamJunction->sIndex[i] > -1)
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunction;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxUpJunctionBC;
      }

      element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxSelf;
    }
  }
  else
  {
    element->computeTempAdvDeriv[0] = &ElementAdvUpwind::inFluxSelf;

    if(element->downstreamElement && !element->downstreamJunction->temperature.isBC && !element->downstreamJunction->tIndex >= 0)
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownNeighbor;
    }
    else if(element->downstreamJunction->tIndex > -1)
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunction;
    }
    else
    {
      element->computeTempAdvDeriv[1] = &ElementAdvUpwind::outFluxDownJunctionBC;
    }

    for(int i = 0 ; i < element->numSolutes; i++)
    {
      element->computeSoluteAdvDeriv[i][0] = &ElementAdvUpwind::inFluxSelf;

      if(element->downstreamElement && !element->downstreamJunction->soluteConcs[i].isBC && !element->downstreamJunction->sIndex[i] >= 0)
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownNeighbor;
      }
      else if(element->downstreamJunction->sIndex[i] > -1)
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunction;
      }
      else
      {
        element->computeSoluteAdvDeriv[i][1] = &ElementAdvUpwind::outFluxDownJunctionBC;
      }
    }
  }
}

double ElementAdvUpwind::inFluxUpJunction(Element *element, double dt, double T[])
{
  return element->rho_cp * element->flow.value *  T[element->upstreamJunction->tIndex];
}

double ElementAdvUpwind::inFluxUpJunctionBC(Element *element, double dt, double T[])
{
  return element->rho_cp * element->flow.value *  element->upstreamJunction->temperature.value;
}

double ElementAdvUpwind::inFluxUpNeighbour(Element *element, double dt, double T[])
{
  return element->rho_cp * element->upstreamElement->flow.value * T[element->upstreamElement->tIndex];
}

double ElementAdvUpwind::inFluxSelf(Element *element, double dt, double T[])
{
  return element->rho_cp * element->flow.value * T[element->tIndex];
}

double ElementAdvUpwind::outFluxSelf(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->flow.value * T[element->tIndex];
}

double ElementAdvUpwind::outFluxDownJunction(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->flow.value * T[element->downstreamJunction->tIndex];
}

double ElementAdvUpwind::outFluxDownJunctionBC(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->flow.value * element->downstreamJunction->temperature.value;
}

double ElementAdvUpwind::outFluxDownNeighbor(Element *element, double dt, double T[])
{
  return -element->rho_cp * element->downstreamElement->flow.value * T[element->downstreamElement->tIndex];
}

double ElementAdvUpwind::inFluxUpJunction(Element *element, double dt, double S[], int soluteIndex)
{
  return element->flow.value * S[element->upstreamJunction->sIndex[soluteIndex]];
}

double ElementAdvUpwind::inFluxUpJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  return element->flow.value * element->upstreamJunction->soluteConcs[soluteIndex].value;
}

double ElementAdvUpwind::inFluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex)
{
  return element->upstreamElement->flow.value * S[element->upstreamElement->sIndex[soluteIndex]];
}

double ElementAdvUpwind::inFluxSelf(Element *element, double dt, double S[], int soluteIndex)
{
  return element->flow.value * S[element->sIndex[soluteIndex]];
}

double ElementAdvUpwind::outFluxSelf(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->flow.value * S[element->sIndex[soluteIndex]];
}

double ElementAdvUpwind::outFluxDownJunction(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->flow.value * S[element->downstreamJunction->sIndex[soluteIndex]];
}

double ElementAdvUpwind::outFluxDownJunctionBC(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->flow.value * element->downstreamJunction->soluteConcs[soluteIndex].value;
}

double ElementAdvUpwind::outFluxDownNeighbor(Element *element, double dt, double S[], int soluteIndex)
{
  return -element->downstreamElement->flow.value * S[element->downstreamElement->sIndex[soluteIndex]];
}




