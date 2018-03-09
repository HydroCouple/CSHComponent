/*!
*  \file    element.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
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
#include "element.h"
#include "elementjunction.h"
#include "stmmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  STMModel *model)
  : id(id),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    length(0.0),
    depth(0.0),
    xSectionArea(0.0),
    width(0.0),
    flow(0.0),
    slope(0.0),
    externalHeatFluxes(0.0),
    radiationFluxes(0.0),
    externalSoluteFluxes(nullptr),
    heatBalance(0.0),
    soluteMassBalance(nullptr),
    pecletNumber(0.0),
    upstreamElement(nullptr),
    downstreamElement(nullptr),
    model(model),
    computeTempAdv(nullptr),
    computeSoluteAdv(nullptr)
{

  initializeSolutes();

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

  computeTempAdv = &Element::computeDTDtUpwind;
  computeSoluteAdv = &Element::computeDSoluteDtUpwind;
}

Element::~Element()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
    delete[] soluteMassBalance;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  //set upstream and downstream elements
  setUpstreamElement();
  setDownStreamElement();

  switch (model->m_advectionMode)
  {
    case STMModel::AdvectionDiscretizationMode::Central:
      {
        computeTempAdv = &Element::computeDTDtCentral;
        computeSoluteAdv = &Element::computeDSoluteDtCentral;
      }
      break;
    case STMModel::AdvectionDiscretizationMode::Hybrid:
      {
        computeTempAdv = &Element::computeDTDtHybrid;
        computeSoluteAdv = &Element::computeDSoluteDtHybrid;
      }
      break;
    case STMModel::AdvectionDiscretizationMode::TVD:
      {
        computeTempAdv = &Element::computeDTDtTVD;
        computeSoluteAdv = &Element::computeDSoluteDtTVD;
      }
      break;
    default:
      {
        computeTempAdv = &Element::computeDTDtUpwind;
        computeSoluteAdv = &Element::computeDSoluteDtUpwind;
      }
      break;
  }
}

void Element::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] soluteMassBalance; soluteMassBalance = nullptr;
    
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    soluteMassBalance = new double[numSolutes]();
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;
  double volume = xSectionArea * length;
  double rho_cp = model->m_waterDensity * model->m_cp;
  double rho_cp_vol = rho_cp * volume;

  //Compute advection
  DTDt += (this->*computeTempAdv)(dt, T);

  //Compute dispersion
  DTDt += computeDTDtDispersion(dt, T);

  //Add external sources
  {
    DTDt += radiationFluxes / depth / rho_cp;
    DTDt += externalHeatFluxes / rho_cp_vol;
  }

  return DTDt;
}

double Element::computeDTDtDispersion(double dt, double T[])
{
  double volume = xSectionArea * length;
  double DTDt = 0.0;

  //If upstream element exists and upstream junction is not a boundary condition
  //use upstream element to compute derivative.
  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {

    DTDt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
            (T[upstreamElement->index] - T[index]) / (volume * (length + upstreamElement->length) / 2.0);

  }
  //Otherwise use upstream junction to compute temperature derivate
  else
  {
    DTDt +=  upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
             (upstreamJunction->temperature.value - T[index]) / (volume * length / 2.0);
  }

  //If downstream element exists and downstream junction element is not a boundary condition
  //use downstream element to compute derivative.
  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    DTDt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
            (T[index] - T[downstreamElement->index]) / (volume * (length + downstreamElement->length) / 2.0);
  }
  //Otherwise use downstream element junction
  else
  {
    DTDt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
            (T[index] - downstreamJunction->temperature.value) / (volume * length / 2.0);
  }

  return DTDt;
}

double Element::computeDTDtUpwind(double dt, double T[])
{
  double volume = xSectionArea * length;
  double incomingFlux = 0.0;
  double outgoingFlux = 0.0;

  //Flow goes from upstream to downstream
  if(flow >= 0)
  {
    if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
    {
      incomingFlux = upstreamFlow * T[upstreamElement->index] / volume;
      outgoingFlux = flow * T[index] / volume;
    }
    else
    {
      incomingFlux = upstreamFlow * upstreamJunction->temperature.value / volume;
      outgoingFlux = flow * T[index] / volume;
    }
  }
  //Opposite Direction
  else
  {
    if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
    {
      incomingFlux = flow * T[index] / volume;
      outgoingFlux = downstreamFlow *  T[downstreamElement->index] / volume ;
    }
    else
    {
      incomingFlux = flow * T[index] / volume;
      outgoingFlux = downstreamFlow * downstreamJunction->temperature.value / volume ;
    }
  }

  return incomingFlux - outgoingFlux;

}

double Element::computeDTDtCentral(double dt, double T[])
{

  double volume = xSectionArea * length;
  double incomingFlux = 0;
  double outgoingFlux = 0;

  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (T[upstreamElement->index] * centerFactor +
                   T[index] * upstreamFactor) / volume;
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->temperature.value / volume;
  }

  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;


    outgoingFlux = downstreamFlow * (T[downstreamElement->index] * downstreamFactor +
                   T[index] * centerFactor) / volume;
  }
  else
  {
    outgoingFlux = downstreamFlow * downstreamJunction->temperature.value / volume;
  }


  return incomingFlux - outgoingFlux;
}

double Element::computeDTDtHybrid(double dt, double T[])
{
  double volume = xSectionArea * length;
  double incomingFlux = 0;
  double outgoingFlux = 0;

  bool hasUpstream = upstreamElement != nullptr && !upstreamJunction->temperature.isBC;
  bool hasDownstream = downstreamElement != nullptr && !downstreamJunction->temperature.isBC;

  //If peclet number is zero use central differencing
  if(hasUpstream && fabs(upstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (T[upstreamElement->index] * centerFactor +
                   T[index] * upstreamFactor) / volume;
  }
  else if(hasUpstream && upstreamPecletNumber >= 2.0)
  {
    incomingFlux = upstreamFlow * T[upstreamElement->index] / volume;
  }
  else if(upstreamPecletNumber <= -2.0)
  {
    incomingFlux = flow * T[index] / volume;
  }
  else if(hasUpstream)
  {
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0;
    double centerFactor = 1.0 / length / 2.0;

    double idwDenomFactor = upstreamFactor + centerFactor;

    upstreamFactor = upstreamFactor / idwDenomFactor;
    centerFactor   = centerFactor / idwDenomFactor;

    upstreamFactor = (1 + (1.0 / upstreamPecletNumber / upstreamFactor)) * upstreamFactor;
    centerFactor   = (1 - (1.0 / upstreamPecletNumber / centerFactor)) * centerFactor;

    incomingFlux =  upstreamFlow * (T[upstreamElement->index] * upstreamFactor +
                    T[index] * centerFactor) / volume;
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->temperature.value / volume;
  }


  if(hasDownstream && fabs(downstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;


    outgoingFlux = downstreamFlow * (T[downstreamElement->index] * downstreamFactor +
                   T[index] * centerFactor) / volume;
  }
  else if(downstreamPecletNumber >= 2.0)
  {
    outgoingFlux = flow * T[index] / volume;
  }
  else if(hasDownstream && downstreamPecletNumber <= -2.0)
  {
    outgoingFlux = downstreamFlow * T[downstreamElement->index] / volume;
  }
  else if(hasDownstream)
  {
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0;
    double centerFactor = 1.0 / length / 2.0;

    double idwDenomFactor = centerFactor + downstreamFactor;

    centerFactor /= idwDenomFactor;
    downstreamFactor   /= idwDenomFactor;

    centerFactor = (1 + (1.0 / downstreamPecletNumber/ centerFactor)) * centerFactor;
    downstreamFactor = (1 - (1.0 / downstreamPecletNumber / downstreamFactor)) * downstreamFactor;

    outgoingFlux =  downstreamFlow * (T[index] * centerFactor +
                                      T[downstreamElement->index] * downstreamFactor) / volume;
  }
  else
  {
    outgoingFlux = downstreamFlow * downstreamJunction->temperature.value / volume;
  }

  double DTDt = incomingFlux - outgoingFlux;

  return DTDt;
}

double Element::computeDTDtTVD(double dt, double T[])
{
  double volume = xSectionArea * length;
  double incomingFlux = 0;
  double outgoingFlux = 0;

  bool hasUpstream = upstreamElement != nullptr && !upstreamJunction->temperature.isBC;
  bool hasDownstream = downstreamElement != nullptr && !downstreamJunction->temperature.isBC;

  double f_pos = 0;
  double f_neg = 0;

  if(hasDownstream)
  {
    f_pos = (T[downstreamElement->index] - T[index]) * downstreamElement->length / 2.0 / ( (downstreamElement->length / 2.0) + (length / 2.0));
  }
  else
  {
    f_pos = downstreamJunction->temperature.value - T[index];
  }


  if(hasUpstream)
  {
    f_neg = (T[index] - T[upstreamElement->index]) * length / 2.0 / ( (upstreamElement->length / 2.0) + (length / 2.0));
  }
  else
  {
    f_neg = T[index] - upstreamJunction->temperature.value;
  }


  //  if(flow >= 0)
  //  {
  //    outgoingFlux = downstreamFlow * (T[index] + max(f_pos, f_neg)) / volume;
  //  }
  //  else
  //  {
  //    outgoingFlux = downstreamFlow * (T[index]- max(f_pos, f_neg)) / volume;
  //  }

  if(flow >= 0)
  {
    incomingFlux = downstreamFlow * (T[index] - min(f_pos, f_neg)) / volume;
  }
  else
  {
    incomingFlux = downstreamFlow * (T[index] + min(f_pos, f_neg)) / volume;
  }


  return incomingFlux - outgoingFlux;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;
  double volume = xSectionArea * length;
  double rho_vol = model->m_waterDensity * volume;

  //Compute advection
  DSoluteDt += (this->*computeSoluteAdv)(dt, S, soluteIndex);

  //Compute dispersion
  DSoluteDt += computeDSoluteDtDispersion(dt, S, soluteIndex);

  //Add external sources
  {
    DSoluteDt += externalSoluteFluxes[soluteIndex] / rho_vol;
  }

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersion(double dt, double S[], int soluteIndex)
{
  double volume = xSectionArea * length;
  double DSoluteDt = 0.0;

  //If upstream element exists and upstream junction is not a boundary condition
  //use upstream element to compute derivative.
  if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
  {

    DSoluteDt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
                 (S[upstreamElement->index] - S[index]) / (volume * (length + upstreamElement->length) / 2.0);

  }
  //Otherwise use upstream junction to compute temperature derivate
  else
  {
    DSoluteDt +=  upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
                  (upstreamJunction->soluteConcs[soluteIndex].value - S[index]) / (volume * length / 2.0);
  }

  //If downstream element exists and downstream junction element is not a boundary condition
  //use downstream element to compute derivative.
  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    DSoluteDt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
                 (S[index] - S[downstreamElement->index]) / (volume * (length + downstreamElement->length) / 2.0);
  }
  //Otherwise use downstream element junction
  else
  {
    DSoluteDt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
                 (S[index] - downstreamJunction->soluteConcs[soluteIndex].value) / (volume * length / 2.0);
  }

  return DSoluteDt;
}

double Element::computeDSoluteDtUpwind(double dt, double S[], int soluteIndex)
{
  double volume = xSectionArea * length;
  double incomingFlux = 0.0;
  double outgoingFlux = 0.0;

  //Flow goes from upstream to downstream
  if(flow >= 0)
  {
    if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
    {
      incomingFlux = upstreamFlow * S[upstreamElement->index] / volume;
      outgoingFlux = flow * S[index] / volume;
    }
    else
    {
      incomingFlux = upstreamFlow * upstreamJunction->soluteConcs[soluteIndex].value / volume;
      outgoingFlux = flow * S[index] / volume;
    }
  }
  //Opposite Direction
  else
  {
    if(downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC)
    {
      incomingFlux = flow * S[index] / volume;
      outgoingFlux = downstreamFlow *  S[downstreamElement->index] / volume ;
    }
    else
    {
      incomingFlux = flow * S[index] / volume;
      outgoingFlux = downstreamFlow * downstreamJunction->soluteConcs[soluteIndex].value / volume ;
    }
  }

  return incomingFlux - outgoingFlux;
}

double Element::computeDSoluteDtCentral(double dt, double S[], int soluteIndex)
{
  double volume = xSectionArea * length;
  double incomingFlux = 0;
  double outgoingFlux = 0;

  if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (S[upstreamElement->index] * centerFactor +
                   S[index] * upstreamFactor) / volume;
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->soluteConcs[soluteIndex].value / volume;
  }

  if(downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;


    outgoingFlux = downstreamFlow * (S[downstreamElement->index] * downstreamFactor +
                   S[index] * centerFactor) / volume;
  }
  else
  {
    outgoingFlux = downstreamFlow * downstreamJunction->soluteConcs[soluteIndex].value / volume;
  }


  return incomingFlux - outgoingFlux;
}

double Element::computeDSoluteDtHybrid(double dt, double S[], int soluteIndex)
{
  double volume = xSectionArea * length;
  double incomingFlux = 0;
  double outgoingFlux = 0;

  bool hasUpstream = upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC;
  bool hasDownstream = downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC;

  //If peclet number is zero use central differencing
  if(hasUpstream && fabs(upstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (S[upstreamElement->index] * centerFactor +
                   S[index] * upstreamFactor) / volume;
  }
  else if(hasUpstream && upstreamPecletNumber >= 2.0)
  {
    incomingFlux = upstreamFlow * S[upstreamElement->index] / volume;
  }
  else if(upstreamPecletNumber <= -2.0)
  {
    incomingFlux = flow * S[index] / volume;
  }
  else if(hasUpstream)
  {
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0;
    double centerFactor = 1.0 / length / 2.0;

    double idwDenomFactor = upstreamFactor + centerFactor;

    upstreamFactor = upstreamFactor / idwDenomFactor;
    centerFactor   = centerFactor / idwDenomFactor;

    upstreamFactor = (1 + (1.0 / upstreamPecletNumber / upstreamFactor)) * upstreamFactor;
    centerFactor   = (1 - (1.0 / upstreamPecletNumber / centerFactor)) * centerFactor;

    incomingFlux =  upstreamFlow * (S[upstreamElement->index] * upstreamFactor +
                    S[index] * centerFactor) / volume;
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->soluteConcs[soluteIndex].value / volume;
  }


  if(hasDownstream && fabs(downstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;


    outgoingFlux = downstreamFlow * (S[downstreamElement->index] * downstreamFactor +
                   S[index] * centerFactor) / volume;
  }
  else if(downstreamPecletNumber >= 2.0)
  {
    outgoingFlux = flow * S[index] / volume;
  }
  else if(hasDownstream && downstreamPecletNumber <= -2.0)
  {
    outgoingFlux = downstreamFlow * S[downstreamElement->index] / volume;
  }
  else if(hasDownstream)
  {
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0;
    double centerFactor = 1.0 / length / 2.0;

    double idwDenomFactor = centerFactor + downstreamFactor;

    centerFactor /= idwDenomFactor;
    downstreamFactor   /= idwDenomFactor;

    centerFactor = (1 + (1.0 / downstreamPecletNumber/ centerFactor)) * centerFactor;
    downstreamFactor = (1 - (1.0 / downstreamPecletNumber / downstreamFactor)) * downstreamFactor;

    outgoingFlux =  downstreamFlow * (S[index] * centerFactor +
                                      S[downstreamElement->index] * downstreamFactor) / volume;
  }
  else
  {
    outgoingFlux = downstreamFlow * downstreamJunction->soluteConcs[soluteIndex].value / volume;
  }

  double DSoluteDt = incomingFlux - outgoingFlux;

  return DSoluteDt;
}

double Element::computeDSoluteDtTVD(double dt, double S[], int soluteIndex)
{

    return 0.0;
}

double Element::computeCourantFactor() const
{
  return  fabs(flow / xSectionArea / length);
}

double Element::computeDispersionFactor() const
{
  return 2.0 * longDispersion.value / (length * length);
}

void Element::computeLongDispersion()
{
  double vel = flow / xSectionArea;

  double fricVel = sqrt(9.81 * depth * slope);
  double dispFischer = 0.11 * vel * vel * width * width / (depth * fricVel);

  double dispNumerical = fabs( vel * length / 2.0);

  longDispersion.value = dispNumerical <= dispFischer ? dispFischer - dispNumerical : 0.0;
}

void Element::computePecletNumbers()
{
  pecletNumber = longDispersion.value ? (flow / xSectionArea) * length / longDispersion.value : flow * 10000 / fabs(flow);
}

void Element::computeUpstreamPeclet()
{
  computeUpstreamFlow();

  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {
    upstreamPecletNumber = upstreamJunction->longDispersion ?
                             upstreamVelocity * ((upstreamElement->length / 2.0) + (length / 2.0)) / (upstreamJunction->longDispersion)
                           : upstreamVelocity * 10000 / fabs(upstreamVelocity);
  }
  else
  {
    upstreamPecletNumber = pecletNumber;
  }
}

void Element::computeDownstreamPeclet()
{
  computeDownstreamFlow();

  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    downstreamPecletNumber = downstreamJunction->longDispersion ?
                               downstreamVelocity * ((downstreamElement->length / 2.0) + (length / 2.0)) / (downstreamJunction->longDispersion)
                             : downstreamVelocity * 10000 / fabs(downstreamVelocity);
  }
  else
  {
    downstreamPecletNumber = pecletNumber;
  }
}

void Element::computeHeatBalance()
{
  //KJ
  heatBalance +=  model->m_waterDensity * model->m_cp * xSectionArea * length * (temperature.value - prevTemperature.value) / 1000.000;
}

void Element::computeSoluteBalance(int soluteIndex)
{
  soluteMassBalance[soluteIndex] += xSectionArea * length * (soluteConcs[soluteIndex].value - prevSoluteConcs[soluteIndex].value);
}

void Element::setUpstreamElement()
{
  upstreamElement = nullptr;
  upstreamElementDirection = 1.0;

  if(upstreamJunction->junctionType == ElementJunction::DoubleElement)
  {
    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        upstreamElementDirection = 1.0;
        return;
      }
    }

    for(Element *element : upstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        upstreamElementDirection = -1.0;
        return;
      }
    }
  }
}

void Element::setDownStreamElement()
{
  downstreamElement = nullptr;
  downstreamElementDirection = 1.0;

  if(downstreamJunction->junctionType == ElementJunction::DoubleElement)
  {
    for(Element *element : downstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        downstreamElementDirection = 1.0;
        return;
      }
    }

    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        downstreamElementDirection = -1.0;
        return;
      }
    }
  }
}

void Element::computeUpstreamFlow()
{
  if(upstreamElement != nullptr)
  {
    upstreamFlow =  ((upstreamElement->flow * upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                     (flow / (length * 0.5))) /
                    ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamVelocity =  (((upstreamElement->flow / upstreamElement->xSectionArea )* upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                         ((flow / xSectionArea) / (length * 0.5))) /
                        ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));
  }
  else
  {
    upstreamFlow = flow;
    upstreamVelocity = upstreamFlow / xSectionArea;
  }
}

void Element::computeDownstreamFlow()
{
  if(downstreamElement != nullptr)
  {
    downstreamFlow = ((downstreamElement->flow * downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                      (flow / (length * 0.5))) /
                     ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamVelocity = (((downstreamElement->flow / downstreamElement->xSectionArea ) * downstreamElementDirection / (downstreamElement->length * 0.5)) +
                          ((flow / xSectionArea) / (length * 0.5))) /
                         ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));
  }
  else
  {
    downstreamFlow = flow;
    downstreamVelocity = downstreamFlow / xSectionArea;
  }
}
