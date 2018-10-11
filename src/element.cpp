/*!
*  \file    element.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
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
#include "cshmodel.h"

#include <math.h>

using namespace std;

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  CSHModel *model)
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
    relativeHumidity(0.0),
    evaporationRate(0.0),
    evaporationHeatFlux(0.0),
    saturationVaporPressureAir(0.0),
    saturationVaporPressureWater(0.0),
    vaporPressureAir(0.0),
    vaporPressureWater(0.0),
    windSpeed(0.0),
    airTemperature(0.0),
    convectionHeatFlux(0.0),
    externalHeatFluxes(0.0),
    radiationFluxes(0.0),
    externalSoluteFluxes(nullptr),
    totalHeatBalance(0.0),
    totalAdvDispHeatBalance(0.0),
    totalEvaporativeHeatFluxesBalance(0.0),
    totalConvectiveHeatFluxesBalance(0.0),
    totalRadiationFluxesHeatBalance(0.0),
    totalExternalHeatFluxesBalance(0.0),
    totalSoluteMassBalance(nullptr),
    totalAdvDispSoluteMassBalance(nullptr),
    totalExternalSoluteFluxesMassBalance(nullptr),
    pecletNumber(0.0),
    upstreamElement(nullptr),
    downstreamElement(nullptr),
    model(model),
    downstreamPecletNumber(1.0),
    downstreamFlow(0.0),
    downstreamVelocity(0.0),
    upstreamPecletNumber(1.0),
    upstreamFlow(0.0),
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
    delete[] totalSoluteMassBalance;
    delete[] totalAdvDispSoluteMassBalance;
    delete[] totalExternalSoluteFluxesMassBalance;
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
    case CSHModel::AdvectionDiscretizationMode::Central:
      {
        computeTempAdv = &Element::computeDTDtCentral;
        computeSoluteAdv = &Element::computeDSoluteDtCentral;
      }
      break;
    case CSHModel::AdvectionDiscretizationMode::Hybrid:
      {
        computeTempAdv = &Element::computeDTDtHybrid;
        computeSoluteAdv = &Element::computeDSoluteDtHybrid;
      }
      break;
    case CSHModel::AdvectionDiscretizationMode::TVD:
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

  totalHeatBalance =  totalAdvDispHeatBalance =
      totalEvaporativeHeatFluxesBalance = totalConvectiveHeatFluxesBalance =
      totalRadiationFluxesHeatBalance = relativeHumidity =
      evaporationRate = saturationVaporPressureAir = saturationVaporPressureWater =
      vaporPressureAir = vaporPressureWater = windSpeed = airTemperature =
      evaporationHeatFlux = convectionHeatFlux = upstreamLongDispersion = downstreamLongDispersion =
      volume = prev_volume = 0.0;

  starting = true;

  for(int i = 0; i < numSolutes; i++)
  {
    totalSoluteMassBalance[i] = 0.0;
    totalAdvDispSoluteMassBalance[i] = 0.0;
    totalExternalSoluteFluxesMassBalance[i] = 0.0;
  }
}

void Element::initializeSolutes()
{
  starting = true;

  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalAdvDispSoluteMassBalance; totalAdvDispSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
  }

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    totalSoluteMassBalance = new double[numSolutes]();
    totalAdvDispSoluteMassBalance = new double[numSolutes]();
    totalExternalSoluteFluxesMassBalance = new double[numSolutes]();
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0.0;

  if(volume > 1e-10)
  {
    //Compute advection
    DTDt += (this->*computeTempAdv)(dt, T);

    //Compute dispersion
    DTDt += computeDTDtDispersion(dt, T);

    //Add external sources
    {
      DTDt += radiationFluxes * top_area / rho_cp_vol;
      DTDt += externalHeatFluxes / rho_cp_vol;
    }

    //Chain rule subtract volume derivative
    {
      DTDt -= (rho_cp * T[index] * (volume - prev_volume)) / (model->m_timeStep * rho_cp_vol);
    }

    if(model->m_useEvaporation)
    {
      DTDt += computeDTDtEvaporation(dt, T);

      if(model->m_useConvection)
      {
        DTDt += computeDTDtConvection(dt, T);
      }
    }
  }


  return DTDt;
}

double Element::computeDTDtUpwind(double dt, double T[])
{
  double incomingFlux = 0.0;
  double outgoingFlux = 0.0;

  //Flow goes from upstream to downstream
  if(flow >= 0)
  {
    if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
    {
      incomingFlux = rho_cp * upstreamElement->flow * upstreamElementDirection * T[upstreamElement->index];
    }
    else
    {
      incomingFlux = rho_cp * flow * upstreamJunction->temperature.value;
    }

    outgoingFlux = rho_cp * flow * T[index];
  }
  //Opposite Direction
  else
  {
    incomingFlux = rho_cp * flow * T[index];

    if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
    {
      outgoingFlux = rho_cp * downstreamElement->flow * downstreamElementDirection * T[downstreamElement->index];
    }
    else
    {
      outgoingFlux = rho_cp * flow * downstreamJunction->temperature.value;
    }
  }

  return (incomingFlux - outgoingFlux) / rho_cp_vol ;
}

double Element::computeDTDtCentral(double dt, double T[])
{
  double incomingFlux = 0;
  double outgoingFlux = 0;

  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = rho_cp * upstreamFlow * (T[upstreamElement->index] * centerFactor +
                   T[index] * upstreamFactor);
  }
  else
  {
    incomingFlux = rho_cp * flow * upstreamJunction->temperature.value / 2.0;
  }

  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;

    outgoingFlux = rho_cp * downstreamFlow * (T[downstreamElement->index] * downstreamFactor +
                   T[index] * centerFactor);
  }
  else
  {
    outgoingFlux = rho_cp * flow * downstreamJunction->temperature.value / 2.0;
  }

  return (incomingFlux - outgoingFlux) / rho_cp_vol;
}

double Element::computeDTDtHybrid(double dt, double T[])
{
  double incomingFlux = 0;
  double outgoingFlux = 0;

  bool hasUpstream = upstreamElement != nullptr && !upstreamJunction->temperature.isBC;
  bool hasDownstream = downstreamElement != nullptr && !downstreamJunction->temperature.isBC;

  //If peclet number is zero use central differencing
  if(hasUpstream && fabs(upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = rho_cp * upstreamFlow * (T[upstreamElement->index] * centerFactor +
                   T[index] * upstreamFactor);
  }
  else if(hasUpstream && upstreamPecletNumber >= 2.0)
  {
    incomingFlux = rho_cp * upstreamElement->flow * upstreamElementDirection * T[upstreamElement->index];
  }
  else if(upstreamPecletNumber <= -2.0)
  {
    incomingFlux = rho_cp * flow * T[index];
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

    incomingFlux = rho_cp * upstreamFlow * (T[upstreamElement->index] * upstreamFactor +
                   T[index] * centerFactor);
  }
  else
  {
    incomingFlux = rho_cp * flow * upstreamJunction->temperature.value;
  }


  if(hasDownstream && fabs(downstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;

    outgoingFlux = rho_cp * downstreamFlow * (T[downstreamElement->index] * downstreamFactor +
                   T[index] * centerFactor);
  }
  else if(downstreamPecletNumber >= 2.0)
  {
    outgoingFlux = rho_cp * flow * T[index];
  }
  else if(hasDownstream && downstreamPecletNumber <= -2.0)
  {
    outgoingFlux = rho_cp * downstreamElement->flow *
                   downstreamElementDirection * T[downstreamElement->index];
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

    outgoingFlux = rho_cp * downstreamFlow * (T[index] * centerFactor +
                                              T[downstreamElement->index] * downstreamFactor);
  }
  else
  {
    outgoingFlux = rho_cp * flow * downstreamJunction->temperature.value;
  }

  double DTDt = (incomingFlux - outgoingFlux)/ rho_cp_vol;

  return DTDt;
}

double Element::computeDTDtTVD(double dt, double T[])
{
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

  //  if(flow >= 0)
  //  {
  //    incomingFlux = downstreamFlow * (T[index] - min(f_pos, f_neg)) / volume;
  //  }
  //  else
  //  {
  //    incomingFlux = downstreamFlow * (T[index] + min(f_pos, f_neg)) / volume;
  //  }

  return incomingFlux - outgoingFlux;
}

double Element::computeDTDtULTIMATE(double dt, double T[])
{

}

double Element::computeDTDtDispersion(double dt, double T[])
{
  double DTDt = 0.0;

  //  printf("l: %f\tx: %f\tv: %f\n" , length, xSectionArea, volume);
  //If upstream element exists and upstream junction is not a boundary condition
  //use upstream element to compute derivative.
  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {

    DTDt -= rho_cp * upstreamLongDispersion * upstreamXSectionArea *
            (T[index] - T[upstreamElement->index]) / ((length + upstreamElement->length) / 2.0);

  }
  //Otherwise use upstream junction to compute temperature derivate
  else
  {
    DTDt -=  rho_cp * upstreamLongDispersion * upstreamXSectionArea *
             (T[index] - upstreamJunction->temperature.value) / (length / 2.0);
  }

  //If downstream element exists and downstream junction element is not a boundary condition
  //use downstream element to compute derivative.
  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    DTDt += rho_cp * downstreamLongDispersion * downstreamXSectionArea *
            (T[downstreamElement->index] - T[index]) / ((length + downstreamElement->length) / 2.0);
  }
  //Otherwise use downstream element junction
  else
  {
    DTDt += rho_cp * downstreamLongDispersion * downstreamXSectionArea *
            (downstreamJunction->temperature.value - T[index]) / (length / 2.0);
  }

  return DTDt / rho_cp_vol;
}

double Element::computeDTDtEvaporation(double dt, double T[])
{
  double currentTemp = T[index];
  double DTDt = 0.0;

  {
    double Le = 1000.0 * (2499.0 - 2.36 * currentTemp);

    saturationVaporPressureAir = 0.61275 * exp(17.27 * airTemperature / (237.3 + airTemperature));

    saturationVaporPressureWater = 0.61275 * exp(17.27 * currentTemp / (237.3 + currentTemp));

    vaporPressureAir = relativeHumidity * saturationVaporPressureAir / 100.0;

    vaporPressureWater = relativeHumidity * saturationVaporPressureWater / 100.0;

    windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);

    //    windFunction = 0.069 + model->m_evapWindFuncCoeffB * fabs(windSpeed);

    evaporationRate = windFunction * (saturationVaporPressureWater - vaporPressureAir);

    evaporationHeatFlux = -Le * evaporationRate * model->m_waterDensity;

    DTDt = evaporationHeatFlux * top_area / rho_cp_vol;
  }

  //  {
  //    DTDt = 0.0;
  //    double esatair = 4.596 * std::exp(17.270 * airTemperature / (237.30 + airTemperature));
  //    double eair = relativeHumidity / 100.0  * esatair ;
  //    //    double fUw =  model->m_evapWindFuncCoeffA +  model->m_evapWindFuncCoeffB * (windSpeed * windSpeed);
  //    double fUw =  19.0 +  0.95 * (windSpeed * windSpeed);
  //    double esat = 4.596 * std::exp(17.270 * currentTemp / (237.30 + currentTemp));
  //    evaporationHeatFlux = -fUw * (esat - eair) * 4.184 * 10000.0  / 86400.0;
  //    DTDt = evaporationHeatFlux * length * width / (model->m_waterDensity * model->m_cp * xSectionArea * length);
  //  }

  return DTDt;
}

double Element::computeDTDtConvection(double dt, double T[])
{
  double currentTemp = T[index];
  double DTDt = 0.0;

  {
    double Le = 1000.0 * (2499.0 - 2.36 * currentTemp);
    convectionHeatFlux = - model->m_waterDensity * Le * windFunction * model->m_bowensCoeff * model->m_pressureRatio * (currentTemp - airTemperature) ;
    DTDt = convectionHeatFlux * width /  (rho_cp * xSectionArea);
  }

  //  {
  //    DTDt = 0.0;
  //    //    double fUw = model->m_evapWindFuncCoeffA +  model->m_evapWindFuncCoeffB * (windSpeed * windSpeed);
  //    double fUw =  19.0 +  0.95 * (windSpeed * windSpeed);
  //    convectionHeatFlux = -0.6 * fUw * (currentTemp - airTemperature) * 4.184 * 10000.0  / 86400.0;
  //    DTDt = convectionHeatFlux * length * width /  (model->m_waterDensity * model->m_cp * xSectionArea * length);
  //  }

  return DTDt;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;

  if(volume > 1e-10)
  {
    //Compute advection
    DSoluteDt += (this->*computeSoluteAdv)(dt, S, soluteIndex);

    //Compute dispersion
    DSoluteDt += computeDSoluteDtDispersion(dt, S, soluteIndex);

    //subtract chain rule volume derivative
    {
      DSoluteDt -= (S[index] * (volume - prev_volume)) / (model->m_timeStep * volume);
    }

    //Add external sources
    {
      DSoluteDt += externalSoluteFluxes[soluteIndex] / volume;
    }
  }

  return DSoluteDt;
}

double Element::computeDSoluteDtUpwind(double dt, double S[], int soluteIndex)
{
  double incomingFlux = 0.0;
  double outgoingFlux = 0.0;

  //Flow goes from upstream to downstream
  if(flow >= 0)
  {
    if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
    {
      incomingFlux = upstreamElement->flow *
                     upstreamElementDirection * S[upstreamElement->index];
    }
    else
    {
      incomingFlux = flow * upstreamJunction->soluteConcs[soluteIndex].value;
    }

    outgoingFlux = flow * S[index];

  }
  //Opposite Direction
  else
  {
    incomingFlux = flow * S[index];

    if(downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC)
    {
      outgoingFlux = downstreamElement->flow *
                     downstreamElementDirection * S[downstreamElement->index];
    }
    else
    {
      outgoingFlux = flow * downstreamJunction->soluteConcs[soluteIndex].value;
    }
  }

  return (incomingFlux - outgoingFlux) / volume;
}

double Element::computeDSoluteDtCentral(double dt, double S[], int soluteIndex)
{
  double incomingFlux = 0;
  double outgoingFlux = 0;

  if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (S[upstreamElement->index] * centerFactor +
                   S[index] * upstreamFactor);
  }
  else
  {
    incomingFlux = flow * upstreamJunction->soluteConcs[soluteIndex].value;
  }

  if(downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;

    outgoingFlux = downstreamFlow * (S[downstreamElement->index] * downstreamFactor +
                   S[index] * centerFactor);
  }
  else
  {
    outgoingFlux = flow * downstreamJunction->soluteConcs[soluteIndex].value;
  }

  return (incomingFlux - outgoingFlux) / volume;
}

double Element::computeDSoluteDtHybrid(double dt, double S[], int soluteIndex)
{
  double incomingFlux = 0;
  double outgoingFlux = 0;

  bool hasUpstream = upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC;
  bool hasDownstream = downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC;

  //If peclet number is zero use central differencing
  if(hasUpstream && fabs(upstreamPecletNumber - 0.0) < std::numeric_limits<double>::epsilon())
  {
    double denom = (1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double upstreamFactor = 1.0 / upstreamElement->length / 2.0 / denom;

    incomingFlux = upstreamFlow * (S[upstreamElement->index] * centerFactor +
                   S[index] * upstreamFactor);
  }
  else if(hasUpstream && upstreamPecletNumber >= 2.0)
  {
    incomingFlux = upstreamElement->flow *
                   upstreamElementDirection * S[upstreamElement->index];
  }
  else if(upstreamPecletNumber <= -2.0)
  {
    incomingFlux = flow * S[index];
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
                    S[index] * centerFactor);
  }
  else
  {
    incomingFlux = flow * upstreamJunction->soluteConcs[soluteIndex].value;
  }


  if(hasDownstream && fabs(downstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
  {
    double denom = (1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0);
    double centerFactor = 1.0 / length / 2.0 / denom;
    double downstreamFactor = 1.0 / downstreamElement->length / 2.0 / denom;


    outgoingFlux = downstreamFlow * (S[downstreamElement->index] * downstreamFactor +
                   S[index] * centerFactor);
  }
  else if(downstreamPecletNumber >= 2.0)
  {
    outgoingFlux = flow * S[index];
  }
  else if(hasDownstream && downstreamPecletNumber <= -2.0)
  {
    outgoingFlux = downstreamElement->flow * downstreamElementDirection * S[downstreamElement->index];
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
                                      S[downstreamElement->index] * downstreamFactor);
  }
  else
  {
    outgoingFlux = flow * downstreamJunction->soluteConcs[soluteIndex].value;
  }

  double DSoluteDt = (incomingFlux - outgoingFlux) / volume;

  return DSoluteDt;
}

double Element::computeDSoluteDtTVD(double dt, double S[], int soluteIndex)
{

  return 0.0;
}

double Element::computeDSoluteDtDispersion(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  //If upstream element exists and upstream junction is not a boundary condition
  //use upstream element to compute derivative.
  if(upstreamElement != nullptr && !upstreamJunction->soluteConcs[soluteIndex].isBC)
  {

    DSoluteDt -= upstreamLongDispersion * upstreamXSectionArea *
                 (S[index] - S[upstreamElement->index]) / ((length + upstreamElement->length) / 2.0);

  }
  //Otherwise use upstream junction to compute temperature derivate
  else
  {
    DSoluteDt -=  upstreamLongDispersion * upstreamXSectionArea *
                  (S[index] - upstreamJunction->soluteConcs[soluteIndex].value) / (length / 2.0);
  }

  //If downstream element exists and downstream junction element is not a boundary condition
  //use downstream element to compute derivative.
  if(downstreamElement != nullptr && !downstreamJunction->soluteConcs[soluteIndex].isBC)
  {
    DSoluteDt += downstreamLongDispersion * downstreamXSectionArea *
                 (S[downstreamElement->index] - S[index]) / ((length + downstreamElement->length) / 2.0);
  }
  //Otherwise use downstream element junction
  else
  {
    DSoluteDt += downstreamLongDispersion * downstreamXSectionArea *
                 (downstreamJunction->soluteConcs[soluteIndex].value - S[index]) / (length / 2.0);
  }

  DSoluteDt = DSoluteDt / volume;

  return DSoluteDt;
}

double Element::computeCourantFactor() const
{
  return  fabs(flow / xSectionArea / length);
}

double Element::computeDispersionFactor() const
{
  return 2.0 * longDispersion.value / (length * length);
}

void Element::computeDerivedHydraulics()
{
  if(starting)
  {
    prev_volume = volume  = xSectionArea * length;
    starting = false;
  }
  else
  {
    prev_volume = volume;
    volume  = xSectionArea * length;
  }

  volume  = xSectionArea * length;
  rho_vol = model->m_waterDensity * volume;
  rho_cp  = model->m_waterDensity * model->m_cp;
  rho_cp_vol = rho_cp * volume;
  top_area = length * width;
}

void Element::computeLongDispersion()
{
  double vel = flow / xSectionArea;
  slope = max(0.000001, fabs(upstreamJunction->z - downstreamJunction->z) / length);

  double fricVel = sqrt(9.81 * depth * slope);
  double dispFischer = (0.011 * vel * vel * width * width) / (depth * fricVel);
  double dispNumerical = fabs( vel * length / 2.0);

  longDispersion.value = dispNumerical <= dispFischer ? dispFischer - dispNumerical : 0.0;
}

void Element::computePecletNumbers()
{
  pecletNumber = longDispersion.value ? (flow / xSectionArea) * length / longDispersion.value : flow * 10000 / fabs(flow);
}

void Element::computeUpstreamPeclet()
{
  if(upstreamElement != nullptr)
  {
    upstreamFlow =  ((upstreamElement->flow * upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                     (flow / (length * 0.5))) /
                    ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamVelocity =  (((upstreamElement->flow / upstreamElement->xSectionArea )* upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                         ((flow / xSectionArea) / (length * 0.5))) /
                        ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamXSectionArea = ((upstreamElement->xSectionArea / (upstreamElement->length * 0.5)) + (xSectionArea / (length * 0.5))) /
                           ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamLongDispersion =  ((upstreamElement->longDispersion.value / (upstreamElement->length * 0.5)) + (longDispersion.value / (length * 0.5))) /
                              ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamPecletNumber = upstreamLongDispersion > 0 ? upstreamVelocity * ((upstreamElement->length / 2.0) + (length / 2.0)) / (upstreamLongDispersion)
                                                      : upstreamVelocity * 10000 / fabs(upstreamVelocity);
  }
  else
  {
    upstreamFlow = flow;
    upstreamVelocity = upstreamFlow / xSectionArea;
    upstreamXSectionArea = xSectionArea;
    upstreamLongDispersion = longDispersion.value;
    upstreamPecletNumber = pecletNumber;
  }

}

void Element::computeDownstreamPeclet()
{
  if(downstreamElement != nullptr)
  {
    downstreamFlow = ((downstreamElement->flow * downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                      (flow / (length * 0.5))) /
                     ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamVelocity = (((downstreamElement->flow / downstreamElement->xSectionArea )* downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                          ((flow / xSectionArea) / (length * 0.5))) /
                         ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamXSectionArea = ((downstreamElement->xSectionArea / (downstreamElement->length * 0.5)) + (xSectionArea / (length * 0.5))) /
                             ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamLongDispersion = ((downstreamElement->longDispersion.value / (downstreamElement->length * 0.5)) + (longDispersion.value / (length * 0.5))) /
                               ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamPecletNumber = downstreamLongDispersion > 0 ? downstreamVelocity * ((downstreamElement->length / 2.0) + (length / 2.0)) / (downstreamLongDispersion)
                                                          : downstreamVelocity * 10000 / fabs(downstreamVelocity);
  }
  else
  {
    downstreamFlow = flow;
    downstreamVelocity = downstreamFlow / xSectionArea;
    downstreamXSectionArea = xSectionArea;
    downstreamLongDispersion = longDispersion.value;
    downstreamPecletNumber = pecletNumber;
  }
}

void Element::computeHeatBalance(double timeStep)
{
  double radiationEnergy = radiationFluxes * top_area * timeStep / 1000.0;
  totalRadiationFluxesHeatBalance += radiationEnergy;

  double externalEnergy = externalHeatFluxes * xSectionArea * length * timeStep / 1000.0;
  totalExternalHeatFluxesBalance += externalEnergy;

  double totalHeatEnergy = model->m_waterDensity * model->m_cp * xSectionArea * length * (temperature.value - prevTemperature.value) / 1000.0;
  totalHeatBalance +=  totalHeatEnergy;

  double totalEvaporationHeat = evaporationHeatFlux * top_area * timeStep / 1000.0;
  totalEvaporativeHeatFluxesBalance += totalEvaporationHeat;

  double totalConvectiveHeat = convectionHeatFlux * top_area * timeStep / 1000.0;
  totalConvectiveHeatFluxesBalance += totalConvectiveHeat;

  double advDispEnergy = totalHeatEnergy - radiationEnergy - externalEnergy - totalEvaporationHeat - totalConvectiveHeat;
  totalAdvDispHeatBalance += advDispEnergy;

}

void Element::computeSoluteBalance(double timeStep, int soluteIndex)
{

  double externalMass = externalSoluteFluxes[soluteIndex] * xSectionArea * length * timeStep;
  totalExternalSoluteFluxesMassBalance[soluteIndex] += externalMass;

  double totalMass = model->m_waterDensity * xSectionArea * length * (soluteConcs[soluteIndex].value - prevSoluteConcs[soluteIndex].value) ;
  totalSoluteMassBalance[soluteIndex] +=  totalMass;

  double advDispMass = totalMass - externalMass;
  totalAdvDispSoluteMassBalance[soluteIndex] += advDispMass;

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
