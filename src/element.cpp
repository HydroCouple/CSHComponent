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
#include "elementadvupwind.h"
#include "elementadvcentral.h"
#include "elementadvhybrid.h"
#include "elementadvtvd.h"

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
    computeTempAdvDeriv(nullptr),
    computeSoluteAdvDeriv(nullptr),
    computeTempDispDeriv(nullptr),
    computeSoluteDispDeriv(nullptr)
{
  starting = true;

  initializeSolutes();

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

  computeTempAdvDeriv  = new ComputeTempAdvDeriv[2]();
  computeTempDispDeriv = new ComputeTempDeriv[2]();
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

    for(int i = 0; i < numSolutes; i++)
    {
      delete[] computeSoluteAdvDeriv[i];
      delete[] computeSoluteDispDeriv[i];
    }

    delete[] computeSoluteAdvDeriv;
    delete[] computeSoluteDispDeriv;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

  delete[] computeTempAdvDeriv;
  delete[] computeTempDispDeriv;
}

void Element::initialize()
{
  //set upstream and downstream elements
  setUpstreamElement();
  setDownStreamElement();

  switch (model->m_advectionMode)
  {
    case CSHModel::Central:
      {
        setAdvectionFuctions = &ElementAdvCentral::setAdvectionFunction;
      }
      break;
    case CSHModel::Hybrid:
      {
        setAdvectionFuctions = &ElementAdvHybrid::setAdvectionFunction;
      }
      break;
   case CSHModel::TVD:
      {
        setAdvectionFuctions = &ElementAdvTVD::setAdvectionFunction;
      }
      break;
    default:
      {
        setAdvectionFuctions = &ElementAdvUpwind::setAdvectionFunction;
      }
      break;
  }

  totalHeatBalance =  0.0;
  totalAdvDispHeatBalance = 0.0;
  totalEvaporativeHeatFluxesBalance = 0.0;
  totalConvectiveHeatFluxesBalance = 0.0;
  totalRadiationFluxesHeatBalance = 0.0;
  relativeHumidity = 0.0;
  evaporationRate = 0.0;
  saturationVaporPressureAir = 0.0;
  saturationVaporPressureWater = 0.0;
  vaporPressureAir = 0.0;
  vaporPressureWater = 0.0;
  windSpeed = 0.0;
  airTemperature = 0.0;
  evaporationHeatFlux = 0.0;
  convectionHeatFlux = 0.0;
  upstreamLongDispersion = 0.0;
  downstreamLongDispersion = 0.0;
  volume = prev_volume = 0.0;

  starting = true;

  dvolume_dt.isBC = false;
  dvolume_dt.value = 0.0;

  for(int i = 0; i < numSolutes; i++)
  {
    totalSoluteMassBalance[i] = 0.0;
    totalAdvDispSoluteMassBalance[i] = 0.0;
    totalExternalSoluteFluxesMassBalance[i] = 0.0;
  }

  setDispersionFunctions();
  (*setAdvectionFuctions)(this);
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

    for(int i = 0; i < numSolutes; i++)
    {
      delete[] computeSoluteAdvDeriv[i];
      delete[] computeSoluteDispDeriv[i];
    }

    delete[] computeSoluteAdvDeriv;
    delete[] computeSoluteDispDeriv;
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

    computeSoluteAdvDeriv  = new ComputeSoluteAdvDeriv*[numSolutes]();
    computeSoluteDispDeriv = new ComputeSoluteDeriv*[numSolutes]();

    for(int i = 0; i < numSolutes; i++)
    {
      computeSoluteAdvDeriv[i] = new ComputeSoluteAdvDeriv[2]();
      computeSoluteDispDeriv[i] = new ComputeSoluteDeriv[2]();
    }
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0.0;

  if(volume > 1e-14)
  {
    //Compute advection
    DTDt += computeDTDtAdv(dt, T);

    //Compute dispersion
    DTDt += computeDTDtDispersion(dt, T);

    //Add external sources
    {
      DTDt += radiationFluxes * top_area / rho_cp_vol;
      DTDt += externalHeatFluxes / rho_cp_vol;
    }

    //Evaporation and convection
    {
      DTDt += evaporationHeatFlux * top_area / rho_cp_vol;
      DTDt += convectionHeatFlux * top_area / rho_cp_vol;
    }

    //Product rule subtract volume derivative
    {
      DTDt -= temperature.value * dvolume_dt.value / volume;
    }
  }


  return DTDt;
}

double Element::computeDTDtAdv(double dt, double T[])
{
  double DTDt = 0.0;

  //Upstream
  DTDt += (*computeTempAdvDeriv[0])(this, dt, T);

  //Downstream
  DTDt += (*computeTempAdvDeriv[1])(this, dt, T);

  DTDt = DTDt / rho_cp_vol;

  return DTDt;
}

double Element::computeDTDtDispersion(double dt, double T[])
{
  double DTDt = 0.0;

  //Upstream
  DTDt += (this->*(computeTempDispDeriv[0]))(dt, T);

  //Downstream
  DTDt += (this->*(computeTempDispDeriv[1]))(dt, T);

  DTDt = DTDt / rho_cp_vol;

  return DTDt;
}

double Element::computeDTDtDispersionUpstreamJunction(double dt, double T[])
{
  double DTDt = upstreamLongDispersion * upstreamXSectionArea * rho_cp *
                (upstreamJunction->temperature.value - T[index]) /
                 (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionUpstreamNeighbour(double dt, double T[])
{
  double DTDt = upstreamLongDispersion * upstreamXSectionArea * rho_cp *
                (T[upstreamElement->index]  - T[index]) /
      ((length / 2.0) + (upstreamElement->length / 2.0));

  return DTDt;
}

double Element::computeDTDtDispersionDownstreamJunction(double dt, double T[])
{
  double DTDt = downstreamLongDispersion * downstreamXSectionArea * rho_cp *
                (downstreamJunction->temperature.value  - T[index]) /
                (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionDownstreamNeighbour(double dt, double T[])
{
  double DTDt = downstreamLongDispersion * downstreamXSectionArea * rho_cp *
                (T[downstreamElement->index]  - T[index]) /
                ((length / 2.0) + (downstreamElement->length / 2.0));

  return DTDt;
}

void Element::computeEvaporation()
{
  if(model->m_useEvaporation)
  {
    double Le = 1000.0 * (2499.0 - 2.36 * temperature.value);

    saturationVaporPressureAir = 0.61275 * exp(17.27 * airTemperature / (237.3 + airTemperature));

    saturationVaporPressureWater = 0.61275 * exp(17.27 * temperature.value / (237.3 + temperature.value));

    vaporPressureAir = relativeHumidity * saturationVaporPressureAir / 100.0;

    vaporPressureWater = relativeHumidity * saturationVaporPressureWater / 100.0;

    windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);

    evaporationRate = windFunction * (saturationVaporPressureWater - vaporPressureAir);

    evaporationHeatFlux = -Le * evaporationRate * model->m_waterDensity;
  }
}

void Element::computeConvection()
{
  if(model->m_useConvection)
  {
    double Le = 1000.0 * (2499.0 - 2.36 * temperature.value);
    windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);
    convectionHeatFlux = - model->m_waterDensity * Le * windFunction * model->m_bowensCoeff * model->m_pressureRatio * (temperature.value - airTemperature) ;
  }
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;

  if(volume > 1e-14)
  {
    //Compute advection
    DSoluteDt += computeDSoluteDtAdv(dt, S, soluteIndex);

    //Compute dispersion
    DSoluteDt += computeDSoluteDtDispersion(dt, S, soluteIndex);

    //First order reaction reaction
    DSoluteDt += model->m_solute_first_order_k[soluteIndex] * soluteConcs[soluteIndex].value;

    //subtract chain rule volume derivative
    {
      DSoluteDt -= (soluteConcs[soluteIndex].value * dvolume_dt.value) / volume;
    }

    //Add external sources
    {
      DSoluteDt += externalSoluteFluxes[soluteIndex] / volume;
    }
  }

  return DSoluteDt;
}

double Element::computeDSoluteDtAdv(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  //Upstream
  DSoluteDt += (*computeSoluteAdvDeriv[soluteIndex][0])(this, dt, S, soluteIndex);

  //Downstream
  DSoluteDt += (*computeSoluteAdvDeriv[soluteIndex][1])(this, dt, S, soluteIndex);

  DSoluteDt = DSoluteDt / volume;

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersion(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  //Upstream
  DSoluteDt += (this->*computeSoluteDispDeriv[soluteIndex][0])(dt, S, soluteIndex);

  //Downstream
  DSoluteDt += (this->*computeSoluteDispDeriv[soluteIndex][1])(dt, S, soluteIndex);


  DSoluteDt = DSoluteDt / volume;

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionUpstreamJunction(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = upstreamLongDispersion * upstreamXSectionArea *
                     (upstreamJunction->soluteConcs[soluteIndex].value - S[index]) /
                     (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionUpstreamNeighbour(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = upstreamLongDispersion * upstreamXSectionArea *
                     (S[upstreamElement->index]  - S[index]) /
      ((length / 2.0) + (upstreamElement->length / 2.0));

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionDownstreamJunction(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = downstreamLongDispersion * downstreamXSectionArea *
                     (downstreamJunction->soluteConcs[soluteIndex].value  - S[index]) /
                     (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionDownstreamNeighbour(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = downstreamLongDispersion * downstreamXSectionArea *
                     (S[downstreamElement->index]  - S[index]) /
                     ((length / 2.0) +  (downstreamElement->length / 2.0));

  return DSoluteDt;
}

void Element::setDispersionFunctions()
{
  if(upstreamJunction->temperature.isBC)
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamJunction;
  }
  else if(upstreamElement)
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamNeighbour;
  }
  else
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamJunction;
  }

  if(downstreamJunction->temperature.isBC)
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamJunction;
  }
  else if(downstreamElement)
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamNeighbour;
  }
  else
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamJunction;
  }

  for(int i = 0; i < numSolutes; i++)
  {
    if(upstreamJunction->soluteConcs[i].isBC)
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamJunction;
    }
    else if(upstreamElement)
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamNeighbour;
    }
    else
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamJunction;
    }

    if(downstreamJunction->soluteConcs[i].isBC)
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamJunction;
    }
    else if(downstreamElement)
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamNeighbour;
    }
    else
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamJunction;
    }
  }
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
  else if(dvolume_dt.isBC == false)
  {
    prev_volume = volume;
    volume  = xSectionArea * length;
    dvolume_dt.value = (volume - prev_volume) / model->m_prevTimeStep;
  }

  rho_vol = model->m_waterDensity * volume;
  rho_cp  = model->m_waterDensity * model->m_cp;
  rho_cp_vol = rho_cp * volume;
  top_area = length * width;

  //  setDispersionFunctions();
  (*setAdvectionFuctions)(this);
}

void Element::computeLongDispersion()
{
  double vel = flow / xSectionArea;
  slope = max(0.00001, fabs(upstreamJunction->z - downstreamJunction->z) / length);

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
                                                      : upstreamVelocity * 10000 / upstreamVelocity;
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
                                                          : downstreamVelocity * 10000 / downstreamVelocity;
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
