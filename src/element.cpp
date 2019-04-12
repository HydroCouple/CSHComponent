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

  sideSlopes = new double[2]();
  externalFlows = 0.0;
}

Element::~Element()
{


  deleteSoluteVariables();

  delete[] sideSlopes;

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

  slope = max(0.00001, fabs(upstreamJunction->z - downstreamJunction->z) / length);

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
  externalFlows = 0.0;

  starting = true;

  dvolume_dt.isBC = false;
  dvolume_dt.value = 0.0;

  xSectionArea = getAofH(depth);
  width = getWofH(depth);
  volume = xSectionArea * length;

  for(int i = 0; i < numSolutes; i++)
  {
    totalSoluteMassBalance[i] = 0.0;
    totalAdvDispSoluteMassBalance[i] = 0.0;
    totalExternalSoluteFluxesMassBalance[i] = 0.0;
  }

  setDispersionFunctions();
  (*setAdvectionFuctions)(this);

  if(model->m_solveHydraulics)
  {
    depth = getHofQ(flow.value);
    xSectionArea = getAofH(depth);
  }
}

void Element::initializeSolutes()
{
  starting = true;

  deleteSoluteVariables();

  if(model->m_solutes.size() > 0)
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
    totalSoluteMassBalance = new double[numSolutes]();
    totalAdvDispSoluteMassBalance = new double[numSolutes]();
    totalExternalSoluteFluxesMassBalance = new double[numSolutes]();
    sIndex = new int[numSolutes]();
    computeSoluteAdvDeriv  = new ComputeSoluteAdvDeriv*[numSolutes]();
    computeSoluteDispDeriv = new ComputeSoluteDeriv*[numSolutes]();

    for(int i = 0; i < numSolutes; i++)
    {
      computeSoluteAdvDeriv[i] = new ComputeSoluteAdvDeriv[2]();
      computeSoluteDispDeriv[i] = new ComputeSoluteDeriv[2]();
    }
  }
}

double Element::computeDADt(double dt, double A[])
{

  double DADt = 0.0;

  if(volume)
  {
    double A1 = getAofQ(upstreamJunction->inflow.value);
    xSectionArea = A[hIndex];

    depth = getHofA(xSectionArea);
    flow.value = getQofH(depth);
    width  = getWofH(depth);
    volume = xSectionArea * length;
    dvolume_dt.value = (volume - prev_volume) / model->m_timeStep;

    DADt += (A1 * upstreamJunction->inflow.value - xSectionArea * flow.value) / volume;
//    DADt -=  xSectionArea * dvolume_dt.value / volume;
    DADt += externalFlows / length;
  }



  return DADt;
}

double Element::getAofH(double H)
{
  double area = H * (bottomWidth + 0.5 * sideSlopes[0] * H + 0.5 * sideSlopes[1] * H); //fix for side slope
//  double area = H * bottomWidth; //fix for side slope

  return area;
}

double Element::getPofH(double H)
{
  double per = bottomWidth + H * sqrt(1 + sideSlopes[0] * sideSlopes[0]) + H * sqrt(1 + sideSlopes[1] * sideSlopes[1]);
//  double per = bottomWidth + H;
  return per;
}

double Element::getWofH(double H)
{
  double w = bottomWidth + sideSlopes[0] * H + sideSlopes[1] * H;
//  double w = bottomWidth;
  return w;
}

double Element::getHofA(double A)
{
  double h = findRoots(depth, A, &Element::getAofH);
//  double h = A / bottomWidth;
  return h;
}

double Element::getQofH(double H)
{
  double area = getAofH(H);
  double per = getPofH(H);
  double flow = pow(area / per, 2.0/3.0) * sqrt(slope) * area / mannings;
  return flow;
}

double Element::getAofQ(double Q)
{
  double area = getAofH(getHofQ(Q));

  return area;
}

double Element::getHofQ(double Q)
{
  double h = findRoots(depth, Q, &Element::getQofH);

  return h;
}

double Element::findRoots(double x, double y, GetXofY function, int maxIters, double derivStepSize, double eps)
{
  double error = 0;
  int iters = 0;

  do
  {

    double fplus = (this->*function)(x + derivStepSize) - y;
    double fminus = (this->*function)(x - derivStepSize) - y;

    double dfx = (fplus - fminus) /  (2 * derivStepSize);

    double fx =  (this->*function)(x) - y;

    double xn = x - fx / dfx;

    error = fabs(xn - x);

    x = xn;

    iters++;

  } while (error > eps && iters > 1 && iters < maxIters);


  return x;
}

void Element::computeHydraulicVariables()
{
  depth = getHofA(xSectionArea);
  flow.value = getQofH(depth);
  width  = getWofH(depth);
  volume = xSectionArea * length;
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0.0;

  if(volume > 1e-18)
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
                (T[upstreamJunction->tIndex] - T[tIndex]) /
      (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionUpstreamJunctionBC(double dt, double T[])
{
  double DTDt = upstreamLongDispersion * upstreamXSectionArea * rho_cp *
                (upstreamJunction->temperature.value - T[tIndex]) /
                (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionUpstreamNeighbour(double dt, double T[])
{
  double DTDt = upstreamLongDispersion * upstreamXSectionArea * rho_cp *
                (T[upstreamElement->tIndex]  - T[tIndex]) /
      ((length / 2.0) + (upstreamElement->length / 2.0));

  return DTDt;
}

double Element::computeDTDtDispersionDownstreamJunction(double dt, double T[])
{
  double DTDt = downstreamLongDispersion * downstreamXSectionArea * rho_cp *
                (T[downstreamJunction->tIndex]  - T[tIndex]) /
      (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionDownstreamJunctionBC(double dt, double T[])
{
  double DTDt = downstreamLongDispersion * downstreamXSectionArea * rho_cp *
                (downstreamJunction->temperature.value  - T[tIndex]) /
                (length / 2.0);

  return DTDt;
}

double Element::computeDTDtDispersionDownstreamNeighbour(double dt, double T[])
{
  double DTDt = downstreamLongDispersion * downstreamXSectionArea * rho_cp *
                (T[downstreamElement->tIndex]  - T[tIndex]) /
      ((length / 2.0) + (downstreamElement->length / 2.0));

  return DTDt;
}

double Element::computeDTDtDispersionSelf(double dt, double T[])
{
  return 0.0;
}

void Element::computeEvaporation()
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

void Element::computeConvection()
{
  double Le = 1000.0 * (2499.0 - 2.36 * temperature.value);
  windFunction = model->m_evapWindFuncCoeffA + model->m_evapWindFuncCoeffB * fabs(windSpeed);
  convectionHeatFlux = - model->m_waterDensity * Le * windFunction * model->m_bowensCoeff * model->m_pressureRatio * (temperature.value - airTemperature);
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;

  if(volume > 1e-18)
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
                     (S[upstreamJunction->sIndex[soluteIndex]] - S[sIndex[soluteIndex]]) /
      (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionUpstreamJunctionBC(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = upstreamLongDispersion * upstreamXSectionArea *
                     (upstreamJunction->soluteConcs[soluteIndex].value - S[sIndex[soluteIndex]]) /
      (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionUpstreamNeighbour(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = upstreamLongDispersion * upstreamXSectionArea *
                     (S[upstreamElement->sIndex[soluteIndex]]  - S[sIndex[soluteIndex]]) /
      ((length / 2.0) + (upstreamElement->length / 2.0));

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionDownstreamJunction(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = downstreamLongDispersion * downstreamXSectionArea *
                     (S[downstreamJunction->sIndex[soluteIndex]]  - S[sIndex[soluteIndex]]) /
      (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionDownstreamJunctionBC(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = downstreamLongDispersion * downstreamXSectionArea *
                     (downstreamJunction->soluteConcs[soluteIndex].value  - S[sIndex[soluteIndex]]) /
      (length / 2.0);

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionDownstreamNeighbour(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = downstreamLongDispersion * downstreamXSectionArea *
                     (S[downstreamElement->sIndex[soluteIndex]]  - S[sIndex[soluteIndex]]) /
      ((length / 2.0) +  (downstreamElement->length / 2.0));

  return DSoluteDt;
}

double Element::computeDSoluteDtDispersionSelf(double dt, double S[], int soluteIndex)
{
  return 0.0;
}

void Element::setDispersionFunctions()
{
  if(upstreamJunction->temperature.isBC)
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamJunctionBC;
  }
  else if(upstreamJunction->tIndex > -1)
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamJunction;
  }
  else if(upstreamElement)
  {
    computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamNeighbour;
  }
  else
  {
    if(flow.value  >= 0)
    {
      computeTempDispDeriv[0] = &Element::computeDTDtDispersionUpstreamJunctionBC;
    }
    else
    {
      computeTempDispDeriv[0] = &Element::computeDTDtDispersionSelf;
    }
  }

  if(downstreamJunction->temperature.isBC)
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamJunctionBC;
  }
  else if(downstreamJunction->tIndex > -1)
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamJunction;
  }
  else if(downstreamElement)
  {
    computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamNeighbour;
  }
  else
  {
    if(flow.value  >= 0)
    {
      computeTempDispDeriv[1] = &Element::computeDTDtDispersionSelf;
    }
    else
    {
      computeTempDispDeriv[1] = &Element::computeDTDtDispersionDownstreamJunctionBC;
    }
  }

  for(int i = 0; i < numSolutes; i++)
  {
    if(upstreamJunction->soluteConcs[i].isBC)
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamJunctionBC;
    }
    else if(upstreamJunction->sIndex[i] > -1)
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamJunction;
    }
    else if(upstreamElement)
    {
      computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamNeighbour;
    }
    else
    {
      if(flow.value  >= 0)
      {
        computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionUpstreamJunctionBC;
      }
      else
      {
        computeSoluteDispDeriv[i][0] = &Element::computeDSoluteDtDispersionSelf;
      }
    }

    if(downstreamJunction->soluteConcs[i].isBC)
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamJunctionBC;
    }
    else if(downstreamJunction->sIndex[i] > -1)
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamJunction;
    }
    else if(downstreamElement)
    {
      computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamNeighbour;
    }
    else
    {
      if(flow.value  >= 0)
      {
        computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionSelf;
      }
      else
      {
        computeSoluteDispDeriv[i][1] = &Element::computeDSoluteDtDispersionDownstreamJunctionBC;
      }
    }
  }
}

double Element::computeCourantFactor() const
{
  return  fabs(flow.value  / xSectionArea / length);
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

  rho_cp  = model->m_waterDensity * model->m_cp;
  rho_vol = model->m_waterDensity * volume;
  rho_cp_vol = rho_cp * volume;
  top_area = length * width;


  if(upstreamElement != nullptr)
  {
    upstreamFlow =  ((upstreamElement->flow.value  * upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                     (flow.value  / (length * 0.5))) /
                    ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamVelocity =  (((upstreamElement->flow.value  / upstreamElement->xSectionArea )* upstreamElementDirection  / (upstreamElement->length * 0.5)) +
                         ((flow.value  / xSectionArea) / (length * 0.5))) /
                        ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    upstreamXSectionArea = ((upstreamElement->xSectionArea / (upstreamElement->length * 0.5)) + (xSectionArea / (length * 0.5))) /
                           ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));


    upstreamCourantNumber = upstreamVelocity * model->m_timeStep / (length / 2.0 + upstreamElement->length / 2.0);

  }
  else
  {
    upstreamFlow = flow.value ;
    upstreamVelocity = upstreamFlow / xSectionArea;
    upstreamXSectionArea = xSectionArea;
    upstreamCourantNumber = upstreamVelocity * model->m_timeStep / (length);
  }

  if(downstreamElement != nullptr)
  {
    downstreamFlow = ((downstreamElement->flow.value  * downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                      (flow.value  / (length * 0.5))) /
                     ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamVelocity = (((downstreamElement->flow.value  / downstreamElement->xSectionArea )* downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                          ((flow.value  / xSectionArea) / (length * 0.5))) /
                         ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamXSectionArea = ((downstreamElement->xSectionArea / (downstreamElement->length * 0.5)) + (xSectionArea / (length * 0.5))) /
                             ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamCourantNumber = downstreamVelocity * model->m_timeStep / (length / 2.0 + downstreamElement->length / 2.0);

  }
  else
  {
    downstreamFlow = flow.value ;
    downstreamVelocity = downstreamFlow / xSectionArea;
    downstreamXSectionArea = xSectionArea;
    downstreamCourantNumber = downstreamVelocity * model->m_timeStep / (length);
  }

  //  setDispersionFunctions();
  (*setAdvectionFuctions)(this);
}

void Element::computeLongDispersion()
{
  double vel = flow.value  / xSectionArea;

  slope = max(0.00001, fabs(upstreamJunction->z - downstreamJunction->z) / length);
  double fricVel = sqrt(9.81 * depth * slope);
  double dispFischer = model->m_computeDispersion  * (0.011 * vel * vel * width * width) / (depth * fricVel);
  double dispNumerical = model->m_computeDispersion * fabs(vel * length / 2.0);

  longDispersion.value = dispNumerical <= dispFischer ? dispFischer - dispNumerical : dispNumerical;
//  longDispersion.value = dispFischer;
}

void Element::computePecletNumbers()
{
  pecletNumber = longDispersion.value ? (flow.value  / xSectionArea) * length / longDispersion.value : flow.value  * 10000 / fabs(flow.value );
}

void Element::computeUpstreamPeclet()
{
  if(upstreamElement != nullptr)
  {

    upstreamLongDispersion =  ((upstreamElement->longDispersion.value / (upstreamElement->length * 0.5)) + (longDispersion.value / (length * 0.5))) /
                              ((1.0 / (upstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    double dispNum =   model->m_computeDispersion * fabs(upstreamVelocity) * (length / 2.0 + upstreamElement->length / 2.0) / 2.0;

    upstreamLongDispersion = dispNum < upstreamLongDispersion ? upstreamLongDispersion - dispNum : dispNum;

    upstreamPecletNumber = upstreamLongDispersion > 0 ? upstreamVelocity * ((upstreamElement->length / 2.0) + (length / 2.0)) / (upstreamLongDispersion)
                                                      : upstreamVelocity * 10000 / upstreamVelocity;
  }
  else
  {
    upstreamLongDispersion = longDispersion.value;

    double dispNum =  model->m_computeDispersion * fabs(upstreamVelocity) * length / 2.0;

    upstreamLongDispersion = dispNum < upstreamLongDispersion ? upstreamLongDispersion - dispNum : dispNum;

    upstreamPecletNumber = pecletNumber;
  }

}

void Element::computeDownstreamPeclet()
{
  if(downstreamElement != nullptr)
  {
    downstreamLongDispersion = ((downstreamElement->longDispersion.value / (downstreamElement->length * 0.5)) + (longDispersion.value / (length * 0.5))) /
                               ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));


    double dispNum = model->m_computeDispersion *  fabs(downstreamVelocity) * (length / 2.0 + downstreamElement->length / 2.0) / 2.0;

    downstreamLongDispersion = dispNum < downstreamLongDispersion ? downstreamLongDispersion - dispNum : dispNum;


    downstreamPecletNumber = downstreamLongDispersion > 0 ? downstreamVelocity * ((downstreamElement->length / 2.0) + (length / 2.0)) / (downstreamLongDispersion)
                                                          : downstreamVelocity * 10000 / downstreamVelocity;
  }
  else
  {
    downstreamLongDispersion = longDispersion.value;

    double dispNum = model->m_computeDispersion * fabs(downstreamVelocity) * length  / 2.0;

    downstreamLongDispersion =  dispNum < downstreamLongDispersion ? downstreamLongDispersion - dispNum : dispNum;

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
    downstreamFlow = ((downstreamElement->flow.value  * downstreamElementDirection  / (downstreamElement->length * 0.5)) +
                      (flow.value  / (length * 0.5))) /
                     ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));

    downstreamVelocity = (((downstreamElement->flow.value  / downstreamElement->xSectionArea ) * downstreamElementDirection / (downstreamElement->length * 0.5)) +
                          ((flow.value  / xSectionArea) / (length * 0.5))) /
                         ((1.0 / (downstreamElement->length * 0.5)) + (1.0 / (length * 0.5)));
  }
  else
  {
    downstreamFlow = flow.value ;
    downstreamVelocity = downstreamFlow / xSectionArea;
  }
}

void Element::deleteSoluteVariables()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
    delete[] totalAdvDispSoluteMassBalance; totalAdvDispSoluteMassBalance = nullptr;
    delete[] totalExternalSoluteFluxesMassBalance; totalExternalSoluteFluxesMassBalance = nullptr;
    delete[] sIndex; sIndex = nullptr;

    for(int i = 0; i < numSolutes; i++)
    {
      delete[] computeSoluteAdvDeriv[i];
      delete[] computeSoluteDispDeriv[i];
    }

    delete[] computeSoluteAdvDeriv; computeSoluteAdvDeriv = nullptr;
    delete[] computeSoluteDispDeriv; computeSoluteDispDeriv = nullptr;
  }
}
