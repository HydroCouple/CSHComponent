/*!
*  \file    elementjunction.cpp
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
#include "elementjunction.h"
#include "element.h"
#include "cshmodel.h"

#include <math.h>

ElementJunction::ElementJunction(const std::string &id, double x, double y, double z, CSHModel *model)
  :id(id), x(x), y(y), z(z),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    model(model)
{
  hIndex = -1;
  tIndex = -1;
  initializeSolutes();
}

ElementJunction::~ElementJunction()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
  }

  while (outgoingElements.size())
  {
    Element *element = *outgoingElements.begin();
    delete element;
  }

  while (incomingElements.size())
  {
    Element *element = *incomingElements.begin();
    delete element;
  }
}

void ElementJunction::initialize()
{
  starting = true;
  volume = prev_volume = 0.0;
}

void ElementJunction::interpTemp()
{
  //IDW interpolation for
  double sum_x = 0;
  double tempTemp = 0;

  if(incomingElements.size() > 0)
  {
    for(Element *element : this->incomingElements)
    {
      tempTemp += element->temperature.value / element->length / 2.0;
      sum_x += 1.0 / element->length / 2.0;
    }
  }
  else if(outgoingElements.size() > 0)
  {
    for(Element *element : this->outgoingElements)
    {
      tempTemp += element->temperature.value / element->length / 2.0;
      sum_x += 1.0 / element->length / 2.0;
    }
  }
  else
  {
    sum_x = 1;
  }


  this->temperature.value = tempTemp / sum_x;
}

void ElementJunction::interpSoluteConcs(int soluteIndex)
{
  //IDW interpolation for
  double sum_S_x = 0;
  double sum_x = 0;

  if(incomingElements.size() > 0)
  {
    for(Element *element : this->incomingElements)
    {
      sum_S_x += element->soluteConcs[soluteIndex].value / element->length / 2.0;
      sum_x += 1.0 / element->length / 2.0;
    }
  }
  else if(outgoingElements.size() > 0)
  {
    for(Element *element : this->outgoingElements)
    {
      sum_S_x += element->soluteConcs[soluteIndex].value / element->length / 2.0;
      sum_x += 1.0 / element->length / 2.0;
    }
  }
  else
  {
    sum_x = 1;
  }

  this->soluteConcs[soluteIndex].value = sum_S_x / sum_x;
}

void ElementJunction::solveHeatContinuity(double dt)
{

  double sumQ = 0.0;
  double sumQT = 0.0;

  for(Element *incomingElement : incomingElements)
  {
    double q = std::max(0.0 , incomingElement->flow.value );
    sumQ += q;
    sumQT += q * incomingElement->temperature.value;
  }

  for(Element *outgoingElement : outgoingElements)
  {
    double q = fabs(std::min(0.0 , outgoingElement->flow.value ));
    sumQ += q;
    sumQT += q * outgoingElement->temperature.value;
  }

  double temp = sumQ ? sumQT / sumQ : 0.0;
  temperature.value = temp;
}

void ElementJunction::solveSoluteContinuity(int soluteIndex, double dt)
{
  double sumQ = 0.0;
  double sumQT = 0.0;

  for(Element *incomingElement : incomingElements)
  {
    double q = std::max(0.0 , incomingElement->flow.value );
    sumQ += q;
    sumQT += q * incomingElement->soluteConcs[soluteIndex].value;
  }

  for(Element *outgoingElement : outgoingElements)
  {
    double q = fabs(std::min(0.0 , outgoingElement->flow.value ));
    sumQ += q;
    sumQT += q * outgoingElement->soluteConcs[soluteIndex].value;
  }

  double temp = sumQ ? sumQT / sumQ : 0.0;
  soluteConcs[soluteIndex].value = temp;
}

double ElementJunction::computeDTDt(double dt, double T[])
{
  double DTDt = 0.0;

  for(Element *element : this->incomingElements)
  {
    DTDt += element->rho_cp * element->flow.value  * T[element->tIndex] / (element->rho_cp * volume);

    DTDt += (element->longDispersion.value * element->xSectionArea * element->rho_cp *
            (T[element->tIndex] - T[tIndex]) /
            (element->length / 2.0)) / (element->rho_cp * volume);
  }

  for(Element *element : this->outgoingElements)
  {
    DTDt -= element->rho_cp * element->flow.value  * T[element->tIndex] / (element->rho_cp * volume);

    DTDt += (element->longDispersion.value * element->xSectionArea * element->rho_cp *
            (T[element->tIndex] - T[tIndex]) /
            (element->length / 2.0))/ (element->rho_cp * volume);
  }

  DTDt -= temperature.value * dvolume_dt / volume;

  return DTDt;
}

double ElementJunction::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(Element *element : this->incomingElements)
  {
    DSoluteDt += element->flow.value  * S[element->sIndex[soluteIndex]] / ( volume);

    DSoluteDt += (element->longDispersion.value * element->xSectionArea *
            (S[element->sIndex[soluteIndex]] - S[sIndex[soluteIndex]]) /
            (element->length / 2.0)) / volume;
  }

  for(Element *element : this->outgoingElements)
  {
    DSoluteDt -= element->flow.value  * S[element->sIndex[soluteIndex]] / ( volume);

    DSoluteDt += (element->longDispersion.value * element->xSectionArea *
            (S[element->sIndex[soluteIndex]] - S[sIndex[soluteIndex]]) /
            (element->length / 2.0)) / volume;
  }

  DSoluteDt -= soluteConcs[soluteIndex].value * dvolume_dt / volume;

  return DSoluteDt;
}

void ElementJunction::copyVariablesToPrev()
{
  prevTemperature.copy(temperature);

  for(int i = 0 ; i < numSolutes; i++)
  {
    prevSoluteConcs[i].copy(soluteConcs[i]);
  }
}

void ElementJunction::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] sIndex; sIndex = nullptr;
  }

  if(model->m_solutes.size() > 0 )
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    sIndex = new int[numSolutes]();

    for(int i = 0 ; i < numSolutes; i++)
    {
      sIndex[i] = -1;
    }
  }
}

void ElementJunction::computeDerivedHydraulics()
{

  if(starting)
  {
    prev_volume = volume = 0;

    for(Element *element : this->outgoingElements)
    {
      volume += element->xSectionArea * element->length / 2.0;
    }

    for(Element *element : this->incomingElements)
    {
      volume += element->xSectionArea * element->length / 2.0;
    }

    prev_volume = volume;

    dvolume_dt = 0.0;

    starting = false;
  }
  else
  {

    prev_volume = volume;
    volume = 0.0;

    for(Element *element : this->outgoingElements)
    {
      volume += element->xSectionArea * element->length / 2.0;
    }

    for(Element *element : this->incomingElements)
    {
      volume += element->xSectionArea * element->length / 2.0;
    }

    dvolume_dt = (volume - prev_volume) / model->m_prevTimeStep;
  }

  if(!inflow.isBC)
  {
    inflow.value = 0;

    for(Element *element : this->incomingElements)
    {
      inflow.value += element->flow.value;
    }
  }
}

void ElementJunction::computeInflow(double Q[])
{
  if(!inflow.isBC)
  {
    inflow.value = 0;

    for(Element *element : this->incomingElements)
    {
      inflow.value += Q[element->hIndex];
    }
  }
}
