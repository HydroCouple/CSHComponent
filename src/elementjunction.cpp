/*!
*  \file    elementjunction.cpp
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
#include "elementjunction.h"
#include "element.h"
#include "stmmodel.h"

#include <math.h>

ElementJunction::ElementJunction(const std::string &id, double x, double y, double z, STMModel *model)
  :id(id), x(x), y(y), z(z),
    heatContinuityIndex(-1),
    soluteContinuityIndexes(nullptr),
    numSolutes(0),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    longDispersion(0.0),
    model(model)
{
  initializeSolutes();
}

ElementJunction::~ElementJunction()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] soluteContinuityIndexes;
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

void ElementJunction::initializeSolutes()
{

  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] soluteContinuityIndexes; soluteContinuityIndexes = nullptr;
  }

  if(model->m_solutes.size() > 0 )
  {
    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    soluteContinuityIndexes = new int[numSolutes];

    for(int i = 0 ; i < numSolutes; i++)
    {
      soluteContinuityIndexes[i] = -1;
    }
  }
}

void ElementJunction::interpLongDispersion()
{
  //IDW interpolation for
  this->longDispersion = 0;
  this->longDispersion_length = 0;
  double sum_x = 0;

  for(Element *element : this->outgoingElements)
  {
    longDispersion += element->longDispersion.value / element->length / 2.0;
    longDispersion_length += element->longDispersion.value / element->length / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  for(Element *element : this->incomingElements)
  {
    longDispersion += element->longDispersion.value / element->length / 2.0;
    longDispersion_length += element->longDispersion.value / element->length / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  this->longDispersion /= sum_x;
  this->longDispersion_length /= sum_x;
}

void ElementJunction::interpTemp()
{
  //IDW interpolation for
  this->temperature.value = 0;
  double sum_x = 0;

  for(Element *element : this->outgoingElements)
  {
    this->temperature.value += element->temperature.value / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  for(Element *element : this->incomingElements)
  {
    this->temperature.value += element->temperature.value / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  this->temperature.value /= sum_x;
}

void ElementJunction::interpXSectionArea()
{
  //IDW interpolation for
  double sum_Area_x = 0;
  double sum_x = 0;

  for(Element *element : this->outgoingElements)
  {
    sum_Area_x += element->xSectionArea / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  for(Element *element : this->incomingElements)
  {
    sum_Area_x += element->xSectionArea / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  xSectionArea = sum_Area_x / sum_x;
}

void ElementJunction::interpSoluteConcs(int soluteIndex)
{
  //IDW interpolation for
  double sum_S_x = 0;
  double sum_x = 0;

  for(Element *element : this->outgoingElements)
  {
    sum_S_x += element->soluteConcs[soluteIndex].value / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  for(Element *element : this->incomingElements)
  {
    sum_S_x += element->soluteConcs[soluteIndex].value / element->length / 2.0;
    sum_x += 1.0 / element->length / 2.0;
  }

  this->soluteConcs[soluteIndex].value = sum_S_x / sum_x;
}

void ElementJunction::solveHeatContinuity(double dt)
{
  double sumQ = 0.0;
  double sumQT = 0.0;

  for(Element *incomingElement : incomingElements)
  {
    double q = std::max(0.0 , incomingElement->flow);
    sumQ += q;
    sumQT += q * incomingElement->temperature.value;
  }

  for(Element *outgoingElement : outgoingElements)
  {
    double q = fabs(std::min(0.0 , outgoingElement->flow));
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
    double q = std::max(0.0 , incomingElement->flow);
    sumQ += q;
    sumQT += q * incomingElement->soluteConcs[soluteIndex].value;
  }

  for(Element *outgoingElement : outgoingElements)
  {
    double q = fabs(std::min(0.0 , outgoingElement->flow));
    sumQ += q;
    sumQT += q * outgoingElement->soluteConcs[soluteIndex].value;
  }

  double temp = sumQ ? sumQT / sumQ : 0.0;
  soluteConcs[soluteIndex].value = temp;
}

void ElementJunction::copyVariablesToPrev()
{
  prevTemperature.copy(temperature);

  for(int i = 0 ; i < numSolutes; i++)
  {
    prevSoluteConcs[i].copy(soluteConcs[i]);
  }
}
