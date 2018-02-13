#include "stdafx.h"
#include "elementjunction.h"
#include "element.h"

ElementJunction::ElementJunction(const std::string &id, double x, double y, double z, STMModel *model)
  :id(id), x(x), y(y), z(z),
   soluteContinuityIndexes(nullptr),
   numSolutes(0),
   soluteConcs(nullptr),
   prevSoluteConcs(nullptr),
   model(model)
{

}

ElementJunction::~ElementJunction()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
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

void ElementJunction::initializeSolutes(int numSolutes)
{
  if(numSolutes > 0 )
  {
    if(soluteConcs)
    {
      delete[] soluteConcs;
      delete[] prevSoluteConcs;
      delete[] externalSoluteFluxes;
      delete[] soluteContinuityIndexes;
    }

    this->numSolutes = numSolutes;
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new Variable[numSolutes];
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

double ElementJunction::computeDTDt(double dt)
{
  double volume = 0;
  double sumGradT = externalHeatFlux.value;

  for(Element *incomingElement : incomingElements)
  {
    volume += incomingElement->length * incomingElement->xSectionArea / 2.0;
    sumGradT += incomingElement->flow * incomingElement->temperature.value;
  }

  for(Element *outgoingElement : outgoingElements)
  {
    volume += outgoingElement->length * outgoingElement->xSectionArea / 2.0;
    sumGradT -= outgoingElement->flow * outgoingElement->temperature.value;
  }

  sumGradT /= volume;

  return sumGradT;
}

double ElementJunction::computeDSoluteDt(double dt, int soluteIndex)
{
  double volume = 0;
  double sumGradSolute = externalSoluteFluxes[soluteIndex].value;

  for(Element *incomingElement : incomingElements)
  {
    volume += incomingElement->length * incomingElement->xSectionArea / 2.0;
    sumGradSolute += incomingElement->flow * incomingElement->soluteConcs[soluteIndex].value;
  }


  for(Element *outgoingElement : outgoingElements)
  {
    volume += outgoingElement->length * outgoingElement->xSectionArea / 2.0;
    sumGradSolute -= outgoingElement->flow * outgoingElement->soluteConcs[soluteIndex].value;
  }

  sumGradSolute /= volume;

  return sumGradSolute;
}

double ElementJunction::computeDTDt(double dt, double T[])
{
  double volume = 0;
  double sumGradT = T[heatContinuityIndex];

  for(Element *incomingElement : incomingElements)
  {
    volume += incomingElement->length * incomingElement->xSectionArea / 2.0;
    sumGradT += incomingElement->flow * incomingElement->temperature.value;
  }


  for(Element *outgoingElement : outgoingElements)
  {
    volume += outgoingElement->length * outgoingElement->xSectionArea / 2.0;
    sumGradT -= outgoingElement->flow * outgoingElement->temperature.value;
  }

  sumGradT /= volume;

  return sumGradT;
}

double ElementJunction::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double volume = 0;
  double sumGradSolute = S[soluteIndex];

  for(Element *incomingElement : incomingElements)
  {
    volume += incomingElement->length * incomingElement->xSectionArea / 2.0;
    sumGradSolute += incomingElement->flow * incomingElement->soluteConcs[soluteIndex].value;
  }


  for(Element *outgoingElement : outgoingElements)
  {
    volume += outgoingElement->length * outgoingElement->xSectionArea / 2.0;
    sumGradSolute -= outgoingElement->flow * outgoingElement->soluteConcs[soluteIndex].value;
  }

  sumGradSolute /= volume;

  return sumGradSolute;
}

void ElementJunction::copyVariablesToPrev()
{
  prevTemperature.copy(temperature);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < numSolutes; i++)
  {
    prevSoluteConcs[i].copy(soluteConcs[i]);
  }
}
