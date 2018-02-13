#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "stmmodel.h"

#include <math.h>

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  STMModel *model)
  : id(id),
    numSolutes(0),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    upstreamElement(nullptr),
    downstreamElement(nullptr),
    model(model)
{
  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);
}


Element::~Element()
{
  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initializeSolutes(int numSolutes)
{
  if(numSolutes > 0 )
  {
    if(soluteConcs)
    {
      delete[] soluteConcs;
      delete[] prevSoluteConcs;
      delete[] externalSoluteFluxes;
    }

    this->numSolutes = numSolutes;
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
  }
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
        upstreamElement = element;
        upstreamElementDirection = 1.0;
        return;
      }
    }

    for(Element *element : upstreamJunction->incomingElements)
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

double Element::computeDTDt(double dt)
{
  double DTdt = 0;
  double volume = xSectionArea * length;
  double rho_cp_vol = model->m_waterDensity * model->m_cp * volume;

  //Advection
  switch (model->m_advectionMode)
  {
    case STMModel::AdvectionDiscretizationMode::Central:
      {

      }
      break;
    case STMModel::AdvectionDiscretizationMode::Hybrid:
      {

      }
      break;
    default:
      {
        //Right Direction
        if(flow >= 0)
        {
          if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
          {

            DTdt += (upstreamFlow * upstreamElement->temperature.value - flow * temperature.value) / volume;
          }
          else
          {
            DTdt += (flow * upstreamJunction->temperature.value - flow * temperature.value) / volume;
          }
        }
        //Wrong Direction
        else
        {
          if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
          {

            DTdt += (flow * temperature.value - downstreamFlow *  downstreamElement->temperature.value) / volume;
          }
          else
          {
            DTdt += ( flow * temperature.value - flow * downstreamElement->temperature.value) / volume;
          }
        }
      }
      break;
  }

  //Dispersion
  {
    if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
    {

      DTdt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
              (upstreamElement->temperature.value - temperature.value) / (volume * (length + upstreamElement->length) / 2.0);
    }
    else
    {
      DTdt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
              (upstreamJunction->temperature.value - temperature.value) / (volume * length / 2.0);
    }

    if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
    {
      DTdt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
              (temperature.value - downstreamElement->temperature.value) / (volume * (length + downstreamElement->length) / 2.0);
    }
    else
    {
      DTdt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
              (temperature.value - downstreamJunction->temperature.value) / (volume * length / 2.0);
    }
  }

  //Add sources
  {
    DTdt += radiationFluxes / depth / rho_cp_vol;
    DTdt += externalHeatFluxes / rho_cp_vol;
  }

  return DTdt;
}

double Element::computeDSoluteDt(double dt, int soluteIndex)
{

}


double Element::computeDTDt(double dt, double T[])
{
  double DTdt = 0;
  double volume = xSectionArea * length;
  double rho_cp_vol = model->m_waterDensity * model->m_cp * volume;

  //Advection
  switch (model->m_advectionMode)
  {
    case STMModel::AdvectionDiscretizationMode::Central:
      {

      }
      break;
    case STMModel::AdvectionDiscretizationMode::Hybrid:
      {

      }
      break;
    default:
      {
        //Right Direction
        if(flow >= 0)
        {
          if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
          {

            DTdt += (upstreamFlow * T[upstreamElement->index] - flow * T[index]) / volume;
          }
          else
          {
            DTdt += (flow * T[upstreamJunction->index] - flow * T[index]) / volume;
          }
        }
        //Wrong Direction
        else
        {
          if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
          {

            DTdt += (flow * T[index] - downstreamFlow *  T[downstreamElement->index]) / volume;
          }
          else
          {
            DTdt += ( flow * T[index] - flow * T[downstreamElement->index]) / volume;
          }
        }
      }
      break;
  }

  //Dispersion
  {
    if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
    {

      DTdt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
              (T[upstreamElement->index] - T[index]) / (volume * (length + upstreamElement->length) / 2.0);
    }
    else
    {
      DTdt += upstreamJunction->longDispersion * upstreamJunction->xSectionArea *
              (upstreamJunction->temperature.value - T[index]) / (volume * length / 2.0);
    }

    if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
    {
      DTdt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
              (T[index] - T[downstreamElement->index]) / (volume * (length + downstreamElement->length) / 2.0);
    }
    else
    {
      DTdt -= downstreamJunction->longDispersion * downstreamJunction->xSectionArea *
              (T[index] - downstreamJunction->temperature.value) / (volume * length / 2.0);
    }
  }

  //Add sources
  {
    DTdt += radiationFluxes / depth / rho_cp_vol;
    DTdt += externalHeatFluxes / rho_cp_vol;
  }

  return DTdt;
}


double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{

}

double Element::computeCourantFactor() const
{
  return length / flow / xSectionArea;
}

void Element::computeLongDispersion()
{
  double vel = flow / xSectionArea;

  double fricVel = sqrt(9.81 * depth * slope);
  double dispFischer = 0.011 * vel * vel * width * width / (depth * fricVel);

  double dispNumerical = fabs( vel * length / 2.0);

  longDispersion.value = dispNumerical <= dispFischer ? dispFischer - dispNumerical : 0.0;
}

void Element::computePecletNumbers()
{
  pecletNumber = (flow / xSectionArea) / (longDispersion.value / length);
}

void Element::computeUpstreamPeclet()
{
  computeUpstreamFlow();
  upstreamPecletNumber = upstreamVelocity / (upstreamJunction->longDispersion_length);
}

void Element::computeDownstreamPeclet()
{
  computeDownstreamFlow();
  downstreamPecletNumber = downstreamVelocity / (downstreamJunction->longDispersion_length);
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
