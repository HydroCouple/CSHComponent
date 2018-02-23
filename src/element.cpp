#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "stmmodel.h"

#include <math.h>

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  STMModel *model)
  : id(id),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    externalSoluteFluxes(nullptr),
    upstreamElement(nullptr),
    downstreamElement(nullptr),
    model(model),
    computeTempAdv(nullptr),
    computeSoluteAdv(nullptr)
{
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
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);

}

void Element::initialize()
{
  if(model->m_solutes.size() > 0 )
  {
    if(soluteConcs)
    {
      delete[] soluteConcs;
      delete[] prevSoluteConcs;
      delete[] externalSoluteFluxes;
    }

    numSolutes = model->m_solutes.size();
    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes];
  }

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
    default:
      {
        computeTempAdv = &Element::computeDTDtUpwind;
        computeSoluteAdv = &Element::computeDSoluteDtUpwind;
      }
      break;
  }
}

double Element::computeDTDt(double dt, double T[])
{
  double DTDt = 0;
  double volume = xSectionArea * length;
  double rho_cp_vol = model->m_waterDensity * model->m_cp * volume;

  //Compute advection
  DTDt += (this->*computeTempAdv)(dt, T);

  //Compute dispersion
  DTDt += computeDTDtDispersion(dt, T);

  //Add external sources
  {
    DTDt += radiationFluxes / depth / rho_cp_vol;
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

    incomingFlux = upstreamFlow * (((T[upstreamElement->index] / upstreamElement->length / 2.0) +
                                   (T[index]/ length / 2.0)) / ((1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0))) / volume;
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->temperature.value / volume;
  }

  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    outgoingFlux = downstreamFlow * (((T[downstreamElement->index] / downstreamElement->length / 2.0) +
                                     (T[index]/ length / 2.0)) / ((1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0))) / volume;
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

  if(upstreamElement != nullptr && !upstreamJunction->temperature.isBC)
  {
    //If peclet number is zero use central differencing
    if(fabs(upstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
    {
      incomingFlux = upstreamFlow * (((T[upstreamElement->index] / upstreamElement->length / 2.0) +
                                     (T[index]/ length / 2.0)) / ((1.0 / upstreamElement->length / 2.0) + (1.0 / length/ 2.0))) / volume;
    }
    else if(upstreamPecletNumber >= 2.0)
    {
      incomingFlux = upstreamFlow * T[upstreamElement->index] / volume;
    }
    else if(upstreamPecletNumber <= -2.0)
    {
      incomingFlux = flow * T[index] / volume;
    }
    else
    {
      double upstreamFactor = 1.0 / upstreamElement->length / 2.0;
      double centerFactor = 1.0 / length / 2.0;

      double idwDenomFactor = upstreamFactor + centerFactor;

      upstreamFactor = upstreamFactor / idwDenomFactor;
      centerFactor   = centerFactor / idwDenomFactor;

      incomingFlux =  upstreamFlow * (((1 + (1.0 / upstreamPecletNumber / upstreamFactor)) * T[upstreamElement->index] * upstreamFactor) +
                      (1 - (1.0 / upstreamPecletNumber / centerFactor)) * T[index] * centerFactor) / volume;
    }
  }
  else
  {
    incomingFlux = upstreamFlow * upstreamJunction->temperature.value / volume;
  }

  if(downstreamElement != nullptr && !downstreamJunction->temperature.isBC)
  {
    if(fabs(downstreamPecletNumber - 0.0) > std::numeric_limits<double>::epsilon() == 0)
    {
      outgoingFlux = downstreamFlow * (((T[downstreamElement->index] / downstreamElement->length / 2.0) +
                                       (T[index]/ length / 2.0)) / ((1.0 / downstreamElement->length / 2.0) + (1.0 / length/ 2.0))) / volume;
    }
    else if(downstreamPecletNumber >= 2.0)
    {
      outgoingFlux = flow * T[index] / volume;
    }
    else if(downstreamPecletNumber <= -2.0)
    {
      outgoingFlux = downstreamFlow * T[downstreamElement->index] / volume;
    }
    else
    {
      double downstreamFactor = 1.0 / downstreamElement->length / 2.0;
      double centerFactor = 1.0 / length / 2.0;

      double idwDenomFactor = centerFactor + downstreamFactor;

      centerFactor /= idwDenomFactor;
      downstreamFactor   /= idwDenomFactor;

      outgoingFlux =  downstreamFlow * (((1 + (1.0 / downstreamPecletNumber/ centerFactor)) *  T[index] * centerFactor ) +
                                        (1 - (1.0 / downstreamPecletNumber / downstreamFactor)) * T[downstreamElement->index] * downstreamFactor  ) / volume;
    }
  }
  else
  {
    outgoingFlux = downstreamFlow * downstreamJunction->temperature.value / volume;
  }


  return incomingFlux - outgoingFlux;
}

double Element::computeDSoluteDt(double dt, double S[], int soluteIndex)
{

}

double Element::computeDSoluteDtDispersion(double dt, double t[])
{

}

double Element::computeDSoluteDtUpwind(double dt, double S[], int soluteIndex)
{

}

double Element::computeDSoluteDtCentral(double dt, double S[], int soluteIndex)
{

}

double Element::computeDSoluteDtHybrid(double dt, double S[], int soluteIndex)
{

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
  double dispFischer = 0.011 * vel * vel * width * width / (depth * fricVel);

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
