/*!
*  \file    nonpointsrctimeseriesbc.cpp
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
*  \todo Test transport on branching networks
*  \warning
*/

#include "nonpointsrctimeseriesbc.h"
#include "element.h"
#include "cshmodel.h"

NonPointSrcTimeSeriesBC::NonPointSrcTimeSeriesBC(Element *startElement, double startElementLFactor,
                                                 Element *endElement, double endElementLFactor,
                                                 VariableType variableType, CSHModel *model)
  :AbstractTimeSeriesBC(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_startElementLFactor(startElementLFactor),
    m_endElementLFactor(endElementLFactor),
    m_variableType(variableType),
    m_soluteIndex(-1)
{

}

NonPointSrcTimeSeriesBC::~NonPointSrcTimeSeriesBC()
{

}

void NonPointSrcTimeSeriesBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
  m_factors.clear();

  for(Element *element : m_profile)
  {
    if(element != m_startElement && element != m_endElement)
    m_factors[element] = element->length;
  }

  m_factors[m_startElement] = m_startElement->length * m_startElementLFactor;
  m_factors[m_endElement]   = m_endElement->length * m_endElementLFactor;
}

void NonPointSrcTimeSeriesBC::prepare()
{

}

void NonPointSrcTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);


  if(found)
  {
    switch (m_variableType)
    {
      case HeatSource:
        for(Element *element : m_profile)
        {
          element->externalHeatFluxes += value * m_factors[element];
        }
        break;
      case FlowSource:
        for(Element *element : m_profile)
        {
          double factor = m_factors[element];

          element->externalHeatFluxes += m_model->m_cp * m_model->m_waterDensity * value *
                                         element->temperature.value * factor;

          for(size_t i = 0; i < m_model->m_solutes.size(); i++)
          {
            element->externalSoluteFluxes[i] += value *
                                                element->soluteConcs[i].value * factor;
          }

        }
        break;
      default:
        for(Element *element : m_profile)
        {
          element->externalSoluteFluxes[m_soluteIndex] += value * m_factors[element];
        }
        break;
    }
  }
}


Element *NonPointSrcTimeSeriesBC::startElement() const
{
  return m_startElement;
}

void NonPointSrcTimeSeriesBC::setStartElement(Element *element)
{
  m_startElement = element;
}


double NonPointSrcTimeSeriesBC::startElementLFactor() const
{
  return m_startElementLFactor;
}


void NonPointSrcTimeSeriesBC::setStartElementLFactor(double factor)
{
  m_startElementLFactor = factor;
}


Element *NonPointSrcTimeSeriesBC::endElement() const
{
  return m_endElement;
}


void NonPointSrcTimeSeriesBC::setEndElement(Element *element)
{
  m_endElement = element;
}


double NonPointSrcTimeSeriesBC::endElementLFactor() const
{
  return m_endElementLFactor;
}


void NonPointSrcTimeSeriesBC::setEndElementLFactor(double factor)
{
  m_endElementLFactor = factor;
}


int NonPointSrcTimeSeriesBC::soluteIndex() const
{
  return m_soluteIndex;
}

void NonPointSrcTimeSeriesBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}
