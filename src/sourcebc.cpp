/*!
*  \file    SourceBC.cpp
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
#include "stdafx.h"
#include "sourcebc.h"
#include "element.h"
#include "cshmodel.h"
#include "core/datacursor.h"
#include "temporal/timeseries.h"

SourceBC::SourceBC(Element *startElement, double startElementLFactor,
                   Element *endElement, double endElementLFactor,
                   VariableType variableType, CSHModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_startElementLFactor(startElementLFactor),
    m_endElementLFactor(endElementLFactor),
    m_variableType(variableType),
    m_soluteIndex(-1),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

SourceBC::~SourceBC()
{
  delete m_dataCursor;
}

void SourceBC::findAssociatedGeometries()
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

void SourceBC::prepare()
{

}

void SourceBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;

  if(m_timeSeries->numColumns() == (int)m_profile.size())
  {

    switch (m_variableType)
    {
      case HeatSource:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
          {
            Element *element = m_profile[i];
            element->externalHeatFluxes += value * m_factors[element];
          }
        }
        break;
      case FlowSource:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
          {
            Element *element = m_profile[i];
            double factor = m_factors[element];

            element->externalHeatFluxes += m_model->m_cp * m_model->m_waterDensity * value *
                                           element->temperature.value * factor;

            for(size_t i = 0; i < m_model->m_solutes.size(); i++)
            {
              element->externalSoluteFluxes[i] += value * element->soluteConcs[i].value * factor;
            }
          }
        }
        break;
      default:
        for(size_t i = 0; i < m_profile.size(); i++)
        {
          if(m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
          {
            Element *element = m_profile[i];
            element->externalSoluteFluxes[m_soluteIndex] += value * m_factors[element];
          }
        }
        break;
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      switch (m_variableType)
      {
        case HeatSource:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];
            element->externalHeatFluxes += value * m_factors[element];
          }
          break;
        case FlowSource:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];

            double factor = m_factors[element];

            element->externalHeatFluxes += m_model->m_cp * m_model->m_waterDensity * value *
                                           element->temperature.value * factor;

            for(size_t i = 0; i < m_model->m_solutes.size(); i++)
            {
              element->externalSoluteFluxes[i] += value * element->soluteConcs[i].value * factor;
            }
          }
          break;
        default:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            Element *element = m_profile[i];
            element->externalSoluteFluxes[m_soluteIndex] += value * m_factors[element];
          }
          break;
      }
    }
  }
}

void SourceBC::clear()
{

}

Element *SourceBC::startElement() const
{
  return m_startElement;
}

void SourceBC::setStartElement(Element *element)
{
  m_startElement = element;
}

double SourceBC::startElementLFactor() const
{
  return m_startElementLFactor;
}

void SourceBC::setStartElementLFactor(double factor)
{
  m_startElementLFactor = factor;
}

Element *SourceBC::endElement() const
{
  return m_endElement;
}

void SourceBC::setEndElement(Element *element)
{
  m_endElement = element;
}

double SourceBC::endElementLFactor() const
{
  return m_endElementLFactor;
}

void SourceBC::setEndElementLFactor(double factor)
{
  m_endElementLFactor = factor;
}

int SourceBC::soluteIndex() const
{
  return m_soluteIndex;
}

void SourceBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

QSharedPointer<TimeSeries> SourceBC::timeSeries() const
{
  return m_timeSeries;
}

void SourceBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}

