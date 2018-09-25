/*!
*  \file    pointsrctimeseriesbc.cpp
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

#include "pointsrctimeseriesbc.h"
#include "element.h"
#include "cshmodel.h"

PointSrcTimeSeriesBC::PointSrcTimeSeriesBC(Element *element, VariableType variableType, CSHModel *model)
  : AbstractTimeSeriesBC(model),
    m_element(element),
    m_variableType(variableType),
    m_soluteIndex(-1)
{

}

PointSrcTimeSeriesBC::~PointSrcTimeSeriesBC()
{

}

void PointSrcTimeSeriesBC::findAssociatedGeometries()
{

}

void PointSrcTimeSeriesBC::prepare()
{

}

void PointSrcTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableType)
    {
      case HeatSource:
        {
          m_element->externalHeatFluxes += value;
        }
        break;
      case FlowSource:
        {
          m_element->externalHeatFluxes +=  m_model->m_cp * m_model->m_waterDensity * value * m_element->temperature.value;

          for(size_t i = 0; i < m_model->m_solutes.size(); i++)
          {
            m_element->externalSoluteFluxes[i] += m_model->m_waterDensity * value * m_element->soluteConcs[i].value;
          }
        }
        break;
      default:
        {
          m_element->externalSoluteFluxes[m_soluteIndex] += value;
        }
        break;
    }
  }
}

Element *PointSrcTimeSeriesBC::element() const
{
  return m_element;
}

void PointSrcTimeSeriesBC::setElement(Element *element)
{
  m_element = element;
}

int PointSrcTimeSeriesBC::soluteIndex() const
{
  return m_soluteIndex;
}

void PointSrcTimeSeriesBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

