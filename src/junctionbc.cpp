/*!
*  \file    JunctionBC.cpp
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
#include "junctionbc.h"
#include "elementjunction.h"
#include "temporal/timeseries.h"
#include "core/datacursor.h"
#include "cshmodel.h"

JunctionBC::JunctionBC(ElementJunction *elementJunction, int variableIndex, CSHModel *model)
  : QObject(model),
    m_elementJunction(elementJunction),
    m_variableIndex(variableIndex),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

JunctionBC::~JunctionBC()
{
  delete  m_dataCursor;
}

void JunctionBC::findAssociatedGeometries()
{

}

void JunctionBC::prepare()
{
  switch (m_variableIndex)
  {
    case -1:
      {
        m_elementJunction->temperature.isBC = true;
      }
      break;
    default:
      {
        m_elementJunction->soluteConcs[m_variableIndex].isBC = true;
      }
      break;
  }
}

void JunctionBC::applyBoundaryConditions(double dateTime)
{
  double value = 0;

  if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
  {
    switch (m_variableIndex)
    {
      case -1:
        {
          m_elementJunction->temperature.value = value;
        }
        break;
      default:
        {
          m_elementJunction->soluteConcs[m_variableIndex].value = value;
        }
        break;
    }
  }
}

void JunctionBC::clear()
{

}

ElementJunction *JunctionBC::elementJunction() const
{
  return m_elementJunction;
}

void JunctionBC::setElementJunction(ElementJunction *elementJunction)
{
  m_elementJunction = elementJunction;
}

QSharedPointer<TimeSeries> JunctionBC::timeSeries() const
{
  return m_timeSeries;
}

void JunctionBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}
