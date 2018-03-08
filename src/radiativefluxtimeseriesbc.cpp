/*!
*  \file    radiativefluxtimeseriesbc.cpp
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
*  \todo Test transport on branching networks
*  \warning
*/


#include "radiativefluxtimeseriesbc.h"
#include "stmmodel.h"
#include "element.h"

RadiativeFluxTimeSeriesBC::RadiativeFluxTimeSeriesBC(Element *startElement, Element *endElement, STMModel *model)
  :AbstractTimeSeriesBC(model),
   m_startElement(startElement),
   m_endElement(endElement)
{

}

RadiativeFluxTimeSeriesBC::~RadiativeFluxTimeSeriesBC()
{

}

void RadiativeFluxTimeSeriesBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void RadiativeFluxTimeSeriesBC::prepare()
{

}

void RadiativeFluxTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    for(Element *element : m_profile)
    {
      element->radiationFluxes += value;
    }
  }
}
