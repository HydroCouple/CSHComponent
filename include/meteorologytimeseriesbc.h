/*!
*  \file    meteorologytimeseriesbc.h
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

#ifndef METEOROLOGYTIMESERIES_H
#define METEOROLOGYTIMESERIES_H

#include "abstracttimeseriesbc.h"

class STMCOMPONENT_EXPORT MeteorologyTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    MeteorologyTimeSeriesBC(Element *element, int variableIndex, STMModel *model);

    virtual ~MeteorologyTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

  private:

    Element *m_element;
    int m_variableIndex;
};


class STMCOMPONENT_EXPORT UniformMeteorologyTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    UniformMeteorologyTimeSeriesBC(Element *startElement, Element *endElement, int variableIndex, STMModel *model);

    virtual ~UniformMeteorologyTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

  private:
    std::list<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_variableIndex;
};

#endif // METEOROLOGYTIMESERIES_H
