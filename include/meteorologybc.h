/*!
*  \file    meteorologytimeseriesbc.h
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
*  \todo
*  \warning
*/

#ifndef METEOROLOGYTIMESERIES_H
#define METEOROLOGYTIMESERIES_H

#include "iboundarycondition.h"
#include "cshcomponent_global.h"

#include <QObject>
#include <QSharedPointer>

struct Element;
class CSHModel;
class DataCursor;
class TimeSeries;


class CSHCOMPONENT_EXPORT MeteorologyBC: public QObject,
    public virtual IBoundaryCondition
{
    Q_OBJECT

  public:

    MeteorologyBC(Element *startElement,
                            Element *endElement,
                            int variableIndex,
                            CSHModel *model);

    virtual ~MeteorologyBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:

    std::vector<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    int m_variableIndex;
    DataCursor *m_dataCursor;
    QSharedPointer<TimeSeries> m_timeSeries;
    CSHModel *m_model;
};

#endif // METEOROLOGYTIMESERIES_H
