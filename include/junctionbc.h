/*!
*  \file    JunctionBC.h
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


#ifndef JUNCTIONBC_H
#define JUNCTIONBC_H

#include "iboundarycondition.h"
#include "cshcomponent_global.h"

#include <QObject>
#include <QSharedPointer>

struct ElementJunction;
class DataCursor;
class CSHModel;
class TimeSeries;

class CSHCOMPONENT_EXPORT JunctionBC: public QObject,
    public virtual IBoundaryCondition
{
    Q_OBJECT

  public:

    JunctionBC(ElementJunction *elementJunction, int variableIndex, CSHModel *model);

    virtual ~JunctionBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override;

    ElementJunction *elementJunction() const;

    void setElementJunction(ElementJunction *elementJunction);

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:

    ElementJunction *m_elementJunction;
    int m_variableIndex;
    DataCursor *m_dataCursor;
    QSharedPointer<TimeSeries> m_timeSeries;
    CSHModel *m_model;
};


#endif // JUNCTIONBC_H
