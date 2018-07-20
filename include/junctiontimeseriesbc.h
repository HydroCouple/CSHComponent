/*!
*  \file    junctiontimeseriesbc.h
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


#ifndef JUNCTIONTIMESERIESBC_H
#define JUNCTIONTIMESERIESBC_H

#include "abstracttimeseriesbc.h"


class STMCOMPONENT_EXPORT JunctionTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    JunctionTimeSeriesBC(ElementJunction *elementJunction, int variableIndex, STMModel *model);

    virtual ~JunctionTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    ElementJunction *elementJunction() const;

    void setElementJunction(ElementJunction *elementJunction);

  private:

    ElementJunction *m_elementJunction;
    int m_variableIndex;
};


#endif // JUNCTIONTIMESERIESBC_H
