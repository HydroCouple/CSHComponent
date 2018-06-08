/*!
*  \file    pointsrctimeseriesbc.h
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

#ifndef POINTSRCTIMESERIESBC_H
#define POINTSRCTIMESERIESBC_H

#include "abstracttimeseriesbc.h"

class STMCOMPONENT_EXPORT PointSrcTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    enum VariableType
    {
      HeatSource,
      FlowSource,
      SoluteSource,
    };

    PointSrcTimeSeriesBC(Element *element, VariableType variableType, STMModel *model);

    virtual ~PointSrcTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *element() const;

    void setElement(Element *element);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:
    Element *m_element;
    VariableType m_variableType;
    int m_soluteIndex;
};


#endif // POINTSRCTIMESERIESBC_H
