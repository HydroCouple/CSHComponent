/*!
 * \file   elementoutput.h
 * \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 * \version   1.0.0
 * \description
 * \license
 * This file and its associated files, and libraries are free software.
 * You can redistribute it and/or modify it under the terms of the
 * Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 3 of the License, or (at your option) any later version.
 * This file and its associated files is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 * \copyright Copyright 2014-2018, Caleb Buahin, All rights reserved.
 * \date 2014-2018
 * \pre
 * \bug
 * \warning
 * \todo
 */

#ifndef ELEMENTOUTPUT_H
#define ELEMENTOUTPUT_H

#include "stmcomponent_global.h"
#include "spatiotemporal/timegeometryoutput.h"

class STMComponent;

class STMCOMPONENT_EXPORT ElementOutput: public TimeGeometryOutputDouble
{

    Q_OBJECT

  public:

    enum VariableType
    {
      Flow,
      XSectionArea,
      Depth,
      TopWidth,
      Temperature,
      SoluteConc
    };

    ElementOutput(const QString &id,
                  Dimension *timeDimension,
                  Dimension *geometryDimension,
                  ValueDefinition *valueDefinition,
                  VariableType variableType,
                  STMComponent *modelComponent);


    virtual ~ElementOutput();

    void updateValues(HydroCouple::IInput *querySpecifier) override;

    void updateValues() override;

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    STMComponent *m_component;
    int m_variableType;
    int m_soluteIndex;
};


#endif // ELEMENTOUTPUT_H
