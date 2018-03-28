/*!
*  \file    nonpointsrctimeseriesbc.h
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

#ifndef NONPOINTSRCTIMESERIESBC_H
#define NONPOINTSRCTIMESERIESBC_H

#include "abstracttimeseriesbc.h"

class STMCOMPONENT_EXPORT NonPointSrcTimeSeriesBC : public AbstractTimeSeriesBC
{
  public:

    NonPointSrcTimeSeriesBC(Element *startElement, double startElementLFactor,
                            Element *endElement, double endElementLFactor,
                            int variableIndex, STMModel *model);

    virtual ~NonPointSrcTimeSeriesBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    Element *startElement() const;

    void setStartElement(Element *element);

    double startElementLFactor() const;

    void setStartElementLFactor(double factor);

    Element *endElement() const;

    void setEndElement(Element *element);

    double endElementLFactor() const;

    void setEndElementLFactor(double factor);

  private:
    std::list<Element*> m_profile;
    Element *m_startElement, *m_endElement;
    double m_startElementLFactor, m_endElementLFactor;
    int m_variableIndex;
};


#endif // NONPOINTSRCTIMESERIESBC_H
