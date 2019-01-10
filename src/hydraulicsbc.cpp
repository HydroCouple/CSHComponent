#include "stdafx.h"
#include "element.h"
#include "cshmodel.h"
#include "hydraulicsbc.h"
#include "temporal/timeseries.h"
#include "core/datacursor.h"

HydraulicsBC::HydraulicsBC(Element *startElement, Element *endElement, int variableIndex, CSHModel *model)
  : QObject(model),
    m_startElement(startElement),
    m_endElement(endElement),
    m_variableIndex(variableIndex),
    m_model(model)
{
  m_dataCursor = new DataCursor();
}

HydraulicsBC::~HydraulicsBC()
{
  delete m_dataCursor;
}

void HydraulicsBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_model->findProfile(m_startElement, m_endElement, m_profile);
}

void HydraulicsBC::prepare()
{

}

void HydraulicsBC::applyBoundaryConditions(double dateTime)
{

  double value = 0;

  if(m_timeSeries->numColumns() == (int)m_profile.size())
  {
    for(size_t i = 0; i < m_profile.size(); i++)
    {
      if( m_timeSeries->interpolate(dateTime, i, m_dataCursor, value))
      {
        switch (m_variableIndex)
        {
          case 1:
            m_profile[i]->depth = value;
            break;
          case 2:
            m_profile[i]->width = value;
            break;
          case 3:
            m_profile[i]->xSectionArea = value;
            break;
          case 4:
            m_profile[i]->flow = value;
            break;
        }
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      switch (m_variableIndex)
      {
        case 1:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            m_profile[i]->depth = value;
          }
          break;
        case 2:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            m_profile[i]->width = value;
          }
          break;
        case 3:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            m_profile[i]->xSectionArea = value;
          }
          break;
        case 4:
          for(size_t i = 0; i < m_profile.size(); i++)
          {
            m_profile[i]->flow = value;
          }
          break;
      }
    }
  }
}

void HydraulicsBC::clear()
{

}

Element *HydraulicsBC::startElement() const
{
  return m_startElement;
}

void HydraulicsBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *HydraulicsBC::endElement() const
{
  return m_endElement;
}

void HydraulicsBC::setEndElement(Element *element)
{
  m_endElement = element;
}

QSharedPointer<TimeSeries> HydraulicsBC::timeSeries() const
{
  return m_timeSeries;
}

void HydraulicsBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}

