
#include "hydraulicstimeseriesbc.h"
#include "element.h"


HydraulicsTimeSeriesBC::HydraulicsTimeSeriesBC(Element *element, int variableIndex, STMModel *model)
  :AbstractTimeSeriesBC(model),
   m_element(element),
   m_variableIndex(variableIndex)
{

}

HydraulicsTimeSeriesBC::~HydraulicsTimeSeriesBC()
{

}

void HydraulicsTimeSeriesBC::findAssociatedGeometries()
{

}

void HydraulicsTimeSeriesBC::prepare()
{

}

void HydraulicsTimeSeriesBC::applyBoundaryConditions(double dateTime)
{
  bool found ;
  double value = interpolate(dateTime, found);

  if(found)
  {
    switch (m_variableIndex)
    {
      case 1:
        m_element->depth = value;
        break;
      case 2:
        m_element->width = value;
        break;
      case 3:
        m_element->xSectionArea = value;
        break;
      case 4:
        m_element->flow = value;
        break;

    }
  }
}
