#include "stdafx.h"
#include "stmcomponent.h"
#include "elementinput.h"
#include "spatial/point.h"
#include "spatial/linestring.h"
#include "stmmodel.h"
#include "element.h"
#include "temporal/timedata.h"

#include <QDebug>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::SpatioTemporal;


ElementInput::ElementInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           VariableType varType,
                           STMComponent *modelComponent)
  : TimeGeometryInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                            valueDefinition, modelComponent),
    m_component(modelComponent),
    m_varType(varType)
{

}

bool ElementInput::setProvider(HydroCouple::IOutput *provider)
{
  m_geometryMapping.clear();
  m_geometryMappingOrientation.clear();

  if(AbstractInput::setProvider(provider) && provider)
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

    if(timeGeometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

        if(lineString->pointCount())
        {
          HCPoint *p1 = lineString->pointInternal(0);
          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

          for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
          {
            if(!mapped[j])
            {
              ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

              IPoint *pp1 = lineStringProvider->point(0);
              IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);

              double deltap1p1 = hypot(p1->x() - pp1->x() , p1->y() - pp1->y());
              double deltap2p2 = hypot(p2->x() - pp2->x() , p2->y() - pp2->y());

              double deltap1p2 = hypot(p1->x() - pp2->x() , p1->y() - pp2->y());
              double deltap2p1 = hypot(p2->x() - pp1->x() , p2->y() - pp1->y());

              if( deltap1p1 < 1e-3 && deltap2p2 < 1e-3)
              {
                m_geometryMapping[i] = j;
                m_geometryMappingOrientation[i] = 1.0;
                mapped[j] = true;
                break;
              }
              else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
              {
                m_geometryMapping[i] = j;
                m_geometryMappingOrientation[i] = -1.0;
                mapped[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    return true;
  }

  return false;
}

bool ElementInput::canConsume(HydroCouple::IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());
  provider()->updateValues(this);
}

void ElementInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();
  ITimeGeometryComponentDataItem *timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider());

  int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
  int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

  double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
  double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

  if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
  {
    double factor = 0.0;

    if(providerCurrentTime > providerPreviousTime)
    {
      double denom = providerCurrentTime - providerPreviousTime;
      double numer =currentTime - providerPreviousTime;
      factor = numer / denom;
    }

    switch (m_varType)
    {
      case Flow:
        {
          for(auto it : m_geometryMapping)
          {
            double value1 = 0;
            double value2 = 0;

            timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
            timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

            double orientation = m_geometryMappingOrientation[it.first];
            value1 *= orientation;
            value2 *= orientation;

            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->flow = value2 + factor *(value1 - value2);

          }
        }
        break;
      case XSectionArea:
        {
          for(auto it : m_geometryMapping)
          {
            double value1 = 0;
            double value2 = 0;

            timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
            timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->xSectionArea = value2 + factor *(value1 - value2);
          }
        }
        break;
      case Depth:
        {
          for(auto it : m_geometryMapping)
          {
            double value1 = 0;
            double value2 = 0;

            timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
            timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->depth = value2 + factor *(value1 - value2);

          }
        }
        break;
      case TopWidth:
        {
          for(auto it : m_geometryMapping)
          {
            double value1 = 0;
            double value2 = 0;

            timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
            timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->width = value2 + factor *(value1 - value2);
          }
        }
        break;
    }
  }
  else
  {
    switch (m_varType)
    {
      case Flow:
        {
          for(auto it : m_geometryMapping)
          {
            double value = 0;
            timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
            value *= m_geometryMappingOrientation[it.first];
            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->flow = value;
          }
        }
        break;
      case XSectionArea:
        {
          for(auto it : m_geometryMapping)
          {
            double value = 0;
            timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->xSectionArea = value;
          }
        }
        break;
      case Depth:
        {
          for(auto it : m_geometryMapping)
          {
            double value = 0;
            timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->depth = value;
          }
        }
        break;
      case TopWidth:
        {
          for(auto it : m_geometryMapping)
          {
            double value = 0;
            timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
            Element *element =  m_component->modelInstance()->getElement(it.first);
            element->width = value;
          }
        }
        break;
    }
  }
}

ElementInput::VariableType ElementInput::variableType() const
{
  return m_varType;
}

void ElementInput::setVariableType(VariableType variableType)
{
  m_varType = variableType;
}

ElementHeatSourceInput::ElementHeatSourceInput(const QString &id,
                                               Dimension *timeDimension,
                                               Dimension *geometryDimension,
                                               ValueDefinition *valueDefinition,
                                               SourceType srcType,
                                               STMComponent *modelComponent)
  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                                 valueDefinition, modelComponent),
    m_component(modelComponent),
    m_srcType(srcType)
{

}

bool ElementHeatSourceInput::addProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

    std::unordered_map<int, int> geometryMapping;

    if(timeGeometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

        if(lineString->pointCount())
        {
          HCPoint *p1 = lineString->pointInternal(0);
          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

          for(int j = 0; j < timeGeometryDataItem->geometryCount() ; j++)
          {
            if(!mapped[j])
            {
              ILineString *lineStringProvider = dynamic_cast<ILineString*>(timeGeometryDataItem->geometry(j));

              IPoint *pp1 = lineStringProvider->point(0);
              IPoint *pp2 = lineStringProvider->point(lineStringProvider->pointCount() - 1);


              if(hypot(p1->x() - pp1->x() , p1->y() - pp1->y()) < 1e-3 && hypot(p2->x() - pp2->x() , p2->y() - pp2->y()) < 1e-3)
              {
                geometryMapping[i] = j;
                mapped[j] = true;
                break;
              }
              else if(hypot(p1->x() - pp2->x() , p1->y() - pp2->y()) < 1e-3 && hypot(p2->x() - pp1->x() , p2->y() - pp1->y()) < 1e-3)
              {
                geometryMapping[i] = j;
                mapped[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    m_geometryMapping[provider] = geometryMapping;

    return true;
  }

  return false;
}

bool ElementHeatSourceInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    m_geometryMapping.erase(provider);
    return true;
  }

  return false;
}

bool ElementHeatSourceInput::canConsume(IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  IGeometryComponentDataItem *geometryDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
          (geometryDataItem->geometryType() == IGeometry::LineString ||
           geometryDataItem->geometryType() == IGeometry::LineStringZ ||
           geometryDataItem->geometryType() == IGeometry::LineStringZM) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementHeatSourceInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

  for(IOutput *provider : m_providers)
  {
    provider->updateValues(this);
  }
}

void ElementHeatSourceInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();

  for(IOutput *provider : m_providers)
  {

    std::unordered_map<int,int> &geometryMapping = m_geometryMapping[provider];

    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeGeometryDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeGeometryDataItem->timeCount() - 2);

      double providerCurrentTime = timeGeometryDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeGeometryDataItem->time(previousTimeIndex)->julianDay();

      if(currentTime >=  providerPreviousTime && currentTime <= providerCurrentTime)
      {
        double factor = 0.0;

        if(providerCurrentTime > providerPreviousTime)
        {
          double denom = providerCurrentTime - providerPreviousTime;
          double numer = currentTime - providerPreviousTime;
          factor = numer / denom;
        }

        switch (m_srcType)
        {
          case RadiativeFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->radiationFluxes += value2 + factor *(value1 - value2);
              }
            }
            break;
          case HeatFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalHeatFluxes += value2 + factor *(value1 - value2);
              }
            }
            break;
        }
      }
      else
      {
        switch (m_srcType)
        {
          case RadiativeFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->radiationFluxes += value;

              }
            }
            break;
          case HeatFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalHeatFluxes += value;
              }
            }
            break;
        }
      }
    }
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)))
    {
      switch (m_srcType)
      {
        case RadiativeFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->radiationFluxes += value;

            }
          }
          break;
        case HeatFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->externalHeatFluxes += value;
            }
          }
          break;
      }
    }
  }
}

ElementHeatSourceInput::SourceType ElementHeatSourceInput::sourceType() const
{
  return m_srcType;
}

void ElementHeatSourceInput::setSourceType(SourceType srcType)
{
  m_srcType = srcType;
}
