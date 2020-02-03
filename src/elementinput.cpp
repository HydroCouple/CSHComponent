#include "stdafx.h"
#include "cshcomponent.h"
#include "elementinput.h"
#include "spatial/point.h"
#include "spatial/linestring.h"
#include "cshmodel.h"
#include "element.h"
#include "temporal/timedata.h"
#include "hydrocouple.h"
#include "hydrocoupletemporal.h"

#include <QDebug>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::Temporal;
using namespace HydroCouple::SpatioTemporal;


ElementInput::ElementInput(const QString &id,
                           Dimension *timeDimension,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           VariableType varType,
                           CSHComponent *modelComponent)
  : TimeGeometryInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                            valueDefinition, modelComponent),
    m_component(modelComponent),
    m_varType(varType)
{

}

ElementInput::~ElementInput()
{

}

bool ElementInput::setProvider(HydroCouple::IOutput *provider)
{
  m_geometryMapping.clear();
  m_geometryMappingOrientation.clear();

  if(AbstractInput::setProvider(provider) && provider)
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
       timeGeometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(timeGeometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));
        Element *element = m_component->modelInstance()->getElement(i);

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

                if(m_varType == VariableType::DVolumeDTime)
                {
                  element->dvolume_dt.isBC = true;
                }

                break;
              }
              else if(deltap1p2 < 1e-3 &&  deltap2p1 < 1e-3)
              {
                m_geometryMapping[i] = j;
                m_geometryMappingOrientation[i] = -1.0;
                mapped[j] = true;

                if(m_varType == VariableType::DVolumeDTime)
                {
                  element->dvolume_dt.isBC = true;
                }

                break;
              }
            }
          }
        }
      }
    }
    else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)) )
    {
      QStringList identifiers = timeIdBasedDataItem->identifiers();
      std::vector<bool> mapped(identifiers.size(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        Element *element = m_component->modelInstance()->getElement(i);

        for(int j = 0; j < identifiers.size(); j++)
        {
          if(!mapped[j])
          {
            QString pId = identifiers[j];

            if(pId.toStdString() == element->id)
            {
              m_geometryMapping[i] = j;
              m_geometryMappingOrientation[i] = 1.0;
              mapped[j] = true;

              if(m_varType == VariableType::DVolumeDTime)
              {
                element->dvolume_dt.isBC = true;
              }

              break;
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
  ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;


  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
     (timeGeometryDataItem->geometryType() == IGeometry::LineString ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZ ||
      timeGeometryDataItem->geometryType() == IGeometry::LineStringZM) &&
     (provider->valueDefinition()->type() == QVariant::Double ||
      provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)) &&
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
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;

  if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider())))
  {
    double currentTime = m_component->modelInstance()->currentDateTime();

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
              element->flow.value = value2 + factor *(value1 - value2);

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
        case STSXSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSXSectionArea = value2 + factor *(value1 - value2);
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
        case STSWidthFraction:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSWidthFraction = value2 + factor *(value1 - value2);
            }
          }
          break;
        case DVolumeDTime:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->dvolume_dt.value = value2 + factor *(value1 - value2);
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
              element->flow.value = value;
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
        case STSXSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSXSectionArea = value;
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
        case STSWidthFraction:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSWidthFraction = value;
            }
          }
          break;
        case DVolumeDTime:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeGeometryDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->dvolume_dt.value = value;
            }
          }
          break;
      }
    }
  }
  else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider())))
  {
    double currentTime = m_component->modelInstance()->currentDateTime();

    int currentTimeIndex = timeIdBasedDataItem->timeCount() - 1;
    int previousTimeIndex = std::max(0 , timeIdBasedDataItem->timeCount() - 2);

    double providerCurrentTime = timeIdBasedDataItem->time(currentTimeIndex)->julianDay();
    double providerPreviousTime = timeIdBasedDataItem->time(previousTimeIndex)->julianDay();

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

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              double orientation = m_geometryMappingOrientation[it.first];
              value1 *= orientation;
              value2 *= orientation;

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->flow.value = value2 + factor *(value1 - value2);

            }
          }
          break;
        case XSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->xSectionArea = value2 + factor *(value1 - value2);
            }
          }
          break;
        case STSXSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSXSectionArea = value2 + factor *(value1 - value2);
            }
          }
          break;
        case Depth:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

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

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->width = value2 + factor *(value1 - value2);
            }
          }
          break;
        case STSWidthFraction:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSWidthFraction = value2 + factor *(value1 - value2);
            }
          }
          break;
        case DVolumeDTime:
          {
            for(auto it : m_geometryMapping)
            {
              double value1 = 0;
              double value2 = 0;

              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
              timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->dvolume_dt.value = value2 + factor *(value1 - value2);
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
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              value *= m_geometryMappingOrientation[it.first];
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->flow.value = value;
            }
          }
          break;
        case XSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->xSectionArea = value;
            }
          }
          break;
        case STSXSectionArea:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSXSectionArea = value;
            }
          }
          break;
        case Depth:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
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
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->width = value;
            }
          }
          break;
        case STSWidthFraction:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->STSWidthFraction = value;
            }
          }
          break;
        case DVolumeDTime:
          {
            for(auto it : m_geometryMapping)
            {
              double value = 0;
              timeIdBasedDataItem->getValue(currentTimeIndex,it.second, & value);
              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->dvolume_dt.value = value;
            }
          }
          break;
      }
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


ElementSourceInput::ElementSourceInput(const QString &id,
                                       Dimension *timeDimension,
                                       Dimension *geometryDimension,
                                       ValueDefinition *valueDefinition,
                                       SourceType srcType,
                                       CSHComponent *modelComponent)
  : TimeGeometryMultiInputDouble(id, IGeometry::LineString, timeDimension, geometryDimension,
                                 valueDefinition, modelComponent),
    m_component(modelComponent),
    m_srcType(srcType),
    m_soluteIndex(0)
{

}

ElementSourceInput::~ElementSourceInput()
{

}

bool ElementSourceInput::addProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
    IIdBasedComponentDataItem *idBasedDataItem = nullptr;

    std::unordered_map<int, int> geometryMapping;

    if((timeGeometryDataItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider)) &&
       timeGeometryDataItem->geometryCount())
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
    else if((geometryDataItem = dynamic_cast<IGeometryComponentDataItem*>(provider)) &&
            geometryDataItem->geometryCount())
    {
      std::vector<bool> mapped(geometryDataItem->geometryCount(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        HCLineString *lineString = dynamic_cast<HCLineString*>(getGeometry(i));

        if(lineString->pointCount())
        {
          HCPoint *p1 = lineString->pointInternal(0);
          HCPoint *p2 = lineString->pointInternal(lineString->pointCount() - 1);

          for(int j = 0; j < geometryDataItem->geometryCount() ; j++)
          {
            if(!mapped[j])
            {
              ILineString *lineStringProvider = dynamic_cast<ILineString*>(geometryDataItem->geometry(j));

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
    else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {
      QStringList identifiers = timeIdBasedDataItem->identifiers();

      std::vector<bool> mapped(identifiers.size(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        Element *element = m_component->modelInstance()->getElement(i);

        for(int j = 0; j < identifiers.size() ; j++)
        {
          if(!mapped[j])
          {
            QString providerId = identifiers[j];

            if(element->id == providerId.toStdString())
            {
              geometryMapping[i] = j;
              mapped[j] = true;
              break;
            }
          }
        }
      }
    }
    else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)))
    {
      QStringList identifiers = idBasedDataItem->identifiers();

      std::vector<bool> mapped(identifiers.size(), false);

      for(int i = 0; i < geometryCount() ; i++)
      {
        Element *element = m_component->modelInstance()->getElement(i);

        for(int j = 0; j < identifiers.size() ; j++)
        {
          if(!mapped[j])
          {
            QString providerId = identifiers[j];

            if(element->id == providerId.toStdString())
            {
              geometryMapping[i] = j;
              mapped[j] = true;
              break;
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

bool ElementSourceInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    m_geometryMapping.erase(provider);
    return true;
  }

  return false;
}

bool ElementSourceInput::canConsume(IOutput *provider, QString &message) const
{
  ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
  IGeometryComponentDataItem *geometryDataItem = nullptr;
  ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
  IIdBasedComponentDataItem *idBasedDataItem = nullptr;

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
  else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }
  else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)) &&
          (provider->valueDefinition()->type() == QVariant::Double ||
           provider->valueDefinition()->type() == QVariant::Int))
  {
    return true;
  }

  message = "Provider must be a LineString";

  return false;
}

void ElementSourceInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();
  int currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_component->modelInstance()->currentDateTime());

  for(IOutput *provider : m_providers)
  {
    provider->updateValues(this);
  }
}

void ElementSourceInput::applyData()
{
  double currentTime = m_component->modelInstance()->currentDateTime();

  for(IOutput *provider : m_providers)
  {

    std::unordered_map<int,int> &geometryMapping = m_geometryMapping[provider];

    ITimeGeometryComponentDataItem *timeGeometryDataItem = nullptr;
    IGeometryComponentDataItem *geometryDataItem = nullptr;
    ITimeIdBasedComponentDataItem *timeIdBasedDataItem = nullptr;
    IIdBasedComponentDataItem *idBasedDataItem = nullptr;

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
                double interpValue = value2 + factor *(value1 - value2);
                element->externalHeatFluxes += interpValue;
              }
            }
            break;
          case SoluteFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeGeometryDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                double interpValue = value2 + factor *(value1 - value2);
                element->externalSoluteFluxes[m_soluteIndex] += interpValue;
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
          case SoluteFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeGeometryDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalSoluteFluxes[m_soluteIndex] += value;
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
        case SoluteFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              geometryDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->externalSoluteFluxes[m_soluteIndex] += value;
            }
          }
          break;
      }
    }
    else if((timeIdBasedDataItem = dynamic_cast<ITimeIdBasedComponentDataItem*>(provider)))
    {
      int currentTimeIndex = timeIdBasedDataItem->timeCount() - 1;
      int previousTimeIndex = std::max(0 , timeIdBasedDataItem->timeCount() - 2);

      double providerCurrentTime = timeIdBasedDataItem->time(currentTimeIndex)->julianDay();
      double providerPreviousTime = timeIdBasedDataItem->time(previousTimeIndex)->julianDay();

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

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

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

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                double interpValue = value2 + factor *(value1 - value2);
                element->externalHeatFluxes += interpValue;
              }
            }
            break;
          case SoluteFlux:
            {
              for(auto it : geometryMapping)
              {
                double value1 = 0;
                double value2 = 0;

                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value1);
                timeIdBasedDataItem->getValue(previousTimeIndex,it.second, &value2);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                double interpValue = value2 + factor *(value1 - value2);
                element->externalSoluteFluxes[m_soluteIndex] += interpValue;
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
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value);

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
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalHeatFluxes += value;
              }
            }
            break;
          case SoluteFlux:
            {
              for(auto it : geometryMapping)
              {
                double value = 0;
                timeIdBasedDataItem->getValue(currentTimeIndex,it.second, &value);

                Element *element =  m_component->modelInstance()->getElement(it.first);
                element->externalSoluteFluxes[m_soluteIndex] += value;
              }
            }
            break;
        }
      }
    }
    else if((idBasedDataItem = dynamic_cast<IIdBasedComponentDataItem*>(provider)))
    {
      switch (m_srcType)
      {
        case RadiativeFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, &value);

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
              idBasedDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->externalHeatFluxes += value;
            }
          }
          break;
        case SoluteFlux:
          {
            for(auto it : geometryMapping)
            {
              double value = 0;
              idBasedDataItem->getValue(it.second, &value);

              Element *element =  m_component->modelInstance()->getElement(it.first);
              element->externalSoluteFluxes[m_soluteIndex] += value;
            }
          }
          break;
      }
    }
  }
}

ElementSourceInput::SourceType ElementSourceInput::sourceType() const
{
  return m_srcType;
}

void ElementSourceInput::setSourceType(SourceType srcType)
{
  m_srcType = srcType;
}

int ElementSourceInput::soluteIndex() const
{
  return m_soluteIndex;
}

void ElementSourceInput::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}
