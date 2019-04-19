/*!
   *  \file    CSHComponent.cpp
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

#include "stdafx.h"
#include "cshcomponent.h"
#include "cshmodel.h"

#include "core/dimension.h"
#include "core/valuedefinition.h"
#include "spatial/linestring.h"
#include "spatial/point.h"
#include "element.h"
#include "elementjunction.h"
#include "core/idbasedargument.h"
#include "progresschecker.h"
#include "elementinput.h"
#include "elementoutput.h"
#include "temporal/timedata.h"
#include "core/abstractoutput.h"
#include "core/unit.h"
#include "core/unitdimensions.h"
#include "elementoutput.h"

using namespace HydroCouple;

CSHComponent::CSHComponent(const QString &id, CSHComponentInfo *modelComponentInfo)
  : AbstractTimeModelComponent(id, modelComponentInfo),
    m_inputFilesArgument(nullptr),
    m_flowInput(nullptr),
    m_xSectionAreaInput(nullptr),
    m_topWidthInput(nullptr),
    m_radiationFluxUnit(nullptr),
    m_heatFluxUnit(nullptr),
    m_temperatureUnit(nullptr),
    m_externalRadiationFluxInput(nullptr),
    m_externalHeatFluxInput(nullptr),
    m_temperatureOutput(nullptr),
    m_modelInstance(nullptr),
    m_parent(nullptr)
{
  m_timeDimension = new Dimension("TimeDimension",this);
  m_geometryDimension = new Dimension("ElementGeometryDimension", this);

  m_radiationFluxUnit = new Unit(this);
  m_radiationFluxUnit->setCaption("Radiation Flux (W/m^2)");
  m_radiationFluxUnit->setConversionFactorToSI(1.0);
  m_radiationFluxUnit->setOffsetToSI(0.0);
  m_radiationFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_radiationFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -3.0);

  m_heatFluxUnit = new Unit(this);
  m_heatFluxUnit->setCaption("Heat Source (W or J/s)");
  m_heatFluxUnit->setConversionFactorToSI(1.0);
  m_heatFluxUnit->setOffsetToSI(0.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Length, 2.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -3.0);

  m_temperatureUnit = new Unit(this);
  m_temperatureUnit->setCaption("Temperature (°C)");
  m_temperatureUnit->setConversionFactorToSI(1.0);
  m_temperatureUnit->setOffsetToSI(273.15);
  m_temperatureUnit->dimensionsInternal()->setPower(HydroCouple::Temperature, 1.0);

  m_soluteFluxUnit = new Unit(this);
  m_soluteFluxUnit->setCaption("Mass Flux (kg/s)");
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -1.0);

  m_soluteUnit = new Unit(this);
  m_soluteUnit->setCaption("Concentration (kg/m^3)");
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Length, -3.0);

  m_soluteConcQuantity = new Quantity(QVariant::Double, m_soluteUnit, this);
  m_soluteConcFluxQuantity = new Quantity(QVariant::Double, m_soluteFluxUnit, this);

  createArguments();
}

CSHComponent::~CSHComponent()
{

  initializeFailureCleanUp();

  while (m_clones.size())
  {
    CSHComponent *clone =  dynamic_cast<CSHComponent*>(m_clones.first());
    removeClone(clone);
    delete clone;
  }

  if(m_parent)
  {
    m_parent->removeClone(this);
    m_parent = nullptr;
  }
}

QList<QString> CSHComponent::validate()
{
  if(isInitialized())
  {
    setStatus(IModelComponent::Validating,"Validating...");

    //check connections

    setStatus(IModelComponent::Valid,"");
  }
  else
  {
    //throw has not been initialized yet.
  }

  return QList<QString>();
}

void CSHComponent::prepare()
{
  if(!isPrepared() && isInitialized() && m_modelInstance)
  {
    for(auto output :  outputsInternal())
    {
      for(auto adaptedOutput : output->adaptedOutputs())
      {
        adaptedOutput->initialize();
      }
    }

    updateOutputValues(QList<HydroCouple::IOutput*>());

    setStatus(IModelComponent::Updated ,"Finished preparing model");
    setPrepared(true);
  }
  else
  {
    setPrepared(false);
    setStatus(IModelComponent::Failed ,"Error occured when preparing model");
  }
}

void CSHComponent::update(const QList<HydroCouple::IOutput*> &requiredOutputs)
{
  if(status() == IModelComponent::Updated)
  {
    setStatus(IModelComponent::Updating);

    double minConsumerTime = std::max(m_modelInstance->currentDateTime(), getMinimumConsumerTime());

    while (m_modelInstance->currentDateTime() <= minConsumerTime &&
           m_modelInstance->currentDateTime() < m_modelInstance->endDateTime())
    {
      m_modelInstance->update();

      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
    }

    updateOutputValues(requiredOutputs);

    currentDateTimeInternal()->setJulianDay(m_modelInstance->currentDateTime());

    if(m_modelInstance->currentDateTime() >=  m_modelInstance->endDateTime())
    {
      setStatus(IModelComponent::Done , "Simulation finished successfully", 100);
    }
    else
    {
      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
      else
      {
        setStatus(IModelComponent::Updated);
      }
    }
  }
}

void CSHComponent::finish()
{
  if(isPrepared())
  {
    setStatus(IModelComponent::Finishing , "CSHComponent with id " + id() + " is being disposed" , 100);

    std::list<std::string> errors;
    m_modelInstance->finalize(errors);
    initializeFailureCleanUp();

    setPrepared(false);
    setInitialized(false);

    setStatus(IModelComponent::Finished , "CSHComponent with id " + id() + " has been disposed" , 100);
    setStatus(IModelComponent::Created , "CSHComponent with id " + id() + " ran successfully and has been re-created" , 100);
  }
}

CSHModel *CSHComponent::modelInstance() const
{
  return m_modelInstance;
}

ICloneableModelComponent *CSHComponent::parent() const
{
  return m_parent;
}

ICloneableModelComponent *CSHComponent::clone()
{
  if(isInitialized())
  {
    CSHComponent *cloneComponent = dynamic_cast<CSHComponent*>(componentInfo()->createComponentInstance());
    cloneComponent->setReferenceDirectory(referenceDirectory());

    IdBasedArgumentString *identifierArg = identifierArgument();
    IdBasedArgumentString *cloneIndentifierArg = cloneComponent->identifierArgument();

    (*cloneIndentifierArg)["Id"] = QString((*identifierArg)["Id"]);
    (*cloneIndentifierArg)["Caption"] = QString((*identifierArg)["Caption"]);
    (*cloneIndentifierArg)["Description"] = QString((*identifierArg)["Description"]);

    QString appendName = "_clone_" + QString::number(m_clones.size()) + "_" + QUuid::createUuid().toString().replace("{","").replace("}","");

    //(*cloneComponent->m_inputFilesArgument)["Input File"] = QString((*m_inputFilesArgument)["Input File"]);

    QString inputFilePath = QString((*m_inputFilesArgument)["Input File"]);
    QFileInfo inputFile = getAbsoluteFilePath(inputFilePath);

    if(inputFile.absoluteDir().exists())
    {
      QString suffix = "." + inputFile.completeSuffix();
      inputFilePath = inputFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      QFile::copy(inputFile.absoluteFilePath(), inputFilePath);
      (*cloneComponent->m_inputFilesArgument)["Input File"] = inputFilePath;
    }

    QString outputNetCDFFilePath = QString((*m_inputFilesArgument)["Output NetCDF File"]);
    QFileInfo outputNetCDFFile = getAbsoluteFilePath(outputNetCDFFilePath);

    if(!outputNetCDFFilePath.isEmpty() && outputNetCDFFile.absoluteDir().exists())
    {
      QString suffix = "." + outputNetCDFFile.completeSuffix();
      outputNetCDFFilePath = outputNetCDFFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      (*cloneComponent->m_inputFilesArgument)["Output NetCDF File"] = outputNetCDFFilePath;
    }

    QString  outputCSVFilePath = QString((*m_inputFilesArgument)["Output CSV File"]);
    QFileInfo outputCSVFile = getAbsoluteFilePath(outputCSVFilePath);

    if(!outputCSVFilePath.isEmpty() && outputCSVFile.absoluteDir().exists())
    {
      QString suffix = "." + outputCSVFile.completeSuffix();
      outputCSVFilePath = outputCSVFile.absoluteFilePath().replace(suffix,"") + appendName + suffix;
      (*cloneComponent->m_inputFilesArgument)["Output CSV File"] = outputCSVFilePath;
    }


    cloneComponent->m_parent = this;
    m_clones.append(cloneComponent);

    emit propertyChanged("Clones");

    cloneComponent->initialize();

    return cloneComponent;
  }

  return nullptr;
}

QList<ICloneableModelComponent*> CSHComponent::clones() const
{
  return m_clones;
}

bool CSHComponent::removeClone(CSHComponent *component)
{
  int removed;

#ifdef USE_OPENMP
#pragma omp critical (CSHComponent)
#endif
  {
    removed = m_clones.removeAll(component);
  }


  if(removed)
  {
    emit propertyChanged("Clones");
  }

  return removed;
}

void CSHComponent::initializeFailureCleanUp()
{
  if(m_modelInstance)
  {
    delete m_modelInstance;
    m_modelInstance = nullptr;
  }

}

void CSHComponent::createArguments()
{
  createInputFileArguments();
}

void CSHComponent::createInputFileArguments()
{
  QStringList fidentifiers;
  fidentifiers.append("Input File");
  fidentifiers.append("Output NetCDF File");
  fidentifiers.append("Output CSV File");

  Quantity *fquantity = Quantity::unitLessValues("InputFilesQuantity", QVariant::String, this);
  fquantity->setDefaultValue("");
  fquantity->setMissingValue("");

  Dimension *dimension = new Dimension("IdDimension","Dimension for identifiers",this);

  m_inputFilesArgument = new IdBasedArgumentString("InputFiles", fidentifiers, dimension, fquantity, this);
  m_inputFilesArgument->setCaption("Model Input Files");
  m_inputFilesArgument->addFileFilter("Input File (*.inp)");
  m_inputFilesArgument->setMatchIdentifiersWhenReading(true);

  addArgument(m_inputFilesArgument);
}

bool CSHComponent::initializeArguments(QString &message)
{
  bool initialized = initializeInputFilesArguments(message);

  if(initialized)
  {
    createGeometries();

    for(AbstractOutput *output : outputsInternal())
      output->updateValues();
  }
  else
  {
    initializeFailureCleanUp();
  }

  return initialized;
}

bool CSHComponent::initializeInputFilesArguments(QString &message)
{
  QString inputFilePath = QString((*m_inputFilesArgument)["Input File"]);
  QFileInfo inputFile = getAbsoluteFilePath(inputFilePath);

  if(inputFile.exists())
  {
    initializeFailureCleanUp();

    m_modelInstance = new CSHModel(this);
    m_modelInstance->setInputFile(inputFile);

    QString netCDFOutput = QString((*m_inputFilesArgument)["Output NetCDF File"]);
    if(!netCDFOutput.isEmpty() && !netCDFOutput.isNull())
      m_modelInstance->setOutputNetCDFFile(QFileInfo(netCDFOutput));

    QString csvOutput = QString((*m_inputFilesArgument)["Output CSV File"]);
    if(!csvOutput.isEmpty() && !csvOutput.isNull())
      m_modelInstance->setOutputCSVFile(QFileInfo(csvOutput));

    std::list<std::string> errors;
    bool initialized = m_modelInstance->initialize(errors);

    for (std::string errorMsg : errors)
    {
      message += "/n" + QString::fromStdString(errorMsg);
    }

    if(initialized)
    {
      timeHorizonInternal()->setJulianDay(m_modelInstance->startDateTime());
      timeHorizonInternal()->setDuration(m_modelInstance->endDateTime() - m_modelInstance->startDateTime());
      currentDateTimeInternal()->setJulianDay(m_modelInstance->startDateTime());
      progressChecker()->reset(m_modelInstance->startDateTime(), m_modelInstance->endDateTime());
    }


    return initialized;
  }
  else
  {
    message = "Input file does not exist: " + inputFile.absoluteFilePath();
    return false;
  }

  return true;
}

void CSHComponent::createGeometries()
{
  m_elementGeometries.clear();
  m_elementJunctionGeometries.clear();

  for(int i = 0; i < m_modelInstance->numElements() ; i++)
  {
    Element *element = m_modelInstance->getElement(i);
    ElementJunction *from = element->upstreamJunction;
    ElementJunction *to   = element->downstreamJunction;

    HCLineString *lineString = new HCLineString(QString::fromStdString(element->id));
    lineString->setMarker(i);
    HCPoint *p1 = new HCPoint(from->x , from->y, QString::fromStdString(from->id), lineString);
    HCPoint *p2 = new HCPoint(to->x , to->y, QString::fromStdString(to->id), lineString);
    lineString->addPoint(p1);
    lineString->addPoint(p2);

    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(from->x , from->y, from->z, QString::fromStdString(from->id), nullptr)));
    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(to->x , to->y, to->z, QString::fromStdString(to->id), nullptr)));

    m_elementGeometries.push_back(QSharedPointer<HCGeometry>(lineString));
  }
}

void CSHComponent::createInputs()
{
  createFlowInput();
  createXSectionAreaInput();
  createDepthInput();
  createTopWidthInput();
  createExternalRadiationFluxInput();
  createExternalHeatFluxInput();
  createDVDtInput();

  for(int i = 0; i < m_modelInstance->numSolutes() ; i++)
  {
    createExternalSoluteFluxInput(i);
  }

  createWaterAgeFluxInput();
}

void CSHComponent::createFlowInput()
{
  Quantity *flowQuantity = Quantity::flowInCMS(this);

  m_flowInput = new ElementInput("ElementFlowInput",
                                 m_timeDimension,
                                 m_geometryDimension,
                                 flowQuantity,
                                 ElementInput::Flow,
                                 this);

  m_flowInput->setCaption("Element Flow (m^3/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_flowInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_flowInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_flowInput);

  m_flowInput->addTime(dt1);
  m_flowInput->addTime(dt2);

  addInput(m_flowInput);
}

void CSHComponent::createXSectionAreaInput()
{
  Quantity *xSectionQuantity = Quantity::areaInSquareMeters(this);

  m_xSectionAreaInput = new ElementInput("ElementXSectionAreaInput",
                                         m_timeDimension,
                                         m_geometryDimension,
                                         xSectionQuantity,
                                         ElementInput::XSectionArea,
                                         this);

  m_xSectionAreaInput->setCaption("Element Cross-Section Area(m^2)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_xSectionAreaInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_xSectionAreaInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_xSectionAreaInput);

  m_xSectionAreaInput->addTime(dt1);
  m_xSectionAreaInput->addTime(dt2);

  addInput(m_xSectionAreaInput);
}

void CSHComponent::createDepthInput()
{
  Quantity *depthQuantity = Quantity::lengthInMeters(this);

  m_depthInput = new ElementInput("ElementDepthInput",
                                  m_timeDimension,
                                  m_geometryDimension,
                                  depthQuantity,
                                  ElementInput::Depth,
                                  this);

  m_depthInput->setCaption("Element Depth (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_depthInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_depthInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_depthInput);

  m_depthInput->addTime(dt1);
  m_depthInput->addTime(dt2);

  addInput(m_depthInput);
}

void CSHComponent::createTopWidthInput()
{
  Quantity *widthQuantity = Quantity::lengthInMeters(this);

  m_topWidthInput = new ElementInput("ElementTopWidthInput",
                                     m_timeDimension,
                                     m_geometryDimension,
                                     widthQuantity,
                                     ElementInput::TopWidth,
                                     this);

  m_topWidthInput->setCaption("Element Top Width (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_topWidthInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_topWidthInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_topWidthInput);

  m_topWidthInput->addTime(dt1);
  m_topWidthInput->addTime(dt2);

  addInput(m_topWidthInput);
}

void CSHComponent::createExternalRadiationFluxInput()
{

  Quantity *radiationQuantity = new Quantity(QVariant::Double, m_radiationFluxUnit, this);

  m_externalRadiationFluxInput = new ElementSourceInput("RadiationFluxInput",
                                                        m_timeDimension,
                                                        m_geometryDimension,
                                                        radiationQuantity,
                                                        ElementSourceInput::RadiativeFlux,
                                                        this);
  m_externalRadiationFluxInput->setCaption("External Radiation Flux (W/m^2)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_externalRadiationFluxInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_externalRadiationFluxInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_externalRadiationFluxInput);

  m_externalRadiationFluxInput->addTime(dt1);
  m_externalRadiationFluxInput->addTime(dt2);

  addInput(m_externalRadiationFluxInput);
}

void CSHComponent::createDVDtInput()
{
  Quantity *flowQuantity = Quantity::lengthInMeters(this);

  m_DVDtInput = new ElementInput("ElementVolumeTimeDerivativeInput",
                                 m_timeDimension,
                                 m_geometryDimension,
                                 flowQuantity,
                                 ElementInput::DVolumeDTime,
                                 this);

  m_DVDtInput->setCaption("Element Volume Time Derivative (m^3/s)");
  m_DVDtInput->setDescription("Element Volume Time Derivative (m^3/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_DVDtInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_DVDtInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_DVDtInput);

  m_DVDtInput->addTime(dt1);
  m_DVDtInput->addTime(dt2);

  addInput(m_DVDtInput);
}

void CSHComponent::createExternalHeatFluxInput()
{
  Quantity *heatFluxQuantity = new Quantity(QVariant::Double, m_heatFluxUnit, this);

  m_externalHeatFluxInput = new ElementSourceInput("HeatFluxInput",
                                                   m_timeDimension,
                                                   m_geometryDimension,
                                                   heatFluxQuantity,
                                                   ElementSourceInput::HeatFlux,
                                                   this);

  m_externalHeatFluxInput->setCaption("External Heat Flux (J/s)");


  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_externalHeatFluxInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_externalHeatFluxInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_externalHeatFluxInput);

  m_externalHeatFluxInput->addTime(dt1);
  m_externalHeatFluxInput->addTime(dt2);

  addInput(m_externalHeatFluxInput);
}

void CSHComponent::createExternalSoluteFluxInput(int soluteIndex)
{
  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementSourceInput *soluteFluxInput  = new ElementSourceInput(soluteName + "Input",
                                                                m_timeDimension,
                                                                m_geometryDimension,
                                                                m_soluteConcFluxQuantity,
                                                                ElementSourceInput::SoluteFlux,
                                                                this);
  soluteFluxInput->setCaption("Element " + soluteName + " Flux (kg/s)");
  soluteFluxInput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  soluteFluxInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteFluxInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteFluxInput);

  soluteFluxInput->addTime(dt1);
  soluteFluxInput->addTime(dt2);

  addInput(soluteFluxInput);
}

void CSHComponent::createWaterAgeFluxInput()
{
  if(m_modelInstance->simulateWaterAge())
  {
    Quantity *timeQuantity = Quantity::unitLessValues("Time/Time", QVariant::Double, this);

    ElementSourceInput *soluteFluxInput  = new ElementSourceInput("WaterAgeFluxInput",
                                                                  m_timeDimension,
                                                                  m_geometryDimension,
                                                                  timeQuantity,
                                                                  ElementSourceInput::SoluteFlux,
                                                                  this);

    soluteFluxInput->setCaption("Element Water Age Flux (days/s)");
    soluteFluxInput->setSoluteIndex(m_modelInstance->numSolutes());

    QList<QSharedPointer<HCGeometry>> geometries;

    for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
    {
      geometries.append(lineString);
    }

    soluteFluxInput->addGeometries(geometries);

    SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteFluxInput);
    SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteFluxInput);

    soluteFluxInput->addTime(dt1);
    soluteFluxInput->addTime(dt2);

    addInput(soluteFluxInput);
  }
}

void CSHComponent::createOutputs()
{
  createTemperatureOutput();
  createFlowOutput();
  createXSectionAreaOutput();
  createDepthOutput();
  createTopWidthOutput();

  for(int i = 0; i < m_modelInstance->numSolutes() ; i++)
  {
    createSoluteConcOutput(i);
  }

  createWaterAgeFluxOutput();
}

void CSHComponent::createTemperatureOutput()
{
  Quantity *temperatureQuantity = new Quantity(QVariant::Double, m_temperatureUnit, this);

  m_temperatureOutput = new ElementOutput("TemperatureOutput",
                                          m_timeDimension,
                                          m_geometryDimension,
                                          temperatureQuantity,
                                          ElementOutput::Temperature,
                                          this);

  m_temperatureOutput->setCaption("Element Temperature (°C)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_temperatureOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, m_temperatureOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_temperatureOutput);

  m_temperatureOutput->addTime(dt1);
  m_temperatureOutput->addTime(dt2);

  addOutput(m_temperatureOutput);

}

void CSHComponent::createSoluteConcOutput(int index)
{

  QString soluteName = QString::fromStdString(m_modelInstance->solute(index));

  ElementOutput *soluteOutput  = new ElementOutput(soluteName + "_Output",
                                                   m_timeDimension,
                                                   m_geometryDimension,
                                                   m_soluteConcQuantity,
                                                   ElementOutput::SoluteConc,
                                                   this);
  soluteOutput->setCaption("Element " + soluteName + " Concentration (kg/m^3)");
  soluteOutput->setSoluteIndex(index);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  soluteOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteOutput);

  soluteOutput->addTime(dt1);
  soluteOutput->addTime(dt2);

  addOutput(soluteOutput);
}

void CSHComponent::createWaterAgeFluxOutput()
{
  if(m_modelInstance->simulateWaterAge())
  {
    Quantity *timeQuantity = Quantity::timeInDays(this);
    ElementOutput *soluteOutput  = new ElementOutput("WaterAgeOutput",
                                                     m_timeDimension,
                                                     m_geometryDimension,
                                                     timeQuantity ,
                                                     ElementOutput::SoluteConc,
                                                     this);
    soluteOutput->setCaption("Element Water Age (days)");
    soluteOutput->setSoluteIndex(m_modelInstance->numSolutes());

    QList<QSharedPointer<HCGeometry>> geometries;

    for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
    {
      geometries.append(lineString);
    }

    soluteOutput->addGeometries(geometries);

    SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteOutput);
    SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteOutput);

    soluteOutput->addTime(dt1);
    soluteOutput->addTime(dt2);

    addOutput(soluteOutput);
  }
}

void CSHComponent::createFlowOutput()
{
  Quantity *flowQuantity = Quantity::flowInCMS(this);

  ElementOutput *flowOutput = new ElementOutput("ElementFlowOutput",
                                   m_timeDimension,
                                 m_geometryDimension,
                                 flowQuantity,
                                 ElementOutput::Flow,
                                 this);

  flowOutput->setCaption("Element Flow (m^3/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  flowOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, flowOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), flowOutput);

  flowOutput->addTime(dt1);
  flowOutput->addTime(dt2);

  addOutput(flowOutput);
}

void CSHComponent::createXSectionAreaOutput()
{
  Quantity *xSectionQuantity = Quantity::areaInSquareMeters(this);

  ElementOutput *xSectionAreaOutput = new ElementOutput("ElementXSectionAreaOutput",
                                         m_timeDimension,
                                         m_geometryDimension,
                                         xSectionQuantity,
                                         ElementOutput::XSectionArea,
                                         this);

  xSectionAreaOutput->setCaption("Element Cross-Section Area (m^2)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  xSectionAreaOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, xSectionAreaOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), xSectionAreaOutput);

  xSectionAreaOutput->addTime(dt1);
  xSectionAreaOutput->addTime(dt2);

  addOutput(xSectionAreaOutput);
}

void CSHComponent::createDepthOutput()
{
  Quantity *depthQuantity = Quantity::lengthInMeters(this);

  ElementOutput *depthOutput = new ElementOutput("ElementDepthOutput",
                                  m_timeDimension,
                                  m_geometryDimension,
                                  depthQuantity,
                                  ElementOutput::Depth,
                                  this);

  depthOutput->setCaption("Element Depth (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  depthOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, depthOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), depthOutput);

  depthOutput->addTime(dt1);
  depthOutput->addTime(dt2);

  addOutput(depthOutput);
}

void CSHComponent::createTopWidthOutput()
{
  Quantity *widthQuantity = Quantity::lengthInMeters(this);

  ElementOutput *topWidthOutput = new ElementOutput("ElementTopWidthOutput",
                                     m_timeDimension,
                                     m_geometryDimension,
                                     widthQuantity,
                                     ElementOutput::TopWidth,
                                     this);

  topWidthOutput->setCaption("Element Top Width (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  topWidthOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, topWidthOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), topWidthOutput);

  topWidthOutput->addTime(dt1);
  topWidthOutput->addTime(dt2);

  addOutput(topWidthOutput);
}
