/*!
*  \file    stmproject.cpp
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
*  \todo
*  \warning
*/

#include "stmmodel.h"
#include "stmcomponent.h"
#include "spatial/point.h"
#include "spatial/network.h"
#include "element.h"
#include "elementjunction.h"
#include "spatial/edge.h"

using namespace std;

STMModel::STMModel(STMComponent *component)
  : QObject(component),
    m_timeStep(0.0001), //seconds
    m_maxTimeStep(0.5), //seconds
    m_minTimeStep(0.001), //seconds
    m_timeStepRelaxationFactor(0.80),
    m_numInitFixedTimeSteps(2),
    m_numCurrentInitFixedTimeSteps(0),
    m_computeDispersion(false),
    m_useAdaptiveTimeStep(true),
    m_converged(false),
    m_advectionMode(AdvectionDiscretizationMode::Upwind),
    m_numHeatElementJunctions(0),
    m_solver(nullptr),
    m_waterDensity(1.0), //kg/m^3
    m_cp(1.0/*for testing*/), //4187.0 J/kg/C
    m_component(component)
{
  m_solver = new ODESolver(1, ODESolver::RKQS);
}

STMModel::~STMModel()
{

  for(Element *element : m_elements)
    delete element;

  m_elements.clear();
  m_elementsById.clear();


  for(ElementJunction *elementJunction : m_elementJunctions)
    delete elementJunction;

  m_elementJunctions.clear();
  m_elementJunctionsById.clear();

  if(m_solver)
    delete m_solver;
}

double STMModel::minTimeStep() const
{
  return m_minTimeStep;
}

void STMModel::setMinTimeStep(double timeStep)
{
  m_minTimeStep = timeStep;
}

double STMModel::maxTimeStep() const
{
  return m_maxTimeStep;
}

void STMModel::setMaxTimeStep(double timeStep)
{
  m_maxTimeStep = timeStep;
}

bool STMModel::useAdaptiveTimeStep() const
{
  return m_useAdaptiveTimeStep;
}

void STMModel::setUseAdaptiveTimeStep(bool use)
{
  m_useAdaptiveTimeStep = use;
}

double STMModel::timeStepRelaxationFactor() const
{
  return m_timeStepRelaxationFactor;
}

void STMModel::setTimeStepRelaxationFactor(double tStepRelaxFactor)
{
  if(tStepRelaxFactor > 0)
    m_timeStepRelaxationFactor = tStepRelaxFactor;
}

double STMModel::currentTimeStep() const
{
  return m_timeStep;
}

double STMModel::startDateTime() const
{
  return m_startDateTime;
}

void STMModel::setStartDateTime(double dateTime)
{
  m_startDateTime = dateTime;
}

double STMModel::endDateTime() const
{
  return m_endDateTime;
}

void STMModel::setEndDateTime(double dateTime)
{
  m_endDateTime = dateTime;
}

double STMModel::outputInterval() const
{
  return m_outputInterval;
}

void STMModel::setOutputInterval(double interval)
{
  m_outputInterval = interval;
}

double STMModel::currentDateTime() const
{
  return m_currentDateTime;
}

ODESolver *STMModel::solver() const
{
  return m_solver;
}

bool STMModel::computeLongDispersion() const
{
  return m_computeDispersion;
}

void STMModel::setComputeLongDispersion(bool calculate)
{
  m_computeDispersion = calculate;

  for(Element *element : m_elements)
  {
    element->longDispersion.isBC = m_computeDispersion == false;
  }
}

STMModel::AdvectionDiscretizationMode STMModel::advectionDiscretizationMode() const
{
  return m_advectionMode;
}

void STMModel::setAdvectionDiscretizationMode(AdvectionDiscretizationMode advectionDiscretizationMode)
{
  m_advectionMode = advectionDiscretizationMode;
}

int STMModel::numSolutes() const
{
  return (int)m_solutes.size();
}

bool STMModel::addSolute(const string &soluteName)
{
  vector<string>::iterator it = find(m_solutes.begin() , m_solutes.end(), soluteName);

  if(it == m_solutes.end())
  {
    m_solutes.push_back(soluteName);
  }

  return false;
}

void STMModel::removeSolute(const string &soluteName)
{
  vector<string>::iterator it = find(m_solutes.begin() , m_solutes.end(), soluteName);

  if(it != m_solutes.end())
  {
    m_solutes.erase(it);
  }
}

vector<string> STMModel::solutes() const
{
  return m_solutes;
}

int STMModel::numElementJunctions() const
{
  return m_elementJunctions.size();
}

ElementJunction *STMModel::addElementJunction(const string &id, double x, double y, double z)
{
  if(m_elementJunctionsById.find(id) == m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = new ElementJunction(id, x, y, z, this);
    eJunction->index = m_elementJunctions.size();
    m_elementJunctions.push_back(eJunction);
    m_elementJunctionsById[id] = eJunction;
    return eJunction;
  }

  return nullptr;
}

void STMModel::deleteElementJunction(const string &id)
{
  std::unordered_map<string,ElementJunction*>::iterator eJIter =  m_elementJunctionsById.find(id) ;

  if(eJIter != m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = eJIter->second;
    m_elementJunctionsById.erase(eJIter);

    std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
    if(it != m_elementJunctions.end())
    {
      m_elementJunctions.erase(it);
    }

    delete eJunction;
  }
}

void STMModel::deleteElementJunction(int index)
{
  ElementJunction *eJunction = m_elementJunctions[index];

  m_elementJunctionsById.erase(eJunction->id);

  std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
  if(it != m_elementJunctions.end())
    m_elementJunctions.erase(it);

  delete eJunction;
}

ElementJunction *STMModel::getElementJunction(const string &id)
{
  return m_elementJunctionsById[id];
}

ElementJunction *STMModel::getElementJunction(int index)
{
  return m_elementJunctions[index];
}

int STMModel::numElements() const
{
  return m_elements.size();
}

Element *STMModel::addElement(const string &id, ElementJunction *upStream, ElementJunction *downStream)
{
  if(upStream && downStream)
  {
    Element *element = new Element(id, upStream, downStream, this);
    element->index = m_elements.size();
    m_elements.push_back(element);
    m_elementsById[id] = element;
    return element;
  }

  return nullptr;
}

void STMModel::deleteElement(const string &id)
{
  unordered_map<string,Element*>::iterator eIter = m_elementsById.find(id);

  if(eIter != m_elementsById.end())
  {
    Element *element = eIter->second;
    m_elementsById.erase(eIter);

    vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);
    if(it != m_elements.end())
      m_elements.erase(it);

    delete element;
  }
}

void STMModel::deleteElement(int index)
{
  Element *element = m_elements[index];
  m_elementJunctionsById.erase(element->id);

  vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);

  if(it != m_elements.end())
    m_elements.erase(it);

  delete element;
}

Element *STMModel::getElement(const string &id)
{
  return m_elementsById[id];
}

Element *STMModel::getElement(int index)
{
  return m_elements[index];
}

QFileInfo STMModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void STMModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

bool STMModel::initialize(list<string> &errors)
{
  bool initialized = initializeInputFiles(errors) &&
                     initializeTimeVariables(errors) &&
                     initializeElements(errors) &&
                     initializeSolver(errors) &&
                     initializeOutputFiles(errors);

  if(initialized)
  {
    applyInitialConditions();
  }

  return initialized;
}

void STMModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    prepareForNextTimeStep();

    m_timeStep = computeTimeStep();

    applyBoundaryConditions(m_currentDateTime + m_timeStep / 86400.0);

    computeLongDispersion();

    //Solve the transport for each element
    {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {
#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          //Solve heat transport first
          solveHeatTransport(m_timeStep);
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(size_t i = 0 ; i < m_solutes.size(); i++)
          {
            solveSoluteTransport(i, m_timeStep);
          }
        }
      }
    }

    //Solve continuity for each junction
    {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {
#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          //Solve heat transport first
          solveJunctionHeatContinuity(m_timeStep);
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(size_t i = 0 ; i < m_solutes.size(); i++)
          {
            solveJunctionSoluteContinuity(i, m_timeStep);
          }
        }
      }
    }

    m_currentDateTime +=  m_timeStep / 86400.0;

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime += m_outputInterval / 86400.0;
    }
  }
}

bool STMModel::finalize(std::list<string> &errors)
{
  closeOutputFiles();
  return true;
}

bool STMModel::initializeInputFiles(std::list<string> &errors)
{
  return true;
}

bool STMModel::initializeTimeVariables(std::list<string> &errors)
{
  if(m_startDateTime >= m_endDateTime)
  {
    errors.push_back("End datetime must be greater than startdatetime");
    return false;
  }

  if( (m_endDateTime - m_startDateTime) *  86400.0 < m_minTimeStep )
  {
    errors.push_back("Make sure timestep is less than the simulation interval");
    return false;
  }

  if(m_minTimeStep <=  0 || m_maxTimeStep <= 0)
  {
    errors.push_back("Make sure time steps are greater 0");
    return false;
  }

  if(m_minTimeStep >= m_maxTimeStep)
  {
    errors.push_back("");
    return false;
  }

  m_numCurrentInitFixedTimeSteps = 0;

  m_currentDateTime = m_startDateTime;
  m_nextOutputTime = m_currentDateTime;

  return true;
}

bool STMModel::initializeElements(std::list<string> &errors)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
    size_t numElements = elementJunction->incomingElements.size() + elementJunction->outgoingElements.size();

    switch (numElements)
    {
      case 0:
        elementJunction->junctionType = ElementJunction::NoElement;
        break;
      case 1:
        elementJunction->junctionType = ElementJunction::SingleElement;
        break;
      case 2:
        elementJunction->junctionType = ElementJunction::DoubleElement;
        break;
      default:
        elementJunction->junctionType = ElementJunction::MultiElement;
        break;
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];
    element->index = i;
    element->initialize();
  }

  //Set tempearture continuity junctions
  m_numHeatElementJunctions = 0;

  //Number of junctions where continuity needs to be enforced.
  m_numSoluteElementJunctions.resize(m_solutes.size(), 0);

  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
    elementJunction->initialize();

    size_t numJunctions = elementJunction->outgoingElements.size() + elementJunction->incomingElements.size();

    if(!elementJunction->temperature.isBC)
    {
      //If more than one junction solve continuity
      if(numJunctions > 2)
      {
        elementJunction->heatContinuityIndex = m_numHeatElementJunctions;
        m_numHeatElementJunctions++;
      }
      else
      {
        elementJunction->heatContinuityIndex = -1;
      }
    }
    else
    {
      //Otherwise don't solve continuity
      elementJunction->heatContinuityIndex = -1;
    }

    //Set solute indexes
    for(size_t j = 0 ; j < m_solutes.size(); j++)
    {
      if(!elementJunction->soluteConcs[j].isBC)
      {
        //If more than one junction solve continuity
        if(numJunctions > 2)
        {
          elementJunction->soluteContinuityIndexes[j] = m_numSoluteElementJunctions[j];
          m_numSoluteElementJunctions[j]++;
        }
        else
        {
          elementJunction->soluteContinuityIndexes[j] = -1;
        }
      }
      else
      {
        elementJunction->soluteContinuityIndexes[j] = -1;
      }
    }
  }

  return true;
}

bool STMModel::initializeSolver(std::list<string> &errors)
{
  m_solver->setSize(m_elements.size());
  m_solver->initialize();

  return true;
}

bool STMModel::initializeOutputFiles(std::list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool STMModel::initializeCSVOutputFile(std::list<string> &errors)
{
  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if(!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.dir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if(!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice( device);
    }

    if(m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      //      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Temperature";
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Temperature, UpstreamJunctionTemperature, DownstreamJunctionTemperature";

      for(size_t i = 0; i < m_solutes.size(); i++)
      {
        m_outputCSVStream << ", " << QString::fromStdString(m_solutes[i]);
      }

      m_outputCSVStream << endl;
    }

    m_outputCSVStream.flush();

    return true;
  }

  return false;
}

bool STMModel::initializeNetCDFOutputFile(std::list<string> &errors)
{
  return true;
}

void STMModel::prepareForNextTimeStep()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

//    if(!elementJunction->temperature.isBC)
    {
      elementJunction->prevTemperature.copy(elementJunction->temperature);
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
//      if(!elementJunction->soluteConcs[j].isBC)
      {
        elementJunction->prevSoluteConcs[j].copy(elementJunction->soluteConcs[j]);
      }
    }
  }


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

//    if(!element->temperature.isBC)
    {
      element->prevTemperature.copy(element->temperature);
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
//      if(!element->soluteConcs[j].isBC)
      {
        element->prevSoluteConcs[j].copy(element->soluteConcs[j]);
      }
    }
  }

}

void STMModel::applyInitialConditions()
{

  //Interpolate nodal temperatures and solute concentrations
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    if(!elementJunction->temperature.isBC)
    {
      elementJunction->interpTemp();
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      if(!elementJunction->soluteConcs[j].isBC)
      {
        elementJunction->interpSoluteConcs(j);
      }
    }
  }

  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;
}

void STMModel::applyBoundaryConditions(double dateTime)
{

}

double STMModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;

  double maxCourantFactor = 0.0;// Δx / v (s^-1)

  if(m_useAdaptiveTimeStep)
  {


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size()  ; i++)
    {
      Element *element = m_elements[i];
      double courantFactor = element->computeCourantFactor();
      //      double dispersionFactor = element->computeDispersionFactor();

      if(courantFactor > maxCourantFactor)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        maxCourantFactor = courantFactor;
      }

      //      if(dispersionFactor > maxCourantFactor)
      //      {
      //#ifdef USE_OPENMP
      //#pragma omp atomic read
      //#endif
      //        maxCourantFactor = dispersionFactor;
      //      }
    }

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;\

    if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
    {
      timeStep = std::min(timeStep, m_minTimeStep);
      m_numCurrentInitFixedTimeSteps++;
    }
  }

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(m_minTimeStep,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  return timeStep;
}

void STMModel::computeLongDispersion()
{
  if(m_computeDispersion)
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0 ; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      element->computeLongDispersion();
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->interpXSectionArea();
    elementJunction->interpLongDispersion();
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computePecletNumbers();
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeUpstreamPeclet();
    element->computeDownstreamPeclet();
  }
}

void STMModel::solveHeatTransport(double timeStep)
{
  //Allocate memory to store inputs and outputs
  double *currentTemperatures = new double[m_elements.size()];
  double *outputTemperatures = new double[m_elements.size()];

  //Set initial input and output values to current values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentTemperatures[element->index] = element->temperature.value;
    outputTemperatures[element->index] = element->temperature.value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;
  m_solver->solve(currentTemperatures, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputTemperatures, &STMModel::computeDTDt, &solverUserData);

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    double outputTemperature = outputTemperatures[element->index];
    element->temperature.value = outputTemperature;
  }

  //Delete allocated memory
  delete[] currentTemperatures;
  delete[] outputTemperatures;
}

void STMModel::solveSoluteTransport(int soluteIndex, double timeStep)
{
  double *currentSoluteConcs = new double[m_elements.size()];
  double *outputSoluteConcs = new double[m_elements.size()];

  //Set initial values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    outputSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;
  m_solver->solve(outputSoluteConcs, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                  outputSoluteConcs, &STMModel::computeDSoluteDt, &solverUserData);


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->soluteConcs[soluteIndex].value = outputSoluteConcs[element->index];
  }


  //Delete allocated memory
  delete[] currentSoluteConcs;
  delete[] outputSoluteConcs;

}

void STMModel::computeDTDt(double t, double y[], double dydt[], void* userData)
{

  SolverUserData *solverUserData = (SolverUserData*) userData;
  STMModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DTDt = element->computeDTDt(dt,y);
    dydt[element->index] = DTDt;
  }
}

void STMModel::computeDSoluteDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  STMModel *modelInstance = solverUserData->model;
  double dt = t - modelInstance->m_currentDateTime *  86400.0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    double DTDt = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
    dydt[element->index] = DTDt;
  }

}

void STMModel::solveJunctionHeatContinuity(double timeStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    if(elementJunction->heatContinuityIndex > -1)
    {
      elementJunction->solveHeatContinuity(timeStep);
    }
    else if(!elementJunction->temperature.isBC)
    {
      elementJunction->interpTemp();
    }
  }
}

void STMModel::solveJunctionSoluteContinuity(int soluteIndex, double timeStep)
{

}

void STMModel::writeOutput()
{
  writeCSVOutput();
  writeNetCDFOutput();
}

void STMModel::writeCSVOutput()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    for(size_t i = 0 ; i < m_elements.size() ; i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        <<  ", " << element->temperature.value << ", " << element->upstreamJunction->temperature.value << ", "
                         << element->downstreamJunction->temperature.value;

      for(size_t j = 0 ; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value;
      }

      m_outputCSVStream << endl;
    }
  }
}

void STMModel::writeNetCDFOutput()
{

}

void STMModel::closeOutputFiles()
{
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeCSVOutputFile();
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      closeOutputNetCDFFile();
    }
  }
}

void STMModel::closeCSVOutputFile()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void STMModel::closeOutputNetCDFFile()
{

}
