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
#include "odesolver.h"

using namespace std;

STMModel::STMModel(STMComponent *component)
  : QObject(component),
    m_component(component)
{
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
  return initializeInputFiles(errors) &&
      initializeTimeVariables(errors) &&
      initializeElements(errors) &&
      initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

void STMModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    prepareForNextTimeStep();

    m_timeStep = computeTimeStep();

    applyBoundaryConditions(m_currentDateTime + m_timeStep / 86400.0);

    computeLongDispersion();

    bool converged = false;
    bool *soluteConverged = new bool[solutes().size()]; std::fill(soluteConverged, soluteConverged + m_solutes.size(), false);

    int iterations = 0;

    while (!converged && iterations < m_maxNumberOfIterations)
    {
      bool tempConverged = false;

      int numSoluteConverged = 0;

      //Solve junction continuity
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {
#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          //Solve heat continuity first
          solveJunctionHeatContinuity(m_timeStep, tempConverged);
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          //Solve for each solute
          for(size_t i = 0 ; i < m_solutes.size(); i++)
          {
            bool soluteConverged = false;
            solveJunctionSoluteContinuity(i, m_timeStep, soluteConverged);

            if(soluteConverged)
            {
#ifdef USE_OPENMP
#pragma omp atomic
#endif
              numSoluteConverged++;
            }
          }
        }
      }

      converged = tempConverged && numSoluteConverged == (int)m_solutes.size();

      //Solve the transport for each element
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {
#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          //Solve heat transport first
          solveElementHeatTransport(m_timeStep);
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
            solveElementSoluteTransport(i, m_timeStep);
          }
        }
      }

      iterations++;
    }

    delete[] soluteConverged;

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

  m_currentDateTime = m_startDateTime;

  return true;
}

bool STMModel::initializeElements(std::list<string> &errors)
{

  //Initialize solutes
  m_numSoluteContinuityElementJunctions.clear();
  m_numSoluteContinuityElementJunctions.resize(m_solutes.size(), 0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
    size_t numElements = elementJunction->incomingElements.size() + elementJunction->outgoingElements.size();

    switch (numElements) {
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
    element->setUpstreamElement();
    element->setDownStreamElement();
  }

  //Set tempearture continuity junctions
  m_numHeatContinuityElementJunctions = 0;

  for(size_t i = 0 ; i < m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;

    size_t numJunctions = elementJunction->outgoingElements.size() + elementJunction->incomingElements.size();

    if(!elementJunction->temperature.isBC)
    {
      //If more than one junction solve continuity
      if(numJunctions > 1)
      {
        elementJunction->heatContinuityIndex = m_numHeatContinuityElementJunctions;
        m_numHeatContinuityElementJunctions++;

        for(size_t j = 0 ; j < m_solutes.size(); j++)
        {
          elementJunction->soluteContinuityIndexes[j] = m_numSoluteContinuityElementJunctions[j];
          m_numSoluteContinuityElementJunctions[j]++;
        }
      }
      else
      {
        //Otherwise if one junction but inflow boundary condition then solve continuity
        if(numJunctions == 1)
        {
          if(elementJunction->externalHeatFlux.isBC)
          {
            elementJunction->heatContinuityIndex = m_numHeatContinuityElementJunctions;
            m_numHeatContinuityElementJunctions++;
          }
          else
          {
            for(size_t j = 0 ; j < m_solutes.size(); j++)
            {
              if(elementJunction->externalSoluteFluxes[j].isBC)
              {
                elementJunction->soluteContinuityIndexes[j] =  m_numSoluteContinuityElementJunctions[j];
                m_numSoluteContinuityElementJunctions[j]++;
              }
              else
              {
                elementJunction->soluteContinuityIndexes[j] = -1;
              }
            }
          }
        }
        //Otherwise don't solve continuity
        else
        {
          elementJunction->heatContinuityIndex = -1;

          for(size_t j = 0 ; j < m_solutes.size(); j++)
          {
            elementJunction->soluteContinuityIndexes[j] = -1;
          }
        }
      }
    }
  }

  return true;
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
      m_outputCSVStream << "ElementId, ElementIndex, Temperature";

      for(size_t i = 0; i < m_solutes.size(); i++)
      {
        m_outputCSVStream << ", " << QString::fromStdString(m_solutes[i]);
      }

      m_outputCSVStream << endl;
    }

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

}

void STMModel::applyInitialConditions()
{
  writeOutput();
}

void STMModel::applyBoundaryConditions(double dateTime)
{

}

double STMModel::computeTimeStep() const
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

      //add dispersion factor

      if(courantFactor > maxCourantFactor)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        maxCourantFactor = courantFactor;
      }
    }

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::min(m_minTimeStep,  m_nextOutputTime - m_currentDateTime);
  }

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

void STMModel::solveJunctionHeatContinuity(double timeStep, bool &converged)
{
  double *currentTemperatures = new double[m_numHeatContinuityElementJunctions];
  double *outputTemperatures = new double[m_numHeatContinuityElementJunctions];

  //Set T_o and calculate initial gradient using euler method.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    if(elementJunction->heatContinuityIndex > -1)
    {
      currentTemperatures[elementJunction->heatContinuityIndex] = elementJunction->temperature.value;
      outputTemperatures[elementJunction->heatContinuityIndex] = elementJunction->temperature.value;
    }
    else if(!elementJunction->temperature.isBC && (elementJunction->outgoingElements.size() || elementJunction->incomingElements.size()))
    {
      elementJunction->interpTemp();
    }
  }

  ODESolver solver(m_numHeatContinuityElementJunctions, ODESolver::RKQS);
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;
  solver.solve(currentTemperatures, m_numHeatContinuityElementJunctions, m_currentDateTime * 86400.0, timeStep, outputTemperatures, &STMModel::computeJunctionDYDt, &solverUserData);


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    if(elementJunction->heatContinuityIndex > -1)
    {
      elementJunction->temperature.value = outputTemperatures[elementJunction->heatContinuityIndex];
    }
  }


  //Check for convergence.
  converged = true;

  for(int i = 0 ; i < m_numHeatContinuityElementJunctions; i++)
  {
    if(fabs(currentTemperatures[i] - outputTemperatures[i]) > m_convergenceTol)
    {
      converged = false;
      break;
    }
  }

  delete[] currentTemperatures;
  delete[] outputTemperatures;
}

void STMModel::solveJunctionSoluteContinuity(int soluteIndex, double timeStep, bool &converged)
{
  int numSoluteContinuityJunctions =  m_numSoluteContinuityElementJunctions[soluteIndex];

  double *currentSoluteConcs = new double[numSoluteContinuityJunctions];
  double *outputSoluteConcs = new double[numSoluteContinuityJunctions];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    int soluteContinuityIndex = elementJunction->soluteContinuityIndexes[soluteIndex];

    if(soluteContinuityIndex > -1)
    {
      currentSoluteConcs[soluteContinuityIndex] = elementJunction->soluteConcs[soluteIndex].value;
      outputSoluteConcs[soluteContinuityIndex] = elementJunction->soluteConcs[soluteIndex].value;
    }
    else if(!elementJunction->soluteConcs[soluteIndex].isBC && (elementJunction->outgoingElements.size() || elementJunction->incomingElements.size()))
    {
      elementJunction->interpSoluteConcs(soluteIndex);
    }
  }

  ODESolver solver(numSoluteContinuityJunctions, ODESolver::RKQS);
  SolverUserData solverUserData; solverUserData.model = this, solverUserData.variableIndex = soluteIndex;
  solver.solve(currentSoluteConcs, numSoluteContinuityJunctions, m_currentDateTime * 86400.0, timeStep, outputSoluteConcs, &STMModel::computeJunctionDYDt, &solverUserData);

  converged = true;

  for(int i = 0 ; i < numSoluteContinuityJunctions; i++)
  {
    if(fabs(currentSoluteConcs[i] - outputSoluteConcs[i]) > m_convergenceTol)
    {
      converged = false;
      break;
    }
  }

  delete[] currentSoluteConcs;
  delete[] outputSoluteConcs;
}

void STMModel::computeJunctionDYDt(double t, double y[], double dydt[], void* userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  STMModel *modelInstance = solverUserData->model;
  int variableIndex = solverUserData->variableIndex;

  double dt = t - modelInstance->m_currentDateTime *  86400.0;

  switch (variableIndex)
  {
    //-1 Is for calculating heat gradient
    case -1:
      {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(size_t i = 0 ; i < modelInstance->m_elementJunctions.size(); i++)
        {
          ElementJunction *elementJunction = modelInstance->m_elementJunctions[i];

          //Calculate continuity only for valid junctions.
          if(elementJunction->heatContinuityIndex  > -1)
          {
            dydt[elementJunction->heatContinuityIndex] = elementJunction->computeDTDt(dt,y);
          }
        }
      }
      break;
      // Everything else is for calculating solute gradient
    default:
      {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(size_t i = 0 ; i < modelInstance->m_elementJunctions.size(); i++)
        {
          ElementJunction *elementJunction = modelInstance->m_elementJunctions[i];
          int soluteContinuityIndex = elementJunction->soluteContinuityIndexes[variableIndex];

          //Calculate continuity only for valid junctions.
          if(soluteContinuityIndex > -1)
          {
            dydt[soluteContinuityIndex] = elementJunction->computeDSoluteDt(dt, y, variableIndex);
          }
        }
      }
      break;
  }

}

void STMModel::solveElementHeatTransport(double timeStep)
{
  double *currentTemperatures = new double[m_elements.size()];
  double *outputTemperatures = new double[m_elements.size()];

  //Set T_o and calculate initial gradient using euler method.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentTemperatures[element->index] = element->temperature.value;
    outputTemperatures[element->index] = element->temperature.value;
  }

  ODESolver solver(m_elements.size(), ODESolver::RKQS);
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;
  solver.solve(currentTemperatures, m_elements.size(), m_currentDateTime * 86400.0, timeStep, outputTemperatures, &STMModel::computeElementDYDt, &solverUserData);


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->temperature.value = outputTemperatures[element->index];
  }

  delete[] currentTemperatures;
  delete[] outputTemperatures;
}

void STMModel::solveElementSoluteTransport(int soluteIndex, double timeStep)
{
  double *currentSoluteConcs = new double[m_elements.size()];
  double *outputSoluteConcs = new double[m_elements.size()];
  double *DSolutedt = new double[m_elements.size()];

  //Set Solute_o and calculate initial gradient using euler method.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    outputSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    DSolutedt[element->index] = element->computeDSoluteDt(timeStep, soluteIndex);
  }

  ODESolver solver(m_elements.size(), ODESolver::RKQS);
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;
  solver.solve(currentSoluteConcs, m_elements.size(), m_currentDateTime * 86400.0, timeStep, outputSoluteConcs, &STMModel::computeElementDYDt, &solverUserData);


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->soluteConcs[soluteIndex].value = outputSoluteConcs[element->index];
  }

  delete[] currentSoluteConcs;
  delete[] outputSoluteConcs;
  delete[] DSolutedt;
}

void STMModel::computeElementDYDt(double t, double y[], double dydt[], void* userData)
{
  SolverUserData *solverUserData = (SolverUserData*)userData;
  STMModel *modelInstance = solverUserData->model;
  int variableIndex = solverUserData->variableIndex;

  double dt = t - modelInstance->m_currentDateTime *  86400.0;

  switch (variableIndex)
  {
    case -1:
      {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
        {
          Element *element = modelInstance->m_elements[i];
          dydt[element->index] = element->computeDTDt(dt,y);
        }
      }
      break;
    default:
      {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
        {
          Element *element = modelInstance->m_elements[i];
          dydt[element->index] = element->computeDSoluteDt(dt,y, variableIndex);
        }
      }
      break;
  }
}

void STMModel::writeOutput()
{
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      writeCSVOutput();
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      writeNetCDFOutput();
    }
  }
}

void STMModel::writeCSVOutput()
{
  if(m_outputCSVStream.device()->isOpen())
  {
    for(size_t i = 0 ; i < m_elements.size() ; i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << QString::fromStdString(element->id) << ", " << element->index << ", " << element->temperature.value ;

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
  }
}

void STMModel::closeOutputNetCDFFile()
{

}
