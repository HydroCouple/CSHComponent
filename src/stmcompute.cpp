/*!
*  \file    stmcompute.cpp
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
#include "element.h"
#include "elementjunction.h"
#include "iboundarycondition.h"

using namespace std;

void STMModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    m_timeStep = computeTimeStep();

    double tempCurrentDateTime = m_currentDateTime + m_timeStep / 86400.0;

    applyBoundaryConditions(tempCurrentDateTime);

    //Retrieve external data from other coupled models

    if(m_retrieveCouplingDataFunction)
    {
      (*m_retrieveCouplingDataFunction)(this, tempCurrentDateTime);
    }

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

    m_prevDateTime = m_currentDateTime;
    m_currentDateTime = tempCurrentDateTime;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime += m_outputInterval / 86400.0;
    }

    if(m_verbose)
    {
      printStatus();
    }
  }
}

void STMModel::prepareForNextTimeStep()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->prevTemperature.copy(elementJunction->temperature);

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      elementJunction->prevSoluteConcs[j].copy(elementJunction->soluteConcs[j]);
    }
  }

  m_minTemp = std::numeric_limits<double>::max();
  m_maxTemp = std::numeric_limits<double>::lowest();

  std::fill(m_maxSolute.begin(), m_maxSolute.end(), m_maxTemp);
  std::fill(m_minSolute.begin(), m_minSolute.end(), m_minTemp);

  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    element->computeHeatBalance();
    m_heatBalance += element->heatBalance;

    element->prevTemperature.copy(element->temperature);

    element->radiationFluxes = 0.0;
    element->externalHeatFluxes = 0.0;

    m_minTemp = min(m_minTemp , element->temperature.value);
    m_maxTemp = max(m_maxTemp , element->temperature.value);

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->computeSoluteBalance(j);
      m_soluteMassBalance[j] += element->soluteMassBalance[j];

      element->prevSoluteConcs[j].copy(element->soluteConcs[j]);

      element->externalSoluteFluxes[j] = 0.0;

      m_minSolute[j] = min(m_minSolute[j] , element->soluteConcs[j].value);
      m_maxSolute[j] = max(m_maxSolute[j] , element->soluteConcs[j].value);
    }
  }
}

void STMModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_heatBalance = 0.0;

  for(size_t i = 0; i < m_solutes.size(); i++)
  {
    m_soluteMassBalance[i] = 0.0;
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->heatBalance = 0.0;

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->soluteMassBalance[j] = 0.0;
    }
  }

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

  applyBoundaryConditions(m_currentDateTime);

  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;
}

void STMModel::applyBoundaryConditions(double dateTime)
{
  //reset external fluxes
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->externalHeatFluxes = 0.0;
    element->radiationFluxes = 0.0;

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->externalSoluteFluxes[j] = 0.0;
    }
  }

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->applyBoundaryConditions(dateTime);
  }
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
      //double dispersionFactor = element->computeDispersionFactor();

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
  m_heatSolver->solve(currentTemperatures, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
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
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;
  m_soluteSolvers[soluteIndex]->solve(outputSoluteConcs, m_elements.size() , m_currentDateTime * 86400.0, timeStep,
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
    dydt[element->index] = element->computeDTDt(dt,y);
  }
}

void STMModel::computeDSoluteDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  STMModel *modelInstance = solverUserData->model;
  double dt = t - modelInstance->m_currentDateTime *  86400.0;

#ifdef USE_OPENMP
  //#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    dydt[element->index] = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
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
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    if(elementJunction->soluteContinuityIndexes[soluteIndex] > -1)
    {
      elementJunction->solveSoluteContinuity(soluteIndex, timeStep);
    }
    else if(!elementJunction->soluteConcs[soluteIndex].isBC)
    {
      elementJunction->interpSoluteConcs(soluteIndex);
    }
  }
}
