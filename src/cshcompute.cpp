/*!
*  \file    CSHcompute.cpp
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
#include "cshmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "iboundarycondition.h"
#include "cshcomponent.h"


#ifdef USE_OPENMP
#include <omp.h>
#endif


using namespace std;

void CSHModel:: update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    applyBoundaryConditions(m_currentDateTime);

    if(m_component)
      m_component->applyInputValues();

    m_prevTimeStep = m_timeStep;

    m_timeStep = computeTimeStep();

    computeDerivedHydraulics();

    computeLongDispersion();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int j = 0; j < 2; j++)
    {
      switch (j)
      {
        case 0:
          {
            solveJunctionHeatContinuity(m_timeStep);
          }
          break;
        case 1:
          {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int i = 0 ; i < (int)m_solutes.size(); i++)
            {
              solveJunctionSoluteContinuity(i, m_timeStep);
            }
          }
          break;
      }
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int j = 0; j < 2; j++)
    {
      switch (j)
      {
        case 0:
          {
            solveHeatTransport(m_timeStep);
          }
          break;
        case 1:
          {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int i = 0 ; i < (int)m_solutes.size(); i++)
            {
              solveSoluteTransport(i, m_timeStep);
            }
          }
          break;
      }
    }

    m_prevDateTime = m_currentDateTime;
    m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime = std::min(m_nextOutputTime + m_outputInterval / 86400.0 , m_endDateTime);
    }

    if(m_verbose)
    {
      printStatus();
    }
  }
}

void CSHModel::prepareForNextTimeStep()
{

  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  for(int i = 0 ; i < (int)m_elementJunctions.size(); i++)
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

    element->computeDVolumeDt();

    element->computeHeatBalance(m_timeStep);
    m_totalHeatBalance += element->totalHeatBalance;
    m_totalRadiationHeatBalance += element->totalRadiationFluxesHeatBalance;
    m_totalAdvDispHeatBalance += element->totalAdvDispHeatBalance;
    m_totalEvaporationHeatBalance += element->totalEvaporativeHeatFluxesBalance;
    m_totalConvectiveHeatBalance += element->totalConvectiveHeatFluxesBalance;
    m_totalExternalHeatFluxBalance += element->totalExternalHeatFluxesBalance;

    element->prevTemperature.copy(element->temperature);

    m_minTemp = min(m_minTemp , element->temperature.value);
    m_maxTemp = max(m_maxTemp , element->temperature.value);

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->computeSoluteBalance(m_timeStep, j);
      m_totalSoluteMassBalance[j] += element->totalSoluteMassBalance[j];
      m_totalAdvDispSoluteMassBalance[j] += element->totalAdvDispSoluteMassBalance[j];
      m_totalExternalSoluteFluxMassBalance[j] += element->totalExternalSoluteFluxesMassBalance[j];

      element->prevSoluteConcs[j].copy(element->soluteConcs[j]);

      m_minSolute[j] = min(m_minSolute[j] , element->soluteConcs[j].value);
      m_maxSolute[j] = max(m_maxSolute[j] , element->soluteConcs[j].value);
    }
  }
}

void CSHModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_totalHeatBalance = 0.0;
  m_totalRadiationHeatBalance = 0.0;
  m_totalEvaporationHeatBalance = 0.0;
  m_totalConvectiveHeatBalance = 0.0;
  m_totalAdvDispHeatBalance = 0.0;
  m_totalExternalHeatFluxBalance = 0.0;

  std::fill(m_totalSoluteMassBalance.begin(), m_totalSoluteMassBalance.end(), 0.0);
  std::fill(m_totalAdvDispSoluteMassBalance.begin(), m_totalAdvDispSoluteMassBalance.end(), 0.0);
  std::fill(m_totalExternalSoluteFluxMassBalance.begin(), m_totalExternalSoluteFluxMassBalance.end(), 0.0);

  //Interpolate nodal temperatures and solute concentrations
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  for(int i = 0 ; i < (int)m_elementJunctions.size(); i++)
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

void CSHModel::applyBoundaryConditions(double dateTime)
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

double CSHModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;
  double maxCourantFactor = 0.0;

  if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
  {
    timeStep = m_minTimeStep;
    m_numCurrentInitFixedTimeSteps++;
  }
  else if(m_useAdaptiveTimeStep)
  {
#ifdef _WIN32
    {

      for(int i = 0 ; i < (int)m_elements.size()  ; i++)
      {
        Element *element = m_elements[i];
        double courantFactor = element->computeCourantFactor() + element->computeDispersionFactor();

        if(!std::isinf(courantFactor) && courantFactor > maxCourantFactor)
        {
          maxCourantFactor = courantFactor;
        }
      }
    }
#else
    {
      //#ifdef USE_OPENMP
      //#pragma omp parallel for
      //#endif
      for(int i = 0 ; i < (int)m_elements.size()  ; i++)
      {
        Element *element = m_elements[i];
        double courantFactor = element->computeCourantFactor() + element->computeDispersionFactor();

        if(!(std::isinf(courantFactor) || std::isnan(courantFactor)) && courantFactor > maxCourantFactor)
        {

          //#ifdef USE_OPENMP
          //#pragma omp atomic read
          //#endif
          maxCourantFactor = courantFactor;

        }
      }
    }
#endif

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;
  }

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(m_minTimeStep,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  return timeStep;
}

void CSHModel::computeDerivedHydraulics()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeDerivedHydraulics();
  }

}

void CSHModel::computeLongDispersion()
{
  if(m_computeDispersion)
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      element->computeLongDispersion();
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computePecletNumbers();
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->computeUpstreamPeclet();
    element->computeDownstreamPeclet();
  }
}

void CSHModel::solveHeatTransport(double timeStep)
{
  //Allocate memory to store inputs and outputs
  m_dheatPrevTime = m_currentDateTime * 86400.0;

  //Set initial input and output values to current values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    m_currTemps[element->index] = element->temperature.value;
    m_outTemps[element->index] = element->temperature.value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;

  if(m_heatSolver->solve(m_currTemps.data(), m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                         m_outTemps.data(), &CSHModel::computeDTDt, &solverUserData))
  {
    printf("CSH Temperature Solver failed \n");
  }

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    double outputTemperature = m_outTemps[element->index];
    element->temperature.value = outputTemperature;
  }

}

void CSHModel::solveSoluteTransport(int soluteIndex, double timeStep)
{
  std::vector<double> &currentSoluteConcs = m_currSoluteConcs[soluteIndex];
  std::vector<double> &outputSoluteConcs = m_outSoluteConcs[soluteIndex];
  m_dsolutePrevTimes[soluteIndex] = m_currentDateTime *  86400.0;

  //Set initial values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    currentSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
    outputSoluteConcs[element->index] = element->soluteConcs[soluteIndex].value;
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;

  if(m_soluteSolvers[soluteIndex]->solve(outputSoluteConcs.data(), m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                                         outputSoluteConcs.data(), &CSHModel::computeDSoluteDt, &solverUserData))
  {
    printf("CSH Solute Solver failed \n");
  }


  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->soluteConcs[soluteIndex].value = outputSoluteConcs[element->index];
  }
}

void CSHModel::computeDTDt(double t, double y[], double dydt[], void* userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  CSHModel *modelInstance = solverUserData->model;
  double dt = t - modelInstance->m_dheatPrevTime;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    dydt[element->index] = element->computeDTDt(dt,y);
  }

  modelInstance->m_dheatPrevTime = t;
}

void CSHModel::computeDSoluteDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  CSHModel *modelInstance = solverUserData->model;
  double dt = t - modelInstance->m_dsolutePrevTimes[solverUserData->variableIndex];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];
    dydt[element->index] = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
  }

  modelInstance->m_dsolutePrevTimes[solverUserData->variableIndex] = t;
}

void CSHModel::solveJunctionHeatContinuity(double timeStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elementJunctions.size(); i++)
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

void CSHModel::solveJunctionSoluteContinuity(int soluteIndex, double timeStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elementJunctions.size(); i++)
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
