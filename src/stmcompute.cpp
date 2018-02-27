#include "stmmodel.h"
#include "element.h"
#include "elementjunction.h"

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

void STMModel::prepareForNextTimeStep()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elementJunctions.size(); i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];

    {
      elementJunction->prevTemperature.copy(elementJunction->temperature);
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
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

    {
      element->prevTemperature.copy(element->temperature);
    }

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
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
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;
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
    double DSoluteDt = element->computeDSoluteDt(dt,y,solverUserData->variableIndex);
    dydt[element->index] = DSoluteDt;
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
