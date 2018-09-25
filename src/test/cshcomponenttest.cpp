/*!
*  \file    CSHComponenttest.cpp
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
*  \todo Test transport on branching networks
*  \warning
*/

#include "stdafx.h"
#include "test/CSHComponenttest.h"
#include "cshmodel.h"
#include "elementjunction.h"
#include "element.h"
#include "variable.h"


void CSHComponentTest::versteegCase1_Upwind()
{
  QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    CSHModel *model = new CSHModel(nullptr);

    model->setInputFile(QFileInfo("../../examples/Versteeg/case1/case1_upwind.inp"));
    model->setOutputNetCDFFile(QFileInfo("../../examples/Versteeg/case1/case1_upwind.nc"));

    //initialize model
    if(model->initialize(errors))
    {
      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void CSHComponentTest::versteegCase1_Central()
{
  QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    CSHModel *model = new CSHModel(nullptr);

    model->setInputFile(QFileInfo("../../examples/Versteeg/case1/case1_central.inp"));
    model->setOutputNetCDFFile(QFileInfo("../../examples/Versteeg/case1/case1_central.nc"));

    //initialize model
    if(model->initialize(errors))
    {
      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void CSHComponentTest::versteegCase1_Hybrid()
{
  QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    CSHModel *model = new CSHModel(nullptr);

    model->setInputFile(QFileInfo("../../examples/Versteeg/case1/case1_hybrid.inp"));
    model->setOutputNetCDFFile(QFileInfo("../../examples/Versteeg/case1/case1_hybrid.nc"));

    //initialize model
    if(model->initialize(errors))
    {
      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void CSHComponentTest::versteegCase2_Upwind()
{
  QBENCHMARK_ONCE
  {
    //    std::list<std::string> errors;

    //    //Stream temperature model instance
    //    CSHModel *model = new CSHModel(nullptr);

    //    //Set advection discretization mode
    //    model->setAdvectionDiscretizationMode(CSHModel::AdvectionDiscretizationMode::Upwind);

    //    //set time and output variables
    //    model->setMaxTimeStep(0.5); //seconds
    //    model->setMinTimeStep(0.0001); //seconds
    //    model->setStartDateTime(0); // Modified Julian DateTime
    //    model->setEndDateTime(1.0 / 1440.0);// Modified Julian DateTime 1 minute
    //    model->setOutputInterval(1.0); //1 second

    //    //Output consider 1 date separate files easier. Cross tab
    //    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_upwind.csv"));

    //    //Set solver type
    //    model->heatSolver()->setSolverType(ODESolver::CVODE_ADAMS);

    //    //Set specific heat to 1.O
    //    model->setSpecificHeatCapacityWater(1.0);

    //    //  Domain spatial discretization
    //    double x , y = 0, z = 0;
    //    int numCells = 5.0;
    //    double dx = 1.0 / numCells;

    //    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    //    for(int i = 0; i <= numCells; i++)
    //    {
    //      //Make it explicit what each variable means
    //      x = i * dx;
    //      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

    //      if(i == 0)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 1.0;
    //      }
    //      else if(i == numCells)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 0.0;
    //      }
    //    }

    //    //Channels
    //    for(int i = 0; i < numCells ; i++)
    //    {
    //      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
    //      element->xSectionArea = 0.1;
    //      element->length = dx;
    //      element->flow = element->xSectionArea * 2.5;
    //      element->depth = 0.1;
    //      element->longDispersion.isBC = true;
    //      element->longDispersion.value = 0.1;
    //    }

    //    //initialize model
    //    if(model->initialize(errors))
    //    {
    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }

    //    //compare results

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void CSHComponentTest::versteegCase2_Central()
{
  QBENCHMARK_ONCE
  {
    //    std::list<std::string> errors;

    //    CSHModel *model = new CSHModel(nullptr);

    //    model->setAdvectionDiscretizationMode(CSHModel::AdvectionDiscretizationMode::Central);

    //    //set time and output variables
    //    model->setMaxTimeStep(0.5); //seconds
    //    model->setMinTimeStep(0.0001); //seconds
    //    model->setStartDateTime(0); // Modified Julian DateTime
    //    model->setEndDateTime(1.0 / 1440.0);//1 minute
    //    model->setOutputInterval(1.0); //1 second

    //    //Output consider 1 date separate files easier. Cross tab
    //    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_central.csv"));

    //    //Set solver type
    //    model->heatSolver()->setSolverType(ODESolver::CVODE_ADAMS);

    //    //Set specific heat to 1.O
    //    model->setSpecificHeatCapacityWater(1.0);

    //    //  Domain spatial discretization
    //    double x , y = 0, z = 0;
    //    int numCells = 5.0;
    //    double dx = 1.0 / numCells;

    //    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    //    for(int i = 0; i <= numCells; i++)
    //    {
    //      //Make it explicit what each variable means
    //      x = i * dx;
    //      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

    //      if(i == 0)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 1.0;
    //      }
    //      else if(i == numCells)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 0.0;
    //      }
    //    }

    //    //Channels
    //    for(int i = 0; i < numCells ; i++)
    //    {
    //      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
    //      element->xSectionArea = 0.1;
    //      element->length = dx;
    //      element->flow = element->xSectionArea * 2.5;
    //      element->depth = 0.1;
    //      element->longDispersion.isBC = true;
    //      element->longDispersion.value = 0.1;
    //    }

    //    //initialize model
    //    if(model->initialize(errors))
    //    {
    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }

    //    //compare results

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void CSHComponentTest::versteegCase2_Hybrid()
{
  QBENCHMARK_ONCE
  {

    //    //Error messages
    //    std::list<std::string> errors;

    //    //Stream temperature model instance
    //    CSHModel *model = new CSHModel(nullptr);

    //    //Set advection discretization mode
    //    model->setAdvectionDiscretizationMode(CSHModel::AdvectionDiscretizationMode::Hybrid);

    //    //set time and output variables
    //    model->setMaxTimeStep(0.5); //seconds
    //    model->setMinTimeStep(0.0001); //seconds
    //    model->setStartDateTime(0); // Modified Julian DateTime
    //    model->setEndDateTime(1.0 / 1440.0);//1 minute
    //    model->setOutputInterval(4.0); //1 second

    //    //Output consider 1 date separate files easier. Cross tab
    //    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_hybrid.csv"));
    //    model->setOutputNetCDFFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_hybrid.nc"));

    //    //Set solver type
    //    model->heatSolver()->setSolverType(ODESolver::CVODE_ADAMS);

    //    //Set specific heat to 1.O
    //    model->setSpecificHeatCapacityWater(1.0);

    //    //Domain spatial discretization
    //    //Element junction coordinates
    //    double x = 0, y = 0, z = 0;
    //    int numCells = 25.0;

    //    //grid size in x direction
    //    double dx = 1.0 / numCells;

    //    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    //    for(int i = 0; i <= numCells; i++)
    //    {
    //      //Make it explicit what each variable means
    //      x = i * dx;
    //      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

    //      //Apply boundary condition at the inlet junction
    //      if(i == 0)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 1.0;
    //      }
    //      //Apply boundary condition at the outlet junction
    //      else if(i == numCells)
    //      {
    //        eJunction->temperature.isBC = true;
    //        eJunction->temperature.value = 0.0;
    //      }
    //    }

    //    //Discretize elements using upstream and downstream junctions created earlier
    //    for(int i = 0; i < numCells ; i++)
    //    {
    //      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
    //      element->xSectionArea = 0.1;
    //      element->length = dx;
    //      element->flow = element->xSectionArea * 2.5;
    //      element->depth = 0.1;
    //      element->longDispersion.isBC = true;
    //      element->longDispersion.value = 0.1;
    //    }

    //    //initialize model
    //    if(model->initialize(errors))
    //    {

    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void CSHComponentTest::green_river_test()
{
  QBENCHMARK_ONCE
  {

    //    //Error messages
    //    std::list<std::string> errors;

    //    //Stream temperature model instance
    //    CSHModel *model = new CSHModel(nullptr);

    //    model->setInputFile(QFileInfo("../../examples/green_river_test/green_river_test.inp"));

    //    //initialize model
    //    if(model->initialize(errors))
    //    {
    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }
    //    else
    //    {
    //      for(std::string error : errors)
    //      {
    //        printf("%s\n", error.c_str());
    //      }
    //    }


    //    //compare results

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void CSHComponentTest::green_river_test1()
{
  QBENCHMARK_ONCE
  {

    //    //Error messages
    //    std::list<std::string> errors;

    //    //Stream temperature model instance
    //    CSHModel *model = new CSHModel(nullptr);

    //    model->setInputFile(QFileInfo("../../examples/green_river_test1/green_river_test1.inp"));

    //    //initialize model
    //    if(model->initialize(errors))
    //    {
    //      //Perform timestep until completion
    //      while (model->currentDateTime() < model->endDateTime())
    //      {
    //        model->update();
    //      }
    //    }
    //    else
    //    {
    //      for(std::string error : errors)
    //      {
    //        printf("%s\n", error.c_str());
    //      }
    //    }

    //    //finalize model
    //    model->finalize(errors);

    //    delete model;
  }
}

void CSHComponentTest::green_river_test2()
{
  QBENCHMARK_ONCE
  {

    //Error messages
    std::list<std::string> errors;

    //Stream temperature model instance
    CSHModel *model = new CSHModel(nullptr);

    model->setInputFile(QFileInfo("../../examples/green_river_test2/green_river_test2.inp"));

    //initialize model
    if(model->initialize(errors))
    {
      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }
    else
    {
      for(std::string error : errors)
      {
        printf("%s\n", error.c_str());
      }
    }

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

