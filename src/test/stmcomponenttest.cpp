#include "stdafx.h"
#include "test/stmcomponenttest.h"
#include "odesolver.h"
#include "stmmodel.h"
#include "elementjunction.h"
#include "element.h"
#include "variable.h"

void STMComponentTest::solveODERK4_Prob1()
{
  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RK4);
    solver.initialize();

    double y = -1.0;
    double y_out = y;
    double t = 0.0;
    double dt = 0.01;
    double maxt = 1.1;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

      double y_anal = problem1(t + dt);
      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("RK4 Problem 1 Error: %1").arg(error).toStdString().c_str());
  }
}

void STMComponentTest::solveODERKQS_Prob1()
{

  ODESolver solver(1, ODESolver::RKQS);
  solver.setRelativeTolerance(1e-5);
  solver.initialize();

  double y = -1.0;
  double y_out = y;
  double t = 0.0;
  double dt = 0.01;
  double maxt = 1.1;

  double error = 0.0;

  while(t + dt < maxt)
  {
    solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

    double y_anal = problem1(t + dt);

    double currError = (y_out - y_anal);
    error += currError * currError;

    t += dt;
    y = y_out;
  }

  error = sqrt(error);

  QVERIFY2( error < 1e-4 , QString("RKQS Problem 1 Error: %1").arg(error).toStdString().c_str());
}

#ifdef USE_CVODE

void STMComponentTest::solveODEAdam_Prob1()
{

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::CVODE_ADAMS);
    solver.setRelativeTolerance(1e-10);
    solver.setAbsoluteTolerance(1e-16);
    solver.setOrder(12);
    solver.initialize();

    double y = -1.0;
    double y_out = y;
    double t = 0.0;
    double dt = 0.01;
    double maxt = 1.1;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

      double y_anal = problem1(t + dt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("CVODE Adams Problem 1 Error: %1").arg(error).toStdString().c_str());
  }
}

void STMComponentTest::solveODEBDF_Prob1()
{

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::CVODE_BDF);
    solver.setRelativeTolerance(1e-10);
    solver.setAbsoluteTolerance(1e-16);
    solver.setOrder(5);
    solver.initialize();

    double y = -1.0;
    double y_out = y;
    double t = 0.0;
    double dt = 0.01;
    double maxt = 1.1;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

      double y_anal = problem1(t + dt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("CVODE BDF Problem 1 Error: %1").arg(error).toStdString().c_str());
  }
}

#endif

void STMComponentTest::solveODERK4_Prob2()
{
  //  printf("RK4\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RK4);
    solver.initialize();

    double y = 3.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 5.0;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = problem2(tdt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("RK4 Problem 2 Error: %1").arg(error).toStdString().c_str());
  }
}

void STMComponentTest::solveODERKQS_Prob2()
{
  //  printf("RKQS\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RKQS);
    solver.setRelativeTolerance(1e-3);
    solver.initialize();

    double y = 3.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 5.0;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = problem2(tdt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("RKQS Problem 2 Error: %1").arg(error).toStdString().c_str());
  }
}

#ifdef USE_CVODE

void STMComponentTest::solveODEAdam_Prob2()
{
  //  printf("RKQS\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::CVODE_ADAMS);
    solver.setRelativeTolerance(1e-10);
    solver.setAbsoluteTolerance(1e-16);
    solver.setOrder(10);
    solver.initialize();

    double y = 3.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 5.0;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = problem2(tdt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("CVODE ADAMS Problem 2 Error: %1").arg(error).toStdString().c_str());
  }
}

void STMComponentTest::solveODEBDF_Prob2()
{
  //  printf("RKQS\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::CVODE_BDF);
    solver.setRelativeTolerance(1e-10);
    solver.setAbsoluteTolerance(1e-16);
    solver.setOrder(5);
    solver.initialize();

    double y = 3.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 5.0;

    double error = 0.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = problem2(tdt);

      double currError = (y_out - y_anal);
      error += currError * currError;

      t += dt;
      y = y_out;
    }

    error = sqrt(error);

    QVERIFY2( error < 1e-4 , QString("CVODE BDF Problem 2 Error: %1").arg(error).toStdString().c_str());
  }
}

#endif

void STMComponentTest::versteegCase1_Upwind()
{
QBENCHMARK_ONCE
  {

    //Error messages
    std::list<std::string> errors;

    //Stream temperature model instance
    STMModel *model = new STMModel(nullptr);

    //Set advection discretization mode
    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Upwind);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);//1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case1/versteegcase1_upwind.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //Domain spatial discretization

    //Element junction coordinates
    double x , y = 0, z = 0;

    int numCells = 5.0;

    //grid size in x direction
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Channels
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 0.1;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

    //initialize model
    if(model->initialize(errors))
    {


      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //compare results

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void STMComponentTest::versteegCase1_Central()
{
QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    STMModel *model = new STMModel(nullptr);

    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Central);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);//1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case1/versteegcase1_central.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //  Domain spatial discretization
    double x , y = 0, z = 0;
    int numCells = 5.0;
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Channels
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 0.1;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

    //initialize model
    if(model->initialize(errors))
    {



      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //compare results

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void STMComponentTest::versteegCase1_Hybrid()
{
QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    //Stream temperature model instance
    STMModel *model = new STMModel(nullptr);

    //Set advection discretization mode
    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Hybrid);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);//1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case1/versteegcase1_hybrid.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //  Domain spatial discretization
    double x , y = 0, z = 0;
    int numCells = 5.0;
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Channels
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 0.1;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

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

void STMComponentTest::versteegCase2_Upwind()
{
  QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    //Stream temperature model instance
    STMModel *model = new STMModel(nullptr);

    //Set advection discretization mode
    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Upwind);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);// Modified Julian DateTime 1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_upwind.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::RKQS);

    //  Domain spatial discretization
    double x , y = 0, z = 0;
    int numCells = 5.0;
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Channels
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 2.5;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

    //initialize model
    if(model->initialize(errors))
    {


      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //compare results

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void STMComponentTest::versteegCase2_Central()
{
  QBENCHMARK_ONCE
  {
    std::list<std::string> errors;

    STMModel *model = new STMModel(nullptr);

    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Central);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);//1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_central.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //  Domain spatial discretization
    double x , y = 0, z = 0;
    int numCells = 5.0;
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Channels
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 2.5;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

    //initialize model
    if(model->initialize(errors))
    {
      //Perform timestep until completion
      while (model->currentDateTime() < model->endDateTime())
      {
        model->update();
      }
    }

    //compare results

    //finalize model
    model->finalize(errors);

    delete model;
  }
}

void STMComponentTest::versteegCase2_Hybrid()
{
  QBENCHMARK_ONCE
  {

    //Error messages
    std::list<std::string> errors;

    //Stream temperature model instance
    STMModel *model = new STMModel(nullptr);

    //Set advection discretization mode
    model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Hybrid);

    //set time and output variables
    model->setMaxTimeStep(0.5); //seconds
    model->setMinTimeStep(0.0001); //seconds
    model->setStartDateTime(0); // Modified Julian DateTime
    model->setEndDateTime(1.0 / 1440.0);//1 minute
    model->setOutputInterval(1.0); //1 second

    //Output consider 1 date separate files easier. Cross tab
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_hybrid.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //Domain spatial discretization

    //Element junction coordinates
    double x = 0, y = 0, z = 0;
    int numCells = 5.0;

    //grid size in x direction
    double dx = 1.0 / numCells;

    //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
    for(int i = 0; i <= numCells; i++)
    {
      //Make it explicit what each variable means
      x = i * dx;
      ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

      //Apply boundary condition at the inlet junction
      if(i == 0)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 1.0;
      }
      //Apply boundary condition at the outlet junction
      else if(i == numCells)
      {
        eJunction->temperature.isBC = true;
        eJunction->temperature.value = 0.0;
      }
    }

    //Discretize elements using upstream and downstream junctions created earlier
    for(int i = 0; i < numCells ; i++)
    {
      Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
      element->xSectionArea = 0.1;
      element->length = dx;
      element->flow = element->xSectionArea * 2.5;
      element->depth = 0.1;
      element->longDispersion.isBC = true;
      element->longDispersion.value = 0.1;
    }

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

void STMComponentTest::derivativeProb1(double t, double y[], double dydt[], void *userData)
{
  dydt[0] = t * pow(y[0],3) / sqrt(1 + t * t);
}

double STMComponentTest::problem1(double t)
{
  return -1.0 / sqrt(3.0 - 2.0 * sqrt(1+t*t));
}

void STMComponentTest::derivativeProb2(double t, double y[], double dydt[], void* userData)
{
  dydt[0] = (3 * t*t + 4 * t - 4)/(2 * y[0] - 4);
}

double STMComponentTest::problem2(double t)
{
  return 2.0 + sqrt(t*t*t + 2.0*t*t - 4.0*t + 2.0);
}
