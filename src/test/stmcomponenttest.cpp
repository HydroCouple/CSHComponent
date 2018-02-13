#include "stdafx.h"
#include "test/stmcomponenttest.h"
#include "odesolver.h"
#include "stmmodel.h"
#include "elementjunction.h"
#include "element.h"
#include "variable.h"

void STMComponentTest::solveODERK4_Prob1()
{
  //  printf("RK4\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RK4);

    double y = 1.0 / 25.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 3.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

      double y_anal = 1.0 / (28.0 - 3 * (t+dt) * (t+dt));

      QVERIFY2( fabs(y_out - y_anal) < 1e-4 , QString("Time: %1 Y_out: %2 Y_anal: %3 Error: %4").arg(t+dt).arg(y_out).arg(y_anal).arg(fabs(y_anal - y_out)).toStdString().c_str());

      t += dt;
      y = y_out;
    }
  }
}


void STMComponentTest::solveODERKQS_Prob1()
{
  //  printf("RKQS\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RKQS);

    double y = 1.0 / 25.0;
    double y_out = y;
    double t = 1.0;
    double dt = 0.01;
    double maxt = 3.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb1, nullptr);

      double y_anal = 1.0 / (28.0 - 3 * (t+dt) * (t+dt));

      QVERIFY2( fabs(y_out - y_anal) < 1e-4 , QString("Time: %1 Y_out: %2 Y_anal: %3 Error: %4").arg(t+dt).arg(y_out).arg(y_anal).arg(fabs(y_anal - y_out)).toStdString().c_str());

      t += dt;
      y = y_out;
    }
  }
}

void STMComponentTest::solveODERK4_Prob2()
{
  //  printf("RK4\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RK4);

    double y = 1.0;
    double y_out = y;
    double t = 0.0;
    double dt = 0.01;
    double maxt = 5.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = -(tdt * tdt * y / (2.0 * sqrt(2.0 - y * y))) + 1.0;
      QVERIFY2( fabs(y_out - y_anal) < 1e-4 , QString("Time: %1 Y_out: %2 Y_anal: %3 Error: %4").arg(t+dt).arg(y_out).arg(y_anal).arg(fabs(y_anal - y_out)).toStdString().c_str());

      t += dt;
      y = y_out;
    }
  }
}


void STMComponentTest::solveODERKQS_Prob2()
{
  //  printf("RKQS\n");

  QBENCHMARK
  {
    ODESolver solver(1, ODESolver::RKQS);

    double y = 1.0;
    double y_out = y;
    double t = 0.0;
    double dt = 0.01;
    double maxt = 5.0;

    while(t + dt < maxt)
    {
      solver.solve(&y, 1, t, dt, &y_out, &STMComponentTest::derivativeProb2, nullptr);

      double tdt = t + dt;
      double y_anal = -(tdt * tdt * y / (2.0 * sqrt(2.0 - y * y))) + 1.0;
      QVERIFY2( fabs(y_out - y_anal) < 1e-4 , QString("Time: %1 Y_out: %2 Y_anal: %3 Error: %4").arg(t+dt).arg(y_out).arg(y_anal).arg(fabs(y_anal - y_out)).toStdString().c_str());

      t += dt;
      y = y_out;
    }
  }
}


void STMComponentTest::versteegCase1_Upwind()
{
  std::list<std::string> errors;

  STMModel *model = new STMModel(nullptr);

  model->setAdvectionDiscretizationMode(STMModel::AdvectionDiscretizationMode::Upwind);

  //set time and output variables
  model->setMaxTimeStep(0.01); //seconds
  model->setMinTimeStep(0.0025); //seconds
  model->setStartDateTime(0); // Modified Julian DateTime
  model->setEndDateTime(0.0416667);//1 hour
  model->setOutputInterval(60); //1 minutes

  //Output consider 1 date separate files easier. Cross tab
  model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case1/versteegcase1.csv"));

  //  Domain spatial discretization
  double x , y = 0, z = 0;
  double dx = 1.0 / 40;

  //Channel Junctions con. Nice diagram illustrating the control volumes. Junctions and elements. Consider changing.
  for(int i = 0; i <= 40; i++)
  {
    //Make it explicit what each variable means
    x = i * dx;
    ElementJunction *eJunction = model->addElementJunction(QString::number(i).toStdString(), x, y, z);

    if(i == 0)
    {
      eJunction->temperature.isBC = true;
      eJunction->temperature.value = 1.0;
    }
    else if(i == 40)
    {
      eJunction->temperature.isBC = true;
      eJunction->temperature.value = 0.0;
    }
  }

  //Channels
  for(int i = 0; i < 40 ; i++)
  {
    Element *element = model->addElement(QString::number(i).toStdString(), model->getElementJunction(i), model->getElementJunction(i+1));
    element->xSectionArea = 0.1;
    element->length = dx;
    element->flow = element->xSectionArea * 2.5;
    element->longDispersion.isBC = true;
    element->longDispersion.value = 0.1;
  }

  //initialize model
  model->initialize(errors);

  //Perform timestep until completion
  while (model->currentDateTime() < model->endDateTime())
  {
    model->update();
  }

  //compare results

  //finalize model
  model->finalize(errors);

  delete model;
}

void STMComponentTest::versteegCase2_Upwind()
{

}

void STMComponentTest::versteegCase1_Central()
{

}

void STMComponentTest::versteegCase2_Central()
{

}

void STMComponentTest::versteegCase1_Hybrid()
{

}

void STMComponentTest::versteegCase2_Hybrid()
{

}

void STMComponentTest::derivativeProb1(double t, double y[], double dydt[], void *userData)
{
  dydt[0] = 6.0 * y[0] * y[0] * t;
}

void STMComponentTest::derivativeProb2(double t, double y[], double dydt[], void* userData)
{
  dydt[0] = -t * y[0] / sqrt( 2.0 - y[0] * y[0]);
}
