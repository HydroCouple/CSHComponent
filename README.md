# STMComponent
[![AUR](https://img.shields.io/aur/license/yaourt.svg)](https://github.com/HydroCouple/STMComponent/blob/master/LICENSE)

A one-dimensional stream heat and solute transport model.

## Input File Specification
---------------------------------------

The STMComponent requires two input files. The first input file defines the spatial descritzation for the model and is defined using the SWMM input file format. The second input file is used to specify the hydrodynamic parameters, initial conditions, boundary conditions, and time varying hydraulic inputs. The conventions for the second file are similar to the SWMM file conventions and are as follows:

```
[OPTIONS]
ADVECTION_MODE UPWIND/CENTRAL/HYBRID/ULTIMATE
START_DATETIME month/day/year hour/minute/second
END_DATETIME month/day/year hour/minute/second
REPORT_INTERVAL seconds
MAX_TIME_STEP seconds
MIN_TIME_STEP seconds
NUM_INITIAL_FIXED_STEPS value
USE_ADAPTIVE_TIME_STEP YES/NO
TIME_STEP_RELAXATION_FACTOR value
TEMP_TRANSPORT YES/NO
WATER_DENSITY value
WATER_SPECIFIC_HEAT_CAPACITY value
NUM_SOLUTES value
VERBOSE YES/NO

[OUTPUTS]
CSV "CSV output file path"
netCDF "NetCDF output file path"

[SOLUTES]
;;SOLUTE_NAME
;;=================================================

[INITIAL_CONDITIONS] <br/>
;;CONDUIT    FLOW    XSECTION_AREA    DEPTH    WIDTH   TEMPERATURE    SOLUTE1    SOLUTE2
;;=======================================================================================

[BOUNDARY_CONDITIONS]
;;JUNCTION  TEMPERATURE SOLUTE1 SOLUTE2
;;=====================================

[POINT_SOURCES]
;;CONDUIT  HEAT_FLUX SOLUTE1_FLUX SOLUTE2_FLUX
;;============================================


[NON_POINT_SOURCES]
;;START_CONDUIT  END_CONDUIT  HEAT_FLUX_PER_UNIT_LENGTH SOLUTE1_FLUX_PER_UNIT_LENGTH
;;==================================================================================

[TIME_VARYING_HYDRAULICS]
csvfile path
```
As with SWMM all lines beginning with ";;" are comment lines and will not be read. The time varying hydraulics input file is optional and can be provided in a csv format with the following conventions:

```
Conduit, DateTime           , Flow, Cross-sectional Area, Depth, Width
C1     , 02/02/2018 02/02/45, 2   , 0.1                 , 0.1  , 0.1
```

## Sample Code
---------------------------------------
``` C++


    /*!
     * Case 2 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the upwind differencing method for advection discretization.
     */

    //Error Messages
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

    //Output
    model->setOutputCSVFile(QFileInfo("../../examples/Versteeg/case2/versteegcase2_upwind.csv"));

    //Set solver type
    model->solver()->setSolverType(ODESolver::CVODE_ADAMS);

    //  Domain spatial discretization
    double x , y = 0, z = 0;
    int numCells = 5.0;
    double dx = 1.0 / numCells;

    //Channel Junctions connections
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

    //finalize model
    model->finalize(errors);
    delete model;

```