/*!
*  \file    stmproject.h
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

#ifndef STMMODEL_H
#define STMMODEL_H

#include "stmcomponent_global.h"
#include "spatial/network.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <QFileInfo>
#include <QTextStream>

class STMComponent;
struct Element;
struct ElementJunction;
class Edge;
class ODESolver;
class STMModel;

struct SolverUserData
{
    STMModel *model = nullptr;
    int variableIndex = -1;
};

class STMCOMPONENT_EXPORT STMModel : public QObject
{
    Q_OBJECT

    friend struct ElementJunction;
    friend struct Element;

  public:

    enum AdvectionDiscretizationMode
    {
      Upwind,
      Central,
      Hybrid
    };

    /*!
     * \brief STMModel - Computational engine for the Stream Temperature Model
     * \param component - Paremt model coupling component if being used for coupling.
     */
    STMModel(STMComponent *component);

    ~STMModel();

    /*!
     * \brief minTimeStep
     * \return
     */
    double minTimeStep() const;

    /*!
     * \brief setMinTimeStep
     * \param timeStep
     */
    void setMinTimeStep(double timeStep);

    /*!
     * \brief maxTimeStep
     * \return
     */
    double maxTimeStep() const;

    /*!
     * \brief setMaxTimeStep
     * \param timeStep
     */
    void setMaxTimeStep(double timeStep);

    /*!
     * \brief useAdaptiveTimeStep
     * \return
     */
    bool useAdaptiveTimeStep() const;

    /*!
     * \brief setUseAdaptiveTimeStep
     * \param use
     */
    void setUseAdaptiveTimeStep(bool use);

    /*!
     * \brief currentTimeStep
     * \return
     */
    double currentTimeStep() const;

    /*!
     * \brief startDateTime
     * \return
     */
    double startDateTime() const;

    /*!
     * \brief setStartDate
     * \param dateTime
     */
    void setStartDateTime(double dateTime);

    /*!
     * \brief endDateTime
     * \return
     */
    double endDateTime() const;

    /*!
     * \brief setEndDateTime
     * \param dateTime
     */
    void setEndDateTime(double dateTime);

    /*!
     * \brief currentDateTime
     * \return
     */
    double currentDateTime() const;

    /*!
     * \brief outputInterval
     * \return
     */
    double outputInterval() const;

    /*!
     * \brief setOutputInterval
     * \param interval
     */
    void setOutputInterval(double interval);

    /*!
     * \brief computeLongDispersion
     * \return
     */
    bool computeLongDispersion() const;

    /*!
     * \brief setCalculateLongDispersion
     * \param calculate
     */
    void setComputeLongDispersion(bool calculate);

    /*!
     * \brief advectionDiscretizationMode
     * \return
     */
    AdvectionDiscretizationMode advectionDiscretizationMode() const;

    /*!
     * \brief setAdvectionDiscretizationMode
     * \param advectionDiscretizationMode
     */
    void setAdvectionDiscretizationMode(AdvectionDiscretizationMode advectionDiscretizationMode);

    /*!
     * \brief numSolutes
     * \return
     */
    int numSolutes() const;

    /*!
     * \brief addSolute
     * \param soluteName
     * \return
     */
    bool addSolute(const std::string &soluteName);

    /*!
     * \brief removeSolute
     * \param soluteName
     */
    void removeSolute(const std::string &soluteName);

    /*!
     * \brief solutes
     * \return
     */
    std::vector<std::string> solutes() const;

    /*!
     * \brief numElementJunctions
     * \return
     */
    int numElementJunctions() const;

    /*!
     * \brief addControlVolumeNode
     * \param id
     * \param x
     * \param y
     * \param z
     * \return
     */
    ElementJunction *addElementJunction(const std::string &id, double x = 0, double y = 0, double z = 0);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(const std::string &id);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(int id);

    /*!
     * \brief getElementJunction
     * \param id
     * \return
     */
    ElementJunction *getElementJunction(const std::string &id) ;

    /*!
     * \brief getElementJunction
     * \param index
     * \return
     */
    ElementJunction *getElementJunction(int index) ;

    /*!
     * \brief numElements
     * \return
     */
    int numElements() const;

    /*!
     * \brief addElement
     * \param id
     * \param fromElement
     * \param toElement
     * \return
     */
    Element *addElement(const std::string &id, ElementJunction *upStream, ElementJunction *downStream);

    /*!
     * \brief removeElement
     * \param id
     */
    void deleteElement(const std::string &id);

    /*!
     * \brief removeElement
     * \param index
     */
    void deleteElement(int index);

    /*!
     * \brief getElement
     * \param id
     * \return
     */
    Element *getElement(const std::string &id);

    /*!
     * \brief getElement
     * \param index
     * \return
     */
    Element *getElement(int index);

    /*!
     * \brief outputCSVFile
     * \return
     */
    QFileInfo outputCSVFile() const;

    /*!
     * \brief setOutputCSVFile
     * \param outputFile
     */
    void setOutputCSVFile(const QFileInfo &outputFile);

    /*!
     * \brief initialize
     * \param errors
     * \return
     */
    bool initialize(std::list<std::string> &errors);

    /*!
     * \brief update
     */
    void update();

    /*!
     * \brief finalize
     * \param errors
     * \return
     */
    bool finalize(std::list<std::string> &errors);

  private:

    /*!
     * \brief initializeInputFiles
     * \param errors
     * \return
     */
    bool initializeInputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeTimeVariables
     * \param errors
     * \return
     */
    bool initializeTimeVariables(std::list<std::string> &errors);

    /*!
     * \brief initializeElements
     * \param errors
     * \return
     */
    bool initializeElements(std::list<std::string> &errors);

    /*!
     * \brief initializeCSVOutputFile
     * \param errors
     * \return
     */
    bool initializeCSVOutputFile(std::list<std::string> &errors);

    /*!
     * \brief initializeNetCDFOutputFile
     * \param errors
     * \return
     */
    bool initializeNetCDFOutputFile(std::list<std::string> &errors);

    /*!
     * \brief prepareForNextTimeStep
     */
    void prepareForNextTimeStep();

    /*!
     * \brief applyInitialConditions
     */
    void applyInitialConditions();

    /*!
     * \brief applyBoundaryConditions
     */
    void applyBoundaryConditions(double dateTime);

    /*!
     * \brief computeTimeStep
     * \return
     */
    double computeTimeStep() const;

    /*!
     * \brief computeLongDispersion
     */
    void computeLongDispersion();

    /*!
     * \brief solveJunctionHeatContinuity
     * \param timeStep
     */
    void solveJunctionHeatContinuity(double timeStep, bool &converged);

    /*!
     * \brief solveJunctionSoluteContinuity
     * \param soluteIndex
     * \param timeStep
     */
    void solveJunctionSoluteContinuity(int soluteIndex, double timeStep, bool &converged);

    /*!
     * \brief computeJunctionDYDt
     * \param model
     * \param variableIndex
     * \param t
     * \param y
     * \param dydt
     */
    static void computeJunctionDYDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief solveElementHeatTransport
     * \param timeStep
     */
    void solveElementHeatTransport(double timeStep);

    /*!
     * \brief solveElementSoluteTransport
     * \param soluteIndex
     * \param timeStep
     */
    void solveElementSoluteTransport(int soluteIndex, double timeStep);

    /*!
     * \brief computeElementDYDt
     * \param model
     * \param variableIndex
     * \param t
     * \param y
     * \param dydt
     */
    static void computeElementDYDt(double t, double y[], double dydt[], void* userData);

    /*!
     * \brief writeOutput
     */
    void writeOutput();

    /*!
     * \brief writeCSVOutput
     */
    void writeCSVOutput();

    /*!
     * \brief writeNetCDFOutput
     */
    void writeNetCDFOutput();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputFiles();

    /*!
     * \brief closeCSVOutputFile
     */
    void closeCSVOutputFile();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputNetCDFFile();

  private:

    //Solute variable names;
    std::vector<std::string> m_solutes; // Names of the solutes.

    //Time variables
    double m_timeStep = 0.001, //seconds
           m_startDateTime, //MJD
           m_endDateTime, //MJD
           m_currentDateTime, //MJD
           m_maxTimeStep = 0.5, //seconds
           m_minTimeStep = 0.001, //seconds
           m_outputInterval, //seconds
           m_nextOutputTime,//MJD
           m_timeStepRelaxationFactor = 0.9; //


    bool m_computeDispersion = false,
         m_useAdaptiveTimeStep = true,
         m_converged = true;

    /*!
     * \brief m_advectionMode
     */
    AdvectionDiscretizationMode m_advectionMode = AdvectionDiscretizationMode::Upwind;

    //Element junctions
    std::vector<ElementJunction*> m_elementJunctions;
    std::unordered_map<std::string, ElementJunction*> m_elementJunctionsById; //added for fast lookup using identifiers instead of indexes.

    //1D Computational elements
    std::vector<Element*> m_elements;
    std::unordered_map<std::string, Element*> m_elementsById; //added for fast lookup using identifiers instead of indexes.

    //Number of junctions where continuity needs to be enforced.
    int m_numHeatContinuityElementJunctions = 0;
    std::vector<int> m_numSoluteContinuityElementJunctions;

    //Convergence tolerance
    double m_convergenceTol = 1e-5;
    int m_maxNumberOfIterations =  2;

    //File input and output
    QFileInfo m_outputCSVFileInfo;
    QTextStream m_outputCSVStream;
    STMComponent *m_component;

    //Global water properties
    double m_waterDensity = 1.0; //kg/m^3
    double m_cp = 1.0;// 4187.0; // J/kg/C
};

#endif // STMMODEL_H
