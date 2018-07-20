/*!
*  \file    stmmodelio.cpp
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

#include "stmmodel.h"
#include "element.h"
#include "elementjunction.h"
#include "junctiontimeseriesbc.h"
#include "radiativefluxtimeseriesbc.h"
#include "hydraulicstimeseriesbc.h"
#include "pointsrctimeseriesbc.h"
#include "nonpointsrctimeseriesbc.h"
#include "meteorologytimeseriesbc.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "timeseries.h"

#include <QDir>
#include <QDate>
#include <cstdlib>
#include <errno.h>

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;


#endif

using namespace std;

bool STMModel::verbose() const
{
  return m_verbose;
}

void STMModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int STMModel::printFrequency() const
{
  return m_printFrequency;
}

void STMModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int STMModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void STMModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo STMModel::inputFile() const
{
  return m_inputFile;
}

void STMModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo STMModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void STMModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

QFileInfo STMModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void STMModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void STMModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("STM TimeStep (s): %f\tDateTime: %f\tTemp (°C) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalHeatBalance: %g (KJ)}", m_timeStep, m_currentDateTime,
           m_heatSolver->getIterations(), m_heatSolver->maxIterations(), m_minTemp, m_maxTemp, m_totalHeatBalance);

    for (size_t j = 0; j < m_solutes.size(); j++)
    {
      std::string &solute = m_solutes[j];
      ODESolver *solver = m_soluteSolvers[j];
      printf("\t%s (kg/m) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalMassBalance: %g (kg)}", solute.c_str(), solver->getIterations(), solver->maxIterations(),
             m_minSolute[j], m_maxSolute[j], m_totalSoluteMassBalance[j]);
    }

    printf("\n");

    m_currentPrintCount = 0;
  }
}

void STMModel::saveAs(const QFileInfo &filePath)
{
  QFileInfo fileInfo;

  if (filePath.isRelative())
  {
    fileInfo = relativePathToAbsolute(filePath);
  }
  else
  {
    fileInfo = filePath;
  }

  QString file = fileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && filePath.absoluteDir().exists() && fileInfo.isFile())
  {


  }
}

bool STMModel::initializeInputFiles(list<string> &errors)
{

  if (QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if (file.open(QIODevice::ReadOnly))
    {
      m_delimiters = QRegExp("(\\,|\\t|\\;|\\s+)");
      int currentFlag = -1;
      m_addedSoluteCount = 0;

      QTextStream streamReader(&file);
      int lineCount = 0;

      while (!streamReader.atEnd())
      {
        QString line = streamReader.readLine().trimmed();
        lineCount++;

        if (!line.isEmpty() && !line.isNull())
        {
          bool readSuccess = true;
          QString error = "";

          auto it = m_inputFileFlags.find(line.toStdString());

          if (it != m_inputFileFlags.cend())
          {
            currentFlag = it->second;
          }
          else if (!QStringRef::compare(QStringRef(&line, 0, 2), ";;"))
          {
            //commment do nothing
          }
          else
          {
            switch (currentFlag)
            {
              case 1:
                readSuccess = readInputFileOptionTag(line, error);
                break;
              case 2:
                readSuccess = readInputFileOutputTag(line, error);
                break;
              case 3:
                readSuccess = readInputFileSolutesTag(line, error);
                break;
              case 4:
                readSuccess = readInputFileElementJunctionsTag(line, error);
                break;
              case 5:
                readSuccess = readInputFileElementsTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileBoundaryConditionsTag(line, error);
                break;
              case 7:
                readSuccess = readInputFilePointSourcesTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileNonPointSourcesTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileUniformHydraulicsTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileNonUniformHydraulicsTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileUniformRadiativeFluxesTag(line, error);
                break;
              case 12:
                readSuccess = readInputFileNonUniformRadiativeFluxesTag(line, error);
                break;
              case 13:
                readSuccess = readInputFileUniformMeteorologyTag(line, error);
                break;
              case 14:
                readSuccess = readInputFileNonUniformMeteorologyTag(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
            break;
          }
        }
      }

      file.close();
    }
  }

  return true;
}

bool STMModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool STMModel::initializeCSVOutputFile(list<string> &errors)
{

  if (m_outputCSVFileInfo.isRelative())
  {
    m_outputCSVFileInfo = relativePathToAbsolute(m_outputCSVFileInfo);
  }

  QString file = m_outputCSVFileInfo.absoluteFilePath();

  if (!file.isEmpty() && !file.isNull() && !m_outputCSVFileInfo.absoluteDir().exists())
  {
    errors.push_back("Output shapefile directory does not exist: " + file.toStdString());
    return false;
  }

  if (!file.isEmpty() && !file.isNull())
  {
    if(m_outputCSVFileInfo.isDir())
      return true;

    if (m_outputCSVStream.device() == nullptr)
    {
      QFile *device = new QFile(file, this);
      m_outputCSVStream.setDevice(device);
    }

    if (m_outputCSVStream.device()->open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_outputCSVStream.setRealNumberPrecision(10);
      m_outputCSVStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_outputCSVStream << "DateTime, ElementId, ElementIndex, x, y, z, Depth, Width, XSectionArea, Dispersion, "
                           "Temperature, TotalAdvDispHeatBalance, TotalExternalHeatFluxesBalance, TotalRadFluxesHeatBalance,"
                           "TotalEvapFluxHeatBalance, TotalConvFluxHeatBalance, TotalHeatBalance";

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        QString soluteName = QString::fromStdString(m_solutes[i]);
        m_outputCSVStream << ", " << soluteName << ", "
                          << "TotalAdvDispMassBalance_" + soluteName << ", "
                          << "TotalExternalFluxMassBalance_" + soluteName << ","
                          << "TotalMassBalance_" + soluteName ;
      }

      m_outputCSVStream << endl;
    }

    m_outputCSVStream.flush();

    return true;
  }

  return false;
}

bool STMModel::initializeNetCDFOutputFile(list<string> &errors)
{

#ifdef USE_NETCDF


  if (m_outputNetCDFFileInfo.isRelative())
  {
    m_outputNetCDFFileInfo = relativePathToAbsolute(m_outputNetCDFFileInfo);
  }

  if (m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() || m_outputNetCDFFileInfo.isDir())
  {
    return true;
  }
  else if (!m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() &&
           !m_outputNetCDFFileInfo.absoluteFilePath().isNull() &&
           !m_outputNetCDFFileInfo.absoluteDir().exists())
  {
    std::string message = "NetCDF output file directory does not exist: " + m_outputNetCDFFileInfo.absoluteFilePath().toStdString();
    errors.push_back(message);
    return false;
  }

  bool returnValue = false;


  closeOutputNetCDFFile();

  try
  {
    m_outNetCDFVariables.clear();

    m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);

    //time variable
    ThreadSafeNcDim timeDim =  m_outputNetCDF->addDim("time");
    ThreadSafeNcVar timeVar =  m_outputNetCDF->addVar("time", NcType::nc_DOUBLE, timeDim);
    timeVar.putAtt("time:long_name", "time");
    timeVar.putAtt("time:units", "days since 1858-11-17 0:0:0");
    timeVar.putAtt("time:calendar", "modified_julian");
    m_outNetCDFVariables["time"] = timeVar;

    //Add Solutes
    ThreadSafeNcDim solutesDim =  m_outputNetCDF->addDim("solutes", m_solutes.size());
    ThreadSafeNcVar solutes =  m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
    solutes.putAtt("solute_names::long_name", "solute names");
    m_outNetCDFVariables["solutes"] = solutes;

    if (m_solutes.size())
    {
      char **soluteNames = new char *[m_solutes.size()];

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        delete[] soluteNames[i];
      }

      delete[] soluteNames;
    }

    //Add element junctions
    ThreadSafeNcDim junctionDim =  m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers =  m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("element_junction_id::long_name", "element junction identifier");
    m_outNetCDFVariables["element_junction_id"] = junctionIdentifiers;

    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_DOUBLE, junctionDim);
    junctionX.putAtt("x:long_name", "projection_x_coordinate");
    junctionX.putAtt("x:units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_DOUBLE, junctionDim);
    junctionY.putAtt("y:long_name", "projection_y_coordinate");
    junctionY.putAtt("y:units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_DOUBLE, junctionDim);
    junctionZ.putAtt("z:long_name", "projection_z_coordinate");
    junctionZ.putAtt("z:units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    double *vertx = new double[m_elementJunctions.size()];
    double *verty = new double[m_elementJunctions.size()];
    double *vertz = new double[m_elementJunctions.size()];
    char **junctionIds = new char *[m_elementJunctions.size()];

    //write other relevant junction attributes here.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elementJunctions.size(); i++)
    {
      ElementJunction *junction = m_elementJunctions[i];

      junctionIds[i] = new char[junction->id.size() + 1];
      strcpy(junctionIds[i], junction->id.c_str());

      vertx[i] = junction->x;
      verty[i] = junction->y;
      vertz[i] = junction->z;
    }

    junctionX.putVar(vertx);
    junctionY.putVar(verty);
    junctionZ.putVar(vertz);
    junctionIdentifiers.putVar(junctionIds);

    delete[] vertx;
    delete[] verty;
    delete[] vertz;

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      delete[] junctionIds[i];
    }

    delete[] junctionIds;

    //Add Elements
    ThreadSafeNcDim elementsDim =  m_outputNetCDF->addDim("elements", m_elements.size());

    ThreadSafeNcVar elementIdentifiers =  m_outputNetCDF->addVar("element_id", NcType::nc_STRING, elementsDim);
    elementIdentifiers.putAtt("element_id::long_name", "element identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;

    ThreadSafeNcVar elementFromJunction =  m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("from_junction:long_name", "upstream junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction =  m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("to_junction:long_name", "downstream junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;

    ThreadSafeNcVar element_x =  m_outputNetCDF->addVar("element_x", NcType::nc_DOUBLE, elementsDim);
    element_x.putAtt("element_x:standard_name", "projection_x_coordinate");
    element_x.putAtt("element_x:long_name", "projection_x_coordinate");
    element_x.putAtt("element_x:units", "m");
    m_outNetCDFVariables["element_x"] = element_x;

    ThreadSafeNcVar element_y =  m_outputNetCDF->addVar("element_y", NcType::nc_DOUBLE, elementsDim);
    element_y.putAtt("element_y:standard_name", "projection_y_coordinate");
    element_y.putAtt("element_y:long_name", "projection_y_coordinate");
    element_y.putAtt("element_y:units", "m");
    m_outNetCDFVariables["element_y"] = element_y;

    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    double *elX = new double[m_elements.size()];
    double *elY = new double[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->index;
      toJunctions[i] = element->downstreamJunction->index;

      elX[i] = element->x;
      elY[i] = element->y;
    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    element_x.putVar(elX);
    element_y.putVar(elY);

    delete[] fromJunctions;
    delete[] toJunctions;
    delete[] elX;
    delete[] elY;

    for (size_t i = 0; i < m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;

    //hydraulics variables
    ThreadSafeNcVar flowVar =  m_outputNetCDF->addVar("flow", "double",
                                                      std::vector<std::string>({"time", "elements"}));
    flowVar.putAtt("flow:long_name", "flow");
    flowVar.putAtt("flow:units", "m^3/s");
    m_outNetCDFVariables["flow"] = flowVar;

    ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    depthVar.putAtt("depth:long_name", "channel flow depth");
    depthVar.putAtt("depth:units", "m");
    m_outNetCDFVariables["depth"] = depthVar;

    ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "double",
                                                       std::vector<std::string>({"time", "elements"}));
    widthVar.putAtt("width:long_name", "channel flow top width");
    widthVar.putAtt("width:units", "m");
    m_outNetCDFVariables["width"] = widthVar;

    ThreadSafeNcVar xsectAreaVar =  m_outputNetCDF->addVar("xsection_area", "double",
                                                           std::vector<std::string>({"time", "elements"}));
    xsectAreaVar.putAtt("xsection_area:long_name", "flow cross-section area");
    xsectAreaVar.putAtt("xsection_area:units", "m^2");
    m_outNetCDFVariables["xsection_area"] = xsectAreaVar;

    ThreadSafeNcVar dispersionVar =  m_outputNetCDF->addVar("dispersion", "double",
                                                            std::vector<std::string>({"time", "elements"}));
    dispersionVar.putAtt("dispersion:long_name", "longitudinal dispersion");
    dispersionVar.putAtt("dispersion:units", "m^2/s");
    m_outNetCDFVariables["dispersion"] = dispersionVar;

    ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("temperature", "double",
                                                             std::vector<std::string>({"time", "elements"}));
    temperatureVar.putAtt("temperature:long_name", "temperature");
    temperatureVar.putAtt("temperature:units", "°C");
    m_outNetCDFVariables["temperature"] = temperatureVar;


    ThreadSafeNcVar totalElementHeatBalanceVar =  m_outputNetCDF->addVar("total_element_heat_balance", "double",
                                                                         std::vector<std::string>({"time", "elements"}));
    totalElementHeatBalanceVar.putAtt("total_element_heat_balance:long_name", "total heat balance");
    totalElementHeatBalanceVar.putAtt("total_element_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_element_heat_balance"] = totalElementHeatBalanceVar;

    ThreadSafeNcVar totalElementAdvDispHeatBalanceVar =  m_outputNetCDF->addVar("total_element_adv_disp_heat_balance", "double",
                                                                                std::vector<std::string>({"time", "elements"}));
    totalElementAdvDispHeatBalanceVar.putAtt("total_element_adv_disp_heat_balance:long_name", "total element advection dispersion heat balance");
    totalElementAdvDispHeatBalanceVar.putAtt("total_element_adv_disp_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_element_adv_disp_heat_balance"] = totalElementAdvDispHeatBalanceVar;

    ThreadSafeNcVar totalElementEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_element_evap_heat_balance", "double",
                                                                             std::vector<std::string>({"time", "elements"}));
    totalElementEvapHeatBalanceVar.putAtt("total_element_evap_heat_balance:long_name", "total element evaporation heat balance");
    totalElementEvapHeatBalanceVar.putAtt("total_element_evap_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_element_evap_heat_balance"] = totalElementEvapHeatBalanceVar;

    ThreadSafeNcVar totalElementConvHeatBalanceVar =  m_outputNetCDF->addVar("total_element_conv_heat_balance", "double",
                                                                             std::vector<std::string>({"time", "elements"}));
    totalElementConvHeatBalanceVar.putAtt("total_element_conv_heat_balance:long_name", "total element convection heat balance");
    totalElementConvHeatBalanceVar.putAtt("total_element_conv_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_element_conv_heat_balance"] = totalElementConvHeatBalanceVar;


    ThreadSafeNcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_element_radiation_flux_heat_balance", "double",
                                                                                      std::vector<std::string>({"time", "elements"}));
    totalElementRadiationFluxHeatBalanceVar.putAtt("total_element_radiation_flux_heat_balance:long_name", "total element radiation flux heat balance");
    totalElementRadiationFluxHeatBalanceVar.putAtt("total_element_radiation_flux_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_element_radiation_flux_heat_balance"] = totalElementRadiationFluxHeatBalanceVar;

    ThreadSafeNcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_element_external_heat_flux_balance", "double",
                                                                                     std::vector<std::string>({"time", "elements"}));
    totalElementExternalHeatFluxBalanceVar.putAtt("total_element_external_heat_flux_balance:long_name", "total element external heat flux balance");
    totalElementExternalHeatFluxBalanceVar.putAtt("total_element_external_heat_flux_balance:units", "KJ");
    m_outNetCDFVariables["total_element_external_heat_flux_balance"] = totalElementExternalHeatFluxBalanceVar;

    ThreadSafeNcVar totalElementSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_element_solute_mass_balance", "double",
                                                                               std::vector<std::string>({"time", "solutes", "elements"}));
    totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:long_name", "total element solute mass balance");
    totalElementSoluteMassBalanceVar.putAtt("total_element_solute_mass_balance:units", "kg");
    m_outNetCDFVariables["total_element_solute_mass_balance"] = totalElementSoluteMassBalanceVar;

    ThreadSafeNcVar totalElementAdvDispSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_element_adv_disp_solute_mass_balance", "double",
                                                                                      std::vector<std::string>({"time", "solutes", "elements"}));
    totalElementAdvDispSoluteMassBalanceVar.putAtt("total_element_adv_disp_solute_mass_balance:long_name", "total element advection dispersion solute mass balance");
    totalElementAdvDispSoluteMassBalanceVar.putAtt("total_element_adv_disp_solute_mass_balance:units", "kg");
    m_outNetCDFVariables["total_element_adv_disp_solute_mass_balance"] = totalElementAdvDispSoluteMassBalanceVar;

    ThreadSafeNcVar totalElementExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_element_external_solute_flux_mass_balance", "double",
                                                                                           std::vector<std::string>({"time", "solutes"}));
    totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
    totalElementExternalSoluteFluxMassBalanceVar.putAtt("total_element_external_solute_flux_mass_balance:units", "kg");
    m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"] = totalElementExternalSoluteFluxMassBalanceVar;

    ThreadSafeNcVar elementEvapHeatFluxVar =  m_outputNetCDF->addVar("element_evap_heat_flux", "double",
                                                                     std::vector<std::string>({"time", "elements"}));

    elementEvapHeatFluxVar.putAtt("element_evap_heat_flux:long_name", "element evaporation heat flux");
    elementEvapHeatFluxVar.putAtt("element_evap_heat_flux:units", "W/m^2");
    m_outNetCDFVariables["element_evap_heat_flux"] = elementEvapHeatFluxVar;

    ThreadSafeNcVar elementConvHeatFluxVar =  m_outputNetCDF->addVar("element_conv_heat_flux", "double",
                                                                     std::vector<std::string>({"time", "elements"}));

    elementConvHeatFluxVar.putAtt("element_conv_heat_flux:long_name", "element convective heat flux");
    elementConvHeatFluxVar.putAtt("element_conv_heat_flux:units", "W/m^2");
    m_outNetCDFVariables["element_conv_heat_flux"] = elementConvHeatFluxVar;


    ThreadSafeNcVar elementRadiationFluxVar =  m_outputNetCDF->addVar("element_radiation_flux", "double",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementRadiationFluxVar.putAtt("element_radiation_flux:long_name", "element radiation flux");
    elementRadiationFluxVar.putAtt("element_radiation_flux:units", "W/m^2");
    m_outNetCDFVariables["element_radiation_flux"] = elementRadiationFluxVar;

    ThreadSafeNcVar elementHeatFluxVar =  m_outputNetCDF->addVar("element_heat_flux", "double",
                                                                 std::vector<std::string>({"time", "elements"}));

    elementHeatFluxVar.putAtt("element_heat_flux:long_name", "element heat flux");
    elementHeatFluxVar.putAtt("element_heat_flux:units", "J/s");
    m_outNetCDFVariables["element_heat_flux"] = elementHeatFluxVar;


    ThreadSafeNcVar elementAirTempVar =  m_outputNetCDF->addVar("element_air_temp", "double",
                                                                std::vector<std::string>({"time", "elements"}));

    elementAirTempVar.putAtt("element_air_temp:long_name", "Air temperature");
    elementAirTempVar.putAtt("element_air_temp:units", "C");
    m_outNetCDFVariables["element_air_temp"] = elementAirTempVar;


    ThreadSafeNcVar elementRHVar =  m_outputNetCDF->addVar("element_relative_humidity", "double",
                                                           std::vector<std::string>({"time", "elements"}));

    elementRHVar.putAtt("element_relative_humidity:long_name", "Relative Humidity");
    elementRHVar.putAtt("element_relative_humidity:units", "%");
    m_outNetCDFVariables["element_relative_humidity"] = elementRHVar;


    ThreadSafeNcVar elementWindSpeedVar =  m_outputNetCDF->addVar("element_wind_speed", "double",
                                                                  std::vector<std::string>({"time", "elements"}));

    elementWindSpeedVar.putAtt("element_wind_speed:long_name", "Relative Humidity");
    elementWindSpeedVar.putAtt("element_wind_speed:units", "m/s");
    m_outNetCDFVariables["element_wind_speed"] = elementWindSpeedVar;

    ThreadSafeNcVar elementVaporPressVar =  m_outputNetCDF->addVar("element_vapor_pressure", "double",
                                                                   std::vector<std::string>({"time", "elements"}));

    elementVaporPressVar.putAtt("element_vapor_pressure:long_name", " VaporPressure");
    elementVaporPressVar.putAtt("element_vapor_pressure:units", "kPa");
    m_outNetCDFVariables["element_vapor_pressure"] = elementVaporPressVar;

    ThreadSafeNcVar elementSatVaporPressVar =  m_outputNetCDF->addVar("element_saturated_vapor_pressure", "double",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementSatVaporPressVar.putAtt("element_saturated_vapor_pressure:long_name", "Saturated Vapor Pressure");
    elementSatVaporPressVar.putAtt("element_saturated_vapor_pressure:units", "kPa");
    m_outNetCDFVariables["element_saturated_vapor_pressure"] = elementSatVaporPressVar;


    ThreadSafeNcVar elementAirVaporPressVar =  m_outputNetCDF->addVar("element_air_vapor_pressure", "double",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementAirVaporPressVar.putAtt("element_air_vapor_pressure:long_name", " VaporPressure");
    elementAirVaporPressVar.putAtt("element_air_vapor_pressure:units", "kPa");
    m_outNetCDFVariables["element_air_vapor_pressure"] = elementAirVaporPressVar;

    ThreadSafeNcVar elementAirSatVaporPressVar =  m_outputNetCDF->addVar("element_air_saturated_vapor_pressure", "double",
                                                                         std::vector<std::string>({"time", "elements"}));

    elementAirSatVaporPressVar.putAtt("element_air_saturated_vapor_pressure:long_name", "Saturated Vapor Pressure");
    elementAirSatVaporPressVar.putAtt("element_air_saturated_vapor_pressure:units", "kPa");
    m_outNetCDFVariables["element_air_saturated_vapor_pressure"] = elementAirSatVaporPressVar;

    ThreadSafeNcVar totalHeatBalanceVar =  m_outputNetCDF->addVar("total_heat_balance", "double",
                                                                  std::vector<std::string>({"time"}));
    totalHeatBalanceVar.putAtt("total_heat_balance:long_name", "total heat balance");
    totalHeatBalanceVar.putAtt("total_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_heat_balance"] = totalHeatBalanceVar;

    ThreadSafeNcVar totalAdvDispHeatBalanceVar =  m_outputNetCDF->addVar("total_adv_disp_heat_balance", "double",
                                                                         std::vector<std::string>({"time"}));
    totalAdvDispHeatBalanceVar.putAtt("total_adv_disp_heat_balance:long_name", "total advection dispersion heat balance");
    totalAdvDispHeatBalanceVar.putAtt("total_adv_disp_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_adv_disp_heat_balance"] = totalAdvDispHeatBalanceVar;


    ThreadSafeNcVar totalEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_evap_heat_balance", "double",
                                                                      std::vector<std::string>({"time"}));
    totalEvapHeatBalanceVar.putAtt("total_evap_heat_balance:long_name", "total evaporation heat balance");
    totalEvapHeatBalanceVar.putAtt("total_evap_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_evap_heat_balance"] = totalEvapHeatBalanceVar;


    ThreadSafeNcVar totalConvHeatBalanceVar =  m_outputNetCDF->addVar("total_conv_heat_balance", "double",
                                                                      std::vector<std::string>({"time"}));
    totalConvHeatBalanceVar.putAtt("total_conv_heat_balance:long_name", "total convection heat balance");
    totalConvHeatBalanceVar.putAtt("total_conv_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_conv_heat_balance"] = totalConvHeatBalanceVar;


    ThreadSafeNcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_radiation_flux_heat_balance", "double",
                                                                               std::vector<std::string>({"time"}));
    totalRadiationFluxHeatBalanceVar.putAtt("total_radiation_flux_heat_balance:long_name", "total radiation flux heat balance");
    totalRadiationFluxHeatBalanceVar.putAtt("total_radiation_flux_heat_balance:units", "KJ");
    m_outNetCDFVariables["total_radiation_flux_heat_balance"] = totalRadiationFluxHeatBalanceVar;

    ThreadSafeNcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_external_heat_flux_balance", "double",
                                                                              std::vector<std::string>({"time"}));
    totalExternalHeatFluxBalanceVar.putAtt("total_external_heat_flux_balance:long_name", "total external heat flux balance");
    totalExternalHeatFluxBalanceVar.putAtt("total_external_heat_flux_balance:units", "KJ");
    m_outNetCDFVariables["total_external_heat_flux_balance"] = totalExternalHeatFluxBalanceVar;


    ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "double",
                                                         std::vector<std::string>({"time", "solutes", "elements"}));
    solutesVar.putAtt("solute_concentration:long_name", "solute concentrations");
    solutesVar.putAtt("solute_concentration:units", "kg/m^3");
    m_outNetCDFVariables["solute_concentration"] = solutesVar;

    ThreadSafeNcVar totalSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_solute_mass_balance", "double",
                                                                        std::vector<std::string>({"time", "solutes"}));
    totalSoluteMassBalanceVar.putAtt("total_solute_mass_balance:long_name", "total solute mass balance");
    totalSoluteMassBalanceVar.putAtt("total_solute_mass_balance:units", "kg");
    m_outNetCDFVariables["total_solute_mass_balance"] = totalSoluteMassBalanceVar;

    ThreadSafeNcVar totalAdvDispSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_adv_disp_solute_mass_balance", "double",
                                                                               std::vector<std::string>({"time", "solutes"}));
    totalAdvDispSoluteMassBalanceVar.putAtt("total_adv_disp_solute_mass_balance:long_name", "total advection dispersion solute mass balance");
    totalAdvDispSoluteMassBalanceVar.putAtt("total_adv_disp_solute_mass_balance:units", "kg");
    m_outNetCDFVariables["total_adv_disp_solute_mass_balance"] = totalAdvDispSoluteMassBalanceVar;

    ThreadSafeNcVar totalExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_external_solute_flux_mass_balance", "double",
                                                                                    std::vector<std::string>({"time", "solutes"}));
    totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:long_name", "total external solute flux mass balance");
    totalExternalSoluteFluxMassBalanceVar.putAtt("total_external_solute_flux_mass_balance:units", "kg");
    m_outNetCDFVariables["total_external_solute_flux_mass_balance"] = totalExternalSoluteFluxMassBalanceVar;

    m_outputNetCDF->sync();

    returnValue = true;

  }
  catch (NcException &e)
  {
    std::string message = std::string(e.what());
    printf("%s\n", e.what());
    errors.push_back(message);
    returnValue = false;
  }


#endif

  return returnValue;
}

bool STMModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  std::string optionsFlag = options[0].toStdString();
  auto it = m_optionsFlags.find(optionsFlag);

  if (it != m_optionsFlags.end())
  {
    int optionsIndex = it->second;

    switch (optionsIndex)
    {
      case 1:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_startDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 2:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_endDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 3:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_outputInterval = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Report interval error";
            return false;
          }
        }
        break;
      case 4:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_maxTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Max time step error";
            return false;
          }
        }
        break;
      case 5:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_minTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Min time step error";
            return false;
          }
        }
        break;
      case 6:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numInitFixedTimeSteps = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of initial time step error";
            return false;
          }
        }
        break;
      case 7:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useAdaptiveTimeStep = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use adaptive time step error";
            return false;
          }
        }
        break;
      case 8:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_timeStepRelaxationFactor = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Time step relaxation factor error";
            return false;
          }
        }
        break;
      case 9:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_advectionFlags.find(code);

            int advectionMode = -1;

            if (it != m_advectionFlags.end())
              advectionMode = it->second;

            switch (advectionMode)
            {
              case 1:
                m_advectionMode = AdvectionDiscretizationMode::Upwind;
                break;
              case 2:
                m_advectionMode = AdvectionDiscretizationMode::Central;
                break;
              case 3:
                m_advectionMode = AdvectionDiscretizationMode::Hybrid;
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Advection mode error";
            return false;
          }
        }
        break;
      case 10:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_computeDispersion = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Compute dispersion error";
            return false;
          }
        }
        break;
      case 11:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int heatSolverMode = -1;

            if (it != m_solverTypeFlags.end())
              heatSolverMode = it->second;

            switch (heatSolverMode)
            {
              case 1:
                m_heatSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_heatSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                m_heatSolver->setSolverType(ODESolver::CVODE_ADAMS);
                break;
              case 4:
                m_heatSolver->setSolverType(ODESolver::CVODE_BDF);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver type error";
            return false;
          }
        }
        break;
      case 12:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver absolute tolerance error";
            return false;
          }
        }
        break;
      case 13:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Temperature solver relative tolerance error";
            return false;
          }
        }
        break;
      case 14:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_waterDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water density error";
            return false;
          }
        }
        break;
      case 15:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_cp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Specific heat capacity of water error";
            return false;
          }
        }
        break;
      case 16:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            int numSolutes = options[1].toInt(&parsed);

            if (parsed)
              setNumSolutes(numSolutes);

            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of solutes error";
            return false;
          }
        }
        break;
      case 17:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Verbose error";
            return false;
          }
        }
        break;
      case 18:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_flushToDiskFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Flush to disk frequency error";
            return false;
          }
        }
        break;
      case 19:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_printFrequency = options[1].toInt(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Print frequency error";
            return false;
          }
        }
        break;
      case 20:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useEvaporation = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use evaporation error";
            return false;
          }
        }
        break;
      case 21:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useConvection = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use convection error";
            return false;
          }
        }
        break;
      case 22:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_evapWindFuncCoeffA = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Evaporation coefficient error";
            return false;
          }
        }
        break;
      case 23:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_evapWindFuncCoeffB = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Evaporation coefficient error";
            return false;
          }
        }
        break;
      case 24:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_bowensCoeff = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Bowens coefficient error";
            return false;
          }
        }
        break;
      case 25:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_pressureRatio = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Pressure ratio";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool STMModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  QString optionsFlag = options[0];

  if (options.size() == 2)
  {
    if (!QString::compare(optionsFlag, "csv", Qt::CaseInsensitive))
    {
      m_outputCSVFileInfo = QFileInfo(options[1]);
    }
    else if (!QString::compare(optionsFlag, "netcdf", Qt::CaseInsensitive))
    {
      m_outputNetCDFFileInfo = QFileInfo(options[1]);
    }
  }
  else
  {
    errorMessage = "Output file error";
    return false;
  }

  return true;
}

bool STMModel::readInputFileSolutesTag(const QString &line, QString &errorMessage)
{
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      std::string solverType = columns[1].toStdString();
      auto it = m_solverTypeFlags.find(solverType);

      if (it != m_solverTypeFlags.end())
      {
        int solverTypeCode = it->second;

        switch (solverTypeCode)
        {
          case 1:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RK4);
            break;
          case 2:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RKQS);
            break;
          case 3:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_ADAMS);
            break;
          case 4:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_BDF);
            break;
          default:
            foundError = true;
            break;
        }

        if (foundError)
        {
          errorMessage = "Solute error";
          return false;
        }

        bool parsed;
        double abs_tol = columns[2].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setAbsoluteTolerance(abs_tol);
        }
        else
        {
          errorMessage = "Solute absolute tolerance error";
          return false;
        }

        double rel_tol = columns[3].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setRelativeTolerance(rel_tol);
        }
        else
        {
          errorMessage = "Solute relative tolerance error";
          return false;
        }
      }
      else
      {
        errorMessage = "Solute error";
        return false;
      }

      m_addedSoluteCount++;
    }
    else
    {
      errorMessage = "Solute error count";
      return false;
    }
  }
  else
  {
    errorMessage = "Solute error";
    return false;
  }

  return true;
}

bool STMModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];

    bool workedX;
    bool workedY;
    bool workedZ;

    double x = columns[1].toDouble(&workedX);

    double y = columns[2].toDouble(&workedY);

    double z = columns[3].toDouble(&workedZ);

    if (workedX && workedY && workedZ)
    {
      addElementJunction(id.toStdString(), x, y, z);
    }
    else
    {
      errorMessage = "Junctions error";
      return false;
    }
  }
  else
  {
    errorMessage = "Junctions error";
    return false;
  }

  return true;
}

bool STMModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() > 9)
  {
    QString id = columns[0];
    QString fromId = columns[1];
    QString toId = columns[2];

    auto fromIt = m_elementJunctionsById.find(fromId.toStdString());
    auto toIt = m_elementJunctionsById.find(toId.toStdString());

    if (fromIt != m_elementJunctionsById.end() &&
        toIt != m_elementJunctionsById.end())
    {
      ElementJunction *ej1 = m_elementJunctionsById[fromId.toStdString()];
      ElementJunction *ej2 = m_elementJunctionsById[toId.toStdString()];

      bool lengthOk;
      double length = columns[3].toDouble(&lengthOk);

      bool depthOk ;
      double depth = columns[4].toDouble(&depthOk);

      bool xsectionAreaOk ;
      double xsectionArea = columns[5].toDouble(&xsectionAreaOk);

      bool widthOk ;
      double width = columns[6].toDouble(&widthOk);

      bool slopeOk ;
      double slope = columns[7].toDouble(&slopeOk);

      bool flowOk ;
      double flow = columns[8].toDouble(&flowOk);

      bool disperseCoeffOk ;
      double disperseCoeff = columns[9].toDouble(&disperseCoeffOk);

      bool tempOk ;
      double temp = columns[10].toDouble(&tempOk);

      if (lengthOk && depthOk && xsectionAreaOk &&
          widthOk && slopeOk && disperseCoeffOk &&
          tempOk && flowOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->length = length;
        element->depth = depth;
        element->xSectionArea = xsectionArea;
        element->width = width;
        element->slope = slope;
        element->longDispersion.value = disperseCoeff;
        element->temperature.value = temp;
        element->flow = flow;

        if (columns.size() > 10)
        {
          for (int i = 11; i < columns.size(); i++)
          {
            bool soluteOk ;
            double solute = columns[i].toDouble(&soluteOk);

            if (soluteOk && i - 11 < (int)m_solutes.size())
            {
              element->soluteConcs[i - 11].value = solute;
            }
            else
            {
              errorMessage = "Wrong initial solute concetration";
              return false;
            }
          }
        }
      }
      else
      {
        errorMessage = "";
        return false;
      }
    }
    else
    {
      errorMessage = "Wrong upstream or downstream junction";
      return false;
    }
  }

  return true;
}

bool STMModel::readInputFileBoundaryConditionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];
    auto it = m_elementJunctionsById.find(id.toStdString());

    if (it != m_elementJunctionsById.end())
    {
      ElementJunction *junction = it->second;

      bool found = false;
      QString type = columns[2];
      QString variable = columns[1].trimmed();

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {
        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          bool valueOk ;
          double value = columns[3].toDouble(&valueOk);

          if (valueOk)
          {
            junction->temperature.isBC = true;
            junction->temperature.value = value;
            found = true;
          }
          else
          {
            errorMessage = "Temperature BC value is invalid";
            return false;
          }
        }
        else
        {
          for (size_t i = 0; i < m_solutes.size(); i++)
          {
            std::string solute = m_solutes[i];

            if (!solute.compare(variable.toStdString()))
            {
              bool ok;
              double value = columns[3].toDouble(&ok);

              if (ok)
              {
                junction->soluteConcs[i].isBC = true;
                junction->soluteConcs[i].value = value;
                found = true;
              }
            }
          }
        }

        if (!found)
        {
          return false;
        }
      }
      else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
      {
        int variableIndex = -2;

        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          variableIndex = -1;
        }
        else
        {
          for (size_t i = 0; i < m_solutes.size(); i++)
          {
            std::string solute = m_solutes[i];

            if (!solute.compare(variable.toStdString()))
            {
              variableIndex = i;
              break;
            }
          }
        }

        if (variableIndex > -2)
        {
          QString filePath = columns[3];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if(TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                JunctionTimeSeriesBC *junctionBC = new JunctionTimeSeriesBC(junction, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  junctionBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(junctionBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified BC filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified BC filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified BC filepath does not exist";
            return false;
          }
        }
        else
        {
          errorMessage = "Variable specified for BC is not valid";
          return false;
        }
      }
    }
    else
    {
      errorMessage = "Boundary condition junction not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFilePointSourcesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];
    auto it = m_elementsById.find(id.toStdString());

    if (it != m_elementsById.end())
    {
      Element *element = it->second;
      QString variable = columns[1].trimmed();
      QString type = columns[2];
      PointSrcTimeSeriesBC::VariableType varType;

      int soluteIndex = -2;

      if (!QString::compare(variable, "HEAT", Qt::CaseInsensitive))
      {
        varType = PointSrcTimeSeriesBC::HeatSource;
        soluteIndex = -1;
      }
      else if (!QString::compare(variable, "FLOW", Qt::CaseInsensitive))
      {
        varType = PointSrcTimeSeriesBC::FlowSource;
        soluteIndex = -1;
      }
      else
      {
        varType = PointSrcTimeSeriesBC::SoluteSource;

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            soluteIndex = i;
            break;
          }
        }
      }

      if (soluteIndex > -2)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[3].toDouble(&valueOk);

          if (valueOk)
          {
            PointSrcTimeSeriesBC *pointSrcTSBC = new PointSrcTimeSeriesBC(element, varType, this);
            pointSrcTSBC->setSoluteIndex(soluteIndex);
            pointSrcTSBC->addValue(m_startDateTime, value);
            pointSrcTSBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(pointSrcTSBC);
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
        {

          QString filePath = columns[3];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                PointSrcTimeSeriesBC *pointSrcTSBC = new PointSrcTimeSeriesBC(element, varType  , this);
                pointSrcTSBC->setSoluteIndex(soluteIndex);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  pointSrcTSBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(pointSrcTSBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Point source is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Point source is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Point source is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Point source is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFileNonPointSourcesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 7)
  {
    QString idFrom = columns[0];
    QString idTo = columns[2];

    bool okStart;
    double startFactor = columns[1].toDouble(&okStart);

    bool okEnd;
    double endFactor = columns[3].toDouble(&okEnd);

    auto itFrom = m_elementsById.find(idFrom.toStdString());
    auto itTo = m_elementsById.find(idTo.toStdString());

    if (okStart && okEnd && itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *elementFrom = itFrom->second;
      Element *elementTo = itTo->second;

      QString variable = columns[4].trimmed();
      QString type = columns[5];
      NonPointSrcTimeSeriesBC::VariableType variableType;

      int soluteIndex = -2;

      if (!QString::compare(variable, "HEAT", Qt::CaseInsensitive))
      {
        soluteIndex = -1;
        variableType = NonPointSrcTimeSeriesBC::HeatSource;
      }
      else if (!QString::compare(variable, "FlOW", Qt::CaseInsensitive))
      {
        soluteIndex = -1;
        variableType = NonPointSrcTimeSeriesBC::FlowSource;
      }
      else
      {
        variableType = NonPointSrcTimeSeriesBC::SoluteSource;

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            soluteIndex = i;
            break;
          }
        }
      }

      if (soluteIndex > -2)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[6].toDouble(&valueOk);

          if (valueOk)
          {
            NonPointSrcTimeSeriesBC *nonPointSrcTSBC = new NonPointSrcTimeSeriesBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            nonPointSrcTSBC->setSoluteIndex(soluteIndex);
            nonPointSrcTSBC->addValue(m_startDateTime, value);
            nonPointSrcTSBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(nonPointSrcTSBC);
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
        {

          QString filePath = columns[6];

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                NonPointSrcTimeSeriesBC *nonPointSrcTSBC = new NonPointSrcTimeSeriesBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
                nonPointSrcTSBC->setSoluteIndex(soluteIndex);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  nonPointSrcTSBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(nonPointSrcTSBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Point source is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Point source is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Point source is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Point source is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Point source is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFileUniformHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_hydraulicVariableFlags.find(variableType.toStdString());

      if (it != m_hydraulicVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                   variableIndex, this);
            uniformHydraulicsBC->addValue(m_startDateTime, value);
            uniformHydraulicsBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(uniformHydraulicsBC);
          }
          else
          {
            errorMessage = "Uniform hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformHydraulicsTimeSeriesBC *uniformHydraulicsBC = new UniformHydraulicsTimeSeriesBC(fromElement, toElement,
                                                                                                       variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformHydraulicsBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformHydraulicsBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform hydraulics filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform hydraulics filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified uniform hydraulics filepath does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform hydraulic is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform hydraulic condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFileNonUniformHydraulicsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 2)
  {
    auto it = m_hydraulicVariableFlags.find(columns[0].toStdString());

    if (it != m_hydraulicVariableFlags.end())
    {
      int variableIndex = it->second;

      QString filePath = columns[1];

      if (!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                HydraulicsTimeSeriesBC *hydraulicsTimeSeries = new HydraulicsTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  hydraulicsTimeSeries->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(hydraulicsTimeSeries);
              }
              else
              {
                errorMessage = "Specified time varying hydraulic file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified time varying hydraulic file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified time varying hydraulic file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying hydraulic file is invalid";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

bool STMModel::readInputFileUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString type = columns[2];
      QString varValue = columns[3].trimmed();

      if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
      {

        bool valueOk;
        double value =  varValue.toDouble(&valueOk);

        if (valueOk)
        {
          UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);
          uniformRadiationFluxBC->addValue(m_startDateTime, value);
          uniformRadiationFluxBC->addValue(m_endDateTime, value);
          m_boundaryConditions.push_back(uniformRadiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "FILE", Qt::CaseInsensitive))
      {
        QString filePath = varValue;

        if (!filePath.isEmpty() && !filePath.isNull())
        {
          QFileInfo fileInfo(filePath);

          if (fileInfo.isRelative())
            fileInfo = relativePathToAbsolute(fileInfo);

          if (QFile::exists(fileInfo.absoluteFilePath()))
          {
            std::map<double, std::vector<double>> timeSeries;
            std::vector<std::string> headers;

            if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
            {
              UniformRadiativeFluxTimeSeriesBC *uniformRadiationFluxBC = new UniformRadiativeFluxTimeSeriesBC(fromElement, toElement, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[0];
                uniformRadiationFluxBC->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(uniformRadiationFluxBC);

              timeSeries.clear();
              headers.clear();
            }
            else
            {
              errorMessage = "Specified BC filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified BC filepath does not exist";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified BC filepath does not exist";
          return false;
        }
      }
    }
    else
    {
      errorMessage = "Boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFileNonUniformRadiativeFluxesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 1)
  {
    QString filePath = columns[0];

    if (!filePath.isEmpty() && !filePath.isNull())
    {
      QFileInfo fileInfo(filePath);

      if (fileInfo.isRelative())
        fileInfo = relativePathToAbsolute(fileInfo);

      if (QFile::exists(fileInfo.absoluteFilePath()))
      {
        std::map<double, std::vector<double>> timeSeries;
        std::vector<std::string> headers;

        if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
        {
          for (size_t i = 0; i < headers.size(); i++)
          {
            auto eit = m_elementsById.find(headers[i]);

            if (eit != m_elementsById.end())
            {
              Element *element = eit->second;
              RadiativeFluxTimeSeriesBC *radiativeFluxTimeSeries = new RadiativeFluxTimeSeriesBC(element, this);

              for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
              {
                double dateTime = it->first;
                double value = it->second[i];
                radiativeFluxTimeSeries->addValue(dateTime, value);
              }

              m_boundaryConditions.push_back(radiativeFluxTimeSeries);
            }
            else
            {
              errorMessage = "Specified time varying radiative flux file is invalid";
              return false;
            }
          }
        }
        else
        {
          errorMessage = "Specified time varying radiative flux file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying radiative flux file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying radiative flux file is invalid";
      return false;
    }

  }
  else
  {
    errorMessage = "Specified time varying radiative flux file is invalid";
    return false;
  }

  return true;
}

bool STMModel::readInputFileUniformMeteorologyTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 5)
  {
    QString fromId = columns[0];
    QString toId = columns[1];

    auto itFrom = m_elementsById.find(fromId.toStdString());
    auto itTo = m_elementsById.find(toId.toStdString());

    if (itFrom != m_elementsById.end() && itTo != m_elementsById.end())
    {
      Element *fromElement = itFrom->second;
      Element *toElement = itTo->second;

      QString variableType = columns[2];
      QString valueType = columns[3];
      QString varValue = columns[4].trimmed();

      auto it = m_meteorologicalVariableFlags.find(variableType.toStdString());

      if (it != m_meteorologicalVariableFlags.end())
      {
        int variableIndex = it->second;

        if (!QString::compare(valueType, "VALUE", Qt::CaseInsensitive))
        {

          bool valueOk;
          double value =  varValue.toDouble(&valueOk);

          if (valueOk)
          {
            UniformMeteorologyTimeSeriesBC *uniformMeteorologyBC = new UniformMeteorologyTimeSeriesBC(fromElement, toElement,
                                                                                                      variableIndex, this);
            uniformMeteorologyBC->addValue(m_startDateTime, value);
            uniformMeteorologyBC->addValue(m_endDateTime, value);
            m_boundaryConditions.push_back(uniformMeteorologyBC);
          }
          else
          {
            errorMessage = "Uniform meteorology value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "FILE", Qt::CaseInsensitive))
        {
          QString filePath = varValue;

          if (!filePath.isEmpty() && !filePath.isNull())
          {
            QFileInfo fileInfo(filePath);

            if (fileInfo.isRelative())
              fileInfo = relativePathToAbsolute(fileInfo);

            if (QFile::exists(fileInfo.absoluteFilePath()))
            {
              std::map<double, std::vector<double>> timeSeries;
              std::vector<std::string> headers;

              if(TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
              {
                UniformMeteorologyTimeSeriesBC *uniformMeteorologyBC = new UniformMeteorologyTimeSeriesBC(fromElement, toElement,
                                                                                                          variableIndex, this);
                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[0];
                  uniformMeteorologyBC->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(uniformMeteorologyBC);

                timeSeries.clear();
                headers.clear();
              }
              else
              {
                errorMessage = "Specified uniform meteorology filepath does not exist";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified uniform meteorology filepath does not exist";
              return false;
            }
          }
          else
          {
            errorMessage = "Specified uniform meteorology filepath does not exist";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Variable specified for uniform meteorology is incorrect";
        return false;
      }
    }
    else
    {
      errorMessage = "Uniform meteorology boundary condition element not found";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool STMModel::readInputFileNonUniformMeteorologyTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 2)
  {
    auto it = m_meteorologicalVariableFlags.find(columns[0].toStdString());

    if (it != m_meteorologicalVariableFlags.end())
    {
      int variableIndex = it->second;

      QString filePath = columns[1];

      if (!filePath.isEmpty() && !filePath.isNull())
      {
        QFileInfo fileInfo(filePath);

        if (fileInfo.isRelative())
          fileInfo = relativePathToAbsolute(fileInfo);

        if (QFile::exists(fileInfo.absoluteFilePath()))
        {
          std::map<double, std::vector<double>> timeSeries;
          std::vector<std::string> headers;

          if (TimeSeries::readTimeSeries(fileInfo, timeSeries, headers))
          {
            for (size_t i = 0; i < headers.size(); i++)
            {
              auto eit = m_elementsById.find(headers[i]);

              if (eit != m_elementsById.end())
              {
                Element *element = eit->second;
                MeteorologyTimeSeriesBC *meteorologyTimeSeries = new MeteorologyTimeSeriesBC(element, variableIndex, this);

                for (auto it = timeSeries.begin(); it != timeSeries.end(); it++)
                {
                  double dateTime = it->first;
                  double value = it->second[i];
                  meteorologyTimeSeries->addValue(dateTime, value);
                }

                m_boundaryConditions.push_back(meteorologyTimeSeries);
              }
              else
              {
                errorMessage = "Specified time varying hydraulic file is invalid";
                return false;
              }
            }
          }
          else
          {
            errorMessage = "Specified time varying hydraulic file is invalid";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified time varying hydraulic file is invalid";
          return false;
        }
      }
      else
      {
        errorMessage = "Specified time varying hydraulic file is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified time varying hydraulic file is invalid";
      return false;
    }
  }
  else
  {
    errorMessage = "Specified time varying hydraulic file is invalid";
    return false;
  }

  return true;
}

void STMModel::writeOutput()
{
  m_currentflushToDiskCount++;

  if (m_currentflushToDiskCount >= m_flushToDiskFrequency)
  {
    m_flushToDisk = true;
    m_currentflushToDiskCount = 0;
  }
  else
  {
    m_flushToDisk = false;
  }

  writeCSVOutput();
  writeNetCDFOutput();
}

void STMModel::writeCSVOutput()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->index
                        << ", " << element->x << ", " << element->y << ", " << element->z
                        << ", " << element->depth
                        << ", " << element->width
                        << ", " << element->xSectionArea
                        << ", " << element->longDispersion.value
                        << ", " << element->temperature.value
                        << ", " << element->totalAdvDispHeatBalance
                        << ", " << element->totalExternalHeatFluxesBalance
                        << ", " << element->totalRadiationFluxesHeatBalance
                        << ", " << element->totalEvaporativeHeatFluxesBalance
                        << ", " << element->totalConvectiveHeatFluxesBalance
                        << ", " << element->totalHeatBalance;

      for (size_t j = 0; j < m_solutes.size(); j++)
      {
        m_outputCSVStream << ", " << element->soluteConcs[j].value
                          << ", " << element->totalAdvDispSoluteMassBalance[j]
                             << ", " << element->totalExternalSoluteFluxesMassBalance[j]
                                << ", " << element->totalSoluteMassBalance[j];
      }

      m_outputCSVStream << endl;
    }

    if (m_flushToDisk)
    {
      m_outputCSVStream.flush();
    }
  }
}

void STMModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF

  if (m_outputNetCDF)
  {

    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();

    //Set current dateTime
    //      ThreadSafeNcVar timeVar =  m_outputNetCDF->getVar("time");
    //      timeVar.putVar(std::vector<size_t>({currentTime}), m_currentDateTime);
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);

    double *flow = new double[m_elements.size()];
    double *depth = new double[m_elements.size()];
    double *width = new double[m_elements.size()];
    double *xsectArea = new double[m_elements.size()];
    double *dispersion = new double[m_elements.size()];
    double *temperature = new double[m_elements.size()];
    double *elementEvapHeatFlux = new double[m_elements.size()];
    double *elementConvHeatFlux = new double[m_elements.size()];
    double *elementRadiationFlux = new double[m_elements.size()];
    double *elementHeatFlux = new double[m_elements.size()];
    double *elementAirTemp =  new double[m_elements.size()];
    double *elementRH =  new double[m_elements.size()];
    double *elementWindSpeed =  new double[m_elements.size()];
    double *elementVaporPress =  new double[m_elements.size()];
    double *elementSatVaporPress =  new double[m_elements.size()];
    double *elementAirVaporPress =  new double[m_elements.size()];
    double *elementAirSatVaporPress =  new double[m_elements.size()];
    double *totalElementHeatBalance = new double[m_elements.size()];
    double *totalElementAdvDispHeatBalance = new double[m_elements.size()];
    double *totalElementEvapHeatBalance = new double[m_elements.size()];
    double *totalElementConvHeatBalance = new double[m_elements.size()];
    double *totalElementRadiationFluxHeatBalance = new double[m_elements.size()];
    double *totalElementExternalHeatFluxBalance = new double[m_elements.size()];
    double *solutes = new double[m_elements.size() * m_solutes.size()];
    double *totalElementSoluteMassBalance = new double[m_elements.size() * m_solutes.size()];
    double *totalElementAdvDispSoluteMassBalance = new double[m_elements.size() * m_solutes.size()];
    double *totalElementExternalSoluteFluxMassBalance = new double[m_elements.size() * m_solutes.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      flow[i] = element->flow;
      depth[i] = element->depth;
      width[i] = element->width;
      xsectArea[i] = element->xSectionArea;
      dispersion[i] = element->longDispersion.value;
      temperature[i] = element->temperature.value;
      elementEvapHeatFlux[i] = element->evaporationHeatFlux;
      elementConvHeatFlux[i] = element->convectionHeatFlux;
      elementRadiationFlux[i] = element->radiationFluxes;
      elementHeatFlux[i] = element->externalHeatFluxes;
      totalElementHeatBalance[i] = element->totalHeatBalance;
      totalElementAdvDispHeatBalance[i] = element->totalAdvDispHeatBalance;
      totalElementRadiationFluxHeatBalance[i] = element->totalRadiationFluxesHeatBalance;
      totalElementExternalHeatFluxBalance[i] = element->totalExternalHeatFluxesBalance;
      totalElementEvapHeatBalance[i] = element->totalEvaporativeHeatFluxesBalance;
      totalElementConvHeatBalance[i] = element->totalConvectiveHeatFluxesBalance;
      elementAirTemp[i] = element->airTemperature;
      elementRH[i] = element->relativeHumidity;
      elementWindSpeed[i] = element->windSpeed;
      elementVaporPress[i] = element->vaporPressureWater;
      elementSatVaporPress[i] = element->saturationVaporPressureWater;
      elementAirVaporPress[i] = element->vaporPressureAir;
      elementAirSatVaporPress[i] = element->saturationVaporPressureAir;

      for (size_t j = 0; j < m_solutes.size(); j++)
      {
        solutes[j * m_elements.size() + i] = element->soluteConcs[j].value;
        totalElementSoluteMassBalance[j * m_elements.size() + i] = element->totalSoluteMassBalance[j];
        totalElementAdvDispSoluteMassBalance[j * m_elements.size() + i] = element->totalAdvDispSoluteMassBalance[j];
        totalElementExternalSoluteFluxMassBalance[j * m_elements.size() + i] = element->totalExternalSoluteFluxesMassBalance[j];
      }
    }

    m_outNetCDFVariables["flow"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), flow);

    m_outNetCDFVariables["depth"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);

    m_outNetCDFVariables["width"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);

    m_outNetCDFVariables["xsection_area"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), xsectArea);

    m_outNetCDFVariables["dispersion"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), dispersion);

    m_outNetCDFVariables["temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);

    m_outNetCDFVariables["element_evap_heat_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementEvapHeatFlux);

    m_outNetCDFVariables["element_conv_heat_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementConvHeatFlux);

    m_outNetCDFVariables["element_radiation_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementRadiationFlux);

    m_outNetCDFVariables["element_heat_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementHeatFlux);

    m_outNetCDFVariables["element_air_temp"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementAirTemp);

    m_outNetCDFVariables["element_relative_humidity"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementRH);

    m_outNetCDFVariables["element_wind_speed"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementWindSpeed);

    m_outNetCDFVariables["element_vapor_pressure"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementVaporPress);

    m_outNetCDFVariables["element_saturated_vapor_pressure"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementSatVaporPress);

    m_outNetCDFVariables["element_air_vapor_pressure"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementAirVaporPress);

    m_outNetCDFVariables["element_air_saturated_vapor_pressure"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementAirSatVaporPress);

    m_outNetCDFVariables["total_element_heat_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementHeatBalance);

    m_outNetCDFVariables["total_element_adv_disp_heat_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementAdvDispHeatBalance);

    m_outNetCDFVariables["total_element_evap_heat_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementEvapHeatBalance);

    m_outNetCDFVariables["total_element_conv_heat_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementConvHeatBalance);

    m_outNetCDFVariables["total_element_radiation_flux_heat_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementRadiationFluxHeatBalance);

    m_outNetCDFVariables["total_element_external_heat_flux_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), totalElementExternalHeatFluxBalance);

    if(m_solutes.size())
    {
      m_outNetCDFVariables["solute_concentration"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), solutes);
    }

    m_outNetCDFVariables["total_element_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementSoluteMassBalance);

    m_outNetCDFVariables["total_element_adv_disp_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementAdvDispSoluteMassBalance);

    m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size()}), totalElementExternalSoluteFluxMassBalance);

    m_outNetCDFVariables["total_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);

    m_outNetCDFVariables["total_adv_disp_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalAdvDispHeatBalance);

    m_outNetCDFVariables["total_evap_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalEvaporationHeatBalance);

    m_outNetCDFVariables["total_conv_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalConvectiveHeatBalance);

    m_outNetCDFVariables["total_radiation_flux_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);

    m_outNetCDFVariables["total_external_heat_flux_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);

    m_outNetCDFVariables["total_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalSoluteMassBalance.data());

    m_outNetCDFVariables["total_adv_disp_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalAdvDispSoluteMassBalance.data());

    m_outNetCDFVariables["total_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_solutes.size()}), m_totalExternalSoluteFluxMassBalance.data());

    delete[] flow;
    delete[] depth;
    delete[] width;
    delete[] xsectArea;
    delete[] dispersion;
    delete[] temperature;
    delete[] elementEvapHeatFlux;
    delete[] elementConvHeatFlux;
    delete[] elementRadiationFlux;
    delete[] elementHeatFlux;
    delete[] totalElementHeatBalance;
    delete[] totalElementAdvDispHeatBalance;
    delete[] totalElementEvapHeatBalance;
    delete[] totalElementConvHeatBalance;
    delete[] totalElementRadiationFluxHeatBalance;
    delete[] totalElementExternalHeatFluxBalance;
    delete[] solutes;
    delete[] totalElementSoluteMassBalance;
    delete[] totalElementAdvDispSoluteMassBalance;
    delete[] totalElementExternalSoluteFluxMassBalance;
    delete[] elementAirTemp ;
    delete[] elementRH ;
    delete[] elementWindSpeed ;
    delete[] elementVaporPress ;
    delete[] elementSatVaporPress ;
    delete[] elementAirVaporPress ;
    delete[] elementAirSatVaporPress ;

    if(m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }
  }
#endif
}

void STMModel::closeOutputFiles()
{
  closeCSVOutputFile();
  closeOutputNetCDFFile();
}

void STMModel::closeCSVOutputFile()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void STMModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF

  if(m_outputNetCDF)
  {
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }

#endif
}

QFileInfo STMModel::relativePathToAbsolute(const QFileInfo &fileInfo)
{
  if (fileInfo.isRelative())
  {
    if (!m_inputFile.filePath().isEmpty() &&
        !m_inputFile.filePath().isNull() &&
        QFile::exists(m_inputFile.absoluteFilePath()))
    {
      QFileInfo absoluteFilePath = m_inputFile.absoluteDir().absoluteFilePath(fileInfo.filePath());

      if (absoluteFilePath.absoluteDir().exists())
      {
        return absoluteFilePath;
      }
    }
  }

  return fileInfo;
}


const unordered_map<string, int> STMModel::m_inputFileFlags({
                                                              {"[OPTIONS]", 1},
                                                              {"[OUTPUTS]", 2},
                                                              {"[SOLUTES]", 3},
                                                              {"[ELEMENTJUNCTIONS]", 4},
                                                              {"[ELEMENTS]", 5},
                                                              {"[BOUNDARY_CONDITIONS]", 6},
                                                              {"[POINT_SOURCES]", 7},
                                                              {"[NON_POINT_SOURCES]", 8},
                                                              {"[UNIFORM_HYDRAULICS]", 9},
                                                              {"[NON_UNIFORM_HYDRAULICS]", 10},
                                                              {"[UNIFORM_RADIATIVE_FLUXES]", 11},
                                                              {"[NON_UNIFORM_RADIATIVE_FLUXES]", 12},
                                                              {"[UNIFORM_METEOROLOGY]", 13},
                                                              {"[NON_UNIFORM_METEOROLOGY]", 14}
                                                            });

const unordered_map<string, int> STMModel::m_optionsFlags({
                                                            {"START_DATETIME", 1},
                                                            {"END_DATETIME", 2},
                                                            {"REPORT_INTERVAL", 3},
                                                            {"MAX_TIME_STEP", 4},
                                                            {"MIN_TIME_STEP", 5},
                                                            {"NUM_INITIAL_FIXED_STEPS", 6},
                                                            {"USE_ADAPTIVE_TIME_STEP", 7},
                                                            {"TIME_STEP_RELAXATION_FACTOR", 8},
                                                            {"ADVECTION_MODE", 9},
                                                            {"COMPUTE_DISPERSION", 10},
                                                            {"TEMP_SOLVER", 11},
                                                            {"TEMP_SOLVER_ABS_TOL", 12},
                                                            {"TEMP_SOLVER_REL_TOL", 13},
                                                            {"WATER_DENSITY", 14},
                                                            {"WATER_SPECIFIC_HEAT_CAPACITY", 15},
                                                            {"NUM_SOLUTES", 16},
                                                            {"VERBOSE", 17},
                                                            {"FLUSH_TO_DISK_FREQ", 18},
                                                            {"PRINT_FREQ", 19},
                                                            {"EVAPORATION", 20},
                                                            {"CONVECTION", 21},
                                                            {"WIND_FUNC_COEFF_A", 22},
                                                            {"WIND_FUNC_COEFF_B", 23},
                                                            {"BOWENS_COEFF", 24},
                                                            {"PRESSURE_RATIO", 25},
                                                            {"TVD_SCHEME", 26},
                                                            {"NETCDF_THREAD_SAFE", 27},

                                                          });

const unordered_map<string, int> STMModel::m_advectionFlags({
                                                              {"UPWIND", 1},
                                                              {"CENTRAL", 2},
                                                              {"HYBRID", 3},
                                                              {"TVD", 4},
                                                            });

const unordered_map<string, int> STMModel::m_solverTypeFlags({{"RK4", 1},
                                                              {"RKQS", 2},
                                                              {"ADAMS", 3},
                                                              {"BDF", 4}});

const unordered_map<string, int> STMModel::m_hydraulicVariableFlags({{"DEPTH", 1},
                                                                     {"WIDTH", 2},
                                                                     {"XSECTION_AREA", 3},
                                                                     {"FLOW", 4}});

const unordered_map<string, int> STMModel::m_meteorologicalVariableFlags({{"RELATIVE_HUMIDITY", 1},
                                                                          {"AIR_TEMPERATURE", 2},
                                                                          {"WIND_SPEED", 3}});

const QRegExp STMModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
