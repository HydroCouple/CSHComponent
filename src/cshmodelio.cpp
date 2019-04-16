/*!
*  \file    CSHModelio.cpp
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
#include "junctionbc.h"
#include "radiativefluxbc.h"
#include "hydraulicsbc.h"
#include "sourcebc.h"
#include "meteorologybc.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "temporal/timeseries.h"

#include <QDir>
#include <QDate>
#include <cstdlib>
#include <errno.h>

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;


#endif

using namespace std;

bool CSHModel::verbose() const
{
  return m_verbose;
}

void CSHModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int CSHModel::printFrequency() const
{
  return m_printFrequency;
}

void CSHModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int CSHModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void CSHModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo CSHModel::inputFile() const
{
  return m_inputFile;
}

void CSHModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo CSHModel::outputCSVFile() const
{
  return m_outputCSVFileInfo;
}

void CSHModel::setOutputCSVFile(const QFileInfo &outputFile)
{
  m_outputCSVFileInfo = outputFile;
}

QFileInfo CSHModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void CSHModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void CSHModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("CSH TimeStep (s): %f\tDateTime: %f\tIters: %i/%i\tTemp (°C) { Min: %f\tMax: %f\tTotalHeatBalance: %g (KJ)}", m_timeStep, m_currentDateTime,
           m_odeSolver->getIterations(), m_odeSolver->maxIterations(), m_minTemp, m_maxTemp, m_totalHeatBalance);

    for (size_t j = 0; j < m_solutes.size(); j++)
    {
      std::string &solute = m_solutes[j];

      if(solute == "WATER_AGE")
      {
        printf("\t%s (days) { Min: %f\tMax: %f}", solute.c_str(), m_minSolute[j], m_maxSolute[j]);
      }
      else
      {
        printf("\t%s (kg/m^3) { Min: %f\tMax: %f\tTotalMassBalance: %g (kg)}", solute.c_str(), m_minSolute[j], m_maxSolute[j], m_totalSoluteMassBalance[j]);
      }
    }

    printf("\n");

    m_currentPrintCount = 0;
  }
}

void CSHModel::saveAs(const QFileInfo &filePath)
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

bool CSHModel::initializeInputFiles(list<string> &errors)
{

  if (QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if (file.open(QIODevice::ReadOnly))
    {

      m_timeSeries.clear();

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

          auto it = m_inputFileFlags.find(line.toUpper().toStdString());

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
              case 12:
                readSuccess = readInputFileElementHydraulicVariablesTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileBoundaryConditionsTag(line, error);
                break;
              case 7:
                readSuccess = readInputFileSourcesTag(line, error);
                break;
              case 8:
                readSuccess = readInputFileHydraulicsTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileRadiativeFluxesTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileMeteorologyTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileTimeSeriesTag(line, error);
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

bool CSHModel::initializeOutputFiles(list<string> &errors)
{
  return initializeCSVOutputFile(errors) &&
      initializeNetCDFOutputFile(errors);
}

bool CSHModel::initializeCSVOutputFile(list<string> &errors)
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

bool CSHModel::initializeNetCDFOutputFile(list<string> &errors)
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
    timeVar.putAtt("long_name", "Time");
    timeVar.putAtt("standard_name", "time");
    //    timeVar.putAtt("units", "days since 1473-01-01 00:00:00 BCE");
    timeVar.putAtt("calendar", "julian");
    m_outNetCDFVariables["time"] = timeVar;

    //Add Solutes
    ThreadSafeNcDim solutesDim =  m_outputNetCDF->addDim("solutes", m_numSolutes);


    //Add element junctions
    ThreadSafeNcDim junctionDim =  m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers =  m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("long_name", "Element Junction Identifiers");
    m_outNetCDFVariables["element_junction_id"] = junctionIdentifiers;

    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_FLOAT, junctionDim);
    junctionX.putAtt("long_name", "Junction X Coordinate");
    junctionX.putAtt("units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_FLOAT, junctionDim);
    junctionY.putAtt("long_name", "Junction Y Coordinate");
    junctionY.putAtt("units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_FLOAT, junctionDim);
    junctionZ.putAtt("long_name", "Junction Z Coordinate");
    junctionZ.putAtt("units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    float *vertx = new float[m_elementJunctions.size()];
    float *verty = new float[m_elementJunctions.size()];
    float *vertz = new float[m_elementJunctions.size()];
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
    elementIdentifiers.putAtt("long_name", "Element Identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;

    ThreadSafeNcVar elementFromJunction =  m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("long_name", "Upstream Junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction =  m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("long_name", "Downstream Junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;

    ThreadSafeNcVar element_x =  m_outputNetCDF->addVar("element_x", NcType::nc_FLOAT, elementsDim);
    element_x.putAtt("standard_name", "projection_x_coordinate");
    element_x.putAtt("long_name", "Element X Coordinate");
    element_x.putAtt("units", "m");
    m_outNetCDFVariables["element_x"] = element_x;

    ThreadSafeNcVar element_y =  m_outputNetCDF->addVar("element_y", NcType::nc_FLOAT, elementsDim);
    element_y.putAtt("standard_name", "projection_y_coordinate");
    element_y.putAtt("long_name", "Element Y Coordinate");
    element_y.putAtt("units", "m");
    m_outNetCDFVariables["element_y"] = element_y;

    //    ThreadSafeNcVar elementsVar =  m_outputNetCDF->addVar("elements", NcType::nc_FLOAT, elementsDim);
    //    elementsVar.putAtt("long_name", "Distance");
    //    elementsVar.putAtt("units", "m");
    //    m_outNetCDFVariables["elements"] = elementsVar;


    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    float *elX = new float[m_elements.size()];
    float *elY = new float[m_elements.size()];
    //    double *els = new double[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->tIndex;
      toJunctions[i] = element->downstreamJunction->tIndex;

      elX[i] = element->x;
      elY[i] = element->y;
      //      els[i] = element->distanceFromUpStreamJunction;
    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    element_x.putVar(elX);
    element_y.putVar(elY);
    //    elementsVar.putVar(els);

    delete[] fromJunctions;
    delete[] toJunctions;
    delete[] elX;
    delete[] elY;
    //    delete[] els;

    for (size_t i = 0; i < m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;


    //hydraulics variables
    ThreadSafeNcVar lengthVar =  m_outputNetCDF->addVar("length", "float",
                                                        std::vector<std::string>({"elements"}));
    lengthVar.putAtt("long_name", "Element Length");
    lengthVar.putAtt("units", "m");
    m_outNetCDFVariables["length"] = lengthVar;

    std::vector<float> lengths(m_elements.size());

    for(size_t i = 0; i < m_elements.size(); i++)
    {
      lengths[i] = m_elements[i]->length;
    }

    lengthVar.putVar(lengths.data());

    //hydraulics variables
    ThreadSafeNcVar flowVar =  m_outputNetCDF->addVar("flow", "float",
                                                      std::vector<std::string>({"time", "elements"}));
    flowVar.putAtt("long_name", "Flow");
    flowVar.putAtt("units", "m^3/s");
    m_outNetCDFVariables["flow"] = flowVar;

    ThreadSafeNcVar depthVar =  m_outputNetCDF->addVar("depth", "float",
                                                       std::vector<std::string>({"time", "elements"}));
    depthVar.putAtt("long_name", "Flow Depth");
    depthVar.putAtt("units", "m");
    m_outNetCDFVariables["depth"] = depthVar;

    ThreadSafeNcVar widthVar =  m_outputNetCDF->addVar("width", "float",
                                                       std::vector<std::string>({"time", "elements"}));
    widthVar.putAtt("long_name", "Flow Top Width");
    widthVar.putAtt("units", "m");
    m_outNetCDFVariables["width"] = widthVar;

    ThreadSafeNcVar xsectAreaVar =  m_outputNetCDF->addVar("xsection_area", "float",
                                                           std::vector<std::string>({"time", "elements"}));
    xsectAreaVar.putAtt("long_name", "Flow Cross-Sectional Area");
    xsectAreaVar.putAtt("units", "m^2");
    m_outNetCDFVariables["xsection_area"] = xsectAreaVar;

    ThreadSafeNcVar dispersionVar =  m_outputNetCDF->addVar("dispersion", "float",
                                                            std::vector<std::string>({"time", "elements"}));
    dispersionVar.putAtt("long_name", "Longitudinal Dispersion");
    dispersionVar.putAtt("units", "m^2/s");
    m_outNetCDFVariables["dispersion"] = dispersionVar;

    ThreadSafeNcVar temperatureVar =  m_outputNetCDF->addVar("temperature", "float",
                                                             std::vector<std::string>({"time", "elements"}));
    temperatureVar.putAtt("long_name", "Temperature");
    temperatureVar.putAtt("units", "°C");
    m_outNetCDFVariables["temperature"] = temperatureVar;

    ThreadSafeNcVar volumeTimeDerivativeVar =  m_outputNetCDF->addVar("volume_time_derivative", "float",
                                                             std::vector<std::string>({"time", "elements"}));
    volumeTimeDerivativeVar.putAtt("long_name", "Volume Time Derivative");
    volumeTimeDerivativeVar.putAtt("units", "m^3/s");
    m_outNetCDFVariables["volume_time_derivative"] = volumeTimeDerivativeVar;

    ThreadSafeNcVar waterAgeVar =  m_outputNetCDF->addVar("water_age", "float",
                                                          std::vector<std::string>({"time", "elements"}));
    waterAgeVar.putAtt("long_name", "Water Age");
    waterAgeVar.putAtt("units", "days");
    m_outNetCDFVariables["water_age"] = waterAgeVar;


    ThreadSafeNcVar totalElementHeatBalanceVar =  m_outputNetCDF->addVar("total_element_heat_balance", "float",
                                                                         std::vector<std::string>({"time", "elements"}));
    totalElementHeatBalanceVar.putAtt("long_name", "Total Element Heat Balance");
    totalElementHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_heat_balance"] = totalElementHeatBalanceVar;

    ThreadSafeNcVar totalElementAdvDispHeatBalanceVar =  m_outputNetCDF->addVar("total_element_adv_disp_heat_balance", "float",
                                                                                std::vector<std::string>({"time", "elements"}));
    totalElementAdvDispHeatBalanceVar.putAtt("long_name", "Total Element Advection Dispersion Heat Balance");
    totalElementAdvDispHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_adv_disp_heat_balance"] = totalElementAdvDispHeatBalanceVar;

    ThreadSafeNcVar totalElementEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_element_evap_heat_balance", "float",
                                                                             std::vector<std::string>({"time", "elements"}));
    totalElementEvapHeatBalanceVar.putAtt("long_name", "Total Element Evaporation Heat Balance");
    totalElementEvapHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_evap_heat_balance"] = totalElementEvapHeatBalanceVar;

    ThreadSafeNcVar totalElementConvHeatBalanceVar =  m_outputNetCDF->addVar("total_element_conv_heat_balance", "float",
                                                                             std::vector<std::string>({"time", "elements"}));
    totalElementConvHeatBalanceVar.putAtt("long_name", "Total Element Convection Heat Balance");
    totalElementConvHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_conv_heat_balance"] = totalElementConvHeatBalanceVar;


    ThreadSafeNcVar totalElementRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_element_radiation_flux_heat_balance", "float",
                                                                                      std::vector<std::string>({"time", "elements"}));
    totalElementRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Element Radiation Flux Heat Balance");
    totalElementRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_radiation_flux_heat_balance"] = totalElementRadiationFluxHeatBalanceVar;

    ThreadSafeNcVar totalElementExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_element_external_heat_flux_balance", "float",
                                                                                     std::vector<std::string>({"time", "elements"}));
    totalElementExternalHeatFluxBalanceVar.putAtt("long_name", "Total Element External Heat Flux Balance");
    totalElementExternalHeatFluxBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_element_external_heat_flux_balance"] = totalElementExternalHeatFluxBalanceVar;

    ThreadSafeNcVar totalElementSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_element_solute_mass_balance", "float",
                                                                               std::vector<std::string>({"time", "solutes", "elements"}));
    totalElementSoluteMassBalanceVar.putAtt("long_name", "Total Element Solute Mass Balance");
    totalElementSoluteMassBalanceVar.putAtt("units", "kg");
    m_outNetCDFVariables["total_element_solute_mass_balance"] = totalElementSoluteMassBalanceVar;

    ThreadSafeNcVar totalElementAdvDispSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_element_adv_disp_solute_mass_balance", "float",
                                                                                      std::vector<std::string>({"time", "solutes", "elements"}));
    totalElementAdvDispSoluteMassBalanceVar.putAtt("long_name", "Total Element Advection Dispersion Solute Mass Balance");
    totalElementAdvDispSoluteMassBalanceVar.putAtt("units", "kg");
    m_outNetCDFVariables["total_element_adv_disp_solute_mass_balance"] = totalElementAdvDispSoluteMassBalanceVar;

    ThreadSafeNcVar totalElementExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_element_external_solute_flux_mass_balance", "float",
                                                                                           std::vector<std::string>({"time", "solutes"}));
    totalElementExternalSoluteFluxMassBalanceVar.putAtt("long_name", "Total External Solute Flux Mass Balance");
    totalElementExternalSoluteFluxMassBalanceVar.putAtt("units", "kg");
    m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"] = totalElementExternalSoluteFluxMassBalanceVar;

    ThreadSafeNcVar elementEvapHeatFluxVar =  m_outputNetCDF->addVar("element_evap_heat_flux", "float",
                                                                     std::vector<std::string>({"time", "elements"}));

    elementEvapHeatFluxVar.putAtt("long_name", "Element Evaporation Heat Flux");
    elementEvapHeatFluxVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["element_evap_heat_flux"] = elementEvapHeatFluxVar;

    ThreadSafeNcVar elementConvHeatFluxVar =  m_outputNetCDF->addVar("element_conv_heat_flux", "float",
                                                                     std::vector<std::string>({"time", "elements"}));

    elementConvHeatFluxVar.putAtt("long_name", "Element Convective Heat Flux");
    elementConvHeatFluxVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["element_conv_heat_flux"] = elementConvHeatFluxVar;


    ThreadSafeNcVar elementRadiationFluxVar =  m_outputNetCDF->addVar("element_radiation_flux", "float",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementRadiationFluxVar.putAtt("long_name", "Element Radiation Flux");
    elementRadiationFluxVar.putAtt("units", "W/m^2");
    m_outNetCDFVariables["element_radiation_flux"] = elementRadiationFluxVar;

    ThreadSafeNcVar elementHeatFluxVar =  m_outputNetCDF->addVar("element_heat_flux", "float",
                                                                 std::vector<std::string>({"time", "elements"}));

    elementHeatFluxVar.putAtt("long_name", "Element Heat Flux");
    elementHeatFluxVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_heat_flux"] = elementHeatFluxVar;


    ThreadSafeNcVar elementAirTempVar =  m_outputNetCDF->addVar("element_air_temp", "float",
                                                                std::vector<std::string>({"time", "elements"}));

    elementAirTempVar.putAtt("long_name", "Air Temperature");
    elementAirTempVar.putAtt("units", "C");
    m_outNetCDFVariables["element_air_temp"] = elementAirTempVar;


    ThreadSafeNcVar elementRHVar =  m_outputNetCDF->addVar("element_relative_humidity", "float",
                                                           std::vector<std::string>({"time", "elements"}));

    elementRHVar.putAtt("long_name", "Relative Humidity");
    elementRHVar.putAtt("units", "%");
    m_outNetCDFVariables["element_relative_humidity"] = elementRHVar;


    ThreadSafeNcVar elementWindSpeedVar =  m_outputNetCDF->addVar("element_wind_speed", "float",
                                                                  std::vector<std::string>({"time", "elements"}));

    elementWindSpeedVar.putAtt("long_name", "Wind Speed");
    elementWindSpeedVar.putAtt("units", "m/s");
    m_outNetCDFVariables["element_wind_speed"] = elementWindSpeedVar;

    ThreadSafeNcVar elementVaporPressVar =  m_outputNetCDF->addVar("element_vapor_pressure", "float",
                                                                   std::vector<std::string>({"time", "elements"}));

    elementVaporPressVar.putAtt("long_name", "Vapor Pressure");
    elementVaporPressVar.putAtt("units", "kPa");
    m_outNetCDFVariables["element_vapor_pressure"] = elementVaporPressVar;

    ThreadSafeNcVar elementSatVaporPressVar =  m_outputNetCDF->addVar("element_saturated_vapor_pressure", "float",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementSatVaporPressVar.putAtt("long_name", "Saturated Vapor Pressure");
    elementSatVaporPressVar.putAtt("units", "kPa");
    m_outNetCDFVariables["element_saturated_vapor_pressure"] = elementSatVaporPressVar;


    ThreadSafeNcVar elementAirVaporPressVar =  m_outputNetCDF->addVar("element_air_vapor_pressure", "float",
                                                                      std::vector<std::string>({"time", "elements"}));

    elementAirVaporPressVar.putAtt("long_name", "Air Vapor Pressure");
    elementAirVaporPressVar.putAtt("units", "kPa");
    m_outNetCDFVariables["element_air_vapor_pressure"] = elementAirVaporPressVar;

    ThreadSafeNcVar elementAirSatVaporPressVar =  m_outputNetCDF->addVar("element_air_saturated_vapor_pressure", "float",
                                                                         std::vector<std::string>({"time", "elements"}));

    elementAirSatVaporPressVar.putAtt("long_name", "Saturated Air Vapor Pressure");
    elementAirSatVaporPressVar.putAtt("units", "kPa");
    m_outNetCDFVariables["element_air_saturated_vapor_pressure"] = elementAirSatVaporPressVar;

    ThreadSafeNcVar totalHeatBalanceVar =  m_outputNetCDF->addVar("total_heat_balance", "float",
                                                                  std::vector<std::string>({"time"}));
    totalHeatBalanceVar.putAtt("long_name", "Total Heat Balance");
    totalHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_heat_balance"] = totalHeatBalanceVar;

    ThreadSafeNcVar totalAdvDispHeatBalanceVar =  m_outputNetCDF->addVar("total_adv_disp_heat_balance", "float",
                                                                         std::vector<std::string>({"time"}));
    totalAdvDispHeatBalanceVar.putAtt("long_name", "Total Advection Dispersion Heat Balance");
    totalAdvDispHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_adv_disp_heat_balance"] = totalAdvDispHeatBalanceVar;


    ThreadSafeNcVar totalEvapHeatBalanceVar =  m_outputNetCDF->addVar("total_evap_heat_balance", "float",
                                                                      std::vector<std::string>({"time"}));
    totalEvapHeatBalanceVar.putAtt("long_name", "Total Evaporation Heat Balance");
    totalEvapHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_evap_heat_balance"] = totalEvapHeatBalanceVar;


    ThreadSafeNcVar totalConvHeatBalanceVar =  m_outputNetCDF->addVar("total_conv_heat_balance", "float",
                                                                      std::vector<std::string>({"time"}));
    totalConvHeatBalanceVar.putAtt("long_name", "Total Convection Heat Balance");
    totalConvHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_conv_heat_balance"] = totalConvHeatBalanceVar;


    ThreadSafeNcVar totalRadiationFluxHeatBalanceVar =  m_outputNetCDF->addVar("total_radiation_flux_heat_balance", "float",
                                                                               std::vector<std::string>({"time"}));
    totalRadiationFluxHeatBalanceVar.putAtt("long_name", "Total Radiation Flux Heat Balance");
    totalRadiationFluxHeatBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_radiation_flux_heat_balance"] = totalRadiationFluxHeatBalanceVar;

    ThreadSafeNcVar totalExternalHeatFluxBalanceVar =  m_outputNetCDF->addVar("total_external_heat_flux_balance", "float",
                                                                              std::vector<std::string>({"time"}));
    totalExternalHeatFluxBalanceVar.putAtt("long_name", "Total External Heat Flux Balance");
    totalExternalHeatFluxBalanceVar.putAtt("units", "KJ");
    m_outNetCDFVariables["total_external_heat_flux_balance"] = totalExternalHeatFluxBalanceVar;


    ThreadSafeNcVar solutes =  m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
    solutes.putAtt("long_name", "Solutes");
    m_outNetCDFVariables["solutes"] = solutes;

    if (m_numSolutes)
    {
      char **soluteNames = new char *[m_numSolutes];

      for (int i = 0; i < m_numSolutes; i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (int i = 0; i < m_numSolutes; i++)
      {
        delete[] soluteNames[i];
      }

      delete[] soluteNames;
    }


    ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "float",
                                                         std::vector<std::string>({"time", "solutes", "elements"}));
    solutesVar.putAtt("long_name", "Solute Concentration");
    solutesVar.putAtt("units", "kg/m^3");
    m_outNetCDFVariables["solute_concentration"] = solutesVar;

    ThreadSafeNcVar totalSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_solute_mass_balance", "float",
                                                                        std::vector<std::string>({"time", "solutes"}));
    totalSoluteMassBalanceVar.putAtt("long_name", "Total Solute Mass Balance");
    totalSoluteMassBalanceVar.putAtt("units", "kg");
    m_outNetCDFVariables["total_solute_mass_balance"] = totalSoluteMassBalanceVar;

    ThreadSafeNcVar totalAdvDispSoluteMassBalanceVar =  m_outputNetCDF->addVar("total_adv_disp_solute_mass_balance", "float",
                                                                               std::vector<std::string>({"time", "solutes"}));
    totalAdvDispSoluteMassBalanceVar.putAtt("long_name", "Total Advection Dispersion Solute Mass Balance");
    totalAdvDispSoluteMassBalanceVar.putAtt("units", "kg");
    m_outNetCDFVariables["total_adv_disp_solute_mass_balance"] = totalAdvDispSoluteMassBalanceVar;

    ThreadSafeNcVar totalExternalSoluteFluxMassBalanceVar =  m_outputNetCDF->addVar("total_external_solute_flux_mass_balance", "float",
                                                                                    std::vector<std::string>({"time", "solutes"}));
    totalExternalSoluteFluxMassBalanceVar.putAtt("long_name", "Total External Solute Flux Mass Balance");
    totalExternalSoluteFluxMassBalanceVar.putAtt("units", "kg");
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

bool CSHModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
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
              case 4:
                m_advectionMode = AdvectionDiscretizationMode::TVD;
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
            if(QString::compare(options[1], "No", Qt::CaseInsensitive))
              m_computeDispersion = 1.0;
            else
              m_computeDispersion = 0.0;
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
      case 29:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int odeSolver = -1;

            if (it != m_solverTypeFlags.end())
              odeSolver = it->second;

            switch (odeSolver)
            {
              case 1:
                m_odeSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_odeSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_ADAMS);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::FUNCTIONAL);
                }
                break;
              case 4:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_BDF);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::NEWTON);
                  m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                }
                break;
              case 5:
                m_odeSolver->setSolverType(ODESolver::EULER);
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
            errorMessage = "Solver type error";
            return false;
          }
        }
        break;
      case 12:
      case 30:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Solver absolute tolerance error";
            return false;
          }
        }
        break;
      case 13:
      case 31:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Solver relative tolerance error";
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
      case 26:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            int fluxLimiter = options[1].toInt(&ok);
            foundError = !ok;

            m_TVDFluxLimiter = (ElementAdvTVD::TVDFluxLimiter)fluxLimiter;
          }
          else
          {
            foundError = true;
          }

          if(foundError)
          {
            errorMessage = "TVD flux limiter error";
            return false;
          }
        }
        break;
      case 28:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_simulateWaterAge = QString::compare(options[1], "No", Qt::CaseInsensitive);

            if(m_simulateWaterAge)
              setNumSolutes(m_numSolutes);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Simulate water age error";
            return false;
          }
        }
        break;
      case 32:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_linearSolverTypeFlags.find(code);

            int linearSolver = -1;

            if (it != m_linearSolverTypeFlags.end())
              linearSolver = it->second;

            switch (linearSolver)
            {
              case 1:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                break;
              case 2:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::FGMRES);
                break;
              case 3:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::Bi_CGStab);
                break;
              case 4:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::TFQMR);
                break;
              case 5:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::PCG);
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
            errorMessage = "Linear solver type error";
            return false;
          }
        }
        break;
      case 33:
        {

          bool foundError = false;

          if (options.size() == 2)
          {
            m_solveHydraulics = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }


          if (foundError)
          {
            errorMessage = "Solve hydraulics tag";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool CSHModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
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

bool CSHModel::readInputFileSolutesTag(const QString &line, QString &errorMessage)
{
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() >= 2)
  {
    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      bool parsed;
      double first_order_k = columns[1].toDouble(&parsed);

      if(parsed)
      {
        m_solute_first_order_k[m_addedSoluteCount] = first_order_k;
      }
      else
      {
        errorMessage = "Invalid solute first order reaction rate";
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

bool CSHModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
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

bool CSHModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
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
        element->width = element->bottomWidth = width;
        element->slope = slope;
        element->longDispersion.value = disperseCoeff;
        element->temperature.value = temp;
        element->flow.value  = element->prevFlow.value = flow;

        if(m_solutes.size() && columns.size() > 10)
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

bool CSHModel::readInputFileElementHydraulicVariablesTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 6)
  {
    QString id = columns[0];


    auto idIt = m_elementsById.find(id.toStdString());

    if(idIt != m_elementsById.end())
    {
      Element *element = m_elementsById[id.toStdString()];

      bool manOk;
      double man = columns[1].toDouble(&manOk);



      Element::XSectType xSectType = Element::RECT;

      QString xsect = columns[2];

      if(!xsect.compare("TRAP", Qt::CaseInsensitive))
      {
        xSectType = Element::TRAP;
      }
      else if(!xsect.compare("IRREGULAR", Qt::CaseInsensitive))
      {

      }

      bool bottomWidthOk ;
      double bottomWidth = columns[3].toDouble(&bottomWidthOk);

      bool slope1Ok ;
      double slope1 = columns[4].toDouble(&slope1Ok);

      bool slope2Ok ;
      double slope2 = columns[5].toDouble(&slope2Ok);

      if (manOk && bottomWidthOk && slope1Ok && slope2Ok)
      {
        element->bottomWidth = bottomWidth;
        element->mannings = man;
        element->sideSlopes[0] = slope1;
        element->sideSlopes[1] = slope2;
      }
      else
      {
        errorMessage = "";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified element not found";
      return false;
    }
  }

  return true;
}

bool CSHModel::readInputFileBoundaryConditionsTag(const QString &line, QString &errorMessage)
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
            JunctionBC *junctionBC = new JunctionBC(junction, -1, this);

            QSharedPointer<TimeSeries> timeSeries = QSharedPointer<TimeSeries>(new TimeSeries(QUuid::createUuid().toString(),1));
            m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
            timeSeries->addRow(m_startDateTime - 0.1, value);
            timeSeries->addRow(m_endDateTime + 0.1, value);
            junctionBC->setTimeSeries(timeSeries);

            m_boundaryConditions.push_back(junctionBC);

            found = true;
          }
          else
          {
            errorMessage = "Temperature BC value is invalid";
            return false;
          }
        }
        else if (!QString::compare(variable, "FLOW", Qt::CaseInsensitive))
        {
          bool valueOk ;
          double value = columns[3].toDouble(&valueOk);

          if (valueOk)
          {
            JunctionBC *junctionBC = new JunctionBC(junction, -2, this);

            QSharedPointer<TimeSeries> timeSeries = QSharedPointer<TimeSeries>(new TimeSeries(QUuid::createUuid().toString(),1));
            m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
            timeSeries->addRow(m_startDateTime - 0.1, value);
            timeSeries->addRow(m_endDateTime + 0.1, value);
            junctionBC->setTimeSeries(timeSeries);

            m_boundaryConditions.push_back(junctionBC);

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
                JunctionBC *junctionBC = new JunctionBC(junction, i, this);

                QSharedPointer<TimeSeries> timeSeries = QSharedPointer<TimeSeries>(new TimeSeries(QUuid::createUuid().toString(),1));
                m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
                timeSeries->addRow(m_startDateTime - 0.1, value);
                timeSeries->addRow(m_endDateTime + 0.1, value);
                junctionBC->setTimeSeries(timeSeries);

                m_boundaryConditions.push_back(junctionBC);

                found = true;
              }
            }

            if(found)
              break;
          }
        }

        if (!found)
        {
          return false;
        }
      }
      else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
      {
        int variableIndex = -3;

        if (!QString::compare(variable, "TEMPERATURE", Qt::CaseInsensitive))
        {
          variableIndex = -1;
        }
        else if (!QString::compare(variable, "FLOW", Qt::CaseInsensitive))
        {
          variableIndex = -2;
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

        if (variableIndex > -3)
        {
          std::string tsId = columns[3].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            JunctionBC *junctionBC = new JunctionBC(junction, variableIndex, this);

            junctionBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(junctionBC);
          }
          else
          {
            errorMessage = "Specified timeseries does not exist";
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

bool CSHModel::readInputFileSourcesTag(const QString &line, QString &errorMessage)
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
      SourceBC::VariableType variableType;

      int soluteIndex = -1;

      if (!QString::compare(variable, "HEAT", Qt::CaseInsensitive))
      {
        variableType = SourceBC::HeatSource;
        soluteIndex = 0;
      }
      else if (!QString::compare(variable, "FlOW", Qt::CaseInsensitive))
      {
        variableType = SourceBC::FlowSource;
        soluteIndex = 0;
      }
      else
      {

        for (size_t i = 0; i < m_solutes.size(); i++)
        {
          std::string solute = m_solutes[i];

          if (!solute.compare(variable.toStdString()))
          {
            soluteIndex = i;
            variableType = SourceBC::SoluteSource;
            break;
          }
        }
      }

      if(soluteIndex > -1)
      {
        if (!QString::compare(type, "VALUE", Qt::CaseInsensitive))
        {
          bool valueOk;
          double value = columns[6].toDouble(&valueOk);

          if (valueOk)
          {
            SourceBC *nonPointSrcTSBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            nonPointSrcTSBC->setSoluteIndex(soluteIndex);

            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            nonPointSrcTSBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(nonPointSrcTSBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
        else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = columns[6].toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            SourceBC *nonPointSrcTSBC = new SourceBC(elementFrom, startFactor, elementTo, endFactor, variableType, this);
            nonPointSrcTSBC->setSoluteIndex(soluteIndex);
            nonPointSrcTSBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(nonPointSrcTSBC);
          }
          else
          {
            errorMessage = "Source is is invalid";
            return false;
          }
        }
      }
      else
      {
        errorMessage = "Source is is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Source is is invalid";
      return false;
    }
  }
  else
  {
    return false;
  }

  return true;
}

bool CSHModel::readInputFileHydraulicsTag(const QString &line, QString &errorMessage)
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
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement,
                                                          variableIndex, this);
            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);

            m_timeSeries[ts->id().toStdString()] = ts;

            hydraulicsBC->setTimeSeries(ts);

            m_boundaryConditions.push_back(hydraulicsBC);
          }
          else
          {
            errorMessage = "Uniform hydraulics value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if (tsIt != m_timeSeries.end())
          {
            HydraulicsBC *hydraulicsBC = new HydraulicsBC(fromElement, toElement,
                                                          variableIndex, this);


            hydraulicsBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(hydraulicsBC);
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

bool CSHModel::readInputFileRadiativeFluxesTag(const QString &line, QString &errorMessage)
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
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);

          QUuid uid = QUuid::createUuid();
          QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

          ts->addRow(m_startDateTime, value);
          ts->addRow(m_endDateTime, value);

          m_timeSeries[ts->id().toStdString()] = ts;

          radiationFluxBC->setTimeSeries(ts);

          m_boundaryConditions.push_back(radiationFluxBC);
        }
        else
        {
          errorMessage = "Radiation BC value is invalid";
          return false;
        }
      }
      else if (!QString::compare(type, "TIMESERIES", Qt::CaseInsensitive))
      {
        std::string tsId = varValue.toStdString();
        auto tsIt = m_timeSeries.find(tsId);

        if(tsIt != m_timeSeries.end())
        {
          RadiativeFluxBC *radiationFluxBC = new RadiativeFluxBC(fromElement, toElement, this);
          radiationFluxBC->setTimeSeries(tsIt->second);
          m_boundaryConditions.push_back(radiationFluxBC);
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

bool CSHModel::readInputFileMeteorologyTag(const QString &line, QString &errorMessage)
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
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement,
                                                             variableIndex, this);
            QUuid uid = QUuid::createUuid();
            QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));
            ts->addRow(m_startDateTime, value);
            ts->addRow(m_endDateTime, value);
            m_timeSeries[ts->id().toStdString()] = ts;

            meteorologyBC->setTimeSeries(ts);
            m_boundaryConditions.push_back(meteorologyBC);
          }
          else
          {
            errorMessage = "Uniform meteorology value is invalid";
            return false;
          }
        }
        else if (!QString::compare(valueType, "TIMESERIES", Qt::CaseInsensitive))
        {
          std::string tsId = varValue.toStdString();
          auto tsIt = m_timeSeries.find(tsId);

          if(tsIt != m_timeSeries.end())
          {
            MeteorologyBC *meteorologyBC = new MeteorologyBC(fromElement, toElement,
                                                             variableIndex, this);

            meteorologyBC->setTimeSeries(tsIt->second);
            m_boundaryConditions.push_back(meteorologyBC);
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

bool CSHModel::readInputFileTimeSeriesTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    QFileInfo fileInfo(options[1].trimmed());

    if (fileInfo.isRelative())
      fileInfo = relativePathToAbsolute(fileInfo);

    if(QFile::exists(fileInfo.absoluteFilePath()))
    {
      QSharedPointer<TimeSeries> timeSeries(TimeSeries::createTimeSeries(options[0], fileInfo, this));

      if(!timeSeries.isNull())
      {
        m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
      }
      else
      {
        errorMessage = "Timeseries specified is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified filepath does not exist";
      return false;
    }
  }
  else
  {
    errorMessage = "TimeSeries must have two columns";
    return false;
  }

  return true;
}

void CSHModel::writeOutput()
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

void CSHModel::writeCSVOutput()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    for (size_t i = 0; i < m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      m_outputCSVStream << m_currentDateTime << ", " << QString::fromStdString(element->id) << ", " << element->tIndex
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

void CSHModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF
  if (m_outputNetCDF)
  {

    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();

    //Set current dateTime
    //      ThreadSafeNcVar timeVar =  m_outputNetCDF->getVar("time");
    //      timeVar.putVar(std::vector<size_t>({currentTime}), m_currentDateTime);
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);


    float *flow = new float[m_elements.size()];
    float *depth = new float[m_elements.size()];
    float *width = new float[m_elements.size()];
    float *xsectArea = new float[m_elements.size()];
    float *dispersion = new float[m_elements.size()];
    float *temperature = new float[m_elements.size()];
    float *dvolumedt = new float[m_elements.size()];
    float *waterAge = new float[m_elements.size()]();
    float *elementEvapHeatFlux = new float[m_elements.size()];
    float *elementConvHeatFlux = new float[m_elements.size()];
    float *elementRadiationFlux = new float[m_elements.size()];
    float *elementHeatFlux = new float[m_elements.size()];
    float *elementAirTemp =  new float[m_elements.size()];
    float *elementRH =  new float[m_elements.size()];
    float *elementWindSpeed =  new float[m_elements.size()];
    float *elementVaporPress =  new float[m_elements.size()];
    float *elementSatVaporPress =  new float[m_elements.size()];
    float *elementAirVaporPress =  new float[m_elements.size()];
    float *elementAirSatVaporPress =  new float[m_elements.size()];
    float *totalElementHeatBalance = new float[m_elements.size()];
    float *totalElementAdvDispHeatBalance = new float[m_elements.size()];
    float *totalElementEvapHeatBalance = new float[m_elements.size()];
    float *totalElementConvHeatBalance = new float[m_elements.size()];
    float *totalElementRadiationFluxHeatBalance = new float[m_elements.size()];
    float *totalElementExternalHeatFluxBalance = new float[m_elements.size()];
    float *solutes = new float[m_elements.size() * m_numSolutes];
    float *totalElementSoluteMassBalance = new float[m_elements.size() * m_numSolutes];
    float *totalElementAdvDispSoluteMassBalance = new float[m_elements.size() * m_numSolutes];
    float *totalElementExternalSoluteFluxMassBalance = new float[m_elements.size() * m_numSolutes];


    if(m_simulateWaterAge)
    {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)m_elements.size(); i++)
      {
        Element *element = m_elements[i];
        waterAge[i] = element->soluteConcs[m_solutes.size() - 1].value;
      }

      m_outNetCDFVariables["water_age"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), waterAge);
    }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      flow[i] = element->flow.value ;
      depth[i] = element->depth;
      width[i] = element->width;
      xsectArea[i] = element->xSectionArea;
      dispersion[i] = element->longDispersion.value;
      temperature[i] = element->temperature.value;
      dvolumedt[i] = element->dvolume_dt.value;
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

      for (int j = 0; j < m_numSolutes; j++)
      {
        solutes[i + j * m_elements.size()] = element->soluteConcs[j].value;
        totalElementSoluteMassBalance[i + j * m_elements.size()] = element->totalSoluteMassBalance[j];
        totalElementAdvDispSoluteMassBalance[i + j * m_elements.size()] = element->totalAdvDispSoluteMassBalance[j];
        totalElementExternalSoluteFluxMassBalance[i + j * m_elements.size()] = element->totalExternalSoluteFluxesMassBalance[j];
      }
    }

    m_outNetCDFVariables["flow"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), flow);

    m_outNetCDFVariables["depth"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), depth);

    m_outNetCDFVariables["width"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), width);

    m_outNetCDFVariables["xsection_area"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), xsectArea);

    m_outNetCDFVariables["dispersion"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), dispersion);

    m_outNetCDFVariables["temperature"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), temperature);

    m_outNetCDFVariables["volume_time_derivative"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), dvolumedt);

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

    m_outNetCDFVariables["total_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalHeatBalance);

    m_outNetCDFVariables["total_adv_disp_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalAdvDispHeatBalance);

    m_outNetCDFVariables["total_evap_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalEvaporationHeatBalance);

    m_outNetCDFVariables["total_conv_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalConvectiveHeatBalance);

    m_outNetCDFVariables["total_radiation_flux_heat_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalRadiationHeatBalance);

    m_outNetCDFVariables["total_external_heat_flux_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalExternalHeatFluxBalance);

    if(m_numSolutes)
    {
      m_outNetCDFVariables["solute_concentration"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size()}), solutes);

      m_outNetCDFVariables["total_element_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size()}), totalElementSoluteMassBalance);

      m_outNetCDFVariables["total_element_adv_disp_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size()}), totalElementAdvDispSoluteMassBalance);

      m_outNetCDFVariables["total_element_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size()}), totalElementExternalSoluteFluxMassBalance);

      m_outNetCDFVariables["total_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numSolutes}), m_totalSoluteMassBalance.data());

      m_outNetCDFVariables["total_adv_disp_solute_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numSolutes}), m_totalAdvDispSoluteMassBalance.data());

      m_outNetCDFVariables["total_external_solute_flux_mass_balance"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numSolutes}), m_totalExternalSoluteFluxMassBalance.data());
    }

    delete[] flow;
    delete[] depth;
    delete[] width;
    delete[] xsectArea;
    delete[] dispersion;
    delete[] temperature;
    delete[] dvolumedt;
    delete[] waterAge;
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

void CSHModel::closeOutputFiles()
{
  closeCSVOutputFile();
  closeOutputNetCDFFile();
}

void CSHModel::closeCSVOutputFile()
{
  if (m_outputCSVStream.device() && m_outputCSVStream.device()->isOpen())
  {
    m_outputCSVStream.flush();
    m_outputCSVStream.device()->close();
    delete m_outputCSVStream.device();
    m_outputCSVStream.setDevice(nullptr);
  }
}

void CSHModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF

  if(m_outputNetCDF)
  {
    m_outputNetCDF->sync();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }

#endif
}

QFileInfo CSHModel::relativePathToAbsolute(const QFileInfo &fileInfo)
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


const unordered_map<string, int> CSHModel::m_inputFileFlags({
                                                              {"[OPTIONS]", 1},
                                                              {"[OUTPUTS]", 2},
                                                              {"[SOLUTES]", 3},
                                                              {"[ELEMENTJUNCTIONS]", 4},
                                                              {"[ELEMENTS]", 5},
                                                              {"[ELEMENT_HYDRAULIC_VARIABLES]", 12},
                                                              {"[BOUNDARY_CONDITIONS]", 6},
                                                              {"[SOURCES]", 7},
                                                              {"[HYDRAULICS]", 8},
                                                              {"[RADIATIVE_FLUXES]", 9},
                                                              {"[METEOROLOGY]", 10},
                                                              {"[TIMESERIES]", 11}
                                                            });

const unordered_map<string, int> CSHModel::m_optionsFlags({
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
                                                            {"TVD_FLUX_LIMITER", 26},
                                                            {"NETCDF_THREAD_SAFE", 27},
                                                            {"SIMULATE_WATER_AGE", 28},
                                                            {"SOLVER", 29},
                                                            {"SOLVER_ABS_TOL", 30},
                                                            {"SOLVER_REL_TOL", 31},
                                                            {"LINEAR_SOLVER", 32},
                                                            {"SOLVE_HYDRAULICS", 33},
                                                          });

const unordered_map<string, int> CSHModel::m_advectionFlags({
                                                              {"UPWIND", 1},
                                                              {"CENTRAL", 2},
                                                              {"HYBRID", 3},
                                                              {"TVD", 4},
                                                            });

const unordered_map<string, int> CSHModel::m_solverTypeFlags({{"RK4", 1},
                                                              {"RKQS", 2},
                                                              {"ADAMS", 3},
                                                              {"BDF", 4},
                                                              {"EULER", 5}
                                                             });

const unordered_map<string, int> CSHModel::m_linearSolverTypeFlags({{"GMRES", 1},
                                                                    {"FGMRES", 2},
                                                                    {"Bi_CGStab", 3},
                                                                    {"TFQMR", 4},
                                                                    {"PCG", 5}
                                                                   });

const unordered_map<string, int> CSHModel::m_hydraulicVariableFlags({{"DEPTH", 1},
                                                                     {"WIDTH", 2},
                                                                     {"XSECTION_AREA", 3},
                                                                     {"FLOW", 4}});

const unordered_map<string, int> CSHModel::m_meteorologicalVariableFlags({{"RELATIVE_HUMIDITY", 1},
                                                                          {"AIR_TEMPERATURE", 2},
                                                                          {"WIND_SPEED", 3}});

const QRegExp CSHModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");
