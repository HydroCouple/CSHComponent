/*!
*  \file    elementjunction.h
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

#ifndef ELEMENJUNCTION_H
#define ELEMENJUNCTION_H

#include "variable.h"
#include "cshcomponent_global.h"

#include <string>
#include <set>

class CSHModel;
class HCVertex;
struct Element;

/*!
 * \brief The ElementJunction struct represents the
 */
struct CSHCOMPONENT_EXPORT ElementJunction
{

    enum JunctionType
    {
      NoElement = 0,
      SingleElement = 1,
      DoubleElement = 2,
      MultiElement = 3
    };

    /*!
     * \brief ElementJunction
     * \param numsolutes - Number of solutes
     * \param model -
     */
    ElementJunction(const std::string &id, double x, double y, double z, CSHModel *model);

    /*!
     * \brief ~ElementJunction - Deletes the ElementJunction and its associated data.
     */
    ~ElementJunction();

    /*!
     * \brief id unique identifier for this junction
     */
    std::string id;

    /*!
     * \brief hIndex
     */
    int hIndex;

    /*!
     * \brief index
     */
    int tIndex;

    /*!
     * \brief sIndex
     */
    int *sIndex;

    /*!
     * \brief x location
     */
    double x;

    /*!
     * \brief y location
     */
    double y;

    /*!
     * \brief z location
     */
    double z;

    /*!
     * \brief numSolutes
     */
    int numSolutes;

    /*!
     * \brief temperature
     */
    Variable temperature;

    /*!
     * \brief prevTemperature
     */
    Variable prevTemperature;

    /*!
     * \brief soluteConcs
     */
    Variable *soluteConcs;

    /*!
     * \brief prevSoluteConcs
     */
    Variable *prevSoluteConcs;

    /*!
     * \brief incomingElements
     */
    std::set<Element*> incomingElements;

    /*!
     * \brief outgoingElements
     */
    std::set<Element*> outgoingElements;

    /*!
     * \brief junctionType
     */
    JunctionType junctionType;

    /*!
     * \brief flow
     */
    Variable inflow;

    /*!
     * \brief model
     */
    CSHModel *model;

    bool starting = true;

    double volume;

    double prev_volume;

    double dvolume_dt;

    void initialize();

    /*!
     * \brief interpTemp
     */
    void interpTemp();

    /*!
     * \brief interpSoluteConcs
     */
    void interpSoluteConcs(int soluteIndex);

    /*!
     * \brief solveHeatContinuity
     * \param dt
     */
    void solveHeatContinuity(double dt);

    /*!
     * \brief solveSoluteContinuity
     * \param soluteIndex
     * \param dt
     * \return
     */
    void solveSoluteContinuity(int soluteIndex, double dt);

    /*!
     * \brief computeDTDt
     * \param dt
     * \param T
     * \return
     */
    double computeDTDt(double dt, double T[]);

    /*!
     * \brief computeDSoluteDt
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDt(double dt, double S[], int soluteIndex);

    /*!
     * \brief copyVariablesToPrev
     */
    void copyVariablesToPrev();

    /*!
     * \brief initializeSolutes
     * \param numSolutes
     */
    void initializeSolutes();

    /*!
     * \brief computeDerivedHydraulics
     */
    void computeDerivedHydraulics();



};

#endif // ELEMENJUNCTION_H
