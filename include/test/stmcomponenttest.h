/*!
*  \file    stmcomponenttest.h
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


#ifndef STMCOMPONENTTEST_H
#define STMCOMPONENTTEST_H

#include <QtTest/QtTest>

class STMComponentTest : public QObject
{
    Q_OBJECT

  private slots:

    /*!
     * \brief solveODERK4_Prob1 Solve ODE problem 1 using RK4
     */
    void solveODERK4_Prob1();

    /*!
     * \brief solveODERK4_Prob1 Solve ODE problem 2 using RK45
     */
    void solveODERKQS_Prob1();

#ifdef USE_CVODE

    /*!
     * \brief solveODEAdam_Prob1
     */
    void solveODEAdam_Prob1();

    /*!
     * \brief solveODEBDF_Prob1
     */
    void solveODEBDF_Prob1();

#endif

    /*!
     * \brief solveODERK4_Prob2 Solve ODE problem 1 using RK4
     */
    void solveODERK4_Prob2();

    /*!
     * \brief solveODERKQS_Prob2 Solve ODE problem 1 using RK45
     */
    void solveODERKQS_Prob2();

#ifdef USE_CVODE

    /*!
     * \brief solveODEAdam_Prob2
     */
    void solveODEAdam_Prob2();

    /*!
     * \brief solveODEBDF_Prob2
     */
    void solveODEBDF_Prob2();

#endif

    /*!
     * \brief versteegCase1_Upwind Solve steady state problem Case 1 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the upwind differencing method for advection discretization.
     */
    void versteegCase1_Upwind();

    /*!
     * \brief versteegCase1_Central Solve steady state problem Case 1 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the central differencing method for advection discretization.
     */
    void versteegCase1_Central();

    /*!
     * \brief versteegCase1_Hybrid Solve steady state problem Case 1 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the hybrid differencing method for advection discretization.
     */
    void versteegCase1_Hybrid();

    /*!
     * \brief versteegCase2_Upwind Solve steady state problem Case 2 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the upwind differencing method for advection discretization.
     */
    void versteegCase2_Upwind();

    /*!
     * \brief versteegCase2_Central Solve steady state problem Case 2 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the central differencing method for advection discretization.
     */
    void versteegCase2_Central();

    /*!
     * \brief versteegCase2_Hybrid Solve steady state problem Case 2 Page 137 (Versteeg, H.K. and W. Malalasekera, 2007.
     * An Introduction to Computational Fluid Dynamics: The Finite Volume Method. Pearson Education Ltd., Harlow, England; New York.)
     * in using the hybrid differencing method for advection discretization.
     */
    void versteegCase2_Hybrid();

    /*!
     * \brief green_river_test
     */
    void green_river_test();

    /*!
     * \brief green_river_test
     */
    void green_river_test1();

    /*!
     * \brief green_river_test1
     */
    void green_river_test2();

    /*!
     * \brief derivativeProb1 Example ODE problem: dy/dt = x * y ^3 / sqrt(1 + x^2); y(0) = -1; y = -1 / sqrt(3 - 2 * sqrt(1+t^2))
     * \param t
     * \param y
     * \param dydt
     * \param userData
     */
    static void derivativeProb1(double t, double y[], double dydt[], void* userData);

    /*!
     * \brief problem1  y = -1 / sqrt(3 - 2 * sqrt(1+t^2))
     * \param t
     * \return
     */
    static double problem1(double t);

    /*!
     * \brief derivativeProb2 Example ODE problem: dy/dt = (3 * t^2 + 4 * t - 4)/(2 * y - 4), y(1) = 3, y = 2 + sqrt(t^3 + 2 * t^2-4t+2)
     * \param t time
     * \param y distance
     * \param dydt
     * \param userData optional user supplied data.
     */
    static void derivativeProb2(double t, double y[], double dydt[], void* userData);

    /*!
     * \brief problem2 y = - (t^2 * y) / (2.0 * sqrt(2.0 - y^2)
     * \param t
     * \return
     */
    static double problem2(double t);

};


#endif // STMCOMPONENTTEST_H
