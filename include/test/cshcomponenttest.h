/*!
*  \file    CSHComponenttest.h
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


#ifndef CSHCOMPONENTTEST_H
#define CSHCOMPONENTTEST_H

#include <QtTest/QtTest>

class CSHComponentTest : public QObject
{
    Q_OBJECT

  private slots:

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

};


#endif // CSHCOMPONENTTEST_H
