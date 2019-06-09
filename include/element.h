/*!
*  \file    element.h
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

#ifndef ELEMENT_H
#define ELEMENT_H


#include "variable.h"
#include "cshcomponent_global.h"

#include <string>

struct Element;
struct ElementJunction;
class CSHModel;


/*!
 *\brief Function pointer to calculate temperature advection to eliminate costly if else function calls
 */
typedef double (Element::*ComputeTempDeriv)(double dt, double T[]);

/*!
 *
 */
typedef double (*ComputeTempAdvDeriv)(Element *element, double dt, double T[]);

/*!
 *\brief Function pointer to calculate solute advection to eliminate costly if else function calls
 */
typedef double (Element::*ComputeSoluteDeriv)(double dt, double S[], int soluteIndex);

typedef double (Element::*GetXofY)(double y);

/*!
 *
 */
typedef double (*ComputeSoluteAdvDeriv)(Element *element, double dt, double S[], int soluteIndex);

/*!
 *
 */
typedef void (*SetAdvectionFuctions)(Element *element);

/*!
 * \brief This struct represents the channel control volume
 */
struct CSHCOMPONENT_EXPORT Element
{
   friend struct ElementJunction;
   friend class ElementAdvUpwind;
   friend class ElementAdvCentral;
   friend class ElementAdvHybrid;
   friend class ElementAdvQUICK;
   friend class ElementAdvULTIMATE;
   friend class ElementAdvTVD;
   friend class CSHModel;

   enum XSectType
   {
     RECT,
     TRAP,
     IRREGULAR
   };

    /*!
    * \brief Element - Creates an instance of the control volume element used to represent a computational
    * element in a reach.
    * \param numSolutes - Number of solutes that are going to be transported in addition to temperature.
    * \param from - The upstream junction of this element.
    * \param to - The downstream junction of this element.
    * \param project
    */
   Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  CSHModel *model);

   /*!
    * \brief ~Element - Destructor for this class.
    */
   ~Element();

   /*!
    * \brief hIndex
    */
   int hIndex;

   /*!
    * \brief index unique identifier for element
    */
   int tIndex;

   /*!
    * \brief sIndex
    */
   int *sIndex;

   /*!
    * \brief id
    */
   std::string id;

   /*!
    * \brief x
    */
   double x;

   /*!
    * \brief y
    */
   double y;

   /*!
    * \brief z
    */
   double z;

   /*!
    * \brief temperature (°C)
    */
   Variable temperature;

   /*!
    * \brief prevTemperature (°C)
    */
   Variable prevTemperature;

   /*!
    * \brief numSolutes
    */
   int numSolutes = 0;

   /*!
    * \brief soluteConcs (kg/m^3)
    */
   Variable *soluteConcs;

   /*!
    * \brief prevSoluteConcs (kg/m^3)
    */
   Variable *prevSoluteConcs;

   /*!
    * \brief longDispersion (m^2/s)
    */
   Variable longDispersion;

   /*!
    * \brief fromJunction
    */
   ElementJunction *upstreamJunction;

   /*!
    * \brief toJunction
    */
   ElementJunction *downstreamJunction;

   /*!
    * \brief length (m)
    */
   double length;

   /*!
    * \brief depth (m)
    */
   double depth;

   /*!
    * \brief xSectionArea (m^2)
    */
   double xSectionArea;

   /*!
    * \brief prevXSectionArea
    */
//   double prevXSectionArea;

   /*!
    * \brief upstreamXSectionArea
    */
   double upstreamXSectionArea;

   /*!
    * \brief downstreamXSectionArea
    */
   double downstreamXSectionArea;

   /*!
    * \brief width (m)
    */
   double width;

   /*!
    * \brief bottomWidth
    */
   double bottomWidth;

   /*!
    * \brief flow (m^3/s)
    */
   Variable flow;

   /*!
    * \brief prevFlow
    */
   Variable prevFlow;

   /*!
    * \brief externalFlows
    */
   double externalFlows;

   /*!
    * \brief slope
    */
   double slope;

   /*!
    * \brief xsectionType
    */
   XSectType xsectionType = RECT;

   /*!
    * \brief mannings
    */
   double mannings = 0.039896;

   /*!
    * \brief sideSlopes
    */
   double *sideSlopes;

   /*!
    * \brief relativeHumidity (%)
    */
   double relativeHumidity;

   /*!
    * \brief windFunction
    */
   double windFunction;

   /*!
    * \brief evaporationRate (m/s)
    */
   double evaporationRate;

   /*!
    * \brief evaporationHeatFlux J/s
    */
   double evaporationHeatFlux;

   /*!
    * \brief saturationVaporPressure (kPa)
    */
   double saturationVaporPressureAir;

   /*!
    * \brief saturationVaporPressureWater (kPa)
    */
   double saturationVaporPressureWater;

   /*!
    * \brief vaporPressure (kPa)
    */
   double vaporPressureAir;

   /*!
    * \brief vaporPressureWater (kPa)
    */
   double vaporPressureWater;

   /*!
    * \brief windVelocity (m/s)
    */
   double windSpeed;

   /*!
    * \brief airTemperature (°C)
    */
   double airTemperature;

   /*!
    * \brief convectionHeatFlux in units of (J/s)
    */
   double convectionHeatFlux;

   /*!
    * \brief externalHeatFluxes in units of (J/s)
    */
   double externalHeatFluxes;

   /*!
    * \brief externalSoluteFluxes (W/m^2)
    */
   double radiationFluxes;

   /*!
    * \brief fluidFrictionHeatFlux
    */
   double fluidFrictionHeatFlux;

   /*!
    * \brief externalSoluteFluxes of the form (kg/s)
    */
   double *externalSoluteFluxes;

   /*!
    * \brief heatBalance (KJ)
    */
   double totalHeatBalance;

   /*!
    * \brief totalAdvDispHeatBalance (KJ)
    */
   double totalAdvDispHeatBalance;

   /*!
    * \brief totalRadiationHeatBalance  (KJ)
    */
   double totalRadiationFluxesHeatBalance;

   /*!
    * \brief totalExternalHeatFluxesBalance  (KJ)
    */
   double totalExternalHeatFluxesBalance;

   /*!
    * \brief totalEvaporativeHeatFluxesBalance  (KJ)
    */
   double totalEvaporativeHeatFluxesBalance;

   /*!
    * \brief totalConvectiveHeatFluxesBalance  (KJ)
    */
   double totalConvectiveHeatFluxesBalance;

   /*!
    * \brief soluteMassBalance  (kg)
    */
   double *totalSoluteMassBalance;

   /*!
    * \brief totalAdvDispSoluteMassBalance  (kg)
    */
   double *totalAdvDispSoluteMassBalance;

   /*!
    * \brief totalExternalSoluteFluxesMassBalance  (kg)
    */
   double *totalExternalSoluteFluxesMassBalance;

   /*!
    * \brief pecletNumber
    */
   double pecletNumber;

   /*!
    * \brief upstreamElement
    */
   Element *upstreamElement;

   /*!
    * \brief upstreamElementDirection
    */
   double upstreamElementDirection;

   /*!
    * \brief downstreamElement
    */
   Element *downstreamElement;

   /*!
    * \brief downstreamElementDirection
    */
   double downstreamElementDirection;

   /*!
    * \brief distanceFromUpStreamJunction
    */
   double distanceFromUpStreamJunction;

   /*!
    * \brief upstreamCourantNumber
    */
   double upstreamCourantNumber;

   /*!
    * \brief downstreamCourantNumber
    */
   double downstreamCourantNumber;

   /*!
    * \brief dvolume_dt
    */
   Variable dvolume_dt;

   /*!
    * \brief model
    */
   CSHModel *model;

   /*!
    * \brief initializeSolutes
    * \param numSolutes
    */
   void initialize();

   /*!
    * \brief computeDADt
    * \param dt
    * \param A
    * \return
    */
   double computeDADt(double dt, double A[]);

   void calculateQfromA(double A[]);

   double getAofH(double H);

   double getPofH(double H);

   double getWofH(double H);

   double getHofA(double A);

   double getQofH(double H);

   double getAofQ(double Q);

   double getHofQ(double Q);

   double findRoots(double x, double y, GetXofY function, int maxIters = 10000, double derivStepSize = 1e-8, double eps = 1e-10);

   void computeHydraulicVariables();

   /*!
    * \brief computeDTDt - Computes the time derivative of temperature based on data generated by the ODE solver.
    * \param dt - The timestep over which to compute the solute gradient.
    * \param T - The temperature array for all elements.
    * \return
    */
   double computeDTDt(double dt, double T[]);

   /*!
    * \brief computeDTDtAdv
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtAdv(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersion
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersion(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionUpstreamJunction
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionUpstreamJunction(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionUpstreamJunctionBC
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionUpstreamJunctionBC(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionUpstreamNeighbour
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionUpstreamNeighbour(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionDownstreamJunction
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionDownstreamJunction(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionDownstreamJunctionBC
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionDownstreamJunctionBC(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionDownstreamNeighbour
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionDownstreamNeighbour(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersionDownstreamNeighbour
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersionSelf(double dt, double T[]);

   /*!
    * \brief computeEvaporation
    * \param dt
    * \param T
    * \return
    */
   void computeEvaporation();

   /*!
    * \brief computeConvection
    * \param dt
    * \param T
    * \return
    */
   void computeConvection();

   /*!
    * \brief computeFluidFrictionHeat
    */
   void computeFluidFrictionHeat();

   /*!
    * \brief computeDSoluteDt
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDt(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersion
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtAdv(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersion
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersion(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDTDtDispersionUpstreamJunction
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionUpstreamJunction(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersionUpstreamJunctionBC
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionUpstreamJunctionBC(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersionUpstreamNeighbour
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionUpstreamNeighbour(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersionDownstreamJunction
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionDownstreamJunction(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersionDownstreamJunctionBC
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionDownstreamJunctionBC(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtDispersionDownstreamNeighbour
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtDispersionDownstreamNeighbour(double dt, double S[], int soluteIndex);

   double computeDSoluteDtDispersionSelf(double dt, double S[], int soluteIndex);

   /*!
    * \brief setDispersionFunctions
    */
   void setDispersionFunctions();

   /*!
    * \brief computeCourantNumber
    * \return
    */
   double computeCourantFactor() const;

   /*!
    * \brief commputeDispersionFactor
    * \return
    */
   double computeDispersionFactor() const;

   /*!
    * \brief computeDerivedHydraulics
    */
   void computeDerivedHydraulics();

   /*!
    * \brief computeLongDispersion
    */
   void computeLongDispersion();

   /*!
    * \brief computePecletNumbers
    */
   void computePecletNumbers();

   /*!
    * \brief calculateUpstreamPeclet
    */
   void computeUpstreamPeclet();

   /*!
    * \brief calculateDownstreamPeclet
    */
   void computeDownstreamPeclet();

   /*!
    * \brief computeHeatBalance
    */
   void computeHeatBalance(double timeStep);

   /*!
    * \brief computeSoluteBalance
    * \param soluteIndex
    */
   void computeSoluteBalance(double timeStep, int soluteIndex);

   /*!
    * \brief initializeSolutes
    */
   void initializeSolutes();

  private:

   /*!
    * \brief setUpstreamElement
    */
   void setUpstreamElement();

   /*!
    * \brief setDownStreamElement
    */
   void setDownStreamElement();

   /*!
    * \brief computeUpstreamFlow
    */
   void computeUpstreamFlow();

   /*!
    * \brief computeDownstreamFlow
    */
   void computeDownstreamFlow();

   /*!
    * \brief deleteSoluteVariables
    */
   void deleteSoluteVariables();

   /*!
    * \brief downstreamDispersion
    */
   double downstreamLongDispersion;

   /*!
    * \brief downStreamPecletNumber
    */
   double downstreamPecletNumber;

   /*!
    * \brief downstreamFlow
    */
   double downstreamFlow;

   /*!
    * \brief downstreamVelocity
    */
   double downstreamVelocity;

   /*!
    * \brief upstreamDispersion
    */
   double upstreamLongDispersion;

   /*!
    * \brief upStreamPecletNumber
    */
   double upstreamPecletNumber;

   /*!
    * \brief upstreamFlow
    */
   double upstreamFlow;

   /*!
    * \brief upstreamVelocity
    */
   double upstreamVelocity;

   /*!
    * \brief rho_cp
    */
   double rho_cp;

   /*!
    * \brief rho_cp_vol
    */
   double rho_cp_vol;

   /*!
    * \brief rho_vol
    */
   double rho_vol;

   /*!
    * \brief top_area
    */
   double top_area;

   /*!
    * \brief volume
    */
   double volume;

   /*!
    * \brief prev_volume
    */
   double prev_volume;

   /*!
    * \brief starting
    */
   bool starting;

   /*!
    * \brief computeTempAdv Pointer to function to compute temperature advection.
    */
   ComputeTempAdvDeriv *computeTempAdvDeriv;

   /*!
    * \brief computeSoluteAdv Pointer to function to compute solute advection.
    */
   ComputeSoluteAdvDeriv **computeSoluteAdvDeriv;

   /*!
    * \brief computeTempDispDeriv
    */
   ComputeTempDeriv *computeTempDispDeriv;

   /*!
    * \brief computeTempSoluteDeriv
    */
   ComputeSoluteDeriv **computeSoluteDispDeriv;

   /*!
    * \brief setAdvectionFuctions
    */
   SetAdvectionFuctions setAdvectionFuctions;

};

#endif // ELEMENT_H
