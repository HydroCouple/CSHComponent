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
typedef double (Element::*ComputeTempAdv)(double dt, double T[]);

/*!
 *\brief Function pointer to calculate solute advection to eliminate costly if else function calls
 */
typedef double (Element::*ComputeSoluteAdv)(double dt, double S[], int soluteIndex);

/*!
 * \brief This struct represents the channel control volume
 */
struct CSHCOMPONENT_EXPORT Element
{
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
    * \brief index unique identifier for element
    */
   int index;


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
    * \brief flow (m^3/s)
    */
   double flow;

   /*!
    * \brief slope
    */
   double slope;


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
    * \brief computeDTDt - Computes the time derivative of temperature based on data generated by the ODE solver.
    * \param dt - The timestep over which to compute the solute gradient.
    * \param T - The temperature array for all elements.
    * \return
    */
   double computeDTDt(double dt, double T[]);

   /*!
    * \brief computeDTDtUpwind
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtUpwind(double dt, double T[]);

   /*!
    * \brief computeDTDtCentral
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtCentral(double dt, double T[]);

   /*!
    * \brief computeDTDtHybrid
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtHybrid(double dt, double T[]);

   /*!
    * \brief computeDTDtUltimate
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtTVD(double dt, double T[]);


   /*!
    * \brief computeTVDLimiter
    * \param r
    * \param limiter
    * \return
    */
   double computeTVDLimiter(double r, int limiter);

   /*!
    * \brief computeDTDtULTIMATE
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtULTIMATE(double dt, double T[]);

   /*!
    * \brief computeDTDtDispersion
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtDispersion(double dt, double T[]);

   /*!
    * \brief computeEvaporation
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtEvaporation(double dt, double T[]);

   /*!
    * \brief computeConvection
    * \param dt
    * \param T
    * \return
    */
   double computeDTDtConvection(double dt, double T[]);

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
   double computeDSoluteDtDispersion(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtUpwind
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtUpwind(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtCentral
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtCentral(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtHybrid
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtHybrid(double dt, double S[], int soluteIndex);

   /*!
    * \brief computeDSoluteDtENO
    * \param dt
    * \param S
    * \param soluteIndex
    * \return
    */
   double computeDSoluteDtTVD(double dt, double S[], int soluteIndex);

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
    * \brief computeDVolumeDt
    */
   void computeDVolumeDt();

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

  private:

   /*!
    * \brief initializeSolutes
    */
   void initializeSolutes();


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
   ComputeTempAdv computeTempAdv;

   /*!
    * \brief computeSoluteAdv Pointer to function to compute solute advection.
    */
   ComputeSoluteAdv computeSoluteAdv;

};

#endif // ELEMENT_H
