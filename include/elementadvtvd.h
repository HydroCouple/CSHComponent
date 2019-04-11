/*!
 *  \file    CSHComponent.h
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
 *  \todo fix flow in opposite direction to account for boundary conditions
 *  \warning
 */

#ifndef ELEMENTADVTVD_H
#define ELEMENTADVTVD_H

struct Element;

class ElementAdvTVD
{
  public:

    enum TVDFluxLimiter
    {
      MIN_MOD = 0,
      SUPERBEE = 1,
      VAN_LEER = 2,
      MUSCL = 3,
      SWEBY = 4,
      VAN_ALBADA = 5,
      QUICK = 6,
      UMIST = 7,
      SOU = 8,
      FROMM = 9,
      ULTIMATE_QUICKEST = 10,
      SUPER_C = 11,
      HYPER_C = 12
    };

    static void setAdvectionFunction(Element *element);


    static double fluxUpNeighbour(Element *element, double dt, double T[]);

    static double fluxUpJunction(Element *element, double dt, double T[]);

    static double fluxUpJunctionBC(Element *element, double dt, double T[]);

    static double fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double T[]);

    static double fluxDownNeighbourUpstreamJunction(Element *element, double dt, double T[]);

    static double fluxDownNeighbourUpstreamJunctionBC(Element *element, double dt, double T[]);


    static double fluxDownNeighbour(Element *element, double dt, double T[]);

    static double fluxDownJunction(Element *element, double dt, double T[]);

    static double fluxDownJunctionBC(Element *element, double dt, double T[]);

    static double fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double T[]);

    static double fluxUpNeighbourDownstreamJunction(Element *element, double dt, double T[]);

    static double fluxUpNeighbourDownstreamJunctionBC(Element *element, double dt, double T[]);


    static double fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpJunction(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpJunctionBC(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbourUpstreamJunction(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbourUpstreamJunctionBC(Element *element, double dt, double S[], int soluteIndex);


    static double fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownJunction(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownJunctionBC(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpNeighbourDownstreamJunction(Element *element,  double dt, double S[], int soluteIndex);

    static double fluxUpNeighbourDownstreamJunctionBC(Element *element,  double dt, double S[], int soluteIndex);


    static double computeTVDLimiter(double r, TVDFluxLimiter limiter, Element *element, int upstream = 0);

};

#endif // ELEMENTADVTVD_H
