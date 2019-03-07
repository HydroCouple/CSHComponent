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

    static double fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double T[]);

    static double fluxDownNeighbourUpstreamJunction(Element *element, double dt, double T[]);


    static double fluxDownNeighbour(Element *element, double dt, double T[]);

    static double fluxDownJunction(Element *element, double dt, double T[]);

    static double fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double T[]);

    static double fluxUpNeighbourDownstreamJunction(Element *element, double dt, double T[]);


    static double fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpJunction(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbourUpstreamNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbourUpstreamJunction(Element *element, double dt, double S[], int soluteIndex);


    static double fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownJunction(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpNeighbourDownstreamNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxUpNeighbourDownstreamJunction(Element *element,  double dt, double S[], int soluteIndex);


    static double computeTVDLimiter(double r, TVDFluxLimiter limiter, Element *element, int upstream = 0);

};

#endif // ELEMENTADVTVD_H
