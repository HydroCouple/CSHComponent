#ifndef ELEMENTADVTVD_H
#define ELEMENTADVTVD_H

struct Element;

class ElementAdvTVD
{
  public:

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


    static double computeTVDLimiter(double r, int limiter);

};

#endif // ELEMENTADVTVD_H
