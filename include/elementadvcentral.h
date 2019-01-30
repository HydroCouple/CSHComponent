#ifndef ELEMENTADVCENTRAL_H
#define ELEMENTADVCENTRAL_H

struct Element;

class ElementAdvCentral
{
  public:

    static void setAdvectionFunction(Element *element);

    static double fluxUpNeighbour(Element *element, double dt, double T[]);

    static double fluxDownNeighbour(Element *element, double dt, double T[]);

    static double fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex);
};


#endif // ELEMENTADVCENTRAL_H
