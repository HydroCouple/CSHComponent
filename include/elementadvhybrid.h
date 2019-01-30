#ifndef ELEMENTADVHYBRID_H
#define ELEMENTADVHYBRID_H

struct Element;

class ElementAdvHybrid
{
  public:

    static void setAdvectionFunction(Element *element);

    static double fluxUpNeighbour(Element *element, double dt, double T[]);

    static double fluxDownNeighbour(Element *element, double dt, double T[]);

    static double fluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double fluxDownNeighbour(Element *element, double dt, double S[], int soluteIndex);
};


#endif // ELEMENTADVHYBRID_H
