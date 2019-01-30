#ifndef ELEMENTADVUPWIND_H
#define ELEMENTADVUPWIND_H

struct Element;

class ElementAdvUpwind
{
  public:

    static void setAdvectionFunction(Element *element);

    static double inFluxUpJunction(Element *element, double dt, double T[]);

    static double inFluxUpNeighbour(Element *element, double dt, double T[]);

    static double inFluxSelf(Element *element, double dt, double T[]);

    static double outFluxSelf(Element *element, double dt, double T[]);

    static double outFluxDownJunction(Element *element, double dt, double T[]);

    static double outFluxDownNeighbor(Element *element, double dt, double T[]);

    static double inFluxUpJunction(Element *element, double dt, double S[], int soluteIndex);

    static double inFluxUpNeighbour(Element *element, double dt, double S[], int soluteIndex);

    static double inFluxSelf(Element *element, double dt, double S[], int soluteIndex);

    static double outFluxSelf(Element *element, double dt, double S[], int soluteIndex);

    static double outFluxDownJunction(Element *element, double dt, double S[], int soluteIndex);

    static double outFluxDownNeighbor(Element *element, double dt, double S[], int soluteIndex);
};

#endif // ELEMENTADVUPWIND_H
