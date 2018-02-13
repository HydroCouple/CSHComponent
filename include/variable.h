#ifndef VARIABLE_H
#define VARIABLE_H

/*!
 * \brief The Variable struct represents a variable value
 */
struct Variable
{

    /*!
     * \brief value - Magnitude of the variable
     */
    double value = 0;

    /*!
     * \brief associatedValue - A related value that may be stored with this variable.
     */
    double associatedValue = 0;

    /*!
     * \brief isBC - Boolean indicating whether this variable is a boundary condition or not.
     */
    bool isBC = false;


    void copy(const Variable & variable)
    {
      value = variable.value;
      associatedValue = variable.associatedValue;
      isBC = variable.isBC;
    }

};



#endif // VARIABLE_H
