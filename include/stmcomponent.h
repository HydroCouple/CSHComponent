/*!
 *  \file    stmcomponent.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 *  This file and its associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2018
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef STMCOMPONENT_H
#define STMCOMPONENT_H

#include "stmcomponent_global.h"
#include "stmcomponentinfo.h"
#include "core/abstractmodelcomponent.h"

class STMCOMPONENT_EXPORT STMComponent : public AbstractModelComponent
{
    Q_OBJECT

  public:

    /*!
     * \brief STMComponent constructor
     * \param id Unique identifier for this component instance.
     * \param modelComponentInfo the parent ModelComponentInfo that generated this component instance.
     */
    STMComponent(const QString &id, STMComponentInfo* modelComponentInfo = nullptr);

    /*!
     * \brief ~STMComponent destructor
     */
    virtual ~STMComponent();

    /*!
     * \brief validate validates this component model instance
     * \return Returns a list of error messages.
     */
    QList<QString> validate() override;

    /*!
     * \brief prepare Prepares the model component instance.
     */
    void prepare() override;

    /*!
     * \brief update
     * \param requiredOutputs
     */
    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    /*!
     * \brief finish
     */
    void finish() override;

  protected:

    /*!
     * \brief intializeFailureCleanUp
     */
    void intializeFailureCleanUp() override;

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    /*!
     * \brief updateOutputValues
     * \param requiredOutputs
     */
    void updateOutputValues(const QList<HydroCouple::IOutput*>& requiredOutputs) override;


  private:




};

#endif //STMCOMPONENT_H
