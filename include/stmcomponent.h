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
#include "temporal/abstracttimemodelcomponent.h"

#include <unordered_map>

class STMModel;
class HCGeometry;
struct Element;
class ElementInput;
class ElementOutput;
class Dimension;
class ElementHeatSourceInput;
class Unit;

class STMCOMPONENT_EXPORT STMComponent : public AbstractTimeModelComponent,
    public virtual HydroCouple::ICloneableModelComponent
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

    /*!
     * \brief modelInstance
     * \return
     */
    STMModel *modelInstance() const;

    /*!
     * \brief parent
     * \return
     */
    HydroCouple::ICloneableModelComponent* parent() const override;

    /*!
     * \brief clone
     * \return
     */
    HydroCouple::ICloneableModelComponent* clone() override;

    /*!
     * \brief clones
     * \return
     */
    QList<HydroCouple::ICloneableModelComponent*> clones() const override;

  protected:

    bool removeClone(STMComponent *component);


    /*!
     * \brief intializeFailureCleanUp
     */
    void initializeFailureCleanUp() override;

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    /*!
     * \brief createInputFileArguments
     */
    void createInputFileArguments();

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief initializeInputFilesArguments
     * \param message
     * \return
     */
    bool initializeInputFilesArguments(QString &message);

    /*!
     * \brief createGeometriesMap
     */
    void createGeometries();

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    /*!
     * \brief createFlowInput
     */
    void createFlowInput();

    /*!
     * \brief createXSectionAreaInput
     */
    void createXSectionAreaInput();

    /*!
     * \brief createDepthInput
     */
    void createDepthInput();

    /*!
     * \brief createTopWidthInput
     */
    void createTopWidthInput();

    /*!
     * \brief createRadiationFluxInput
     */
    void createExternalRadiationFluxInput();

    /*!
     * \brief createExternalHeatFluxInput
     */
    void createExternalHeatFluxInput();

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    /*!
     * \brief createTemperatureOutput
     */
    void createTemperatureOutput();

  private:

    IdBasedArgumentString *m_inputFilesArgument;

    ElementInput *m_flowInput,
                 *m_xSectionAreaInput,
                 *m_depthInput,
                 *m_topWidthInput;


    Unit *m_radiationFluxUnit,
         *m_heatFluxUnit,
         *m_temperatureUnit;

    ElementHeatSourceInput *m_externalRadiationFluxInput,
                           *m_externalHeatFluxInput;

    ElementOutput *m_temperatureOutput;

    Dimension *m_timeDimension,
              *m_geometryDimension;

    std::vector<QSharedPointer<HCGeometry>> m_elementGeometries;
    std::vector<QSharedPointer<HCGeometry>> m_elementJunctionGeometries;
    STMModel *m_modelInstance;

    STMComponent *m_parent;
    QList<HydroCouple::ICloneableModelComponent*> m_clones;
};

#endif //STMCOMPONENT_H
