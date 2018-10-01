/*!
 *  \file    CSHComponentinfo.h
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
 *  \todo
 *  \warning
 */

#ifndef CSHCOMPONENTINFO_H
#define CSHCOMPONENTINFO_H

#include "cshcomponent_global.h"
#include "core/abstractmodelcomponentinfo.h"


class CSHComponent_EXPORT CSHComponentInfo : public AbstractModelComponentInfo
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "CSHComponentInfo")

  public:

    CSHComponentInfo(QObject *parent = nullptr);

    virtual ~CSHComponentInfo();

    HydroCouple::IModelComponent* createComponentInstance() override;
};


Q_DECLARE_METATYPE(CSHComponentInfo*)

#endif //CSHCOMPONENTINFO_H
