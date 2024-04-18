/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef nlstructuralsurfaceelement_h
#define nlstructuralsurfaceelement_h

#include "sm/Elements/nlstructuralelement.h"

///@name Input fields for NLStructuralElement
//@{surfaceEvaldNdx
#define _IFT_NLStructuralSurfaceElement_nlgeoflag "nlgeo"
//@}

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Abstract base class for "structural" surace finite elements with geometrical nonlinearities.
 */
class NLStructuralSurfaceElement : public NLStructuralElement
{
protected:

public:

    NLStructuralSurfaceElement( int n, Domain *d );

    virtual ~NLStructuralSurfaceElement() {}

    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;

    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    virtual FloatArray giveNormal( GaussPoint *gp ) const = 0;

    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override
    {
        return false;
    };
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override
    {
        return 0;
    };

    /*void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;*/

};
} // end namespace oofem
#endif // nlstructuralelement_h
