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

#ifndef lineplanestrainsurface_h
#define lineplanestrainsurface_h

#include "sm/Elements/structural2dsurfaceelement.h"
#include "floatmatrixf.h"
// #include "sm/ErrorEstimators/huertaerrorestimator.h"
// #include "zznodalrecoverymodel.h"
// #include "sprnodalrecoverymodel.h"
// #include "nodalaveragingrecoverymodel.h"
// #include "spatiallocalizer.h"

#define _IFT_LinePlaneStrainSurface_Name "lineplanestrainsurface"

namespace oofem {
// class FEI3dQuadLin;
class FEI1dLin;

class LinePlaneStrainSurface : public Structural2DSurfaceElement
{
protected:
    static FEI1dLin interpolation;

    /// Local vertex coordinates
    std ::vector<FloatArray> lc;
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;


public:
    LinePlaneStrainSurface( int n, Domain *d );
    virtual ~LinePlaneStrainSurface()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LinePlaneStrainSurface_Name; }
    const char *giveClassName() const override { return "LinePlaneStrainSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 4; }
    int computeNumberOfGlobalDofs() override { return 4; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;


    MaterialMode giveMaterialMode() override { return _3dMat; }
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    //double computeVolumeAround( GaussPoint *gp ) override;

    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    // void postInitialize() override;
protected:
    void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;


protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 2; }


}; // end namespace oofem
} // namespace oofem
#endif