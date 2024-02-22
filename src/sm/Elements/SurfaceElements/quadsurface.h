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

#ifndef quadsurface_h
#define quadsurface_h

#include "sm/Elements/structural3dsurfaceelement.h"
#include "floatmatrixf.h"
//#include "sm/ErrorEstimators/huertaerrorestimator.h"
//#include "zznodalrecoverymodel.h"
//#include "sprnodalrecoverymodel.h"
//#include "nodalaveragingrecoverymodel.h"
//#include "spatiallocalizer.h"

#define _IFT_QuadSurface_Name "quadsurface"

namespace oofem {
//class FEI3dQuadLin;
class FEI2dQuadLin;

class QuadSurface : public Structural3DSurfaceElement
{
protected:
    static FEI2dQuadLin interpolation;

    /// Local vertex coordinates
    std ::vector<FloatArray> lc;
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;


public:
    QuadSurface( int n, Domain *d );
    virtual ~QuadSurface() {}
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QuadSurface_Name; }
    const char *giveClassName() const override { return "QuadSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;
    //void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer) override;
    //void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int li, int ui ) override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 12; }
    int computeNumberOfGlobalDofs() override { return 12; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

    //void postInitialize() override;


protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }
    


}; // end namespace oofem
} // namespace oofem
#endif 