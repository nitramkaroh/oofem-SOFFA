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

#ifndef triangularsurfacequadratic_h
#define triangularsurfacequadratic_h

#include "sm/Elements/structural3dsurfaceelement.h"
#include "floatmatrixf.h"


#define _IFT_TriangularSurfaceQuadratic_Name "triangularsurfacequadratic"

namespace oofem {
// class FEI3dQuadLin;
class FEI2dTrQuad;

class TriangularSurfaceQuadratic : public Structural3DSurfaceElement
{
protected:
    static FEI2dTrQuad interpolation;

    /// Local vertex coordinates
    std ::vector<FloatArray> lc;
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;


public:
    TriangularSurfaceQuadratic( int n, Domain *d );
    virtual ~TriangularSurfaceQuadratic()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TriangularSurfaceQuadratic_Name; }
    const char *giveClassName() const override { return "TriangularSurfaceQuadratic"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 18; }
    int computeNumberOfGlobalDofs() override { return 18; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

    int computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords ) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }


}; // end namespace oofem
} // namespace oofem
#endif