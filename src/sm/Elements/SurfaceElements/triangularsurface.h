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

#ifndef triangularsurface_h
#define triangularsurface_h

#include "sm/Elements/structural3dsurfaceelement.h"
#include "floatmatrixf.h"


#define _IFT_TriangularSurface_Name "triangularsurface"

namespace oofem {
// class FEI3dQuadLin;
class FEI2dTrLin;

class TriangularSurface : public Structural3DSurfaceElement
{
protected:
    static FEI2dTrLin interpolation;

    /// Local vertex coordinates
    std ::vector<FloatArray> lc;
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;


public:
    TriangularSurface( int n, Domain *d );
    virtual ~TriangularSurface()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TriangularSurface_Name; }
    const char *giveClassName() const override { return "TriangularSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 9; }
    int computeNumberOfGlobalDofs() override { return 9; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 3; }


}; // end namespace oofem
} // namespace oofem
#endif