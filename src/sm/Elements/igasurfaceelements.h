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

#ifndef igasurfaceelements_h
#define igasurfaceelements_h

#include "iga/iga.h"
#include "iga/feibspline.h"
#include "iga/feinurbs.h"
#include "iga/feitspline.h"
//#include "sm/Elements/PlaneStress/planestresselementevaluator.h"

#include "sm/Elements/SurfaceElements/space2dsurfaceelementevaluator.h"
//#include "sm/CrossSections/simplesurfacecrosssection.h"

#include "sm/Elements/3D/space3delementevaluator.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "mathfem.h"

//#define _IFT_BsplinePlaneStressElement_Name "bsplineplanestresselement"
//#define _IFT_NURBSPlaneStressElement_Name "nurbsplanestresselement"
//#define _IFT_TSplinePlaneStressElement_Name "tsplineplanestresselement"
//#define _IFT_NURBSSpace3dElement_Name "nurbs3delement"
#define _IFT_BsplinePlaneStrainSurfaceElement_Name "bsplineplanestrainsurfaceelement"
#define _IFT_BsplineAxisymSurfaceElement_Name "bsplineaxisymsurfaceelement"
namespace oofem {

class BsplinePlaneStrainSurfaceElement : public IGAElement, public Space2dStructuralSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplinePlaneStrainSurfaceElement( int n, Domain *aDomain );

    void initializeFrom( InputRecord &ir ) override;
    int checkConsistency() override;

    void giveCharacteristicMatrix( FloatMatrix &answer, CharType mtrx, TimeStep *tStep ) override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveCharacteristicMatrix( answer, mtrx, tStep );
    }
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t ) override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveCharacteristicVector( answer, type, mode, t );
    }

    FEInterpolation *giveInterpolation() const override { return const_cast<BSpline2dLineInterpolation *>( &this->interpolation ); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveDofManDofIDMask( inode, answer );
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState( TimeStep *tStep ) override { Space2dStructuralSurfaceElementEvaluator ::updateInternalState( tStep ); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_BsplinePlaneStrainSurfaceElement_Name; }
    const char *giveClassName() const override { return "BsplinePlaneStrainSurfaceElement"; }

    MaterialMode giveMaterialMode() override { return _3dMat; }

    FloatArray giveNormal( GaussPoint *gp ) const override;

    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;

protected:
    int giveNsd() override { return 1; }
};


//// Axisym

class BsplineAxisymSurfaceElement : public IGAElement, public AxisymSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplineAxisymSurfaceElement( int n, Domain *aDomain );

    void initializeFrom( InputRecord &ir ) override;
    int checkConsistency() override;

    void giveCharacteristicMatrix( FloatMatrix &answer, CharType mtrx, TimeStep *tStep ) override
    {
        AxisymSurfaceElementEvaluator ::giveCharacteristicMatrix( answer, mtrx, tStep );
    }
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t ) override
    {
        AxisymSurfaceElementEvaluator ::giveCharacteristicVector( answer, type, mode, t );
    }

    FEInterpolation *giveInterpolation() const override { return const_cast<BSpline2dLineInterpolation *>( &this->interpolation ); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override
    {
        AxisymSurfaceElementEvaluator ::giveDofManDofIDMask( inode, answer );
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState( TimeStep *tStep ) override { AxisymSurfaceElementEvaluator ::updateInternalState( tStep ); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_BsplineAxisymSurfaceElement_Name; }
    const char *giveClassName() const override { return "BsplineAxisymSurfaceElement_Name"; }

    MaterialMode giveMaterialMode() override { return _3dMat; }

    FloatArray giveNormal( GaussPoint *gp ) const override;

    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;

protected:
    int giveNsd() override { return 1; }
};

//class SecondGradientInterface 
//{
//}


} // end namespace oofem
#endif // igaelements_h
