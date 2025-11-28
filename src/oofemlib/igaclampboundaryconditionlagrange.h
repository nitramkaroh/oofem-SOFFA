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

#ifndef igaclampboundaryconditionlagrange_h
#define igaclampboundaryconditionlagrange_h

// #include "floatarray.h"
#include "activebc.h"
#include "iga/iga.h"
#include "boundarycondition.h"


///@name Input fields for _IFT_IGAClampBoundaryCondition
//@{
#define _IFT_IGAClampBoundaryConditionLagrange_Name "igaclampboundaryconditionlagrange"
#define _IFT_IGAClampBoundaryConditionLagrange_ElementNumber "numel"
#define _IFT_IGAClampBoundaryConditionLagrange_LocalPosition "locpos"
#define _IFT_IGAClampBoundaryConditionLagrange_BcKnotSpan "bcknotspan"


//@}

namespace oofem {
// class Domain;
// class SparseMtrx;
// class TimeStep;
// class DofManager;
// class GaussPoint;
// class UnknownNumberingScheme;
// class FloatMatrix;
// class IntegrationRule;
// class Node;

class OOFEM_EXPORT IGAClampBoundaryConditionLagrange : public ActiveBoundaryCondition
// class OOFEM_EXPORT IGAClampBoundaryCondition : public BoundaryCondition
{
protected:
    int elementNumber; // surface element where the zero derivative is prescribed
    FloatArray localPosition; // local position where the zero derivative is prescribed
    IntArray BcKnotSpan;
    std ::unique_ptr<DofManager> lmdm;

public:
    /// Constructor.
    IGAClampBoundaryConditionLagrange( int n, Domain *d );
    //: ActiveBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~IGAClampBoundaryConditionLagrange() {};

    void initializeFrom( InputRecord &ir ) override;

    void assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0,
        void *lock = nullptr ) override;

    void assembleVector( FloatArray &answer, TimeStep *tStep,
        CharType type, ValueModeType mode,
        const UnknownNumberingScheme &s, FloatArray *eNorms = NULL,
        void *lock = nullptr ) override;


    const char *giveClassName() const override { return "IGAClampBoundaryConditionLagrange"; }
    const char *giveInputRecordName() const override { return _IFT_IGAClampBoundaryConditionLagrange_Name; }

    int giveLocalCodeNumbers( IntArray &answer, Element *elem, IntArray &span );

     int giveNumberOfInternalDofManagers() override { return 1; }
     DofManager *giveInternalDofManager( int i ) override { return this->lmdm.get(); }


    // double computeTangentFromContact( FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep );
    // void computeGap( double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep );

    // void computeNormalMatrixAt( FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep );


    // void computeExternalForcesFromContact( FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep );

     void giveLocationArrays( std ::vector<IntArray> &rows, std ::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s ) override;
     void giveLagrangianMultiplierLocationArray( const UnknownNumberingScheme &r_s, IntArray &answer );
};
} // end namespace oofem
#endif // node2nodecontact_h
