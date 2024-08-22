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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#pragma once

#include "activebc.h"
#include "element.h"
#include "inputrecord.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
///@name Input fields for HardMagneticBoundaryCondition
//@{
#define _IFT_HardMagneticBoundaryCondition_Name "hardmagneticboundarycondition"
#define _IFT_HardMagneticBoundaryCondition_mu_0 "mu_0"
#define _IFT_HardMagneticBoundaryCondition_b_ext "b_ext"
#define _IFT_HardMagneticBoundaryCondition_ltf "ltf" //load time function for applied field
//@}

namespace oofem {
/**
 * This class implements a hard magnetic boundary condition
 *
 * @author Ondrej Faltus
 *
 * @note
 */
    
class HardMagneticBoundaryCondition : public ActiveBoundaryCondition {

protected:

    double mu0; //vacuum permeability
    Tensor1_3d b_ext; //external magnetic field in free space
  //FloatArrayF<9> sigma_star; //precomputed free space stress
    Tensor2_3d maxwell_stress;
    int ltf_index; //index of load time function for applied load

public:
    /// <summary>
    ///  Base Constructor
    /// </summary>
    /// <param name="n"></param>
    /// <param name="d"></param>
    HardMagneticBoundaryCondition( int n, Domain *d ) :
        ActiveBoundaryCondition( n, d ) {}
    /// Destructor.
    virtual ~HardMagneticBoundaryCondition() {}


    virtual void initializeFrom(InputRecord &ir) override;

        /**
     * Assembles B.C. contributions to specified matrix.
     * @param[in,out] answer Matrix to assemble to.
     * @param tStep Active time step.
     * @param tStep Active time step.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     * @param scale Scaling factor.
     * @param omp_lock optional OMP lock to ensure correct update of answer
     */
    virtual void assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0, void *lock = nullptr ) override;

    /**
     * Assembles B.C. contributions to specified vector.
     * @param[in,out] answer Vector to assemble to.
     * @param tStep Active time step.
     * @param type Type of matrix to assemble.
     * @param mode Mode of value.
     * @param s Numbering scheme.
     * @param eNorms Norms for each dofid.
     * @param omp_lock optional OMP lock to ensure correct update of answer
     * @return Equivalent of the sum of the element norm (squared) of assembled vector.
     */
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms=nullptr, void*lock=nullptr) override;

    /**
     * Gives a list of location arrays that will be assembled.
     * This should only be used to construct zero structure in sparse matrices.
     * The rows and columns location arrays returned in tuples (stored in vector),
     * allowing to efficiently assemble and allocate off-diagonal blocks.
     * The nonzero entries are assembled and allocated for entries at (rows[i], cols[i]) positions.
     * @param rows List of location arrays for r_s.
     * @param cols List of location arrays for c_s.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     */
    virtual void giveLocationArrays( std ::vector<IntArray> &rows, std ::vector<IntArray> &cols, CharType type,
        const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s ) override;

    const char *giveInputRecordName() const override { return _IFT_HardMagneticBoundaryCondition_Name; }
    const char *giveClassName() const override { return "HardMagneticBoundaryCondition"; }

protected:
    /**
     * Helper function for computing the contributions to the load vector.
     */
    void computeLoadVectorFromElement(FloatArray &answer, Element *e, int side, TimeStep *tStep, ValueModeType mode);
    /**
     * Helper function for computing the tangent (@f$ K = \frac{\mathrm{d}F}{\mathrm{d}u} @f$)
     */
    void computeTangentFromElement(FloatMatrix &answer, Element *e, int side, TimeStep *tStep);

private:

    void evaluateFreeSpaceStress();
};
} // end namespace oofem
