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

#ifndef bifurcationinterface_h
#define bifurcationinterface_h

//#include "nummet.h"
//#include "convergedreason.h"
//#include "linsystsolvertype.h"
//#include "sparsemtrxtype.h"
#include "sparsemtrx.h"
#include "oofemcfg.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include <cstdio>
//#include "inputrecord.C"


namespace oofem {
//class EngngModel;
//class SparseMtrx;
//class FloatArray;

/**
 * Interface for bifurcation and continuation analysis
 */
class OOFEM_EXPORT BifurcationInterface 
{
protected:
    bool bifurcation          = false;
    double alpha              = 10.;
    bool choleskyBif          = false;
    bool deflationBifurcation = false;
    bool foundLimitPoint      = false;
    bool postBifurcationLineSearchSolver = false;
    bool LineSearchState                 = false;
    int numberOfFoundUnstableSolutions  = 0;
    double p                                = 1.; // deflation norm exponent

    // for deflation bifurcation
    FloatArray x0_Defl, dx_Defl;

    // Saved eigenvectors and eigenvalues
    FloatArray eigenvalues; // stored eigenvalues
    FloatMatrix evectors; // stored eigenvectors
    FloatArray Xeigs; // computed at X

    // Definiteness of the tangent stiffness at the last converged solution.
    // Filled in by checkPD, reported in the .out file by the solver printState.
    bool pdCheckPerformed = false; // true if checkPD has been run for the current step
    bool pdLastResult     = true; // result of the last checkPD call
    int pdNumNegPivots    = 0; // number of negative pivots of the LDL^T factorization (inertia)
    double pdMinPivot     = 0.; // smallest diagonal entry of D; a pivot, not an eigenvalue
    int pdMatrixSize      = 0; // size of the checked matrix, identifies a condensed one

public:
    /**
     * Constructor.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    BifurcationInterface(){};
    virtual ~BifurcationInterface(){};
    //void initializeFrom( InputRecord &ir ){};
    
    virtual bool checkPD( SparseMtrx &A ) = 0;
    virtual void CholeskyUpdate( SparseMtrx &Ae, FloatArray &b, FloatArray &x ) { OOFEM_LOG_INFO( "Solver cannot perform Cholesky bifurcation\n" ); }; 
    virtual int computeEigenValuesVectors( SparseMtrx &A, FloatArray &evaluesFA, FloatMatrix &evectorsFM )
    {
        OOFEM_LOG_INFO( "Solver cannot compute eigenvectors\n" );
        return 0;
    };
    virtual int computeEigenValues( SparseMtrx &A, FloatArray &evaluesFA )
    {
        OOFEM_LOG_INFO( "Solver cannot compute eigenvalues\n" );
        return 0;
    };
    virtual bool canCholeskyBifurcation() { return false; };
    virtual bool canComputeEigenVectors() { return false; };

    void performBifurcation( SparseMtrx &Ae, FloatArray &b, FloatArray &x );
    void setEigenValuesVectors( SparseMtrx &A, FloatArray &Xeigs );
    void setBifurcation( bool doBif ) { this->bifurcation = doBif; }
    void setCholesky( bool doChol ) { this->choleskyBif = doChol; }
    void setDeflation( bool doDeflation ) { this->deflationBifurcation = doDeflation; }
    void setX0Defl( FloatArray &x0 ) { this->x0_Defl = x0; }
    void setFoundLimitPoint( bool setLP ) { this->foundLimitPoint = setLP; }
    bool giveFoundLimitPoint() const { return this->foundLimitPoint; }
    void setPostBifurcationLineSearchSolver( bool setPBLS ) { this->postBifurcationLineSearchSolver = setPBLS; }
    bool givePostBifurcationLineSearchSolver() const { return this->postBifurcationLineSearchSolver; }
    void compute_dx_defl( FloatArray &xcurr ){this->dx_Defl = xcurr - this->x0_Defl;}
    void setLineSearchState( bool lsstate ) { this->LineSearchState = lsstate; }
    bool giveLineSearchState() const { return this->LineSearchState; }
    FloatArray &giveX0Defl() { return this->x0_Defl; }
    bool getBifurcation() const { return this->bifurcation; }
    void setAlpha( double alphaNew ) { this->alpha = alphaNew; }
    double giveAlpha() const { return this->alpha; }
    double give_p() const { return this->p; }
    void storeEigenValuesVectors( FloatArray &evaluesFA, FloatMatrix &evectorsFM, FloatArray &Xeigs );
    //{
    //    this->eigenvalues = evaluesFA;
    //    this->evectors    = evectorsFM;
    //    this->Xeigs       = Xeigs;}
    FloatMatrix &getEigenVectors() { return this->evectors; };
    FloatArray &getEigenValues() { return this->eigenvalues; };
    FloatArray &getXeigs() { return this->Xeigs; };

    void incrementNumFoundSols() { this->numberOfFoundUnstableSolutions++; };
    void nullNumFoundSols() { this->numberOfFoundUnstableSolutions = 0; };
    int giveNumFoundSols() const { return this->numberOfFoundUnstableSolutions; }

    /**
     * Records the definiteness of the checked matrix, so that it can be reported
     * after the step has been completed. Called by checkPD implementations.
     * @param isPD True if the checked matrix is positive definite.
     * @param numNegPivots Number of negative pivots of the LDL^T factorization. By
     * Sylvester's law of inertia this equals the number of negative eigenvalues.
     * @param minPivot Smallest diagonal entry of D. Its sign is meaningful, its
     * magnitude is a pivot and must not be interpreted as the smallest eigenvalue.
     * @param matrixSize Size of the checked matrix. Reported so that a check performed
     * on a statically condensed matrix can be told apart from one on the full tangent.
     */
    void storePDStatus( bool isPD, int numNegPivots, double minPivot, int matrixSize )
    {
        this->pdCheckPerformed = true;
        this->pdLastResult     = isPD;
        this->pdNumNegPivots   = numNegPivots;
        this->pdMinPivot       = minPivot;
        this->pdMatrixSize     = matrixSize;
    }
    /**
     * Discards the recorded definiteness. Has to be called at the beginning of every
     * solve, so that the verdict of the previous step (or of a previous attempt of the
     * same step, when the time step is reduced) is never reported for the current one.
     */
    void resetPDStatus()
    {
        this->pdCheckPerformed = false;
        this->pdLastResult     = true;
        this->pdNumNegPivots   = 0;
        this->pdMinPivot       = 0.;
        this->pdMatrixSize     = 0;
    }
    /// Returns true if a definiteness check has been recorded for the current step.
    bool givePDCheckPerformed() const { return this->pdCheckPerformed; }
    /// Returns the result of the last recorded definiteness check.
    bool giveLastPDResult() const { return this->pdLastResult; }
    /// Returns the number of negative pivots (equal to the number of negative eigenvalues).
    int giveNumNegPivots() const { return this->pdNumNegPivots; }
    /// Returns the smallest diagonal entry of D of the last recorded check.
    double giveMinPivot() const { return this->pdMinPivot; }
    /// Returns the size of the matrix of the last recorded check.
    int givePDMatrixSize() const { return this->pdMatrixSize; }
    /**
     * Prints the recorded definiteness of the converged solution to the given stream.
     * Prints nothing when no check has been recorded for the current step, so that a
     * step which did not converge is not decorated with a stale verdict.
     */
    void printPDStatus( FILE *outputStream ) const;
};
} // end namespace oofem
#endif // bifurcationinterface_h
