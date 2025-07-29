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
    double alpha              = 10;
    bool choleskyBif          = false;
    bool deflationBifurcation = false;
    bool foundLimitPoint      = false;
    bool postBifurcationLineSearchSolver = false;
    bool LineSearchState                 = false;
    int numberOfFoundUnstableSolutions  = 0;

    // for deflation bifurcation
    FloatArray x0_Defl, dx_Defl;

    // Saved eigenvectors and eigenvalues
    FloatArray eigenvalues; // stored eigenvalues
    FloatMatrix evectors; // stored eigenvectors
    FloatArray Xeigs; // computed at X

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
    void storeEigenValuesVectors( FloatArray &evaluesFA, FloatMatrix &evectorsFM, FloatArray &Xeigs ) {
        this->eigenvalues = evaluesFA;
        this->evectors    = evectorsFM;
        this->Xeigs       = Xeigs;}
    FloatMatrix &getEigenVectors() { return this->evectors; };
    FloatArray &getEigenValues() { return this->eigenvalues; };
    FloatArray &getXeigs() { return this->Xeigs; };

    void incrementNumFoundSols() { this->numberOfFoundUnstableSolutions++; };
    void nullNumFoundSols() { this->numberOfFoundUnstableSolutions = 0; };
    int giveNumFoundSols() const { return this->numberOfFoundUnstableSolutions; }
};
} // end namespace oofem
#endif // bifurcationinterface_h
