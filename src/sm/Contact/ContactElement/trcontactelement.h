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

#ifndef trcontactelement_h
#define trcontactelement_h

#include "sm/Elements/PlaneStrain/trplanestrain.h"

#define _IFT_TrContactElement_Name "trcontactelement"
#define _IFT_TrContactElement_normalPenalty "normalpenalty"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements an triangular three-node  plane-
 * strain elasticity finite element. Each node has 2 degrees of freedom.
 */
class TrContactElement : public TrPlaneStrain
{
protected:
    static FEI2dTrLin interp;
    double area;
    FloatMatrix um;

public:
    TrContactElement(int n, Domain * d);
    virtual ~TrContactElement() { }

    FEInterpolation *giveInterpolation() const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TrContactElement_Name; }
    const char *giveClassName() const override { return "TrContactElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_1;}

    void postInitialize() override;

    // OOFEM-style interface methods (signatures as you specified)
    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep,
                                  int useUpdatedGpRecord) override;

    void computeStiffnessMatrix(FloatMatrix &answer,
                                MatResponseMode rMode,
                                TimeStep *tStep) override;


  
protected:
    double giveXiPrev() const { return xiPrev; }
    void   setXiPrev(double v) { xiPrev = v; }

    bool   isXiInitialized() const { return xiInitialized; }
    void   setXiInitialized(bool b) { xiInitialized = b; }

    void   giveDisplacementMatrix(FloatMatrix &u) { u = this->um; }
    void   setDisplacementMatrix(FloatMatrix &un) { this->um = un; }

  
    double xiPrev       = 0.0;
    bool   xiInitialized = false;

private:
      // ---- material parameters (were d(1..4) in MATLAB) ----
    double normalPenalty_;     // d(1)
    double tangentialPenalty_; // d(2)
    double muStick_;           // d(3)
    double muSlide_;           // d(4)
    FloatMatrix xl; // initial geometry of the element
    // Small helper struct carrying all contact state for one Gauss point
    struct ContactState {
        bool active = false;
        double xi = 0.0;
        double dxi = 0.0;
        double penetr = 0.0;
        double m_11 = 0.0;
        double m11 = 0.0;
        double detm = 0.0;
        double T_trial = 0.0;
        double T_real = 0.0;
        double phi = 0.0;   // friction function
        double dA = 0.0;
        FloatArray rhoxi;   // size 2
        FloatArray normal;  // size 2
        FloatArray tau;     // size 2
        FloatMatrix a;      // 2x6
        FloatMatrix da;     // 2x6
    };

    // ---- core driver for contact state ----
    void computeContactState(ContactState &cs,
                             const FloatMatrix &xl,   // 2x3 reference coords
                             const FloatMatrix &u,    // 2x3 displacements
                             TimeStep *tStep,
                             bool updateHistory);

    // ---- geometry & constitutive helpers ----
    double project(const FloatMatrix &xupd) const;
    void   shape(double xi, FloatMatrix &shape) const;
    void   pos(double xi, FloatMatrix &a, FloatMatrix &da) const;
    void   surfGeom(const FloatMatrix &xupd, const FloatMatrix &da,
                    FloatArray &rhoxi, double &m_11, double &m11,
                    double &detm, double &dA, FloatArray &normal, FloatArray &tau) const;
    double penetr(const FloatMatrix &xupd, const FloatArray &normal,
                  double xi) const;
    double tangentTractionTrial(double dxi,
                                double m_11) const;
    void computeDisplacementMatrix(FloatMatrix &u, TimeStep *tStep) const;

    void initializeFrom(InputRecord &ir) override;
    void updateYourself(TimeStep *tStep) override {;}
    void printOutputAt(FILE *file, TimeStep *tStep) override{;}
    

  

};
} // end namespace oofem
#endif // trplanstrain_h
