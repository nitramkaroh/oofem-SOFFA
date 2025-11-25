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

#include "sm/Contact/ContactElement/trcontactelement.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(TrContactElement);

FEI2dTrLin TrContactElement :: interp(1, 2);

TrContactElement :: TrContactElement(int n, Domain *aDomain) :
    TrPlaneStrain(n, aDomain)
    // Constructor.
{
    numberOfDofMans  = 3;
    area = -1;
    numberOfGaussPoints = 2;
    um.resize(2,3);
}


FEInterpolation *TrContactElement :: giveInterpolation() const { return & interp; }


void
TrContactElement :: initializeFrom(InputRecord &ir)
{
    TrPlaneStrain::initializeFrom(ir);
    IR_GIVE_FIELD(ir, normalPenalty_, _IFT_TrContactElement_normalPenalty);
    IR_GIVE_FIELD( ir, cheatFactor_, _IFT_TrContactElement_cheatFactor );
    this->tangentialPenalty_ = 0;
    this->muStick_ = 0;           
    this->muSlide_ = 0;
    xl.resize(2,3);
    for (int i = 1; i <= 3; ++i) {
      const auto &node = this->giveCellGeometryWrapper()->giveVertexCoordinates(i);
      xl.at(1,i) = node.at(1);
      xl.at(2,i) = node.at(2);
    }

}

  

void
TrContactElement :: postInitialize()
{
    numberOfGaussPoints = 1;
    PlaneStrainElement :: postInitialize();
}




// -----------------------------------------------------------------------------
// Helper: shape functions and derivatives (shape109)
// -----------------------------------------------------------------------------

void
TrContactElement::shape(double xi, FloatMatrix &shape) const
{
    shape.resize(2, 2);
    // row 1: N1, N2
    shape.at(1, 1) = 0.5 * (1.0 - xi);
    shape.at(1, 2) = 0.5 * (1.0 + xi);
    // row 2: dN1/dxi, dN2/dxi
    shape.at(2, 1) = -0.5;
    shape.at(2, 2) =  0.5;
}

// -----------------------------------------------------------------------------
// Helper: build a and da matrices (pos109)
// -----------------------------------------------------------------------------

void
TrContactElement::pos(double xi, FloatMatrix &a, FloatMatrix &da) const
{
    FloatMatrix s;
    this->shape(xi, s);

    a.resize(2, 6);
    da.resize(2, 6);

    const double N1  = s.at(1, 1);
    const double N2  = s.at(1, 2);
    const double dN1 = s.at(2, 1);
    const double dN2 = s.at(2, 2);

    // a matrix (2x6)
    a.at(1, 1) = -N1;  a.at(1, 2) = 0.0;
    a.at(1, 3) = -N2;  a.at(1, 4) = 0.0;
    a.at(1, 5) =  1.0; a.at(1, 6) = 0.0;

    a.at(2, 1) = 0.0;  a.at(2, 2) = -N1;
    a.at(2, 3) = 0.0;  a.at(2, 4) = -N2;
    a.at(2, 5) = 0.0;  a.at(2, 6) =  1.0;

    // da matrix (2x6)
    da.at(1, 1) = -dN1; da.at(1, 2) = 0.0;
    da.at(1, 3) = -dN2; da.at(1, 4) = 0.0;
    da.at(1, 5) =  0.0; da.at(1, 6) = 0.0;

    da.at(2, 1) = 0.0;  da.at(2, 2) = -dN1;
    da.at(2, 3) = 0.0;  da.at(2, 4) = -dN2;
    da.at(2, 5) = 0.0;  da.at(2, 6) = 0.0;
}

// -----------------------------------------------------------------------------
// Helper: projection of slave node on master segment (project109)
// -----------------------------------------------------------------------------

double
TrContactElement::project(const FloatMatrix &xupd) const
{
    // xupd is 2x3: columns = nodes 1,2,3

    double x1 = xupd.at(1, 1);
    double y1 = xupd.at(2, 1);
    double x2 = xupd.at(1, 2);
    double y2 = xupd.at(2, 2);
    double x3 = xupd.at(1, 3);
    double y3 = xupd.at(2, 3);

    double dx = x2 - x1;
    double dy = y2 - y1;

    double num = 2.0 * (x3 * dx + y3 * dy)
                 - (x2 * x2 + y2 * y2)
                 + (x1 * x1 + y1 * y1);
    double den = dx * dx + dy * dy;

    return num / den;
}

// -----------------------------------------------------------------------------
// Helper: surface geometry (surfgeom109)
// -----------------------------------------------------------------------------

void
TrContactElement::surfGeom(const FloatMatrix &xupd,  // 2x3 updated positions
                           const FloatMatrix &da,    // 2x6 matrix from pos()
                           FloatArray &rhoxi,        // out: size 2
                           double &m_11,             // out
                           double &m11,              // out
                           double &detm,             // out
			   double &dA,
                           FloatArray &normal,       // out: size 2
                           FloatArray &tau) const    // out: size 2
{
    rhoxi.resize(2);
    normal.resize(2);
    tau.resize(2);

    rhoxi.zero();
    normal.zero();
    tau.zero();

    // 1) rhoxi = - da * xupd (with xupd seen as 6x1 via column-major order)
    //
    // MATLAB used linear indexing xupd(j) for a 2x3 matrix. In column-major:
    //   j=1 -> (1,1)
    //   j=2 -> (2,1)
    //   j=3 -> (1,2)
    //   j=4 -> (2,2)
    //   j=5 -> (1,3)
    //   j=6 -> (2,3)
    //
    for (int i = 1; i <= 2; ++i) {
        double sum = 0.0;
        for (int j = 1; j <= 6; ++j) {
            int row = ((j - 1) % 2) + 1;   // 1 or 2
            int col = ((j - 1) / 2) + 1;   // 1..3
            double xval = xupd.at(row, col);
            sum -= da.at(i, j) * xval;
        }
        rhoxi.at(i) = sum;
    }

    const double r1 = rhoxi.at(1);
    const double r2 = rhoxi.at(2);

    m_11 = r1 * r1 + r2 * r2;
    detm = m_11;

    if (m_11 <= 0.0) {
        m11 = 0.0;
        normal.zero();
        tau.zero();
        return;
    }

    m11 = 1.0 / m_11;
    const double sqrtDetm = std::sqrt(detm);

    // 3) normal vector (2D)
    //
    // MATLAB:
    // norm(1) = -(1/(2*sqrt(detm)))*(xupd(4) - xupd(2));
    // norm(2) = -(1/(2*sqrt(detm)))*(xupd(1) - xupd(3));
    //
    // xupd(1) -> (1,1)
    // xupd(2) -> (2,1)
    // xupd(3) -> (1,2)
    // xupd(4) -> (2,2)
    double x11 = xupd.at(1, 1);
    double y11 = xupd.at(2, 1);
    double x12 = xupd.at(1, 2);
    double y12 = xupd.at(2, 2);

    normal.at(1) = -(1.0 / (2.0 * sqrtDetm)) * (y12 - y11);
    normal.at(2) = -(1.0 / (2.0 * sqrtDetm)) * (x11 - x12);

    // 4) tangential direction tau = rhoxi / sqrt(detm)
    tau.at(1) = r1 / sqrtDetm;
    tau.at(2) = r2 / sqrtDetm;

    dA = sqrtDetm;
}

// -----------------------------------------------------------------------------
// Helper: normal penetration (gap) along normal direction
// g_n = n · (x_S - x_M(xi))
// -----------------------------------------------------------------------------

double
TrContactElement::penetr(const FloatMatrix &xupd,
                         const FloatArray &normal,
                         double xi) const
{
    // Slave node = node 3
    double xSx = xupd.at(1, 3);
    double xSy = xupd.at(2, 3);

    // Master nodes = 1 and 2
    double x1x = xupd.at(1, 1);
    double x1y = xupd.at(2, 1);
    double x2x = xupd.at(1, 2);
    double x2y = xupd.at(2, 2);

    // Shape functions at xi
    double N1 = 0.5 * (1.0 - xi);
    double N2 = 0.5 * (1.0 + xi);

    // Projection point on master
    double xMx = N1 * x1x + N2 * x2x;
    double xMy = N1 * x1y + N2 * x2y;

    // Gap vector
    double gx = xSx - xMx;
    double gy = xSy - xMy;

    // Normal gap
    double gn = normal.at(1) * gx + normal.at(2) * gy;

    return gn;
}

// -----------------------------------------------------------------------------
// Helper: tangential trial traction (tangenttr109)
// T_trial = -kT * m_11 * dxi, where kT = tangentialPenalty_
// -----------------------------------------------------------------------------

double
TrContactElement::tangentTractionTrial(double dxi, double m_11) const
{
    return -tangentialPenalty_ * m_11 * dxi;
}

// -----------------------------------------------------------------------------
// Core: compute contact state (geometry, penetration, friction state)
// Translated from main logic in stiff109.m (without assembling matrices)
// -----------------------------------------------------------------------------

void
TrContactElement::computeContactState(ContactState &cs,
                                      const FloatMatrix &xl,   // 2x3 reference coords
                                      const FloatMatrix &u,    // 2x3 displacements
                                      TimeStep * /*tStep*/,
                                      bool updateHistory)
{
    cs.active = false;
    cs.rhoxi.resize(2);
    cs.normal.resize(2);
    cs.tau.resize(2);
    cs.a.resize(2, 6);
    cs.da.resize(2, 6);

    cs.rhoxi.zero();
    cs.normal.zero();
    cs.tau.zero();
    cs.a.zero();
    cs.da.zero();

    // 1) updated coordinates xupd = xl + u
    FloatMatrix xupd = xl;
    xupd.add(u);

    // 2) projection
    cs.xi = this->project(xupd);

    // Projection outside the element -> no contact
    if (std::abs(cs.xi) > 1. + cheatFactor_) {
        cs.active = false;
        return;
    }

    // 3) position matrices a, da
    this->pos(cs.xi, cs.a, cs.da);

    // 4) surface geometry (rhoxi, m_11, m11, detm, normal, tau)
    this->surfGeom(xupd, cs.da,
                   cs.rhoxi, cs.m_11, cs.m11,
                   cs.detm, cs.dA, cs.normal, cs.tau);

    // 5) penetration
    cs.penetr = this->penetr(xupd, cs.normal, cs.xi);

    // Positive gap -> no contact
    if (cs.penetr > 0.0) {
        cs.active = false;
        return;
    }

    cs.active = true;

   
    if (this->isXiInitialized()) {
        cs.dxi = 0.0;
        if (updateHistory) {
            this->setXiPrev(cs.xi);
            this->setXiInitialized(true);
        }
    } else {
        double xiPrev = this->giveXiPrev();
        cs.dxi = cs.xi - xiPrev;
        if (updateHistory) {
            this->setXiPrev(cs.xi);
        }
    }

    // 7) trial tangential traction
    cs.T_trial = this->tangentTractionTrial(cs.dxi, cs.m_11);

    // 8) yield function for friction:
    //     phi = |T_trial|*sqrt(m11) - muStick_*normalPenalty_*|penetr|
    cs.phi = std::sqrt(cs.T_trial * cs.T_trial * cs.m11)
             - muStick_ * normalPenalty_ * std::abs(cs.penetr);

    if (cs.phi <= 0.0) {
        // sticking
        cs.T_real = cs.T_trial;
    } else {
        // sliding
        double signT = (cs.T_trial >= 0.0) ? 1.0 : -1.0;
        cs.T_real = signT * muSlide_ * normalPenalty_
                    * std::abs(cs.penetr) * std::sqrt(cs.m_11);
    }
}

// -----------------------------------------------------------------------------
// Internal forces vector: p = - kN*penetr * a^T*n - T_real*sqrt(m11) * a^T*tau
// -----------------------------------------------------------------------------

void
TrContactElement::giveInternalForcesVector(FloatArray &answer,
                                           TimeStep *tStep,
                                           int useUpdatedGpRecord)
{
    answer.resize(this->giveNumberOfDofs());
    answer.zero();

    // 1) gather geometry
    FloatMatrix u(2, 3);    // displacements
    this->computeDisplacementMatrix(u, tStep);
    this->setDisplacementMatrix(u);
    // 2) compute contact state (may update history)
    ContactState cs;
    bool updateHistory = (useUpdatedGpRecord != 0);
    this->computeContactState(cs, xl, u, tStep, updateHistory);

    if (!cs.active) {
        return; // no contribution
    }

    // 3) assemble internal forces
    FloatArray aTn, aTt;
    aTn.beTProductOf(cs.a, cs.normal); // a^T * n  -> 6x1
    aTt.beTProductOf(cs.a, cs.tau);    // a^T * tau -> 6x1

    double factorN = normalPenalty_ * cs.penetr;
    double factorT = cs.T_real * std::sqrt(cs.m11);
    double dA = 1; //cs.dA;

    answer.add(factorN * dA, aTn);
    answer.add(factorT * dA, aTt);
}

// -----------------------------------------------------------------------------
// Tangent stiffness matrix: s = sN + sT
// This follows the structure of stiff109.m (normal + tangential parts).
// -----------------------------------------------------------------------------

void
TrContactElement::computeStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode rMode,
                                         TimeStep *tStep)
{
    answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
    answer.zero();

    if (rMode != TangentStiffness) {
        // Implement other modes as needed.
        return;
    }

    // 1) geometry
    FloatMatrix u(2, 3);    // displacements
    this->giveDisplacementMatrix(u);
    // 2) contact state WITHOUT updating history
    ContactState cs;
    this->computeContactState(cs, xl, u, tStep, /*updateHistory=*/false);

    if (!cs.active) {
        return;
    }

    const double kN = normalPenalty_;
    const double kT = tangentialPenalty_;

    FloatMatrix sNm, sNr, sN;
    FloatMatrix sT1, sT2, sT;

    sNm.resize(6,6); sNm.zero();
    sNr.resize(6,6); sNr.zero();
    sN.resize(6,6);  sN.zero();

    sT1.resize(6,6); sT1.zero();
    sT2.resize(6,6); sT2.zero();
    sT.resize(6,6);  sT.zero();

    const FloatArray &n  = cs.normal;
    const FloatArray &t  = cs.tau;
    const FloatArray &rh = cs.rhoxi;
    const double dA = 1.0;//cs.dA;

    // --- Normal part ---
    // sNm = kN * a^T * n * n^T * a
    FloatMatrix nn(2,2);
    nn.beDyadicProductOf(n, n);   // n * n^T

    FloatMatrix tmp6x2;
    tmp6x2.beTProductOf(cs.a, nn);   // a^T * nn  (6x2)

    sNm.beProductOf(tmp6x2, cs.a);   // (6x2)*(2x6) = 6x6
    sNm.times(kN);

    // sNr = kN * penetr * m11 * da^T * n * rh^T * a
    FloatMatrix nr(2,2);
    nr.beDyadicProductOf(n, rh);     // n * rh^T

    FloatMatrix tmp6x2b;
    tmp6x2b.beTProductOf(cs.da, nr);  // da^T * nr (6x2)

    sNr.beProductOf(tmp6x2b, cs.a);   // (6x2)*(2x6)
    sNr.times(kN * cs.penetr * cs.m11);

    sN.add(sNm);
    sN.add(sNr);
    FloatMatrix sNrT;
    sNrT.beTranspositionOf(sNr);
    sN.add(sNrT);

    // --- Tangential part ---
    if (cs.phi <= 0.0) {
        // sticking
        // sT1 = kT * a^T * tau * tau^T * a
        FloatMatrix tt(2,2);
        tt.beDyadicProductOf(t, t); // tau * tau^T

        tmp6x2.beTProductOf(cs.a, tt);  // a^T * tt
        sT1.beProductOf(tmp6x2, cs.a);  // (6x2)*(2x6)
        sT1.times(kT);

        // sT2 = T_real * m11 * a^T * tau * tau^T * da
        tmp6x2.beTProductOf(cs.a, tt);  // a^T * tt
        sT2.beProductOf(tmp6x2, cs.da); // (6x2)*(2x6)
        sT2.times(cs.T_real * cs.m11);

        sT.add(-1.0, sT1);
        sT.add( 1.0, sT2);
        FloatMatrix sT2T;
        sT2T.beTranspositionOf(sT2);
        sT.add(1.0, sT2T);
    } else {
        // sliding
        const double signT = (cs.T_trial >= 0.0) ? 1.0 : -1.0;

        // sT1 = muSlide_*kN*sign(T_trial) * a^T * tau * norm^T * a
        FloatMatrix tn(2,2);
        tn.beDyadicProductOf(t, n);     // tau * n^T

        tmp6x2.beTProductOf(cs.a, tn);  // a^T * tn
        sT1.beProductOf(tmp6x2, cs.a);  // (6x2)*(2x6)
        sT1.times(muSlide_ * kN * signT);

        // sT2 = sign(T_trial)*muSlide_*kN*|penetr|*sqrt(m11) * da^T * tau * tau^T * a
        FloatMatrix tt(2,2);
        tt.beDyadicProductOf(t, t);      // tau * tau^T

        tmp6x2.beTProductOf(cs.da, tt);  // da^T * tt
        sT2.beProductOf(tmp6x2, cs.a);   // (6x2)*(2x6)
        sT2.times(signT * muSlide_ * kN * std::abs(cs.penetr) * std::sqrt(cs.m11));

        sT.add(-1.0, sT1);
        sT.add( 1.0, sT2);
        FloatMatrix sT2T;
        sT2T.beTranspositionOf(sT2);
        sT.add(1.0, sT2T);
    }

    // full tangent
    answer.add(dA, sN);
    answer.add(dA, sT);
}


void TrContactElement::computeDisplacementMatrix(FloatMatrix & u, TimeStep * tStep) const
{
  FloatArray un;
  for (int i = 1; i <= 3; ++i) {
      this->giveNode(i)->giveUnknownVector(un, { D_u, D_v}, VM_Total, tStep);
      u.at(1,i) = un.at(1);
      u.at(2,i) = un.at(2);
    }

}




} // end namespace oofem
