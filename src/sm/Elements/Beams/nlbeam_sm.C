#include "../sm/Elements/Beams/nlbeam_sm.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"

#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"
#include "function.h"
#include "domain.h"


namespace oofem {
REGISTER_Element(NlBeam_SM);


NlBeam_SM :: NlBeam_SM(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    this->internalForces.resize(3);
    this->fb_init.resize(3);

}

void
NlBeam_SM :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
  StructuralElement :: updateYourself(tStep);
  this->fb_init = this->internalForces;
}


void
NlBeam_SM :: initForNewStep()
// Updates the receiver at end of step.
{
  Element :: initForNewStep();
  this->internalForces =   this->fb_init;
}


void
NlBeam_SM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}




void
NlBeam_SM :: initializeFrom(InputRecord &ir)
{
    // first call parent
    NLStructuralElement :: initializeFrom(ir);
    this->computeLength();
    // Numerical parameters
    // 1. number of segments for numerical integration along the beam, default value 100
    IR_GIVE_OPTIONAL_FIELD(ir, NIP, _IFT_NlBeam_SM_NIP);
    IR_GIVE_FIELD(ir, EA, _IFT_NlBeam_SM_EA);
    IR_GIVE_FIELD(ir, EI, _IFT_NlBeam_SM_EI);
    IR_GIVE_FIELD(ir, GAs, _IFT_NlBeam_SM_GAs);
    //
    IR_GIVE_OPTIONAL_FIELD(ir, this->num, _IFT_NlBeam_SM_num);
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM_Beam_MaxIteration);
    // initialize initial position vector
    FloatArray a = {this->giveNode(1)->giveCoordinate(1), this->giveNode(1)->giveCoordinate(3), 0, this->giveNode(2)->giveCoordinate(1), this->giveNode(2)->giveCoordinate(3), 0};
    X = FloatArrayF<6>(a);
    // set size of the fields
    this->x.resize(NIP+1);
    this->z.resize(NIP+1);
    this->phi.resize(NIP+1);
    //compute alpha
    this->alpha = 0;
    // calculate the loading -- change the name
    double fx = 0, fz = 0;
    ////
    IR_GIVE_OPTIONAL_FIELD(ir, fx, _IFT_NlBeam_SM_fx);
    IR_GIVE_OPTIONAL_FIELD(ir, fx_ltf, _IFT_NlBeam_SM_fxLTF);
    ////
    IR_GIVE_OPTIONAL_FIELD(ir, fz, _IFT_NlBeam_SM_fz);
    IR_GIVE_OPTIONAL_FIELD(ir, fz_ltf, _IFT_NlBeam_SM_fzLTF);
    ////
    M = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, M, _IFT_NlBeam_SM_m);
    IR_GIVE_OPTIONAL_FIELD(ir, m_ltf, _IFT_NlBeam_SM_mLTF);
    //
    Fx.resize(NIP+1);
    Fz.resize(NIP+1);
    this->computeLoadingVariables(fx, fz);

    /*    this->vN.resize(2*NIP+1);
    this->vQ.resize(NIP+1);
    this->vM.resize(NIP+1);
    */
    /*
      
    this->s.resize(NIP+1);
    this->x.resize(NIP+1);
    this->z.resize(NIP+1);
    this->phi.resize(NIP+1);
    this->kappa.resize(NIP+1);

    */
}


double
NlBeam_SM :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( beamLength == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        beamLength  = sqrt(dx * dx + dy * dy);
    }

    return beamLength;
}



  
void
NlBeam_SM :: computeLoadingVariables(double fx, double fz)
{
  // assuming only load with constant intesity
  double dxi = beamLength/NIP;
  Fx.at(1) = 0.5 * fx * dxi;
  Fz.at(1) = 0.5 * fz * dxi;
  
  for (int i=2; i <= NIP; i++) {
    Fx.at(i) = Fx.at(i-1) + fx * dxi;
    Fz.at(i) = Fz.at(i-1) + fz * dxi;
  }
  Fx.at(NIP+1) = Fx.at(NIP) + 0.5 * fx * dxi;
  Fz.at(NIP+1) = Fz.at(NIP) + 0.5 * fz * dxi;
}


double 
NlBeam_SM :: computeLoadAt(FloatArray &Px, FloatArray &Pz, TimeStep *tStep)
{
  double m;
  if(fx_ltf != 0) {
    Px = Fx * this->giveDomain()->giveFunction(fx_ltf)->evaluate(tStep, VM_Total);
  } else {
    Px = Fx;
  }
  
  if(fz_ltf != 0) {
    Pz = Fz * this->giveDomain()->giveFunction(fz_ltf)->evaluate(tStep, VM_Total);
  } else {
    Pz = Fz;
  }

  if(m_ltf != 0) {
    m = M * this->giveDomain()->giveFunction(m_ltf)->evaluate(tStep, VM_Total);
  } else {
    m = M;
  }

  return m;

}

  
  

  std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_SM :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
 {
   FloatArray Px, Pz;
   double m = this->computeLoadAt(Px, Pz, tStep);
   // initial configuration
   x.at(1) = xa.at(1);
   z.at(1) = xa.at(2);
   phi.at(1) = xa.at(3) - alpha;
   //
   double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
   //this->vM.at(1) = -Mab;
   //
   FloatArray dx(4), dz(4), dphi(3);
   dphi  = {0,0,0,1};
   //
   FloatArray dN(4), dQ(4), dM(4);
   dM = {0 , 0, -1, 0};
   //
   double dxi = beamLength/NIP;
   //
   double M = -Mab;
   double Mp = 0;
   FloatArray dMp;
   for (int i=2; i <= NIP+1; i++) {
     //this->s.at(i) = this->s.at(i-1) + ds;
     //
     double kappa = computeCurvatureFromMoment(M);
     double dMdKappa = computeDerMomentFromCurvature(kappa);
     auto dkappa = dM / dMdKappa;
     //
     double phi_mid = phi.at(i-1) + 0.5 * kappa * dxi;
     auto dphi_mid = dphi + 0.5 * dkappa * dxi;
     //
     double c = cos(phi_mid);
     double s = sin(phi_mid);
     double N_mid = - c * (Xab + Px.at(i)) + s * (Zab + Pz.at(i));
     double Q_mid = - s * (Xab + Px.at(i)) - c * (Zab + Pz.at(i));
     //
     /*
     vN.at(2*(i-1)) = N_mid;
     vQ.at(2*(i-1)) = Q_mid;
     */
     //
     auto dN_mid = - Q_mid * dphi_mid;
     dN_mid.at(1) += - c;
     dN_mid.at(2) +=   s;
     //
     auto dQ_mid = N_mid * dphi_mid;
     dQ_mid.at(1) += - s;
     dQ_mid.at(2) +=   c;
     //
     double lambda_mid = 1 + N_mid / EA;
     auto dLambda_mid = dN_mid / EA;
     //
     double gamma_mid = Q_mid / GAs;
     auto dGamma_mid = dQ_mid / GAs;
     //
     double delta_x = ( c * lambda_mid + s * gamma_mid ) * dxi;
     FloatArray ddx;
     ddx = ( c * dLambda_mid + s * dGamma_mid - ( s * lambda_mid - c * gamma_mid) *dphi_mid ) * dxi;
     //
     double delta_z = ( c * gamma_mid - s * lambda_mid ) * dxi;
     FloatArray ddz;
     ddz =( c * dGamma_mid - s * dLambda_mid - ( s * gamma_mid + c * lambda_mid) * dphi_mid ) * dxi;
     //
     x.at(i) = x.at(i-1) + delta_x;
     dx += ddx;
     //
     z.at(i) = z.at(i-1) + delta_z;
     dz += ddz;
     //
     Mp += - m * dxi + Px.at(i-1) * delta_z - Pz.at(i-1) * delta_x;
     dMp += Px.at(i-1) * ddz - Pz.at(i-1) * ddx;
     //
     M = -Mab + Xab *  (z.at(i) - z.at(1)) - Zab * (x.at(i) - x.at(1)) + Mp;
     //
     dM = dMp + FloatArray({z.at(i)-z.at(1) + Xab * dz.at(1) - Zab * dx.at(1), -x.at(i) + x.at(1) + Xab * dz.at(2) - Zab * dx.at(2), -1. + Xab * dz.at(3) - Zab * dx.at(3), Xab * dz.at(4) - Zab * dx.at(4)});
     //
     //     vM.at(i) = M;
     //
     //
     kappa = computeCurvatureFromMoment(M);
     dMdKappa = computeDerMomentFromCurvature(kappa);
     dkappa = dM / dMdKappa;
     phi.at(i) = phi_mid + 0.5 * kappa * dxi;
     dphi = dphi_mid + 0.5 * dxi * dkappa;
   }
   // right end displacements
   FloatArrayF<3> xb;
   xb.at(1) = x.at(NIP+1); xb.at(2) = z.at(NIP+1); xb.at(3) = phi.at(NIP+1);
   // jacobi matrix G
   FloatMatrix G(3,3);
   G.copySubVectorRow(FloatArray({dx.at(1), dx.at(2), dx.at(3)}), 1,1);
   G.copySubVectorRow(FloatArray({dz.at(1), dz.at(2), dz.at(3)}), 2,1);
   G.copySubVectorRow(FloatArray({dphi.at(1), dphi.at(2), dphi.at(3)}), 3,1);
   // 3rd column of Gr
   FloatArrayF<3> Gr;
   Gr.at(1) = dx.at(4);
   Gr.at(2) = dz.at(4);
   Gr.at(3) = dphi.at(4);
   return {xb, FloatMatrixF<3,3> (G), Gr, Mp, FloatArrayF<4> (dMp)};

 }

 /*
 Find forces and moment at the left end that lead to given displacements and rotations
 at both ends.
  The input displacements and output forces are taken with respect to the global coordinate system.
 Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
 */
std::tuple<bool,FloatArrayF<3>, double>
 NlBeam_SM :: findLeftEndForces(const FloatArrayF<6> &x, TimeStep *tStep)
 {
   FloatArray dforces;
   int iter = 0;
   double tolerance = beam_tol;
   FloatArray fa;
   fa = this->internalForces;
   FloatArray xa = x[{0,1,2}];
   FloatArray xb_target = x[{3,4,5}];
   
   if(xb_target.computeNorm() > 1.e-10) {
     tolerance *= xb_target.computeNorm();
   }
   auto [xb_shoot, G, Gr, Mp, dMp] = this->integrateAlongBeam(fa, xa, tStep);
   double Mpr = Mp;

   //compute the jacobi matrix numerically
   /*   FloatArray fnum, dMn(3);
   FloatMatrix Gn(3,3);
   double pert = 1.e-6;
   for(int i = 1; i <= 3; i++) {
     fnum = this->internalForces;
     fnum.at(i) += pert;
     auto [xb_pert, G, vGr, Mp_pert, dMp] = integrateAlongBeam(fnum, xa, tStep);
     for(int j = 1; j <= 3; j++ ) {
       Gn.at(j,i) = (xb_pert.at(j) - xb_shoot.at(j)) ;
     }
     dMn.at(i) = Mp_pert - Mp;
   }
   Gn.times(1./pert);
   */
   
   ////
   FloatMatrix jacobi(G);
   auto res = xb_target-xb_shoot; 
   double error = res.computeNorm();
   //
   while ((iter==0) || (error>tolerance && iter < beam_maxit && error ==error) ){
     iter++;
     jacobi.solveForRhs(res, dforces);
     fa += dforces;
     auto [xb_shoot, G, Gr, Mp, dMp] = this->integrateAlongBeam(fa, xa, tStep);
     Mpr = Mp;
     jacobi = G;
     res = xb_target - xb_shoot; 
     error = res.computeNorm();
   }
   
   if (iter >= beam_maxit || error != error) {
     // what to do ???
     return {false, FloatArrayF<3>(fb_init), Mpr};
   }
   return {true, FloatArrayF<3>(fa), Mpr};
   
 }

 /*
 Find end forces and moments that lead to given displacements and rotations
 at both ends.
 The input displacements and output forces are taken with respect to the global coordinate system.
 */
void
NlBeam_SM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
   // solution vector
   answer.resize(6);
   FloatArray u;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   auto y = X + u;
   // only the first three entries of f are computed
   auto [converged, fa, Mp] = this->findLeftEndForces(y, tStep);
   this->internalForces = fa;
   // internal forces
   answer.at(1) = this->internalForces.at(1);
   answer.at(2) = this->internalForces.at(2);
   answer.at(3) = this->internalForces.at(3);  
   //
   double c1 =  y.at(5) - y.at(2);
   double c2 =  y.at(4) - y.at(1);
   FloatArray Px, Pz;
   this->computeLoadAt(Px, Pz, tStep);
     //
   answer.at(4) = -answer.at(1) - Px.at(NIP + 1);
   answer.at(5) = -answer.at(2) - Pz.at(NIP + 1);
   answer.at(6) = c1 * answer.at(1) - c2 * answer.at(2) - answer.at(3) + Mp;
}


 /*
 Evaluate the tangent stiffness matrix based on given displacements at both ends
 and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
 */
 void
 NlBeam_SM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
 {
   answer.resize(6,6);
   // solution vector
   FloatArray u;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   auto y = X + u;
   FloatArray xa= FloatArray({y.at(1), y.at(2), y.at(3)});
   auto [xb_shoot, G, vGr, Mp, dMp] = this->integrateAlongBeam(this->internalForces, xa, tStep);   
   FloatMatrix Ginv;
   Ginv.beInverseOf(G);
   //
   FloatMatrix Gr = {{1,0,0}, {0,1,0}, {vGr.at(1), vGr.at(2), vGr.at(3)}};
   auto K11 = -1. * ( Ginv * Gr );   
   //ading sub matrix starting at 1 1
   answer.setSubMatrix(K11, 1, 1);
   //ading sub matrix starting at 1 4
   answer.setSubMatrix(Ginv, 1, 4);
   //ading sub vector starting at row 4, col 1
   FloatArray K1 = {answer.at(1,1),answer.at(1,2),answer.at(1,3),answer.at(1,4),answer.at(1,5),answer.at(1,6)};
   answer.addSubVectorRow(-1. * K1, 4, 1);
   //ading sub vector starting at row 5, col 1
   FloatArray K2 = {answer.at(2,1),answer.at(2,2),answer.at(2,3),answer.at(2,4),answer.at(2,5),answer.at(2,6)};
   answer.addSubVectorRow(-1. * K2, 5, 1);
   // construct the sixth row
   FloatArray K3 =  {answer.at(3,1),answer.at(3,2),answer.at(3,3),answer.at(3,4),answer.at(3,5),answer.at(3,6)};
   FloatArray K6 =   (y.at(5) - y.at(2) + dMp.at(1)) * K1 + (-y.at(4) + y.at(1) + dMp.at(2)) * K2 + (dMp.at(3)-1.) * K3;
   K6 += FloatArray({this->internalForces.at(2), -this->internalForces.at(1), 0, -this->internalForces.at(2), this->internalForces.at(1), 0});
   answer.addSubVectorRow(K6, 6, 1);
   //
   answer.at(6,3) += dMp.at(4);
   /*
   answer.at(6,1) =  answer.at(1,6);
   answer.at(6,2) =  answer.at(2,6);
   answer.at(6,3) =  answer.at(3,6);
   answer.at(6,4) =  answer.at(4,6);
   answer.at(6,5) =  answer.at(5,6);
   */

   if(num == true) {
     FloatMatrix jacobi(G);
     FloatMatrix K(6,6);
     this->computeStiffnessMatrix_num(K, rMode, tStep);
     FloatMatrix Ktest(answer);
     Ktest -= K;
    
     FloatMatrix aT, k23;
     k23 = K2-K3;
     aT.beTranspositionOf(answer);
     FloatMatrix aaT(answer);
     aaT -= aT;
     
     //compute the jacobi matrix numerically
     FloatArray fnum, dMn(3);
     FloatMatrix Gn(3,3);
     double pert = 1.e-6;
     for(int i = 1; i <= 3; i++) {
       fnum = this->internalForces;
       fnum.at(i) += pert;
       auto [xb_pert, G, vGr, Mp_pert, dMp] = integrateAlongBeam(fnum, xa, tStep);
       for(int j = 1; j <= 3; j++ ) {
	 Gn.at(j,i) = (xb_pert.at(j) - xb_shoot.at(j)) ;
       }
       dMn.at(i) = Mp_pert - Mp;
     }
     
     FloatArray Grn(3);
     auto xa_num = xa;
     xa_num.at(3) += pert;
     fnum = this->internalForces;
     auto [xb_pert2, G2, vGr2, Mp_pert2, dMp2] = integrateAlongBeam(this->internalForces, xa_num, tStep);
     for(int j = 1; j <= 3; j++ ) {
       Grn.at(j) = (xb_pert2.at(j) - xb_shoot.at(j)) ;
     }   
     double dMn4 = Mp_pert2 - Mp;
     
     Gn.times(1./pert);
     dMn.times(1./pert);
     Grn.times(1./pert);
     dMn4 /= pert;
     // Numerical stiffness
     auto [xb_shoot3, G3, vGr3, Mp3, dMp3] = this->integrateAlongBeam(this->internalForces, xa, tStep);   
     //   FloatMatrix K(6,6);
     //
     this->computeStiffnessMatrix_num(K, rMode, tStep);
     int test = 1;
   }
}
  


 void
 NlBeam_SM :: computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
 {
   // solution vector
   FloatArray u, iF, iFn, iP, iPn;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   this->giveInternalForcesVector_from_u(iF, tStep, u);

   double eps = 1.e-8;
   answer.resize(6,6);
    
   FloatArray du;
   /////////////
   for(int i = 1; i <= 6; i++) {
     du = u;
     du.at(i) += eps;
     this->giveInternalForcesVector_from_u(iFn, tStep, du);
     ////////////
     answer.at(1,i) = iFn.at(1) - iF.at(1);
     answer.at(2,i) = iFn.at(2) - iF.at(2);
     answer.at(3,i) = iFn.at(3) - iF.at(3);
     answer.at(4,i) = iFn.at(4) - iF.at(4);
     answer.at(5,i) = iFn.at(5) - iF.at(5);
     answer.at(6,i) = iFn.at(6) - iF.at(6);
   }
   answer.times(1./eps);
 }


void
NlBeam_SM :: giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u)
{
   // solution vector
   answer.resize(6);
   auto y = X + u;
   // only the first three entries of f are computed
   auto [converged, fa, Mp] = this->findLeftEndForces(y, tStep);
   // internal forces
   answer.at(1) = fa.at(1);
   answer.at(2) = fa.at(2);
   answer.at(3) = fa.at(3);  
   //
   double c1 =  y.at(5) - y.at(2);
   double c2 =  y.at(4) - y.at(1);
   FloatArray Px, Pz;
   this->computeLoadAt(Px, Pz, tStep);
     //
   answer.at(4) = -answer.at(1) - Px.at(NIP + 1);
   answer.at(5) = -answer.at(2) - Pz.at(NIP + 1);
   answer.at(6) = c1 * answer.at(1) - c2 * answer.at(2) - answer.at(3)+ Mp;
}



 

double
NlBeam_SM :: computeMomentFromCurvature(double kappa)
{
   return EI*kappa;
}

double
NlBeam_SM :: computeDerMomentFromCurvature(double kappa)
{
  return EI;
}

double
NlBeam_SM :: computeCurvatureFromMoment(double M)
{
  return M/EI;
}




Interface *NlBeam_SM :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}

void
NlBeam_SM :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
  
    vtkPieces.resize(1);
 
    int numCells = this->NIP;
    const int numCellNodes  = 2; // linear line
    int nNodes = numCells * numCellNodes;
    //
    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);
    //
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
    Node *nodeA = this->giveNode(1);
    int nodeNum = 1;
    double dxi = beamLength/NIP;

    FloatArray nodeCoords(3);
    IntArray connectivity(2);
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
      for (int iNode = 1; iNode <= numCellNodes; iNode++) {
	double L = (iElement-1) * dxi + (iNode -1) * dxi;
	nodeCoords.at(1) = nodeA->giveCoordinate(1) + L * cos(alpha);
	nodeCoords.at(3) = nodeA->giveCoordinate(3) + L * sin(alpha);
	nodeCoords.at(2) = 0;
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	nodeNum++;
	connectivity.at(iNode) = val++;
      }
      vtkPieces.at(0).setConnectivity(iElement, connectivity);
      offset += 2;
      vtkPieces.at(0).setOffset(iElement, offset);
      vtkPieces.at(0).setCellType(iElement, 3);
    }


    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  for ( int nN = 1; nN <= nNodes; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }
	    double L = (iNode -1) * dxi;
	    FloatArray X(3);
	    X.at(1) = nodeA->giveCoordinate(1) + L * cos(alpha);
	    X.at(3) = nodeA->giveCoordinate(3) + L * sin(alpha);
	    //	    
	    FloatArray u;
	    u = {this->x.at(iNode) - X.at(1),0,this->z.at(iNode) - X.at(3)};
	    vtkPieces.at(0).setPrimaryVarInNode(utype, nN, u);
	  }
        }
    }

    /////////////////
    /*    int nC = cellVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfCellVarsToExport(nC, numCells);
    for ( int i = 1; i <= nC; i++ ) {
      InternalStateType istype = ( InternalStateType ) cellVarsToExport.at(i);
      if ( istype == IST_BeamForceMomentTensor ) {
	for ( int iC = 1; iC <= numCells; iC++ ) {
	  FloatArray bfmt(6);
	  bfmt.at(1) = vN.at(iC);
	  bfmt.at(2) = vQ.at(iC);
	  bfmt.at(3) = vM.at(iC);
	  vtkPieces.at(0).setCellVar(istype, iC, bfmt);
	}
      }   
    }
    */
  
}

  /*
void
NlBeam_SM :: printOutputAt(FILE *file, TimeStep *tStep)
{


  FILE *FID;
  

  char fext[100];
  sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
  std :: string fileName, functionname, temp;
  fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
  size_t foundDot;
  foundDot = fileName.rfind(".");
  
  while (foundDot != std :: string :: npos) {
    fileName.replace(foundDot, 1, "_");
    foundDot = fileName.rfind(".");   
  }
  
  fileName += fext;
  temp = fileName;
  size_t backslash = temp.rfind("/");
  if (backslash != std :: string :: npos ) {
    functionname = temp.substr(backslash+1, std :: string :: npos);
  } else {
    functionname = temp;
  }


  
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }



  //transform displacements to global coordinate system
  FloatArray uab, ug(NIP+1), wg(NIP+1), phig(NIP+1), u_l, u_g;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ug.at(1) = uab.at(1);
  wg.at(1) = uab.at(2);
  phig.at(1) = uab.at(3);
  double L = 0;
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    L = L + dx;
    FloatArray l;
    FloatMatrix T;
    this->construct_T(T, uab.at(3));
    this->construct_l(l, uab.at(3), L);
    u_l = {this->u.at(i), this->w.at(i), 0};
    u_l.subtract(l);
    u_g.beTProductOf(T, u_l);
    ug.at(i) = u_g.at(1) + ug.at(1);
    wg.at(i) = u_g.at(2) + wg.at(1);
    phig.at(i) = this->phi.at(i) + phig.at(1);
    
  }

  
  

  fprintf( FID, " function [x z u w phi N M]= %s \n", functionname.c_str() );

   fprintf(FID, "x=[");
   for ( double val: s ) {
     fprintf( FID, "%f,", val * cosAlpha );
   }   
   fprintf(FID, "];\n");
   
   fprintf(FID, "z=[");
   for ( double val: s ) {
     fprintf( FID, "%f,", val * sinAlpha );
   }   
   fprintf(FID, "];\n");

   
   
   fprintf(FID, "u=[");
   for ( double val: ug ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "w=[");
   for ( double val: wg ) {
     fprintf( FID, "%f,", val );
   }     
   fprintf(FID, "];\n");
   
   fprintf(FID, "phi=[");
   for ( double val: phig ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "N=[");
   for ( double val: vN ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "M=[");
   for ( double val: vM ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");

   
   fprintf(FID, "end\n");

   fclose(FID);

}

  */

FILE *
NlBeam_SM :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
    std :: string fileName;
    fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
    fileName += fext;
    fileName += ".m";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    return answer;
}
  

} // end namespace oofem
