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
REGISTER_Element(NlBeam_Reissner);
REGISTER_Element(NlBeam_Ziegler);



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
    IR_GIVE_FIELD(ir, EA, "ea");
    IR_GIVE_FIELD(ir, EI, "ei");
    IR_GIVE_FIELD(ir, GAs, "gas");
    //
    IR_GIVE_OPTIONAL_FIELD(ir, this->num, _IFT_NlBeam_SM_num);
    IR_GIVE_OPTIONAL_FIELD(ir, this->sym, _IFT_NlBeam_SM_sym);
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
    this->alpha = this->computeAlpha();
    // Applied magnetic field
    double BaA = 0, mu0 = 1, theta = 0,rho = 0, Br = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, BaA, "baa");
    IR_GIVE_OPTIONAL_FIELD(ir, mu0, "mu0");
    IR_GIVE_OPTIONAL_FIELD(ir, theta, "theta");
    Btildax = (BaA/mu0)*cos(theta);
    Btildaz = (BaA/mu0)*sin(theta);
    // Residual magnetic field
    IR_GIVE_OPTIONAL_FIELD(ir, Br, "br");
    IR_GIVE_OPTIONAL_FIELD(ir, rho, "rho");
    auto Brx = Br*cos(rho);
    auto Brz = Br*sin(rho);
     IR_GIVE_OPTIONAL_FIELD(ir, Ba_ltf, "ba_ltf");
    //
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
    FBrx.resize(NIP+1);
    FBrz.resize(NIP+1);
    this->computeLoadingVariables(fx, fz, Brx, Brz);
    ///
    this->vN.resize(NIP+1);
    this->vQ.resize(NIP+1);
    this->vM.resize(NIP+1);
    ///    
    this->s.resize(NIP+1);
    this->x.resize(NIP+1);
    this->z.resize(NIP+1);
    //
    this->phi.resize(NIP+1);
    this->vKappa.resize(NIP+1);
    
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



double
NlBeam_SM :: computeAlpha()
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    xA     = nodeA->giveCoordinate(1);
    xB     = nodeB->giveCoordinate(1);
    yA     = nodeA->giveCoordinate(3);
    yB     = nodeB->giveCoordinate(3);
    auto    alpha  = atan2(yB - yA, xB - xA);
    return alpha;
}



  
  
void
NlBeam_SM :: computeLoadingVariables(double fx, double fz, double Brx, double Brz)
{
  // assuming only load with constant intesity
  double dxi = beamLength/NIP;
  Fx.at(1) = 0.5 * fx * dxi;
  Fz.at(1) = 0.5 * fz * dxi;

  FBrx.at(1) = 0.5 * Brx * dxi;
  FBrz.at(1) = 0.5 * Brz * dxi;
  
  for (int i=2; i <= NIP; i++) {
    Fx.at(i) = Fx.at(i-1) + fx * dxi;
    Fz.at(i) = Fz.at(i-1) + fz * dxi;

    FBrx.at(i) = Brx;//FBrx.at(i-1) + Brx * dxi;
    FBrz.at(i) = Brz;//FBrz.at(i-1) + Brz * dxi;
  }

  Fx.at(NIP+1) = Fx.at(NIP) + 0.5 * fx * dxi;
  Fz.at(NIP+1) = Fz.at(NIP) + 0.5 * fz * dxi;

  FBrx.at(NIP+1) = FBrx.at(NIP);// + 0.5 * Brx * dxi;
  FBrz.at(NIP+1) = FBrz.at(NIP);// + 0.5 * Brz * dxi;
}


double 
NlBeam_SM :: computeLoadAt(FloatArray &Px, FloatArray &Pz, double &tBax, double &tBaz, TimeStep *tStep)
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

  if(Ba_ltf != 0) {
    tBax = Btildax * this->giveDomain()->giveFunction(Ba_ltf)->evaluate(tStep, VM_Total);
    tBaz = Btildaz * this->giveDomain()->giveFunction(Ba_ltf)->evaluate(tStep, VM_Total);
  } else {
    tBax = Btildax;
    tBaz = Btildaz;
  }


  return m;

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
   
   if(num == 1) {
     //compute the jacobi matrix numerically
     FloatArray fnum, dMn(3);
     FloatMatrix Gn(3,3);
     double pert = 1.e-7;
     for(int i = 1; i <= 3; i++) {
       fnum = this->internalForces;
       fnum.at(i) += pert;
       auto [xb_pert, Gpert, vGr, Mp_pert, dMp] = integrateAlongBeam(fnum, xa, tStep);
       
       for(int j = 1; j <= 3; j++ ) {
	 Gn.at(j,i) = (xb_pert.at(j) - xb_shoot.at(j)) ;
       }
       dMn.at(i) = Mp_pert - Mp;
     }
     FloatMatrix tG(G);
     Gn.times(1./pert);
   }
     
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
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
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
   if(sym) {
     FloatMatrix at;
     at.beTranspositionOf(answer);
     answer += at;
     answer *= 0.5;
   }
}
  


 void
 NlBeam_SM :: computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
 {
   // solution vector
   FloatArray u, iF, iFn, iP, iPn;
   this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
   this->giveInternalForcesVector_from_u(iF, tStep, u);

   double eps = 1.e-6;
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
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
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

 
void
NlBeam_SM :: doMatlabOutput(std::string fileName)
{

  FILE *FID;
    fileName += "_element" + std::to_string(this->giveNumber());
  auto temp = fileName;
  std :: string functionname;
  size_t backslash = temp.rfind("/");
  if (backslash != std :: string :: npos ) {
    functionname = temp.substr(backslash+1, std :: string :: npos);
  } else {
    functionname = temp;
  }  

  //
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }
  //
  fprintf( FID, " function [x z phi kappa N Q M]= %s \n", functionname.c_str() );
  //
  fprintf(FID, "x=[");
  for ( double val: this->x ) {
    fprintf( FID, "%5.8f,", val);
  }   
  fprintf(FID, "];\n");
  //
  fprintf(FID, "z=[");
  for ( double val: this->z ) {
    fprintf( FID, "%5.8f,", val);
  }   
  fprintf(FID, "];\n");
  //      
  fprintf(FID, "phi=[");
  for ( double val: this->phi ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
    //      
  fprintf(FID, "kappa=[");
  for ( double val: this->vKappa ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  
  fprintf(FID, "N=[");
  for ( double val: this->vN ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "Q=[");
  for ( double val: this->vQ ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "M=[");
  for ( double val: this->vM ) {
    fprintf( FID, "%5.8f,", val );
  }
  fprintf(FID, "];\n");
  //
  fprintf(FID, "end\n");
  //
  fclose(FID);
}






  /////////////////////////Reissner//////////////////////////////////////////

  std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Reissner :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
 {
   FloatArray Px, Pz;
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
   // initial configuration
   this->x.at(1) = xa.at(1);
   this->z.at(1) = xa.at(2);
   this->phi.at(1) = xa.at(3) - alpha;
   //
   double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
   // store variables for postprocessing
   this->s.at(1) = 0.;
   this->vM.at(1) = -Mab;
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
   double Mm = 0;
   FloatArray dMm;
   for (int i=2; i <= NIP+1; i++) {
     this->s.at(i) = this->s.at(i-1) + dxi;
     //
     double kappa = computeCurvatureFromMoment(M);
     double dMdKappa = computeDerMomentFromCurvature(kappa);
     auto dkappa = dM / dMdKappa;
     //
     double phi_mid = this->phi.at(i-1) + 0.5 * kappa * dxi;
     auto dphi_mid = dphi + 0.5 * dkappa * dxi;
     //
     double c = cos(phi_mid);
     double s = sin(phi_mid);
     double N_mid = - c * (Xab + Px.at(i-1)) + s * (Zab + Pz.at(i-1)) + FBrx.at(i) * ( c * tBax - s * tBaz);
     double Q_mid = - s * (Xab + Px.at(i-1)) - c * (Zab + Pz.at(i-1)) + FBrx.at(i) * ( c * tBaz + s * tBax);
     // store variables for postprocessing
     this->vN.at(i-1) = N_mid;
     this->vQ.at(i-1) = Q_mid;
     //
     auto dN_mid = - Q_mid * dphi_mid;
     dN_mid.at(1) += - c;
     dN_mid.at(2) +=   s;
     //
     auto dQ_mid = N_mid * dphi_mid;
     dQ_mid.at(1) += - s;
     dQ_mid.at(2) += - c;
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
     this->x.at(i) = x.at(i-1) + delta_x;
     dx += ddx;
     //
     this->z.at(i) = z.at(i-1) + delta_z;
     dz += ddz;
     //
     Mp += - m * dxi + Px.at(i-1) * delta_z - Pz.at(i-1) * delta_x;
     dMp += Px.at(i-1) * ddz - Pz.at(i-1) * ddx;
     //
     Mm +=  FBrx.at(i) * ( tBaz * delta_x - tBax*delta_z ) + FBrz.at(i) * ( s * tBaz - c * tBax ) * dxi;
     dMm += FBrx.at(i) * ( tBaz * ddx - tBax * ddz ) + FBrz.at(i) * ( c * tBaz + s * tBax ) * dxi * dphi_mid;
     //
     M = -Mab + Xab *  (z.at(i) - z.at(1)) - Zab * (x.at(i) - x.at(1)) + Mp + Mm;
     //
     dM = dMp + dMm + FloatArray({z.at(i)-z.at(1) + Xab * dz.at(1) - Zab * dx.at(1), -x.at(i) + x.at(1) + Xab * dz.at(2) - Zab * dx.at(2), -1. + Xab * dz.at(3) - Zab * dx.at(3), Xab * dz.at(4) - Zab * dx.at(4)});
     //
     //
     kappa = computeCurvatureFromMoment(M);
     dMdKappa = computeDerMomentFromCurvature(kappa);
     dkappa = dM / dMdKappa;
     this->phi.at(i) = phi_mid + 0.5 * kappa * dxi;
     dphi = dphi_mid + 0.5 * dxi * dkappa;
     // store variables for postprocessing
     this->vKappa.at(i) = kappa;
     this->vM.at(i) = M;
     //
     
   }
   // right end displacements
   FloatArrayF<3> xb;
   xb.at(1) = x.at(NIP+1); xb.at(2) = z.at(NIP+1); xb.at(3) = phi.at(NIP+1) + alpha;
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
   return {xb, FloatMatrixF<3,3> (G), Gr, Mp+Mm, FloatArrayF<4> (dMp+dMm)};

 }





  /////////////////////////Ziegler//////////////////////////////////////////
  std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
NlBeam_Ziegler :: integrateAlongBeam(const FloatArray &fa, const FloatArray &xa, TimeStep *tStep)
 {
   FloatArray Px, Pz;
   double tBax, tBaz;
   double m = this->computeLoadAt(Px, Pz, tBax, tBaz,  tStep);
   // initial configuration
   this->x.at(1) = xa.at(1);
   this->z.at(1) = xa.at(2);
   this->phi.at(1) = xa.at(3) - alpha;
   //
   double Xab = fa.at(1), Zab = fa.at(2), Mab = fa.at(3);
   // store variables for postprocessing
   this->s.at(1) = 0.;
   this->vM.at(1) = -Mab;
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
   double Mm = 0;
   FloatArray dMm;
   //
   auto chi_mid = 0.;
   for (int i=2; i <= NIP+1; i++) {
     this->s.at(i) = this->s.at(i-1) + dxi;
     //
     double kappa = computeCurvatureFromMoment(M);
     double dMdKappa = computeDerMomentFromCurvature(kappa);
     auto dkappa = dM / dMdKappa;
     //
     double phi_mid = phi.at(i-1) + 0.5 * kappa * dxi;
     auto dphi_mid = dphi + 0.5 * dkappa * dxi;
     //
     auto iter = 0;
     double c, s, N_mid, Q_mid, lambda_mid;
     while (iter <= chi_iter) {
       iter++;
       c = cos(phi_mid - chi_mid);
       s = sin(phi_mid - chi_mid);
       N_mid = - c * (Xab + Px.at(i-1)) + s * (Zab + Pz.at(i-1));
       Q_mid = - s * (Xab + Px.at(i-1)) - c * (Zab + Pz.at(i-1));
       //
       lambda_mid = 1 + N_mid / EA;
       //
       auto F = GAs * chi_mid - lambda_mid * Q_mid;
       if(fabs(F) < chi_tol) {
	 break;
       }
       auto dF = GAs + lambda_mid * N_mid - Q_mid * Q_mid / EA;
       chi_mid -= F/dF;
       if(iter == chi_iter){
	 //@todo: cut the step...
	 OOFEM_ERROR("Maximum number of iteraton for calculation of chi reached");
       }
     }
     // store variables for postprocessing
     this->vN.at(i-1) = N_mid;
     this->vQ.at(i-1) = Q_mid;
     // linearization
     auto denom  = (GAs + lambda_mid * N_mid - Q_mid * Q_mid / EA);
     auto dchi_mid = (lambda_mid * N_mid - Q_mid * Q_mid / EA) * dphi_mid;
     dchi_mid.at(1) += (-c * Q_mid / EA - s * lambda_mid);
     dchi_mid.at(2) += (s * Q_mid / EA - c * lambda_mid);
     dchi_mid /= denom;
     //
     auto dN_mid = - Q_mid * (dphi_mid - dchi_mid);
     dN_mid.at(1) += - c;
     dN_mid.at(2) +=   s;
     // not needed ?
     /*     auto dQ_mid = N_mid * (dphi_mid - dchi_mid);
     dQ_mid.at(1) += - s;
     dQ_mid.at(2) += - c;
     */
     //
     auto dLambda_mid = dN_mid / EA;
     //
     double delta_x = c * lambda_mid  * dxi;
     FloatArray ddx;
     ddx = ( c * dLambda_mid - s * lambda_mid * ( dphi_mid - dchi_mid) ) * dxi;
     //
     double delta_z = - s * lambda_mid * dxi;
     FloatArray ddz;
     ddz = (- s * dLambda_mid - c * lambda_mid * ( dphi_mid - dchi_mid ) ) * dxi;
     //
     this->x.at(i) = this->x.at(i-1) + delta_x;
     dx += ddx;
     //
     this->z.at(i) = this->z.at(i-1) + delta_z;
     dz += ddz;
     //
     Mp += - m * dxi + Px.at(i-1) * delta_z - Pz.at(i-1) * delta_x;
     dMp += Px.at(i-1) * ddz - Pz.at(i-1) * ddx;
     //
     M = -Mab + Xab *  (z.at(i) - z.at(1)) - Zab * (x.at(i) - x.at(1)) + Mp;
     //
     dM = dMp + FloatArray({z.at(i)-z.at(1) + Xab * dz.at(1) - Zab * dx.at(1), -x.at(i) + x.at(1) + Xab * dz.at(2) - Zab * dx.at(2), -1. + Xab * dz.at(3) - Zab * dx.at(3), Xab * dz.at(4) - Zab * dx.at(4)});
     //
     //
     kappa = computeCurvatureFromMoment(M);
     dMdKappa = computeDerMomentFromCurvature(kappa);
     dkappa = dM / dMdKappa;
     this->phi.at(i) = phi_mid + 0.5 * kappa * dxi;
     dphi = dphi_mid + 0.5 * dxi * dkappa;
     // store variables for postprocessing
     vM.at(i) = M;
     vKappa.at(i) = kappa;
   }
   // right end displacements
   FloatArrayF<3> xb;
   xb.at(1) = x.at(NIP+1); xb.at(2) = z.at(NIP+1); xb.at(3) = phi.at(NIP+1) + alpha;
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
   return {xb, FloatMatrixF<3,3> (G), Gr, Mp+Mm, FloatArrayF<4> (dMp+dMm)};
 }

  

} // end namespace oofem
