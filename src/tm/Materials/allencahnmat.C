#include "tm/Materials/allencahnmat.h"
#include "classfactory.h"
#include "element.h"
#include "gausspoint.h"
#include "inputrecord.h"
#include "materialmode.h"
#include "node.h"
#include "timestep.h"
#include <cmath>

namespace oofem {

REGISTER_Material(AllenCahnMaterial);

AllenCahnMaterial::AllenCahnMaterial(int n, Domain *d) : TransportMaterial(n, d) { }

void
AllenCahnMaterial::initializeFrom(InputRecord &ir)
{
    Material::initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, conductivity, _IFT_AllenCahnMaterial_k);
    normalConductivity = conductivity;
    IR_GIVE_OPTIONAL_FIELD(ir, normalConductivity, _IFT_AllenCahnMaterial_normalConductivity);
    IR_GIVE_OPTIONAL_FIELD(ir, capacity, _IFT_AllenCahnMaterial_c);
    IR_GIVE_OPTIONAL_FIELD(ir, mobility, _IFT_AllenCahnMaterial_mobility);
    IR_GIVE_OPTIONAL_FIELD(ir, surfaceBeta, _IFT_AllenCahnMaterial_surfaceBeta);
    IR_GIVE_OPTIONAL_FIELD(ir, hyperbolic, _IFT_AllenCahnMaterial_hyperbolic);
    IR_GIVE_OPTIONAL_FIELD(ir, density, _IFT_AllenCahnMaterial_rho);
    IR_GIVE_OPTIONAL_FIELD(ir, surfaceCapacity, _IFT_AllenCahnMaterial_surfaceCapacity);
    IR_GIVE_OPTIONAL_FIELD(ir, manufactured, _IFT_AllenCahnMaterial_manufactured);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsBase, _IFT_AllenCahnMaterial_mmsBase);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsAmplitude, _IFT_AllenCahnMaterial_mmsAmplitude);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsRate, _IFT_AllenCahnMaterial_mmsRate);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsBeta, _IFT_AllenCahnMaterial_mmsBeta);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsLx, _IFT_AllenCahnMaterial_mmsLx);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsLy, _IFT_AllenCahnMaterial_mmsLy);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsTemporal, _IFT_AllenCahnMaterial_mmsTemporal);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsSpatial, _IFT_AllenCahnMaterial_mmsSpatial);
    if ( mmsLx <= 0.0 || mmsLy <= 0.0 ) {
        OOFEM_ERROR("MMS domain lengths must be positive");
    }
    if ( surfaceBeta <= 0.0 ) {
        OOFEM_ERROR("Allen-Cahn surfacebeta must be positive");
    }
    if ( density <= 0.0 || surfaceCapacity <= 0.0 ) {
        OOFEM_ERROR("Allen-Cahn rho and surfacecapacity must be positive");
    }
    if ( mmsTemporal != 0 && mmsTemporal != 1 ) {
        OOFEM_ERROR("Unknown manufactured Allen-Cahn temporal profile %d", mmsTemporal);
    }
    if ( mmsSpatial < 0 || mmsSpatial > 3 ) {
        OOFEM_ERROR("Unknown manufactured Allen-Cahn spatial profile %d", mmsSpatial);
    }
}

double
AllenCahnMaterial::computeManufacturedSource(GaussPoint *gp, TimeStep *tStep) const
{
    if ( !manufactured ) {
        return 0.0;
    }

    // Exact fields on [0,mmsLx] x [0,mmsLy]:
    // 0: mixed cosine/sine profile; 1: sin(pi*x/Lx)sin(pi*y/Ly);
    // 2: (x/Lx)(y/Ly); 3: x*y on a polygonal disk.
    // Its trace is smooth around the rectangle corners, and its normal derivative
    // is nonzero on the edges when beta != 0. This option is only for MMS checks.
    const auto &coords = gp->giveGlobalCoordinates();
    const double x = coords.at(1);
    const double y = coords.at(2);
    const double time = tStep->giveTargetTime();
    const double b = mmsTemporal == 0 ? mmsAmplitude * std::exp(-mmsRate * time)
                                      : mmsAmplitude * std::sin(mmsRate * time);
    const double beta = mmsBeta;
    const double xi = M_PI * x / mmsLx;
    const double eta = M_PI * y / mmsLy;
    const double c = mmsSpatial == 0
        ? std::cos(xi) * std::cos(eta) + beta * std::sin(xi) * std::sin(eta)
        : ( mmsSpatial == 1 ? std::sin(xi) * std::sin(eta)
            : ( mmsSpatial == 2 ? (x / mmsLx) * (y / mmsLy) : x * y ) );
    const double u = mmsBase + b * c;
    const double u_t = mmsTemporal == 0 ? -mmsRate * b * c
                                        : mmsAmplitude * mmsRate * std::cos(mmsRate * time) * c;
    const double u_tt = mmsTemporal == 0 ? mmsRate * mmsRate * b * c
                                         : -mmsAmplitude * mmsRate * mmsRate * std::sin(mmsRate * time) * c;

    if ( gp->giveElement()->giveGeometryType() == EGT_line_1 ) {
        if ( mmsSpatial == 3 ) {
            Element *edge = gp->giveElement();
            const FloatArray &x1 = edge->giveNode(1)->giveCoordinates();
            const FloatArray &x2 = edge->giveNode(2)->giveCoordinates();
            const double dx = x2.at(1) - x1.at(1);
            const double dy = x2.at(2) - x1.at(2);
            const double length = std::sqrt(dx * dx + dy * dy);
            const double tx = dx / length;
            const double ty = dy / length;
            double nx = ty;
            double ny = -tx;
            const double xm = 0.5 * (x1.at(1) + x2.at(1));
            const double ym = 0.5 * (x1.at(2) + x2.at(2));
            if ( nx * xm + ny * ym < 0.0 ) {
                nx = -nx;
                ny = -ny;
            }
            const double surfaceLaplacian = 2.0 * b * tx * ty;
            const double normalDerivative = b * (y * nx + x * ny);
            return capacity * u_t - conductivity * surfaceLaplacian
                - surfaceBeta * mobility * (u - u * u * u)
                + normalConductivity * normalDerivative;
        }
        // On each straight edge, the sine profiles have a tangential Laplacian;
        // the polynomial x*y profile is affine along each edge, hence Delta_Gamma u=0.
        const bool horizontalEdge = std::abs(y) < 1.e-8 || std::abs(y - mmsLy) < 1.e-8;
        const bool verticalEdge = std::abs(x) < 1.e-8 || std::abs(x - mmsLx) < 1.e-8;
        if ( !horizontalEdge && !verticalEdge ) {
            OOFEM_ERROR("Manufactured Allen-Cahn boundary source requires boundary elements on the rectangle defined by mmslx and mmsly");
        }
        const double tangentLength = horizontalEdge ? mmsLx : mmsLy;
        const double surfaceLaplacian = mmsSpatial == 2 ? 0.0
            : -M_PI * M_PI * b * c / ( tangentLength * tangentLength );
        const double normalLength = horizontalEdge ? mmsLy : mmsLx;
        const double normalCoordinate = horizontalEdge ? x / mmsLx : y / mmsLy;
        double normalDerivative;
        if ( mmsSpatial == 2 ) {
            const double tangentCoordinate = horizontalEdge ? x / mmsLx : y / mmsLy;
            const double outwardSign = normalCoordinate < 0.5 ? -1.0 : 1.0;
            normalDerivative = outwardSign * b * tangentCoordinate / normalLength;
        } else {
            const double normalAmplitude = mmsSpatial == 0 ? beta : 1.0;
            normalDerivative = -normalAmplitude * M_PI * b
                * std::sin(M_PI * normalCoordinate) / normalLength;
        }
        return capacity * u_t - conductivity * surfaceLaplacian - surfaceBeta * mobility * ( u - u * u * u )
            + normalConductivity * normalDerivative;
    }

    const double bulkLaplacian = mmsSpatial >= 2 ? 0.0
        : -M_PI * M_PI * b * c
            * ( 1.0 / ( mmsLx * mmsLx ) + 1.0 / ( mmsLy * mmsLy ) );
    const double inertial = hyperbolic ? density * u_tt : 0.0;
    return inertial + capacity * u_t - conductivity * bulkLaplacian - mobility * ( u - u * u * u );
}

FloatArrayF<3>
AllenCahnMaterial::computeFlux3D(const FloatArrayF<3> &grad, double field,
                                 GaussPoint *gp, TimeStep *) const
{
    auto status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    status->setTempField(field);
    status->setTempGradient(grad);
    auto flux = -conductivity * grad;
    status->setTempFlux(flux);
    return flux;
}

FloatMatrixF<3,3>
AllenCahnMaterial::computeTangent3D(MatResponseMode, GaussPoint *, TimeStep *) const
{
    return conductivity * eye<3>();
}

double
AllenCahnMaterial::giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp,
                                          TimeStep *) const
{
    if ( mode == Capacity ) {
        if ( hyperbolic ) {
            return gp->giveElement()->giveGeometryType() == EGT_line_1 ? surfaceCapacity : density;
        }
        return capacity;
    } else if ( mode == IntSource ) {
        // The source vector is placed on the right hand side. Its tangent is
        // therefore moved to the left hand side with a minus sign.
        const double u = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempField();
        const double beta = gp->giveElement()->giveGeometryType() == EGT_line_1 ? surfaceBeta : 1.0;
        return beta * mobility * ( 3.0 * u * u - 1.0 );
    }

    OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode));
    return 0.0;
}

void
AllenCahnMaterial::computeInternalSourceVector(FloatArray &val, GaussPoint *gp,
                                               TimeStep *tStep, ValueModeType) const
{
    const double u = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempField();
    const double beta = gp->giveElement()->giveGeometryType() == EGT_line_1 ? surfaceBeta : 1.0;
    val = FloatArray{ beta * mobility * ( u - u * u * u ) + this->computeManufacturedSource(gp, tStep) };
}

} // namespace oofem
