#include "tm/Materials/allencahnmat.h"
#include "classfactory.h"
#include "element.h"
#include "gausspoint.h"
#include "inputrecord.h"
#include "materialmode.h"
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
    IR_GIVE_OPTIONAL_FIELD(ir, manufactured, _IFT_AllenCahnMaterial_manufactured);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsBase, _IFT_AllenCahnMaterial_mmsBase);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsAmplitude, _IFT_AllenCahnMaterial_mmsAmplitude);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsRate, _IFT_AllenCahnMaterial_mmsRate);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsBeta, _IFT_AllenCahnMaterial_mmsBeta);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsLx, _IFT_AllenCahnMaterial_mmsLx);
    IR_GIVE_OPTIONAL_FIELD(ir, mmsLy, _IFT_AllenCahnMaterial_mmsLy);
    if ( mmsLx <= 0.0 || mmsLy <= 0.0 ) {
        OOFEM_ERROR("MMS domain lengths must be positive");
    }
}

double
AllenCahnMaterial::computeManufacturedSource(GaussPoint *gp, TimeStep *tStep) const
{
    if ( !manufactured ) {
        return 0.0;
    }

    // Exact field on [0,mmsLx] x [0,mmsLy]:
    // u*=mmsBase+b(t)[cos(pi*x/Lx)cos(pi*y/Ly)+beta sin(pi*x/Lx)sin(pi*y/Ly)].
    // Its trace is smooth around the rectangle corners, and its normal derivative
    // is nonzero on the edges when beta != 0. This option is only for MMS checks.
    const auto &coords = gp->giveGlobalCoordinates();
    const double x = coords.at(1);
    const double y = coords.at(2);
    const double time = tStep->giveTargetTime();
    const double b = mmsAmplitude * std::exp(-mmsRate * time);
    const double beta = mmsBeta;
    const double xi = M_PI * x / mmsLx;
    const double eta = M_PI * y / mmsLy;
    const double c = std::cos(xi) * std::cos(eta)
        + beta * std::sin(xi) * std::sin(eta);
    const double u = mmsBase + b * c;
    const double u_t = -mmsRate * b * c;

    if ( gp->giveElement()->giveGeometryType() == EGT_line_1 ) {
        // On each straight edge, the trace has Delta_Gamma u = -pi^2 b c.
        const bool horizontalEdge = std::abs(y) < 1.e-8 || std::abs(y - mmsLy) < 1.e-8;
        const bool verticalEdge = std::abs(x) < 1.e-8 || std::abs(x - mmsLx) < 1.e-8;
        if ( !horizontalEdge && !verticalEdge ) {
            OOFEM_ERROR("Manufactured Allen-Cahn boundary source requires boundary elements on the rectangle defined by mmslx and mmsly");
        }
        const double tangentLength = horizontalEdge ? mmsLx : mmsLy;
        const double surfaceLaplacian = -M_PI * M_PI * b * c / ( tangentLength * tangentLength );
        const double normalLength = horizontalEdge ? mmsLy : mmsLx;
        const double normalCoordinate = horizontalEdge ? x / mmsLx : y / mmsLy;
        const double normalDerivative = -beta * M_PI * b * std::sin(M_PI * normalCoordinate) / normalLength;
        return capacity * u_t - conductivity * surfaceLaplacian - mobility * ( u - u * u * u )
            + normalConductivity * normalDerivative;
    }

    const double bulkLaplacian = -M_PI * M_PI * b * c
        * ( 1.0 / ( mmsLx * mmsLx ) + 1.0 / ( mmsLy * mmsLy ) );
    return capacity * u_t - conductivity * bulkLaplacian - mobility * ( u - u * u * u );
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
        return capacity;
    } else if ( mode == IntSource ) {
        // The source vector is placed on the right hand side. Its tangent is
        // therefore moved to the left hand side with a minus sign.
        const double u = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempField();
        return mobility * ( 3.0 * u * u - 1.0 );
    }

    OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode));
    return 0.0;
}

void
AllenCahnMaterial::computeInternalSourceVector(FloatArray &val, GaussPoint *gp,
                                               TimeStep *tStep, ValueModeType) const
{
    const double u = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempField();
    val = FloatArray{ mobility * ( u - u * u * u ) + this->computeManufacturedSource(gp, tStep) };
}

} // namespace oofem
