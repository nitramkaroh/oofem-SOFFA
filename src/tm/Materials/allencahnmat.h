#ifndef allencahnmat_h
#define allencahnmat_h

#include "tm/Materials/transportmaterial.h"
#include "floatmatrixf.h"

#define _IFT_AllenCahnMaterial_Name "allencahn"
#define _IFT_AllenCahnMaterial_k "k"
#define _IFT_AllenCahnMaterial_c "c"
#define _IFT_AllenCahnMaterial_mobility "mobility"
#define _IFT_AllenCahnMaterial_normalConductivity "knormal"
#define _IFT_AllenCahnMaterial_manufactured "manufactured"
#define _IFT_AllenCahnMaterial_mmsBase "mmsbase"
#define _IFT_AllenCahnMaterial_mmsAmplitude "mmsamplitude"
#define _IFT_AllenCahnMaterial_mmsRate "mmsrate"
#define _IFT_AllenCahnMaterial_mmsBeta "mmsbeta"
#define _IFT_AllenCahnMaterial_mmsLx "mmslx"
#define _IFT_AllenCahnMaterial_mmsLy "mmsly"

namespace oofem {

/**
 * Isotropic Allen-Cahn material for the bulk equation
 *
 *     c u_t - k Delta u = mobility (u - u^3).
 *
 * The defaults c = k = mobility = 1 give the equation
 * u_t - Delta u + u^3 - u = 0.
 */
class AllenCahnMaterial : public TransportMaterial
{
protected:
    double conductivity = 1.0;
    // Coefficient multiplying the bulk normal flux in the manufactured
    // dynamic-boundary source. For a boundary material this is the bulk k.
    double normalConductivity = 1.0;
    double capacity = 1.0;
    double mobility = 1.0;
    bool manufactured = false;
    double mmsBase = 0.2;
    double mmsAmplitude = 0.1;
    double mmsRate = 1.0;
    double mmsBeta = 0.25;
    double mmsLx = 1.0;
    double mmsLy = 1.0;

    double computeManufacturedSource(GaussPoint *gp, TimeStep *tStep) const;

public:
    AllenCahnMaterial(int n, Domain *d);

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field,
                                 GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp,
                                       TimeStep *tStep) const override;
    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp,
                                   TimeStep *tStep) const override;
    bool hasInternalSource() const override { return true; }
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp,
                                     TimeStep *tStep, ValueModeType mode) const override;

    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_AllenCahnMaterial_Name; }
    const char *giveClassName() const override { return "AllenCahnMaterial"; }
};

} // namespace oofem

#endif // allencahnmat_h
