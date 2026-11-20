#ifndef arclengthmaterialinterface_h
#define arclengthmaterialinterface_h

#include "interface.h"
//#include "interfacetype.h"


namespace oofem {
/**
 * Interface of a material driven by the material arc-length method (solver matcalm of StaticStructural): one
 * material parameter p is scaled by the load level lambda of the solver, p = lambda p1, where p1 is the input value.
 * The solver sets lambda in every material with such a parameter (setLambda). The material load vector, i.e. the
 * derivative of the internal forces with respect to lambda, is assembled with setReturnAll( false ): the material
 * then evaluates only the terms of p, with p = p1, and all other terms as zero (the internal forces are linear in p);
 * setReturnAll( true ) restores the full response with p = lambda p1.
 */
class ArcLengthMaterialInterface: public Interface
{
protected:
    bool hasALP = false; ///< The material has a parameter driven by the arc-length method.
    bool returnAll = true; ///< true: the full response (p = lambda p1); false: only the terms of p, with p = p1 (d/dlambda).
    double lambdaCALM = 0.; ///< Load level lambda of the arc-length method.

public:
    ArcLengthMaterialInterface() :
        Interface()
    {};
    ~ArcLengthMaterialInterface() {};

    /// Sets the load level lambda of the arc-length method.
    void setLambda(double lambdaSet) {
        this->lambdaCALM = lambdaSet;
    }

    /// true: evaluate the full response; false: only the part proportional to the arc-length parameter (see above).
    void setReturnAll( bool val ) {
        this->returnAll = val;
    }

    /// Whether the material has a parameter driven by the arc-length method (set from its input).
    virtual bool hasArcLengthParameter() const
    {
        return this->hasALP;
    };

};
} // namespace oofem
#endif