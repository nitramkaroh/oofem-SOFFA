#ifndef arclengthmaterialinterface_h
#define arclengthmaterialinterface_h

#include "interface.h"
//#include "interfacetype.h"


namespace oofem {
class ArcLengthMaterialInterface: public Interface
{
protected:
    bool hasALP = false;
    bool returnAll = true; // parameter determining if all material parameters should be used or only the arclength one
    double lambdaCALM = 0.; // calm loading parameter

public:
    ArcLengthMaterialInterface() :
        Interface()
    {};
    ~ArcLengthMaterialInterface() {};

    void setLambda(double lambdaSet) {
        this->lambdaCALM = lambdaSet;
    }

    void setReturnAll( bool val ) {
        this->returnAll = val;
    }

    virtual bool hasArcLengthParameter() const
    {
        return this->hasALP;
    };

};
} // namespace oofem
#endif