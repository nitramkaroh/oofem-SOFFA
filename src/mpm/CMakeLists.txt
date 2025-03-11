set (mpm_element
    termlibrary.C
    termlibrary2.C
    termlibrary3.C
    prototype.C
    prototype2.C
    up.C
    tm.C
    )

set (magnetoelasticity
    termlibrary4.C
    magnetoelasticelement.C
    MagnetoElasticity/Materials/hardmagneticmooneyrivlincompressiblematerial.C
    MagnetoElasticity/Materials/magnetoelasticmaterialstatus.C
    MagnetoElasticity/Materials/hardmagneticmooneyrivlincompressiblematerialprestrain.C
    MagnetoElasticity/CrossSections/magnetoelasticcrosssection.C
    )


    set (mpm_problems
    mpmproblem.C
    )

set (mpm
    ${mpm_element}
    ${mpm_problems}
    ${magnetoelasticity}
    )

add_library (mpm
    ${LIB_TYPE}
    ${mpm}    
    )






add_dependencies(mpm version)
