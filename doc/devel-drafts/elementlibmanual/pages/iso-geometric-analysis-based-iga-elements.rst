Iso Geometric Analysis based (IGA) elements
===========================================

The following record describes the common part of IGA element record:

\recentry**IGAElement**:elemparam:`num{in}`  \recentry\ :elemparam:`mat{in}`
:elemparam:`crossSect{in}` :elemparam:`nodes{ia}`    \recentry\
:elemparam:`knotvectoru{ra}` :elemparam:`knotvectorv{ra}` :elemparam:`knotvectorw{ra}`
\recentry\ :optelemparam:`knotmultiplicityu{ia}` :optelemparam:`knotmultiplicityv{ia}`
\recentry\ :optelemparam:`knotmultiplicityw{ia}`  \recentry\ :elemparam:`degree{ia}`
:elemparam:`nip{ia}`  \recentry\ :optelemparam:`nlgeo{in}` :optelemparam:`numexp{in}`
\recentry⟨:optelemparam:`partitions{ia}`\ ⟩ ⟨:optelemparam:`remote`\ ⟩

One IGA element is a whole patch. The :param:`knotvectoru`, :param:`knotvectorv`, and
:param:`knotvectorw` parameters specify knot vectors in individual parametric directions,
considering only distinct knots. Open knot vector is always assumed, so the multiplicity
of the first and last knot should be equal to :math:`p+1`, where :math:`p` is polynomial
degree in coresponding direction (determined by :param:`degree` parameter, see further).
The knot multiplicity can be set using optional parameters :param:`knotmultiplicityu`,
:param:`knotmultiplicityv`, and :param:`knotmultiplicityw`. By default, the open knot
vector is assumed and multiplicity of internal knots is assumed to be equal to one.
Note, that total number of knots in particular direction (including multiplicity) must
be equal to number of control points in this direction increased by degree in this
direction plus 1.  The degree of approximation for each parametric direction is
determined from :param:`degree` array, dimension of which is equal to number of
parametric dimensions of the patch (one for the curves of the surface elements below).
In case of elements with BSpline or Nurbs interpolation, the nodes forming the
rectangular array of control points of the element are ordered in a such way, that
u-index is changing most quickly, and w-index (or v-index in case of 2d problems) most
slowly. In case of elements with T-spline interpolation, the nodes forming the T-mesh of
the element are ordered arbitrarily.

The parameter :param:`nip` gives the number of integration points in each nonzero knot
span (a square number for two parametric dimensions, e.g. 9 for 3 x 3 points). The
parameter :param:`nlgeo` selects the small strain formulation (0, default) or the large
strain formulation in terms of the deformation gradient and the first Piola-Kirchhoff
stress (1); other values are rejected. The parameter :param:`numexp` (at least 2,
default 2) is the number of points per knot span and parametric direction in which the
elements supporting the VTK export (see below) are sampled.

The supported **IGAElement** values are following:

**Keyword**: :param:`bsplineplanestresselement`  **Parameters**: None.

**Keyword**: :param:`nurbsplanestresselement`  **Parameters**: None.

**Keyword**: :param:`nurbs3delement`  **Parameters**: None.

**Keyword**: :param:`tsplineplanestresselement`  **Parameters**:
:elemparam:`localindexknotvectoru{ia}` :elemparam:`localindexknotvectorv{ia}`
:elemparam:`localindexknotvectorw{ia}` :elemparam:`weights{ra}`  The parameters
:param:`localindexknotvectoru`, :param:`localindexknotvectorv`, and
:param:`localindexknotvectorw` defined by the (one based) indices to global knot
vectors (given by :param:`knotvectoru`, :param:`knotvectorv`, and :param:`knotvectorw`
parameters, i.e. to their distinct knots) specify the local knot vectors for each
control point of T-mesh (node) in the same order as the nodes have been specified for
the element. The local knot vector in a particular direction has :math:`p+2` entries,
where the :math:`p` is the polynomial degree in that direction. :param:`weights` holds
one weight per control point.

**Keyword**: :param:`bsplineplanestrainelement`  **Parameters**: None.
B-spline patch in plane strain (material mode PlaneStrain, two dofs
:math:`u, v` per control point). The volume of the integration points is multiplied by
the thickness of the cross section.

**Keyword**: :param:`bsplineaxisymelement`  **Parameters**: None.
Axisymmetric B-spline patch: :math:`x` is the radius and :math:`y` the axis of rotation,
as for the other axisymmetric elements. The hoop stretch is :math:`1 + u/r`, the
material is evaluated in the full 3d mode and the volume is taken per radian
(:math:`r \, \mathrm{d}r \, \mathrm{d}z`, the cross section thickness is not used).
Only the large strain formulation is implemented (:param:`nlgeo` 1).

IGA surface elements
--------------------

The surface elements carry a surface energy (surface tension, a surface membrane, and
for second gradient materials also bending) on the boundary of a two-dimensional
B-spline patch. Each of them is a planar B-spline curve: :param:`degree` has a single
entry, :param:`knotvectoru` and :param:`nodes` describe the curve, and to couple it with
the bulk it has to use the control points (nodes) of the patch boundary with the same
degree and knots. The cross section is a :param:`simplesurfacecs` with a surface
material. Only the large strain formulation is available, so :param:`nlgeo` has to be 1.

**Keyword**: :param:`bsplineplanestrainsurfaceelement`  **Parameters**: None.
Surface of a plane strain body; the surface energy is per unit length of the curve (and
unit thickness, the cross section thickness is not used, so use thickness 1 in the
bulk).

**Keyword**: :param:`bsplineaxisymsurfaceelement`  **Parameters**: None.
Surface of revolution of the curve (:math:`x` radius, :math:`y` axis); measured per
radian, with the hoop stretch :math:`1 + u/r`.

The normal of the curve is :math:`\mathbf{n} = (t_y, -t_x)/|\mathbf{t}|` with the
tangent :math:`\mathbf{t} = \partial \mathbf{X}/\partial \xi`, so it points to the
right of the direction given by the order of the nodes (outwards for a boundary
traversed counterclockwise). Isotropic surface energies do not depend on this
orientation; the sign of the curvature, and so a prescribed spontaneous curvature of a
bending material, does.

The end of a surface curve can be clamped by :param:`igaclampboundarycondition` or
:param:`igaclampboundaryconditionlagrange`, see the boundary conditions in the input
manual.

Second gradient materials
-------------------------

The elements :param:`bsplineplanestrainelement`, :param:`bsplineaxisymelement`,
:param:`bsplineplanestrainsurfaceelement` and :param:`bsplineaxisymsurfaceelement`
support materials depending also on the gradient of the deformation gradient (e.g.
:param:`gradientsimopistermat` and :param:`gradientmooneyrivlin` in the bulk, bending of
:param:`surfacetensionmat2` on the surfaces). Such a material is detected automatically
from the cross section of the element; an element that does not support the second
gradient stops with an error. The second gradient needs :math:`C^1` continuity of the
interpolation, i.e. degree at least 2 without repeated internal knots; only the degree
is checked.

Export to VTK
-------------

The composite VTK export (:param:`vtkxml`) is implemented for
:param:`bsplineplanestrainelement`, :param:`bsplineaxisymelement` and both surface
elements: every knot span is sampled in :param:`numexp` points per parametric
direction. Gauss point values (:param:`ipvars`) are exported at the true positions of
the integration points.

Tests/Examples: `tests/regression/sm/ex-bspline-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-01.in>`_,
`tests/regression/sm/ex-bspline-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-02.in>`_,
`tests/regression/sm/ex3d-nurbs-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-01.in>`_,
`tests/regression/sm/ex3d-nurbs-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-02.in>`_,
`tests/regression/sm/iga-surface-reactions.in`,
`tests/regression/sm/iga-element-material.in`,
`tests/regression/sm/iga-tspline-nip.in`.
