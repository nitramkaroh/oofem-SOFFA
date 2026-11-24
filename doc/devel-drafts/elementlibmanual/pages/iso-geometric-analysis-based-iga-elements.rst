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
The knot values are divided by the length of the knot vector but not shifted, so the
parameter range is :math:`[0, 1]` when the first knot is 0 (which also gives the scale of
the clamp parameter :param:`locpos`).
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

**Keyword**: :param:`nurbsplanestresselement`  **Parameters**: :elemparam:`weights{ra}`
One weight per control point, finite and positive.

**Keyword**: :param:`nurbs3delement`  **Parameters**: :elemparam:`weights{ra}`
One weight per control point, finite and positive.

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
one weight per control point (finite and positive).

**Keyword**: :param:`bsplineplanestrainelement`  **Parameters**: None.
B-spline patch in plane strain (material mode PlaneStrain, two dofs
:math:`u, v` per control point). The volume of the integration points is multiplied by
the thickness of the cross section.

**Keyword**: :param:`bsplineaxisymelement`  **Parameters**: None.
Axisymmetric B-spline patch: :math:`x` is the radius and :math:`y` the axis of rotation,
as for the other axisymmetric elements. The hoop stretch is :math:`1 + u/r`, the
material is evaluated in the full 3d mode and the volume is taken per radian
(:math:`r \, \mathrm{d}r \, \mathrm{d}z`, the cross section thickness is not used).
Only the large strain formulation is implemented, so :param:`nlgeo` has to be 1.

IGA surface elements
--------------------

The surface elements carry a surface energy (surface tension, a surface membrane, and
for second gradient materials also bending). The first two are planar B-spline curves
on the boundary of a two-dimensional B-spline patch: :param:`degree` has a single
entry, :param:`knotvectoru` and :param:`nodes` describe the curve, and to couple it with
the bulk it has to use the control points (nodes) of the patch boundary with the same
degree and knots. The other two are a B-spline and a NURBS surface in space. The cross section is a
:param:`simplesurfacecs` with a surface material. Only the large strain formulation is
available, so :param:`nlgeo` has to be 1.

**Keyword**: :param:`bsplineplanestrainsurfaceelement`  **Parameters**: None.
Surface of a plane strain body; the surface energy is per unit length of the curve (and
unit thickness, the cross section thickness is not used, so use thickness 1 in the
bulk).

**Keyword**: :param:`bsplineaxisymsurfaceelement`  **Parameters**: None.
Surface of revolution of the curve (:math:`x` radius, :math:`y` axis); measured per
radian, with the hoop stretch :math:`1 + u/r`.

**Keyword**: :param:`bspline3dsurfaceelement`  **Parameters**: None.
B-spline surface in space: two parametric directions (:param:`degree` with two entries,
each at least 1, :param:`knotvectoru` and :param:`knotvectorv`), control points with
three coordinates and the degrees of freedom :math:`u, v, w`; the surface energy is per
unit area of the reference surface. Its normal is
:math:`\mathbf{N} = \mathbf{A}_1 \times \mathbf{A}_2 / |\mathbf{A}_1 \times \mathbf{A}_2|`
with the tangents :math:`\mathbf{A}_1 = \partial \mathbf{X}/\partial \xi` and
:math:`\mathbf{A}_2 = \partial \mathbf{X}/\partial \eta`. A bending material (second
gradient) needs degree at least 2 and :math:`C^1` continuity in both directions (see
below); the anisotropic surface material and the clamp are rejected. As with the other
IGA elements, the loads are applied as nodal loads in the control points.

**Keyword**: :param:`nurbs3dsurfaceelement`  **Parameters**: :elemparam:`weights{ra}`
NURBS surface in space: as :param:`bspline3dsurfaceelement`, with one weight per control
point (finite and positive, in the order of the nodes; only their ratios matter, they are
divided by the largest one); the geometry and the displacements are interpolated by the
rational basis functions. They represent conic sections exactly, e.g. the circles of a
surface of revolution: a quarter circle of radius :math:`r` in one knot span of degree 2
has the control points :math:`(r, 0)`, :math:`(r, r)`, :math:`(0, r)` with the weights
:math:`1, 1/\sqrt{2}, 1` (the VTK export samples this exact surface and connects the
samples by flat cells). A surface of revolution can thus be compared with
:param:`bsplineaxisymsurfaceelement` on the same generatrix: for an axisymmetric
deformation the generalized reactions of the two models agree up to the Gauss quadrature of
the arc length (regression test ``iga-nurbs-axisym-checks.sh``). A bending material needs
:math:`C^1` continuity (see below), which neither the usual full circle of degree 2 with
double knots at the quarters nor a closed surface across its seam has; model a quarter or a
half with symmetry conditions instead. For bending these include the slope: at a coordinate
plane of symmetry met perpendicularly by the surface (as by the quarter circle), the control
points next to the edge get the displacements of the edge points in the directions of the
plane (simple slave dofs, ``doftype`` and ``mastermask`` of the nodes), besides the zero
normal displacement of the edge.

The surface gradient uses the reciprocal tangents
:math:`\mathbf{A}^1 = (\mathbf{A}_2 \times \mathbf{N})/J` and
:math:`\mathbf{A}^2 = (\mathbf{N} \times \mathbf{A}_1)/J`, where
:math:`J = |\mathbf{A}_1 \times \mathbf{A}_2|`. This avoids cancellation in the
inverse metric on strongly skew patches. A zero or non-finite reference area is rejected.
The second gradient is built from the tangential Hessian of the basis functions,
:math:`(N_{,\alpha\beta} - \mathbf{X}_{,\alpha\beta} \cdot \nabla_s N)\, \mathbf{A}^\alpha \otimes \mathbf{A}^\beta`
(the Christoffel symbols of the reference surface included). Its normal part is left out,
as on the curves: the bending material uses the second gradient only through the current
normal and tangential tensors, which the normal part does not reach.

The normal of the curves is :math:`\mathbf{n} = (t_y, -t_x)/|\mathbf{t}|` with the
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
:param:`bsplineplanestrainsurfaceelement`, :param:`bsplineaxisymsurfaceelement`,
:param:`bspline3dsurfaceelement` and :param:`nurbs3dsurfaceelement`
support materials depending also on the gradient of the deformation gradient (e.g.
:param:`gradientsimopistermat` and :param:`gradientmooneyrivlin` in the bulk, bending of
:param:`surfacetensionmat2` on the surfaces). Such a material is detected automatically
from the cross section of the element; an element that does not support the second
gradient stops with an error. The second gradient needs :math:`C^1` continuity of the
interpolation, i.e. degree at least 2 and internal knot multiplicities at most degree
minus 1 in every parametric direction; the six elements check it and stop with an error
otherwise. It is implemented only for :param:`nlgeo` 1.

Export to VTK
-------------

The composite VTK export (:param:`vtkxml`) is implemented for the B-spline and NURBS
elements (not for :param:`tsplineplanestresselement`): every knot span is sampled in
:param:`numexp` points per parametric direction, connected by lines on the curves, by
quadrilaterals on the patches with two parametric dimensions (planar or in space) and by
hexahedra on :param:`nurbs3delement`. A hexahedron with a negative volume (a left-handed
parametrization) is written with the reversed vertex order. The export checks the
Jacobian of every hexahedron at its corners, Gauss points and centre: if it changes sign
or vanishes inside, the export points are too coarse for the patch (e.g. a strongly
twisted one) and a warning asks for a larger :param:`numexp`. Exported are the positions, the primary variable
DisplacementVector and,
for the surface elements with a bending material, the internal variable
IST_CurvatureTensor (the current curvature :math:`\mathbf{K} = -\mathbf{n} \cdot \nabla\mathbf{F}`,
on the curves and on the surface in space evaluated at the export point). Other primary and internal variables and cell variables stop
with an error. Gauss point values (:param:`ipvars`) are exported at the true positions
of the integration points.

Tests/Examples: `tests/regression/sm/ex-bspline-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-01.in>`_,
`tests/regression/sm/ex-bspline-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-02.in>`_,
`tests/regression/sm/ex3d-nurbs-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-01.in>`_,
`tests/regression/sm/ex3d-nurbs-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-02.in>`_,
`tests/regression/sm/iga-surface-reactions.in`,
`tests/regression/sm/iga-surface3d-curved.in`,
`tests/regression/sm/iga-surface3d-quad.in`,
`tests/regression/sm/iga-surface3d-extruded.in`,
`tests/regression/sm/iga-surface3d-skew.in`,
`tests/regression/sm/iga-surface3d-nonlinear.in`,
`tests/regression/sm/iga-surface3d-bending.in`,
`tests/regression/sm/iga-surface3d-bending-extruded.in`,
`tests/regression/sm/iga-surface-bending-curve.in`,
`tests/regression/sm/iga-element-material.in`,
`tests/regression/sm/iga-tspline-nip.in`.
