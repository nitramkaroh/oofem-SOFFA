.. _SurfaceElements:

Surface elements
================

The surface elements add a surface energy (surface tension, a surface membrane) to a body
modelled by ordinary bulk elements. A surface element has no volume of its own: it is
defined on the nodes of the bulk elements on the boundary, so that the surface deforms
together with the bulk, and it takes a surface material (see the surface materials in the
Material Library Manual) from a :param:`simplesurfacecs` cross section. The elements on
this page interpolate the displacements by Lagrange polynomials; the IGA surface elements
:param:`bsplineplanestrainsurfaceelement` and :param:`bsplineaxisymsurfaceelement` are
described with the IGA elements (:doc:`iso-geometric-analysis-based-iga-elements`).

Kinematics
----------

Only the large strain formulation is implemented, so :param:`nlgeo` has to be 1 (see
:doc:`geometric-nonlinear-analysis`); with :param:`nlgeo` 0 the computation stops with an
error. Each element works in an orthonormal frame
:math:`(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)` in which :math:`\mathbf{e}_3` is the unit
normal of the reference surface. The surface deformation gradient

.. math::

   \mathbf{F} = \mathbf{I}_s + \nabla_s \mathbf{u}, \qquad
   \mathbf{I}_s = \mathbf{I} - \mathbf{e}_3 \otimes \mathbf{e}_3,

with the surface gradient :math:`\nabla_s \mathbf{u}` of the displacement, has nine
components in this frame (in the order 11, 22, 33, 23, 13, 12, 32, 31, 21) and
:math:`\mathbf{F} \mathbf{e}_3 = \mathbf{0}`. The material receives :math:`\mathbf{F}` and
the normal :math:`\mathbf{N} = \mathbf{e}_3` (for a prestrained :param:`triangularsurface`
see below) and returns the first Piola-Kirchhoff stress; the surface energy is integrated
over the reference surface.

* The planar elements :param:`lineplanestrainsurface` and :param:`axilinesurface` are
  straight segments in the :math:`xy` plane. :math:`\mathbf{e}_1` is the unit tangent from
  node 1 to node 2, :math:`\mathbf{e}_2` the out-of-plane direction (the circumferential
  direction of the axisymmetric element) and :math:`\mathbf{e}_3 = (-e_{1y}, e_{1x})` the
  in-plane normal, which points to the left of the direction from node 1 to node 2. The
  out-of-plane stretch is 1 in plane strain and :math:`r/R = 1 + u_r/R` in axisymmetry, where
  :math:`x` is the radius and :math:`y` the axis of rotation. The energy is per unit length
  of the segment in plane strain (and unit thickness: the cross section thickness is not
  used, so use thickness 1 in the bulk) and per radian in axisymmetry.

* The spatial elements :param:`quadsurface`, :param:`quadsurfacequadratic`,
  :param:`triangularsurface` and :param:`triangularsurfacequadratic` are flat.
  :math:`\mathbf{e}_1` is the unit vector from node 1 to node 2,
  :math:`\mathbf{e}_3 = \mathbf{e}_1 \times (\mathbf{X}_3 - \mathbf{X}_1)` normalized (so the
  normal follows the order of the nodes by the right-hand rule) and
  :math:`\mathbf{e}_2 = \mathbf{e}_3 \times \mathbf{e}_1`. The geometry is given by the
  coordinates of the nodes in the plane of :math:`\mathbf{e}_1, \mathbf{e}_2`, so all nodes
  have to lie in one plane. The energy is per unit area.

The unknowns are the displacements of the nodes in the global coordinate system; the
elements rotate them to their frame. The isotropic surface energies do not depend on the
orientation of the normal.

These elements are only :math:`C^0` continuous and provide no second gradient of the
deformation, so they reject a material with bending stiffness (:param:`surfacetensionmat2`,
even with a zero bending modulus); bending is available on the IGA surface elements of
degree at least 2.

Planar surface elements
-----------------------

.. table:: lineplanestrainsurface element summary
   :name: lineplanestrainsurfacesummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lineplanestrainsurface                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Straight two-node surface segment of a plane strain body                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in 2 points by default, which the parameter :param:`NIP` changes.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | ``tests/regression/sm/surface-membrane-parameters.in``,                                      |
   |                          | ``tests/regression/sm/surface-arclength-gamma.in``,                                          |
   |                          | ``tests/regression/sm/iga-surface-reactions.in``                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: axilinesurface element summary
   :name: axilinesurfacesummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | axilinesurface                                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Straight two-node segment of an axisymmetric surface (surface of revolution)                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in 2 points by default, which the parameter :param:`NIP` changes.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | ``tests/regression/sm/iga-surface-reactions.in``                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+

Spatial surface elements
------------------------

The node numbering of the spatial elements is the same as that of the plane stress elements
with the same interpolation (see :ref:`Planestress2dfig`, :ref:`qplanstrssfig`,
:ref:`TrPlanestressfig` and :ref:`qtrplanstressfig`).

.. table:: quadsurface element summary
   :name: quadsurfacesummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | quadsurface                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Flat four-node quadrilateral surface element in 3D                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Bilinear approximation of displacements and geometry.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in 4 points by default, which the parameter :param:`NIP` changes.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | \-                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: quadsurfacequadratic element summary
   :name: quadsurfacequadraticsummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | quadsurfacequadratic                                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Flat eight-node quadrilateral surface element in 3D                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic (serendipity) approximation of displacements and geometry.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in 16 points by default, which the parameter :param:`NIP` changes.         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | \-                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: triangularsurface element summary
   :name: triangularsurfacesummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | triangularsurface                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Flat three-node triangular surface element in 3D, optionally prestrained                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`prestrain{in}` :optelemparam:`f0{ra}`                 |
   |                          | :optelemparam:`normal{ra}` :optelemparam:`anisodir{ra}`                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   |                          |                                                                                              |
   |                          | :param:`prestrain`: 1 switches the prestrain on (default 0), :param:`f0` (required then) is  |
   |                          | the prestrain deformation gradient, nine components in the global coordinate system in the   |
   |                          | order 11, 22, 33, 23, 13, 12, 32, 31, 21, and :param:`normal` (three components, global      |
   |                          | coordinates, only with the prestrain) the normal of the original configuration, which is     |
   |                          | only checked; see below.                                                                     |
   |                          |                                                                                              |
   |                          | :param:`anisodir`: direction for the anisotropic surface material (three components, global  |
   |                          | coordinates), used instead of the one the material constructs; see the material.             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in one point by default, which the parameter :param:`NIP` changes.         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | ``tests/regression/sm/surface-prestrain-normal.in``                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: triangularsurfacequadratic element summary
   :name: triangularsurfacequadraticsummary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | triangularsurfacequadratic                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Flat six-node triangular surface element in 3D                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: number of integration points.                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in 4 points by default, which the parameter :param:`NIP` changes.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | A :param:`simplesurfacecs` cross section with a surface material; the thickness is not used. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 1 (0 stops with an error).                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | \-                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+

Prestrain of triangularsurface
------------------------------

With :param:`prestrain` 1 the mesh is not the stress-free configuration of the surface but
the image of an original configuration under the constant deformation gradient
:math:`\mathbf{F}_0` (:param:`f0`, one for the element). The material then works with
:math:`\mathbf{F}\mathbf{F}_0` and the normal :math:`\mathbf{N}` of the original
configuration, and the energy is :math:`W(\mathbf{F}\mathbf{F}_0)/J_0` per unit area of the
mesh, i.e. :math:`W` per unit area of the original configuration. Only the part of
:math:`\mathbf{F}_0` which maps onto the element plane counts: with
:math:`\mathbf{F}_{0s} = \mathbf{I}_s \mathbf{F}_0` in the frame of the element,

.. math::

   J_0 = |\mathrm{cof}(\mathbf{F}_{0s})^T \mathbf{e}_3|, \qquad
   \mathbf{N} = \mathrm{cof}(\mathbf{F}_{0s})^T \mathbf{e}_3 / J_0,

so that :math:`\mathbf{F}\mathbf{F}_0\mathbf{N} = \mathbf{0}`; :math:`J_0` has to be
positive. The original configuration may thus be tilted out of the element plane. The
optional :param:`normal` is not used in the computation: the element stops with an error
when it is not :math:`\pm\mathbf{N}`, which checks that :param:`f0` describes the intended
original configuration.
