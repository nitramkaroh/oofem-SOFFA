.. _SurfaceMaterials:

Surface materials
=================

Hyperelastic materials of the surface elements (the Lagrange surface elements and the IGA
surface elements, see the Element Library Manual), which get them from the cross section
:param:`simplesurfacecs`. They are formulated for large deformations only, so the elements
have to use :param:`nlgeo` 1.

The element hands over the surface deformation gradient :math:`\mathbf{F}` (nine components
in the order 11, 22, 33, 23, 13, 12, 32, 31, 21) and the unit normal :math:`\mathbf{N}` of
the reference surface, both in the same orthonormal frame and with
:math:`\mathbf{F}\mathbf{N} = \mathbf{0}`. With the surface identity
:math:`\mathbf{I}_s = \mathbf{I} - \mathbf{N} \otimes \mathbf{N}`, the area ratio
:math:`J = |\mathrm{cof}(\mathbf{F}) \mathbf{N}|` and the norm
:math:`|\mathbf{F}| = \sqrt{\mathbf{F} : \mathbf{F}}`, the energy :math:`W` below is per unit
area of the reference surface and the material returns the first Piola-Kirchhoff stress
:math:`\mathbf{P} = \partial W / \partial \mathbf{F}`. The surface Cauchy stress
:math:`\mathbf{P}\mathbf{F}^T / J` is stored as the Cauchy stress tensor of the integration
point, in the frame of :math:`\mathbf{F}`.

Membrane
--------

The energy of the material :param:`surfacetensionmat` is

.. math::
   :label: surfaceMembraneEnergy

   W(\mathbf{F}) = \gamma J + \alpha |\mathbf{F}| + \frac{1}{2} \alpha_s (J - 1)^2
   + \mu \left( \frac{1}{2} |\mathbf{F}|^2 - 1 - \ln J \right) + \frac{\lambda}{2} (\ln J)^2 .

The terms :math:`\gamma J` (surface tension) and :math:`\alpha |\mathbf{F}|` form the
surface-polyconvex model of M. Horák, M. Šmejkal, M. Kružík, Surface-polyconvex models for
soft elastic solids, arXiv:2503.24294 (2025), eq. (90); they are positively one-homogeneous in
:math:`\mathbf{F}`, which the existence theory of that model needs. The area term with
:math:`\alpha_s` and the neo-Hookean membrane with :math:`\mu, \lambda` are not, so they lie
outside that theory. This material has no bending stiffness and rejects a nonzero
:param:`omega1`; bending needs :param:`surfacetensionmat2`.

.. table:: Surface membrane - summary.
   :name: SurfaceTensionMat_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic surface membrane (surface tension and surface elasticity)                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`surfacetensionmat` :elemparam:`num{in}` :elemparam:`d{rn}`                        |
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}` :optelemparam:`gaminit{rn}`         |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :optelemparam:`arclenghtpar{in}`   |
   |                    | :optelemparam:`checktangent{rn}`                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`gamma1` surface tension :math:`\gamma` (also :param:`gamma`, the older name)       |
   |                    | - :param:`gaminit` surface tension :math:`\gamma_0` added when :param:`gamma_ltf` is given,  |
   |                    |   :math:`\gamma(t) = \gamma_0 + \gamma_1 f(t)`                                               |
   |                    | - :param:`alpha1` coefficient :math:`\alpha` of :math:`|\mathbf{F}|` (also :param:`alpha`,   |
   |                    |   the older name)                                                                            |
   |                    | - :param:`alphas` coefficient :math:`\alpha_s` of the area term                              |
   |                    | - :param:`mu`, :param:`lambda` the neo-Hookean coefficients :math:`\mu`, :math:`\lambda`     |
   |                    | - :param:`gamma_ltf`, :param:`alpha_ltf`, :param:`alphas_ltf`, :param:`mu_ltf`,              |
   |                    |   :param:`lambda_ltf` time functions of the parameters, see below                            |
   |                    | - :param:`arclenghtpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d` are optional and zero by default.                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | All surface elements (Lagrange and IGA).                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | ``tests/regression/sm/surface-membrane-parameters.in``,                                      |
   |                    | ``tests/regression/sm/surface-arclength-gamma.in``,                                          |
   |                    | ``tests/regression/sm/surface-prestrain-normal.in``,                                         |
   |                    | ``tests/regression/sm/iga-surface-reactions.in``                                             |
   +--------------------+----------------------------------------------------------------------------------------------+

Membrane with bending stiffness
-------------------------------

The material :param:`surfacetensionmat2` adds the bending energy to the membrane
:eq:`surfaceMembraneEnergy`:

.. math::

   W(\mathbf{F}, \nabla \mathbf{F}) = W_{membrane}(\mathbf{F})
   + \frac{1}{2} \omega J (H - H_0)^2, \qquad
   H = \frac{1}{2} \mathrm{tr}\, \mathbf{b}, \qquad
   \mathbf{b} = -\mathbf{F}^{+T} \mathbf{K} \mathbf{F}^{+}, \qquad
   K_{jk} = -n_i \frac{\partial F_{ij}}{\partial X_k},

where :math:`\mathbf{n}` is the unit normal of the deformed surface (oriented by
:math:`\mathbf{N}`) and :math:`\mathbf{F}^{+}` the pseudoinverse of :math:`\mathbf{F}`.
:math:`H` is the mean curvature of the deformed surface and :math:`H_0` the spontaneous
curvature given by :param:`useh0`: 0 (default) :math:`H_0 = 0`, 1 the mean curvature of the
reference surface (the reference shape is then free of bending stress), 2 the value
:param:`h0`. The sign of the curvature, and so that of :param:`h0`, follows the orientation
of the normal. The bending term needs the second gradient of the deformation, so only the
IGA surface elements of degree at least 2 accept this material; all other surface elements
reject it, even with :math:`\omega = 0`.

.. table:: Surface membrane with bending stiffness - summary.
   :name: SurfaceTensionMat2_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic surface membrane with bending stiffness                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`surfacetensionmat2` :elemparam:`num{in}` :elemparam:`d{rn}`                       |
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}` :optelemparam:`gaminit{rn}`         |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :optelemparam:`omega1{rn}`         |
   |                    | :optelemparam:`omega_ltf{in}` :optelemparam:`useh0{in}` :optelemparam:`h0{rn}`               |
   |                    | :optelemparam:`arclenghtpar{in}` :optelemparam:`checktangent{rn}`                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - the parameters of :param:`surfacetensionmat`                                               |
   |                    | - :param:`omega1` bending modulus :math:`\omega`, :param:`omega_ltf` its time function       |
   |                    | - :param:`useh0` choice of the spontaneous curvature (0, 1 or 2, default 0)                  |
   |                    | - :param:`h0` spontaneous curvature :math:`H_0` for :param:`useh0` 2                         |
   |                    | - :param:`arclenghtpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d` are optional and zero by default.                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | :param:`bsplineplanestrainsurfaceelement` and :param:`bsplineaxisymsurfaceelement` of degree |
   |                    | at least 2.                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | ``tests/regression/sm/surface-membrane-parameters.in``,                                      |
   |                    | ``tests/regression/sm/surface-arclength-gamma.in``                                           |
   +--------------------+----------------------------------------------------------------------------------------------+

Anisotropic membrane
--------------------

The material :param:`anisosurfacetensionmat` adds two fibre terms to the membrane
:eq:`surfaceMembraneEnergy` (Horák, Šmejkal, Kružík 2025):

.. math::

   W(\mathbf{F}) = W_{membrane}(\mathbf{F}) + \eta |\mathbf{F}\mathbf{a}|
   + \beta \sqrt{|\mathbf{C}\mathbf{a}|}, \qquad \mathbf{C} = \mathbf{F}^T \mathbf{F},

with the fibre direction :math:`\mathbf{a}` in the reference configuration. The material
takes a vector :math:`\mathbf{d}`, either the :param:`anisodir` of the element or

.. math::

   \mathbf{d} = \frac{((\mathbf{X} - \mathbf{p}_1) \times (\mathbf{X} - \mathbf{p}_2)) \times \mathbf{n}_e}
   {|((\mathbf{X} - \mathbf{p}_1) \times (\mathbf{X} - \mathbf{p}_2)) \times \mathbf{n}_e|}

from the position :math:`\mathbf{X}` of the integration point, the points
:math:`\mathbf{p}_1` (:param:`p_1`) and :math:`\mathbf{p}_2` (:param:`p_2`) and the normal
:math:`\mathbf{n}_e` of the element, and uses
:math:`\mathbf{a} = \mathbf{N} \times \mathbf{d} + (\mathbf{N} \cdot \mathbf{d}) \mathbf{N}`.
Without prestrain, the constructed :math:`\mathbf{a}` is the unit tangential part of
:math:`(\mathbf{X} - \mathbf{p}_1) \times (\mathbf{X} - \mathbf{p}_2)`: for
:math:`\mathbf{p}_1, \mathbf{p}_2` on an axis of revolution, the circumferential direction
around this axis. The direction is built in 3D and rotated with the frame of the element,
so only the spatial Lagrange surface elements are supported; other elements stop with an
error.

.. table:: Anisotropic surface membrane - summary.
   :name: AnisoSurfaceTensionMat_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Surface membrane with one family of fibres                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`anisosurfacetensionmat` :elemparam:`num{in}` :elemparam:`d{rn}`                   |
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}` :optelemparam:`gaminit{rn}`         |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :elemparam:`p_1{ra}`               |
   |                    | :elemparam:`p_2{ra}` :optelemparam:`eta{rn}` :optelemparam:`eta_ltf{in}`                     |
   |                    | :optelemparam:`beta{rn}` :optelemparam:`beta_ltf{in}` :optelemparam:`arclenghtpar{in}`       |
   |                    | :optelemparam:`checktangent{rn}`                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - the parameters of :param:`surfacetensionmat` (:param:`arclenghtpar` only 0 or 1)           |
   |                    | - :param:`p_1`, :param:`p_2` two different points (three components each, required) which    |
   |                    |   define the fibre direction                                                                 |
   |                    | - :param:`eta`, :param:`beta` the fibre coefficients :math:`\eta`, :math:`\beta`,            |
   |                    |   :param:`eta_ltf`, :param:`beta_ltf` their time functions                                   |
   |                    | - :param:`arclenghtpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d`, :param:`p_1` and :param:`p_2` are optional and zero by     |
   |                    | default.                                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | :param:`quadsurface`, :param:`quadsurfacequadratic`, :param:`triangularsurface`,             |
   |                    | :param:`triangularsurfacequadratic`.                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | \-                                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+

Time dependence of the parameters
---------------------------------

Each of the parameters :param:`gamma1`, :param:`alpha1`, :param:`alphas`, :param:`mu`,
:param:`lambda`, :param:`omega1`, :param:`eta` and :param:`beta` may have its own time
function, given by the parameter of the same name with the suffix ``_ltf`` (e.g.
:param:`mu_ltf` for :param:`mu`); the value :math:`p` of the parameter is then
:math:`p\,f(t)` with the value :math:`f(t)` of the time function, otherwise it is constant.
The surface tension is the exception: with :param:`gamma_ltf` it is
:math:`\gamma(t) = \gamma_0 + \gamma_1 f(t)` with :math:`\gamma_0` = :param:`gaminit` and
:math:`\gamma_1` = :param:`gamma1`.

Material arc-length method
--------------------------

With :param:`arclenghtpar` (spelled so in the input) 1 or 2, the material is driven by the
material arc-length method (solver ``matcalm`` of the static structural analysis), whose load
level :math:`\lambda` multiplies the surface tension, :math:`\gamma = \lambda \gamma_1`
(:param:`arclenghtpar` 1), or the bending modulus, :math:`\omega = \lambda \omega_1`
(:param:`arclenghtpar` 2, only :param:`surfacetensionmat2`). All other parameters keep their
values; the time functions and :param:`gaminit` are not used then. :param:`arclenghtpar` 0
(default) switches it off.

Check of the tangent stiffness
------------------------------

:param:`checktangent` > 0 compares, at every evaluation of the tangent stiffness, the
derivative :math:`\partial\mathbf{P}/\partial\mathbf{F}` (for :param:`surfacetensionmat2`
also the derivatives of the second order stress) with central finite differences and stops
with an error when the relative difference exceeds the given tolerance; a negative value
only writes the differences to the log, 0 (default) switches the check off. It is meant for
testing and makes the computation considerably slower.
