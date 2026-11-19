.. _SurfaceMaterials:

Surface materials
=================

Hyperelastic materials of the surface elements (the Lagrange surface elements and the IGA
surface elements, see the Element Library Manual), which get them from the cross section
:param:`simplesurfacecs`. They are formulated for large deformations only, so the elements
have to use :param:`nlgeo` 1.

The registered name :param:`isolesurf` is a placeholder for an unfinished linear elastic
surface material and is rejected on input. It does not provide a supported small-strain
alternative to these materials.

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
soft elastic solids, arXiv:2503.24294 (2025), eq. (90): as a function of the independent
variables :math:`\mathbf{F} \times \mathbf{N}` and :math:`\mathrm{cof}(\mathbf{F}) \mathbf{N}`
they are convex and positively one-homogeneous (surface polyconvexity), which the existence
theory of that model needs (as a function of :math:`\mathbf{F}` alone, :math:`J` is
homogeneous of degree 2). The area term with
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
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}`                                     |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :optelemparam:`arclengthpar{in}`   |
   |                    | :optelemparam:`checktangent{rn}`                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`gamma1` surface tension :math:`\gamma` (also :param:`gamma`, the older name)       |
   |                    | - :param:`alpha1` coefficient :math:`\alpha` of :math:`|\mathbf{F}|` (also :param:`alpha`,   |
   |                    |   the older name)                                                                            |
   |                    | - :param:`alphas` coefficient :math:`\alpha_s` of the area term                              |
   |                    | - :param:`mu`, :param:`lambda` the neo-Hookean coefficients :math:`\mu`, :math:`\lambda`     |
   |                    | - :param:`gamma_ltf`, :param:`alpha_ltf`, :param:`alphas_ltf`, :param:`mu_ltf`,              |
   |                    |   :param:`lambda_ltf` time functions of the parameters, see below                            |
   |                    | - :param:`arclengthpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d` are optional and zero by default.                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | All surface elements (Lagrange and IGA).                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | ``tests/sm/surface-membrane-parameters.in``,                                                 |
   |                    | ``tests/sm/surface-arclength-gamma.in``,                                                     |
   |                    | ``tests/sm/surface-prestrain-normal.in``,                                                    |
   |                    | ``tests/sm/iga-surface-reactions.in``                                                        |
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
IGA surface elements of degree at least 2 with at least :math:`C^1` continuity accept this
material. Each interior knot must have multiplicity at most the degree minus one. All other
surface elements reject it, even with :math:`\omega = 0`. The planar IGA curve uses the normal
obtained by rotating its oriented tangent clockwise; the Lagrange planar line uses the
opposite rotation. A change of orientation must also change the sign of a prescribed
nonzero :param:`h0` to preserve the same preferred curvature.

.. table:: Surface membrane with bending stiffness - summary.
   :name: SurfaceTensionMat2_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic surface membrane with bending stiffness                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`surfacetensionmat2` :elemparam:`num{in}` :elemparam:`d{rn}`                       |
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}`                                     |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :optelemparam:`omega1{rn}`         |
   |                    | :optelemparam:`omega_ltf{in}` :optelemparam:`useh0{in}` :optelemparam:`h0{rn}`               |
   |                    | :optelemparam:`arclengthpar{in}` :optelemparam:`checktangent{rn}`                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - the parameters of :param:`surfacetensionmat`                                               |
   |                    | - :param:`omega1` bending modulus :math:`\omega`, :param:`omega_ltf` its time function       |
   |                    | - :param:`useh0` choice of the spontaneous curvature (0, 1 or 2, default 0)                  |
   |                    | - :param:`h0` spontaneous curvature :math:`H_0` for :param:`useh0` 2                         |
   |                    | - :param:`arclengthpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d` are optional and zero by default.                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | :param:`bsplineplanestrainsurfaceelement` and :param:`bsplineaxisymsurfaceelement` of degree |
   |                    | at least 2.                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | ``tests/sm/surface-membrane-parameters.in``,                                                 |
   |                    | ``tests/sm/surface-arclength-gamma.in``                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+

Anisotropic membrane
--------------------

The material :param:`anisosurfacetensionmat` adds two fibre terms to the membrane
:eq:`surfaceMembraneEnergy` (Horák, Šmejkal, Kružík 2025):

.. math::

   W(\mathbf{F}) = W_{membrane}(\mathbf{F}) + \eta |\mathbf{F}\mathbf{a}|
   + \beta \sqrt{|\mathbf{C}\mathbf{a}|}, \qquad \mathbf{C} = \mathbf{F}^T \mathbf{F},

with the unit fibre direction :math:`\mathbf{a}` in the plane of the element (the energy is even
in :math:`\mathbf{a}`, so its sign does not matter). The direction is given in global coordinates
in one of two ways. The material rotates it into the frame of the element, keeps its part in the
element plane (the surface identity) and normalizes it:

- :param:`anisodir` of the element (only :param:`triangularsurface`) is a general direction for
  any fibre layout, e.g. computed for each element by a script.
- Otherwise :param:`p_1` and :param:`p_2` define an axis, a construction made for fibres wound
  around it. The vector

  .. math::

     \mathbf{m} = (\mathbf{X} - \mathbf{p}_1) \times (\mathbf{X} - \mathbf{p}_2)

  is normal to the plane through the integration point :math:`\mathbf{X}` and the axis, so it
  points along the circle through :math:`\mathbf{X}` around the axis. By default
  :math:`\mathbf{a}` is the normalized part of :math:`\mathbf{m}` in the element plane, fibres
  along the parallels; :param:`meridional` 1 rotates it by 90 degrees in that plane, fibres
  along the meridians.

On a surface of revolution about the axis (a sphere with its centre on the axis, a cylinder,
a cone) :math:`\mathbf{m}` is tangent, so the two options are exactly the parallels and the
meridians. On a mesh of flat elements this holds only approximately, and on other surfaces
:math:`\mathbf{m}` has a normal part: the meridional fibres then lie in the plane through the
axis and :math:`\mathbf{X}`, the default ones are perpendicular to them in the surface. On the
axis, or where the circle is normal to the element, there is no direction and the analysis stops
with an error; so does an :param:`anisodir` normal to the element.

The material is made for elements without prestrain: the direction is taken on the mesh, which is
the reference configuration only without prestrain, so a prestrained element stops with an error.
The direction is built in 3D and rotated with the frame of the element, so only the spatial
Lagrange surface elements are supported; other elements stop with an error.

.. table:: Anisotropic surface membrane - summary.
   :name: AnisoSurfaceTensionMat_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Surface membrane with one family of fibres                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`anisosurfacetensionmat` :elemparam:`num{in}` :elemparam:`d{rn}`                   |
   |                    | :optelemparam:`gamma1{rn}` :optelemparam:`gamma_ltf{in}`                                     |
   |                    | :optelemparam:`alpha1{rn}` :optelemparam:`alpha_ltf{in}` :optelemparam:`alphas{rn}`          |
   |                    | :optelemparam:`alphas_ltf{in}` :optelemparam:`mu{rn}` :optelemparam:`mu_ltf{in}`             |
   |                    | :optelemparam:`lambda{rn}` :optelemparam:`lambda_ltf{in}` :optelemparam:`p_1{ra}`            |
   |                    | :optelemparam:`p_2{ra}` :optelemparam:`meridional{in}` :optelemparam:`eta{rn}`               |
   |                    | :optelemparam:`eta_ltf{in}` :optelemparam:`beta{rn}` :optelemparam:`beta_ltf{in}`            |
   |                    | :optelemparam:`arclengthpar{in}` :optelemparam:`checktangent{rn}`                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - the parameters of :param:`surfacetensionmat` (:param:`arclengthpar` only 0 or 1)           |
   |                    | - :param:`p_1`, :param:`p_2` two different points on the axis of the fibres (three           |
   |                    |   components each, given together); required unless every element gives :param:`anisodir`    |
   |                    | - :param:`meridional` 1 puts the fibres along the meridians of the axis instead of the       |
   |                    |   parallels (default 0; needs :param:`p_1` and :param:`p_2`)                                 |
   |                    | - :param:`eta`, :param:`beta` the fibre coefficients :math:`\eta`, :math:`\beta`,            |
   |                    |   :param:`eta_ltf`, :param:`beta_ltf` their time functions                                   |
   |                    | - :param:`arclengthpar` load parameter of the material arc-length method, see below          |
   |                    | - :param:`checktangent` tolerance of the check of the tangent stiffness, see below           |
   |                    |                                                                                              |
   |                    | All parameters except :param:`d` are optional; :param:`p_1` and :param:`p_2` have no         |
   |                    | default, the others are zero by default.                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, through the surface elements and :param:`simplesurfacecs`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Elements           | :param:`quadsurface`, :param:`quadsurfacequadratic`, :param:`triangularsurface`,             |
   |                    | :param:`triangularsurfacequadratic`.                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | ``tests/sm/surface-tangent-check.in``, ``tests/sm/surface-input-checks.sh``                  |
   +--------------------+----------------------------------------------------------------------------------------------+

Time dependence of the parameters
---------------------------------

The time-function keywords are :param:`gamma_ltf` for :param:`gamma1`, :param:`alpha_ltf`
for :param:`alpha1`, :param:`alphas_ltf` for :param:`alphas`, :param:`mu_ltf` for :param:`mu`,
:param:`lambda_ltf` for :param:`lambda`, :param:`omega_ltf` for :param:`omega1`,
:param:`eta_ltf` for :param:`eta` and :param:`beta_ltf` for :param:`beta`.
The value :math:`p` of a parameter with a time function is then
:math:`p\,f(t)` with the value :math:`f(t)` of the time function, otherwise it is constant.
An initial value belongs to the time function: for example,
:math:`\gamma(t) = \gamma_0 + \gamma_1 t` is :param:`gamma1` 1 with the time function
:math:`f(t) = \gamma_0 + \gamma_1 t` (``PiecewiseLinFunction`` or ``UsrDefLTF``).

Material arc-length method
--------------------------

The input keyword is :param:`arclengthpar`. The earlier misspelling ``arclenghtpar`` is
rejected; older input files must replace it with ``arclengthpar``.
With :param:`arclengthpar` 1 or 2, the material is driven by the
material arc-length method (solver ``matcalm`` of the static structural analysis), whose load
level :math:`\lambda` multiplies the surface tension, :math:`\gamma = \lambda \gamma_1`
(:param:`arclengthpar` 1), or the bending modulus, :math:`\omega = \lambda \omega_1`
(:param:`arclengthpar` 2, only :param:`surfacetensionmat2`). All other parameters keep their
values; the time functions are not used then. :param:`arclengthpar` 0
(default) switches it off.

All materials with a nonzero :param:`arclengthpar` share the same global load level.
The time functions of materials without this setting are still evaluated normally.
External loads, including loads marked ``reference``, remain fixed in ``matcalm``;
the same convention is used for the equilibrium residual and the reported reactions.
The ordinary ``calm`` solver instead scales its external reference loads by its load level.

The material load vector is the derivative of the internal forces with respect to the
load level. It generally depends on the current geometry. With ``Psi`` zero, the step
length is measured in the displacement increments; ``hpc`` selects the degrees of freedom
used for this measure. For nonzero ``Psi``, its parameter contribution uses the material
load vector, also when ``hpc`` is given.

Continuing in :math:`\omega` changes the bending stiffness, whereas continuing in
:math:`\gamma` changes the surface tension. The former is useful for stiffness studies,
but needs a nonzero initial material load direction. In particular, when :math:`H=H_0`
the bending contribution and its load derivative vanish, so a stress-free reference
configuration is not a valid starting point for bending-only continuation. Such a zero
direction is rejected rather than producing an undefined predictor.

A failed step is not accepted as an equilibrium: after exhausting the internal restarts,
the solver restores its starting state and returns failure. The static structural model
then uses its configured time-step reduction strategy, or stops if no such strategy is active.

Check of the tangent stiffness
------------------------------

:param:`checktangent` > 0 compares, at every evaluation of the tangent stiffness, the
derivative :math:`\partial\mathbf{P}/\partial\mathbf{F}` (for :param:`surfacetensionmat2`
also the derivatives of the second order stress) with central finite differences and stops
with an error when the relative difference exceeds the given tolerance; a negative value
only writes the differences to the log, 0 (default) switches the check off. It is meant for
testing and makes the computation considerably slower. Non-finite compared values stop
the check with an error, including when a negative tolerance requests logging only.
