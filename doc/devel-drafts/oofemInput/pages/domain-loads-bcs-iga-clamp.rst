.. To be inserted into doc/oofemInput/pages/domain-loads-bcs.rst, section "Constraints and couplings".

IGAClampBoundaryCondition, IGAClampBoundaryConditionLagrange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clamp a planar IGA curve, i.e. a :param:`bsplineplanestrainsurfaceelement` or a
:param:`bsplineaxisymsurfaceelement` (see the IGA elements in the element library
manual), at one value :math:`\xi_c` of its parameter. The condition

.. math::

   \mathbf{n} \cdot \frac{\partial \mathbf{u}}{\partial \xi} (\xi_c) = 0

with the normal :math:`\mathbf{n}` of the reference curve keeps the direction of the
tangent at :math:`\xi_c` (the curve may stretch there but not rotate). It is linear in the
displacements and exact for large deformations. Typically it clamps an end of a surface
curve whose position is prescribed by ordinary boundary conditions.

:param:`IGAClampBoundaryCondition` adds the condition by a penalty,
:param:`IGAClampBoundaryConditionLagrange` by a Lagrange multiplier, which is an
additional degree of freedom introduced by this boundary condition.

.. record::

   :descitem:`IGAClampBoundaryCondition` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`numel{in}`
   :elemparam:`locpos{rn}` :elemparam:`bcknotspan{in}` :optelemparam:`p{rn}`

.. record::

   :descitem:`IGAClampBoundaryConditionLagrange` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`numel{in}`
   :elemparam:`locpos{rn}` :elemparam:`bcknotspan{in}`

**Parameters**

:param:`numel{in}`
    Number of the clamped IGA curve element (degree at least 1).

:param:`locpos{rn}`
    The curve parameter :math:`\xi_c` (not a physical coordinate). The element divides
    the knot values of :param:`knotvectoru` by their range (last minus first knot), and
    :param:`locpos` refers to these scaled values: for knots starting at 0 the parameter
    runs from 0 to 1, so 1 is the end of the curve also for the knots 0, 2.

:param:`bcknotspan{in}`
    Index of the knot span containing :param:`locpos`, counted from zero in the full
    knot vector including the repeated end knots: the span is
    :math:`[\xi_k, \xi_{k+1}]` with :math:`k` = :param:`bcknotspan`, it must have a
    nonzero length, and :param:`locpos` may lie at either of its ends (so an internal
    knot can be reached from both sides). For degree :math:`p` the first span has the
    index :math:`p`; e.g. the curve with the knots 0, 0.5, 1 and degree 2 (knot vector
    0 0 0 0.5 1 1 1) has the spans 2 and 3. The basis functions and the normal are
    evaluated in this span.

:optparam:`p{rn}`
    Penalty factor of :param:`IGAClampBoundaryCondition`, default :math:`10^6`. The
    condition uses the derivative with respect to the curve parameter, so the effective
    stiffness also depends on the parametrization.

:param:`loadTimeFunction{in}`
    Required by the record; the clamp itself does not depend on time.

The location is checked when the input is read (element type and degree, the knot span
and that :param:`locpos` is a finite number inside it).

Example (curve element 2 with the knots 0 and 1 of degree 2, clamped at its end):

.. code-block:: none

   igaclampboundarycondition 3 loadTimeFunction 1 numel 2 locpos 1.0 bcknotspan 2 p 100
