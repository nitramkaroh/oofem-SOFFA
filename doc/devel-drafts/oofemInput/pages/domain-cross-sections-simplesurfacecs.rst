.. To be inserted into doc/oofemInput/pages/domain-cross-sections.rst, after "Integral cross section with constant properties".

.. _SimpleSurfaceCS:

Surface cross section
---------------------

The cross section of the surface elements (the Lagrange surface elements and the IGA surface
elements, see the Element Library Manual), which carry a surface energy.

.. record::

   :descitem:`SimpleSurfaceCS` :elemparam:`num{in}` :optelemparam:`material{in}`
   :optelemparam:`set{in}`

**Parameters**

:optparam:`material{in}`
    Number of the material record, which has to be a surface material
    (:param:`surfacetensionmat`, :param:`surfacetensionmat2` or
    :param:`anisosurfacetensionmat`, see the Material Library Manual). Without it the
    material of the element (:param:`mat`) is used.

The record accepts the parameters of :ref:`SimpleCS`, but the surface elements integrate over
the reference surface and use none of the geometric properties (in particular not the
thickness).
