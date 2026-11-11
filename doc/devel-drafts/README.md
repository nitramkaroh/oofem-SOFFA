# Documentation drafts for the upstream port

Upstream OOFEM `devel` has its manuals as Sphinx pages (`doc/*/pages/*.rst`); this
branch still has the old LaTeX manuals. The pages here are written in the `devel`
format and are meant to be copied there when the IGA and surface energy work is ported
(they are not built on this branch).

| Draft | Target in `devel` |
|---|---|
| `elementlibmanual/pages/iso-geometric-analysis-based-iga-elements.rst` | replaces the page of the same name |
| `oofemInput/pages/domain-loads-bcs-iga-clamp.rst` | section for `doc/oofemInput/pages/domain-loads-bcs.rst`, "Constraints and couplings" |

The Lagrange surface elements and the surface materials will follow once their input
is settled.
