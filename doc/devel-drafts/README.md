# Documentation drafts for the upstream port

Upstream OOFEM `devel` has its manuals as Sphinx pages (`doc/*/pages/*.rst`); this
branch still has the old LaTeX manuals. The pages here are written in the `devel`
format and are meant to be copied there when the IGA and surface energy work is ported
(they are not built on this branch).

| Draft | Target in `devel` |
|---|---|
| `elementlibmanual/pages/iso-geometric-analysis-based-iga-elements.rst` | replaces the page of the same name |
| `oofemInput/pages/domain-loads-bcs-iga-clamp.rst` | section for `doc/oofemInput/pages/domain-loads-bcs.rst`, "Constraints and couplings" |
| `elementlibmanual/pages/surface-elements.rst` | new page, in the toctree of `doc/elementlibmanual/pages/sec-elements-for-structural-analysis-sm-module.rst` after `iso-geometric-analysis-based-iga-elements` |
| `matlibmanual/pages/surface-materials.rst` | new page, in the toctree of `doc/matlibmanual/pages/sec-material-models-for-structural-analysis.rst` after `large-strain-master-material` |
| `oofemInput/pages/domain-cross-sections-simplesurfacecs.rst` | section for `doc/oofemInput/pages/domain-cross-sections.rst`, after "Integral cross section with constant properties" |

The tests are referred to as `tests/regression/sm/...`, where they go in `devel`; on this
branch they are in `tests/sm/`.
