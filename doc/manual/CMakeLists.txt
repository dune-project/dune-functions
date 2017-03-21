set(SOURCES
  ../../examples/stokes-taylorhood.cc
)

set(IMAGES
  gfx/driven_cavity.pdf
  gfx/driven_cavity_result.png
  gfx/taylor_hood_tree.pdf
  gfx/febasis_interface_schematic.pdf
  gfx/taylor-hood-matrix-interleaved.pdf
  gfx/taylor-hood-matrix-lexicographic.pdf
)

if(LATEX_USABLE)
    dune_add_latex_document(dune-functions-manual.tex
      BIBFILES dune-functions-manual.bib
      INPUTS ${SOURCES}
      IMAGES ${IMAGES})
endif(LATEX_USABLE)
