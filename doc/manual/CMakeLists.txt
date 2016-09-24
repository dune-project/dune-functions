set(SOURCES
  ../../examples/stokes-taylorhood.cc
)

set(IMAGES
  gfx/driven_cavity.pdf
  gfx/driven_cavity_result.png
  gfx/taylor_hood_tree.pdf
)

if(LATEX_USABLE)
    dune_add_latex_document(dune-functions-manual.tex
      INPUTS ${SOURCES}
      IMAGES ${IMAGES})
endif(LATEX_USABLE)
