set(SOURCES
  ../../examples/interpolation.cc
  ../../examples/stokes-taylorhood.cc
)

set(IMAGES
  gfx/driven_cavity.pdf
  gfx/driven_cavity_result.png
  gfx/febasis_interface_schematic.pdf
)

if(LATEX_USABLE)
    dune_add_latex_document(dune-functions-bases.tex
      BIBFILES dune-functions-bases.bib
      INPUTS ${SOURCES}
      IMAGES ${IMAGES})
endif(LATEX_USABLE)
