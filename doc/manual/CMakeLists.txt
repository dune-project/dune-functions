if(LATEX_USABLE)

  set(BASES_SOURCES
    ../../examples/interpolation.cc
    ../../examples/stokes-taylorhood.cc
  )

  set(BASES_IMAGES
    gfx/driven_cavity.pdf
    gfx/driven_cavity_result.png
    gfx/febasis_interface_schematic.pdf
  )

  dune_add_latex_document(dune-functions-bases.tex
    BIBFILES dune-functions-bases.bib
    INPUTS ${BASES_SOURCES}
    IMAGES ${BASES_IMAGES})

  set(FUNCTIONS_SOURCES
    src/integration-test.cc
    gfx/timings_gcc.pgf
    gfx/timings_gcc_pgo.pgf
    gfx/timings_clang.pgf
  )

  set(FUNCTIONS_IMAGES
    gfx/functions_interface_schematic.pdf
  )

  dune_add_latex_document(dune-functions-functions.tex
    BIBFILES dune-functions-functions.bib
    INPUTS ${FUNCTIONS_SOURCES}
    IMAGES ${FUNCTIONS_IMAGES})

endif(LATEX_USABLE)
