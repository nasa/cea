# cmake/CompilerFlags.cmake
function(project_enable_fortran_std target)
  # GNU Fortran
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008>
  )

  # Intel classic ifort (Windows uses /, Linux often accepts - or -stand)
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:
      $<IF:$<PLATFORM_ID:Windows>,/stand:f08,-stand f08>
    >
  )

  # IntelLLVM ifx
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:
      $<IF:$<PLATFORM_ID:Windows>,/stand:f08,-stand f08>
    >
  )
endfunction()


