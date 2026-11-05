! Stand-in for the MPI Fortran header, used only by the FORD docs build
! (fpm.toml: include). The compiled code uses the real mpif.h.
! FORD fails to parse an included file with no statements, hence the line below.
      implicit none
