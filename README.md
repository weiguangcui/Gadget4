

GADGET-4
========


Documentation
=============

To build a local copy of the documentation:

- run `doxygen`

- open `doxygen/html/index.html` in a browser


Developers
==========

See the file DEVELOPERS.


LIMITS
======


* Maximum number of particles per MPI-rank is restricted to 2^31 ~ 2 Billion

* Maximum number of groups/subhalos per MPI-rank is restricted to 2^31 ~ 2 Billion

* The total number of particles, and the total number of groups/subhalos, can however
  be much larger than this if a sufficiently large number of MPI ranks is used. 
  
* Maximum particle number in a single group or subhalo is for the default setting  2^31 ~ 2 billion, 
  but this can be enlarged if needed by setting the option FOF_ALLOW_HUGE_GROUPLENGTH
  
* The maximum number of subhalos per group is restricted nevertheless to 2^31 ~ 2 billion.
  
* If the legacy output formats 1/2 are used, each block in the output is restricted to 2 Gbyte in 
  size or less. This also means that the maximum number of particles per single output file, and the
  maximum number of groups/subhalos per single group catalogue file can be at most of order 2 Billion,
  normally substantially less than that. While this restriction can be circumvented by splitting output
  over enough files, it provides one further motivation to use the HDF5 formay!
  
* Maximum number of lightcones: 256

* If more MPI-ranks than PMGRID are used, not all MPI tasks get anything to do in FFT calculations, which severely 
  compromises achievable work-load and memory balance. In this case, the column-based FFT algorithm is
  a way out (it incurs higher cost however, as twice as many transpose operations are needed)







