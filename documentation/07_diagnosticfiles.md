
Diagnostic outputs
==================


[TOC]

Aside from the verbose stdout log-file that is produced when the code
is run (take a detailed look!), various files in the output directory
provide useful metrics about the progress of a simulation, its memory
usage, work-load balance, and cpu consumption. Here, these files are
described in turn. (Make sure that your browser window is set wide
enough to avoid line-breaks in some of the wide tables.)


stdout                                            {#stdout}
======

After starting the code, the stdout will report a welcome message that
includes the git version hash-key and the date this corresponds
to. This uniquely identifies the last update/pull of the code from the
version control system. The output then lists all the compile time
options that were set. Note that these can be easily copied into an
empty `Config.sh` file to configure the source with the same
configuration. The code then reports the number of MPI ranks used, as
well as the detected node configuration and the pinning settings that
are applied, if any. This is followed by a detailed examination of the
memory available on the execution hosts (unless this is
disabled). After that, the parameterfile settings that are used for
the run are listed. Again, note that these can be easily pasted into
an empty parameter file to reproduce the settings used here, if
needed. The stdout-file thus always contains all the information
needed to reproduce a given run. Before the actual simulation begins,
some other information that is sometimes useful is reported as well,
such as the unit system that is used and the sizes of the most
important data structures used by the code.

During a run, the code outputs frequent log-message what is currently
done by GADGET-4. Usually, the corresponding lines begin with a
capitalized key-word that identifies the corresponding code part. For
example, lines beginning with "DOMAIN:" refer to information issued by
the domain decomposition, lines starting with "PM-PERIODIC:" will be identified with the periodic FFT-based computation of the
long-range gravitational force. It can be useful to filter the file
with grep for one of these phrases to get a clearer picture of what is
happening in a particular code part.

In case a problem should occur that forces the code to stop, you will
typically see a line starting with "Code termination on task" in the
output file, followed with information on which task, in which
function, and in which line number of a particular source file the
crash was triggered. This is followed by some information about the
nature of the problem. For example, such a message could look like
this:

~~~~~~~~~~~~~
Code termination on task=48, function restart(), file src/io/restart.c, line 207: RESTART: Restart file '../../L1/output//restartfiles/restart.48' not found.
~~~~~~~~~~~~~

All of this together hopefully provides a good clue about what may
have happened, and how the issue can be fixed.



info.txt                                          {#info}
========

The file with name `info.txt` in the output directory just contains a
list of all the timesteps. The last entry always holds the timestep
that is currently processed. For a running simulation, the command

~~~~~~~~~~~~~ 
tail info.txt
~~~~~~~~~~~~~

will thus inform you about the current progress of the
simulation. Typical output in this file looks like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sync-Point 26566, Time: 0.789584, Redshift: 0.26649, Systemstep: 5.84556e-05, Dloga: 7.40361e-05, Nsync-grv:   35730685, Nsync-hyd:          0

Sync-Point 26567, Time: 0.789642, Redshift: 0.266396, Systemstep: 5.84599e-05, Dloga: 7.40361e-05, Nsync-grv:  521313377, Nsync-hyd:          0

Sync-Point 26568, Time: 0.789701, Redshift: 0.266302, Systemstep: 5.84642e-05, Dloga: 7.40361e-05, Nsync-grv:   35760584, Nsync-hyd:          0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first number just counts and identifies all the timesteps in terms
of the synchronization points of the timestep hierarchy, which are the
places where force computations occur. The values given after
Time/Redshift are the current simulation times (time is the scale
factor in cosmological simulations). The Systemstep-values give the
time difference to the preceeding step, which is equal to by how much
the simulation as a whole has been advanced since the last
synchronization point. For cosmological integrations, this is
supplemented with the systemstep in logarithmic form in terms of the
ln of the scale factor.  Finally, the number of collisionless
particles that are in sync with the current system time (i.e. these
are the ones that have finished their timestep there and start a new
one) is reported. This is also done separately for the SPH particles.


timebins.txt                                          {#timebins}
============

The file `timebins.txt` in the output directory provides a much more
detailed account of the distribution of particles onto the timestep
hierarchy. A typical entry, created for every synchronization point,
of this file looks as follows:


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sync-Point 16521, Time: 0.934985, Redshift: 0.0695359, Systemstep: 6.92201e-05, Dloga: 7.40361e-05
Occupied timebins: gravity         sph          dt              cumul-grav   cumul-sph A D    avg-time  cpu-frac
    bin=17        12287592           0     0.001184577701         27292287           0           21.46     21.3%
    bin=16         8012442           0     0.000592288851         15004695           0           16.42     16.3%
 X  bin=15         4946259           0     0.000296144425          6992253           0 < *       13.71     27.2%
 X  bin=14         1809726           0     0.000148072213          2045994           0            7.22     28.6%
 X  bin=13          236268           0     0.000074036106           236268           0            0.84      6.7%
               ------------------------
Total active:      6992253           0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first line corresponds to information also included in
`info.txt`. This is followed with a table that shows the distribution
of all particles (listed unter "gravity") onto different timebins. The
different possible timestep sizes are identified via the timebin
number and the corresponding timestep size, which is listed in column
"dt". Also given is the distribution of gaseous SPH particles onto the
timebins. The columns "cumul-grav" and "cumul-sph" list the cumulative
numbers of particles in this timebin and below in the categories of
gravity and SPH calculations. The "X" symbols in the beginning mark
the timebins that are synchronized with the current system time. All
particles that are included in these timebins need a force calculation
at the current system time, hence the cumulative numbers reported for
the top timebin of this set, which is marked with a "<" sign in the
"A" column give an indication of the amount of computational work
required for this step. The value reported under "avg-time" represents
an estimated execution time of this current step, obtained by
averaging around five past executions of this timebin when it was
equally marked with a "<" sign. Note that sometimes these values can
be distorted a bit when the code had to do a special operation during
one of these last averaging steps, like computing a group catalogue,
or writing restart files. Also note that lower timebins need to be
executed more frequently than higher timebins. For example, there will
be twice as many executions of timebin 15 than of timebin 16 for every
execution of timebin 17. The last column represents the resulting
distribution of consumption of total CPU-time when the different
execution frequency of the steps is taken into account, based on the
values reported under `avg-time` . While this estimate is approximate
in nature, it gives a good idea about what cost is incurred in the
total simulation run-time due to the presence of certain timebins. In
particular, a situation where the lowest occupied timebins consume the
dominant fraction of the CPU time despite the fact that these are only
thinly populated normally indicates that the simulation does not run
very efficiently and only slowly makes progress.



cpu.txt                                           {#cpustat}
=======

In the file cpu.txt, you get detailed statistics about the total CPU
consumption measured in various parts of the code while it is
running. At each timestep, a table is added to this file, which
roughly looks like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 328, Time: 0.0485236, CPUs: 2496, HighestActiveTimeBin: 19
                          diff                 cumulative
total                   260.48  100.0%   63526.74  100.0%
  treegrav              213.41   81.9%   52556.59   82.7%
    treebuild             3.01    1.2%     931.09    1.5%
      insert              1.56    0.6%     524.64    0.8%
      branches            0.34    0.1%     113.39    0.2%
      toplevel            0.19    0.1%      72.99    0.1%
    treeforce           210.39   80.8%   51622.08   81.3%
      treewalk          174.82   67.1%   42617.69   67.1%
      treeimbalance      27.99   10.7%    6866.13   10.8%
      treefetch           0.53    0.2%     105.94    0.2%
      treestack           7.05    2.7%    2032.32    3.2%
  pm_grav                31.68   12.2%    6395.72   10.1%
  ngbtreebuild            0.00    0.0%       0.00    0.0%
  ngbtreevelupdate        0.11    0.0%      25.42    0.0%
  ngbtreehsmlupdate       0.00    0.0%       0.21    0.0%
  sph                     0.00    0.0%       0.00    0.0%
    density               0.00    0.0%       0.00    0.0%
      densitywalk         0.00    0.0%       0.00    0.0%
      densityfetch        0.00    0.0%       0.00    0.0%
      densimbalance       0.00    0.0%       0.00    0.0%
    hydro                 0.00    0.0%       0.00    0.0%
      hydrowalk           0.00    0.0%       0.00    0.0%
      hydrofetch          0.00    0.0%       0.00    0.0%
      hydroimbalance      0.00    0.0%       0.00    0.0%
  domain                 11.95    4.6%    2177.10    3.4%
  peano                   1.24    0.5%     291.91    0.5%
  drift/kicks             1.44    0.6%     312.75    0.5%
  timeline                0.07    0.0%      18.62    0.0%
  treetimesteps           0.00    0.0%       0.00    0.0%
  i/o                     0.00    0.0%      71.90    0.1%
  logs                    0.18    0.1%      43.06    0.1%
  fof                     0.00    0.0%       0.00    0.0%
    fofwalk               0.00    0.0%       0.00    0.0%
    fofimbal              0.00    0.0%       0.00    0.0%
  restart                 0.00    0.0%    1506.94    2.4%
  misc                    0.40    0.2%     126.52    0.2%
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two columns of measurements are provided. The entries under "diff"
measure the time difference in seconds relative to the last time
step. This total elapsed time reported in the row "total" is then
subdivided onto different code parts, as labelled. If the indention
level is increased, a further subdivision of the corresponding code
part is provided in terms of the subsequent indented items.  The
relative fractions of these various code parts are also reported as a
percentage of the total time of the step.

In addition, the values reported under "cumulative" give the same
analysis for the cumulative elapsed time since the start of the
simulation. Again, the first number gives absolute elapsed times, so
that the number reported under "total" is the total consumed time
since the start of the simulation. The count continues across
restarts, i.e. when the simulation is completed, this number gives a
faithful account of the total time (in seconds) needed to bring the
simulation to completion. Note that to get the full CPU-time
consumption in core-hours, this number has still to be multiplied with
the number of cores occupied (and to be divided by 3600 to convert
from seconds to hours, of course).

The data contained in the `cpu.txt` file is additionally output as a
column-separated file `cpu.csv` , where all the numbers for one
timestep appear in a one line entry. This can be used to more easily
make plots of the CPU time consumption in different code parts, if
this is desired.


domain.txt                                    {#domainbalance}
==========

The log file `domain.txt` receives a new entry whenever a new domain
decomposition is carried out by the code. A typical output for one
step may look like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOMAIN BALANCE, Sync-Point 31488, Time: 0.994978
Timebins:       Gravity       Hydro  cumulative      grav-balance       hydro-balance
 |bin=17     5882690219           0 10204440916  m  1.106 | 1.000       0.000 | 0.000
 |bin=16     2353095591           0  4321750697  m  1.525 | 1.000       0.000 | 0.000
>|bin=15     1425860069           0  1968655106  m  1.031 | 1.016  *    0.000 | 0.000
 |bin=14      503495925           0   542795037  m  1.070 | 1.031  *    0.000 | 0.000
 |bin=13       39299112           0    39299112  m  2.247 | 1.040  *    0.000 | 0.000
-------------------------------------------------------------------------------------
BALANCE,  LOAD:   1.003       0.000       1.003  WORK:      1.087               0.000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the first line indicates the current system time of the code, and
the following table reports all timebins and their occupancy with
gravity-only and hydrodynamic particles, as well as the total
cumulative occupancy up to the given timebin. The currently highest
synchronized timebin is marked with a "<" sign. The current timestep
distribution and the settings of the code (in particular the parameter
`ActivePartFracForNewDomainDecomp`), allow the code to tell when the
next domain decomposition will be done. In particular, in the above
example, this will not happen when timebins 14 or 13 are the highest
active timebin. Therefore, the current domain decomposition has
attempted to balance timebins 15, 14, and 13 simultaneously, which is
indicated by the `*` marker sign. Only the timebins marked by `*` need
to be balanced by the current domain decomposition. In this balancing,
the code attempts to reach a good balance for the particle load, the
gravity work-load, and the hydrodynamic work-load. For the latter two,
the algorithm takes into account that several steps need to be
executed until the next domain decomposition will occur, and that the
goal should be to minimize the overall execution time until then. This
for example means that a larger relative imbalance for bin 13 may be
acceptable if the absolute time for executing this step is reasonably
short.

The result of the domain decomposition is reported under
"grav-balance" and "hydro-balance", respectively. The first number
gives the work-imbalance if this step would be executed alone. This
imbalance factor is defined as the maximum (estimated) execution time
among the MPI-ranks divided by the average of the execution
times. Because the total work required per step is invariant under the
domain decomposition, this should be approximately independent of the
domain decomposition, so the goal is to push the maximum execution
time as close to the average as possible. The imbalance factor gives
the relative slowdown due to an imperfect work-load balance.

The values reported under grav-balance and hydro-balance inform about
the relative success of the domain decomposition to reach the desired
balance among multiple different quantities. The first number gives
the imbalance factor if only the current step would have to be
executed. This is however only the full story if another domain
decomposition will be carried out in the next step immediately. If
this is not the case, several steps need to be averaged appropriately,
and the relative slow-down of the residual imbalance in the present
step is reported as the second number, while the overall imbalance
over all simultaneously balanced steps is reported behind "WORK"
(first number in the case of gravity.) The second numbers in the table
above should add up to this number, as they give for every timebin
involved in the balancing the relative contribution they are
responsible for in this overall imbalance. Effectively, the "WORK"
factor tells how much faster the code may run if the domain
decomposition could be carried out perfectly for every involved
step. Note that timebins that are occupied with particles but do not
need to be balanced by the current domain decomposition may have a
large intrinsic imbalance, but this doesn't affect the run-time
behaviour at all, hence the second number is reported as a 1.0 in this
case.

Similarly, the hydrodynamic imbalance is reported in a second pair of
columns, yielding a further overall imbalance factor that is reported
as second number after "WORK".

Finally, there is an overall memory-load imbalance factor reported
behind "LOAD", which is also subdivided into gravity and
SPH-particles. This is another important metric as GADGET-4 attempts
to balance the work-load without allowing for significant memory
imbalance (because often simulations can be memory-bound). This is in
practice achieved by trying to push down all values reported under
LOAD and WORK simultaneously. The categories memory-load and work-load
in gravity and SPH is given equal weight in this context.


balance.txt                                    {#balance}
===========

The file `balance.txt` provides another quick look at the performance
and execution pattern of the code. It is important to view this in a
terminal window that is set wide enough to avoid extra line
wrapping. In this file, each step of the code is reported with a
single line, giving step number, total execution time, and the number
of active gravity and hydro particles. This is followed with a block
of characters of a total length representing the full execution time
of the step. Different code parts are represented by different
characters, and are filling this block with their corresponding
character in proportion to the time spent in this code part. A piece
of the the resulting output may then for example look something like
this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step=    761  sec=     9.406 Nsync-grv=  27292287 Nsync-hyd=         0  :r************************************************************************ttt|||||||||||||||||||||||||||||||+++++"??7---
Step=    762  sec=     0.373 Nsync-grv=     15778 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::::::::rrrDDDDDd.................*ttA--------------------------------------------------
Step=    763  sec=     8.025 Nsync-grv=  26662053 Nsync-hyd=         0  ::.**********************************************************************************tttt|||||||||||||||||++++++"??F----
Step=    764  sec=     0.371 Nsync-grv=     16226 Nsync-hyd=         0  :::::::::::::::::::::::::::::::::::::::::rrDDDDDD..................ttt7-------------------------------------------------
Step=    765  sec=     9.385 Nsync-grv=  27292287 Nsync-hyd=         0  :r************************************************************************ttt|||||||||||||||||||||||||||||||+++++"??F---
Step=    766  sec=     0.403 Nsync-grv=     16748 Nsync-hyd=         0  :::::::::::::::::::::::::::::::::::::rrrDDDDD................*tt0-------------------------------------------------------
Step=    767  sec=     8.049 Nsync-grv=  26662117 Nsync-hyd=         0  ::.**********************************************************************************ttt=|||||||||||||||||+++++""?KF----
Step=    768  sec=     0.419 Nsync-grv=     17202 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::::::rrDDDDD................tttF-------------------------------------------------------
Step=    769  sec=     9.399 Nsync-grv=  27292287 Nsync-hyd=         0  :r************************************************************************ttt|||||||||||||||||||||||||||||||+++++"??F---
Step=    770  sec=     0.407 Nsync-grv=     17695 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::::::rDDDDDd................tttF-------------------------------------------------------
Step=    771  sec=     8.016 Nsync-grv=  26662178 Nsync-hyd=         0  ::.**********************************************************************************tttt|||||||||||||||||++++++"??F----
Step=    772  sec=     0.430 Nsync-grv=     18152 Nsync-hyd=         0  :::::::::::::::::::::::::::::::::::rrDDDDDE...............tttF----------------------------------------------------------
Step=    773  sec=     9.517 Nsync-grv=  27292287 Nsync-hyd=         0  :r***********************************************************************ttt|||||||||||||||||||||||||||||||+++++"??7----
Step=    774  sec=     0.435 Nsync-grv=     18726 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::rrrDDDDD...............ttt7-----------------------------------------------------------
Step=    775  sec=     8.062 Nsync-grv=  26662000 Nsync-hyd=         0  ::.*********************************************************************************tttt||||||||||||||||||+++++"o?7-----
Step=    776  sec=     0.506 Nsync-grv=     19217 Nsync-hyd=         0  :::::::::::::::::::::::::::::rrrrrrrrDDDdE..............tt0-------------------------------------------------------------
Step=    777  sec=    19.647 Nsync-grv=  27292287 Nsync-hyd=         0  :::D*****************************ttttttt|||||||||||||||++++++"?0333333333333333555555555D6666666666666666666666F--------
Step=    778  sec=     0.445 Nsync-grv=     19759 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::rrDDDDD...............tttF------------------------------------------------------------
Step=    779  sec=     8.101 Nsync-grv=  26662013 Nsync-hyd=         0  ::.*********************************************************************************tttt|||||||||||||||||++++++"??F-----
Step=    780  sec=     0.373 Nsync-grv=     20232 Nsync-hyd=         0  :::::::::::::::::::::::::::::::::::::::rrrrDDDDDd..................ttttF------------------------------------------------
Step=    781  sec=     9.538 Nsync-grv=  27292287 Nsync-hyd=         0  :r**********************************************************************ttt=|||||||||||||||||||||||||||||||+++++"?KF----
Step=    782  sec=     0.366 Nsync-grv=     20811 Nsync-hyd=         0  ::::::::::::::::::::::::::::::::::::::::::rrDDDDDD..................*ttt7-----------------------------------------------
Step=    783  sec=     8.053 Nsync-grv=  26662006 Nsync-hyd=         0  ::.**********************************************************************************ttt||||||||||||||||||++++++"?)F----
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A key to the different symbols used for different code parts is
included in the beginning of the file `balance.txt`. The idea of this
output is to allow a quick graphical analysis of the execution patterns
of the code, in particular also to allow visual identification of
sudden changes of it. For example, normally the appearance of
additional timestep bins, or the dominance of certain code parts can
be readily inferred from this graphical text output. Likewise, the
occurrence of things like group finding or light-cone file output
tends to show up. Also imbalances in certain code parts are reported
separately by different symbols, so this can also be a way to tell
whether imbalances in certain places are particularly strong.

One should be cautious, however, to avoid over-interpreting this
graphical text output.  Because short and long steps are all stretched
to the same width in their corresponding output lines, short timesteps
(which may be comparatively unimportant for the total CPU budget) tend
to be overrepresented in this graphical representation. Note that by
filtering out certain timebins, this effect can be avoided and the
variation in the execution metrics of this particular step as the
simulation progresses can be monitored.


memory.txt                                                  {#memory}
==========

Another interesting diagnostic information about the simulation code
is contained in the file `memory.txt`. There, each time a new
high-watermark is reached in the total memory consumption on any of
the MPI-ranks, a new entry in the form of an extended table is
produced. An example for this table is reproduced below.

One of the most important numbers is the value reported behind
"Largest Allocation Without Generic". This is the minimum amount of
memory the code needed to run in the present configuration during this
timestep, excluding communication buffers. The latter are flexible in
size and will automatically adjust to the amount of free memory left
according to the `MaxMemSize` parameter.

For a more detailed view, the table tells the name (normally identical
to the variable name in the source code used to refer to the buffer)
of each allocated memory block. This is followed by its size in
MBytes, and the cumulative size of all allocated blocks up to this
point. In addition, the function, file name, and line number where this
particular block was allocated is given as well.

GADGET-4 organizes its internal memory handling in the form of a stack
in order to avoid memory fragmentation. The flag reported under "F"
shows whether the block has been explicitly allocated as movable (so
that previous blocks may be freed or resized), or whether this hasn't
been done by the source code. The number reported for "Task" is just
the MPI-rank on which this maximum allocation had occurred.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MEMORY:  Largest Allocation = 1542.21 Mbyte  |  Largest Allocation Without Generic = 301.181 Mbyte

-------------------------- Allocated Memory Blocks---- ( Step     5086 )------------------
Task    Nr F                  Variable      MBytes   Cumulative  Function|File|Linenumber
------------------------------------------------------------------------------------------
  12     0 0                               Exportflag      0.0002       0.0002  allocate_comm_tables()|src/data/allocate.c|37
  12     1 0                              Exportindex      0.0002       0.0005  allocate_comm_tables()|src/data/allocate.c|39
  12     2 0                          Exportnodecount      0.0002       0.0007  allocate_comm_tables()|src/data/allocate.c|40
  12     3 0                                     Send      0.0005       0.0012  allocate_comm_tables()|src/data/allocate.c|42
  12     4 0                                     Recv      0.0005       0.0017  allocate_comm_tables()|src/data/allocate.c|43
  12     5 0                               Send_count      0.0002       0.0020  allocate_comm_tables()|src/data/allocate.c|45
  12     6 0                              Send_offset      0.0002       0.0022  allocate_comm_tables()|src/data/allocate.c|46
  12     7 0                               Recv_count      0.0002       0.0024  allocate_comm_tables()|src/data/allocate.c|47
  12     8 0                              Recv_offset      0.0002       0.0027  allocate_comm_tables()|src/data/allocate.c|48
  12     9 0                         Send_count_nodes      0.0002       0.0029  allocate_comm_tables()|src/data/allocate.c|50
  12    10 0                        Send_offset_nodes      0.0002       0.0032  allocate_comm_tables()|src/data/allocate.c|51
  12    11 0                         Recv_count_nodes      0.0002       0.0034  allocate_comm_tables()|src/data/allocate.c|52
  12    12 0                        Recv_offset_nodes      0.0002       0.0037  allocate_comm_tables()|src/data/allocate.c|53
  12    13 0                          Tree.Send_count      0.0002       0.0039  allocate_comm_tables()|src/data/allocate.c|55
  12    14 0                         Tree.Send_offset      0.0002       0.0042  allocate_comm_tables()|src/data/allocate.c|56
  12    15 0                          Tree.Recv_count      0.0002       0.0044  allocate_comm_tables()|src/data/allocate.c|57
  12    16 0                         Tree.Recv_offset      0.0002       0.0046  allocate_comm_tables()|src/data/allocate.c|58
  12    17 1                                IO_Fields      0.0035       0.0082  init_field()|src/io/io_fields.c|91
  12    18 1                                IO_Fields      0.0123       0.0205  init_field()|src/io/io_fields.c|91
  12    19 0                             slab_to_task      0.0020       0.0225  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|43
  12    20 0                         slabs_x_per_task      0.0002       0.0227  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|58
  12    21 0                     first_slab_x_of_task      0.0002       0.0229  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|61
  12    22 0                         slabs_y_per_task      0.0002       0.0232  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|64
  12    23 0                     first_slab_y_of_task      0.0002       0.0234  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|67
  12    24 0                             slab_to_task      0.0020       0.0254  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|43
  12    25 0                         slabs_x_per_task      0.0002       0.0256  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|58
  12    26 0                     first_slab_x_of_task      0.0002       0.0259  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|61
  12    27 0                         slabs_y_per_task      0.0002       0.0261  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|64
  12    28 0                     first_slab_y_of_task      0.0002       0.0264  my_slab_based_fft_init()|src/pm/pm_mpi_fft.c|67
  12    29 0                                kernel[1]      8.0312       8.0576  pm_init_nonperiodic()|src/pm/pm_nonperiodic.c|482
  12    30 1                      def->ntype_in_files      0.0015       8.0591  read_files_driver()|src/io/io_fields.c|301
  12    31 1                                        P     98.3903     106.4494  allocate_memory()|src/data/allocate.c|71
  12    32 1                                     SphP      0.0001     106.4495  allocate_memory()|src/data/allocate.c|72
  12    33 1                  NextActiveParticleHydro      0.0001     106.4495  timebins_allocate()|src/time_integration/timestep.c|501
  12    34 1                       NextInTimeBinHydro      0.0001     106.4496  timebins_allocate()|src/time_integration/timestep.c|504
  12    35 1                       PrevInTimeBinHydro      0.0001     106.4496  timebins_allocate()|src/time_integration/timestep.c|507
  12    36 1                NextActiveParticleGravity      3.6441     110.0938  timebins_allocate()|src/time_integration/timestep.c|501
  12    37 1                     NextInTimeBinGravity      3.6441     113.7379  timebins_allocate()|src/time_integration/timestep.c|504
  12    38 1                     PrevInTimeBinGravity      3.6441     117.3820  timebins_allocate()|src/time_integration/timestep.c|507
  12    39 1                             D->StartList      0.0020     117.3839  domain_allocate()|src/domain/domain.cc|183
  12    40 1                               D->EndList      0.0020     117.3859  domain_allocate()|src/domain/domain.cc|184
  12    41 1                    D->FirstTopleafOfTask      0.0002     117.3861  domain_allocate()|src/domain/domain.cc|185
  12    42 1                      D->NumTopleafOfTask      0.0002     117.3864  domain_allocate()|src/domain/domain.cc|186
  12    43 1                              D->TopNodes      0.0234     117.4098  domain_allocate()|src/domain/domain.cc|187
  12    44 1                            D->TaskOfLeaf      0.0103     117.4200  domain_allocate()|src/domain/domain.cc|188
  12    45 1                       D->ListOfTopleaves      0.0103     117.4303  domain_decomposition()|src/domain/domain.cc|135
  12    46 1                                       PS     65.5936     183.0239  fof_fof()|src/fof/fof.cc|181
  12    47 1                                    Group      0.0610     183.0848  fof_fof()|src/fof/fof.cc|293
  12    48 1                                 SubGroup      2.0081     185.0930  subfind()|src/subfind/subfind.cc|319
  12    49 1                             D->StartList      0.0001     185.0931  domain_allocate()|src/domain/domain.cc|183
  12    50 1                               D->EndList      0.0001     185.0932  domain_allocate()|src/domain/domain.cc|184
  12    51 1                    D->FirstTopleafOfTask      0.0001     185.0933  domain_allocate()|src/domain/domain.cc|185
  12    52 1                      D->NumTopleafOfTask      0.0001     185.0933  domain_allocate()|src/domain/domain.cc|186
  12    53 1                              D->TopNodes      0.0011     185.0944  domain_allocate()|src/domain/domain.cc|187
  12    54 1                            D->TaskOfLeaf      0.0005     185.0949  domain_allocate()|src/domain/domain.cc|188
  12    55 1                       D->ListOfTopleaves      0.0005     185.0954  domain_decomposition()|src/domain/domain.cc|135
  12    56 1                                IndexList      2.5806     187.6760  subfind_processing()|src/subfind/subfind_processing.cc|166
  12    57 1                             D->StartList      0.0001     187.6761  domain_allocate()|src/domain/domain.cc|183
  12    58 1                               D->EndList      0.0001     187.6761  domain_allocate()|src/domain/domain.cc|184
  12    59 1                    D->FirstTopleafOfTask      0.0001     187.6762  domain_allocate()|src/domain/domain.cc|185
  12    60 1                      D->NumTopleafOfTask      0.0001     187.6763  domain_allocate()|src/domain/domain.cc|186
  12    61 1                              D->TopNodes      0.0004     187.6766  domain_allocate()|src/domain/domain.cc|187
  12    62 1                            D->TaskOfLeaf      0.0002     187.6768  domain_allocate()|src/domain/domain.cc|188
  12    63 1                       D->ListOfTopleaves      0.0002     187.6770  domain_decomposition()|src/domain/domain.cc|135
  12    64 1                                       sd      8.5197     196.1967  subfind_process_single_group()|src/subfind/subfind_processing.cc|236
  12    65 1                                     Head      2.1299     198.3266  subfind_process_single_group()|src/subfind/subfind_processing.cc|278
  12    66 1                                     Next      2.1299     200.4565  subfind_process_single_group()|src/subfind/subfind_processing.cc|279
  12    67 1                                     Tail      2.1299     202.5865  subfind_process_single_group()|src/subfind/subfind_processing.cc|280
  12    68 1                                      Len      1.0650     203.6515  subfind_process_single_group()|src/subfind/subfind_processing.cc|281
  12    69 1                          coll_candidates      0.2130     203.8645  subfind_process_single_group()|src/subfind/subfind_processing.cc|291
  12    70 1                             D->StartList      0.0001     203.8646  domain_allocate()|src/domain/domain.cc|183
  12    71 1                               D->EndList      0.0001     203.8647  domain_allocate()|src/domain/domain.cc|184
  12    72 1                    D->FirstTopleafOfTask      0.0001     203.8648  domain_allocate()|src/domain/domain.cc|185
  12    73 1                      D->NumTopleafOfTask      0.0001     203.8649  domain_allocate()|src/domain/domain.cc|186
  12    74 1                              D->TopNodes      0.0011     203.8660  domain_allocate()|src/domain/domain.cc|187
  12    75 1                            D->TaskOfLeaf      0.0005     203.8665  domain_allocate()|src/domain/domain.cc|188
  12    76 1                       D->ListOfTopleaves      0.0005     203.8669  domain_decomposition()|src/domain/domain.cc|135
  12    77 1                              unbind_list      1.1685     205.0354  subfind_process_single_group()|src/subfind/subfind_processing.cc|659
  12    78 0                                     dold      1.1685     206.2039  subfind_unbind()|src/subfind/subfind_unbind.cc|56
  12    79 0                                   potold      2.3369     208.5407  subfind_unbind()|src/subfind/subfind_unbind.cc|57
  12    80 1                           Tree.NodeLevel      0.0001     208.5408  force_treeallocate()|src/forcetree/forcetree.c|1444
  12    81 1                         Tree.NodeSibling      0.0005     208.5413  force_treeallocate()|src/forcetree/forcetree.c|1445
  12    82 1                           Tree.NodeIndex      0.0005     208.5418  force_treeallocate()|src/forcetree/forcetree.c|1446
  12    83 1                           Tree.Task_list      2.1603     210.7021  force_treeallocate()|src/forcetree/forcetree.c|1447
  12    84 1                           Tree.Node_list      2.1603     212.8625  force_treeallocate()|src/forcetree/forcetree.c|1448
  12    85 1                               Tree.Nodes     35.7231     248.5856  force_treeallocate()|src/forcetree/forcetree.c|1450
  12    86 1                              Tree.Points      0.0001     248.5857  force_treebuild_construct()|src/forcetree/forcetree.c|310
  12    87 1                            Tree.Nextnode      3.6446     252.2303  force_treebuild_construct()|src/forcetree/forcetree.c|311
  12    88 1                              Tree.Father      3.6441     255.8744  force_treebuild_construct()|src/forcetree/forcetree.c|312
  12    89 1                                 PartList   1019.9025    1275.7769  src/subfind/subfind_treepotential.cc|16generic_alloc_partlist_nodelist_ngblist()|src/subfind/../mpi_utils/generic_comm_h|84
  12    90 1                                  Ngblist      2.1603    1277.9373  src/subfind/subfind_treepotential.cc|16generic_alloc_partlist_nodelist_ngblist()|src/subfind/../mpi_utils/generic_comm_h|88
  12    91 1                                   DataIn     13.8409    1291.7781  src/subfind/subfind_treepotential.cc|16generic_exchange()|src/subfind/../mpi_utils/generic_comm_h|417
  12    92 1                               NodeDataIn     29.4413    1321.2195  src/subfind/subfind_treepotential.cc|16generic_exchange()|src/subfind/../mpi_utils/generic_comm_h|418
  12    93 1                                  DataOut      2.3068    1323.5263  src/subfind/subfind_treepotential.cc|16generic_exchange()|src/subfind/../mpi_utils/generic_comm_h|419
  12    94 1                                  DataGet      6.2735    1329.7998  src/subfind/subfind_treepotential.cc|16generic_multiple_phases()|src/subfind/../mpi_utils/generic_comm_h|318
  12    95 1                              NodeDataGet     19.1392    1348.9390  src/subfind/subfind_treepotential.cc|16generic_multiple_phases()|src/subfind/../mpi_utils/generic_comm_h|319
  12    96 1                               DataResult      1.0456    1349.9846  src/subfind/subfind_treepotential.cc|16generic_multiple_phases()|src/subfind/../mpi_utils/generic_comm_h|320
------------------------------------------------------------------------------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



timings.txt                                                 {#timings}
===========

The file timings.txt contains detailed performance statistics of the
gravitational tree algorithm for each timestep, which is usually (but
not always) the main sink of computational time in a simulation. A
typical output for a certain step may look like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step(*): 364, t: 0.050221, dt: 4.79423e-05, highest active timebin: 19  (lowest active: 17, highest occupied: 19)
Nf=8589934592  timebin=19  total-Nf=2087354121396
   work-load balance: 1.15054   part/sec: raw=19746.3, effective=17162.6     ia/part: avg=635.949   (593.204|42.7446)
   maximum number of nodes: 636059, filled: 0.943865
   NumForeignNodes: max=395371 avg=258155 fill=0.201039
   NumForeignPoints: max=1.51499e+06 avg=973285 fill=0.0961538  cycles=42
   avg times: <all>=208.996  <tree>=174.285  <wait>=26.7201  <fetch>=0.593317  <stack>=7.38049  (lastpm=29.0059) sec
   total interaction cost: 5.46276e+12  (imbalance=1.13168)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first line of the block generated for each step informs about the
number of the current timestep, the current simulation time, and the
system timestep itself (i.e. the time difference to the last force
computation). Also the highest active time at this step and the range
of occupied timebins is reported. Then, the number behind `Nf` is the
number of gravitational forces that are computed in this invokation of
the tree code, whereas the number behind `total-Nf` gives the total
number of such force calculations since the beginning of the
simulation.

The line starting with `work-load balance` gives the actual work-load
balance measured for the summed execution times of the tree walks. The
next number, for `part/sec`, measures the raw force speed in terms of
tree-force computations per processor per second. The first number
basically gives the speed that would be achieved for perfect work-load
balance, while the actually achieved average effective force speed
will always be lower in practice due to work-load imbalance, and this
number is given after the label `effective`. The number reported for
`ia/part` gives the average number of particle-node interactions
required to compute the force for each of the active particles, and
this should be roughly anti-proportional to the raw calculational
speed. The following two numbers in parathenses give the average
number of particle-particle and particle-node interactions that were
computed per particle.

The next line reports the maximum number of tree nodes that were used
among the MPI-ranks, while the number behind `filled` gives this
quantity normalised to the number of allocated tree nodes, and hence
indicates the degree to which the tree storage was filled. The next
two lines report how many nodes and particles were imported from other
nodes by each MPI rank, both in terms of the maximum that occurred,
and the average values. Again, the values behind 'fill' give the
maximum degree to which the buffer storage for these two components
were filled. The value behind cycles indicates the maximum number of
times an MPI process needed to call the routine that fetches foreign
nodes or particles before it could continue.

The line beginning with `avg times` reports the average execution
times of different parts in the tree calculation. For `all` the total
execution time is reported, for `tree` the time the code carries out
actual tree walks and gravity calculations, and for `wait` the time
lost because some processes finish before others and then need to wait
until everybody is done (this reflects the work-load imbalance). The
time reported for `fetch` is the time MPI ranks needed to wait for the
arrival of data requested from foreign nodes, while `stack` measures
further bookkeeping time to organize the importing of data in the
first place. If the PM-algorithm is used, a number in parenthesis
gives the execution time of the most recent PM-force calculation.

Finally, the last line reports the total cost measure the code
computes for the work done in this step, and the imbalance therein. It
is this cost measure that the code tries to balance in the domain
decomposition. This is thus the imbalance the code expects to be there
based on the domain decomposition it has done, whereas the one
reported for `work-load balance` is the one measured based on the
actual execution time.


density.txt                                                 {#density}
===========

The file `density.txt` contains detailed performance statistics of the
SPH density calculation for each timestep. This is very simular in
structure and information content to the `timings.txt` output for the
gravitational tree walks. A typical output for a certain step may look
like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step: 404, t: 0.0940046, dt: 6.98469e-05, highest active timebin: 17  (lowest active: 17, highest occupied: 19)
Nf= 16777216  highest active timebin=19  total-Nf=3875640800
   work-load balance: 1.18233   part/sec: raw=224112, effective=189550
   maximum number of nodes: 11024, filled: 0.835658
   NumForeignNodes: max=8058 avg=3646.1 fill=0.00688725
   NumForeignPoints: max=99570 avg=52765 fill=0.010639  cycles=8
   avg times: <all>=1.24922  <tree>=0.959757  <wait>=0.168166  <fetch>=0.0849678  sec

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just like for the tree-gravity, the `work-load balance` gives the
ratio between the maximum execution time for SPH density loops
relative to the average among all MPI ranks. The numbers behind
`part/sec` and `effective` give the number of SPH density computations
per second completed per particle, ignoring work-load imbalance or
including it, respectively. The numbers reported for imported nodes
and particles have the same meaning as for the gravity tree, except
that they here refer to the neighbor tree, and imported particles are
exclusively SPH particles (type 0).

Finally, the last line reports the average times spent in different
parts of the calculation. Note that `tree` refers here to walking of
the neighbor tree to find neighbours and doing all SPH calculations on
them.


hydro.txt                                                 {#hydro}
=========

There is also a log file informing in detail about the performance of
the SPH hydrodynamical force calculations. Its structure and
information content follows very closely the `density.txt` file for
the SPH density computation, we therefore refrain from inlining an
example output here.



energy.txt                                                {#egystat}
==========

In the file `energy.txt`, the code gives some statistics about the
total energy of the system. In regular intervals (specified by
`TimeBetStatistics`), the code computes the total kinetic, thermal and
potential energy of the system, and it then adds one line to this
file. Each of these lines contains 28 numbers (if `NTYPES=6` is used),
which you may process by some analysis script. The first number in the
line is the output time, followed by the total internal energy of the
system (will be 0 if no gas physics is included), the total potential
energy, and the total kinetic energy.

The next 18 numbers are the internal energy, potential energy, and
kinetic energy of the six particle types. Finally, the last six
numbers give the total mass in these components.

Note that while frequent outputs of the energy quantities allow a
check of energy conservation in Newtonian dynamics, this is more
difficult for cosmological integrations, where the Layzer-Irvine
equation is needed. (Note that the softening needs to be fixed in
comoving coordinates for it.) We remark that it is in practice not
easy to obtain a precise value of the peculiar potential energy at
high redshift (should be exactly zero for a homogeneous particle
distribution). Also, the cosmic energy integration is a differential
equation, so a test of conservation of energy in an expanding cosmos
is less straightforward that one may think.




sfr.txt                                    {#sfrlog}
=======

This is only present in simulations with cooling and star
formation. It can then be used to obtain a simple global overview of the
total star formation rate in the simulation.


