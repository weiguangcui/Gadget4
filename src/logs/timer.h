/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file timer.h
 *
 *  \brief defines the CPU consumption timers
 */

/*!
 *  All CPU timers are defined here. Additionally the macros to use these timers are defined.
 */

#if defined(TIMER_ENUM) || defined(TIMER_STRUCT)

/*! \def TIMER_CREATE(name, desc, par, symba, symbb)
 * \brief creates a new CPU timer
 *
 * \param name name used in the code to reference this timer
 * \param desc description string used in output files
 * \param parent parent of this timer to build a tree-like hierarchy of timers
 * \param symba character used for active time in balance.txt
 * \param symbb character used for imbalance in balance.txt
 *
 */

#if defined(TIMER_ENUM)

#define TIMER_CREATE(name, desc, parent, symba, symbb) name,

#else

#define TIMER_CREATE(name, desc, par, symba, symbb) \
  Timer_data[name].parent = par;                    \
  strcpy(Timer_data[name].shortname, #name);        \
  strcpy(Timer_data[name].longname, (desc));        \
  Timer_data[name].symb      = (symba);             \
  Timer_data[name].symbImbal = (symbb);
#endif

/*add your counter here, they must appear in the right order*/

TIMER_CREATE(CPU_ALL, "total", CPU_ROOT, 0, 0) /*!< root timer, everything should be below this timer */
TIMER_CREATE(CPU_TREE, "treegrav", CPU_ALL, 0, 0)
TIMER_CREATE(CPU_TREEBUILD, "treebuild", CPU_TREE, 0, 0)
TIMER_CREATE(CPU_TREEBUILD_INSERT, "insert", CPU_TREEBUILD, ':', 'r')
TIMER_CREATE(CPU_TREEBUILD_BRANCHES, "branches", CPU_TREEBUILD, 'D', 'd')
TIMER_CREATE(CPU_TREEBUILD_TOPLEVEL, "toplevel", CPU_TREEBUILD, 'E', 'e')
TIMER_CREATE(CPU_TREEFORCE, "treeforce", CPU_TREE, 0, 0)
TIMER_CREATE(CPU_TREEWALK, "treewalk", CPU_TREEFORCE, '*', 't')
TIMER_CREATE(CPU_TREEIMBALANCE, "treeimbalance", CPU_TREEFORCE, 'x', 'X')
TIMER_CREATE(CPU_TREEFETCH, "treefetch", CPU_TREEFORCE, '&', 'e')
TIMER_CREATE(CPU_TREESTACK, "treestack", CPU_TREEFORCE, '#', 'a')
#ifdef ALLOW_DIRECT_SUMMATION
TIMER_CREATE(CPU_TREEDIRECT, "treedirect", CPU_TREE, 'r', '2')
#endif
#ifdef PMGRID
TIMER_CREATE(CPU_PM_GRAVITY, "pm_grav", CPU_ALL, '|', 'n')
#endif
TIMER_CREATE(CPU_NGBTREEUPDATEVEL, "ngbtreevelupdate", CPU_ALL, 'B', 'b')
TIMER_CREATE(CPU_NGBTREEUPDATEMAXHSML, "ngbtreehsmlupdate", CPU_ALL, 'h', 'H')
TIMER_CREATE(CPU_SPH, "sph", CPU_ALL, 0, 0)
TIMER_CREATE(CPU_DENSITY, "density", CPU_SPH, 'd', 'D')
TIMER_CREATE(CPU_DENSWALK, "densitywalk", CPU_DENSITY, 'e', 'E')
TIMER_CREATE(CPU_DENSFETCH, "densityfetch", CPU_DENSITY, 'f', 'F')
TIMER_CREATE(CPU_DENSIMBALANCE, "densimbalance", CPU_DENSITY, 'c', 'C')
TIMER_CREATE(CPU_HYDRO, "hydro", CPU_SPH, 'J', 'j')
TIMER_CREATE(CPU_HYDROWALK, "hydrowalk", CPU_HYDRO, '9', 'H')
TIMER_CREATE(CPU_HYDROFETCH, "hydrofetch", CPU_HYDRO, '0', 'K')
TIMER_CREATE(CPU_HYDROIMBALANCE, "hydroimbalance", CPU_HYDRO, '[', 'y')
TIMER_CREATE(CPU_DOMAIN, "domain", CPU_ALL, '+', 'l')
TIMER_CREATE(CPU_PEANO, "peano", CPU_ALL, '"', 'o')
TIMER_CREATE(CPU_DRIFTS, "drift/kicks", CPU_ALL, '?', 'A')
TIMER_CREATE(CPU_TIMELINE, "timeline", CPU_ALL, ')', 'I')
TIMER_CREATE(CPU_TREE_TIMESTEPS, "treetimesteps", CPU_ALL, '0', 'W')
TIMER_CREATE(CPU_SNAPSHOT, "i/o", CPU_ALL, '3', 'B')
TIMER_CREATE(CPU_LOGS, "logs", CPU_ALL, 'K', 'k')
#if defined(COOLING) || defined(STARFORMATION)
TIMER_CREATE(CPU_COOLING_SFR, "sfrcool", CPU_ALL, '1', 'T')
#endif
#ifdef FOF
TIMER_CREATE(CPU_FOF, "fof", CPU_ALL, '5', 'D')
TIMER_CREATE(CPU_FOFWALK, "fofwalk", CPU_FOF, '6', 'G')
TIMER_CREATE(CPU_FOFIMBAL, "fofimbal", CPU_FOF, '7', 'H')
#endif
#ifdef SUBFIND
TIMER_CREATE(CPU_SUBFIND, "subfind", CPU_ALL, '6', 'E')
#endif
#ifdef NGENIC
TIMER_CREATE(CPU_NGENIC, "ngenic", CPU_ALL, 'n', 'N')
#endif
TIMER_CREATE(CPU_RESTART, "restart", CPU_ALL, 'C', 'c')
#ifdef FORCETEST
TIMER_CREATE(CPU_FORCETEST, "forcetest", CPU_ALL, 't', 'T')
#endif
TIMER_CREATE(CPU_MISC, "misc", CPU_ALL, '7', 'F')
TIMER_CREATE(CPU_LAST, "LAST", CPU_NONE, ' ', ' ') /*!<last item, do not use! */

#undef TIMER_CREATE

#endif
