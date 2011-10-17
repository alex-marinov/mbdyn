/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * Copyright 1999-2011 Giuseppe Quaranta <quaranta@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * This copyright statement applies to the MPI related code, which was
 * merged from files schur.h/schur.cc
 */

/*
 *
 * Copyright (C) 2003-2011
 * Giuseppe Quaranta	<quaranta@aero.polimi.it>
 *
 */

/* metodo per la soluzione del modello */

#ifndef SOLVER_IMPL_H
#define SOLVER_IMPL_H

#include <limits>

#if defined(HAVE_SIGNAL) && defined(HAVE_SIGNAL_H)
#include <signal.h>
#endif /* HAVE_SIGNAL && HAVE_SIGNAL_H */

#ifdef USE_MPI
#include "mbcomm.h"
#ifdef USE_EXTERNAL
#include "external.h"
#endif /* USE_EXTERNAL */
#endif /* USE_MPI */


#if defined(HAVE_SYS_MMAN_H)
#include <sys/mman.h>
#endif /* HAVE_SYS_MMAN_H */

#if 0
#ifdef HAVE_SIGNAL
extern volatile sig_atomic_t mbdyn_keep_going;
extern __sighandler_t mbdyn_sh_term;
extern __sighandler_t mbdyn_sh_int;
extern __sighandler_t mbdyn_sh_hup;
extern __sighandler_t mbdyn_sh_pipe;

extern "C" void mbdyn_really_exit_handler(int signum);
extern "C" void mbdyn_modify_last_iteration_handler(int signum);
extern "C" void mbdyn_modify_final_time_handler(int signum);
#endif /* HAVE_SIGNAL */
#endif

extern "C" void mbdyn_signal_init(int pre);

extern int mbdyn_reserve_stack(unsigned long size);

/* Parametri locali */
static const doublereal dDefaultDerivativesCoefficient = 1.e-6;
static const integer iDefaultDummyStepsNumber = 0;
static const doublereal dDefaultDummyStepsRatio = 1.e-3;
static const integer iDefaultIterationsBeforeAssembly = 2;
static const integer iDefaultIterativeSolversMaxSteps = 100;
static const integer iDefaultPreconditionerSteps = 20;
static const doublereal dDefaultTol = 1.e-6;
static const doublereal defaultIterativeEtaMax = 0.9;
static const doublereal defaultIterativeTau = 1.e-7;

static const integer iDefaultMaxIterations = 1;
static const doublereal dDefaultMinTimeStep = -1.;
static const doublereal dDefaultMaxTimeStep = std::numeric_limits<doublereal>::max();
static const doublereal dDefaultDummyStepsTolerance = dDefaultTol;

#endif /* ! SOLVER_IMPL_H */
