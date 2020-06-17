/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2020
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

/** derived from @file inifcns.cpp
 *
 *  Implementation of GiNaC's initially known functions. */

/*
 *  GiNaC Copyright (C) 1999-2020 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * Implementation of abs() function that evaluates to 0 when the real-valued arg is 0
 *
 * NOTE: will be removed as soon as GiNaC relaxes the derivative of abs().
 * 
 * 2020 Pierangelo Masarati <pierangelo.masarati@polimi.it>
 */

#include "ginac/ginac.h"

namespace GiNaC {

/** Absolute value with 0-valued derivative in when real-valued arg is 0. */
DECLARE_FUNCTION_1P(mbdyn_abs)

//////////
// absolute value with correct derivative when the real-valued arg is 0
//////////

static ex mbdyn_abs_evalf(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return mbdyn_abs(ex_to<numeric>(arg));
	
	return mbdyn_abs(arg).hold();
}

static ex mbdyn_abs_eval(const ex & arg)
{
	if (is_exactly_a<numeric>(arg))
		return abs(ex_to<numeric>(arg));

	if (arg.info(info_flags::nonnegative))
		return arg;

	if (arg.info(info_flags::negative) || (-arg).info(info_flags::nonnegative))
		return -arg;

	if (is_ex_the_function(arg, mbdyn_abs))
		return arg;

	if (is_ex_the_function(arg, exp))
		return exp(arg.op(0).real_part());

	if (is_exactly_a<power>(arg)) {
		const ex& base = arg.op(0);
		const ex& exponent = arg.op(1);
		if (base.info(info_flags::positive) || exponent.info(info_flags::real))
			return pow(mbdyn_abs(base), exponent.real_part());
	}

	if (is_ex_the_function(arg, conjugate_function))
		return mbdyn_abs(arg.op(0));

	if (is_ex_the_function(arg, step))
		return arg;

	return mbdyn_abs(arg).hold();
}

static ex mbdyn_abs_expand(const ex & arg, unsigned options)
{
	if ((options & expand_options::expand_transcendental)
		&& is_exactly_a<mul>(arg)) {
		exvector prodseq;
		prodseq.reserve(arg.nops());
		for (const_iterator i = arg.begin(); i != arg.end(); ++i) {
			if (options & expand_options::expand_function_args)
				prodseq.push_back(mbdyn_abs(i->expand(options)));
			else
				prodseq.push_back(mbdyn_abs(*i));
		}
		return dynallocate<mul>(prodseq).setflag(status_flags::expanded);
	}

	if (options & expand_options::expand_function_args)
		return mbdyn_abs(arg.expand(options)).hold();
	else
		return mbdyn_abs(arg).hold();
}

static ex mbdyn_abs_expl_derivative(const ex & arg, const symbol & s)
{
	ex diff_arg = arg.diff(s);
	if (arg.info(info_flags::real))
		return diff_arg*(2*step(arg) - 1);
	return (diff_arg*arg.conjugate()+arg*diff_arg.conjugate())/2/mbdyn_abs(arg);
}

static void mbdyn_abs_print_latex(const ex & arg, const print_context & c)
{
	c.s << "{|"; arg.print(c); c.s << "|}";
}

static void mbdyn_abs_print_csrc_float(const ex & arg, const print_context & c)
{
	c.s << "fabs("; arg.print(c); c.s << ")";
}

static ex mbdyn_abs_conjugate(const ex & arg)
{
	return mbdyn_abs(arg).hold();
}

static ex mbdyn_abs_real_part(const ex & arg)
{
	return mbdyn_abs(arg).hold();
}

static ex mbdyn_abs_imag_part(const ex& arg)
{
	return 0;
}

static ex mbdyn_abs_power(const ex & arg, const ex & exp)
{
	if ((is_a<numeric>(exp) && ex_to<numeric>(exp).is_even()) || exp.info(info_flags::even)) {
		if (arg.info(info_flags::real) || arg.is_equal(arg.conjugate()))
			return pow(arg, exp);
		else
			return pow(arg, exp/2) * pow(arg.conjugate(), exp/2);
	} else
		return power(mbdyn_abs(arg), exp).hold();
}

bool mbdyn_abs_info(const ex & arg, unsigned inf)
{
	switch (inf) {
		case info_flags::integer:
		case info_flags::even:
		case info_flags::odd:
		case info_flags::prime:
			return arg.info(inf);
		case info_flags::nonnegint:
			return arg.info(info_flags::integer);
		case info_flags::nonnegative:
		case info_flags::real:
			return true;
		case info_flags::negative:
			return false;
		case info_flags::positive:
			return arg.info(info_flags::positive) || arg.info(info_flags::negative);
		case info_flags::has_indices: {
			if (arg.info(info_flags::has_indices))
				return true;
			else
				return false;
		}
	}
	return false;
}

REGISTER_FUNCTION(mbdyn_abs, eval_func(mbdyn_abs_eval).
                       evalf_func(mbdyn_abs_evalf).
                       expand_func(mbdyn_abs_expand).
                       expl_derivative_func(mbdyn_abs_expl_derivative).
                       info_func(mbdyn_abs_info).
                       print_func<print_latex>(mbdyn_abs_print_latex).
                       print_func<print_csrc_float>(mbdyn_abs_print_csrc_float).
                       print_func<print_csrc_double>(mbdyn_abs_print_csrc_float).
                       conjugate_func(mbdyn_abs_conjugate).
                       real_part_func(mbdyn_abs_real_part).
                       imag_part_func(mbdyn_abs_imag_part).
                       power_func(mbdyn_abs_power));

} // namespace GiNaC
