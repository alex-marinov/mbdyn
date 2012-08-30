/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2012) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <mbconfig.h>

#ifdef USE_OCTAVE

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/dynamic-ld.h>
#include <octave/oct-map.h>
#include <octave/oct-stream.h>
#include <octave/ov-base-scalar.h>

#include <iostream>

using namespace std;

#include "octave_object.h"

namespace oct {

octave_object::method_function* octave_object::class_object::lookup_method(const std::string& method_name)const
{
	typedef method_table_t::const_iterator iter_t;

	iter_t p = method_table.find(method_name);

	if ( p != method_table.end() )
		return p->second;

	if ( parent_object != NULL )
	{
		return parent_object->lookup_method(method_name);
	}
	else
	{
		cerr << "octave_object::class_object::lookup_method(\"" << method_name << "\") failed!" << endl;
		const class_object* class_obj = this;

		cerr << "dump of octave_object::class_object::method_table:" << endl;

		while ( class_obj != NULL )
		{
			for ( iter_t p = class_obj->method_table.begin(); p != class_obj->method_table.end(); ++p )
				cerr << " " << p->first << "->" << p->second <<  endl;
			class_obj = class_obj->parent_object;
		}
	}

	return NULL;
}

octave_object::class_object octave_object::dispatch_class_object(NULL,NULL);

octave_object::octave_object()
{

}

octave_object::~octave_object()
{

}

static bool any_arg_is_magic_colon (const octave_value_list& args)
{
  int nargin = args.length ();

  for (int i = 0; i < nargin; i++)
    if (args(i).is_magic_colon ())
      return true;

  return false;
}

#if TRACE_SUBSREF == 1

static void print_args(const octave_value_list& args)
{
  int nargin = args.length ();

  for (int i = 0; i < nargin; i++)
  {
    args(i).print(octave_stdout);

    if( i < nargin - 1 )
    	octave_stdout << ',';
  }
}

#endif

octave_object::method_function* octave_object::lookup_method(const std::string& method_name)
{
	return get_class_object()->lookup_method(method_name);
}

octave_value_list octave_object::subsref (const std::string& type,
		const std::list<octave_value_list>& idx,
		int nargout)
{
	#if TRACE_SUBSREF == 1
	{
		octave_stdout << endl;
		octave_stdout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << endl;
		octave_stdout << "type=" << type << endl;

		typedef std::list<octave_value_list>::const_iterator iter_t;

		int i = 0;

		for ( iter_t it = idx.begin(); it != idx.end(); ++it, ++i )
		{
			octave_stdout << "idx(" << i << ")=";
			print_args(*it);
			octave_stdout << endl;
		}
	}
	#endif

	octave_value_list retval;
	size_t skip = 1;

	if ( type.length() < 1 )
	{
		panic_impossible();
	}

	const char operator_char = type[0];

	switch( operator_char )
	{
	case '.':
	{
		octave_value_list method_name_list = idx.front();

		if ( method_name_list.length () != 1 || !method_name_list(0).is_string() )
			panic_impossible();

		const std::string method_name = method_name_list(0).string_value();

		method_function *method_pfn = lookup_method(method_name);

		if ( method_pfn == NULL )
		{
			error("octave_object: class \"%s\" has no member \"%s\"", class_name().c_str(), method_name.c_str());
			return retval;
		}

		if ( idx.size() < 2 || type[1] != '(')
		{
			retval = (*method_pfn)(this,octave_value_list(),nargout);
		}
		else
		{
			skip = 2;
			std::list<octave_value_list>::const_iterator pidx = idx.begin();
			pidx++;

			if (any_arg_is_magic_colon (*pidx))
			{
				error("octave object: invalid use of colon in method argument list");
				return retval;
			}
			else
			{
				retval = (*method_pfn)(this,*pidx,nargout);
			}
		}
	}
	break;
	case '(':
	{
		if (idx.size() < 1)
		{
			error("octave object: invalid number of indices");
			return retval;
		}

		retval = (*this)(idx.front());
	}
	break;
	case '{':
	{
		std::string nm = type_name();
		error("%s cannot be indexed with %c", nm.c_str(), operator_char);
		return retval;

	}
	break;

	default:
		panic_impossible();
	}

	if ( !error_state && idx.size() > skip )
	{
		retval = retval(0).next_subsref(type, idx, skip);
	}

	return retval;
}

octave_value octave_object::operator()(const octave_value_list& idx) const
{
	std::string nm = type_name();
	error("%s cannot be indexed with %c", nm.c_str(), '(');
	return octave_value();
}

octave_value octave_object::subsref(const std::string& type,
	  	  	  	  	  	  	 const std::list<octave_value_list>& idx)
{
	octave_value_list ret_val = subsref(type,idx,1);

	if ( ret_val.length() < 1 )
		return octave_value();

	return ret_val(0);
}

void error(const char* fmt, ...)
{
	va_list va;
	va_start(va, fmt);
	::verror(fmt, va);
	va_end(va);
}

} // namespace

#endif // USE_OCTAVE
