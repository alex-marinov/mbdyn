/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
        Copyright (C) 2011(-2014) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef OCTAVE_OBJECT_H_
#define OCTAVE_OBJECT_H_

#include <mbconfig.h>

#ifdef USE_OCTAVE

#include <octave/oct.h>
#include <cassert>
#include <typeinfo>
#include <stdarg.h>
#include <stdexcept>
#include <string>
#include <map>

namespace oct {

class octave_object: public octave_base_value
{
protected:
	typedef octave_value_list method_function(octave_object*, const octave_value_list&, int);
	typedef void setup_class_object_function();
	class class_object
	{
	public:
		explicit class_object(class_object* parent_obj,setup_class_object_function* setup_class_obj_func)
			:parent_object(parent_obj)
		{
			if ( NULL != setup_class_obj_func )
				(*setup_class_obj_func)();
		}

		method_function*& operator[](const std::string& method_name)
		{
			return method_table[method_name];
		}

		size_t size()const
		{
			return method_table.size();
		}
		method_function* lookup_method(const std::string& method_name)const;
	private:
		typedef std::map<std::string,method_function*> method_table_t;
		method_table_t method_table;
		class_object* const parent_object;
	};
protected:
	octave_object();
	virtual ~octave_object();
protected:
	static class_object dispatch_class_object;
private:
	virtual octave_value operator()(const octave_value_list& idx) const;
	virtual const class_object* get_class_object()=0;
	virtual octave_value_list subsref(const std::string& type,
			const std::list<octave_value_list>& idx,
			int nargout);
	virtual octave_value subsref(const std::string& type,
			const std::list<octave_value_list>& idx);
	method_function* lookup_method(const std::string& method_name);

	virtual bool is_constant(void)const{ return true; }
	virtual bool is_defined(void)const{ return true; }
};

template <typename T>
class octave_object_ptr: public octave_value
{
public:
	octave_object_ptr()
	{

	}

	explicit octave_object_ptr(T* p)
		:octave_value(p)
	{

	}

	explicit octave_object_ptr(octave_base_value* p)
		:octave_value(p)
	{
		check_type(*p);
	}

	octave_object_ptr(const octave_object_ptr& ptr)
		:octave_value(ptr)
	{
	}

	octave_object_ptr(const octave_value& rhs)
		:octave_value(rhs)
	{
		check_type(rhs.get_rep());
	}

	template <typename U>
	octave_object_ptr(const octave_object_ptr<U>& rhs)
		:octave_value(rhs)
	{
		check_type(rhs.get_rep());
	}

	octave_object_ptr& operator=(const octave_object_ptr& rhs)
	{
		octave_value::operator =(rhs);
		return *this;
	}

	octave_object_ptr& operator=(const octave_value& rhs)
	{
		check_type(rhs.get_rep());

		octave_value::operator=(rhs);
		return *this;
	}

	template <typename U>
	octave_object_ptr& operator=(const octave_object_ptr<U>& rhs)
	{
		check_type(rhs.get_rep());

		octave_value::operator=(rhs);
		return *this;
	}

	T& get_rep()const
	{
		return static_cast<T&>(const_cast<octave_base_value&>(octave_value::get_rep()));
	}

	T* operator->()const
	{
		return &get_rep();
	}

private:
	static void check_type(const octave_base_value& ref) throw(std::bad_cast)
	{
		if ( NULL == dynamic_cast<const T*>(&ref) )
			throw std::bad_cast();
	}
};

extern void error(const char* fmt, ...)
#ifdef __GNUC__
	__attribute__((format (printf, 1, 2)))
#else
	#warning "printf format will not be checked!"
#endif
;

#define BEGIN_METHOD_TABLE_DECLARE() \
	virtual const class_object* get_class_object(); \
	static class_object dispatch_class_object; \
	static void setup_class_object();

#define METHOD_DECLARE(method_name) \
   octave_value_list method_name(const octave_value_list&, int); \
   static octave_value_list method_name(octave_object* object__,const octave_value_list&,int);

#define END_METHOD_TABLE_DECLARE()

#define METHOD_DEFINE(class_name,method_name,args,nargout) \
		octave_value_list class_name::method_name(octave_object *object__, \
          	  	  	  	  	  	  	  	  		  const octave_value_list& args, \
          	  	  	  	  	  	  	  	  		  int nargout) \
        { \
			return dynamic_cast<class_name&>(*object__).method_name(args,nargout); \
		} \
		octave_value_list class_name::method_name(const octave_value_list& args, int nargout)

#define BEGIN_METHOD_TABLE(class_name,base_class_name) \
   octave_object::class_object class_name::dispatch_class_object(&base_class_name::dispatch_class_object,&class_name::setup_class_object); \
   const octave_object::class_object* class_name::get_class_object() \
   { \
	   return &dispatch_class_object; \
   } \
   void class_name::setup_class_object() \
   { \
	   if ( dispatch_class_object.size() == 0 ) \
	   {

#define METHOD_DISPATCH(class_name,method_name) \
			dispatch_class_object[#method_name] = &class_name::method_name;

#define END_METHOD_TABLE() \
	   } else assert(false); \
	}

} // namespace

#endif // USE_OCTAVE

#endif /* OCTAVE_OBJECT_H_ */
