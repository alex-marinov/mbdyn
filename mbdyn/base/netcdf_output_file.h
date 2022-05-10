/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * ACKNOWLEDGEMENTS:
 * Support for output with NetCDF is based on a contribution
 * by Patrick Rix <patrick.rix@online.de>
 */

/* gestore dell'output */

#ifndef NETCDF_OUTPUT_FILE_H
#define NETCDF_OUTPUT_FILE_H

#if defined(USE_NETCDF)
#include <string>
#include <netcdf>


#include "binaryoutput.h"
#include "filename.h"



class NetCDFOutput: public BinaryOutput, public FileName {
private:
	netCDF::NcFile *m_pBinFile;

	netCDF::NcDim DimTime;
	netCDF::NcDim DimV1;
	netCDF::NcDim DimV3;

	std::vector<netCDF::NcVar> BinaryVars;
	typedef std::vector<netCDF::NcDim> NcDimVec;
	std::map<MBDynOutType, netCDF::NcType> TypeMap = { 
			{MBDynOutType::OutInt, netCDF::NcType::nc_INT}, 
			{MBDynOutType::OutDouble, netCDF::NcType::nc_DOUBLE}, 
			{MBDynOutType::OutInt, netCDF::NcType::nc_CHAR}, 
		};

	netCDF::NcDim CreateDim(const std::string& name, integer size = -1)
	{
		ASSERT(m_pBinFile != 0);

		netCDF::NcDim dim;
		if (size == -1) {
			dim = m_pBinFile->addDim(name);  // .c_str is useless here
		} else {
			dim = m_pBinFile->addDim(name, size);
		}

		return dim;
	}
public:
	virtual ~NetCDFOutput(void) {
	};

	virtual void Open(const int format);

	bool IsOpen() const {
		return !m_pBinFile->isNull();
	}

	size_t
	CreateVar(const std::string& name, const MBDynOutType& type,
		const AttrValVec& attrs, const NcDimVec& dims) {
		netCDF::NcVar var;

		var = m_pBinFile->addVar(name, netCDF::NcType(TypeMap[type]), dims);
		for (AttrValVec::const_iterator i = attrs.begin(); i != attrs.end(); ++i) {
			var.putAtt(i->attr, i->val);
		}
		BinaryVars.push_back(var);

		return BinaryVars.size()-1;
	}

	template <class T>
	size_t CreateVar(const std::string& name,
		const MBUnits::Dimensions phys_dim, const std::string& description)
	{
		AttrValVec attrs(3);
		NcDimVec dims{DimTime};

		//attrs[0] = AttrVal("units", units);
		attrs[0] = AttrVal("units", MBUnits::Units[phys_dim]);
		attrs[2] = AttrVal("description", description);
	
		if (typeid(T) == typeid(integer)) {
			attrs[1] = AttrVal("type", "integer");
			return CreateVar(name, MBDynOutType::OutInt, attrs, dims); // put this one inside if because couldn't figure out how to assign type inside if while declaring it outside of if.. alternative would be to create a pointer and change CreateVar function to copy type parameter instead of passing by reference (because type cannot be delete in CreateVar since addVar passes it by reference to netcdf)...
		} else if (typeid(T) == typeid(doublereal)) {
			attrs[1] = AttrVal("type", "doublereal");
			return CreateVar(name, MBDynOutType::OutDouble, attrs, dims); // see comment above
		} else if (typeid(T) == typeid(Vec3)) {
			attrs[1] = AttrVal("type", "Vec3");
			dims.resize(2);
			dims[1] = DimV3;
			return CreateVar(name, MBDynOutType::OutDouble, attrs, dims); // see comment above
		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	size_t
	CreateVar(const std::string& name, const std::string& type) {
		AttrValVec attrs(1);
		attrs[0] = AttrVal("type", type);

		NcDimVec dims{DimV1};
		return CreateVar(name, MBDynOutType::OutChar, attrs, dims);
	};

	size_t
	CreateRotationVar(const std::string& name_prefix,
		const std::string& name_postfix,
		OrientationDescription od,
		const std::string& description) {
		NcDimVec dim;
		AttrValVec attrs;
		std::string name(name_prefix);

		switch (od) {
		case ORIENTATION_MATRIX:
			dim.resize(3);
			dim[0] = DimTime;
			dim[1] = DimV3;
			dim[2] = DimV3;

			attrs.resize(3);
			attrs[0] = AttrVal("units", "-");
			attrs[1] = AttrVal("type", "Mat3x3");
			attrs[2] = AttrVal("description",
				description + " orientation matrix "
				"(R11, R21, R31, R12, R22, R32, R13, R23, R33)");

			name += "R";
			break;

		case ORIENTATION_VECTOR:
		case UNKNOWN_ORIENTATION_DESCRIPTION:  // only relevant to displacement nodes
			dim.resize(2);
			dim[0] = DimTime;
			dim[1] = DimV3;

			attrs.resize(3);
			attrs[0] = AttrVal("units", "radian");
			attrs[1] = AttrVal("type", "Vec3");
			attrs[2] = AttrVal("description",
				description + " orientation vector "
				"(Phi_X, Phi_Y, Phi_Z)");

			name += "Phi";
			break;

		case EULER_123:
		case EULER_313:
		case EULER_321:
		{
			dim.resize(2);
			dim[0] = DimTime;
			dim[1] = DimV3;

			attrs.resize(3);
			attrs[0] = AttrVal("units", "deg");
			attrs[1] = AttrVal("type", "Vec3");

			std::string etype;
			switch (od) {
			case EULER_123:
				etype = "123";
				break;

			case EULER_313:
				etype = "313";
				break;

			case EULER_321:
				etype = "321";
				break;

			default:
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			attrs[2] = AttrVal("description",
				description + " orientation Euler angles (" + etype + ") "
				"(E_X, E_Y, E_Z)");

			name += "E";
			} break;

		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		name += name_postfix;
		return CreateVar(name, MBDynOutType::OutDouble, attrs, dim);
	};
	
	virtual void WriteNcVar(size_t Var_Var, const Mat3x3& pGetVar) {
		BinaryVars[Var_Var].putVar(Start1x3x3, Count1x3x3, pGetVar.pGetMat());
	}
	virtual void WriteNcVar(size_t Var_Var, const Mat3x3& pGetVar,
		const size_t& ncStart) 
	{
		std::vector<size_t> Start1x3x3Tmp = Start1x3x3;
		Start1x3x3Tmp[0] = ncStart;
		BinaryVars[Var_Var].putVar(Start1x3x3Tmp, Count1x3x3, pGetVar.pGetMat());
	}
	virtual void WriteNcVar(size_t Var_Var, const Vec3& pGetVar) {
		BinaryVars[Var_Var].putVar(Start1x3, Count1x3, pGetVar.pGetVec());
	}
	virtual void WriteNcVar(size_t& Var_Var, const Vec3& pGetVar,
		const size_t& ncStart) 
	{
		std::vector<size_t> Start1x3Tmp = Start1x3;
		Start1x3Tmp[0] = ncStart;
		BinaryVars[Var_Var].putVar(Start1x3Tmp, Count1x3, pGetVar.pGetVec());
	}

};

#else /* !defined USE_NETCDF */
class NetCDFOutput {
#error
};
#endif /* !defined USE_NETCDF */

#endif /* NETCDF_OUTPUT_FILE_H */
