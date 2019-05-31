/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
 * Author: Andrea Zanoni  <andrea.zanoni@polimi.it>
 *
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */


#include <iostream>
#include <netcdf>
#include <exception>
#include <vector>

#define OUTPRECISION 8

static void help(char* pname)
{
	std::cerr << "usage: " << pname 
		<< " <NetCDF_file> <MBDyn_variable>" << std::endl;
	std::cerr << "    Reads the MBDyn NetCDF output file <NetCDF_file> and prints" << std::endl;
	std::cerr << "    The variable <MBDyn_variable> to stdout, in tab separated format." << std::endl << std::endl;
	std::cerr << "    1D variables are printed in a single column." << std::endl;
	std::cerr << "    2D Vec3 are printed as rows of an N x 3 matrix." << std::endl;
	std::cerr << "    3D Mat3x3 values are printed as row-wise expansions in" << std::endl;
	std::cerr << "    an N x 6 matrix" << std::endl << std::endl;
	std::cerr << "    N is the number of timesteps" << std::endl << std::endl;
	std::cerr << "    If the option \"-h\", \"-H\" or \"--help\" is given, this" << std::endl;
	std::cerr << "    message is printed." << std::endl;
}

class exWrongArgC: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Wrong number of arguments. Expecting 2.";
	}
} exwargc;

class exNcFileInvalid: public std::exception
{
	virtual const char* what() const throw()
	{
		return "The NetCDF file is not readable.";
	}
} exncfiv;

class exNcDimInvalid: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Unsupported number of dimensions.";
	}
} exncdimiv;

class exNcTypeInvalid: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Unsupported variable type.";
	}
} exnctyiv;

class exArgInvalid: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Invalid argument.";
	}
} exariv;

void printDoublerealNcVar(const netCDF::NcVar var)
{
	std::vector<netCDF::NcDim> dims = var.getDims();
	std::cout.precision(OUTPRECISION);
	switch (dims.size())
	{
		case 1:
		{
			double varIn[dims[0].getSize()];
			var.getVar(varIn);
			for (size_t ii = 0; ii < dims[0].getSize(); ii++)
			{
				std::cout << std::fixed << varIn[ii] << std::endl;
			}
		}
		break;
		case 2:
		{
			const size_t SZX = dims[0].getSize();
			const size_t SZY = dims[1].getSize();
			double varIn[SZX][SZY];

			var.getVar(varIn);
			for (size_t ii = 0; ii < SZX; ii++)
			{
				std::cout << std::fixed << varIn[ii][0];
				for (size_t jj = 1; jj < SZY; jj++)
				{
					std::cout << '\t' << std::fixed << varIn[ii][jj];
				}
				std::cout << std::fixed << std::endl;
			}

		}
		break;
		default:
		throw exncdimiv;
	}
}

void printIntegerNcVar(const netCDF::NcVar var)
{
	std::vector<netCDF::NcDim> dims = var.getDims();
	switch (dims.size())
	{
		case 1:
		{
			int varIn[dims[0].getSize()];
			var.getVar(varIn);
			for (size_t ii = 0; ii < dims[0].getSize(); ii++)
			{
				std::cout << std::fixed << varIn[ii] << std::endl;
			}

		}
		break;
		default:
		throw exncdimiv;
	}
}

void printVec3NcVar(const netCDF::NcVar var)
{
	std::vector<netCDF::NcDim> dims = var.getDims();
	std::cout.precision(OUTPRECISION);
	
	const size_t SZX = dims[0].getSize();
	const size_t SZY = dims[1].getSize();
	double varIn[SZX][SZY];

	var.getVar(varIn);
	for (size_t ii = 0; ii < SZX; ii++)
	{
		std::cout << std::fixed << varIn[ii][0];
		for (size_t jj = 1; jj < SZY; jj++)
		{
			std::cout << '\t' << std::fixed << varIn[ii][jj];
		}
		std::cout << std::fixed << std::endl;
	}
}

void printMat3x3NcVar(const netCDF::NcVar var)
{
	std::vector<netCDF::NcDim> dims = var.getDims();
	std::cout.precision(OUTPRECISION);
	
	const size_t SZX = dims[0].getSize();
	const size_t SZY = dims[1].getSize();
	const size_t SZZ = dims[2].getSize();
	double varIn[SZX][SZY][SZZ];

	var.getVar(varIn);
	for (size_t ii = 0; ii < SZX; ii++)
	{
		std::cout << std::fixed << varIn[ii][0][0] << '\t'
		<< std::fixed << varIn[ii][0][1] << '\t'
		<< std::fixed << varIn[ii][0][2] << '\t'
		<< std::fixed << varIn[ii][1][0] << '\t'
		<< std::fixed << varIn[ii][1][1] << '\t'
		<< std::fixed << varIn[ii][1][2] << '\t'
		<< std::fixed << varIn[ii][2][0] << '\t'
		<< std::fixed << varIn[ii][2][1] << '\t'
		<< std::fixed << varIn[ii][2][2]
		<< std::endl;
	}
}

void printNcVar(const netCDF::NcVar var)
{
	netCDF::NcVarAtt varTypeAtt = var.getAtt("type");
	char varType[varTypeAtt.getAttLength()];
	varTypeAtt.getValues(varType);

	// Sometimes (e.g. with doublereal) the value returned by NcVarAtt::getVaules() is
	// doesn't have the same length as given by NcVarAtt::getAttLength()
	if (std::string(varType).find("doublereal") != std::string::npos)
	{
		printDoublerealNcVar(var);
	} 
	else if (std::string(varType).find("Vec3") != std::string::npos)
	{
		printVec3NcVar(var);
	}
	else if (std::string(varType).find("Mat3x3") != std::string::npos)
	{
		printMat3x3NcVar(var);
	} 
	else if (std::string(varType).find("integer") != std::string::npos)
	{
		printIntegerNcVar(var);
	}
	else
	{
		throw exnctyiv;
	}
}

int main(int argc, char** argv)
{
	try
	{
		switch (argc)
		{
			case 1:
			{
				help(argv[0]);
				return 0;
			}
			break;
			case 2:
			{
				if (!(!strcasecmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) { 
					help(argv[0]);
					throw exariv;
					return 1;
				} else {
					help(argv[0]);
					return 0;
				}	
			}
			break;
			case 3:
			{
				netCDF::NcFile ncfile(argv[1], netCDF::NcFile::read);

				if (ncfile.isNull()) {
					throw exncfiv;
				}	

				// Retrieve the variable that we need
				netCDF::NcVar var = ncfile.getVar(argv[2]);

				// Print its contents to stdout
				printNcVar(var);

				return 0;
			}
			break;
			default:
			{
				help(argv[0]);
				throw exwargc;
			}
		}


	}
	catch (std::exception& e)
	{
		std::cerr << "An exception occurred: " 
			<< e.what() << " Aborting..." << std::endl;
		return 1;
	}
}
