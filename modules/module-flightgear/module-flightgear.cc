/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

//#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <fstream>

#include "userelem.h"
#include "socketstream_out_elem.h"
#include "socketstreamdrive.h"
#include "bufmod.h"
#include "drive_.h"

#include "module-flightgear.h"

FGNetFDMWordSet fGNetFDMWordSet;
FGNetCtrlsWordSet fGNetCtrlsWordSet;

FieldsDescriptionFG fieldsDescriptionFGNetFDM;
FieldsDescriptionFG fieldsDescriptionFGNetCtrls;
int FGNetFDMCurVersion;
int FGNetCtrlsCurVersion;

extern "C"
int module_init(const char *module_name, void *pdm, void *php){

#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif
	/*Stream output element content type for Flight Gear*/
	StreamOutputContentTypeReader *rf1 = new FlightGearStreamOutputReader;
	if (!SetStreamOutputContentType("FlightGear", rf1)) {
		delete rf1;
		return -1;
	}

	/*file drive content type for Flight Gear*/
	FileDriveContentTypeReader *rf2 = new FlightGearFileDriveReader;
	if (!SetFileDriveContentType("FlightGear", rf2)) {
		delete rf2;
		return -1;
	}

	/*file drive caller content type for Flight Gear*/
	FileDriveCallerTypeReader *rf3 = new FlightGearFileDriveCallerTypeReader;
	if (!setFileDriveCallerType("FlightGear", rf3)) {
		delete rf3;
		return -1;
	}

  buildFieldsDescriptionFG();

	return 0;
}

void buildFieldsDescriptionFG(void){
	//FGStructuresLoaded = true;
	readFGStructuresFromFile("FG_structures/FGNetFDM_structure.txt", "FG_NET_FDM_VERSION", fieldsDescriptionFGNetFDM, FGNetFDMCurVersion);
	readFGStructuresFromFile("FG_structures/FGNetCtrls_structure.txt", "FG_NET_CTRLS_VERSION", fieldsDescriptionFGNetCtrls, FGNetCtrlsCurVersion);
}

void readFGStructuresFromFile(const char *filePath, std::string currentVersionLabel, FieldsDescriptionFG &fieldsDescriptionFG, int &currentVersion){

	std::string fieldNameDelimiter = ":";//after name field this character must occur
	std::string fieldInfoDelimiter = "= ";//before field info (type, size, offset) this character must occur
	std::string fieldInfoSeparator = ",";//field info (type, size, offset) are comma separated
	std::string openSquareBracket = "[";//used to identify fields with multiple entities (array)

	std::ifstream infile(filePath);

	/*checks if file is correctly opened*/
	if(!infile.is_open()){
	  silent_cerr("Error while loading of " << filePath << ": file not found" << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::string line; /*the file is read by lines: every line is temporary stored inside of this object*/
	std::getline(infile, line);/*gets first line of the file, where field 'version' should be stored*/

	/*checks if file is empty*/
	if(infile.eof()){
	  silent_cerr("Error while reading " << filePath << ": file is empty"<< std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*reads current version*/

	if(line.find(currentVersionLabel) == std::string::npos){/*checks if version is present and it is the first field*/
	  silent_cerr("Error while reading " << filePath << ": " << currentVersionLabel << " not found at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::size_t found = line.find(fieldInfoDelimiter);/*finds delimiter before version value*/
	if(found == std::string::npos){/*checks if delimiter is present before version value*/
	  silent_cerr("Error while reading " << filePath << ": " << "delimiter '" << fieldInfoDelimiter << "' not found after " << currentVersionLabel << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	/*acquisition of the version value*/

	std::string value = std::string();
	for(int i = 0 ; isdigit(line.at(found + fieldInfoDelimiter.size() + i)) ; i++){
	  char digit = line.at(found + fieldInfoDelimiter.size() + i);
	  std::stringstream ss;
	  ss << digit;
	  std::string tmp;
	  ss >> tmp;
	  value.append(tmp);
	}

	if(value.size() == 0){/*checks if version value is present in the file*/
	  silent_cerr("Error while reading " << filePath << ": " << " missing value of " << currentVersionLabel << " at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
	  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::stringstream tmpss1;
	tmpss1 << value;
	tmpss1 >> currentVersion; //stores the acquired value inside of FGNetFDMCurVersion/FGNetCtrlsCurVersion

	/*reads all of the other fields*/

	int fieldPosition = 1; /*used to register every field with its proper position*/
	do{
	  std::getline(infile, line);
	  if(!infile.eof()){

				/*reads field name*/

	      std::string fieldName; /*stores the name of field read from the file*/
	      found = line.find(fieldNameDelimiter); /*searches for name delimiter in order to acquire field name*/
	      if(found == std::string::npos){/*checks if delimiter is missing*/
	        silent_cerr("Error while reading " << filePath << ": " << "delimiter '" << fieldNameDelimiter << "' not found at line '" << line.substr(0,line.size()-1) << "'" << std::endl);
	        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	      }
	      fieldName = line.substr(0, found); /*collects only the part of the current line containing the field name*/
	      found = line.find(fieldInfoDelimiter);/*ex: eng_state: type, size, offset = unsigned int [4], 16, 124 finds where '= ' is located*/
	      if(found == std::string::npos){
	        silent_cerr("Error while reading " << filePath << ": " << "delimiter '" << fieldInfoDelimiter << "' not found for field '" << fieldName << "'" << std::endl);
	        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	      }
	      std::string fieldInfoRawString = line.substr(line.find(fieldInfoDelimiter) + fieldInfoDelimiter.size(), line.size()); //ex: "unsigned int [4], 16, 124"

				/*reads field type*/

	      found = fieldInfoRawString.find(openSquareBracket);
	      int fieldNumOfElem; //store how many elements that field has (ex array fields will have fieldNumOfElem > 1)
				std::string typeFullName;
	      if(found != std::string::npos){/*field is an array: '[' found inside of the current line*/
	        /*reads the number inside of []*/
	        value = std::string();
	        for(int i = 0 ; isdigit(fieldInfoRawString.at(found + openSquareBracket.size() + i)) ; i++){
	          char digit = fieldInfoRawString.at(found + openSquareBracket.size() + i);
	          std::stringstream ss;
	          ss << digit;
	          std::string tmp;
	          ss >> tmp;
	          value.append(tmp);
	        }
	        if(value.size() == 0){/*checks if [] are empty*/
	          silent_cerr("Error while reading " << filePath << ": " << " empty brackets for field '" << fieldName << "' at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
	          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	        }
	        std::stringstream tmpss2;
	        tmpss2 << value;
	        tmpss2 >> fieldNumOfElem;

	        typeFullName = fieldInfoRawString.substr(0, found-1); //ex: gets "unsigned int" from current line*/
	        found = fieldInfoRawString.find(fieldInfoSeparator); //finds the first comma
	      }else{
	        fieldNumOfElem = 1; /*the current field is a single-variable field*/
	        found = fieldInfoRawString.find(fieldInfoSeparator); //finds the first comma
					if(found != std::string::npos){
	        	typeFullName = fieldInfoRawString.substr(0, found);
					}
	      }

				/*checks if both info separators ',' are missing*/
				if(found == std::string::npos){
					silent_cerr("Error while reading " << filePath << ": " << "both info separators '" <<fieldInfoSeparator << "' missing for field '" << fieldName << "' at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/*checks if one info separator ',' are missing*/
				std::string tmp = fieldInfoRawString.substr(found+1, fieldInfoRawString.size());
				if(tmp.find(fieldInfoSeparator) == std::string::npos){
					silent_cerr("Error while reading " << filePath << ": " << "missing info separator '" <<fieldInfoSeparator << "' for field '" << fieldName << "' at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

	      std::string sizeAndOffsetRawString = fieldInfoRawString.substr(found+1, fieldInfoRawString.size()); /*ex: got " 16, 124" */
	      found = sizeAndOffsetRawString.find(fieldInfoSeparator); /*gets second comma in the current line (between size and offset)*/

	      /*reads offset*/
	      int offset;
	      value = std::string();
	     for(int i = 0 ; isdigit(sizeAndOffsetRawString.at(found+2+i)) && found+2+i < sizeAndOffsetRawString.size() ; i++){
	        char digit = sizeAndOffsetRawString.at(found+2+i);
	        std::stringstream ss;
	        ss << digit;
	        std::string tmp;
	        ss >> tmp;
	        value.append(tmp);
	      }

	      if(value.size() == 0){/*checks if offset is missing*/
	        silent_cerr("Error while reading " << filePath << ": " << " missing offset for field '" << fieldName << "' at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
	        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	      }

	      std::stringstream tmpss3;
	      tmpss3 << value;
	      tmpss3 >> offset;

	      /*once all information has been acquired, field is registered inside of fieldsDescriptionFG*/

	      for(int i = 0 ; i < fieldNumOfElem ; i++){//creates 'fieldNumOfElem' entries for the current field
	        std::string label = std::string();
	        std::string finalFieldName = std::string(fieldName);

	        if(fieldNumOfElem > 1){
	          label = std::string("_");
	          std::string num;
	          std::stringstream convert;
	          convert << i+1;
	          num = convert.str();

	          if(fieldNumOfElem < 10){
	            label.append("0");
						}

	          label.append(num);
	        }
	        finalFieldName.append(label);

	        if(typeFullName.compare("unsigned int") == 0){
	          fieldsDescriptionFG.insert(FieldsDescriptionFG::value_type(finalFieldName,{typeid(uint32_t),sizeof(uint32_t),fieldPosition,offset+sizeof(uint32_t)*i}));
					}else if(typeFullName.compare("int") == 0){
	            fieldsDescriptionFG.insert(FieldsDescriptionFG::value_type(finalFieldName,{typeid(int32_t),sizeof(int32_t),fieldPosition,offset+sizeof(int32_t)*i}));
	        }else if(typeFullName.compare("float") == 0){
	          fieldsDescriptionFG.insert(FieldsDescriptionFG::value_type(finalFieldName,{typeid(float),sizeof(float),fieldPosition,offset+sizeof(float)*i}));
	        }else if(typeFullName.compare("double") == 0){
	          fieldsDescriptionFG.insert(FieldsDescriptionFG::value_type(finalFieldName,{typeid(double),sizeof(double),fieldPosition,offset+sizeof(double)*i}));
					}else{
						silent_cerr("Error while reading " << filePath << ": " << " unknown type for field '" << fieldName << "' at line '" << line.substr(0, line.size()-1) << "'" << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}

					fieldPosition++;
	      }
	  }else if(fieldPosition == 1){
			silent_cerr("Error while reading " << filePath << ": " << " file does not contain any field description" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}while(!infile.eof());

	infile.close();
}

void printOptionsOnTextFile(const char * fileName, FieldsDescriptionFG &fieldsDescriptionFG){
	std::ofstream outFile;
	outFile.open(fileName);

	if(!outFile.is_open()){
		silent_cerr("Error while writing on file: " << fileName << ", impossible to open it" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for(int i = 0 ; i < (int)fieldsDescriptionFG.size() ; i++){
		for(FieldsDescriptionFG::iterator it = fieldsDescriptionFG.begin() ; it != fieldsDescriptionFG.end() ; it++){
			if(it->second.position == i+1){
			  outFile << it->first << std::endl;
			}
		}
	}

  outFile.close();
	}

void buildFGBufCast(std::vector<BufCast *>& data, FieldsDescriptionFG &fieldsDescriptionFGInUse){
		TypeMap_t swapmap;
		SwapMapInit(swapmap);

		//every type must be reversed for FG communication! (FG expects big-endian format data)
		if(bIsLittleEndian()){
		  for(TypeMap_t::iterator i = swapmap.begin() ; i != swapmap.end() ; i++){
				i->second = true;
			}
		}

		for (FieldsDescriptionFG::iterator it = fieldsDescriptionFGInUse.begin() ; it != fieldsDescriptionFGInUse.end() ; it++){
			data[it->second.position-1] = buildOneFGBufCast(it->second.offset, swapmap, it->second.type);
		}
	}

BufCast *buildOneFGBufCast(size_t& offset, TypeMap_t& swapmap, const std::type_info &fieldType){
		BufCast *pBC(0);

		if (fieldType == typeid(int8_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(int8_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<int8_t>(offset);
			} else {
				pBC = new TBufCast<int8_t>(offset);
			}

		} else if (fieldType == typeid(uint8_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(uint8_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<uint8_t>(offset);
			} else {
				pBC = new TBufCast<uint8_t>(offset);
			}

		} else if (fieldType == typeid(int16_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(int16_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<int16_t>(offset);
			} else {
				pBC = new TBufCast<int16_t>(offset);
			}

		} else if (fieldType == typeid(uint16_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(uint16_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<uint16_t>(offset);
			} else {
				pBC = new TBufCast<uint16_t>(offset);
			}

		} else if (fieldType == typeid(int32_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(int32_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<int32_t>(offset);
			} else {
				pBC = new TBufCast<int32_t>(offset);
			}

		} else if (fieldType == typeid(uint32_t)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(uint32_t).name());
			if (i->second) {
				pBC = new TBufCastHToN<uint32_t>(offset);
			} else {
				pBC = new TBufCast<uint32_t>(offset);
			}

		} else if (fieldType == typeid(float)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(float).name());
			if (i->second) {
				pBC = new TBufCastHToN<float>(offset);
			} else {
				pBC = new TBufCast<float>(offset);
			}

		} else if (fieldType == typeid(double)) {
			TypeMap_t::const_iterator i = swapmap.find(typeid(double).name());
			if (i->second) {
				pBC = new TBufCastHToN<double>(offset);
			} else {
				pBC = new TBufCast<double>(offset);
			}

		}

		return pBC;
	}
