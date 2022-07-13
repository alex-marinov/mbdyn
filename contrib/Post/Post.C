/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
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
 * along with this program; If not, see <https://www.gnu.org/licenses/>.
 *
 */
/* Copyright Marco Morandini */

#include<cstdlib>
#include<fstream>
// #include<iostream>
// #include<sstream>
// #include<iosfwd>
#include<string>
#include<sstream>
#include<map>
#include<list>
#include<vector>
#include<limits>
#include<algorithm>

#undef yyFlexLexer
#define yyFlexLexer MbdynpostFlexLexer
#include <FlexLexer.h>
#undef yyFlexLexer

std::string output_frequency_string;
int output_frequency = 0;

int step = 0;
double current_time = 0.0;
double tstep = 0.;
int niter = 0;
double reserr = 0.;
double solerr = 0.;
int solconv = 0;
bool out_flag;


#include "Post.hh"

void skipws(std::istream * file) {
	typedef std::streambuf::traits_type traits_type;
	const int eof = traits_type::eof();
	int c;
	bool testdelim;// = (c == ' ' || c == '\t');
	c = file->get();
	while (traits_type::eq_int_type(c, ' ') 
		|| traits_type::eq_int_type(c, '\t') ) {
		c = file->get();
	}
	testdelim = ( traits_type::eq_int_type(c, ' ') 
		|| traits_type::eq_int_type(c, '\t') );
	while (!testdelim) {
		c = file->get();
		testdelim = ( traits_type::eq_int_type(c, ' ' )
			|| traits_type::eq_int_type(c, '\t' )
			|| traits_type::eq_int_type(c, '\n' )
			|| traits_type::eq_int_type(c, eof ) 
		);
	}
	if (c == '\n') {
		file->unget();
		//buf->sputback((char)c);
	}


// 	std::streambuf * buf = file->rdbuf();
// 	typedef std::streambuf::traits_type traits_type;
// 	const int eof = traits_type::eof();
// 	int c = buf->sgetc();
// 	bool testdelim;// = (c == ' ' || c == '\t');
// 	while (traits_type::eq_int_type(c, ' ') 
// 		|| traits_type::eq_int_type(c, '\t') ) {
// 		c = buf->sbumpc();
// 	}
// 	testdelim = ( traits_type::eq_int_type(c, ' ') 
// 		|| traits_type::eq_int_type(c, '\t') );
// 	while (!testdelim) {
// 		c = buf->sbumpc();
// 		testdelim = ( traits_type::eq_int_type(c, ' ' )
// 			|| traits_type::eq_int_type(c, '\t' )
// 			|| traits_type::eq_int_type(c, '\n' )
// 			|| traits_type::eq_int_type(c, eof ) 
// 		);
// 	}
// 	if (c == '\n') {
// 		buf->sungetc();
// 		//buf->sputback((char)c);
// 	}
// // 	while (testdelim) {
// // 		c = buf->sbumpc();
// // 		testdelim = (c == ws || c == tab);
// // 	}
// // 	//buf->sputback(c);
// // 	buf->sungetc();
}

void skipline(std::istream * file) {
	std::streambuf * buf = file->rdbuf();
// 	const int nl = '\n';
// 	const int eof = std::streambuf::traits_type::eof();
// 	int c = buf->sbumpc();
// 	bool testdelim = (c == '\n' || c == eof);
// 	while (!testdelim) {
// 		c = buf->sbumpc();
// 		testdelim = (c == '\n' || c == eof);
// 	}

	std::streambuf::int_type eof = std::streambuf::traits_type::eof();
	typedef std::streambuf::traits_type traits_type;
	std::streambuf::int_type c = buf->sgetc();
	while (!traits_type::eq_int_type(c, eof)
			 && !traits_type::eq_int_type(c, '\n')) {
		c = buf->snextc();
	}
	buf->sbumpc();
}


struct OutputElement {
	int label;
	int num_outputs;
	std::map<int,int> output_cols;
	std::vector<double> outputs;
	OutputElement(const int i) : label(i), num_outputs(0) {};
	bool operator == (const int i) {return i == label;};
};

struct OutputFile {
	std::string name;
	int num_elements;
	std::istream * file;
	std::list<OutputElement> elements;	//ordinati per posizione in output
	std::map<int, std::list<OutputElement>::iterator> element_positions; //<position in the file, position in output> 
						//only for elements stored in std::list<OutputElement> elements
	OutputFile(std::string s) : name(s), num_elements(0), file(0) {};
};

std::list<OutputFile> files;


bool ParseLog(std::istream& in, std::ostream& out) {
	MbdynpostFlexLexer flex(&in, &out);
	//flex.set_debug(1);
	int type;
	while ((type = flex.yylex()) != EOF_TOK) {
		//std::cerr << "tipo token " << type << " " << flex.YYText() <<
		//std::endl << "-----------------------------" << std::endl;
		switch (type) {
			case OUTPUT_FREQUENCY_TOK:
				in >> output_frequency_string;
				if (output_frequency_string != std::string("custom")) {
					output_frequency = std::stoi(output_frequency_string);
				}
					
				return true;
				break;
			default:
				break;
		}
	}
	return false;
}
bool ParseOut(MbdynpostFlexLexer &flex, std::istream& in) {
	int type;
	while ((type = flex.yylex()) != EOF_TOK) {
		switch (type) {
			case STEP_TOK:
				in >> step >> current_time >> tstep >> niter >> reserr >> solerr >> solconv >> out_flag;
				if (out_flag) {
					return true;
				}
				break;
			default:
				break;
		}
	}
	return false;
}

bool ParseCommands(std::istream& in, std::ostream& out) {
	MbdynpostFlexLexer flex(&in, &out);
// 	std::string a; in >> a;
// 	std::cerr << a << " xxxx\n";
	//flex.set_debug(1);
	int type;
	do {
		if ((type = flex.yylex()) != FILE_EXTENSION_TOK) {
			std::cerr << "Unrecognized post_description token \""
				<< flex.YYText() << "\"\n"
				<< "Expected a file extension\n";
			return false;
		}
		files.push_back(OutputFile(flex.YYText()));
		if ((type = flex.yylex()) != DUEPUNTI_TOK) {
			std::cerr << "Unrecognized post_description token \""
				<< flex.YYText() << "\"\n"
				<< "Expected a \":\"\n";
			return false;
		}
		do {
			if ((type = flex.yylex()) != LABEL_TOK) {
				std::cerr << "Unrecognized post_description token \""
					<< flex.YYText() << "\"\n"
					<< "Expected an element label\n";
				return false;
			}
			OutputFile &fl(*files.rbegin());
			fl.num_elements++;
			int i = atoi(flex.YYText());
			auto el_it = std::find(fl.elements.begin(), fl.elements.end(), i);
			if (el_it == fl.elements.end()){
				fl.elements.push_back(OutputElement(i));
				el_it = fl.elements.end();
				el_it--;
			} else {
				std::cerr << "Error, request for the ouptu of element " << i << std::endl;
				std::cerr << "Some output for same element has already been requested." << std::endl;
				std::cerr << "you should ask for an element only once," << std::endl;
				std::cerr << "separating the different requested columns with commas" << std::endl;
				std::cerr << std::endl;
				std::cerr << "See Post --help for details" << std::endl;
				exit(1);
			}
			OutputElement &el(*el_it);
			if ((type = flex.yylex()) != DUEPUNTI_TOK) {
				std::cerr << "Unrecognized post_description token \""
					<< flex.YYText() << "\"\n"
					<< "Expected a \":\"\n";
				return false;
			}
			do {
				if ((type = flex.yylex()) != LABEL_TOK) {
					std::cerr << "Unrecognized post_description token \""
						<< flex.YYText() << "\"\n"
						<< "Expected an integer\n";
					return false;
				}
				int i = atoi(flex.YYText());
				el.output_cols[i] = el.num_outputs;
				el.num_outputs++;
			} while ((type = flex.yylex()) == COMMA_TOK);
			el.outputs.resize(el.num_outputs);
		} while (type == DUEPUNTI_TOK);
	} while (type == DUEDUEPUNTI_TOK);
	if ((type = flex.yylex()) != EOF_TOK) {
		std::cerr << "Unrecognized post_description token \""
			<< flex.YYText() << "\"\n"
			<< "Expected the end of commands\n";
			return false;
	}
		
// 		std::cerr << "tipo token " << type << " " << flex.YYText() <<
// 		std::endl << "-----------------------------" << std::endl;
// 		switch (type) {
// 			case FILE_EXTENSION_TOK:
// 				std::cerr << "Estensione " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			case LABEL_TOK:
// 				std::cerr << "Label " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			case COLUMN_TOK:
// 				std::cerr << "Colonna " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			case COMMA_TOK:
// 				std::cerr << "Comma " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			case DUEPUNTI_TOK:
// 				std::cerr << "Duepunti " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			case DUEDUEPUNTI_TOK:
// 				std::cerr << "Dueduepunti " << flex.YYText() << 
// 					" = " << flex.YYLeng() << "\n";
// 				break;
// 			default:
// 				std::cerr << "Unrecognized post_description token \""
// 					<< flex.YYText() << "\"\n";
// 				return false;
// 				break;
// 		}
// 	}
	return true;
}




#include <argp.h>
const char *argp_program_version = "MbdynPost-0.1";
const char *argp_program_bug_address = "<marco.morandini@polimi.it>";
static char args_doc[] = "post_description input_file_basename";
static char doc[] = "MbdynPost -- another Mbdyn result parser\n"
"\n"
" MBDyn (C) is a multibody analysis code. \n"
" http://www.mbdyn.org\n"
"\n"
" Copyright (C) 1996-2017\n"
"\n"
" Pierangelo Masarati	<pierangelo.masarati@polimi.it>\n"
" Paolo Mantegazza	<paolo.mantegazza@polimi.it>\n"
"\n"
" Changing this copyright notice is forbidden.\n"
"\n"
" This program is free software; you can redistribute it and/or modify\n"
" it under the terms of the GNU General Public License as published by\n"
" the Free Software Foundation (version 2 of the License).\n"
"\n"
" This program is distributed in the hope that it will be useful,\n"
" but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" GNU General Public License for more details.\n"
"\n"
" You should have received a copy of the GNU General Public License\n"
" along with this program; if not, write to the Free Software\n"
" Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n"
"\n"
"\nCopyright 2017 Marco Morandini\n"
"\n"
"\nIt is assumed that the file input_file_basename.log exists.\n"
"\nThe post_description argument has the following syntax:\n"
"\tpost_description: file_output[::file_output]*\n"
"\tfile output:      file_suffix:el_output[:el_output]*\n"
"\tel_output:        el_number:col_number[,col_number]*\n\n"
;
static struct argp_option options[] = {
       {"output", 'o', "FILEOUT", 0, "Output to FILEOUT.dat instead of standard output" },
       {"skip",   's', "SKIP", 0, "Output every SKIP time steps" },
       {"begin",  'b', "t0", 0, "Output from t = t0" },
       {"end",    'e', "t1", 0, "Output up to t = t1" },
       { 0 }
     };
     
struct arguments {
       char *args[2];                /* input file basename */
       char *output_file;
       int skip;
       float t0;
       float t1;
       arguments() : output_file(0), skip(0), t0(0.), t1(0.) {};
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
	/* Get the INPUT argument from `argp_parse', which we
	   know is a pointer to our arguments structure. */
	struct arguments *argu = (arguments*)state->input;

	switch (key) {
		case 'o':
			argu->output_file = arg;
		break;

		case 's':
			argu->skip = atoi(arg);
		break;

		case 'b':
			argu->t0 = atof(arg);
		break;

		case 'e':
			argu->t1 = atof(arg);
		break;

		case ARGP_KEY_ARG:
			if (state->arg_num >= 2) {
				/* Too many arguments. */
				std::cerr << "too many\n";
				argp_usage (state);
			}

			argu->args[state->arg_num] = arg;
			break;

		case ARGP_KEY_END:
			if (state->arg_num < 2) {
				/* Not enough arguments. */
				std::cerr << "too less\n";
				argp_usage (state);
			}
		break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp_s = { options, parse_opt, args_doc, doc };

struct allocate_file {
	std::string basename;
	bool success;
	allocate_file(std::string bs) : basename(bs), success(true) {};
	void operator()(OutputFile& it) {
		std::string namefile(basename);
		namefile += ".";
		namefile += it.name;
		it.name = namefile;
		std::ifstream * file;
		file = new std::ifstream(namefile.c_str());
		if (file == 0 || !file->is_open()) {
			file = 0;
			std::cerr << "Unable to open result file " << namefile << "\n";
		} else {
			int label, first_label, num_labels;
			*file >> label;
			first_label = label;
			num_labels = 0;
			std::map<int, int> element_labels;
			if (!file->eof()) {
				do {
					element_labels[label] = num_labels;
					num_labels++;
					std::string line;
					std::getline(*file, line);
					*file >> label;
				} while (!file->eof() && (label != first_label));
			}
			it.num_elements = num_labels;
			file->seekg(0, std::ios::beg);
			for (std::list<OutputElement>::iterator ite = it.elements.begin();
				ite != it.elements.end(); ite++) {
				std::map<int, int>::const_iterator itep;
				if ((itep = element_labels.find(ite->label)) != element_labels.end()) {
					it.element_positions[itep->second] = ite;
				} else {
					std::cerr << "Requested output for inexistent element " << namefile
						<< ":" << ite->label << "\n";
					success = false;
				};
			}
		}
		it.file = file;
	}
};
struct skip_time_step {
	mutable bool success;
	mutable bool last;
	skip_time_step(bool l) : success(true), last(l) {};
	void operator()(OutputFile& it, bool last = false) const {
		std::istream * file = it.file;
		std::string line;
		std::ios::pos_type pos = file->tellg();
		for (int i=0; i < it.num_elements; i++) {
			std::getline(*file, line);
		}
		if (file->eof() && !last) {
			std::cerr << "Unexpected end of file skipping " << it.name << "; repositioning stream\n";
			file->seekg(pos);
			success = false;
		}
	}
};

struct read_time_step {
	mutable bool success;
	mutable bool last;
	read_time_step(bool l) : success(true), last(l) {};
	void operator()(OutputFile& it, bool last = false) const {
		std::istream * file = it.file;
		std::string line;
		int cur_element, prev_element;
		cur_element = 0; prev_element = -1;
		std::ios::pos_type pos = file->tellg();
		for (std::map<int, std::list<OutputElement>::iterator>::const_iterator ite = it.element_positions.begin();
			ite != it.element_positions.end(); ite++) {
			cur_element = ite->first;
			for (int i=1; i < cur_element - prev_element; i++) {
				file->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
				//file->ignore(10000,'\n');
				//skipline(file);
// 				std::getline(*file, line);
			}
			std::list<OutputElement>::iterator elit = ite->second;
			std::map<int,int>::iterator colit;
			int cur_col, prev_col;
			prev_col = 0;
			for (colit = elit->output_cols.begin(); colit != elit->output_cols.end(); colit++) {
				cur_col = colit->first;
				for (int i = 1; i < cur_col - prev_col; i++) {
					skipws(file);
					//*file >> x;
				}
				*file >> elit->outputs[colit->second];
				//elit->outputs[colit->second] = x;
				prev_col = cur_col;
			}

			file->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
			//file->ignore(10000,'\n');
			//skipline(file);
// 			std::getline(*file, line);
			prev_element = cur_element;
		}
		for (int i=1; i < it.num_elements - prev_element; i++) {
			file->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
			//file->ignore(10000,'\n');
			//skipline(file);
// 			std::getline(*file, line);
		}
		if (file->eof() && !last) {
			//std::cerr << "Unexpected end of file reading " << it.name << "; repositioning stream\n";
			file->seekg(pos);
			success = false;
		}
	}
};
 
struct write_time_step {
	mutable bool success;
	mutable std::ostream * file;
	write_time_step(std::ostream * f) : success(true), file(f) {};
	void operator()(OutputFile& it) const {
		for (std::list<OutputElement>::iterator elit = it.elements.begin();
			elit != it.elements.end(); elit++) {
			for (int i = 0; i<elit->num_outputs; i++) {
				*file << elit->outputs[i] << " ";
			}
		}
	}
};
 
struct close_file {
	void operator()(OutputFile& it) const {
		if (it.file != 0) {
			static_cast<std::ifstream*>(it.file)->close();
		}
	}
};
 
int main(int argc, char *argv[]) {
	std::cout.sync_with_stdio(false);
	std::cerr.sync_with_stdio(false);
	struct arguments argu;
	
	argp_parse (&argp_s, argc, argv, 0, 0, &argu);
	std::ostream *out;
	std::istream *commands;
	std::istream *in_log;
	std::istream *in_out;
	if (argu.output_file != 0) {
		std::string namefile(argu.output_file);
		namefile += ".dat";
		std::cerr << "out: " << namefile << "\n";
		out = new std::ofstream(namefile.c_str());
		if (out == 0) {
			std::cerr << "Unable to open output file " << namefile << "\n";
			return 1;
		}
	} else {
		out = &std::cout;
	}
	commands = new std::istringstream(argu.args[0]);
	if (!ParseCommands(*commands, std::cerr)) {
		return 1;
	}
	std::string basename(argu.args[1]);
	std::string namefile(basename);
	namefile += ".log";
	in_log = new std::ifstream(namefile.c_str());
	if (in_log == 0 || !static_cast<std::ifstream*>(in_log)->is_open()) {
		std::cerr << "Unable to open log file " << namefile << "\n";
		return 1;
	}
	namefile = basename;
	namefile += ".out";
	in_out = new std::ifstream(namefile.c_str());
	if (in_out == 0 || !static_cast<std::ifstream*>(in_out)->is_open()) {
		std::cerr << "Unable to open out file " << namefile << "\n";
		return 1;
	}
	
	if (!ParseLog(*in_log, std::cerr)) {
		return 1;
	}

	MbdynpostFlexLexer flexout(in_out, out);

	bool status;
	status = std::for_each(files.begin(), files.end(), allocate_file(basename)).success;
	if (!ParseOut(flexout, *in_out)) {
		return 1;
	}
	//da eliminare
	do {
		status = std::for_each(files.begin(), files.end(), read_time_step(false)).success;
		if (status && ((argu.t1 == 0.) || (current_time <= argu.t1) ) ) {
			if (current_time >= argu.t0) {
				*out << current_time << " ";
				status = std::for_each(files.begin(), files.end(), write_time_step(out)).success;
				*out << "\n";
				for (int i=0; i < argu.skip; i++) {
					if (!std::for_each(files.begin(), files.end(), skip_time_step(false)).success) {
						break;
					}
				}
				for (int i=1; i < (argu.skip + 1) * output_frequency; i++) {
					if (!ParseOut(flexout, *in_out)) {
						break;
					}
				}
			} else {
				for (int i=1; i < output_frequency; i++) {
					if (!ParseOut(flexout, *in_out)) {
						break;
					}
				}
			}
		} else {
			break;
		}
	} while (ParseOut(flexout, *in_out));
	std::for_each(files.begin(), files.end(), close_file());
//	WriteTimeData(*out, *in_out, *in_mov, argu.skip, strnode_num);
	if (argu.output_file) {
		static_cast<std::ofstream*>(out)->close();
	}
	return 0;
}

