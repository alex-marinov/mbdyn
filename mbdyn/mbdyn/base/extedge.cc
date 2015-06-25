/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2015
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

#include "dataman.h"
#include "extedge.h"
#include "except.h"
#include "solver.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <unistd.h>
#include <cerrno>
#include <cstring>

int
mbedge_goto_eol(std::istream& fin, char *buf, size_t bufsiz)
{
	size_t i;

	for (i = 0; i < bufsiz; i++) {
		buf[i] = fin.get();

		if (!fin) {
			return -1;
		}

		if (buf[i] == '\n') {
			buf[i] = '\0';
			break;
		}
	}

	if (i == bufsiz) {
		return -1;
	}

	return 0;
}

char *
mbedge_eat_sep(char *buf, size_t& buflen)
{
	size_t	len = buflen;

	ASSERT(buflen > 0);

	while (len > 0) {
		if (buf[0] == '\n' || buf[0] == '\0' || (buf[0] != ',' && !std::isspace(buf[0]))) {
			buflen = len;
			return buf;
		}

		len--;
		buf++;
	}

	return 0;
}

char *
mbedge_eat_field(char *buf, size_t& buflen, const char *val)
{
	ASSERT(buflen > 0);

	size_t vallen = strlen(val);
	if (buflen < vallen) {
		return 0;
	}

	if (strncmp(buf, val, vallen) == 0) {
		buflen -= vallen;
		return &buf[vallen];
	}

	return 0;
}

#if 0
char *
mbedge_eat_field(char *buf, size_t& buflen)
{
	size_t len = buflen;

	ASSERT(buflen > 0);

	while (len > 0) {
		if (buf[0] == ',' || std::isspace(buf[0])) {
			buflen = len;
			return buf;
		}

		len--;
		buf++;
	}

	return 0;
}
#endif

/* ExtFileHandlerEDGE - begin */

ExtFileHandlerEDGE::ExtFileHandlerEDGE(std::string& fflagname,
	std::string& fdataname, mbsleep_t SleepTime, std::streamsize Precision)
: ExtFileHandlerBase(SleepTime, Precision),
fflagname(fflagname), fdataname(fdataname),
bReadForces(true)
{
	NO_OP;
}

ExtFileHandlerEDGE::~ExtFileHandlerEDGE(void)
{
	int cnt = -1;
	switch (CheckFlag(cnt)) {
	case EDGE_READ_READY:
	case EDGE_GOTO_NEXT_STEP:
		SendFlag(EDGE_QUIT);
		break;

	default:
		break;
	}
}

static const char *sEDGEcmd2str[] = {
	"INITIALIZING",
	"BUSY",
	"READ_READY",
	"MBDYN_WRITE_DONE",
	"GOTO_NEXT_STEP",
	"QUIT",
	0
};

const char *
ExtFileHandlerEDGE::EDGEcmd2str(int cmd) const
{
	if (cmd < 0 || cmd >= EDGE_LAST) {
		return "UNKNOWN";
	}

	return sEDGEcmd2str[cmd];
}

ExtFileHandlerEDGE::EDGEcmd
ExtFileHandlerEDGE::CheckFlag(int& cnt)
{
	int cmd;
	unsigned lineno = 0;

	cnt++;

	infile.open(fflagname.c_str());
	if (!infile && cnt == -1) {
		return EDGE_QUIT;
	}

	// Needs tuning?
	const int max_retries = 10;
	int retrying = 0;
retry:;

	if (SleepTime > 0) {
		// WARNING: loops forever
		// add optional, configurable limit?
		for (; !infile; cnt++) {
			silent_cout("flag file \"" << fflagname.c_str() << "\" "
				"missing, try #" << cnt << "; "
				"sleeping " << SleepTime << " s" << std::endl); 
			if (mbdyn_stop_at_end_of_iteration()) {
				cmd = EDGE_QUIT;
				goto done;
			}

			mbsleep(&SleepTime);

			infile.clear();
			infile.open(fflagname.c_str());
		}

	} else if (!infile) {
		// error
		silent_cerr("flag file \"" << fflagname.c_str() << "\" "
			"missing" << std::endl); 
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	UPDATE,N,0,0,1
	FLAG,I,1,1,0
	<flag value>
	*/

	// old syntax without this block
	while (true) {
		char buf[BUFSIZ], *p;
		if (mbedge_goto_eol(infile, buf, sizeof(buf))) {
			if (!infile) {
				break;
			}
		}

		lineno++;

		if (buf[0] == '*') {
			continue;
		}

		if (strncasecmp(buf, "UPDATE", STRLENOF("UPDATE")) == 0) {
			size_t buflen = sizeof(buf) - STRLENOF("UPDATE");
			p = &buf[0] + STRLENOF("UPDATE");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "N");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"N\" "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << &buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (p[0] != '\0' && p[0] != '\n') {
				silent_cerr("ExtRigidForceEDGE: no line terminator "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (mbedge_goto_eol(infile, buf, sizeof(buf))) {
				silent_cerr("ExtRigidForceEDGE: unable to get \"FLAG\" line "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			lineno++;
				
			if (strncasecmp(buf, "FLAG", STRLENOF("FLAG")) != 0) {
				silent_cerr("ExtRigidForceEDGE: \"FLAG\" line expected "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			buflen = sizeof(buf) - STRLENOF("FLAG");
			p = &buf[0] + STRLENOF("FLAG");

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "I");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"I\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "1");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"1\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_field(p, buflen, "0");
			if (p == 0) {
				silent_cerr("ExtFileHandlerEDGE: unable to skip field \"0\" "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			p = mbedge_eat_sep(p, buflen);
			if (p == 0) {
				silent_cerr("ExtRigidForceEDGE: unable to skip separator "
					"at line=" << lineno << ", \"" << buf[sizeof(buf) - buflen] << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			if (p[0] != '\0' && p[0] != '\n') {
				silent_cerr("ExtRigidForceEDGE: no line terminator "
					"at line=" << lineno << ", \"" << p << "\"" << std::endl);
				if (retrying < max_retries) {
					retrying++;
					goto retry;
				}
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			lineno++;
				
			break;
		}
	}

	infile >> cmd;

done:;
	silent_cerr("<<< checking flag file \"" << fflagname.c_str() << "\": "
		"cmd=" << cmd << " (" << EDGEcmd2str(cmd) << ")"
		<< std::endl);

	infile.close();
	infile.clear();

	return EDGEcmd(cmd);
}

void
ExtFileHandlerEDGE::SendFlag(EDGEcmd cmd)
{
	// write to a temporary, unique file and then rename to fflagname
#ifdef HAVE_MKSTEMP
	char ftmpname[] = "mbedgeXXXXXX";
	int filedes = mkstemp(ftmpname);
	FILE *fd = fdopen(filedes, "w");
#else // ! HAVE_MKSTEMP
	std::string ftmpn(fflagname + ".tmp");
	const char *ftmpname = ftmpn.c_str();
	FILE *fd = fopen(ftmpname, "w");
#endif // ! HAVE_MKSTEMP

	fprintf(fd, "UPDATE,N,0,0,1\n");
	fprintf(fd, "FLAG,I,1,1,0\n");
	fprintf(fd, "%d", int(cmd));
	fclose(fd);
retry:;
	if (rename(ftmpname, fflagname.c_str()) == -1) {
		switch (errno) {
		case EBUSY: {
			mbsleep_t timeout = mbsleep_init(1);
			// TODO: configurable?
			mbsleep(&timeout);

			if (mbdyn_stop_at_end_of_iteration()) {
				// ultimately give up
				unlink(ftmpname);
				return;
			}

			goto retry;
			}

		default: {
			int save_errno = errno;
			silent_cerr("unable to rename flag file "
				"\"" << fdataname.c_str() << "\" "
				"(errno=" << save_errno << ": "
				<< strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	silent_cerr(">>> sending flag file \"" << fflagname.c_str() << "\": "
		"cmd=" << cmd << " (" << EDGEcmd2str(cmd) << ")"
		<< std::endl);

}

bool
ExtFileHandlerEDGE::Prepare_pre(void)
{
	return true;
}

void
ExtFileHandlerEDGE::Prepare_post(bool ok)
{
	NO_OP;
}

void
ExtFileHandlerEDGE::AfterPredict(void)
{
	bReadForces = true;
}

bool
ExtFileHandlerEDGE::Send_pre(SendWhen when)
{
	int cnt = 0;
	EDGEcmd cmd = EDGE_UNKNOWN;

	if (!bReadForces) {
		goto bad;
	}

	// open data file for writing
retry:;
	cmd = CheckFlag(cnt);
	switch (cmd) {
	case EDGE_INITIALIZING:
	case EDGE_BUSY:
		if (SleepTime > 0) {
			silent_cout("flag file \"" << fflagname.c_str() << "\": "
				"cmd=" << cmd << " (" << EDGEcmd2str(cmd) << ")"
				" try #" << cnt << "; "
				"sleeping " << SleepTime << " s" << std::endl); 

			mbsleep(&SleepTime);
		}

		if (mbdyn_stop_at_end_of_iteration()) {
			bReadForces = false;
			goto bad;
		}

		goto retry;

	case EDGE_READ_READY:
	case EDGE_GOTO_NEXT_STEP:
		outfile.open(fdataname.c_str());
		if (!outfile) {
			int save_errno = errno;
			silent_cerr("unable to open data file "
				"\"" << fdataname.c_str() << "\" "
				"for output (" << errno << ": "
				<< strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		outfile.setf(std::ios::scientific);
		outfile.precision(Precision);
		break;

	default:
bad:;
		outfile.setstate(std::ios_base::badbit);
		break;
	}

	return outfile.good();
}

void
ExtFileHandlerEDGE::Send_post(SendWhen when)
{
	if (!bReadForces) {
		return;
	}

	// close data file after writing
	if (outfile.is_open()) {
		outfile.close();
	}

	if (when == SEND_AFTER_CONVERGENCE) {
#if 0
		// This stops EDGE's subiterations
		// and advances to next step
		outfile.open("update.ainp");
		outfile <<
			"UPDATE,N,0,0,1\n"
			"IBREAK,I,1,1,0\n"
			"5\n";
		outfile.close();
#endif

		bReadForces = true;
	}

	SendFlag(EDGE_MBDYN_WRITE_DONE);
}

bool
ExtFileHandlerEDGE::Recv_pre(void)
{
	for (int cnt = 0; ;) {

		EDGEcmd cmd = CheckFlag(cnt);
		switch (cmd) {
		case EDGE_READ_READY:
			goto done;

		case EDGE_GOTO_NEXT_STEP:
			// FIXME: EDGE is done; do not read forces,
			// keep using old
			bReadForces = false;
			goto done;

		case EDGE_QUIT:
			silent_cout("EDGE requested end of simulation"
				<< std::endl);
			mbdyn_set_stop_at_end_of_time_step();
			bReadForces = false;
			goto done;

		default:
			break;
		}

		if (mbdyn_stop_at_end_of_iteration()) {
			bReadForces = false;
			goto done;
		}

		// WARNING: loops forever
		// add optional, configurable limit?
		if (SleepTime > 0) {
			silent_cout("flag file \"" << fflagname.c_str() << "\": "
				"cmd=" << cmd << " (" << EDGEcmd2str(cmd) << ")"
				" try #" << cnt << "; "
				"sleeping " << SleepTime << " s" << std::endl); 

			mbsleep(&SleepTime);
		} else {
			silent_cout("flag file \"" << fflagname.c_str() << "\": "
				"cmd=" << cmd << " (" << EDGEcmd2str(cmd) << ")"
				" try #" << cnt << std::endl); 
		}
	}

done:;
	if (bReadForces) {
		infile.open(fdataname.c_str());
		if (!infile) {
			silent_cerr("unable to open data file "
				"\"" << fdataname.c_str() << "\" "
				"for input" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		infile.setstate(std::ios_base::badbit);
	}

	return infile.good();
}

bool
ExtFileHandlerEDGE::Recv_post(void)
{
	if (infile.is_open()) {
		infile.close();
	}

	return !bReadForces;
}

std::ostream *
ExtFileHandlerEDGE::GetOutStream(void)
{
	return &outfile;
}

std::istream *
ExtFileHandlerEDGE::GetInStream(void)
{
	return &infile;
}

/* ExtFileHandlerEDGE - end */

ExtFileHandlerBase *
ReadExtFileHandlerEDGE(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH = 0;

	const char *s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForceEDGE(" << uLabel << "): "
			"unable to get flag file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::string fflagname = s;

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForceEDGE(" << uLabel << "): "
			"unable to get data file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::string fdataname = s;

	mbsleep_t SleepTime = mbsleep_init(1);
	std::streamsize Precision = 0;
	ReadExtFileParams(pDM, HP, uLabel, SleepTime, Precision);

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtFileHandlerEDGE,
		ExtFileHandlerEDGE(fflagname, fdataname, SleepTime, Precision));

	return pEFH;
}

