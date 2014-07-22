/* $Header$ */
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstring>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "sock.h"
#include "s2s.h"

#include "solman.h"
#include "linsol.h"

/* I/O classes (stream, using cin/cout) or socket,
 * using MBDyn native socket communication */

class BasicIO {
public:
	virtual ~BasicIO(void) { NO_OP; };
	virtual int ReadMeasures(s2s_t& s2s) = 0;
	virtual int SendControls(s2s_t& s2s) = 0;
};

class SocketBasicIO : public BasicIO {
public:
	int ReadMeasures(s2s_t& s2s);
	int SendControls(s2s_t& s2s);
};

static SocketBasicIO	socket_basic_IO;

class StreamBasicIO : public BasicIO {
public:
	int ReadMeasures(s2s_t& s2s);
	int SendControls(s2s_t& s2s);
};

static StreamBasicIO	stream_basic_IO;

int
SocketBasicIO::ReadMeasures(s2s_t& s2s)
{
	int len = recv(s2s.sock, (char *)&s2s.dbuf[0], sizeof(double)*s2s.nChannels, 0);

	switch (len) {
	case -1: {
		int		save_errno = errno;
		const char	*err_msg = strerror(save_errno);
		
		silent_cerr("recv(" << s2s.sock << ",\"" << s2s.buf << "\") "
			"failed (" << save_errno << ": " << err_msg << ")"
			<< std::endl);
		return -1;
	}

	case 0:
		return 0;

	default:
		if ((unsigned)len < sizeof(double)*s2s.nChannels) {
			silent_cerr("recv(" << s2s.sock << " \"" << s2s.buf << "\") "
				"returned only partial results"
				<< std::endl);
			return -1;
		}
		return 1;
	}
}

int
SocketBasicIO::SendControls(s2s_t& s2s)
{
	int	len;

	len = send(s2s.sock, (char *)&s2s.dbuf[0], sizeof(double)*s2s.nChannels, 0);

	if ( len < 0 ) {
		return len;
	}

	return len/sizeof(double);
}

int
StreamBasicIO::ReadMeasures(s2s_t& s2s)
{
	for (int i = 0; i < s2s.nChannels; i++) {
		std::cin >> s2s.dbuf[i];
		if (!std::cin) {
			if (i != 0) {
				return -1;
			}
			return 0;
		}
	}

	return 1;
}

int
StreamBasicIO::SendControls(s2s_t& s2s)
{
	static const char	*sep = " ";

	for (int i = 0; i < s2s.nChannels - 1; i++) {
		std::cout << s2s.dbuf[i] << sep;
	}
	std::cout << s2s.dbuf[s2s.nChannels - 1] << std::endl;

	return s2s.nChannels;
}

class IO {
private:
	s2s_t	s2s_measures;
	BasicIO	&measures;
	s2s_t	s2s_controls;
	BasicIO	&controls;

public:
	IO(void) : measures(stream_basic_IO), controls(stream_basic_IO) {};
	IO(int argc, char *argv[]);
	~IO(void) {};

	int nMeasures(void) const { return s2s_measures.nChannels; };
	int nControls(void) const { return s2s_controls.nChannels; };

	int Parse(int argc, char *argv[]);
	int Setup(int argc, char *argv[]);

	const std::vector<double> &Measures(void) const { return s2s_measures.dbuf; };
	std::vector<double> &Controls(void) { return s2s_controls.dbuf; };

	int ReadMeasures(void) { return measures.ReadMeasures(s2s_measures); };
	int SendControls(void) { return controls.SendControls(s2s_controls); };
};

IO::IO(int argc, char *argv[])
: measures(stream_basic_IO), controls(stream_basic_IO)
{
	if (Setup(argc, argv) == EXIT_FAILURE) {
		throw;
	}
}

int
IO::Parse(int argc, char *argv[])
{
	while (true) {
		int	opt = getopt(argc, argv, "c:f:H:m:");
		if (opt == EOF) {
			break;
		}

		s2s_t	*s2s = 0;

		switch (opt) {
		case 'H':
			if (strncasecmp(optarg, "measures:", STRLENOF("measures:")) == 0) {
				optarg += STRLENOF("measures:");
				s2s = &s2s_measures;

			} else if (strncasecmp(optarg, "controls:", STRLENOF("controls:")) == 0) {
				optarg += STRLENOF("controls:");
				s2s = &s2s_controls;

			} else {
				silent_cerr("unknown value \"" << optarg << "\" for -H switch" << std::endl);
				throw;
			}

			if (strncasecmp(optarg, "inet:", STRLENOF("inet:")) == 0) {
				optarg += STRLENOF("inet:");
				s2s->host = optarg;

			} else if (strncasecmp(optarg, "path:", STRLENOF("path:")) == 0) {
				optarg += STRLENOF("path:");
				s2s->path = optarg;

			} else if (strcasecmp(optarg, "stdin") == 0) {
				if (s2s != &s2s_measures) {
					silent_cerr("invalid value \"" << optarg << "\" for -H switch" << std::endl);
					throw;
				}
				break;

			} else if (strcasecmp(optarg, "stdout") == 0) {
				if (s2s != &s2s_controls) {
					silent_cerr("invalid value \"" << optarg << "\" for -H switch" << std::endl);
					throw;
				}
				break;

			} else {
				silent_cerr("unknown value \"" << optarg << "\" for -H switch" << std::endl);
				throw;
			}

			break;

		case 'c': {
			char	*next;

			s2s_controls.nChannels = strtoul(optarg, &next, 10);

			if (next == NULL || next[0] != '\0') {
				silent_cerr("invalid value \"" << optarg << "\" for -c switch" << std::endl);
				throw;
			}
			} break;

		case 'm': {
			char	*next;

			s2s_measures.nChannels = strtoul(optarg, &next, 10);

			if (next == NULL || next[0] != '\0') {
				silent_cerr("invalid value \"" << optarg << "\" for -m switch" << std::endl);
				throw;
			}
			} break;
		}
	}

	return 0;
}

int
IO::Setup(int argc, char *argv[])
{
	try {
		Parse(argc, argv);

	} catch (...) {
		return EXIT_FAILURE;
	}

	if (s2s_measures.path || s2s_measures.host) {
		try {
			s2s_measures.prepare();

		} catch (...) {
			return EXIT_FAILURE;
		}
	}

	if (s2s_controls.path || s2s_controls.host) {
		try {
			s2s_controls.prepare();

		} catch (...) {
			return EXIT_FAILURE;
		}
	}

	if (s2s_measures.sock == -1) {
		measures = stream_basic_IO;

		if (s2s_measures.nChannels == 0) {
			std::cin.getline(s2s_measures.buf, sizeof(s2s_measures.buf));

			std::istringstream	str(s2s_measures.buf);

			for (;; s2s_measures.nChannels++) {
				double	d;

				str >> d;

				if (!str) {
					break;
				}
			
				s2s_measures.dbuf.insert(s2s_measures.dbuf.end(), d);
			}

		} else {
			s2s_measures.dbuf.resize(s2s_measures.nChannels);
		}

	} else {
		measures = socket_basic_IO;

		if (s2s_measures.nChannels == 0) {
			silent_cerr("number of measures required in socket mode" << std::endl);
			return EXIT_FAILURE;
		}

		s2s_measures.dbuf.resize(s2s_measures.nChannels);
	}

	if (s2s_controls.sock == -1) {
		controls = stream_basic_IO;

	} else {
		controls = socket_basic_IO;
	}

	if (s2s_controls.nChannels == 0) {
		silent_cerr("number of controls required" << std::endl);
		return EXIT_FAILURE;
	}

	s2s_controls.dbuf.resize(s2s_controls.nChannels);


	return EXIT_SUCCESS;
	
}

class ConvergenceCheck {
protected:
	double tol;

public:
	ConvergenceCheck(int n, double t) : tol(t) {};
	virtual ~ConvergenceCheck(void) {};
	
	virtual bool Check(const std::vector<double> &measures) = 0;
};

class FFDConvergenceCheck : public ConvergenceCheck {
private:
	double rho;
	std::vector<double> prevMeasures;
	std::vector<double> prevDiff;
	std::vector<double> diff;

public:
	FFDConvergenceCheck(int n, double t, double r);
	~FFDConvergenceCheck(void) {};

	bool Check(const std::vector<double> &measures);
};

FFDConvergenceCheck::FFDConvergenceCheck(int n, double t, double r)
: ConvergenceCheck(n, t),
rho(r),
prevMeasures(n, 0.),
prevDiff(n, 0.),
diff(n, 0.)
{
	return;
}

bool
FFDConvergenceCheck::Check(const std::vector<double> &measures)
{
	double	d = 0.;

	for (unsigned i = 0; i < measures.size(); i++) {
		/*
			diff = rho * (measures - prevMeasures) + (1 - rho) * prevDiff
		 */
		diff[i] = (1 - rho)*prevDiff[i];
		prevDiff[i] = measures[i] - prevMeasures[i];
		diff[i] += rho*prevDiff[i];
		prevMeasures[i] = measures[i];

		d += diff[i]*diff[i];
	}

	if (std::sqrt(d) < tol) {
		return true;
	}

	return false;
}

class TrimEval {
protected:
	IO			&io;
	ConvergenceCheck	&cc;

public:
	TrimEval(IO &io, ConvergenceCheck &cc) : io(io), cc(cc) {};
	~TrimEval(void) {};

	int FuncEval(std::vector<double> &X, std::vector<double> &F, int nSteps, int maxSteps);
};

int
TrimEval::FuncEval(std::vector<double> &X, std::vector<double> &F, int nSteps, int maxSteps)
{
	if (X.size() != io.Controls().size()) {
		return -1;
	}

	std::vector<double>	Xref(X.size());
	std::vector<double>	deltaX(X.size());

	for (unsigned int i = 0; i < X.size(); i++) {
		Xref[i] = io.Controls()[i];
		deltaX[i] = X[i] - Xref[i];
	}

	int	currStep;
	for (currStep = 0; currStep < nSteps; currStep++) {
		for (unsigned int i = 0; i < X.size(); i++) {
			io.Controls()[i] = Xref[i] + (deltaX[i]*currStep)/nSteps;
		}

		io.SendControls();
		io.ReadMeasures();
	}

	for (; currStep < maxSteps; currStep++) {
		if (cc.Check(io.Measures())) {
			return 1;
		}

		io.SendControls();
		io.ReadMeasures();
	}

	return 0;
}

class NRTrim {
private:
	IO			&io;
	ConvergenceCheck	*pcc;

public:
	NRTrim(IO &io, int n, double t);
	~NRTrim(void) {};

	int DoTrim(const std::vector<double> expectedMeasures,
			const std::vector<double> incrControls,
			std::vector<double> controls);
};

NRTrim::NRTrim(IO &io, int n, double t)
: io(io),
pcc(0)
{
	pcc = new FFDConvergenceCheck(n, t, .99);
}

int
NRTrim::DoTrim(const std::vector<double> expectedMeasures,
		const std::vector<double> incrControls,
		std::vector<double> controls)
{
	TrimEval		te(io, *pcc);
	std::vector<double>	Xref(incrControls.size()),
				X(incrControls.size());
	std::vector<double>	Fref(expectedMeasures.size()),
				F(expectedMeasures.size(), 0.);
	LinSol			LS;
	SolutionManager		*psm(LS.GetSolutionManager(expectedMeasures.size()));
	MatrixHandler		*pM = psm->pMatHdl();
	VectorHandler		*pR = psm->pResHdl();
	VectorHandler		*pX = psm->pSolHdl();

	if (expectedMeasures.size() != incrControls.size()) {
		return -1;
	}

	for (unsigned int i = 0; i < incrControls.size(); i++) {
		Xref[i] = io.Controls()[i];
	}

	for (unsigned int i = 0; i < expectedMeasures.size(); i++) {
		Fref[i] = io.Measures()[i];
	}

	for (unsigned int c = 0; c < incrControls.size(); c++) {
		for (unsigned int i = 0; i < incrControls.size(); i++) {
			X[i] = Xref[i];
		}
		X[c] += incrControls[c];

		te.FuncEval(X, F, 10, 1000);

		for (unsigned int i = 0; i < expectedMeasures.size(); i++) {
			pM->PutCoef(i + 1, c + 1, F[i] - Fref[i]);
		}
	}

	while (!pcc->Check(io.Measures())) {
		for (unsigned int i = 0; i < expectedMeasures.size(); i++) {
			pR->PutCoef(i + 1, expectedMeasures[i] - Fref[i]);
		}

		psm->Solve();

		for (unsigned int i = 0; i < incrControls.size(); i++) {
			X[i] = Xref[i] + pX->operator()(i + 1);
			Xref[i] = X[i];
		}

		te.FuncEval(X, F, 10, 1000);

		for (unsigned int i = 0; i < expectedMeasures.size(); i++) {
			Fref[i] = F[i];
		}
	}

	return 0;
}

int
main(int argc, char *argv[])
{
	IO			io(argc, argv);
	if (io.nControls() != 4) {
		silent_cerr("need 4 controls" << std::endl);
		throw;
	}

	NRTrim			trim(io, 4, 1.e-6);

	std::vector<double>	F(4, 0.);
	std::vector<double>	I(4, 0.);
	std::vector<double>	X(4, 0.);

	F[0] = .6;
	F[1] = .0;
	F[2] = .0;
	F[3] = 260.;

	I[0] = .5/180*M_PI;
	I[1] = .5/180*M_PI;
	I[2] = .5/180*M_PI;
	I[3] = .1/180*M_PI;

	trim.DoTrim(F, I, X);
}

