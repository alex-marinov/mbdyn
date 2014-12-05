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

#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "dataman.h"
#include "filedrv.h"

class JoystickDrive : public FileDrive {
private:
	int m_fd;

	std::vector<bool> m_b_reset;
	integer m_nButtons;
	std::vector<doublereal> m_lc_scale;
	integer m_nLC;

	doublereal *m_pdLC, *m_pdB;

	// TODO: echo?
	// TODO: self-describing?

	bool get_one(void);
	void init(void);

public:
	JoystickDrive(unsigned int uL, const DriveHandler* pDH,
		const std::string& sFileName,
		integer nButtons, const std::vector<doublereal>& lc_scale,
		const std::vector<doublereal>& v0);

	virtual ~JoystickDrive(void);

	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void ServePending(const doublereal& t);
};

JoystickDrive::JoystickDrive(unsigned int uL, const DriveHandler* pDH,
	const std::string& sFileName,
	integer nButtons, const std::vector<doublereal>& lc_scale,
	const std::vector<doublereal>& v0)
: FileDrive(uL, pDH, sFileName, nButtons + lc_scale.size(), v0),
m_fd(-1),
m_b_reset(nButtons), m_nButtons(nButtons),
m_lc_scale(lc_scale), m_nLC(m_lc_scale.size()),
m_pdLC(pdVal + 1), m_pdB(m_pdLC + m_nLC)
{
}

JoystickDrive::~JoystickDrive(void)
{
	if (m_fd != -1) {
		close(m_fd);
	}
	m_fd = -1;
}

static int
fd_set_blocking(int fd, bool bBlocking)
{
	int flags, rc;

	flags = fcntl(fd, F_GETFL, 0);
	if (flags == -1) {
		return 0;
	}

	if (bBlocking) {
		flags &= ~O_NONBLOCK;

	} else {
		flags |= O_NONBLOCK;
	}

	rc = fcntl(fd, F_SETFL, flags);

	return (rc != -1);
}

bool
JoystickDrive::get_one(void)
{
	fd_set readfds;
	FD_ZERO(&readfds);
	FD_SET(m_fd, &readfds);
	struct timeval tv = { 0, 0 };
	int rc = select(m_fd + 1, &readfds, NULL, NULL, &tv);

	// std::cerr << "select=" << rc << std::endl;

	switch (rc) {
	case -1: {
		int save_errno = errno;
		char *err_msg = strerror(save_errno);

		silent_cout("JoystickDrive(" << uLabel << ", " << sFileName << "): select failed"
			<< " (" << save_errno << ": " << err_msg << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	case 0:
		return true;

	default:
		if (!FD_ISSET(m_fd, &readfds)) {
			silent_cout("JoystickDrive"
				"(" << sFileName << "): "
				"socket " << m_fd << " reset"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	char buf[8];
	fd_set_blocking(m_fd, true);
	ssize_t n2 = read(m_fd, (void *)&buf[0], sizeof(buf));
	fd_set_blocking(m_fd, false);
	if (n2 == -1) {
		// int save_errno = errno;
		// std::cerr << "recv=" << save_errno << " " << strerror(save_errno) << std::endl;
	}

	// uint32_t cnt = *((uint32_t *)&buf[0]);
	uint8_t type = (uint8_t)buf[6];
	uint8_t idx = (uint8_t)buf[7];
	int16_t value = *((int16_t *)&buf[4]);

	// std::cerr << "n2=" << n2 << " cnt=" << cnt << " value=" << int(value) << " type=" << unsigned(type) << " idx=" << unsigned(idx) << std::endl;

	switch (type) {
	case 1:
		// buttons
		if (idx < m_nButtons) {
			/*
			 * when a button is pressed, an event
			 * with non-zero value is recorded;
			 * when the button is released, an event
			 * with zero value is recorded.
			 *
			 * when we see an event with non-zero value,
			 * we set the value to 1., and clear
			 * the "reset" flag; when we see an event with
			 * zero value, we set the reset flag.
			 *
			 * if the button was not previously set,
			 * and it is pressed between two calls to
			 * ServePending(), at the end it is set
			 * and not scheduled for clear.
			 *
			 * if the button was set, and it is 
			 * released between two calls to ServePending(),
			 * it is scheduled for clear.
			 *
			 * if the button is pressed and released
			 * multiple times between two calls,
			 * it remains set to 1. and scheduled for clear.
			 */
			if (value) {
				m_pdB[idx] = 1.;
				m_b_reset[idx] = false;

			} else {
				m_b_reset[idx] = true;
			}
		}
		// else ignore
		break;

	case 2:
		// linear controls

		/*
		 * linear controls vary either between 0 and UINT16_MAX
		 * or between -INT16_MAX and INT16_MAX.
		 *
		 * FIXME: currently, we assume -INT16_MAX and INT16_MAX.
		 *
		 * we simply record the last value, and scale it.
		 */
		if (idx < m_nLC) {
			m_pdLC[idx] = (m_lc_scale[idx]*value)/INT16_MAX;
		}
		break;

	default:
		// ignore
		break;
	}

	return false;
}

void
JoystickDrive::init(void)
{
	m_fd = open(sFileName.c_str(), O_RDONLY);
	if (m_fd == -1) {
		int save_errno = errno;
		silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << ")::init(): "
			"unable to open device <" << sFileName << "> "
			"(" << save_errno << ": " << strerror(save_errno) << ")"
		<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// this probably contains the description of the channels...
	integer nB = 0, nLC = 0;
	uint8_t iBmax = 0, iLCmax = 0;
	for (int i; i < m_nLC + m_nButtons; i++) {
		char buf[8];
		ssize_t n = read(m_fd, (void *)&buf[0], sizeof(buf));

		// std::cerr << "read idx=" << i << " n=" << n << std::endl;

		if (n == -1) {
			int save_errno = errno;
			silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << ")::init(): "
				"read failed (" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		uint8_t type = (uint8_t)buf[6];
		uint8_t idx = (uint8_t)buf[7];
		// int16_t value = *((int16_t *)&buf[4]);
		if ((type & 0x7F) == 1) {
			nB++;
			iBmax = std::max(iBmax, idx);

		} else if ((type & 0x7F) == 2) {
			nLC++;
			iLCmax = std::max(iLCmax, idx);

		} else {
			silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << "): "
				"warning, unknown type " << int(type & 0x7F) << ", ignored" << std::endl); }

		// std::cerr << "    type=" << uint(type) << " idx=" << uint(idx) << " value=" << value << std::endl;
	}

	bool bFail(false);
	if (nB != m_nButtons) {
		silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << "): "
			"inconsistent number of buttons: expected " << m_nButtons << ", got " << nB << std::endl);
		bFail = true;
	}

	if (iBmax >= m_nButtons) {
		silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << "): "
			"inconsistent largest button index: expected " << m_nButtons - 1 << ", got " << unsigned(iBmax) << std::endl);
		bFail = true;
	}

	if (nLC != m_nLC) {
		silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << "): "
			"inconsistent number of linear controls: expected " << m_nLC << ", got " << nLC << std::endl);
		bFail = true;
	}

	if (iLCmax >= m_nLC) {
		silent_cerr("JoystickDrive(" << uLabel << ", " << sFileName << "): "
			"inconsistent largest linear control index: expected " << m_nLC - 1 << ", got " << unsigned(iLCmax) << std::endl);
		bFail = true;
	}

	fd_set_blocking(m_fd, false);
}

std::ostream&
JoystickDrive::Restart(std::ostream& out) const
{
	return out << "# not implemented yet" << std::endl;
}

void
JoystickDrive::ServePending(const doublereal& t)
{
	if (m_fd == -1) {
		init();
	}

	// reset if needed
	for (integer idx = 0; idx < m_nButtons; idx++) {
		if (m_b_reset[idx]) {
			m_pdB[idx] = 0.;
			m_b_reset[idx] = false;
		}
	}

	// flush buffer
	while (!get_one()) { NO_OP; };

	// std::cerr << "JoystickDrive(" << sFileName << ")" << std::endl;
	// for (int i = 0; i < m_nLC + m_nButtons; i++) {
		// std::cerr << "    V[" << i << "]=" << pdVal[i] << std::endl;
	// }
}

struct JoystickDR : public DriveRead {
	virtual Drive *
	Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP);
};

Drive *
JoystickDR::Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	joystick						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"	file: <label> , joystick ,					\n"
"		<file> ,						\n"
"		<number_of_buttons> ,					\n"
"		<number_of_linear_controls> ,				\n"
"			[ scale, <scale_factor> [ , ... ] ] ;		\n"
"									\n"
"	<file> ::= \"/dev/input/js0\"					\n"
"									\n"
"	<scale_factor>'s default to 1.0					\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	const char *s = HP.GetFileName();
	if (s == 0) {
		// error
	}
	std::string sFileName(s);

	integer nButtons = HP.GetInt();
	if (nButtons < 0) {
		// error
	}

	integer nLC = HP.GetInt();
	if (nLC <= 0) {
		// error
	}

	std::vector<doublereal> lc_scale(nLC);
	if (HP.IsKeyWord("scale")) {
		for (integer idx = 0; idx < nLC; idx++) {
			lc_scale[idx] = HP.GetReal();
		}

	} else {
		for (integer idx = 0; idx < nLC; idx++) {
			lc_scale[idx] = 1.;
		}
	}

	const std::vector<doublereal> v0;

	return new JoystickDrive(uLabel, pDM->pGetDrvHdl(), sFileName, nButtons, lc_scale, v0);
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	DriveRead *rf = new JoystickDR;

	if (!SetDriveData("joystick", rf)) {
		delete rf;

		silent_cerr("JoystickDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

