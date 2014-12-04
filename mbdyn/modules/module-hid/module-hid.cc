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

#include <cmath>
#include <cfloat>
#include <cstdio>
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

	/* Scrive il contributo del Drive al file di restart */
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
m_pdLC(pdVal), m_pdB(pdVal + m_nLC)
{
}

JoystickDrive::~JoystickDrive(void)
{
	if (m_fd != -1) {
		close(m_fd);
	}
	m_fd = -1;
}

bool
JoystickDrive::get_one(void)
{
#if 0
	char buf[4*2];
	read(m_fd, (void *)&buf[0], sizeof(buf));
	short value = (buf[5] << 8) + buf[4];
	unsigned short key = (buf[7] << 8) + buf[6];

	unsigned short type = key & 0xF;
	unsigned short idx = (key & 0xF00) >> 16;
#endif

#if 1
	char buf[4];
	size_t n = read(m_fd, (void *)&buf[0], sizeof(buf));
	if (n == 0) {
		// done
		return true;

	} else if (n != 4) {
		// error
	}
	// discard by now

	int16_t value;
	n = read(m_fd, (void *)&value, sizeof(value));
	if (n != 1) {
		// error
	}

	uint8_t idx, type;
	n = read(m_fd, (void *)&type, sizeof(type));
	if (n != 1) {
		// error
	}

	n = read(m_fd, (void *)&idx, 1);
	if (n != 1) {
		// error
	}
#endif

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
	m_fd = open(sFileName.c_str(), O_RDONLY|O_NONBLOCK);
	if (m_fd == -1) {
		int save_errno = errno;
		silent_cerr("JoystickDrive(" << uLabel << "): "
			"unable to open device <" << sFileName << "> "
			"(" << save_errno << ": " << strerror(save_errno) << ")"
		<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// this probably contains the description of the channels...
	char buf[14*4];
	ssize_t n = read(m_fd, (void *)&buf[0], sizeof(buf));
	if (n != sizeof(buf)) {
		// error?
	}
}

/* Scrive il contributo del Drive al file di restart */
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
}


/* prototype of the functional object: reads a drive */
struct JoystickDR : public DriveRead {
	virtual Drive *
	Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP);
};

Drive *
JoystickDR::Read(unsigned uLabel, const DataManager* pDM, MBDynParser& HP)
{
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

