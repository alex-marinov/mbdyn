# Copyright (C) 2007 Pierangelo Masarati <masarati@aero.polimi.it>
#
# This file can be freely used and modified to be run with MBDyn
# provided this copyright notice is not altered nor removed.
#
# THIS MODEL IS PROVIDED AS IS, WITHOUT ANY WARRANTY OF ACCURACY
# OR FITNESS FOR ANY PURPOSE

# MBDyn model of the CART -  Controls Advanced Research Turbine

joint: CURR_BLADE + BLADE_PITCH, total joint,
	TEETER,
		position, reference, CURR_BLADE + BLADE_PITCH, null,
		position orientation, reference, CURR_BLADE + BLADE_PITCH, eye,
		rotation orientation, reference, CURR_BLADE + BLADE_PITCH, eye,
	CURR_BLADE + 0,
		position, reference, CURR_BLADE + BLADE_PITCH, null,
		position orientation, reference, CURR_BLADE + BLADE_PITCH, eye,
		rotation orientation, reference, CURR_BLADE + BLADE_ORIGIN, eye,
	position constraint,
		active, active, active,
			null,
	orientation constraint,
		active, active, active,
			1., 0., 0.,
				#file, 1, 2;
				reference, BLADE_PITCH;

# vim:ft=mbd
