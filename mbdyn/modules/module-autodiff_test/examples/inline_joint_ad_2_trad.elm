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
 AUTHOR: Reinhard Resch <r.resch@secop.com>
        Copyright (C) 2011(-2014) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

    joint: joint_id_slider1, in line,
        node_id_ground,
            position, reference, ref_id_slider1, null,
            orientation, reference, ref_id_slider1, 3, 1., 0., 0.,
                                       2, 0., 1., 0.,
        node_id_slider,
            offset, reference, ref_id_slider1, dX1, dY1, dZ1;

    joint: joint_id_slider2, in line,
        node_id_ground,
            position, reference, ref_id_slider2, null,
            orientation, reference, ref_id_slider2, 3, 1., 0., 0.,
                                       2, 0., 1., 0.,
        node_id_slider,
            offset, reference, ref_id_slider2, dX2, dY2, dZ2;
