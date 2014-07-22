#{

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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2014) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

    user defined: joint_id_slider1, 
# enable the next line to use the octave version of the element
        octave, "InLineFriction", embed octave, yes, octave search path, "..",
# enable the next line to use the C++ version of the element
        # inline friction,
        node1,
            node_id_ground,
            offset, reference, ref_id_slider1, null,
            orientation, reference, ref_id_slider1, eye,
        node2,
            node_id_slider,
            offset, reference, ref_id_slider1, dX1, dY1, dZ1,
        coulomb friction coefficient, muc,
        static friction coefficient, mus,
        sliding velocity coefficient, vs,
        sliding velocity exponent, i,
        micro slip displacement, delta,
        initial stiction state, zdelta0,
        initial stiction derivative, zP0,
        viscous friction coefficient, kv,
        stiction state equation scale, s;

    user defined: joint_id_slider2, 
# enable the next line to use the octave version of the element
        octave, "InLineFriction", embed octave, yes, octave search path, "..",
# enable the next line to use the C++ version of the element
        # inline friction,
        node1,
            node_id_ground,
            offset, reference, ref_id_slider2, null,
            orientation, reference, ref_id_slider2, eye,
        node2,
            node_id_slider,
            offset, reference, ref_id_slider2, dX2, dY2, dZ2,
        coulomb friction coefficient, muc,
        static friction coefficient, mus,
        sliding velocity coefficient, vs,
        sliding velocity exponent, i,
        micro slip displacement, delta,
        initial stiction state, zdelta0,
        initial stiction derivative, zP0,        
        viscous friction coefficient, kv,
        stiction state equation scale, s;

    bind: joint_id_slider1, user defined, prm_node_id_lambda1_1, string, "lambda1";
    bind: joint_id_slider1, user defined, prm_node_id_lambda2_1, string, "lambda2";
    bind: joint_id_slider1, user defined, prm_node_id_z_1,       string, "z";
    bind: joint_id_slider1, user defined, prm_node_id_zP_1,      string, "zP";
    bind: joint_id_slider1, user defined, prm_node_id_tau_1,     string, "tau";

    bind: joint_id_slider2, user defined, prm_node_id_lambda1_2, string, "lambda1";
    bind: joint_id_slider2, user defined, prm_node_id_lambda2_2, string, "lambda2";
    bind: joint_id_slider2, user defined, prm_node_id_z_2,       string, "z";
    bind: joint_id_slider2, user defined, prm_node_id_zP_2,      string, "zP";
    bind: joint_id_slider2, user defined, prm_node_id_tau_2,     string, "tau";

#}
# /*
# fprintf(stderr,"inline_friction2.elm: embed octave is enabled ...\n");
# warning("error","Octave:divide-by-zero");
# */
