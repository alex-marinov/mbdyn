#!/bin/sh
# the next line restarts using wish \
	exec wish "$0" "$@"

# $Header$
# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2014
#    
# Pierangelo Masarati  <masarati@aero.polimi.it>
#      
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
#          
# Changing this copyright notice is forbidden.
#            
# This program is free software; you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License). 
# 
#                 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#                      
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

proc CreatePart { model name id x y z e1 e2 e3 } {
	global marker_id

	incr marker_id

	puts "!"
	puts "part create rigid_body name_and_position  &"
        puts "   part_name = .$model.$name  &"
        puts "   adams_id = $id  &"
        puts "   location = $x, $y, $z  &"
        puts "   orientation = $e1, $e2, $e3"
        puts "!"
        puts "marker create  &"
        puts "   marker_name = .$model.$name.cm  &"
        puts "   adams_id = $marker_id  &"
        puts "   location = $x, $y, $z  &"
        puts "   orientation = $e1, $e2, $e3"
        puts "!"
        puts "geometry create shape ellipsoid  &"
        puts "   ellipsoid_name = .$model.$name.sphere  &"
        puts "   center_marker = .$model.$name.cm  &"
        puts "   x_scale_factor = 0.1414213562  &"
        puts "   y_scale_factor = 0.1414213562  &"
        puts "   z_scale_factor = 0.1414213562"
	puts "!"
	puts "part attributes  &"
        puts "   part_name = .$model.$name  &"
        puts "   color = MAIZE  &"
        puts "   name_visibility = off"
}

set marker_id 1

set model [gets stdin]

while { 1 } {
	set name [gets stdin]
	if {[eof stdin]} {
		exit
	}

	set data [gets stdin]

	CreatePart $model $name [lindex $data 0] \
		[lindex $data 1] \
                [lindex $data 2] \
                [lindex $data 3] \
                [lindex $data 4] \
                [lindex $data 5] \
                [lindex $data 6]

}
