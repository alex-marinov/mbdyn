#!/bin/sh
# the next line restarts using wish \
	exec wish "$0" "$@"


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
