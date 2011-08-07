/*
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2011) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

// FIXME: this is a dirty trick needed to compile
#define real mbdyn_real_type
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */

#include <matvec3.h>
#include <matvecexp.h>
#include <Rot.hh>
#undef real

#define BOUNDS_CHECKING

#include <octave/oct.h>
//
DEFUN_DLD(rotation_matrix_to_rotation_vector,args,nargout,"[phi] = rotation_matrix_to_rotation_vector(R)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 3 )
	{
		error("R must be a 3x3 matrix!");
		return retval;
	}

	Matrix R_oct = args(0).matrix_value();

	Mat3x3 R;

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R(i+1,j+1) = R_oct(i,j);

	Vec3 phi(RotManip::VecRot(R));

	ColumnVector phi_oct(3);

	for ( int i = 0; i < 3; ++i )
		phi_oct(i) = phi(i+1);

	retval.append(octave_value(phi_oct));

	return retval;
}

DEFUN_DLD(rotation_vector_to_rotation_matrix,args,nargout,"[R] = rotation_vector_to_rotation_matrix(phi)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 1 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	ColumnVector phi_oct = args(0).matrix_value().column(0);

	Vec3 phi;

	for ( int i = 0; i < 3; ++i )
		phi(i+1) = phi_oct(i);

	Mat3x3 R(RotManip::Rot(phi));

	Matrix R_oct(3,3);

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R_oct(i,j) = R(i+1,j+1);

	retval.append(octave_value(R_oct));

	return retval;
}

DEFUN_DLD(rotation_matrix_to_euler123,args,nargout,"[phi] = rotation_matrix_to_euler123(R)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 3 )
	{
		error("R must be a 3x3 matrix!");
		return retval;
	}

	Matrix R_oct = args(0).matrix_value();

	Mat3x3 R;

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R(i+1,j+1) = R_oct(i,j);

	Vec3 phi(MatR2EulerAngles123(R));

	ColumnVector phi_oct(3);

	for ( int i = 0; i < 3; ++i )
		phi_oct(i) = phi(i+1);

	retval.append(octave_value(phi_oct));

	return retval;
}

DEFUN_DLD(euler123_to_rotation_matrix,args,nargout,"[R] = euler123_to_rotation_matrix(phi)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 1 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	ColumnVector phi_oct = args(0).matrix_value().column(0);

	Vec3 phi;

	for ( int i = 0; i < 3; ++i )
		phi(i+1) = phi_oct(i);

	const Mat3x3 R(EulerAngles123_2MatR(phi));

	Matrix R_oct(3,3);

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R_oct(i,j) = R(i+1,j+1);

	retval.append(octave_value(R_oct));

	return retval;
}

DEFUN_DLD(rotation_matrix_to_euler313,args,nargout,"[phi] = rotation_matrix_to_euler313(R)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 3 )
	{
		error("R must be a 3x3 matrix!");
		return retval;
	}

	Matrix R_oct = args(0).matrix_value();

	Mat3x3 R;

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R(i+1,j+1) = R_oct(i,j);

	Vec3 phi(MatR2EulerAngles313(R));

	ColumnVector phi_oct(3);

	for ( int i = 0; i < 3; ++i )
		phi_oct(i) = phi(i+1);

	retval.append(octave_value(phi_oct));

	return retval;
}

DEFUN_DLD(euler313_to_rotation_matrix,args,nargout,"[R] = euler313_to_rotation_matrix(phi)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 1 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	ColumnVector phi_oct = args(0).matrix_value().column(0);

	Vec3 phi;

	for ( int i = 0; i < 3; ++i )
		phi(i+1) = phi_oct(i);

	const Mat3x3 R(EulerAngles313_2MatR(phi));

	Matrix R_oct(3,3);

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R_oct(i,j) = R(i+1,j+1);

	retval.append(octave_value(R_oct));

	return retval;
}

DEFUN_DLD(rotation_matrix_to_euler321,args,nargout,"[phi] = rotation_matrix_to_euler321(R)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 3 )
	{
		error("R must be a 3x3 matrix!");
		return retval;
	}

	Matrix R_oct = args(0).matrix_value();

	Mat3x3 R;

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R(i+1,j+1) = R_oct(i,j);

	Vec3 phi(MatR2EulerAngles321(R));

	ColumnVector phi_oct(3);

	for ( int i = 0; i < 3; ++i )
		phi_oct(i) = phi(i+1);

	retval.append(octave_value(phi_oct));

	return retval;
}

DEFUN_DLD(euler321_to_rotation_matrix,args,nargout,"[R] = euler321_to_rotation_matrix(phi)\n")
{
	octave_value_list retval;

	int nargin = args.length();

	if ( nargin < 1 )
	{
		print_usage();
		return retval;
	}

	if ( !args(0).is_matrix_type() || args(0).rows() != 3 || args(0).columns() != 1 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	ColumnVector phi_oct = args(0).matrix_value().column(0);

	Vec3 phi;

	for ( int i = 0; i < 3; ++i )
		phi(i+1) = phi_oct(i);

	const Mat3x3 R(EulerAngles321_2MatR(phi));

	Matrix R_oct(3,3);

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			R_oct(i,j) = R(i+1,j+1);

	retval.append(octave_value(R_oct));

	return retval;
}
/*
%!test
%! for i = 1:100
%!    phi = pi * ( 2 * rand() - 1 );
%!    n = ( 2 * rand(3,1) - 1 );
%!    n /= norm(n);
%!    Phi = n * phi;
%! 	  R = rotation_vector_to_rotation_matrix(Phi);
%!    assert(R * R.',eye(3),sqrt(eps));
%!    assert(R.' * R,eye(3),sqrt(eps));
%!	  Phi2 = rotation_matrix_to_rotation_vector(R);
%! 	  assert(Phi2,Phi,sqrt(eps)*norm(Phi));
%! endfor

%!test
%! for i = 1:100
%!    phi = pi/2 * ( 2 * rand() - 1 );
%!    n = ( 2 * rand(3,1) - 1 );
%!    n /= norm(n);
%!    Phi = n * phi;
%! 	  R = euler123_to_rotation_matrix(Phi);
%!    assert(R * R.',eye(3),sqrt(eps));
%!    assert(R.' * R,eye(3),sqrt(eps));
%!	  Phi2 = rotation_matrix_to_euler123(R);
%! 	  assert(Phi2,Phi,sqrt(eps)*norm(Phi));
%! endfor

%!test
%!    Phi = [ -pi/4, 0, 0,		pi/3,   0;
%!			   pi/2, 0, pi/2,	pi/10,	0;
%!			   pi/3, 0, 0,		0,		pi/2 ];
%!
%!    for i=1:columns(Phi)
%! 	  	R = euler313_to_rotation_matrix(Phi(:,i));
%!    	assert(R * R.',eye(3),sqrt(eps));
%!    	assert(R.' * R,eye(3),sqrt(eps));
%!	  	Phi2 = rotation_matrix_to_euler313(R);
%! 	  	assert(fmod(2*pi+Phi2,2*pi),fmod(2*pi+Phi(:,i),2*pi),sqrt(eps)*norm(Phi));
%!    endfor

%!test
%! for i = 1:100
%!    phi = pi/2 * ( 2 * rand() - 1 );
%!    n = ( 2 * rand(3,1) - 1 );
%!    n /= norm(n);
%!    Phi = n * phi;
%! 	  R = euler321_to_rotation_matrix(Phi);
%!    assert(R * R.',eye(3),sqrt(eps));
%!    assert(R.' * R,eye(3),sqrt(eps));
%!	  Phi2 = rotation_matrix_to_euler321(R);
%! 	  assert(Phi2,Phi,sqrt(eps)*norm(Phi));
%! endfor
*/
