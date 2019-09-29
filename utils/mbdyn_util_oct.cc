/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 AUTHOR: Reinhard Resch <r.resch@a1.net>
        Copyright (C) 2011(-2019) all rights reserved.

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
        error("R must be a 3x3xN matrix!");
		return retval;
	}

    const NDArray R_oct = args(0).array_value();
    const octave_idx_type N = R_oct.ndims() >= 3 ? R_oct.dim3() : 1;
    Matrix phi_oct(3, N);

	Mat3x3 R;

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R(i + 1, j + 1) = R_oct(i, j, k);

	Vec3 phi(RotManip::VecRot(R));

	for ( int i = 0; i < 3; ++i )
            phi_oct(i, k) = phi(i + 1);
    }

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

    if ( !args(0).is_matrix_type() || args(0).rows() != 3 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	Vec3 phi;
    const Matrix phi_oct = args(0).matrix_value();
    const octave_idx_type N = phi_oct.columns();
    NDArray R_oct(dim_vector(3, 3, N));

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
            phi(i + 1) = phi_oct(i, k);

	Mat3x3 R(RotManip::Rot(phi));

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R_oct(i, j, k) = R(i + 1, j + 1);
    }

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

    const NDArray R_oct = args(0).array_value();
    const octave_idx_type N = R_oct.ndims() >= 3 ? R_oct.dim3() : 1;
    Matrix phi_oct(3, N);        
	Mat3x3 R;

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R(i + 1, j + 1) = R_oct(i, j, k);

	Vec3 phi(MatR2EulerAngles123(R));

	for ( int i = 0; i < 3; ++i )
            phi_oct(i, k) = phi(i + 1);
    }

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

    if ( !args(0).is_matrix_type() || args(0).rows() != 3)
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	Vec3 phi;
    const Matrix phi_oct = args(0).matrix_value();
    const octave_idx_type N = phi_oct.columns();
    NDArray R_oct(dim_vector(3, 3, N));

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
            phi(i + 1) = phi_oct(i, k);

	const Mat3x3 R(EulerAngles123_2MatR(phi));

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R_oct(i, j, k) = R(i + 1, j + 1);
    }

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

    const NDArray R_oct = args(0).array_value();
    const octave_idx_type N = R_oct.ndims() >= 3 ? R_oct.dim3() : 1;
    Matrix phi_oct(3, N);        
	Mat3x3 R;

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R(i + 1, j + 1) = R_oct(i, j, k);

	Vec3 phi(MatR2EulerAngles313(R));

	for ( int i = 0; i < 3; ++i )
            phi_oct(i, k) = phi(i + 1);
    }

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

    if ( !args(0).is_matrix_type() || args(0).rows() != 3)
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	Vec3 phi;
    const Matrix phi_oct = args(0).matrix_value();
    const octave_idx_type N = phi_oct.columns();
    NDArray R_oct(dim_vector(3,3, N));

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
            phi(i + 1) = phi_oct(i, k);

	const Mat3x3 R(EulerAngles313_2MatR(phi));

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R_oct(i, j, k) = R(i + 1, j + 1);
    }

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

    const NDArray R_oct = args(0).array_value();
    const octave_idx_type N = R_oct.ndims() >= 3 ? R_oct.dim3() : 1;
    Matrix phi_oct(3, N);
	Mat3x3 R;

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R(i + 1, j + 1) = R_oct(i, j, k);

	Vec3 phi(MatR2EulerAngles321(R));

	for ( int i = 0; i < 3; ++i )
            phi_oct(i, k) = phi(i + 1);
    }

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

    if ( !args(0).is_matrix_type() || args(0).rows() != 3 )
	{
		error("Phi must be a 3x1 matrix!");
		return retval;
	}

	Vec3 phi;
    const Matrix phi_oct = args(0).matrix_value();
    const octave_idx_type N = phi_oct.columns();
    
    NDArray R_oct(dim_vector(3, 3, N));

    for (octave_idx_type k = 0; k < N; ++k) {
	for ( int i = 0; i < 3; ++i )
            phi(i + 1) = phi_oct(i, k);

	const Mat3x3 R(EulerAngles321_2MatR(phi));

	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
                R_oct(i, j, k) = R(i + 1, j + 1);
    }

	retval.append(octave_value(R_oct));

	return retval;
}

/*
%!test
%! func_R = {@euler313_to_rotation_matrix, @euler321_to_rotation_matrix, @euler123_to_rotation_matrix, @rotation_vector_to_rotation_matrix};
%! func_Phi = {@rotation_matrix_to_euler313, @rotation_matrix_to_euler321, @rotation_matrix_to_euler123, @rotation_matrix_to_rotation_vector};
%!    Phi = [ -pi/4, 0, 0,		pi/3,   0;
%!			   pi/2, 0, pi/2,	pi/10,	0;
%!			   pi/3, 0, 0,		0,		pi/2 ];
%! for k=1:numel(func_R)
%! R = feval(func_R{k}, Phi);
%! Phi2 = feval(func_Phi{k}, R);
%!    for i=1:columns(Phi)
%!   assert(R(:, :, i) * R(:, :, i).',eye(3),sqrt(eps));
%!   assert(R(:, :, i).' * R(:, :, i),eye(3),sqrt(eps));
%!   assert(mod(2*pi+Phi2(:, i),2*pi),mod(2*pi+Phi(:,i),2*pi),sqrt(eps)*norm(Phi(:,i)));
%! endfor
%! for i=1:columns(Phi)
%!   R = feval(func_R{k}, Phi(:, i));
%!   Phi2 = feval(func_Phi{k}, R);
%!    	assert(R * R.',eye(3),sqrt(eps));
%!    	assert(R.' * R,eye(3),sqrt(eps));
%!   assert(mod(2*pi+Phi2,2*pi),mod(2*pi+Phi(:,i),2*pi),sqrt(eps)*norm(Phi(:,i)));
%! endfor
%! if (k != 1)
%!   assert(feval(func_Phi{k}, eye(3)), zeros(3, 1));
%! endif
%! assert(feval(func_R{k}, zeros(3, 1)), eye(3));
%!    endfor

%!test
%! func_R = {@euler321_to_rotation_matrix, @euler123_to_rotation_matrix, @rotation_vector_to_rotation_matrix};
%! func_Phi = {@rotation_matrix_to_euler321, @rotation_matrix_to_euler123, @rotation_matrix_to_rotation_vector};
%! for k=1:numel(func_R)
%! for i = 1:100
%!    phi = 0.5 * pi * ( 2 * rand() - 1 );
%!    n = ( 2 * rand(3,1) - 1 );
%!    n /= norm(n);
%!    Phi = n * phi;
%!    R = feval(func_R{k}, Phi);
%!    assert(R * R.',eye(3),sqrt(eps));
%!    assert(R.' * R,eye(3),sqrt(eps));
%!    Phi2 = feval(func_Phi{k}, R);
%!    assert(mod(pi + Phi2, 2 * pi) - pi, mod(pi + Phi, 2 * pi) - pi, sqrt(eps) * norm(Phi));
%! endfor
%! endfor
*/
