## MBDyn (C) is a multibody analysis code.
## http://www.mbdyn.org
##
## Copyright (C) 1996-2020
##
## Pierangelo Masarati	<masarati@aero.polimi.it>
## Paolo Mantegazza	<mantegazza@aero.polimi.it>
##
## Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
## via La Masa, 34 - 20156 Milano, Italy
## http://www.aero.polimi.it
##
## Changing this copyright notice is forbidden.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation (version 2 of the License).
##
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## AUTHOR: Reinhard Resch <octave-user@a1.net>
## Copyright (C) 2020(-2020) all rights reserved.

## The copyright of this code is transferred
## to Pierangelo Masarati and Paolo Mantegazza
## for use in the software MBDyn as described
## in the GNU Public License version 2.1

## Usage:
## MBDyn installation:
##   configure --enable-octave --enable-autodiff --with-static-modules
##
## Gmsh installation:
##   See http://www.gmsh.info/
##
## Octave package installation:
##   octave --eval 'pkg install -forge nurbs'
##   for pkg in mboct-octave-pkg mboct-numerical-pkg mboct-mbdyn-pkg mboct-fem-pkg; do
##     git clone https://github.com/octave-user/${pkg}.git
##     make -C ${pkg} install_local
##   done

clear all;
close all;
pkg load mboct-fem-pkg;

param.N1 = 4;
param.N2 = 100;
param.r2 = 1.5;
param.r1 = 0.1;
param.rc = 0;
param.rs = 0.01;
param.g = 9.81;
param.m1 = 1.2;
param.s1 = 1e4;
param.mu = 0.5;
param.sigma0 = 1e3;
param.sigma1 = 0;

fd = -1;
fname = [];

unwind_protect
  opts.mesh.promote_elem = {};
  opts.mesh.order = 1;
  opts.mesh.dim = 2;
  
  unwind_protect
    [fd, fname] = mkstemp(fullfile(tempdir(), "triangular_contact2_XXXXXX"));
    
    if (fd == -1)
      error("faild to create temporary file");
    endif
    
    fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
    fputs(fd, "Point(1) = {0, 0, 0};\n");
    fputs(fd, "Point(2) = {0, 0, -r2};\n");
    fputs(fd, "Point(3) = {r2, 0, 0};\n");
    fputs(fd, "Circle(1) = {2, 1, 3};\n");
    fputs(fd, "surf[] = Extrude{{0, 0, 1}, {0, 0, 0}, 2 * Pi}{ Curve{1}; };\n");
    fputs(fd, "Physical Surface(\"target\", 2) = {surf[]};\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
    fd = -1;
  end_unwind_protect

  opts.mesh.element_size = 2 * param.r2 * pi / param.N2;
  mesh_data(2).mesh= fem_pre_mesh_unstruct_create(fname, param, opts);

  unwind_protect
    fd = fopen(fname, "wt");
    
    if (fd == -1)
      error("faild to create temporary file");
    endif
    
    fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
    fputs(fd, "Point(1) = {r1, -r1, r1};\n");
    fputs(fd, "Point(2) = {-r1, -r1, r1};\n");
    fputs(fd, "Point(3) = {-r1, -r1, -r1};\n");
    fputs(fd, "Point(4) = {r1, -r1, -r1};\n");
    fputs(fd, "Line(1) = {1, 2};\n");
    fputs(fd, "Line(2) = {2, 3};\n");
    fputs(fd, "Line(3) = {3, 4};\n");
    fputs(fd, "Line(4) = {4, 1};\n");
    fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
    fputs(fd, "Plane Surface(6) = {5};\n");
    fputs(fd, "vol[] = Extrude{0, 2 * r1, 0}{ Surface{6}; };\n");
    fputs(fd, "Physical Surface(\"contact\", 1) = {vol[1]};\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
    fd = -1;
  end_unwind_protect

  opts.mesh.element_size = 2 * param.r1 * pi / param.N1;
  mesh_data(1).mesh = fem_pre_mesh_unstruct_create(fname, param, opts);

  mbdyn_pre_write_param_file([fname, ".set"], param);

  putenv("TRIANGULAR_CONTACT_SET", [fname, ".set"]);
  
  fd = -1;
 
  unwind_protect
    fd = fopen([fname, ".elm"], "wt");

    if (fd == -1)
      error("failed to open file \"%s.elm\"", fname);
    endif

    putenv("TRIANGULAR_CONTACT_ELEM", [fname, ".elm"]);
    
    fputs(fd, "user defined: elem_id_contact, triangular contact,\n");
    fputs(fd, "target node, node_id_target,\n");
    fputs(fd, "penalty function, \"penalty\",\n");
    fputs(fd, "search radius, rs,\n");
    
    if (param.mu > 0)
      fputs(fd, "friction model, lugre,\n");
      fputs(fd, "method, implicit euler,\n");
      fputs(fd, "coulomb friction coefficient, mu,\n");
      fputs(fd, "micro slip stiffness, sigma0,\n");
      fputs(fd, "micro slip damping, sigma1,\n");
    endif
    
    fprintf(fd, "number of target vertices, %d,\n", rows(mesh_data(2).mesh.nodes));
    fprintf(fd, "reference, global, %.16e, %.16e, %.16e,\n", mesh_data(2).mesh.nodes(:, 1:3).');
    fprintf(fd, "number of target faces, %d,\n", rows(mesh_data(2).mesh.elements.tria3));
    fprintf(fd, "%d, %d, %d,\n", mesh_data(2).mesh.elements.tria3(:, end:-1:1).');
    fprintf(fd, "number of contact nodes, %d,\n", 1);
    fprintf(fd, "node_id_contact, number of contact vertices, %d", rows(mesh_data(1).mesh.nodes));
    fprintf(fd, ",\noffset, reference, ref_id_contact, %.16e, %.16e, %.16e, radius, rc", mesh_data(1).mesh.nodes(:, 1:3).');
    fputs(fd, ";\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  options.mbdyn_command = "mbdyn";
  options.output_file = [fname, "_mbd"];

  for i=1:numel(mesh_data)
    mesh_data(i).dof_map.totdof = 3 * rows(mesh_data(i).mesh.nodes);
    mesh_data(i).dof_map.ndof = [reshape(1:mesh_data(i).dof_map.totdof, rows(mesh_data(i).mesh.nodes), 3), zeros(rows(mesh_data(i).mesh.nodes), 3)];
    mesh_data(i).mesh.material_data = struct("C",[],"rho",[])([]);
  endfor
  
  [mesh, dof_map] = fem_post_mesh_merge(mesh_data);

  mbdyn_solver_run("triangular_contact2.mbd", options);

  log_dat = mbdyn_post_load_log(options.output_file);

  r1 = log_dat.vars.r1;
  r2 = log_dat.vars.r2;
  m1 = log_dat.vars.m1;
  J1 = log_dat.vars.J1;
  g = log_dat.vars.g;

  [t, traj, def, vel, acc, node_id] = mbdyn_post_load_output_struct(options.output_file);
  [drive_id, drive_data] = mbdyn_post_load_output_drv(options.output_file);

  sol.def = zeros(rows(mesh.nodes), 6, numel(t));

  Rref = eye(3);
  Xref = zeros(3, 1);
  for i=1:2
    X0 = traj{i}(:, 1:3).';
    Phi = traj{i}(:, 4:6).';
    R0 = euler123_to_rotation_matrix(Phi);
    l0 = mesh_data(i).mesh.nodes(:, 1:3).';
    idx = (1:columns(l0)) + dof_map.submesh.offset.nodes(i);
    for j=1:3
      sol.def(idx, j, :) = repmat(X0(j, :) - Xref(j), numel(idx), 1);
      for k=1:3
	sol.def(idx, j, :) += reshape((R0(j, k, :) - Rref(j, k)) .* l0(k, :), ...
				      numel(idx), 1, columns(X0));
      endfor
    endfor
  endfor

  opt_post.elem_types = {"tria3"};
  opt_post.skin_only = false;

  fem_post_sol_external(mesh, sol, opt_post);
unwind_protect_cleanup
  if (numel(fname))
    fn = dir([fname, "*"]);
    for i=1:numel(fn)
      unlink(fullfile(fn(i).folder, fn(i).name));
    endfor
    unlink(fname);
  endif
end_unwind_protect
