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
## Octave package installation:
##   octave --eval 'pkg install -forge nurbs'
##   for pkg in mboct-octave-pkg mboct-numerical-pkg mboct-mbdyn-pkg; do
##     git clone https://github.com/octave-user/${pkg}.git
##     make -C ${pkg} install_local
##   done

pkg load mboct-mbdyn-pkg;

clear all;
close all;

param.N = 1000;
param.g = 9.81;
param.h2 = 1.5;
param.r2 = 2.5;
param.r1 = 0.25;
param.m1 = 1.2;
param.s1 = 1e6;
param.mu = 0.2;
param.n = 5;
param.sigma0 = 10;
param.sigma1 = 1;
param.alpha = 0.9;

options.output_file = "";

unwind_protect
  phi = (0:(2 * pi / param.N):(2 * pi * (param.N - 1) / param.N))(:);

  n=[zeros(1, 3);
     param.r2 * cos(phi), param.r2 * sin(phi), repmat(param.h2, numel(phi), 1)];

  v1 = 2:param.N + 1;
  v2 = [3:param.N + 1, 2];

  v=[v1(:), v2(:), ones(numel(v1), 1)];

  options.output_file = tempname();

  mbdyn_pre_write_param_file([options.output_file, ".set"], param);

  putenv("TRIANGULAR_CONTACT_SET", [options.output_file, ".set"]);

  fd = -1;

  unwind_protect
    fd = fopen([options.output_file, ".elm"], "wt");

    if (fd == -1)
      error("failed to open file \"%s.elm\"", options.output_file);
    endif

    putenv("TRIANGULAR_CONTACT_ELEM", [options.output_file, ".elm"]);

    fputs(fd, "user defined: elem_id_contact, triangular surface contact, help,\n");
    fputs(fd, "target node, node_id_cone,\n");
    fputs(fd, "penalty function, \"penalty\",\n");
    fputs(fd, "search radius, r1,\n");
    if (param.mu > 0)
      fputs(fd, "friction model, lugre,\n");
      fputs(fd, "method, implicit euler,\n");
      fputs(fd, "coulomb friction coefficient, mu,\n");
      fputs(fd, "micro slip stiffness, sigma0,\n");
      fputs(fd, "micro slip damping, sigma1,\n");
    endif
    fprintf(fd, "number of target vertices, %d,\n", rows(n));
    fprintf(fd, "reference, global, %.16e, %.16e, %.16e,\n", n.');
    fprintf(fd, "number of target faces, %d,\n", rows(v));
    fprintf(fd, "%d, %d, %d,\n", v.');
    fputs(fd, "number of contact nodes, 1,\n");
    fputs(fd, "node_id_sphere, number of contact vertices, 1, offset, reference, ref_id_sphere, null, radius, r1;\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  options.mbdyn_command = "mbdyn";
  mbdyn_solver_run("triangular_contact.mbd", options);

  log_dat = mbdyn_post_load_log(options.output_file);

  N = log_dat.vars.N;
  r1 = log_dat.vars.r1;
  r2 = log_dat.vars.r2;
  h2 = log_dat.vars.h2;
  m1 = log_dat.vars.m1;
  J1 = log_dat.vars.J1;
  g = log_dat.vars.g;
  Phi = log_dat.vars.Phi;

  if (log_dat.vars.omega0 == 0)
    output_file = options.output_file;
  else
    output_file = [options.output_file, "_rel"];
    ext = {".log", ".drv", ".out"};
    for i=1:numel(ext)
      err = symlink([options.output_file, ext{i}], [output_file, ext{i}]);
      if (err ~= 0)
	error("failed to create symbolic link");
      endif
    endfor
    mbdyn_post_abs_to_rel(log_dat.vars.node_id_cone, options.output_file, output_file, false);
  endif

  [t, traj, def, vel, acc, node_id] = mbdyn_post_load_output_struct(output_file);
  [drive_id, drive_data] = mbdyn_post_load_output_drv(output_file);

  figure("visible", "off");
  plot3(traj{log_dat.vars.node_idx_sphere_center}(:,1), traj{log_dat.vars.node_idx_sphere_center}(:,2), traj{log_dat.vars.node_idx_sphere_center}(:,3));

  patch("vertices", n, "faces", v, "facevertexcdata", repmat([1, 0, 0], rows(n), 1), "facecolor", "interp");
  xlabel("x [m]");
  ylabel("y [m]");
  zlabel("z [m]");
  grid on;
  grid minor on;
  daspect(ones(1, 3));
  title("trajectory");

  figure("visible", "off");
  hold on;
  plot(traj{log_dat.vars.node_idx_sphere_center}(:, 1), traj{log_dat.vars.node_idx_sphere_center}(:, 3), "-;trajectory;1");
  plot([-r2, 0, r2], [h2, 0, h2], "-;surface;0");
  xlabel("x [m]");
  ylabel("z [m]");
  grid on;
  grid minor on;
  title("trajectory xz-plane");

  figure("visible", "off");
  hold on;
  plot(t, drive_data{1}, "-;Wkin;1");
  plot(t, drive_data{2}, "-;Wpot;2");
  plot(t, drive_data{1} + drive_data{2}, "-;Wkin+Wpot;0");
  grid on;
  grid minor on;
  xlabel("t [s]");
  ylabel("W [J]");
  title("energy");

  R = euler123_to_rotation_matrix([0; pi / 2 - Phi; 0]);

  figure("visible", "off");
  hold on;
  plot(t, vel{log_dat.vars.node_idx_sphere_center}(:, 5) * r1, "-;omega * r1;1");
  plot(t, vel{log_dat.vars.node_idx_sphere_center}(:, 1:3) * R(:, 1), "-;v;2");
  plot(t, m1 * g * cos(Phi) / (J1 / r1^2 + m1) * t, "-;v1;0");
  xlabel("t [s]");
  ylabel("v [m/s]");
  grid on;
  grid minor on;
  title("velocity");
  figure_list();
unwind_protect_cleanup
  if (~isempty(options.output_file))
    fn = dir([options.output_file, "*"]);

    for i=1:numel(fn)
      unlink(fullfile(fn(i).folder, fn(i).name));
    endfor
  endif
end_unwind_protect
