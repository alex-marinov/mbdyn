// $Header$
// MBDyn (C) is a multibody analysis code. 
// http://www.mbdyn.org
// 
// Copyright (C) 1996-2013
// 
// Pierangelo Masarati	<masarati@aero.polimi.it>
// Paolo Mantegazza	<mantegazza@aero.polimi.it>
// 
// Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
// via La Masa, 34 - 20156 Milano, Italy
// http://www.aero.polimi.it
// 
// Changing this copyright notice is forbidden.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation (version 2 of the License).
// 
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// Author: Tommaso Solcia <tommaso.solcia@mail.polimi.it>

function [x,y,typ] = SOCKREAD(job,arg1,arg2)

x=[];y=[];typ=[];
select job

case 'plot' then
    standard_draw(arg1)
case 'getinputs' then
    x = []; y = []; typ = [];
case 'getoutputs' then
    [x,y,typ] = standard_outputs(arg1);
case 'getorigin' then
  [x,y]=standard_origin(arg1);
case 'set' then
  x=arg1;
  graphics = arg1.graphics;
  exprs = graphics.exprs;
  model = arg1.model;
  while %t do
    [ok, hostname, port, N_CH, exprs] = getvalue('Set MBDyn parameters', ..
             ['Hostname';'PORT';'N channels'], ..
             list('str', 1, 'vec', 1, 'vec', 1), ..
             exprs);
    if ~ok then break
    end
    model.ipar = [port; N_CH];
    model.out = ones(N_CH, 1);
    model.opar = list(int8(asciimat(hostname)));
    graphics.exprs = exprs;
    x.graphics = graphics;
    x.model = model;
    break
  end

case 'define' then
  port = 10011
  N_CH = 2
  hostname = "127.0.0.1"
  model = scicos_model();
  model.sim = list('sockread',4);
  model.out = ones(N_CH, 1);
  model.ipar = [port; N_CH];
  model.opar = list(int8(asciimat(hostname)))
  model.blocktype = 'c';
  model.dep_ut = [%F; %T];
  exprs = [hostname; string([port; N_CH])];
  gr_i=	['x=orig(1),y=orig(2),w=sz(1),h=sz(2)';
            'txt=[''MBDyn'';''READ'']';
            'xstringb(x+0.1*w, y+0.1*h, txt, 0.7*w, 0.75*h, ''fill'')';
        ];
  x=standard_define([3 3],model,exprs,gr_i);
end
endfunction
