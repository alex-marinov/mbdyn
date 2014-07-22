% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2014
% 
% Pierangelo Masarati   <masarati@aero.polimi.it>
% Paolo Mantegazza      <mantegazza@aero.polimi.it>
% 
% Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
% via La Masa, 34 - 20156 Milano, Italy
% http://www.aero.polimi.it
% 
% Changing this copyright notice is forbidden.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 2 of the License).
% 
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% Author: Giuseppe Quaranta <quaranta@aero.polimi.it>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_mls_interface([inputfilename])
%	inputfilename [optional argument] string that contains the name
% 	of the file created by the MBDyn external mapping force element 
%	using the "echo" option
%       Additional rows must be either added in the comment section of the file
%       or can be enforced using the appropriate options of the "echo" keyword
%	if input values different from default are required. In detail
%
%	# surface: ...
%		used to indicate the name of the file that contains 
%		the external element grid data. Default: 'surface.dat' 
%	
%	# order: used to indicate the polynomial order of the base used by MLS
%		The value must be in the range [1,3]. Default is 2
%
%	# basenode: number of nodes included in the local MLS support
%		shold be higher than 0
%
%	# weight: continuity order of the RBF weight function
%		The value must be in the range [-2,inf]
%		-1 means constant weight == 1 in the support
%		-2 means a special weight function == 0 @ the query point 
%               default 2
%
%	# output: name of the file where the interface matrix is saved in sparse format

% $Header$
function create_mls_interface(varargin)


if (nargin == 0)
	disp('Interactive Mode');
	filename = input('Input filename: ', 's');
else
	filename = varargin{1};
end

% safe defaults
extfile = 'surface.dat';
outputfile = 'output.dat';
labFlag = 0;
RefNodeId = -1;
MBNodeN = 0;
MBNodeList = []; 
IOrder = 2;
IBaseNodN = 0;
for i = 0:2,
    IBaseNodN = IBaseNodN + nchoosek(3 + i - 1, i);
end
IWeight = 3; 

fid = fopen(filename, 'r');
if (fid < 0)
	error(sprintf('unable to open file "%s"', filename));
end

tline = fgetl(fid);
while ~(feof(fid))
	[token,rem] = strtok(tline);
	if strcmp(token(1), '#') 
		[token,rem] = strtok(rem);
		switch token
	
			case 'surface:'
				[token,rem] = strtok(rem);
				extfile = token;

			case 'labels:'
				[token,rem] = strtok(rem);
				if strcmp(token, 'on')
					labFlag = 1;
				end

			case 'reference:'
				[token,rem] = strtok(rem);
				RefNodeId = str2double(token);
		
			case 'points:' 
				[token,rem] = strtok(rem);
				MBNodeN = str2double(token);
				MBNodeList = zeros(MBNodeN, 3);
			
			case 'order:'	
				[token,rem] = strtok(rem);
				IOrder = str2double(token);
				if (IOrder < 1) || (IOrder > 3)
					error('the polynomial base order must be between 1 and 3');
				end

			case 'basenode:'
				[token,rem] = strtok(rem);
				IBaseNodN = str2double(token);
				if (IBaseNodN <= 0)
					error('the number of nodes should be a positive number');
				end
			
			case 'weight:'
				[token,rem] = strtok(rem);
				IWeight = str2double(token);
				if (IWeight < -2) 	
					error('the weight order must be an integer between -2 and inf');
				end

			case 'output:'
				[outputfile,rem] = strtok(rem);

			otherwise 
				% ignore
		end
		tline = fgetl(fid);
	else
		if (MBNodeN == 0)
			error(sprintf('no node number found before the node position list in file "%s"', filename));
		end
		if labFlag  
			Lab = 	sscanf(tline, '%d');
			MBNodeList(1,:) = (sscanf(tline, '%g', 3))';
			for i = 2 : MBNodeN
				if feof(fid)
					error(sprintf('unexpected end of file "%s"', filename));
				end
				Lab = 	fscanf(fid, '%d');
				MBNodeList(i,:) = (fscanf(fid, '%g', 3))';
			end
		else
            MBNodeList(1,:) = (sscanf(tline, '%g', 3))';
			for i = 2 : MBNodeN
				if feof(fid)
					error(sprintf('unexpected end of file "%s"', filename));
				end                
				MBNodeList(i,:) = (fscanf(fid, '%g', 3))';
			end
		end
 		fclose(fid);
		break;
	end
end

% read extfile
if isempty(extfile)
	error('no external surface file is given');
end

fid = fopen(extfile, 'r');
if (fid < 0)
	error(sprintf('unable to open file "%s"', extfile));
end

ExtNodeN = [];
while ~(feof(fid))
	tline = fgetl(fid);
	[token, rem] = strtok(tline);
	if (not(isempty(token))),
		if (strcmp(token, '#')),
			% octave's save -text: '# rows: <nrows>'
			[token, rem] = strtok(rem);
			if (strcmp(token, 'rows:')),
				n = sscanf(rem, '%d');
				if (n > 0),
					ExtNodeN = n;
				end
			end
		else
			% legacy format: first non-comment row
			if (isempty(ExtNodeN)),
				ExtNodeN = sscanf(tline, '%d');

			else
				ExtNodeList = zeros(ExtNodeN, 3);
				ExtNodeList(1, :) = (sscanf(tline, '%g', 3))';
				for i = 2:ExtNodeN 
					if feof(fid)
						error(sprintf('unexpected end of file "%s"', extfile));
						fclose(fid);
						return;
					end
					ExtNodeList(i,:) = (fscanf(fid, '%g', 3))';
				end
				break;	
			end
		end
	end
end
fclose(fid);

I = f_regression_derivatives_kdtree(ExtNodeList,MBNodeList,[],IOrder,IBaseNodN,IWeight);
[r,c,s] = find(I);
fid = fopen(outputfile, 'w');

if (fid < 0)
	error(sprintf('unable to open file "%s"', outputfile));
	return;
end

fprintf(fid, '# Moving Least Square interface matrix\n');
fprintf(fid, '# Date %s\n',  date);
fprintf(fid, '# MBDyn Structural file "%s" nodes=%d\n', filename, MBNodeN);
fprintf(fid, '# External file "%s" nodes=%d\n', extfile, ExtNodeN);
if (labFlag == 1),
	fprintf(fid, '# labels: on\n');
else
	fprintf(fid, '# labels: off\n');
end
if (RefNodeId > -1),
	fprintf(fid, '# reference: %d\n', RefNodeId);
end
fprintf(fid, '# order: %d\n', IOrder);
fprintf(fid, '# basenode: %d\n', IBaseNodN);
if (IWeight == inf),
	fprintf(fid, '# weight: inf\n');
else
	fprintf(fid, '# weight: %d\n', IWeight);
end
fprintf(fid, '# output: %s\n', outputfile);
fprintf(fid, '# mapping matrix: %dx%d (%d non-zero coefficients, %f%%)\n', 3*ExtNodeN, 3*MBNodeN, 3*length(r), 100*3*length(r)/(3*ExtNodeN*3*MBNodeN));
for i = 1 : length(r)
	fprintf(fid, '%d %d %.16g\n', 3 * (r(i) - 1) + 1, 3 * (c(i) - 1) + 1, s(i));
	fprintf(fid, '%d %d %.16g\n', 3 * (r(i) - 1) + 2, 3 * (c(i) - 1) + 2, s(i));
	fprintf(fid, '%d %d %.16g\n', 3 * (r(i) - 1) + 3, 3 * (c(i) - 1) + 3, s(i));
end
fclose(fid);

disp('Interface computation completed successfully'); 
return;

