#
#MBDyn (C) is a multibody analysis code.
#http://www.mbdyn.org
#
#Copyright (C) 1996-2017
#
#Pierangelo Masarati	<pierangelo.masarati@polimi.it>
#Paolo Mantegazza	<paolo.mantegazza@polimi.it>
#
#Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
#via La Masa, 34 - 20156 Milano, Italy
#http://www.aero.polimi.it
#
#Changing this copyright notice is forbidden.
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation (version 2 of the License).
#
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#COPYRIGHT (C) 2016
#
#Marco Morandini <marco.morandini@polimi.it>
#Mattia Alioli   <mattia.alioli@polimi.it>
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.


from __future__ import print_function, division
#import subprocess
import re
import sys
if sys.version_info[0] < 3:
        import __builtin__ as builtins
else:
        import builtins
from MBDynLib import *

script, filename = sys.argv
preprocessing_include = re.compile("[\s]*include:")

nodes = []
bodies = []
joints = []
shells = []
beams = []

class PreProc:
    def __init__(self):
        self.lin = None
        self.lins = []
        self.preprocessing = False
        self.preprocessing_start_token = re.compile("#beginpreprocess")
        self.preprocessing_end_token = re.compile("#endpreprocess")
        self.proprocessingdirectivefile = None
        self.proprocessingdirectiveline = None

    def print_start(self, filename, linnumb):
        assert not self.preprocessing, (
                    'Error, #beginpreprocess: found '
                    'at ' + filename + ':' + str(linnumb) +
                    ' while inside a'
                    ' #preprocess directive\n'
                    'issued at ' + self.proprocessingdirectivefile + ':' + str(self.proprocessingdirectiveline))
        self.proprocessingdirectivefile = filename
        self.proprocessingdirectiveline = linnumb
        self.lins = []        
        print('#python #beginpreprocess at ' + filename + ':' + str(linnumb) + '\n')
            
    def print_end(self, preproc, filename, linnumb):
        assert preproc, (
                 'Error, #endpreprocess  found '
                 'at ' + filename + ':' + str(linnumb) +
                 ' without previous'
                 ' #beginpreprocess directive\n')
        self.proprocessingdirectivefile = None
        self.proprocessingdirectiveline = None
        self.lins = '\n'.join(self.lins)
        print('#python #endpreprocess at ' + filename + ':' + str(linnumb) + '\n')

    def print_p(self, preproc, linx):
        self.lin = linx
        if preproc:
             print('#python ' + str(len(self.lins) + 1) + ' # ' + self.lin.replace('\n', ' ').replace('\r', ''))
             self.lins.append(self.lin.replace('\n', ' ').replace('\r', ''))
        else:
            print(self.lin.replace('\n', ' ').replace('\r', ''))

    def compile(self, linx, filename, linnumb, glb):
        if self.preprocessing_start_token.match(linx):
            self.print_start(filename, linnumb)
            self.preprocessing = True
        elif self.preprocessing_end_token.match(linx):
            self.print_end(self.preprocessing, filename, linnumb)
            self.preprocessing = False
            exec(self.lins, glb) #FIXME
        else:
            self.print_p(self.preprocessing, linx)

def PreprocessMBDynFile(pp, fn):
    with open(fn) as fp:
        lineno = 1
        for ln in fp:
            if preprocessing_include.match(ln):
                assert not pp.preprocessing, (
                        'Error, include: found '
                        'at ' + fp.name + ':' + str(lineno) +
                        ' while inside a'
                        ' #preprocess directive\n'
                        'issued at ' + 
                        pp.proprocessingdirectivefile + ':' + 
                        str(pp.proprocessingdirectiveline))
                includedfilename = re.findall('"([^"]*)"', ln)[0]
                print('\n#')
                print('# include: \"' + includedfilename + '\"')
                print('#\n')
                PreprocessMBDynFile(pp, includedfilename)
            else:
                pp.compile(ln, fn, lineno, globals())
            lineno += 1

pp = PreProc()
PreprocessMBDynFile(pp, filename)


print('\n')
print('# vim:ft=mbd')
print('\n')




