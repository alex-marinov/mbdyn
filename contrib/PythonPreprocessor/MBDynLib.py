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
import sys
from numbers import Number, Integral
import pdb

if sys.version_info[0] < 3:
        import __builtin__ as builtins
else:
        import builtins

declared_ConstMBVars = {}
declared_IfndefMBVars = {}
declared_MBVars = {}

MBDynLib_simplify = True

def get_value(x):
    if isinstance(x, expression):
        return x.__get__()
    else:
        return x

def simplify_null_element_multiplication(l, r):
    if MBDynLib_simplify:
        if l == 0 or r == 0:
            return True
        else:
            return False
    else:
        return False

def simplify_null_element_division(l, r):
    assert get_value(r) != 0, (
        'Error, division by zero: \'' + srt(l) + ' / ' + str(r) + 
        '\'\n')
    if MBDynLib_simplify:
        if l == 0:
            return True
        else:
            return False
    else:
        return False

def simplify_neutral_element(l, r, op, ne):
    if MBDynLib_simplify:
        #if get_value(l) == ne:
        if l == ne:
            return r
        #elif get_value(r) == ne:
        elif r == ne:
            return l
        else:
            return op(l, r)
    else:
        return op(l, r)

class expression:
    def __init__(self):
        pass
    def __neg__(self):
        return negative(self)
    def __add__(self, other):
            return simplify_neutral_element(self, other, addition, 0) #addition(self, other)
    def __sub__(self, other):
            return simplify_neutral_element(self, other, subtraction, 0) #subtraction(self, other)
    def __pow__(self, other):
            return power(self, other)
    def __mul__(self, other):
            if simplify_null_element_multiplication(self, other):
                return 0
            else:
                return simplify_neutral_element(self, other, multiplication, 1) #multiplication(self, other)
    def __truediv__(self, other):
            if simplify_null_element_division(self, other):
                return 0
            else:
                return division(self, other)
    def __radd__(self, other):
            return simplify_neutral_element(other, self, addition, 0) #addition(other, self)
    def __rsub__(self, other):
            return simplify_neutral_element(other, self, subtraction, 0) #subtraction(other, self)
    def __rmul__(self, other):
            if simplify_null_element_multiplication(other, self):
                return 0
            else:
                return simplify_neutral_element(other, self, multiplication, 1) #multiplication(other, self)
    def __rtruediv__(self, other):
            if simplify_null_element_division(other, self):
                return 0
            else:
                return division(other, self)

class negative(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return -get_value(self.left)
    def __str__(self):
        ls = str(self.left)
        if isinstance(self.left, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.right, expression):
            ls = '(' + str(self.left) +')'
        return '-' + ls

import math
class sin(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return math.sin(get_value(self.left))
    def __str__(self):
        ls = str(self.left)
        return 'sin(' + ls + ')'

class cos(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return math.cos(get_value(self.left))
    def __str__(self):
        ls = str(self.left)
        return 'cos(' + ls + ')'

class tan(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.tan(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'tan(' + ls + ')'

class asin(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.asin(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'asin(' + ls + ')'

class acos(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.acos(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'acos(' + ls + ')'

class sqrt(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.sqrt(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'sqrt(' + ls + ')'

class terminal_expression(expression):
    def __init__(self, value):
        expression.__init__(self)
        self.value = value
    def __get___(self):
        return self.value
    def __str__(self):
        return str(value)

class binary_expression(expression):
    def __init__(self, left, right):
        expression.__init__(self)
        self.left = left
        self.right = right
    def __trunc__(self):
        y = self.__get__()
        assert isinstance(y, int), (
                'Error, __trunc__  required for expression \n\'' + 
                str(self) + 
                '\'\nof type ' + str(type(y)) +
                ' \n')
        return y
    def __index__(self):
        return self.__trunc__()

class atan2(binary_expression):
	def __init__(self, left, right):
		binary_expression.__init__(self, left, right)
	def __get__(self):
		return math.atan2(get_value(self.left), get_value(self.right))
	def __str__(self):
		ls = str(self.left)
		rs = str(self.right)
		return 'atan2(' + ls + ', ' + rs + ')'

class addition(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) + get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        return ls + ' + ' + rs
            
class subtraction(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) - get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        return ls + ' - ' + rs
            
class multiplication(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) * get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, addition) or isinstance(self.left, subtraction):
            ls = '(' + ls + ')'
        if isinstance(self.right, addition) or isinstance(self.right, subtraction):
            rs = '(' + rs + ')'
        return ls + ' * ' + rs
            
class division(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) / get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, addition) or isinstance(self.left, subtraction):
            ls = '(' + ls + ')'
        if isinstance(self.right, terminal_expression) or isinstance(self.right, MBVar) or isinstance(self.right, power):
            pass
        elif isinstance(self.right, expression):
            rs = '(' + rs + ')'
        return ls + ' / ' + rs

class power(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return pow(get_value(self.left), get_value(self.right))
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.left, expression):
            ls = '(' + ls + ')'
        if isinstance(self.right, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.right, expression):
            rs = '(' + rs + ')'
        return ls + ' ^ ' + rs

class MBVar(terminal_expression):
    def __init__(self, name, var_type, expression):
        assert(name)
        self.name = name
        self.var_type = var_type
        self.expression = expression
        #self.do_declare = do_declare
        #if self.do_declare:
        assert (name in declared_ConstMBVars) == False, (
            '\n-------------------\nERROR:' + 
            ' re-defining an already declared const variable:\n\t' + 
            var_type + ' ' + name + 
            '\n-------------------\n')
        assert (name in declared_IfndefMBVars) == False, (
            '\n-------------------\nERROR:' + 
            ' re-defining an already declared ifndef variable:\n\t' + 
            var_type + ' ' + name + 
            '\n-------------------\n')
        self.declare()
    def __get__(self):
        return get_value(self.expression)
    def __trunc__(self):
        y = self.__get__()
        assert isinstance(y, int), (
                'Error, __trunc__  required for expression \n\'' + 
                str(self) + 
                '\'\nof type ' + str(type(y)) +
                ' \n')
        return self.expression.__trunc__()
    def __index__(self):
        return self.expression.__trunc__()
    def __str__(self):
        return str(self.name)
    def declare(self):
        if self.name in declared_MBVars:
            assert declared_MBVars[self.name].var_type == self.var_type, (
                '\n-------------------\nERROR:' + 
                ' re-defining an already declared variable of type ' + str(declared_MBVars[self.name].var_type) + '\n' + 
                'with different type ' + str(self.var_type) +
                '\n-------------------\n')
            print('set: ' + self.name + ' = ' + str(self.expression) + ';')
        else:
            declared_MBVars[self.name] = self
            print('set: ' + self.var_type + ' ' + self.name + ' = ' + str(self.expression) + ';')
        #globals()[self.name] = self    
        #__builtins__[self.name] = self    
        setattr(builtins, self.name, self)


class ConstMBVar(MBVar):
    def __init__(self, name, var_type, value):
        MBVar.__init__(self, name, 'const ' + var_type, value)
    def declare(self):
        #assert self.do_declare == True, (
        #    '\n-------------------\nERROR:' +
        #    ' declaring either temporary '
        #    'or already declared variable:\n\t' + 
        #    self.var_type + ' ' + self.name + 
        #    '\n-------------------\n')
        MBVar.declare(self)
        #self.do_declare = False
        declared_ConstMBVars[self.name] = self


class IfndefMBVar(MBVar):
    def __init__(self, name, var_type, value):
        if name in declared_MBVars:
            pass
        else:
            MBVar.__init__(self, name, 'const ' + var_type, value)

class null:
    def __str__(self):
        s = 'null'
        return s

class eye:
    def __str__(self):
        s = 'eye'
        return s

class Reference:
    def __init__(self, idx, pos, orient, vel, angvel):
        assert isinstance(pos, Position), (
            '\n-------------------\nERROR:' + 
            ' the position of a reference must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(orient, Position), (
            '\n-------------------\nERROR:' + 
            ' the orientation of a reference must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(vel, Position), (
            '\n-------------------\nERROR:' + 
            ' the velocity of a reference must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(angvel, Position), (
            '\n-------------------\nERROR:' + 
            ' the angulare velocity of a reference must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        self.idx = idx
        self.position = pos
        self.orientation = orient
        self.velocity = vel
        self.angular_velocity = angvel
    def __str__(self):
        s = 'reference: '
        s = s + str(self.idx) + ', \n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity) + ';\n'
        return s

class Position:
    def __init__(self, ref, rel_pos):
        self.reference = ref
        if isinstance(rel_pos, list):
            self.relative_position = rel_pos
        else:
            self.relative_position = [rel_pos]
    def __str__(self):
        s = ''
        if self.reference != '':
            s = 'reference, ' + str(self.reference) + ', '
        s = s + ', '.join(str(i) for i in self.relative_position)
        return s
    def isnull(self):
        return (self.reference == '') and isinstance(self.relative_position[0], null)
    def iseye(self):
        return (self.reference == '') and isinstance(self.relative_position[0], eye)

class Node:
    def __init__(self, idx, pos, orient, vel, angular_vel, node_type = 'dynamic',
            scale = 'default', output = 'yes'):
        assert isinstance(pos, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial position of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(orient, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial orientation of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(vel, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial velocity of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(angular_vel, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial angular velocity of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert node_type in ('dynamic', 'static',), (
            '\n-------------------\nERROR:' + 
            ' unrecognised or unsupported node type;' + 
            '\n-------------------\n')
        self.idx = idx
        self.position = pos
        self.orientation = orient
        self.velocity = vel
        self.angular_velocity = angular_vel
        self.node_type = node_type
        self.scale = scale
        self.output = output
    def __str__(self):
        s = 'structural: ' + str(self.idx) + ', ' + str(self.node_type) + ',\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity)
        if self.scale != 'default':
            s = s + ',\n\tscale, ' + str(self.scale)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DynamicNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'dynamic')

class StaticNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'static')

class DisplacementNode():
    def __init__(self, idx, pos, vel, node_type = 'dynamic',
            scale = 'default', output = 'yes'):
        self.idx = idx
        self.position = pos
        self.velocity = vel
        self.node_type = node_type
        self.scale = scale
        self.output = output
    def __str__(self):
        s = 'structural: ' + str(self.idx) + ', ' + str(self.node_type) + ' displacement,\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.velocity)
        if self.scale != 'default':
            s = s + ',\n\t scale, ' + str(self.scale)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DynamicDisplacementNode(DisplacementNode):
    def __init__(self, idx, pos, vel):
        DisplacementNode.__init__(self, idx, pos, vel, 'dynamic')

class StaticDisplacementNode(DisplacementNode):
    def __init__(self, idx, pos, vel):
        DisplacementNode.__init__(self, idx, pos, vel, 'static')

class PointMass:
    def __init__(self, idx, node, mass, output = 'yes'):
        self.idx = idx
        self.node = node
        self.mass = mass
        self.output = output
    def __str__(self):
        s = 'body: ' + str(self.idx) + ', ' + str(self.node) + ', ' + str(self.mass)
        if self.output != 'yes':
            s = s + ', output, ' + str(self.output)
        s = s + ';\n'
        return s

class Body:
    def __init__(self, idx, node, mass, position, inertial_matrix, inertial = null,
            output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' +
            ' in defining a body, the center of mass relative position ' + 
            ' mass must be an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(inertial_matrix, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a body, the inertial matrix' + 
            ' must be a list;' + 
            '\n-------------------\n')
        self.idx = idx
        self.node = node
        self.mass = mass
        self.position = position
        self.inertial_matrix = inertial_matrix
        self.inertial = inertial
        self.output = output
    def __str__(self):
        s = 'body: ' + str(self.idx) + ', ' + str(self.node) + ',\n'
        s = s + '\t' + str(self.mass) + ',\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.inertial_matrix) 
        if self.inertial != null:
            s = s + ',\n'
            if isinstance(self.inertial, list):
                s = s + ', '.join(str(i) for i in self.inertial_matrix)
            else:
                s = s + ', ' + self.inertial_matrix
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralForce:
    def __init__(self, idx, node, ftype, position, force_drive, 
            force_orientation = [], moment_orientation = [],
            moment_drive = [], output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural force, the relative arm must be' +
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert ftype in {'absolute', 'follower', 'total'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural force: ' + str(ftype) + 
            '\n-------------------\n')
        if ftype == 'total':
            assert isinstance(force_orientation, Position), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total force, the force orientation ' +
                ' must be an instance of the Position class;' + 
                '\n-------------------\n')
            assert isinstance(moment_orientation, Position), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total force, the moment orientation ' +
                ' must be an instance of the Position class;' + 
                '\n-------------------\n')
        self.idx = idx
        self.node = node
        self.ftype = ftype
        self.position = position
        self.force_drive = force_drive
        self.force_orientation = force_orientation
        self.moment_orientation = moment_orientation
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'force: ' + str(self.idx) + ', ' + self.ftype
        s = s + ',\n\t' + str(self.node)
        s = s + ',\n\t\tposition, ' + str(self.position)
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation)
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation)
            s = s + ',\n\t\tforce, ' + ', '.join(str(i) for i in self.force_drive)
            s = s + ',\n\t\tmoment, ' + ', '.join(str(i) for i in self.moment_drive)
        else: # ftype = { absolute|follower }
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.force_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralInternalForce:
    def __init__(self, idx, nodes, ftype, positions, force_drive, 
            force_orientation = [], moment_orientation = [],
            moment_drive = [], output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal force with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions) , (
            '\n-------------------\nERROR:' + 
            ' in defining a structural internal force all relative arms ' +
            ' must be instances of the Position class;' + 
            '\n-------------------\n')
        assert ftype in {'absolute', 'follower', 'total'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural internal force: ' + str(ftype) + 
            '\n-------------------\n')
        if ftype == 'total':
            assert all(isinstance(pos, Position) for pos in force_orientation), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total internal force all the ' +
                ' force orientations must be instances of the Position class;' + 
                '\n-------------------\n')
            assert all(isinstance(pos, Position) for pos in moment_orientation), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total internal force all the ' +
                ' moment orientations must be instances of the Position class;' + 
                '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.ftype = ftype
        self.positions = positions
        self.force_drive = force_drive
        self.force_orientation = force_orientation
        self.moment_orientation = moment_orientation
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'force: ' + str(self.idx) + ', ' + self.ftype + ' internal'
        s = s + ',\n\t' + str(self.nodes[0])
        s = s + ',\n\t\tposition, ' + str(self.positions[0])
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation[0])
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation[0])
        s = s + ',\n\t' + str(self.nodes[1])
        s = s + ',\n\t\tposition, ' + str(self.positions[1])
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation[1])
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation[1])
            s = s + ',\n\t\tforce, ' + ', '.join(str(i) for i in self.force_drive)
            s = s + ',\n\t\tmoment, ' + ', '.join(str(i) for i in self.moment_drive)
        else: # ftype = { absolute|follower }
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.force_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralCouple:
    def __init__(self, idx, node, ctype, position, moment_drive, output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural couple, the relative arm must be' +
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert ctype in {'absolute', 'follower'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural couple: ' + str(ctype) + 
            ';\n-------------------\n')
        self.idx = idx
        self.node = node
        self.ctype = ctype
        self.position = position
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'couple: ' + str(self.idx) + ', ' + self.ctype
        s = s + ',\n\t' + str(self.node)
        if len(self.position):
            s = s + ',\n\t\tposition, ' + str(self.position)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.moment_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralInternalCouple:
    def __init__(self, idx, nodes, ctype, positions, moment_drive, output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal couple with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(positions) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal couple with ' + str(len(positions)) + 
            ' relative positions (!= 2);' +
            '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural internal couple all the relative positions ' +
            ' must be instances of the Position class;' + 
            '\n-------------------\n')
        assert ctype in {'absolute', 'follower'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural internal couple: ' + str(ctype) + 
            '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.ctype = ctype
        self.positions = positions
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'couple: ' + str(self.idx) + ', ' + self.ctype + ' inernal'
        s = s + ',\n\t' + str(self.nodes[0])
        s = s + ',\n\t\tposition, ' + str(self.position[0])
        s = s + ',\n\t' + str(self.nodes[1])
        s = s + ',\n\t\tposition, ' + str(self.position[1])
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.moment_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class Clamp:
    def __init__(self, idx, node, pos = Position('', 'node'), 
            orient = Position('', 'node'), output = 'yes'):
        self.idx = idx
        self.node = node
        self.position = pos
        self.orientation = orient
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', clamp, ' + str(self.node) + ',\n'
        s = s + '\tposition, ' + str(self.position) + ',\n'
        s = s + '\torientation, ' + str(self.orientation)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class TotalJoint:
    def __init__(self, idx, nodes, positions, \
            position_orientations, rotation_orientations, \
            position_constraints, orientation_constraints, \
            position_drive, orientation_drive,
            output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' nodes must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative positions must be given in a list' + 
            '\n-------------------\n')    
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(position_orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative position orientations must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == len(position_orientations), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions_orientations)) + ' position orientations;\n' +
            '\n-------------------\n')
        assert isinstance(rotation_orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative rotation orientations must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == len(rotation_orientations), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(rotation_orientations)) + ' rotation orientations;\n' +
            '\n-------------------\n')
        assert isinstance(position_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' position constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(position_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' +
            str(len(position_constrains)) + ' position constraints;\n' +
            '\n-------------------\n')
        assert isinstance(orientation_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' orientation constraints must be given as a list;' + 
            '\n-------------------\n')    
        assert len(orientation_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' +
            str(len(orientation_constrains)) + ' orientation constraints;\n' +
            '\n-------------------\n')
        assert all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint all offsets must be instances of ' + 
            ' the class Position;\n' +
            '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.position_orientations = position_orientations
        self.rotation_orientations = rotation_orientations
        self.position_constraints = position_constraints
        self.orientation_constraints = orientation_constraints
        self.position_drive = position_drive
        self.orientation_drive = orientation_drive
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', total joint'
        for (node, pos, pos_or, rot_or) in zip(self.nodes, self.positions,
                self.position_orientations, self.rotation_orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s = s + ',\n\t\tposition, ' + str(pos)
            if not(pos_or.iseye()):
                s = s + ',\n\t\tposition orientation, ' + str(pos_or)
            if not(rot_or.iseye()):
                s = s + ',\n\t\trotation orientation, ' + str(rot_or)
        if sum(self.position_constraints):
            s = s + ',\n\tposition constraint, '\
                    + ', '.join(str(pc) for pc in self.position_constraints)
            if isinstance(self.position_drive, list):
                s = s + ',\n\t\t' + ', '.join(str(i) for i in self.position_drive)
            else:
                s = s + ',\n\t\t' + str(self.position_drive)
        if sum(self.orientation_constraints):
            s = s + ',\n\torientation constraint, '\
                    + ', '.join(str(oc) for oc in self.orientation_constraints)
            if isinstance(self.orientation_drive, list):
                s = s + ',\n\t\t' + ', '.join(str(i) for i in self.orientation_drive)
            else:
                s = s + ',\n\t\t', + str(self.orientation_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class TotalPinJoint:
    def __init__(self, idx, node, 
            positions, position_orientations, rotation_orientations, 
            position_constraints, orientation_constraints, 
            position_drive, orientation_drive,
            output = 'yes'):
        if not isinstance(positions, list):
            positions = [positions]
        assert (len(positions) in [1, 2]) and all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative positions must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        if not isinstance(position_orientations, list):
            position_orientations = [position_orientations]
        assert ((len(position_orientations) in [1, 2]) and all([isinstance(pos, Position) for pos in position_orientations])), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative position orientations must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        if not isinstance(rotation_orientations, list):
            rotation_orientations = [rotation_orientations]
        assert ((len(rotation_orientations) in [1, 2]) and all([isinstance(pos, Position) for pos in rotation_orientations])), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative rotation orientations must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        assert isinstance(position_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' position constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(position_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(position_constrains)) + 
            ' position constraints;' + '\n-------------------\n')
        assert isinstance(orientation_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' orientation constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(orientation_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(orientation_constrains)) + 
            ' orientation constraints;' + '\n-------------------\n')
        self.idx = idx
        self.node = node
        self.positions = positions
        self.position_orientations = position_orientations
        self.rotation_orientations = rotation_orientations
        self.position_constraints = position_constraints
        self.orientation_constraints = orientation_constraints
        self.position_drive = position_drive
        self.orientation_drive = orientation_drive
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', total pin joint'
        s = s + ',\n\t' + str(self.node)
        if not(self.positions[0].isnull()):
            s = s + ',\n\t\tposition, ' + str(self.positions[0])
        if not(self.position_orientations[0].iseye()):
            s = s + ',\n\t\tposition orientation, ' + str(self.position_orientations[0])
        if not(self.rotation_orientations[0].iseye()):
            s = s + ',\n\t\trotation orientation, ' + str(self.rotation_orientations[0])
        if len(self.positions) == 2 and not(self.positions[1].isnull()):
            s = s + ',\n\t# GROUND'
            s = s + '\n\t\tposition, ' + str(self.positions[1])
        if len(self.position_orientations) == 2 and not(self.position_orientations[1].iseye()):
            s = s + ',\n\t\tposition orientation, ' + str(self.position_orientations[1])
        if len(self.rotation_orientations) == 2 and not(self.rotation_orientations[1].iseye()):
            s = s + ',\n\t\trotation orientation, ' + str(self.rotation_orientations[1])
        if sum(self.position_constraints):
            s = s + ',\n\tposition constraint, '\
                    + ', '.join(str(pc) for pc in self.position_constraints)
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.position_drive)
        if sum(self.orientation_constraints):
            s = s + ',\n\torientation constraint, '\
                    + ', '.join(str(oc) for oc in self.orientation_constraints)
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.orientation_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class JointRegularization:
    def __init__(self, idx, coefficients):
        assert (isinstance(coefficients, list) and len(coefficients) >= 1) or (isinstance(coefficients, Number)), (
            '\n-------------------\nERROR:' + 
            ' joint regularization needs at least one' +
            ' coefficient ' + '\n-------------------\n')
        self.idx = idx
        self.coefficients = coefficients
    def __str__(self):
        s = 'joint regularization: ' + str(self.idx) + ", tikhonov"
        if isinstance(self.coefficients, list):
            s = s + 'list, ' + ', '.join(str(co) for co in self.coefficients)
        else:
            s = s + ', ' + str(self.coefficients)
        s = s + ';\n'
        return s


class Rod:
    def __init__(self, idx, nodes, positions, const_law, length = 'from nodes', 
            output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a rod with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a rod with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        if not isinstance(positions, list):
            positions = [positions]
        assert all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a rod all offsets must be instances of ' + 
            ' the class Position;\n' +
            '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.const_law = const_law
        self.length = length
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', rod'
        for (node, position) in zip(self.nodes, self.positions):
            s = s + ',\n\t' + str(node)
            if not(position.isnull()):
                s = s + ',\n\t\tposition, ' + str(position)
        s = s + ',\n\t' + str(self.length) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.const_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DeformableDiaplacement:
    def __init__(self, idx, nodes, positions, orientations, const_law, output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a deformable displacement joint, the' +
            ' nodes must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a deformable displacement joint with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a deformable displacement joint, the' +
            ' relative positions must be given in a list' + 
            '\n-------------------\n')    
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a deformable displacement joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a deformable displacement joint, the' +
            ' relative position orientations must be given in a list' + 
            '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.constitutive_law = const_law
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', deformable displacement'
        for (node, pos, orient) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s + s + ',\n\t\tposition, ' + str(pos)
            if not(pos_or.iseye()):
                s + s + ',\n\t\torientation, ' + str(orient)
        s = s + '\n\t'
        if isinstance(self.constitutive_law, str):
            s = s + self.constitutive_law
        else:
            s = s + ', '.join(str(i) for i in self.constitutive_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DeformableHinge:
    def __init__(self, idx, nodes, positions, orientations, const_law, output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a deformable hinge, the' +
            ' nodes must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a deformable hinge with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a displacement hinge, the' +
            ' relative positions must be given in a list' + 
            '\n-------------------\n')    
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a deformable hinge with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a deformable hinge, the' +
            ' relative position orientations must be given in a list' + 
            '\n-------------------\n')
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.constitutive_law = const_law
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', deformable hinge'
        for (node, pos, orient) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s + s + ',\n\t\tposition, ' + str(pos)
            if not(orient.iseye()):
                s + s + ',\n\t\torientation, ' + str(orient)
        s = s + ',\n\t'
        if isinstance(self.constitutive_law, str):
            s = s + self.constitutive_law
        else:
            s = s + ', '.join(str(i) for i in self.constitutive_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class Shell:
    def __init__(self, shell_type, idx, nodes, const_law, output = 'yes'):
        self.shell_type = shell_type
        self.idx = idx
        self.nodes = nodes
        if isinstance(const_law, list):
            self.const_law = const_law
        else:
            self.const_law = [const_law]
        self.output = output
    def __str__(self):
        s = str(self.shell_type) + ': ' + str(self.idx) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.nodes) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.const_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s
        
class Beam:
    def __init__(self, idx, nodes, positions, orientations, const_laws_orientations,
            const_laws, output = 'yes'):
        assert len(nodes) == 3 or len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert len(nodes) == len(orientations), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(orientations)) + ' relative orientations;\n' +
            '\n-------------------\n')
        assert len(const_laws_orientations) == len(const_laws), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(const_laws)) +
            ' coonstitutive laws and ' + str(len(const_laws_orientations)) + ' constitutive law orientations;' +
            '\n-------------------\n')
        if len(nodes) == 2:
            self.beam_type = 'beam2'
        else:
            self.beam_type = 'beam3'
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.const_laws_orientations = const_laws_orientations
        self.const_laws = const_laws
        self.output = output
    def __str__(self):
        s = str(self.beam_type) + ': ' + str(self.idx)
        for (node, position, orientation) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node) + ',\n\t\tposition, ' + str(position) + ',\n\t\torientation, ' + str(orientation)
        for (cl_or, cl) in zip(self.const_laws_orientations, self.const_laws):
            s = s + ',\n\t' + str(cl_or) + ',\n\t' 
            if isinstance(cl, str):
                s = s + cl 
            else: 
                s  = s + ', '.join(str(i) for i in cl)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class AerodynamicBody:
    def __init__(self, idx, node, 
            position, orientation, span,
            chord, aero_center, b_c_point, twist, integration_points,
            induced_velocity = [], tip_loss = [], control = [], 
            airfoil_data = [], unsteady = [], 
            jacobian = 'no', custom_output = [], output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the' + 
            ' relative surface offset must be an instance of the' + 
            ' Position class;' + '\n-------------------\n')
        assert isinstance(orientation, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' relative surface orientation must be an instance of the' 
            ' Position class;' + '\n-------------------\n')
        assert isinstance(span, (Number)) or (isinstance(span, MBVar) and (span.var_type in ('real', 'const real'))), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body, the' + 
            ' surface span must be numeric' + 
            '\n-------------------\n')
        assert (isinstance(integration_points, Integral) and (integration_points > 0)) \
                or (isinstance(integration_points, MBVar) and integration_points.var_type in ('integer', 'const integer')), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with ' + str(integration_points) +
            ' integration_points' + '\n-------------------\n')
        assert (induced_velocity == []) or isinstance(induced_velocity, (Integral, MBVar)), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' induced velocity elment tag must be an integer or MBVar;' 
            '\n-------------------\n')
        assert not(len(unsteady)) or ((len(unsteady) > 0)*'bielawa' == 'bielawa'), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with unrecognised unsteady flag'
            '\n-------------------\n')
        assert (jacobian in {'yes', 'no'}) or isinstance(jacobian, bool), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with unrecognised jacobian flag'
            '\n-------------------\n')
        self.idx = idx
        self.node = node
        self.position = position
        self.orientation = orientation
        self.span = span
        self.chord = chord
        self.aero_center = aero_center
        self.b_c_point = b_c_point
        self.twist = twist
        self.integration_points = integration_points
        self.induced_velocity = induced_velocity
        self.tip_loss = tip_loss
        self.control = control
        self.airfoil_data = airfoil_data
        self.unsteady = unsteady
        self.jacobian = jacobian
        self.custom_output = custom_output
        self.output = output
    def __str__(self):
        s = 'aerodynamic body: ' + str(self.idx)
        s = s + ',\n\t ' + str(self.node)
        if self.induced_velocity:
            s = s + ',\n\t\tinduced velocity ' + str(self.induced_velocity)
        s = s + ',\n\t\t' + str(self.position)
        s = s + ',\n\t\t' + str(self.orientation)
        s = s + ',\n\t\t' + str(self.span)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.chord)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.aero_center)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.b_c_point)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.twist)
        if len(self.tip_loss):
            s = s + ',\n\t\ttip loss, ' + ', '.join(str(i) for i in self.tip_loss)
        s = s + '\n\t\t' + str(self.integration_points)
        if len(self.control):
            s = s + ',\n\t\tcontrol, ' + ', '.join(str(i) for i in self.control)
        if len(self.airfoil_data):
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.airfoil_data)
        if len(self.unsteady):
            s = s + ',\n\t\tunsteady, ' + str(self.unsteady)
        if len(self.jacobian):
            s = s + ',\n\t\tjacobian, ' + str(self.jacobian)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        if len(self.custom_output):
            s = s + ',\n\tcustom output, ' + ', '.join(str(i) for i in self.custom_output)
        s = s + ';\n'
        return s


class AerodynamicBeam:
    def __init__(self, idx, beam, 
            positions, orientations,
            chord, aero_center, b_c_point, twist, integration_points, 
            induced_velocity = [], tip_loss = [], control = [], 
            airfoil_data = [], unsteady = [], 
            jacobian = 'no', custom_output = [], output = 'yes'):
        assert len(positions) in {2,3}, (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(len(positions)) +
            ' relative surface offsets (not in [2,3])' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions), (
            ' in defining an aerodynamic beam the' + 
            ' relative surface offsets must be instances of the' + 
            ' Position class;' + '\n-------------------\n')
        assert len(orientations) in {2,3}, (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(len(orientations)) +
            ' relative surface orientations (not in [2,3])' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in orientations), (
            ' in defining an aerodynamic beam the' + 
            ' relative surface orientations must be instances of the' + 
            ' Position class;' + '\n-------------------\n')
        assert len(positions) == len(orientations), (
            '\n-------------------\nERROR:' + 
            ' definining an aerodynamic beam with ' + str(len(positions)) + 
            ' relative surface offsets and ' + str(len(orientations)) + 
            ' relative surface orientations' + '\n-------------------\n')
        assert (isinstance(integration_points, Integral) and (integration_points > 0))\
                or (isinstance(integration_points, MBVar) and integration_points.var_type in ('integer', 'const integer')), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(integration_points) +
            ' integration_points' + '\n-------------------\n')
        assert (induced_velocity == []) or isinstance(induced_velocity, (Integral, MBVar)), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' induced velocity elment tag must be an integer or an MBVar;' 
            '\n-------------------\n')
        assert not(len(unsteady)) or ((len(unsteady) > 0)*'bielawa' == 'bielawa'), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with unrecognised unsteady flag'
            '\n-------------------\n')
        assert (jacobian in {'yes', 'no'}) or isinstance(jacobian, bool), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with unrecognised jacobian flag'
            '\n-------------------\n')
        self.idx = idx
        self.beam = beam
        self.positions = positions
        self.orientations = orientations
        self.chord = chord
        self.aero_center = aero_center
        self.b_c_point = b_c_point
        self.twist = twist
        self.integration_points = integration_points
        self.induced_velocity = induced_velocity
        self.tip_loss = tip_loss
        self.control = control
        self.airfoil_data = airfoil_data
        self.unsteady = unsteady
        self.jacobian = jacobian
        self.custom_output = custom_output
        self.output = output
    def __str__(self):
        s = 'aerodynamic beam' + str(len(self.positions)) + ': ' + str(self.idx)
        s = s + ',\n\t ' + str(self.beam)
        if self.induced_velocity:
            s = s + ',\n\t\tinduced velocity ' + str(self.induced_velocity)
        for (pos, ori) in zip(self.positions, self.orientations):
            s = s + ',\n\t\t' + str(pos)
            s = s + ',\n\t\t' + str(ori)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.chord)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.aero_center)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.b_c_point)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.twist)
        if len(self.tip_loss):
            s = s + ',\n\t\ttip loss, ' + ', '.join(str(i) for i in self.tip_loss)
        s = s + ',\n\t\t' + str(self.integration_points)
        if len(self.control):
            s = s + ',\n\t\tcontrol, ' + ', '.join(str(i) for i in self.control)
        if len(self.airfoil_data):
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.airfoil_data)
        if len(self.unsteady):
            s = s + ',\n\t\tunsteady, ' + str(self.unsteady)
        if self.jacobian == 'yes':
            s = s + ',\n\t\tjacobian, ' + str(self.jacobian)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        if len(self.custom_output):
            s = s + ',\n\tcustom output, ' + ', '.join(str(i) for i in self.custom_output)
        s = s + ';\n'
        return s
