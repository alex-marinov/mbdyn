#
#MBDyn (C) is a multibody analysis code.
#http://www.mbdyn.org
#
#Copyright (C) 1996-2016
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
if sys.version_info[0] < 3:
        import __builtin__ as builtins
else:
        import builtins

declared_ConstMBVars = {}
declared_MBVars = {}

MBDynLib_simplfy = True

def get_value(x):
    if isinstance(x, expression):
        return x.__get__()
    else:
        return x

def simplify_null_element_multiplication(l, r):
    if MBDynLib_simplfy:
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
    if MBDynLib_simplfy:
        if l == 0:
            return True
        else:
            return False
    else:
        return False

def simplify_neutral_element(l, r, op, ne):
    if MBDynLib_simplfy:
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


class null:
    def __str__(self):
        s = 'null'
        return s

class eye:
    def __str__(self):
        s = 'eye'
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


class Node:
    def __init__(self, idx, pos, orient, vel, angular_vel, node_type = 'dynamic'):
        self.idx = idx
        self.position = pos
        self.orientation = orient
        self.velocity = vel
        self.angular_velocity = angular_vel
        self.node_type = node_type
    def __str__(self):
        s = 'structural: ' + str(self.idx) + ', ' + str(self.node_type) + ',\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity) + ';\n'
        return s

class DynamicNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'dynamic')

class StaticNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'static')

class Body:
    def __init__(self, idx, node, mass, position, inertial_matrix, inertial = null):
        self.idx = idx
        self.node = node
        self.mass = mass
        self.position = position
        self.inertial_matrix = inertial_matrix
        self.inertial = inertial
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
        s = s + ';\n'
        return s

class Clamp:
    def __init__(self, idx, node, pos = Position('', 'node'), orient = Position('', 'node')):
        self.idx = idx
        self.node = node
        self.position = pos
        self.orientation = orient
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', clamp, ' + str(self.node) + ',\n'
        s = s + '\tposition, ' + str(self.position) + ',\n'
        s = s + '\torientation, ' + str(self.orientation) + ';\n'
        return s

class Shell:
    def __init__(self, shell_type, idx, nodes, const_law):
        self.shell_type = shell_type
        self.idx = idx
        self.nodes = nodes
        if isinstance(const_law, list):
            self.const_law = const_law
        else:
            self.const_law = [const_law]
    def __str__(self):
        s = str(self.shell_type) + ': ' + str(self.idx) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.nodes) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.const_law) + ';\n'
        return s
        
class Beam:
    def __init__(self, idx, nodes, positions, orientations, const_laws_orientations, const_laws):
        assert len(nodes) == 3 or len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative orientations;\n' +
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
        s = s + ';\n'
        return s
