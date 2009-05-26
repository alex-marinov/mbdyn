#!BPY

__author__ = "G. Douglas Baldwin, douglasbaldwin AT verizon.net"
__url__ = ["http://www.baldwintechnology.com"]
__version__ = "0.3.0"
__bpydoc__ = """\
Description:

This script and its companion script, mbdyn_gui.py, provide a design environment for
creating, running, and importing the results of an MBDyn multibody dynamic model.

Usage:

Must download MBDyn from http://www.aero.polimi.it/~mbdyn/, compile, and make the 
executable "mbdyn" accessible by Blender.  MBDyn input file is saved in a Blender Text file, 
and can be manually editied before pressing the Run button to execute MBDyn.

"""

# --------------------------------------------------------------------------
# Blender MBDyn
# Copyright (C) 2008, 2009 G. Douglas Baldwin - http://www.baldwintechnology.com
# --------------------------------------------------------------------------
# ***** BEGIN GPL LICENSE BLOCK *****
#
#    This file is part of Blender MBDyn.
#
#    Blender MBDyn is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Blender MBDyn is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
# ***** END GPL LICENCE BLOCK *****
# -------------------------------------------------------------------------- 

import Blender, bpy
from Blender import Draw, Object, Scene, Mathutils, Window, Modifier
from Blender.Mathutils import Vector
from copy import copy
from math import sqrt

class Common(object):

	def makecopy(self):
		newcopy = copy(self)
		try:
			newcopy.objects = copy(self.objects)
		except:
			pass
		newcopy.links = copy(self.links)
		for link in newcopy.links:
			if link:
				link.users += 1
		return newcopy

	def name_check(self, database, preserve=False):
		names = [entity.name for entity in database if entity != self]
		hold = self.name
		if self.name in names:
			if '.001' <= self.name[-4:] and self.name[-4:] <= '.999':
				self.name = self.name[:-4]
			if self.name in names:
				self.name = self.name+'.'+str(1).zfill(3)
			qty = 1
			while self.name in names:
				qty += 1
				self.name = self.name[:-4]+'.'+str(qty).zfill(3)
				if qty >=999: raise ValueError, name
			if preserve:
				for entity in database:
					if entity != self and entity.name == hold:
						hold = self.name
						entity.name = self.name
						self.name = hold

	def permit_deletion(self, linking, string, delete, single):
		if linking: base = 1
		else: base = 0
		if not self.users and not linking:
			string += [('Delete',	delete, 'Delete this item')]
		elif self.users == 1 and not linking:
			string += [(str(self.users+base)+' user(s)', single, 'Duplicate')]
		elif self.users:
			string += [(str(self.users+base)+' user(s)', single, 'Make single user')]

	def finalize(self, nval, delete, single, entities):
		self.name = nval.val
		if self not in entities:
			entities.append(self)
			self.name_check(entities, False)
		elif delete.val:
			for link in self.links:
				if link:
					link.users -= 1
			entities.remove(self)
		elif single.val:
			self.name_check(entities, True)
			new = self.makecopy()
			new.name_check(entities, False)
			new.users = 0
			entities.append(new)
			scn = Scene.GetCurrent()
			try:
				if self.objects:
					old_obs = [ob for ob in scn.objects.selected if ob in self.objects]
					scn.objects.selected = old_obs
					Object.Duplicate()
					for old_ob, new_ob in zip(old_obs, scn.objects.selected):
						for i, ob in enumerate(new.objects):
							if ob == old_ob:
								new.objects[i] = new_ob
			except:
				pass
			return new
		else:
			self.name_check(entities, True)
		return self

	def select(self, index, clas, clas_types=[], head=None, exclude=None):
		if clas != 'Object' and self.links[index]:
			self.links[index].users -= 1
		menu = Menu(self.database, [clas], clas_types, head, exclude)
		sel = -1
		while sel == -1:
			sel = Draw.PupTreeMenu(menu.string)
		result = menu.action(sel, linking=True)
		if result:
			if clas == 'Object':
				self.objects[index] = result
			else:
				self.links[index] = result
				result.users += 1

	def check_loop(self, head):
		if self == head:
			return True
		looped = False
		for link in self.links:
			if link and link.check_loop(head):
				looped = True
		return looped

	def rotationMatrix_write(self, rot, text, pad):
		text.write(
		pad+', '.join([str(rot[0][j]) for j in range(3)])+',\n'+
		pad+', '.join([str(rot[1][j]) for j in range(3)])+',\n'+
		pad+', '.join([str(rot[2][j]) for j in range(3)]))

	def persistantPupMenu(self, string):
		sel = -1
		while sel == -1:
			sel = Draw.PupMenu(string)

class MBDyn(Common):

	entity_classes = ['Element', 'Drive', 'Driver', 'Friction', 'Shape', 'Function', 'NS_Node', 'Constitutive', 'Matrix']

	def __init__(self):
		self.NS_Node = []
		self.Element = []
		self.Constitutive = []
		self.Drive = []
		self.Driver = []
		self.Shape = []
		self.Function = []
		self.Friction = []
		self.Matrix = []
		self.Frame = []
		self.filename = ''
		self.defaults()

	def defaults(self):
		self._integrator = 'initial value'
		self._t0 = 0.
		self._tN = 0.
		self._dt = 1.e-3
		self._maxI = 10
		self._tol = 1.e-6
		self._dTol = 2.
		self._dC = 1.e-3
		self._fps = 30
		self._IPOs = 1
		self._port = 5500
		self._hostname = '127.0.0.1'
		self._reinit = 0
		self._posixRT = 0
		self._change_filename = 1

	def modify(self):
		integrator = Draw.Create(self._integrator)
		t0 = Draw.Create(self._t0)
		tN = Draw.Create(self._tN)
		dt_1e3 = Draw.Create(self._dt * 1.e3)
		maxI = Draw.Create(self._maxI)
		tol_1e6 = Draw.Create(self._tol * 1.e6)
		dTol = Draw.Create(self._dTol)
		dC_1e3 = Draw.Create(self._dC * 1.e3)
		fps = Draw.Create(self._fps)
		IPOs = Draw.Create(self._IPOs)
		port = Draw.Create(self._port)
		hostname = Draw.Create(self._hostname)
		reinit = Draw.Create(self._reinit)
		posixRT = Draw.Create(self._posixRT)
		try:
			change_filename = Draw.Create(self._change_filename)
		except:
			change_filename = Draw.Create(0)
			

		string = [
		('', integrator, 0, 12, 'integrator'),
		('t0', t0, 0.0, 1.0e6, 'initial time'),
		('tN', tN, 0.0, 1.0e6, 'final time (tN==0. => Use Blender animation range)'),
		('dt*1.e3', dt_1e3, 1.0e-3, 1.0e3, 'time step'),
		('maxI', maxI, 1, 100, 'max iterations'),
		('tol*1.e6', tol_1e6, 1.0e-3, 1.0e3, 'tolerance'),
		('dTol', dTol, 1.0e-3, 1.0e3, 'deriatives tolerance'),
		('dC*1.e3', dC_1e3, 1.0e-3, 1.0e3, 'derivatives coefficient'),
		('fps', fps, 1, 1000, 'frames per second: used to format output'),
		('IPOs', IPOs, 'Create new IPs, else overwrite existing IPOs'),
		('port', port, 1024, 49151, 'First sequential socket port number'),
		('Host name: ',	hostname, 0, 30, 'Name or IP of Game Blender Console'),
		('PosixRT',	posixRT, 'When using File Driver in Game Blender, run in POSIX Real Time Mode (WARNING: controls system clock)'),
		('Change Filename', change_filename, 'Change MBDyn filename'),
		('defaults', reinit, 'Revert to default values')]
		Draw.PupBlock('Parameters', string)
		self._integrator = integrator.val
		self._t0 = t0.val
		self._tN = tN.val
		self._dt = dt_1e3.val * 1.e-3
		self._maxI = maxI.val
		self._tol = tol_1e6.val * 1.e-6
		self._dTol = dTol.val
		self._dC = dC_1e3.val * 1.e-3
		self._fps = fps.val
		self._IPOs = IPOs.val
		self._port = port.val
		self._hostname = hostname.val
		self._reinit = reinit.val
		self._posixRT = posixRT.val
		if self.filename:
			self._change_filename = change_filename.val
		else:
			self._change_filename = 1

	def write(self, text):

		rigids = []
		self.rigid_dict = {}
		for element in self.Element:
			if element.type == 'Rigid':
				if self.rigid_dict.has_key(element.objects[0]):
					print ('Model Error: Object '+element.objects[0].name+
						' is linked to more than one Rigid Joint')
					return
				self.rigid_dict[element.objects[0]] = element.objects[1]
				rigids.append(element.objects[0])
		nodes = set([])
		for clas in MBDyn.entity_classes:
			eval('self.'+clas+'.sort(key = lambda x: x.name)')
		self.structural_dynamic_nodes = set([])
		self.structural_static_nodes = set([])
		for clas in ['Element', 'Drive']:
			for entity in eval('self.'+clas):
				try:
					if entity.objects[0] not in rigids:
						nodes |= set([entity.objects[0]])
						if entity.type in Element.structural_dynamic:
							self.structural_dynamic_nodes |= set([entity.objects[0]])
						elif entity.type in Element.structural_static:
							self.structural_static_nodes |= set([entity.objects[0]])
					elif entity.objects[0] in rigids:
						ob = self.rigid_dict[entity.objects[0]]
						nodes |= set([ob])
						if entity.type in Element.structural_dynamic:
							self.structural_dynamic_nodes |= set([ob])
						elif entity.type in Element.structural_static:
							self.structural_static_nodes |= set([ob])
				except:
					pass
		nodes -= set([None])
		self.structural_static_nodes -= set([None])
		self.structural_dynamic_nodes -= set([None])
		self.structural_static_nodes -= self.structural_dynamic_nodes
		structural_node_count = len(self.structural_static_nodes | self.structural_dynamic_nodes)
		self.Node = list(nodes)
		self.Node.sort(key = lambda x: x.name)
		self.Frame.sort(key = lambda x: x.objects[0].name)
		frame_dict = {}
		for i, frame in enumerate(self.Frame):
			for ob in frame.objects[1:]:
				frame_dict[ob] = str(i)

		joint_count = 0
		force_count = 0
		rigid_body_count = 0
		air_properties = False
		gravity = False
		aerodynamic_element_count = 0
		rotor_count = 0
		genel_count = 0
		beam_count = 0
		for element in self.Element:
			cat = element.category()
			if cat == 'joint':  joint_count += 1
			elif cat == 'force': force_count += 1
			elif cat == 'rigid body': rigid_body_count += 1
			elif cat == 'air properties': air_properties = True
			elif cat == 'aerodynamic element': aerodynamic_element_count += 1
			elif cat == 'gravity': gravity = True
			elif cat == 'rotor': rotor_count += 1
			elif cat == 'genel': genel_count += 1
			elif cat == 'beam': beam_count += 1
		file_driver_count = 0
		streaming_file_driver_count = 0
		for driver in self.Driver:
			driver.columns = []
		for drive in self.Drive:
			if drive.type == 'File drive' and drive.users:
				drive.links[0].columns.append(drive)
		for driver in self.Driver:
			if driver.columns:
				file_driver_count += 1
				if driver._args[0]:
					streaming_file_driver_count += 1
		electric_node_count = 0
		abstract_node_count = 0
		hydraulic_node_count = 0
		for ns_node in self.NS_Node:
			if ns_node.type == 'Electric': electric_node_count += 1
			elif ns_node.type == 'Abstract': abstract_node_count += 1
			elif ns_node.type == 'Hydraulic': hydraulic_node_count += 1

		finalTime = self._tN
		if finalTime == 0.:
			finalTime = self._t0 + float(Blender.Get('endframe') - Blender.Get('staframe')
				)/float(self._fps)

		text.write(
		'\n/* Label Indexes\n')
		if self.Frame:
			text.write('\nReference:\n')
			for i, frame in enumerate(self.Frame):
				text.write('\t'+str(i)+'\t- '+frame.objects[0].name+'\n')
		for clas in ['Node', 'NS_Node', 'Drive', 'Driver', 'Element']:
			if eval('self.'+clas):
				text.write('\n'+clas+'s:\n')
				for i, entity in enumerate(eval('self.'+clas)):
					text.write('\t'+str(i)+'\t- '+entity.name+' ('+entity.type)
					if entity.users:
						text.write(', users='+str(entity.users))
					text.write(')\n')
		if streaming_file_driver_count:
			text.write('\t'+str(len(self.Element))+'\t- Stream motion output (Output)\n')
		text.write('\n*/\n\n')

		text.write(
		'begin: data;\n'+
		'\tintegrator: '+self._integrator+';\n'+
		'end: data;\n\n'+
		'begin: '+self._integrator+';\n'+
		'\tinitial time: '+str(self._t0)+';\n')
		if streaming_file_driver_count:
			text.write('\tfinal time: forever;\n')
		else:
			text.write('\tfinal time: '+str(finalTime)+';\n')
		text.write(
		'\ttime step: '+str(self._dt)+';\n'+
		'\tmax iterations: '+str(self._maxI)+';\n'+
		'\ttolerance: '+str(self._tol)+';\n'+
		'\tderivatives tolerance: '+str(self._dTol)+';\n'+
		'\tderivatives coefficient: '+str(self._dC)+';\n')
		if self._posixRT and streaming_file_driver_count:
			text.write('\trealtime: POSIX, mode, period, time step, '+
			str(int(10e9*self._dt))+';\n')
#		text.write('\tlinear solver: naive, pivot factor, 1.e-9;\n') # for compatibility with mbdyn 1.3.7
		text.write(
		'end: '+str(self._integrator)+';\n\n'+
		'begin: control data;\n'+
		'\tdefault orientation: orientation matrix;\n'+
		'\toutput meter: meter, 0., forever, steps, '+
			str(int(round(1./(self._fps * self._dt))))+';\n')

		if structural_node_count:
			text.write('\tstructural nodes: '+str(structural_node_count)+';\n')
		if electric_node_count:
			text.write('\telectric nodes: '+str(electric_node_count)+';\n')
		if abstract_node_count:
			text.write('\tabstract nodes: '+str(abstract_node_count)+';\n')
		if hydraulic_node_count:
			text.write('\thydraulic nodes: '+str(hydraulic_node_count)+';\n')
		if joint_count:
			text.write('\tjoints: '+str(joint_count)+';\n')
		if force_count:
			text.write('\tforces: '+str(force_count)+';\n')
		if genel_count:
			text.write('\tgenels: '+str(genel_count)+';\n')
		if beam_count:
			text.write('\tbeams: '+str(beam_count)+';\n')
		if rigid_body_count:
			text.write('\trigid bodies: '+str(rigid_body_count)+';\n')
		if air_properties:
			text.write('\tair properties;\n')
		if gravity:
			text.write('\tgravity;\n')
		if aerodynamic_element_count:
			text.write('\taerodynamic elements: '+str(aerodynamic_element_count)+';\n')
		if rotor_count:
			text.write('\trotors: '+str(rotor_count)+';\n')
		if file_driver_count:
			text.write('\tfile drivers: '+str(file_driver_count)+';\n')
		if streaming_file_driver_count:
			text.write('\toutput elements: 1;\n'+
			'\tdefault output: none;\n')

		text.write('end: control data;\n')

		if self.Frame:
			text.write('\n')
		for i, frame in enumerate(self.Frame):
			label = str(i)
			parent_label = 'global'
			for ip, parent_frame in enumerate(self.Frame):
				if frame.objects[0] in parent_frame.objects[1:]:
					parent_label = str(ip)
					break
			rot = frame.objects[0].getMatrix().toQuat().toMatrix().transpose()
			if frame.links[0]._args[0]:
				localV = Mathutils.Vector(0., 0., 0.)
			else:
				localV = Mathutils.Vector([frame.links[0]._args[i] for i in range(1,4)])
				if frame.links[0]._args[4]:
					localV *= frame.links[0]._args[5]
			globalV = rot*localV
			if frame.links[1]._args[0]:
				localO = Mathutils.Vector(0., 0., 0.)
			else:
				localO = Mathutils.Vector([frame.links[1]._args[i] for i in range(1,4)])
				if frame.links[1]._args[4]:
					localO *= frame.links[1]._args[5]
			globalO = rot*localO

			text.write('reference: '+label+',\n'+
			'\treference, global, '+
			str(frame.objects[0].LocX)+', '+str(frame.objects[0].LocY)+', '+
				str(frame.objects[0].LocZ)+',\n'+
				'\treference, global, matr,\n')

			self.rotationMatrix_write(rot, text, '\t\t')

			text.write(',\n\treference, '+parent_label+', '+
			str(globalV.x)+', '+str(globalV.y)+', '+str(globalV.z)+',\n'+
			'\treference, '+parent_label+', '+
			str(globalO.x)+', '+str(globalO.y)+', '+str(globalO.z)+';\n')

		text.write('\nbegin: nodes;\n')

		for i, node in enumerate(self.Node):
			if node in self.structural_dynamic_nodes:
				text.write('\tstructural: '+str(i)+', dynamic,\n')
			elif node in self.structural_static_nodes:
				text.write('\tstructural: '+str(i)+', static,\n')
			else:
				print 'Program Error: Node '+node.name+' not in any node category'
				break
			rot = node.getMatrix().toQuat().toMatrix().transpose()
			if frame_dict.has_key(node):
				frame_label = frame_dict[node]
			else:
				frame_label = 'global'

			text.write(
				'\t\treference, global, '+
				str(node.LocX)+', '+str(node.LocY)+', '+str(node.LocZ)+',\n'+
				'\t\treference, global, matr,\n')
			self.rotationMatrix_write(rot, text, '\t'*3)
			text.write(',\n'+
				'\t\treference, '+frame_label+', null,\n'+
				'\t\treference, '+frame_label+', null;\n'
				)

		for i, ns_node in enumerate(self.NS_Node):
			if ns_node.type == 'Electric':
				text.write('\telectric: '+str(i)+', value, '+str(ns_node._args[0]))
				if ns_node._args[1]: text.write(', derivative, '+str(ns_node._args[2]))
				text.write(';\n')
			if ns_node.type == 'Abstract':
				text.write('\tabstract: '+str(i)+', value, '+str(ns_node._args[0]))
				if ns_node._args[1]: text.write(', differential, '+str(ns_node._args[2]))
				text.write(';\n')
			if ns_node.type == 'Hydraulic':
				text.write('\thydraulic: '+str(i)+', value, '+str(ns_node._args[0])+';\n')

		text.write('end: nodes;\n')

		if file_driver_count:
			text.write('\nbegin: drivers;\n')
			for driver in self.Driver:
				if driver.users:
					driver.write(text)
			text.write('end: drivers;\n')

		if self.Function:
			text.write('\n')
		for function in self.Function:
			function.written = False
		for function in self.Function:
			function.write(text)

		text.write('\nbegin: elements;\n')
		for element_type in [
			'Body',
			'Beam',
			'Rotor',
			'Aerodynamic',
			'Joint',
			'Force',
			'GENEL',
			'Air properties',
			'Gravity',
			'Driven']:
			for element in self.Element:
				if element.type == element_type:
					element.write(text)
		if streaming_file_driver_count:
			text.write(
			'\tstream motion output: '+str(len(self.Element))+',\n'+
			'\t\tstream name, "MAILBX",\n'+
			'\t\tcreate, no,\n'+
			'\t\t\tport, '+str(self._port+len(self.Driver))+',\n'
			'\t\t\thost, "'+self._hostname+'",\n'+
			'\t\tnon blocking,\n'+
			'\t\toutput every, 1,\n'+
			'\t\toutput flags, position, orientation matrix,\n')
			string = ''
			for i in range(len(self.Node)):
				string += str(i)+','
			text.write('\t\t'+string[:-1]+';\n')
		text.write('end: elements;\n')
		if streaming_file_driver_count:
			return 1

	def duplicate(self):
		obs = Object.GetSelected()
		elements = [element for element in self.Element if 
			element.objects and element.objects[0] in obs]
		for element in self.Element:
			if element.type == 'Driven' and element.links[1] in elements:
				elements.append(element)
		new_elements = []
		for element in elements:
			new_element = element.makecopy() 
			new_element.name_check(self.Element, False)
			new_element.users = 0
			new_elements.append(new_element)
			self.Element.append(new_element)
		new_obs = []
		scn = Scene.GetCurrent()
		scn.objects.selected = []
		for ob in obs:
			ob.sel = 1
			Object.Duplicate()
			new_obs += Object.GetSelected()
			scn.objects.selected = []
		scn.objects.selected = new_obs
		for new_element in new_elements:
			for i, ob in enumerate(new_element.objects):
				if ob in obs:
					new_element.objects[i] = new_obs[obs.index(ob)]
			for i, link in enumerate(new_element.links):
				if link in elements:
					link.users -= 1
					new_element.links[i] = new_elements[elements.index(link)]
					new_element.links[i].users += 1
		frames = [frame for frame in self.Frame if frame.objects[0] in obs]
		new_frames = []
		for frame in frames:
			new_frame = frame.makecopy() 
			new_frame.name_check(self.Frame, False)
			new_frames.append(new_frame)
			self.Frame.append(new_frame)
		for new_frame in new_frames:
			for i, ob in enumerate(new_frame.objects):
				if ob in obs:
					new_frame.objects[i] = new_obs[obs.index(ob)]
		Window.Redraw()
					
class Menu(Common):

	def __init__(self, database, classes=MBDyn.entity_classes+['Frame'], clas_types=[],
			head=None, exclude=None):
		self.database = database
		self.classes = classes
		self.clas_types = clas_types
		self.head = head
		self._event = 0
		self.events = {-1:'no menu pick'}
		self.eventsDict = {'no menu pick': -1}
		self._groups = {}
		self.string = []
		self._selected = Object.GetSelected()

		scn = Scene.GetCurrent()
		for ob in scn.objects:
			if ob.name in ['Edit', 'Make']:
				ret = Draw.PupMenu('Warning: Click here to rename object '+ob.name+' to be __'+ob.name)
				if not ret:
					return
				ob.name = '__'+ob.name

		for ob in self._selected:
			if ob.type != 'Mesh':
				Draw.PupMenu('Error: Must select only Blender Mesh objects. |'+
					'Your selection included a '+ob.type+' object.')
				return

		if clas_types:
			for typ in clas_types:
				try:
					subtypes = eval(self.classes[0]+'.'+typ)
				except:
					subtypes = []
				new = []
				for subtype in subtypes:
					new += [(subtype, [])]
				if new:
					self._assignevents(new)
					self._groups['New '+typ] = self._event
				hold = [(entity.name, []) for entity in eval('self.database.'+self.classes[0]) if
					entity.type == typ]
				hold.sort()
				self.string += [(typ, [('New', new)]+hold)]
				self._assignevents(self.string)
				self._groups[typ] = self._event
			return

		for clas in self.classes:

			if clas == 'Constitutive':
				types = eval(clas+'.types')
				clas_menu = []
				for typ in types:
					typ_menu = [('New', [(subtype,[]) for subtype in eval(clas+'.'+typ)])]
					self._assignevents(typ_menu)
					self._groups['New '+typ] = self._event
					hold = [(entity.name, []) for entity in self.database.Constitutive if entity.type == typ]
					hold.sort()
					typ_menu += hold
					self._assignevents(typ_menu)
					self._groups[typ] = self._event
					clas_menu += [(typ, typ_menu)]
				self.string += [(clas, clas_menu)]
				self._assignevents(self.string)
				self._groups[clas] = self._event

			elif clas == 'Matrix':
				matrices = []
				for typ in Matrix.types:
					hold = [(entity.name, []) for entity in self.database.Matrix if entity.type == typ]
					hold.sort()
					matrices += [(typ, [('New', [])]+hold)]
					self._assignevents(matrices)
					self._groups[typ] = self._event
				self.string += [('Matrix', matrices)]
				self._groups['Matrix'] = self._event

			elif clas == 'Frame':
				if not self._selected:
					break
				selected = copy(self._selected)
				looped = Frame.check_loop(selected, self.database.Frame)
				parent = [frame for frame in self.database.Frame if
					frame.objects[0] == self._selected[0]]
				hold = []
				if parent and not looped:
					hold = [('Edit', [])]
				if not parent and not looped:
					hold = [('Make', [])]
				child = [(frame.objects[0].name, []) for frame in self.database.Frame if
					self._selected[0] in frame.objects[1:]]
				if hold or child:
					self.string += [('Frame', hold + child)]
				obs = []
				for element in self.database.Element:
					obs += element.objects
				if not parent and not child and self._selected[0] not in obs:
					self.string += [('Delete', [])]
				self._assignevents(self.string)
				self._groups[clas] = self._event

			elif clas == 'Object':
				self.string += [('Object',[(ob.name, []) for ob in scn.objects
					if ob.type == 'Mesh' and ob != head])]
				self._assignevents(self.string)
				self._groups['Object'] = self._event

			elif clas == 'All Elements':
				self.string += [('Elements',[(e.name, []) for e in self.database.Element
					if (e != head and e.type != exclude)])]
				self._assignevents(self.string)
				self._groups['All Elements'] = self._event

			elif clas == 'Element':
				types = eval(clas+'.types')
				new_menu = []
				for typ in types:
					try:
						new_menu += [(typ, [(subtype,[]) for subtype in eval(clas+'.'+typ)])]
					except:
						new_menu += [(typ,[])]
				rigids = [e.objects[0] for e in self.database.Element if e.type == 'Rigid']
				if self._selected and self._selected[0] not in rigids:
					if not len(self._selected) == 2 or self._selected[1] not in rigids:
						new_menu += [('Rigid', [])]
				self._assignevents(new_menu)
				self._groups['New '+clas] = self._event
				hold = []
				if self._selected:
					primary = [(element.name, []) for element in self.database.Element if 
						element.objects and self._selected[0] == element.objects[0] and 
						not element.check_loop(head) and element.type != exclude] 
				else:
					primary = [(element.name, []) for element in self.database.Element if
						not element.objects and element.type != exclude]
				linked = [(element.name, []) for element in self.database.Element if 
					element.objects and self._selected and
					self._selected[0] in element.objects[1:] and element.type != exclude]
				primary.sort()
				linked.sort()
				if primary: hold += [('Primary', primary)]
				if linked: hold += [('Linked', linked)]
				self.string += [(clas, [('New', new_menu)]+hold)]
				self._assignevents(self.string)
				self._groups[clas] = self._event

			else:
				types = eval(clas+'.types')
				new_menu = [(typ,[]) for typ in types]
				self._assignevents(new_menu)
				self._groups['New '+clas] = self._event
				hold = [(entity.name, []) for entity in eval('self.database.'+clas) if
					not entity.check_loop(head)]
				hold.sort()
				self.string += [(clas, [('New', new_menu)]+hold)]
				self._assignevents(self.string)
				self._groups[clas] = self._event

	def _assignevents(self, menu):
		for i, item in enumerate(menu):
			if type(item[1]) is list and item[1] != []:
				self._assignevents(item[1])
			elif type(item[1]) is list:
				self.events[self._event] = item[0]
				self.eventsDict[item[0]] = self._event
				menu[i] = (item[0], self._event)
				self._event += 1

	def action(self, sel, linking=False):
		scn = Scene.GetCurrent()
		for typ in self.clas_types:
			if sel < self._groups[typ]:
				if self.events[sel] == 'New':
					return eval(self.classes[0]+'(typ, self.database, linking)')
				for entity in eval('self.database.'+self.classes[0]):
					if entity.name == self.events[sel]:
						hold = entity.modify(linking)
						return hold

		if 'Object' in self.classes and sel < self._groups['Object']:
			for ob in scn.objects:
				if ob.name == self.events[sel]:
					return ob

		if 'All Elements' in self.classes and sel < self._groups['All Elements']:
			for e in self.database.Element:
				if e.name == self.events[sel]:
					return e

		for clas in [clas for clas in self.classes if clas not in ['Matrix', 'Frame']]:
			if self._groups.has_key('New '+clas) and sel < self._groups['New '+clas]:
				return eval(clas+'(self.events[sel], self.database, linking)')
			elif self._groups.has_key(clas) and sel < self._groups[clas]:
				for entity in eval('self.database.'+clas):
					if entity.name == self.events[sel]:
						return entity.modify(linking, self.head)

		if 'Constitutive' in self.classes:
			for typ in Constitutive.types:
				if self._groups.has_key('New '+typ) and sel < self._groups['New '+typ]:
					return eval(clas+'(self.events[sel], self.database, linking)')
				elif self._groups.has_key(typ) and sel < self._groups[typ]:
					for entity in eval('self.database.'+clas):
						if entity.type == typ and entity.name == self.events[sel]:
							return entity.modify(linking, self.head)

		if 'Matrix' in self.classes and sel < self._groups['Matrix']:
			for typ in Matrix.types:
				if sel < self._groups[typ]:
					if self.events[sel] == 'New':
						return Matrix(typ, self.database, linking)
					for entity in self.database.Matrix:
						if entity.name == self.events[sel]:
							return entity.modify(linking)

		if self._groups.has_key('Frame') and sel < self._groups['Frame']:
			if self.events[sel] == 'Make':
				return Frame(self._selected, self.database)
			elif self.events[sel] == 'Edit':
				for frame in self.database.Frame:
					if frame.objects[0] == self._selected[0]:
						return frame.modify(self._selected)
			elif self.events[sel] == 'Delete':
				scn.objects.unlink(self._selected[0])
				Window.Redraw()
				return
			else:
				for frame in self.database.Frame:
					if frame.objects[0].name == self.events[sel]:
						self._selected = [frame.objects[0]]
						return frame.modify(self._selected)

class Frame(Common):

	def __init__(self, selected, database):
		self.database = database
		self.objects = selected
		self.links = [None]*2
		self.modify(selected)

	def modify(self, selected):
		names = []
		toggles = []
		for i in range(2):
			if self.links[i]:
				names.append(self.links[i].name)
				toggles.append(Draw.Create(0))
			else:
				names.append('select vector')
				toggles.append(Draw.Create(1))
		release = (Draw.Create(0))
		string = [
		(names[0],	toggles[0], 'Select frame linear velocity vector'),
		(names[1],	toggles[1], 'Select frame angular velocity vector')]
		if self in self.database.Frame:
			string += [
		('Release',	release, 'Release children to Global Reference Frame')]
		else:
			self.database.Frame.append(self)
		if not Draw.PupBlock('Frame', string): return
		if self not in self.database.Frame:
			self.database.Frame.append(self)
		if toggles[0].val and not release.val:
			self.persistantPupMenu('Select frame linear velocity vector:')
			self.select(0, 'Matrix', ['3x1'])
		if toggles[1].val and not release.val:
			self.persistantPupMenu('Select frame angular velocity vector:')
			self.select(1, 'Matrix', ['3x1'])
		for ob in set(selected) - set(self.objects):
			for frame in self.database.Frame:
				if ob in frame.objects[1:]:
					frame.objects.remove(ob)
			self.objects.append(ob)
		for ob in self.objects[1:]:
			ob.sel = 1
		self.objects[0].sel = 1
		if release.val:
			self.database.Frame.remove(self)
		else:
			self.draw()

	@classmethod
	def check_loop(self, selected, frames):
		looped = False
		if selected[0] in selected[1:]: return True
		for frame in frames:
			if selected[0] in frame.objects[1:]:
				selected[0] = frame.objects[0]
				if Frame.check_loop(selected, frames):
					looped = True
		return looped

	def draw(self):
		editmode = Window.EditMode()    # are we in edit mode?  If so ...
		if editmode: Window.EditMode(0) # leave edit mode before getting the mesh

		obj = self.objects[0]

		saveData = obj.getData(mesh=True)
		exists = False
		for mesh in bpy.data.meshes:
			if mesh.name == '_'+obj.name:
				exists = True
				me = mesh
		if not exists:
			me = bpy.data.meshes.new('_'+obj.name)
		me.verts = None
		me.verts.extend([(1.,0.,0.),(0.,2.,0.),(-1.,0.,0.),(0.,-2.,0.),(0.,0.,3.)]) 
		me.faces.extend([(3,2,1,0),(0,1,4),(1,2,4),(2,3,4),(3,0,4)])
		for f in me.faces: f.smooth = 0
		for e in me.edges: e.crease = 255
			
		obj.link(me)	
		mods = obj.modifiers
		hasSubsurf = False
		for mod in mods:
			if mod.name == 'Subsurf': hasSubsurf = True
		if not hasSubsurf:
			mod = mods.append(Modifier.Type.SUBSURF)
			mod[Modifier.Settings.LEVELS] = 3     # set subsurf subdivision levels to 3
		obj.makeDisplayList()
		if editmode: Window.EditMode(1)  # optional, just being nice
		Window.Redraw()

class NS_Node(Common):

	types = [
	'Electric',
	'Abstract',
	'Hydraulic',
	'Parameter''']

	def __init__(self, nodeType, database, linking=False):
		self.type = nodeType
		self.name = nodeType
		self.database = database
		self.users = 0
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'Abstract':
			self._args[1] = 0
		self.name_check(self.database.NS_Node)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Electric':
			string = [
			('Name: ', nval, 0, 30),
			('Initial value:',		args[0], -9.9e10, 9.9e10),
			('Differential IV:',	args[1], -9.9e10, 9.9e10, 'Differential initial value')]
			if self in self.database.NS_Node:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == 'Abstract':
			string = [
			('Name: ', nval, 0, 30),
			('Initial value:',	args[0], -9.9e10, 9.9e10),
			('Differential',		args[1], 'Diffferential (or Algebraic)'),
			('Derivative IV:',	args[2], -9.9e10, 9.9e10, 'Derivative initial value')]
			if self in self.database.NS_Node:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		elif self.type == 'Hydraulic':
			string = [
			('Name: ', nval, 0, 30),
			('Initial value:',	args[0], -9.9e10, 9.9e10)]
			if self in self.database.NS_Node:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
#		elif self.type == 'Parameter':
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.NS_Node)

class Element(Common):

	types = [
	'Aerodynamic',
	'Beam',
	'Body',
	'Force',
	'GENEL',
	'Joint',
	'Rotor',
	'Air properties',
	'Gravity',
	'Driven']
#	'Output']

	Aerodynamic = [
	'Aerodynamic body',
	'Aerodynamic beam2',
	'Aerodynamic beam3']

	Beam = [
	'Beam segment',
	'3-node beam']

	Force = [
	'Abstract force',
	'Structural force',
	'Structural internal force',
	'Structural couple',
	'Structural internal couple']

	GENEL = [
	'Swashplate']

	Joint = [
	'Axial rotation',
	'Clamp',
	'Distance',
	'Deformable displacement joint',
	'Deformable hinge',
	'Deformable joint',
	'In line',
	'In plane',
	'Revolute hinge',
	'Rod',
	'Spherical hinge',
	'Total joint']

	rigid_body = [
	'Body']

	structural_static = ['Aerodynamic', 'Joint', 'Rotor', 'Beam', 'Force']

	structural_dynamic = rigid_body

	def category(self):
		if self.type == 'Joint':
			return 'joint'
		elif self.type in Element.rigid_body:
			return 'rigid body'
		elif self.type == 'Air properties':
			return 'air properties'
		elif self.type == 'Aerodynamic':
			return 'aerodynamic element'
		elif self.type == 'Beam':
			if self.subtype == 'Beam segment':
				for element in self.database.Element:
					if (element.type == 'Beam' and element.subtype == '3-node beam' and
						self in element.links):
							return ''
			return 'beam'
		elif self.type == 'Gravity':
			return 'gravity'
		elif self.type == 'Rotor':
			return 'rotor'
		elif self.type == 'Force':
			return 'force'
		elif self.type == 'GENEL':
			return 'genel'

	def __init__(self, elemType, database, linking=False):
		if elemType in Element.Joint:
			self.type = 'Joint'
			self.subtype = elemType
		elif elemType in Element.Force:
			self.type = 'Force'
			self.subtype = elemType
		elif elemType in Element.GENEL:
			self.type = 'GENEL'
			self.subtype = elemType
		elif elemType in Element.Beam:
			self.type = 'Beam'
			self.subtype = elemType
		elif elemType in Element.Aerodynamic:
			self.type = 'Aerodynamic'
			self.subtype = elemType
		else:
			self.type = elemType
		self.name = elemType
		self.database = database
		self.users = 0
		self.objects = []
		self.links = []
		if self.type == 'Rotor':
			self._args = [1]*2 + [1, 0, 0, 0, 0] + [1.]*2 + [0]*3 + [0.]*3 + [0]*2 + [1] + [1.]*4
			self.objects = [None]*3
			self.links = [None]
		elif self.type == 'Aerodynamic':
			if self.subtype == 'Aerodynamic body':
				self._args = [1, 0, 0, 1., 1, 1, 1, 1, 1]
				self.objects = [None]
				self.links = [None]*5
			elif self.subtype in ['Aerodynamic beam2', 'Aerodynamic beam3']:
				beam_subtype = 'Beam segment'
				self.objects = [None]*2
				if self.subtype == 'Aerodynamic beam3':
					beam_subtype = '3-node beam'
					self.objects = [None]*3
				self._args = [0, 0, 1, 1, 1, 1, 1]
				self.links = [None]*6
				self.persistantPupMenu('Select '+beam_subtype+':')
				candidates = []
				candidate_dict = {}
				for i, element in enumerate(self.database.Element):
					if (element.type == 'Beam' and element.subtype == beam_subtype and 
						element.objects[0] == Object.GetSelected()[0] and not element.users):
						candidates += [(element.name, i)]
						candidate_dict[i] = element
				if candidates:
					result = Draw.PupTreeMenu([(beam_subtype+':', candidates)])
					if result != -1:
						self.links[0] = candidate_dict[result]
					else:
						return
				else:
					Draw.PupMenu('Error:  No unused '+beam_subtype+'(s) assigned to active object.')
					return
				scn = Scene.GetCurrent()
				scn.objects.selected = []
		elif self.type == 'Body':
			self._args = [1, 1., 1]
			self.objects = [None]
			self.links = [None]
		elif self.type == 'Beam':
			if self.subtype == 'Beam segment':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			if self.subtype == '3-node beam':
				self._args = [1, 1]
				self.objects = [None]*3
				self.links = [None]*2
		elif self.type == 'Force':
			if self.subtype == 'Abstract force':
				self._args = [1, 1]
				self.links = [None]*2
			if self.subtype == 'Structural force':
				self._args = [1, 1, 1]
				self.objects = [None]
				self.links = [None]
			if self.subtype == 'Structural internal force':
				self._args = [1, 1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			if self.subtype == 'Structural couple':
				self._args = [1, 1, 1]
				self.objects = [None]
				self.links = [None]
			if self.subtype == 'Structural internal couple':
				self._args = [1, 1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
		elif self.type == 'GENEL':
			if self.subtype == 'Swashplate':
				self._args = [1, 0, 0., 0.]*3 + [1]*3 + [0] + [0.]*3
				self.links = [None]*6
		elif self.type == 'Joint':
			if self.subtype == 'Axial rotation':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Clamp':
				self._args = [1]
				self.objects = [None]
			elif self.subtype == 'Deformable displacement joint':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Deformable hinge':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Deformable joint':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Distance':
				self._args = [1, 1, 1, 0]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'In line':
				self._args = [1, 1]
				self.objects = [None]*2
			elif self.subtype == 'In plane':
				self._args = [1, 1]
				self.objects = [None]*2
			elif self.subtype == 'Revolute hinge':
				self._args = [1, 1, 0, 0., 0, 0., 0, 0., 0]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Rod':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Spherical hinge':
				self._args = [1, 1]
				self.objects = [None]*2
			elif self.subtype == 'Total joint':
				self._args = [1]*15
				self.objects = [None]*2
				self.links = [None]*6
		elif self.type == 'Gravity':
			for element in self.database.Element:
				if element.type == 'Gravity':
					Draw.PupMenu("Error: Gravity element named '"+element.name+"' is already defined.")
					return
			self._args = [1, 1]
			self.links = [None]*2
		elif self.type == 'Air properties':
			for element in self.database.Element:
				if element.type == 'Air properties':
					Draw.PupMenu("Error: Air properties element named '"+element.name+"' is already defined.")
					return
			self._args = [1, 0, 0., 0., 1, 1, 0, 0, 0, 0., 0]
			self.links = [None]*5
		elif self.type == 'Driven':
			self._args = [1,1]
			self.links = [None]*2
		elif self.type == 'Rigid':
			self._args = [1, 1]
			self.objects = [None]*2
		sel = Object.GetSelected()
		if sel and self.objects:
			if len(sel) == 2 and len(self.objects) >= 2:
				self.objects[:2] = sel
			else:
				self.objects[0] = sel[0]
		self.name_check(self.database.Element)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Rotor':
			obj_name = ['select node']*3
			link_name = ['select drive']
			for i,j,ob in zip([0,1,2],[0,1,10],self.objects):
				if ob:
					obj_name[i] = ob.name
					args[j] = Draw.Create(0)
			if self.links[0]: link_name[0] = self.links[0].name
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			(obj_name[0],		 	args[0], 'Change craft node'),
			(obj_name[1],		 	args[1], 'Change rotor node'),
			('none',				args[2], 'Inflow model (select only one)'),
			('uniform',				args[3], 'Inflow model (select only one)'),
			('glauert',				args[4], 'Inflow model (select only one)'),
			('mangler',				args[5], 'Inflow model (select only one)'),
			('dynamic inflow',		args[6], 'Inflow model (select only one)'),
			('Ref omega:',			args[7], 1., 9.9e10),
			('Ref radius:',			args[8], 1., 100.),
			('Ground node',			args[9], 'Enable use of ground node'),
			(obj_name[2],		 	args[10], 'Change ground node'),
			('Initial values',		args[11], 'Enable specification of vel (dyn inflow only)'),
			('Avg vel:',			args[12], 1., 1.e2, 'Nominal inflow'),
			('Cosine vel:',			args[13], -1.e2, 1.e2, 'Fore-aft inflow'),
			('Sine vel:',			args[14], -1.e2, 1.e2, 'Lateral inflow'),
			('Delay',				args[15], 'Enable memory factor (all but dyn inflow)'),
			(link_name[0],		 	args[16], 'Change memory factor drive'),
			('Max iterations:',		args[17], 1, 1e2),
			('Tolerance:',			args[18], 1., 1.e2, 'Difference in nominal induced velocity (m/s)'),
			('eta:',				args[19], 0., 1., 'Fraction of induced velocity used'),
			('Hover cf:',			args[20], 0., 2., 'Hover correction factor'),
			('Forward cf:',			args[21], 0., 2., 'Forward flight correction factor')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[2:7])
			self._args = [arg.val for arg in args]
			if not delete.val:
				title = ['Select craft node:', 'Select rotor node:']
				for i in range(2):
					if self._args[i] or not self.objects[i]:
						self.persistantPupMenu(title[i])
						self.select(i, 'Object')
				if self._args[10] or (self._args[9] and not self.objects[2]):
					self.persistantPupMenu('Select ground node:')
					self.select(2, 'Object')
				if self._args[16] or (self._args[15] and not self.links[0]):
					self.persistantPupMenu('Select memory factor drive:')
					self.select(0, 'Drive')
			for i in [0, 1, 10, 16]:
				self._args[i] = 0
		elif self.type == 'Aerodynamic':
			if self.subtype == 'Aerodynamic body':
				obj_name = ['select node']
				link_name = ['select element'] + ['select shape']*4
				if self.objects[0]:
					obj_name[0] = self.objects[0].name
					args[0] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				for i, shape in enumerate(self.links[1:5]):
					if shape:
						link_name[i+1] = shape.name
						args[i+4] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change aerodynamic body node'),
				('Rotor',				args[1], 'Enable a rotor element'),
				(link_name[0],		 	args[2], 'Change rotor element'),
				('Surface span:',		args[3], 0., 1.e2),
				(link_name[1],		 	args[4], 'Change chord shape'),
				(link_name[2],		 	args[5], 'Change aerodynamic center shape'),
				(link_name[3],		 	args[6], 'Change boundary condition points shape'),
				(link_name[4],		 	args[7], 'Change surface twist shape'),
				('Integration points:',	args[8], 1, 1e2)]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0]:
						self.persistantPupMenu('Select aerodynamic body node:')
						self.select(0, 'Object')
					if self._args[2] or (self._args[1] and not self.links[0]):
						self.persistantPupMenu('Select rotor element:')
						self.select(0, 'Element', ['Rotor'], head=head)
					title = ['Select chord shape:', 'Select aerodynamic center shape:', 
						'Select boundary condition points shape', 'Select surface twist shape']
					for i in range(4):
						if self._args[i+4]:
							self.persistantPupMenu(title[i])
							self.select(i+1, 'Shape')
				for i in [0, 2, 4, 5, 6, 7]:
					self._args[i] = 0
			elif self.subtype in ['Aerodynamic beam2', 'Aerodynamic beam3']:
				link_name = ['select element']*2 + ['select shape']*4
				for i, link in enumerate(self.links[1:6]):
					if link:
						link_name[i+1] = link.name
						args[i+1] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				('beam='+self.links[0].name),
				('Rotor',				args[0], 'Enable a rotor element'),
				(link_name[1],		 	args[1], 'Change rotor element'),
				(link_name[2],		 	args[2], 'Change chord shape'),
				(link_name[3],		 	args[3], 'Change aerodynamic center shape'),
				(link_name[4],		 	args[4], 'Change boundary condition points shape'),
				(link_name[5],		 	args[5], 'Change surface twist shape'),
				('Integration points:',	args[6], 1, 1e2)]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					scn = Scene.GetCurrent()
					rotmat = Mathutils.RotationMatrix(90., 4, 'y')
					for i, ob in enumerate(self.objects):
						if not ob:
							scn.objects.selected = [self.links[0].objects[i]]
							Object.Duplicate()
							self.objects[i] = scn.objects.selected[0]
							scn.objects.selected = [self.objects[i], self.links[0].objects[i]]
							for element in self.database.Element:
								if element.type == 'Rigid' and element.objects[0] == self.links[0].objects[i]:
									scn.objects.selected = [self.objects[i], element.objects[1]]
							self.objects[i].setMatrix(rotmat * self.objects[i].getMatrix("localspace"))
							self.objects[i].setLocation(self.objects[i].getLocation())
							self.objects[i].sel = 1
							Element('Rigid', self.database)
					if self._args[1] or (self._args[0] and not self.links[1]):
						self.persistantPupMenu('Select rotor element:')
						self.select(1, 'Element', ['Rotor'], head=head)
					title = ['Select chord shape:', 'Select aerodynamic center shape:', 
						'Select boundary condition points shape:', 'Select surface twist shape:']
					for i in range(4):
						if self._args[i+2]:
							self.persistantPupMenu(title[i])
							self.select(i+2, 'Shape')
				for i in [1, 2, 3, 4, 5]:
					self._args[i] = 0
				if self not in self.database.Element:
					self.links[0].users += 1
		elif self.type == 'Body':
			obj_name = ['select node']
			link_name = ['select matrix']
			if self.objects[0]:
				obj_name[0] = self.objects[0].name
				args[0] = Draw.Create(0)
			if self.links[0]:
				link_name[0] = self.links[0].name
				args[2] = Draw.Create(0)
			string = [
			('Name: ', nval, 0, 30),
			(obj_name[0],		 	args[0], 'Change body node'),
			('Mass:',				args[1], 0.0001, 9.9e10),
			(link_name[0],		 	args[2], 'Change inertia matrix')]
			if self in self.database.Element:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[0]:
					self.persistantPupMenu('Select body node:')
					self.select(0, 'Object')
				if self._args[2]:
					self.persistantPupMenu('Select inertia matrix:')
					self.select(0, 'Matrix', ['3x3'])
			for i in [0, 2]:
				self._args[i] = 0
		elif self.type == 'Beam':
			if self.subtype == 'Beam segment':
				obj_name = ['select node']*2
				link_name = ['select constitutive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change first end'),
				(obj_name[1],	args[1], 'Change second end'),
				(link_name[0],	args[2], 'Change constitutive law')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select first end:')
						self.select(0, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select second end:')
						self.select(1, 'Object', head=self.objects[0])
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select constitutive law:')
						self.select(0, 'Constitutive', ['Const_6D'])
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == '3-node beam':
				link_name = ['select 1st segment', 'select 2nd segment']
				for i, link in enumerate(self.links):
					if link:
						link_name[i] = link.name
						args[i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
# Probably best to prohibit edit, and should now eliminate args[0:2]
#				(link_name[0],	args[0], 'Change 1st segment: '+link_name[0]),
#				(link_name[1],	args[1], 'Change 2nd segment: '+link_name[1])]
				(link_name[0]),
				(link_name[1])]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.links[0]:
						self.persistantPupMenu('Select 1st segment:')
						if self.links[0]:
							self.links[0].users -= 1
							self.links[0] = None
						candidates = []
						candidate_dict = {}
						for i, element in enumerate(self.database.Element):
							if (element.type == 'Beam' and element.subtype == 'Beam segment' and 
								element.objects[0] == Object.GetSelected()[0] and not element.users):
								candidates += [(element.name, i)]
								candidate_dict[i] = element
						if candidates:
							result = Draw.PupTreeMenu([('1st beam segment:', candidates)])
							if result != -1:
								self.links[0] = candidate_dict[result]
								self.links[0].users += 1
								self.objects[1] = self.links[0].objects[1]
							else:
								return
						else:
							Draw.PupMenu('Error:  No Beam Segment(s) assigned to active object.')
							return
					if self._args[1] or not self.links[1] or self._args[0]:
						self.persistantPupMenu('Select 2nd segment:')
						if self.links[1]:
							self.links[1].users -= 1
							self.links[1] = None
						candidates = []
						candidate_dict = {}
						for i, element in enumerate(self.database.Element):
							if (element.type == 'Beam' and element.subtype == 'Beam segment' and 
								element.objects[0] == self.objects[1] and not element.users):
								candidates += [(element.name, i)]
								candidate_dict[i] = element
						if candidates:
							result = Draw.PupTreeMenu([('2nd beam segment:', candidates)])
							if result != -1:
								self.links[1] = candidate_dict[result]
								self.links[1].users += 1
							else:
								return
						else:
							Draw.PupMenu('Error:  No 2nd Beam Segment(s) contiguous with 1st Beam.')
							self.links[0].users -= 1
							self.links[0] = None
							return
				for i, link in enumerate(self.links):
					if link:
						self.objects[i+1] = link.objects[1]
				for i in [0, 1]:
					self._args[i] = 0
		elif self.type == 'Force':
			if self.subtype == 'Abstract force':
				link_name = ['select abstract node', 'select drive']
				for i in range(2):
					if self.links[i]:
						link_name[i] = self.links[i].name
						args[i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(link_name[0],	args[0], 'Change abstract node'),
				(link_name[1],	args[1], 'Change force drive')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.links[0]:
						self.persistantPupMenu('Select abstract node:')
						self.select(0, 'NS_Node', clas_types=['Abstract'])
					if self._args[1] or not self.links[1]:
						self.persistantPupMenu('Select force drive')
						self.select(1, 'Drive')
				for i in [0, 1]:
					self._args[i] = 0
			elif self.subtype == 'Structural force':
				obj_name = ['select node']
				link_name = ['select drive']
				if self.objects[0]:
					obj_name[0] = self.objects[0].name
					args[0] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[1] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change forced node'),
				(link_name[0],	args[1], 'Change force drive'),
				('Follower',	args[2], 'Force orientation follows its node, else fixed in global frame')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select forced node:')
						self.select(0, 'Object')
					if self._args[1] or not self.links[0]:
						self.persistantPupMenu('Select force drive')
						self.select(0, 'Drive')
				for i in [0, 1]:
					self._args[i] = 0
			elif self.subtype == 'Structural internal force':
				obj_name = ['select node']*2
				link_name = ['select drive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change force producing node'),
				(obj_name[1],	args[1], 'Change attached node'),
				(link_name[0],	args[2], 'Change force producing drive'),
				('Follower',	args[3], 'Force orientation follows its node, else fixed in global frame')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select force producing node:')
						self.select(0, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select attached node:')
						self.select(1, 'Object', head=self.objects[0])
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select force producing drive')
						self.select(0, 'Drive')
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Structural couple':
				obj_name = ['select node']
				link_name = ['select drive']
				if self.objects[0]:
					obj_name[0] = self.objects[0].name
					args[0] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[1] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change forced node'),
				(link_name[0],	args[1], 'Change force drive'),
				('Follower',	args[2], 'Force orientation follows its node, else fixed in global frame')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select forced node:')
						self.select(0, 'Object')
					if self._args[1] or not self.links[0]:
						self.persistantPupMenu('Select force drive')
						self.select(0, 'Drive')
				for i in [0, 1]:
					self._args[i] = 0
			elif self.subtype == 'Structural internal couple':
				obj_name = ['select node']*2
				link_name = ['select drive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change force producing node'),
				(obj_name[1],	args[1], 'Change attached node'),
				(link_name[0],	args[2], 'Change force producing drive'),
				('Follower',	args[3], 'Force orientation follows its node, else fixed in global frame')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select force producing node:')
						self.select(0, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select attached node:')
						self.select(1, 'Object', head=self.objects[0])
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select force producing drive')
						self.select(0, 'Drive')
				for i in [0, 1, 2]:
					self._args[i] = 0
		elif self.type == 'GENEL':
			if self.subtype == 'Swashplate':
				link_name = ['select abstract node']*6
				for i, link in enumerate(self.links):
					if link: link_name[i] = self.links[i].name
				string = [
				('Name: ', nval, 0, 30),
				(link_name[0],		 	args[0], 'Change collective node'),
				('collective lmts', 	args[1], 'Set collective limits'),
				('min',					args[2], -1.e2, 1.e2, 'Minimum collective'),
				('max',					args[3], -1.e2, 1.e2, 'Maximum collective'),
				(link_name[1],		 	args[4], 'Change fore/aft node'),
				('fore/aft limits', 	args[5], 'Set fore/aft limits'),
				('min',					args[6], -1.e2, 1.e2, 'Minimum fore/aft'),
				('max',					args[7], -1.e2, 1.e2, 'Maximum fore/aft'),
				(link_name[2],		 	args[8], 'Change lateral node'),
				('Lateral limits', 		args[9], 'Set lateral limits'),
				('min',					args[10], -1.e2, 1.e2, 'Minimum lateral'),
				('max',					args[11], -1.e2, 1.e2, 'Maximum lateral'),
				(link_name[3],		 	args[12], 'Change actuator 1 abstract node'),
				(link_name[4],		 	args[13], 'Change actuator 2 abstract node'),
				(link_name[5],		 	args[14], 'Change actuator 3 abstract node'),
				('Enable',				args[15], 'Use the following coefs/factors'),
				('Dynamic',				args[16], -1.e2, 1.e2, 'Dynamic coefficient'),
				('Cyclic',				args[17], -1.e2, 1.e2, 'Cyclic factor'),
				('Collective',			args[18], -1.e2, 1.e2, 'Collective factor')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					pup = ['Select collective node',
					'Select fore/aft node',
					'Select lateral node',
					'Select actuator 1',
					'Select actuator 2',
					'Select actuator 3']
					for i, j in enumerate([0, 4, 8, 12, 13, 14]):
						if self._args[j] or not self.links[i]:
							self.persistantPupMenu(pup[i])
							self.select(i, 'NS_Node', clas_types=['Abstract'])
				for j in [0, 4, 8, 12, 13, 14]:
					self._args[j] = 0
		elif self.type == 'Joint':
			if self.subtype == 'Axial rotation':
				obj_name = ['select node']*2
				link_name = ['select drive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],			 	args[0], 'Change hinge node'),
				(obj_name[1],			 	args[1], 'Change attached node'),
				(link_name[0],			 	args[2], 'Change angular velocity drive')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select hinge node:')
						self.select(0, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select attached node:')
						self.select(1, 'Object', head=self.objects[0])
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select angular velocity drive')
						self.select(0, 'Drive')
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Clamp':
				obj_name = ['select node']
				if self.objects[0]:
					obj_name[0] = self.objects[0].name
					args[0] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change clamp node')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select clamp node:']
					if self._args[0] == 1 or self.objects[0] == None:
						self.persistantPupMenu(title[0])
						self.select(0, 'Object')
				for i in [0]:
					self._args[i] = 0
			elif self.subtype == 'Deformable displacement joint':
				obj_name = ['select node']*2
				link_name = ['select constitutive law']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
						link_name[0] = self.links[0].name
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change joint node'),
				(obj_name[1],		 	args[1], 'Change attached node'),
				(link_name[0],		 	args[2], 'Change constitutive law')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select joint node:', 'Select attached node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select constitutive law:')
						self.select(0, 'Constitutive', ['Const_3D'])
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Deformable hinge':
				obj_name = ['select node']*2
				link_name = ['select constitutive law']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
						link_name[0] = self.links[0].name
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change hinge node'),
				(obj_name[1],		 	args[1], 'Change attached node'),
				(link_name[0],		 	args[2], 'Change constitutive law')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select hinge node:', 'Select attached node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select constitutive law:')
						self.select(0, 'Constitutive', ['Const_3D'])
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Deformable joint':
				obj_name = ['select node']*2
				link_name = ['select constitutive law']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
						link_name[0] = self.links[0].name
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change joint node'),
				(obj_name[1],		 	args[1], 'Change attached node'),
				(link_name[0],		 	args[2], 'Change constitutive law')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select joint node:', 'Select attached node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[2] or not self.links[0]:
						self.persistantPupMenu('Select constitutive law:')
						self.select(0, 'Constitutive', ['Const_6D'])
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Distance':
				obj_name = ['select node']*2
				link_name = ['select drive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0], 	args[0], 'Change first node'),
				(obj_name[1], 	args[1], 'Change second node'),
				('From nodes',	args[2], 'Constant distance from initial node positons'),
				(link_name[0], 	args[3], 'Change distance drive (active if "From nodes" is deselected)')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select first node:', 'Select second node:']
					for i in range(2):
						if self._args[i] or not self.objects[i]:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[3] or (not self.links[0] and not self._args[2]):
						self.persistantPupMenu('Select distance drive')
						self.select(0, 'Drive')
				for i in [0, 1, 3]:
					self._args[i] = 0
			elif self.subtype == 'In line':
				obj_name = ['select node']*2
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change line node'),
				(obj_name[1],		 	args[1], 'Change point node')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select line node:', 'Select point node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
				for i in [0, 1]:
					self._args[i] = 0
			elif self.subtype == 'In plane':
				obj_name = ['select node']*2
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change plane node'),
				(obj_name[1],		 	args[1], 'Change point node')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select plane node:', 'Select point node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
				for i in [0, 1]:
					self._args[i] = 0
			elif self.subtype == 'Revolute hinge':
				obj_name = ['select node']*2
				link_name = ['select model']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
						link_name[0] = self.links[0].name
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],		 	args[0], 'Change hinge node'),
				(obj_name[1],		 	args[1], 'Change attached node'),
				('Enable theta',		args[2], 'Enable use of initial theta'),
				('Initial theta:',		args[3], 0., 9.9e10),
				('Enable friction',		args[4], 'Enable friction model'),
				('Avg radius:',			args[5], 0., 9.9e10, 'Used in friction model'),
				('Enable preload',		args[6], 'Enable preload'),
				('Preload:',			args[7], 0., 9.9e10, 'Used in friction model'),
				(link_name[0],		 	args[8], 'Change friction model')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select hinge node:', 'Select attached node:']
					for i in range(2):
						if self._args[i] or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[8] or (self._args[4] and not self.links[0]):
						self.persistantPupMenu('Select friction model:')
						self.select(0, 'Friction')
				for i in [0, 1, 8]:
					self._args[i] = 0
			elif self.subtype == 'Rod':
				obj_name = ['select node']*2
				link_name = ['select constitutive']
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				if self.links[0]:
					link_name[0] = self.links[0].name
					args[2] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0], 	args[0], 'Change first node'),
				(obj_name[1], 	args[1], 'Change second node'),
				(link_name[0], 	args[2], 'Change rod constitutive')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select first node:', 'Select second node:']
					for i in range(2):
						if self._args[i] or not self.objects[i]:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
					if self._args[2]:
						self.persistantPupMenu('Select rod constitutive')
						self.select(0, 'Constitutive', ['Const_1D'])
				for i in [0, 1, 2]:
					self._args[i] = 0
			elif self.subtype == 'Spherical hinge':
				obj_name = ['select node']*2
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],	args[0], 'Change hinge node'),
				(obj_name[1],	args[1], 'Change attached node')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					title = ['Select hinge node:', 'Select attached node:']
					for i in range(2):
						if self._args[i] == 1 or self.objects[i] == None:
							self.persistantPupMenu(title[i])
							self.select(i, 'Object')
				for i in [0, 1]:
					self._args[i] = 0
			if self.subtype == 'Total joint':
				obj_name = ['select node']*2
				link_name = ['select drive']*6
				for i, node in enumerate(self.objects):
					if node:
						obj_name[i] = node.name
						args[i] = Draw.Create(0)
				for i, link in enumerate(self.links):
					if link:
						link_name[i] = link.name
						args[4+2*i] = Draw.Create(0)
				string = [
				('Name: ', nval, 0, 30),
				(obj_name[0],			 	args[0], 'Change total joint node'),
				(obj_name[1],			 	args[1], 'Change attached node'),
				('Displace/Rotate',			args[2], 'Displace (or rotate) first'),
				('Active', 					args[3], 'Displacement-X drive is active'),
				(link_name[0],			 	args[4], 'Change displacement-X drive'),
				('Active', 					args[5], 'Displacement-Y drive is active'),
				(link_name[1],			 	args[6], 'Change displacement-Y drive'),
				('Active', 					args[7], 'Displacement-Z drive is active'),
				(link_name[2],			 	args[8], 'Change displacement-Z drive'),
				('Active', 					args[9], 'Angular displacement-X drive is active'),
				(link_name[3],			 	args[10], 'Change angular displacement-X drive'),
				('Active', 					args[11], 'Angular displacement-Y drive is active'),
				(link_name[4],			 	args[12], 'Change angular displacement-Y drive'),
				('Active', 					args[13], 'Angular displacement-Z drive is active'),
				(link_name[5],			 	args[14], 'Change angular displacement-Z drive')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type+': '+self.subtype, string): return
				self._args = [arg.val for arg in args]
				if not delete.val:
					if self._args[0] or not self.objects[0]:
						self.persistantPupMenu('Select total joint node:')
						self.select(0, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select attached node:')
						self.select(1, 'Object', head=self.objects[0])
					select_drive = ['Select displacement-X drive',
									'Select displacement-Y drive',
									'Select displacement-Z drive',
									'Select angular displacement-X drive',
									'Select angular displacement-Y drive',
									'Select angular displacement-Z drive']
					for i in range(6):
						if self._args[3+2*i] and (self._args[4+2*i]or not self.links[i]):
							self.persistantPupMenu(select_drive[i])
							self.select(i, 'Drive')
				for i in [0, 1, 4, 6, 8, 10, 12, 14]:
					self._args[i] = 0
		elif self.type == 'Gravity':
			link_name = ['select vector', 'select drive']
			for i, link in enumerate(self.links):
				if link:
					link_name[i] = link.name
					args[i] = Draw.Create(0)
			string = [
			('Name: ', nval, 0, 30),
			(link_name[0],	args[0], 'Change gravity vector'),
			(link_name[1],	args[1], 'Change gravity drive')]
			if self in self.database.Element:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select gravity vector:')
					self.select(0, 'Matrix', ['3x1'])
				if self._args[1] or not self.links[1]:
					self.persistantPupMenu('Select gravity drive:')
					self.select(1, 'Drive')
			for i in [0, 1]:
				self._args[i] = 0
		elif self.type == 'Air properties':
			link_name = ['select vector', 'select drive']+['select vector']*2+['select drive']
			for i, j, link in zip(range(5), [4,5,7,8,10], self.links):
				if link:
					link_name[i] = link.name
					args[j] = Draw.Create(0)
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('SI',					args[0], 'Standard atmosphere model (select only one)'),
			('British',				args[1], 'Standard atmosphere model (select only one)'),
			('Temp dev:',			args[2], -1.e3, 1.e3, 'Temperature deviation'),
			('Ref Alt:',			args[3], -1.e4, 1.e6, 'Reference altitude'),
			(link_name[0],		 	args[4], 'Change air velocity vector'),
			(link_name[1],		 	args[5], 'Change air velocity drive'),
			('Enable gust',			args[6], 'Enable gust model'),
			(link_name[2],		 	args[7], 'Change front direction vector'),
			(link_name[3],		 	args[8], 'Change perturbation direction vector'),
			('Front speed:',		args[9], -1.e3, 1.e3, 'Magnitude of front velocity'),
			(link_name[4],		 	args[10], 'Drive representing front profile')]
				if self in self.database.Element:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:2])
			self._args = [arg.val for arg in args]
			if not delete.val:
				title = ['Select front direction vector:', 'Select perturbation direction vector']
				if self._args[4] or not self.links[0]:
					self.persistantPupMenu('Select velocity vector:')
					self.select(0, 'Matrix', ['3x1'])
				if self._args[5] or not self.links[1]:
					self.persistantPupMenu('Select velocity drive:')
					self.select(1, 'Drive')
				for i in range(2):
					if self._args[i+7] or (not self.links[i+2] and self._args[6]):
						self.persistantPupMenu(title[i])
						self.select(i+2, 'Matrix', ['3x1'])
				if self._args[10] or (not self.links[4] and self._args[6]):
					self.persistantPupMenu('Select front profile drive:')
					self.select(4, 'Drive')
			for i in [4, 5, 7, 8, 10]:
				self._args[i] = 0
		elif self.type == 'Driven':
			link_name = ['select drive', 'select element']
			for i in range(2):
				if self.links[i]:
					link_name[i] = self.links[i].name
					args[i] = Draw.Create(0)
			string = [
			('Name: ', nval, 0, 30),
			(link_name[0],		 	args[0], 'Change drive'),
			(link_name[1],		 	args[1], 'Change driven element')]
			if self in self.database.Element:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[0] == 1:
					self.persistantPupMenu('Select drive:')
					self.select(0, 'Drive')
				if self._args[1] == 1:
					self.persistantPupMenu('Select driven element:')
					if not self.links[1] and Object.GetSelected():
						self.select(1, 'Element', head=self, exclude='Rigid')
					else:
						self.select(1, 'All Elements', head=self, exclude='Rigid')
			for i in [0, 1]:
				self._args[i] = 0
		elif self.type == 'Rigid':
			obj_name = ['select node']*2
			for i, node in enumerate(self.objects):
				if node:
					obj_name[i] = node.name
					args[i] = Draw.Create(0)
			string = [
			('Name: ', nval, 0, 30),
			(obj_name[0],	args[0], 'Change first node'),
			(obj_name[1],	args[1], 'Change second node')]
			if self in self.database.Element:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				title = ['Select first node:', 'Select second node:']
				for i in range(2):
					if self._args[i] == 1 or self.objects[i] == None:
						self.persistantPupMenu(title[i])
						self.select(i, 'Object')
			for i in [0, 1]:
				self._args[i] = 0
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		result = self.finalize(nval, delete, single, self.database.Element)
		scn = Scene.GetCurrent()
		scn.objects.selected = self.objects
		if self.objects:
			self.objects[0].sel = 1
		self.draw()
		return result

	def write(self, text):
		if self.type == 'Rotor':
			text.write(
			'\trotor: '+str(self.database.Element.index(self))+', '+
			str(self.database.Node.index(self.objects[0]))+', '+
			str(self.database.Node.index(self.objects[1]))+', induced velocity, ')
			if self._args[2]:
				text.write('no;\n')
				return
			if self._args[3]:
				text.write('uniform,\n')
			elif self._args[4]:
				text.write('glauert,\n')
			elif self._args[5]:
				text.write('mangler,\n')
			elif self._args[6]:
				text.write('dynamic inflow,\n')
			text.write('\t\t'+str(self._args[7:9]).strip('[]'))
			if self._args[9]:
				text.write(',\n\t\tground, '+str(self.database.Node.index(self.objects[2])))
			if self._args[11] and self._args[6]:
				text.write(',\n\t\tinitial value, '+str(self._args[12:15]).strip('[]'))
			if self._args[15] and not self._args[6]:
				text.write(',\n\t\tdelay, '+self.links[0].string(self))
			text.write(',\n\t\tmax iterations, '+str(self._args[17])+
			',\n\t\ttolerance, '+str(self._args[18])+
			',\n\t\teta, '+str(self._args[19])+
			',\n\t\tcorrection, '+str(self._args[20:22]).strip('[]')+';\n')
		elif self.type == 'Aerodynamic':
			if self.subtype == 'Aerodynamic body':
				text.write('\taerodynamic body: '+str(self.database.Element.index(self))+',\n')
				self.write_node(text, 0, position=False, orientation=False)
				if self._args[1]:
					text.write('\t\t\trotor, '+str(self.database.Element.index(self.links[0]))+',\n')
				self.write_node(text, 0, node=False, p_label='', o_label='')
				text.write(',\n\t\t'+str(self._args[3])+',\n')
				for shape in self.links[1:4]:
					text.write(shape.string()+',\n')
				text.write(self.links[4].string(3.14159265358979323846/180.)+',\n')
				text.write('\t\t'+str(self._args[8])+';\n')
			elif self.subtype == 'Aerodynamic beam2' or self.subtype == 'Aerodynamic beam3':
				text.write('\t'+self.subtype+': '+str(self.database.Element.index(self))+', '+
				str(self.database.Element.index(self.links[0]))+',\n')
				if self._args[0]:
					text.write(', rotor, '+str(self.database.Element.index(self.links[1])))
				for i in range(len(self.links[0].objects)):
					self.write_node(text, i, node=False)
					text.write(',\n')
				for shape in self.links[2:5]:
					text.write(shape.string()+',\n')
				text.write(self.links[5].string(3.14159265358979323846/180.)+',\n')
				text.write('\t\t'+str(self._args[6])+';\n')
		elif self.type == 'Body':
			text.write('\tbody: '+str(self.database.Element.index(self))+',\n')
			self.write_node(text, 0, position=False, orientation=False)
			text.write('\t\t\t'+str(self._args[1])+',\n')
			self.write_node(text, 0, node=False, orientation=False, p_label='')
			text.write(', '+self.links[0].string())
			self.write_node(text, 0, node=False, position=False, o_label='inertial')
			text.write(';\n')
		elif self.type == 'Beam':
			if self.subtype == 'Beam segment':
				for element in self.database.Element:
					if (element.type == 'Beam' and element.subtype == '3-node beam' and
						self in element.links):
							return
				text.write('\tbeam2: '+str(self.database.Element.index(self))+',\n')
				for i in range(len(self.objects)):
					self.write_node(text, i, p_label='position', o_label='orientation')
					text.write(',\n')
				text.write('\t\tfrom nodes, '+self.links[0].string()+';\n')
			if self.subtype == '3-node beam':
				for i, link in enumerate(self.links):
					if link:
						self.objects[i+1] = link.objects[1]
				text.write('\tbeam3: '+str(self.database.Element.index(self))+',\n')
				for i in range(len(self.objects)):
					self.write_node(text, i, p_label='position', o_label='orientation')
					text.write(',\n')
				text.write('\t\tfrom nodes, '+self.links[0].links[0].string())
				text.write(',\n\t\tfrom nodes, '+self.links[1].links[0].string()+';\n')
		elif self.type == 'Force':
			if self.subtype == 'Abstract force':
				text.write('\tforce: '+str(self.database.Element.index(self))+', abstract, '+
				str(self.database.NS_Node.index(self.links[0]))+', abstract,\n\t\t'+
				self.links[1].string(self)+';\n')
			elif self.subtype == 'Structural force':
				rot_0, globalV_0, Node_0 = self.rigid_offset(0)
				rotT_0 = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				relative_dir = rot_0*rotT_0*Vector([0., 0., 1.])
				relative_arm_0 = rot_0*globalV_0
				string = '\tforce: '+str(self.database.Element.index(self))+', '
				if self._args[2]:
					string += 'follower'
					relative_dir = rot_0*rotT_0*Vector([0., 0., 1.])
				else:
					string += 'absolute'
					relative_dir = rotT_0*Vector([0., 0., 1.])
				text.write(string+
				',\n\t\t'+str(Node_0)+
				',\n\t\t\tposition, '+str(relative_arm_0[0])+', '+str(relative_arm_0[1])+', '+str(relative_arm_0[2])+
				',\n\t\t\t'+str(relative_dir[0])+', '+str(relative_dir[1])+', '+str(relative_dir[2])+
				',\n\t\t'+self.links[0].string(self)+';\n')
			elif self.subtype == 'Structural internal force':
				rot_0, globalV_0, Node_0 = self.rigid_offset(0)
				rotT_0 = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				relative_arm_0 = rot_0*globalV_0
				rot_1, globalV_1, Node_1 = self.rigid_offset(1)
				relative_arm_1 = rot_1*globalV_1
				string = '\tforce: '+str(self.database.Element.index(self))+', '
				if self._args[3]:
					string += 'follower internal'
					relative_dir = rot_0*rotT_0*Vector([0., 0., 1.])
				else:
					string += 'absolute internal'
					relative_dir = rotT_0*Vector([0., 0., 1.])
				text.write(string+
				',\n\t\t'+str(Node_0)+
				',\n\t\t\t'+str(relative_dir[0])+', '+str(relative_dir[1])+', '+str(relative_dir[2])+
				',\n\t\t\t'+str(relative_arm_0[0])+', '+str(relative_arm_0[1])+', '+str(relative_arm_0[2])+
				',\n\t\t'+str(Node_1)+
				',\n\t\t\t'+str(relative_arm_1[0])+', '+str(relative_arm_1[1])+', '+str(relative_arm_1[2])+
				',\n\t\t'+self.links[0].string(self)+';\n')
			elif self.subtype == 'Structural couple':
				rot_0, globalV_0, Node_0 = self.rigid_offset(0)
				rotT_0 = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				string = '\tcouple: '+str(self.database.Element.index(self))+', '
				if self._args[2]:
					string += 'follower'
					relative_dir = rot_0*rotT_0*Vector([0., 0., 1.])
				else:
					string += 'absolute'
					relative_dir = rotT_0*Vector([0., 0., 1.])
				text.write(string+
				',\n\t\t'+str(Node_0)+
				',\n\t\t\t'+str(relative_dir[0])+', '+str(relative_dir[1])+', '+str(relative_dir[2])+
				',\n\t\t'+self.links[0].string(self)+';\n')
			elif self.subtype == 'Structural internal couple':
				rot_0, globalV_0, Node_0 = self.rigid_offset(0)
				rotT_0 = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				rot_1, globalV_1, Node_1 = self.rigid_offset(1)
				string = '\tcouple: '+str(self.database.Element.index(self))+', '
				if self._args[3]:
					string += 'follower internal'
					relative_dir = rot_0*rotT_0*Vector([0., 0., 1.])
				else:
					string += 'absolute internal'
					relative_dir = rotT_0*Vector([0., 0., 1.])
				text.write(string+
				',\n\t\t'+str(Node_0)+
				',\n\t\t\t'+str(relative_dir[0])+', '+str(relative_dir[1])+', '+str(relative_dir[2])+
				',\n\t\t'+str(Node_1)+
				',\n\t\t'+self.links[0].string(self)+';\n')
		elif self.type == 'GENEL':
			if self.subtype == 'Swashplate':
				text.write(
				'\tgenel: '+str(self.database.Element.index(self))+', swashplate')
				for i in range(3):
					text.write(',\n\t\t'+str(self.database.NS_Node.index(self.links[i])))
					if self._args[4*i+1]:
						text.write(', limits, '+str(self._args[4*i+2])+', '+str(self._args[4*i+3]))
				text.write(',\n\t\t'+str(self.database.NS_Node.index(self.links[3]))+
					', '+str(self.database.NS_Node.index(self.links[4]))+
					', '+str(self.database.NS_Node.index(self.links[5])))
				if self._args[15]:
					text.write(',\n\t\t'+str(self._args[16])+', '+str(self._args[17])+', '+str(self._args[18]))
				text.write(';\n')
		elif self.type == 'Joint':
			if self.subtype == 'Axial rotation':
				self.write_hinge(text, 'axial rotation')
				text.write(',\n\t\t'+self.links[0].string(self)+';\n')
			elif self.subtype == 'Clamp':
				text.write(
				'\tjoint: '+str(self.database.Element.index(self))+', clamp,\n'+
				'\t\t'+str(self.database.Node.index(self.objects[0]))+', node, node;\n')
			elif self.subtype == 'Deformable displacement joint':
				self.write_hinge(text, 'deformable displacement joint')
				text.write(',\n\t\t'+self.links[0].string()+';\n')
			elif self.subtype == 'Deformable hinge':
				self.write_hinge(text, 'deformable hinge', V1=False, V2=False)
				text.write(',\n\t\t'+self.links[0].string()+';\n')
			elif self.subtype == 'Deformable joint':
				self.write_hinge(text, 'deformable joint')
				text.write(',\n\t\t'+self.links[0].string()+';\n')
			elif self.subtype == 'Distance':
				text.write('\tjoint: '+str(self.database.Element.index(self))+', distance,\n')
				for i in range(2):
					self.write_node(text, i, orientation=False, p_label='position')
					text.write(',\n')
				if self._args[2]:
					text.write('\t\tfrom nodes;\n')
				else:
					text.write('\t\t'+self.links[0].string(self)+';\n')
			elif self.subtype == 'In line':
				rot0, globalV0, iNode0 = self.rigid_offset(0)
				localV0 = rot0*globalV0
				rot_1, globalV_1, Node_1 = self.rigid_offset(1)
				to_point = rot_1*(globalV_1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
				rot = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				text.write('\tjoint: '+str(self.database.Element.index(self))+', inline,\n')
				self.write_node(text, 0)
				rot_1, globalV_1, Node_1 = self.rigid_offset(1)
				to_point = rot_1*(globalV_1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
				text.write(',\n\t\t'+str(Node_1))
				text.write(',\n\t\t\toffset, '+str(to_point[0])+', '+str(to_point[1])+', '+str(to_point[2])+';\n')
			elif self.subtype == 'In plane':
				rot0, globalV0, iNode0 = self.rigid_offset(0)
				localV0 = rot0*globalV0
				rot1, globalV1, iNode1 = self.rigid_offset(1)
				to_point = rot1*(globalV1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
				rot = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				normal = rot*rot0*Vector(0., 0., 1.)
				text.write(
				'\tjoint: '+str(self.database.Element.index(self))+', inplane,\n'+
				'\t\t'+str(iNode0))
				text.write(',\n\t\t\t'+str(localV0[0])+', '+str(localV0[1])+', '+str(localV0[2]))
				text.write(',\n\t\t\t'+str(normal[0])+', '+str(normal[1])+', '+str(normal[2]))
				text.write(',\n\t\t'+str(iNode1))
				text.write(',\n\t\t\toffset, '+str(to_point[0])+', '+str(to_point[1])+', '+str(to_point[2])+';\n')
			elif self.subtype == 'Revolute hinge':
				self.write_hinge(text, 'revolute hinge')
				if self._args[2]:
					text.write(',\n\t\tinitial theta, '+str(self._args[3]))
				if self._args[4]:
					text.write(',\n\t\tfriction, '+str(self._args[5]))
					if self._args[6]:
						text.write(',\n\t\t\tpreload, '+str(self._args[7]))
					text.write(',\n\t\t\t'+self.links[0].string())
				text.write(';\n')
			elif self.subtype == 'Rod':
				text.write('\tjoint: '+str(self.database.Element.index(self))+', rod,\n')
				for i in range(2):
					self.write_node(text, i, orientation=False, p_label='position')
					text.write(',\n')
				text.write('\t\tfrom nodes,\n\t\t'+self.links[0].string()+';\n')
			elif self.subtype == 'Spherical hinge':
				self.write_hinge(text, 'spherical hinge')
				text.write(';\n')
			elif self.subtype == 'Total joint':
				rot_0, globalV_0, Node_0 = self.rigid_offset(0)
				localV_0 = rot_0*globalV_0
				rot_1, globalV_1, Node_1 = self.rigid_offset(1)
				to_joint = rot_1*(globalV_1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
				rot = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
				if Node_1 == self.objects[1]:
					rot_position = rot
				else:
					rot_position = self.objects[1].getMatrix().toQuat().toMatrix().transpose()
				text.write('\tjoint: '+str(self.database.Element.index(self))+', total joint')
				if not self._args[2]:
					text.write(',\n\t\t'+str(Node_0))
					text.write(', position, '+str(localV_0[0])+', '+str(localV_0[1])+', '+str(localV_0[2]))
					text.write(',\n\t\t\tposition orientation, matr,\n')
					self.rotationMatrix_write(rot_0*rot_position, text, '\t\t\t\t')
					text.write(',\n\t\t\trotation orientation, matr,\n')
					self.rotationMatrix_write(rot_0*rot, text, '\t\t\t\t')
				text.write(',\n\t\t'+str(Node_1))
				text.write(', position, '+str(to_joint[0])+', '+str(to_joint[1])+', '+str(to_joint[2]))
				text.write(',\n\t\t\tposition orientation, matr,\n')
				self.rotationMatrix_write(rot_1*rot_position, text, '\t\t\t\t')
				text.write(',\n\t\t\trotation orientation, matr,\n')
				self.rotationMatrix_write(rot_1*rot, text, '\t\t\t\t')
				if self._args[2]:
					text.write(',\n\t\t'+str(Node_0))
					text.write(', position, '+str(localV_0[0])+', '+str(localV_0[1])+', '+str(localV_0[2]))
					text.write(',\n\t\t\tposition orientation, matr,\n')
					self.rotationMatrix_write(rot_0*rot_position, text, '\t\t\t\t')
					text.write(',\n\t\t\trotation orientation, matr,\n')
					self.rotationMatrix_write(rot_0*rot, text, '\t\t\t\t')
				text.write(',\n\t\t\tposition constraint')
				for i in range(3):
					if self._args[3+2*i]:
						text.write(', active')
					else:
						text.write(', inactive')
				text.write(', component')
				for i in range(3):
					if self._args[3+2*i]:
						text.write(',\n\t\t\t\t'+self.links[i].string())
					else:
						text.write(',\n\t\t\t\tinactive')
				text.write(',\n\t\t\torientation constraint')
				for i in range(3):
					if self._args[9+2*i]:
						text.write(', active')
					else:
						text.write(', inactive')
				text.write(', component')
				for i in range(3):
					if self._args[9+2*i]:
						text.write(',\n\t\t\t\t'+self.links[3+i].string())
					else:
						text.write(',\n\t\t\t\tinactive')
				text.write(';\n')
		elif self.type == 'Gravity':
			text.write('\tgravity: '+self.links[0].string()+', '+self.links[1].string(self)+';\n')
		elif self.type == 'Air properties':
			if self._args[0]:
				text.write('\tair properties: std, SI,\n')
			elif self._args[0]:
				text.write('\tair properties: std, British,\n')
			text.write('\t\ttemperature deviation, '+str(self._args[2])+
				',\n\t\treference altitude, '+str(self._args[3])+','+
				self.links[0].string()+', '+self.links[1].string(self))
			if self._args[5]:
				text.write(',\n\t\tgust, front 1D,\n\t\t'+
					self.links[2].string()+',\n\t\t\t'+
					self.links[3].string()+',\n\t\t\t'+
					str(self._args[9])+', '+self.links[4].string(self))
			text.write(';\n')
		elif self.type == 'Driven':
			text.write('\tdriven: '+str(self.database.Element.index(self.links[1]))+', '+
			self.links[0].string(self)+',\n'+
			'\t\texisting: '+self.links[1].type+', '+str(self.database.Element.index(self.links[1]))+';\n')

	def write_hinge(self, text, name, V1=True, V2=True, M1=True, M2=True):
		rot_0, globalV_0, Node_0 = self.rigid_offset(0)
		localV_0 = rot_0*globalV_0
		rot_1, globalV_1, Node_1 = self.rigid_offset(1)
		to_hinge = rot_1*(globalV_1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
		rotT = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
		text.write(
		'\tjoint: '+str(self.database.Element.index(self))+', '+name+',\n'+
		'\t\t'+str(Node_0))
		if V1:
			text.write(', '+str(localV_0[0])+', '+str(localV_0[1])+', '+str(localV_0[2]))
		if M1:
			text.write(',\n\t\t\thinge, matr,\n')
			self.rotationMatrix_write(rot_0*rotT, text, '\t\t\t\t')
		text.write(', \n\t\t'+str(Node_1))
		if V2:
			text.write(', '+str(to_hinge[0])+', '+str(to_hinge[1])+', '+str(to_hinge[2]))
		if M2:
			text.write(',\n\t\t\thinge, matr,\n')
			self.rotationMatrix_write(rot_1*rotT, text, '\t\t\t\t')

	def write_node(self, text, i, node=True, position=True, orientation=True, p_label='', o_label=''):
		rot_i, globalV_i, Node_i = self.rigid_offset(i)
		localV_i = rot_i*globalV_i
		rotT = self.objects[i].getMatrix().toQuat().toMatrix().transpose()
		if node:
			text.write('\t\t'+str(Node_i)+',\n')
		if position:
			text.write('\t\t\t')
			if p_label:
				text.write(p_label+', ')
			text.write(str(localV_i[0])+', '+str(localV_i[1])+', '+str(localV_i[2]))
		if orientation:
			text.write(',\n\t\t\t')
			if o_label:
				text.write(o_label+', ')
			text.write('matr,\n')
			self.rotationMatrix_write(rot_i*rotT, text, '\t\t\t\t')

	def rigid_offset(self, i):
		if self.objects[i] in self.database.Node:
			ob = self.objects[i]
		elif self.database.rigid_dict.has_key(self.objects[i]):
			ob = self.database.rigid_dict[self.objects[i]]
		else:
			print 'Model Error: Object '+self.objects[i].name+' is not associated with a Node'
		rot = ob.getMatrix().toQuat().toMatrix()
		globalV = Vector(self.objects[i].loc) - Vector(ob.loc)
		return rot, globalV, self.database.Node.index(ob)

	def draw(self):
		if not self.objects:
			return
		editmode = Window.EditMode()    # are we in edit mode?  If so ...
		if editmode: Window.EditMode(0) # leave edit mode before getting the mesh

		obj = self.objects[0]
		saveData = obj.getData(mesh=True)
		exists = False
		for mesh in bpy.data.meshes:
			if mesh.name == '_'+obj.name:
				exists = True
				me = mesh
		if not exists:
			me = bpy.data.meshes.new('_'+obj.name)

		if self.type == 'Aerodynamic':
			if self.subtype == 'Aerodynamic body':
				b = self._args[3]
				ords = [[-.5*b, 0., .5*b]]+[[]]*3
				t = .12
				for i, j in enumerate([0, 1, 3]):
					if self.links[j+1].type == 'const':
						ords[i+1] = [self.links[j+1]._args[0]]*3
					if self.links[j+1].type == 'linear':
						ords[i+1] = [self.links[j+1]._args[0], 0., self.links[j+1]._args[1]]
						ords[i+1][1] = (ords[i+1][0] + ords[i+1][2])/2.
					if self.links[j+1].type == 'piecewise linear':
						N = self.links[j+1]._args[0]
						ords[i+1] = [self.links[j+1]._series[1], 0., self.links[j+1]._series[2*N-1]]
						for i in range(N):
							if self.links[j+1]._series[2*i] > 0:
								break
						ords[i+1][1] = (self.links[j+1]._series[2*i-1] + self.links[j+1]._series[2*i+1])/2.
					if self.links[j+1].type == 'parabolic':
						ords[i+1] = self.links[j+1]._args[:3]
			else:
				ords = [[]]*4
				t = .12
				for i, j in enumerate([1, 2, 4]):
					if self.links[j+1].type == 'const':
						ords[i+1] = [self.links[j+1]._args[0]]*3
					if self.links[j+1].type == 'linear':
						ords[i+1] = [self.links[j+1]._args[0], 0., self.links[j+1]._args[1]]
						ords[i+1][1] = (ords[i+1][0] + ords[i+1][2])/2.
					if self.links[j+1].type == 'piecewise linear':
						N = self.links[j+1]._args[0]
						ords[i+1] = [self.links[j+1]._series[1], 0., self.links[j+1]._series[2*N-1]]
						for i in range(N):
							if self.links[j+1]._series[2*i] > 0:
								break
						ords[i+1][1] = (self.links[j+1]._series[2*i-1] + self.links[j+1]._series[2*i+1])/2.
					if self.links[j+1].type == 'parabolic':
						ords[i+1] = self.links[j+1]._args[:3]
				b = (ords[1][0] + ords[1][1] + ords[1][2])/3.
				ords[0] = [-.5*b, 0., .5*b]
			coords = []
			xLE, yLE = LeadingEdge = (0.292, 0.850)
			xTE, yTE = TrailingEdge = (-0.792, 0.086)
			for z, c, ac, a in zip(ords[0], ords[1], ords[2], ords[3]):
				acV = Vector(ac, 0., 0.)
				rot = Mathutils.RotationMatrix(-a, 3, 'z')
				for ky in [-1., 1.]:
					coords.append(rot*Vector(c*xTE, t*c*yTE*ky, z)+acV)
					coords.append(rot*Vector(c*xLE, t*c*yLE*ky, z)+acV)
			faces = [[5,4,0,1],[7,5,1,3],[6,7,3,2],[4,6,2,0]]
			faces += [[4 + face[i] for i in range(4)] for face in faces] 
			faces += [[2,3,1,0], [11,10,8,9]]
			me.verts = None
			me.verts.extend(coords)          # add vertices to mesh
			me.faces.extend(faces)           # add faces to the mesh (also adds edges)
			for i in me.findEdges([(0,1),(0,2),(2,3),(3,1),	(8,9),(8,10),(10,11),(11,9)]):
				me.edges[i].crease = 255
			for i in me.findEdges([(4,5),(4,6),(6,7),(7,5)]):
				me.edges[i].crease = 127
			for ob in self.objects:
				ob.link(me)
#			for i in range(8):
#				me.faces[i].smooth = 1
		elif self.type == 'Body':
			mass = self._args[1]
			shell = []
			coords = []
			if self.links[0]._args[4]:
				II = [1.]*3
			else:
				II = [self.links[0]._series[4*i] for i in range(3)]
			for I in II:
				shell.append(0.5*sqrt(I/mass))
			for z in [-1., 1.]:
				for y in [-1., 1.]:
					for x in [-1., 1.]:
						coords.append(Vector(x*shell[0],y*shell[1],z*shell[2]))
			faces = [[1,0,2,3],[4,5,7,6],[0,1,5,4],[1,3,7,5],[3,2,6,7],[2,0,4,6]]
			me.verts = None
			me.verts.extend(coords)          # add vertices to mesh
			me.faces.extend(faces)           # add faces to the mesh (also adds edges)
			for e in me.edges: e.crease = 47
#			for f in me.faces: f.smooth = 1
		elif self.type == 'Beam':
			coords = []
			shell = [0.3, 0.1, 0.2]
			for z in [-1., 1.]:
				for y in [-1., 1.]:
					for x in [-1., 1.]:
						coords.append(Vector(x*shell[0],y*shell[1],z*shell[2]))
			faces = [[1,0,2,3],[4,5,7,6],[0,1,5,4],[1,3,7,5],[3,2,6,7],[2,0,4,6]]
			me.verts = None
			me.verts.extend(coords)          # add vertices to mesh
			me.faces.extend(faces)           # add faces to the mesh (also adds edges)
			for e in me.edges: e.crease = 255
			for ob in self.objects:
				ob.link(me)
#			for f in me.faces: f.smooth = 1
		elif self.type == 'Joint':
			if self.subtype in ['Revolute hinge', 'Axial rotation']:
				coords = []
				scale = .5
				for z in [-1., 1.]:
					for y in [-1., 1.]:
						for x in [-1., 1.]:
							coords.append(Vector(scale*x,scale*y,scale*z))
				faces = [[1,0,2,3],[4,5,7,6],[0,1,5,4],[1,3,7,5],[3,2,6,7],[2,0,4,6]]
				me.verts = None
				me.verts.extend(coords)          # add vertices to mesh
				me.faces.extend(faces)           # add faces to the mesh (also adds edges)
				for i in me.findEdges([(0,1),(0,2),(2,3),(3,1),	(4,5),(4,6),(6,7),(7,5)]):
					me.edges[i].crease = 255
#				for i in range(2,6): me.faces[i].smooth = 1
			elif self.subtype == 'Spherical hinge':
				coords = []
				scale = .5
				for z in [-1., 1.]:
					for y in [-1., 1.]:
						for x in [-1., 1.]:
							coords.append(Vector(scale*x,scale*y,scale*z))
				faces = [[1,0,2,3],[4,5,7,6],[0,1,5,4],[1,3,7,5],[3,2,6,7],[2,0,4,6]]
				me.verts = None
				me.verts.extend(coords)          # add vertices to mesh
				me.faces.extend(faces)           # add faces to the mesh (also adds edges)
#				for f in me.faces: f.smooth = 1
			elif self.subtype == 'Clamp':
				coords = []
				scale = .5
				for y in [-1., 1.]:
					for x in [-1., 1.]:
						coords.append(Vector(scale*x,scale*y,scale*(-1.)))
				coords.append(Vector(0.,0.,0.))
				faces = [[2,3,1,0],[0,1,4],[1,3,4],[3,2,4],[2,0,4]]
				me.verts = None
				me.verts.extend(coords)          # add vertices to mesh
				me.faces.extend(faces)           # add faces to the mesh (also adds edges)
#				for f in me.faces: f.smooth = 1
				for i in range(4,8): me.edges[i].crease = 255
			else:
				me = saveData
		elif self.type == 'Rigid':
			me.verts = None
			me.verts.extend([(.333,0.,0.),(0.,.666,0.),(-.333,0.,0.),(0.,-.666,0.),(0.,0.,1.)]) 
			me.faces.extend([(3,2,1,0),(0,1,4),(1,2,4),(2,3,4),(3,0,4)])
			for f in me.faces: f.smooth = 0
			for e in me.edges: e.crease = 255
		else:
			me = saveData
			
		obj.link(me)	
		mods = obj.modifiers
		hasSubsurf = False
		for mod in mods:
			if mod.name == 'Subsurf': hasSubsurf = True
		if not hasSubsurf:
			mod = mods.append(Modifier.Type.SUBSURF)
			mod[Modifier.Settings.LEVELS] = 3     # set subsurf subdivision levels to 3
		obj.makeDisplayList()
		if editmode: Window.EditMode(1)  # optional, just being nice
		Window.Redraw()

class Constitutive(Common):

	types = [
	'Const_1D',
	'Const_3D',
	'Const_6D']

	Const_1D = [
	'Linear elastic(1D)',
	'Linear elastic generic(1D)',
	'Cubic elastic generic(1D)',
	'Log elastic(1D)',
#	'Linear elastic bi-stop generic(1D)',
	'Double linear elastic(1D)',
	'Isotropic hardening elastic(1D)',
	'Scalar function elastic(1D)',
	'Scalar function elastic orthotropic(1D)',
	'nlsf elastic(1D)',
	'nlp elastic(1D)',
	'Linear viscous(1D)',
	'Linear viscous generic(1D)',
	'nlsf viscous(1D)',
	'nlp viscous(1D)',
	'Linear viscoelastic(1D)',
	'Linear viscoelastic generic(1D)',
	'Cubic viscoelastic generic(1D)',
	'Double linear viscoelastic(1D)',
	'Turbulent viscoelastic(1D)',
#	'Linear viscoelastic bi-stop generic(1D)',
	'GRAALL damper(1D)',
	'shock absorber(1D)',
	'symbolic elastic(1D)',
	'symbolic viscous(1D)',
	'symbolic viscoelastic(1D)',
	'ann elastic(1D)',
	'ann viscoelastic(1D)',
	'nlsf viscoelastic(1D)',
	'nlp viscoelastic(1D)']

	Const_3D = [
	'Linear elastic(3D)',
	'Linear elastic generic(3D)',
	'Cubic elastic generic(3D)',
#	'Linear elastic bi-stop generic(3D)',
	'Double linear elastic(3D)',
	'Isotropic hardening elastic(3D)',
	'Scalar function elastic(3D)',
	'Scalar function elastic orthotropic(3D)',
	'nlsf elastic(3D)',
	'nlp elastic(3D)',
	'Linear viscous(3D)',
	'Linear viscous generic(3D)',
	'nlsf viscous(3D)',
	'nlp viscous(3D)',
	'Linear viscoelastic(3D)',
	'Linear viscoelastic generic(3D)',
	'Cubic viscoelastic generic(3D)',
	'Double linear viscoelastic(3D)',
#	'Linear viscoelastic bi-stop generic(3D)',
	'GRAALL damper(3D)',
	'symbolic elastic(3D)',
	'symbolic viscous(3D)',
	'symbolic viscoelastic(3D)',
	'ann elastic(3D)',
	'ann viscoelastic(3D)',
	'nlsf viscoelastic(3D)',
	'nlp viscoelastic(3D)',
	'invariant angular(3D)']

	Const_6D = [
	'Linear elastic(6D)',
	'Linear elastic generic(6D)',
	'Linear elastic generic axial torsion coupling(6D)',
#	'Linear elastic bi-stop generic(6D)',
	'Isotropic hardening elastic(6D)',
	'Scalar function elastic(6D)',
	'Scalar function elastic orthotropic(6D)',
	'nlsf elastic(6D)',
	'nlp elastic(6D)',
	'Linear viscous(6D)',
	'Linear viscous generic(6D)',
	'nlsf viscous(6D)',
	'nlp viscous(6D)',
	'Linear viscoelastic(6D)',
	'Linear viscoelastic generic(6D)',
	'Linear viscoelastic generic axial torsion coupling(6D)',
#	'Linear viscoelastic bi-stop generic(6D)',
	'GRAALL damper(6D)',
	'symbolic elastic(6D)',
	'symbolic viscous(6D)',
	'symbolic viscoelastic(6D)',
	'ann elastic(6D)',
	'ann viscoelastic(6D)',
	'nlsf viscoelastic(6D)',
	'nlp viscoelastic(6D)']

	def __init__(self, constType, database, linking=False):
		if constType in Constitutive.Const_1D:
			self.type = 'Const_1D'
			self.eligible = ['1x1']
			self.dimension = 1
		elif constType in Constitutive.Const_3D:
			self.type = 'Const_3D'
			self.eligible = ['3x3']
			self.dimension = 3
		elif constType in Constitutive.Const_6D:
			self.type = 'Const_6D'
			self.eligible = ['6x6']
			self.dimension = 6
		self.subtype = constType[:-4]
		self.name = self.subtype
		self.database = database
		self.users = 0
		self.links = []
		if self.subtype == 'Linear elastic':
			self._args = [1.]
		elif self.subtype == 'Linear elastic generic':
			self._args = [1]
			self.links = [None]
		elif self.subtype == 'Linear elastic generic axial torsion coupling':
			self.eligible = ['6x1']
			self._args = [1, 0.]
			self.links = [None]
		elif self.subtype == 'Cubic elastic generic':
			if self.dimension == 3:
				self.eligible = ['3x1']
			self._args = [1, 1, 1]
			self.links = [None]*3
		elif self.subtype == 'Log elastic':
			self._args = [0.]
		elif self.subtype == 'Linear elastic bi-stop generic':
			self._args = [1, 0, 0, 1, 1]
			self.links = [None]*3
		elif self.subtype == 'Double linear elastic':
			self._args = [0.]*4
		elif self.subtype == 'Isotropic hardening elastic':
			self._args = [0., 0., 0, 0.]
		elif self.subtype == 'Scalar function elastic':
			self._args = [1]
			self.links = [None]
		elif self.subtype == 'Scalar function elastic orthotropic':
			self._args = [1]*self.dimension
			self.links = [None]*self.dimension
		elif self.subtype == 'nlsf elastic' or self.subtype == 'nlp elastic':
			dim = self.dimension
			self._args = [1]+[1]*dim
			self.links = [None]*(1+dim)
		elif self.subtype == 'Linear viscous':
			self._args = [0.]
		elif self.subtype == 'Linear viscous generic':
			self._args = [1]
			self.links = [None]
		elif self.subtype == 'nlsf viscous' or self.subtype == 'nlp viscous':
			dim = self.dimension
			self._args = [1]+[1]*dim
			self.links = [None]*(1+dim)
		elif self.subtype == 'Linear viscoelastic':
			self._args = [0., 0., 0, 0.]
		elif self.subtype == 'Linear viscoelastic generic':
			self._args = [1, 1, 0, 0.]
			self.links = [None]*2
		elif self.subtype == 'Linear viscoelastic generic axial torsion coupling':
			self.eligible = ['6x1']
			self._args = [1, 1, 0, 0., 0.]
			self.links = [None]*2
		elif self.subtype == 'Cubic viscoelastic generic':
			self._args = [1]*4
			self.links = [None]*4
		elif self.subtype == 'Double linear viscoelastic':
			self._args = [0.]*5 + [0] + [0.]
		elif self.subtype == 'Turbulent viscoelastic':
			self._args = [0., 0., 0, 0., 0, 0.]
		elif self.subtype == 'Linear viscoelastic bi-stop generic':
			self._args = [1, 1, 0, 0, 1, 1]
			self.links = [None]*4
			"""
		elif self.subtype == 'GRAALL damper':
		elif self.subtype == 'shock absorber':
		elif self.subtype == 'symbolic elastic':
		elif self.subtype == 'symbolic viscous':
		elif self.subtype == 'symbolic viscoelastic':
		elif self.subtype == 'ann elastic':
		elif self.subtype == 'ann viscoelastic':
			"""
		elif self.subtype == 'nlsf viscoelastic' or self.subtype == 'nlp viscoelastic':
			dim = self.dimension
			self._args = [1]+[1]*dim+[1]+[0]*2+[1]*dim
			self.links = [None]*(2+2*dim)

#		elif self.subtype == 'invariant angular ':

		self.name_check(self.database.Constitutive)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.subtype == 'Linear elastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',	args[0], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Linear elastic generic':
			tmp = ['select matrix']
			for i in range(1):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change stiffness matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', self.eligible)
				self._args = [0]
		elif self.subtype == 'Linear elastic generic axial torsion coupling':
			tmp = ['select matrix']
			for i in range(1):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change stiffness matrix'),
			('Coef:',	args[1], 0., 9.9e10, 'Coupling coefficient')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', self.eligible)
				self._args[0] = 0
		elif self.subtype == 'Cubic elastic generic':
			tmp = ['select']*3
			for i in range(3):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Stiffness matrix 1'),
			(tmp[1],	args[1], 'Stiffness matrix 2'),
			(tmp[2], 	args[2], 'Stiffness matrix 3')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				for i in range(3):
					if self._args[i] or not self.links[i]:
						self.persistantPupMenu('Select S'+str(i+1))
						self.select(i, 'Matrix', self.eligible)
				self._args = [0, 0, 0]
		elif self.subtype == 'Log elastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',	args[0], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Linear elastic bi-stop generic':
			tmp = ['select matrix'] + ['select drive']*2
			for i in range(3):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],			args[0], 'Change stiffness matrix'),
			('Initial state',	args[1], 'Prescribe initial state'),
			('Active',			args[2], 'Select initial state as inactive or active'),
			(tmp[1],			args[3], 'Change activitating condition'),
			(tmp[2], 			args[4], 'Change deactivating condition')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', self.eligible)
				if args[3].val or not self.links[1]:
					self.persistantPupMenu('Select activating condition:')
					self.select(1, 'Drive')
				if args[4].val or not self.links[2]:
					self.persistantPupMenu('Select deactivating condition:')
					self.select(2, 'Drive')
				self._args[0] = 0
				self._args[3:5] = [0]*2
		elif self.subtype == 'Double linear elastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness 1:',	args[0], 0., 9.9e10),
			('Upper strain:',	args[1], 0., 9.9e10),
			('Lower Strain:',	args[2], 0., 9.9e10),
			('Stiffness 2:',	args[3], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Isotropic hardening elastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',			args[0], 0., 9.9e10),
			('Reference strain:',	args[1], 0., 9.9e10),
			('Linear stiffness',	args[2], 'Use linear stiffness parameter'),
			('Linear stiffness:',	args[3], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Scalar function elastic':
			tmp = 'select'
			if self.links[0]:
				tmp = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp, args[0], 'Scalar function')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select scalar function')
					self.select(0, 'Function')
				self._args[0] = 0
		elif self.subtype == 'Scalar function elastic orthotropic':
			string = [('Name: ', nval, 0, 30)]
			for i in range(self.dimension):
				if not self.links[i]:
					string += [('select  f-'+str(i+1)+':', args[i], 'Select function')]
				else:
					string += [(self.links[i].name, args[i],'Modify f-'+str(i+1))]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				for i, function in enumerate(self.links):
					if self._args[i] or not function:
						self.persistantPupMenu('Select  f-'+str(i+1)+':')
						self.select(i, 'Function')
			self._args = [0]*self.dimension
		elif self.subtype == 'nlsf elastic' or self.subtype == 'nlp elastic':
			dim = self.dimension
			tmp = []
			for i, link in enumerate(self.links):
				if self.links[i]:
					tmp.append(self.links[i].name)
				else:
					tmp.append('select')
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],		args[0], 'Stiffness matrix')]
			for i in range(1, 1+dim):
				string += [(tmp[i], args[i], 'Stiffness function-'+str(i))]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', self.eligible)
				for i in range(1, 1+dim):
					if self._args[i] or not self.links[i]:
						self.persistantPupMenu('Select stiffness function-'+str(i))
						self.select(i, 'Function')
			self._args[0:2+dim] = [0]*(2+dim)
		elif self.subtype == 'Linear viscous':
			string = [
			('Name: ', nval, 0, 30),
			('Viscosity:',	args[0], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Linear viscous generic':
			tmp = ['select matrix']
			for i in range(1):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change viscosity matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select viscosity matrix:')
					self.select(0, 'Matrix', self.eligible)
				self._args = [0]
		elif self.subtype == 'nlsf viscous' or self.subtype == 'nlp viscous':
			dim = self.dimension
			tmp = []
			for i, link in enumerate(self.links):
				if self.links[i]:
					tmp.append(self.links[i].name)
				else:
					tmp.append('select')
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],		args[0], 'Viscosity matrix')]
			for i in range(1, 1+dim):
				string += [(tmp[i], args[i], 'Viscosity function-'+str(i))]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(0, 'Matrix', self.eligible)
				for i in range(1, 1+dim):
					if self._args[i] or not self.links[i]:
						self.persistantPupMenu('Select viscosity function-'+str(i))
						self.select(i, 'Function')
			self._args[0:2+dim] = [0]*(2+dim)
		elif self.subtype == 'Linear viscoelastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',	args[0], 0., 9.9e10),
			('Viscosity:',	args[1], 0., 9.9e10),
			('Proportional',args[2], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Linear viscoelastic generic':
			tmp = ['select']*2
			for i in range(2):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],		args[0], 'Stiffness matrix'),
			(tmp[1], 		args[1], 'Viscosity matrix'),
			('Proportional',args[2], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', self.eligible)
				if not self._args[2] and (self._args[1] or not self.links[1]):
					self.persistantPupMenu('Select viscosity matrix')
					self.select(1, 'Matrix', self.eligible)
			self._args[0:2] = [0, 0]
		elif self.subtype == 'Linear viscoelastic generic axial torsion coupling':
			tmp = ['select']*2
			for i in range(2):
				if self.links[i]:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],		args[0], 'Stiffness matrix'),
			(tmp[1], 		args[1], 'Viscosity matrix'),
			('Proportional',args[2], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3], 0., 9.9e10),
			('Coupling:',	args[4], 0., 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', self.eligible)
				if self._args[1] or not self.links[1]:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(1, 'Matrix', self.eligible)
			self._args[0:2] = [0, 0]
		elif self.subtype == 'Cubic viscoelastic generic':
			tmp = ['select']*4
			for i, link in enumerate(self.links):
				if link:
					tmp[i] = link.name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Stiffness matrix 1'),
			(tmp[1],	args[1], 'Stiffness matrix 2'),
			(tmp[2], 	args[2], 'Stiffness matrix 3'),
			(tmp[3], 	args[3], 'Viscosity matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			eligible_S_matrix = self.eligible
			if self.type == 'Const_3D':
				eligible_S_matrix = ['3x1']
			if not delete.val:
				self._args = [arg.val for arg in args]
				for i in range(3):
					if self._args[i] or not self.links[i]:
						self.persistantPupMenu('Select S'+str(i+1))
						self.select(i, 'Matrix', eligible_S_matrix)
				if self._args[3] or not self.links[3]:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(3, 'Matrix', self.eligible)
				self._args = [0]*4
		elif self.subtype == 'Double linear viscoelastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness-1:',	args[0], -9.9e10, 9.9e10),
			('Upper strain:',	args[1], -9.9e10, 9.9e10),
			('Lower Strain:',	args[2], -9.9e10, 9.9e10),
			('Stiffness-2:',	args[3], -9.9e10, 9.9e10),
			('Viscosity:',		args[4], -9.9e10, 9.9e10),
			('Second damping',	args[5], 'Use viscosity-2 for outside lower/upper strain range'),
			('Viscosity-2:',	args[6], -9.9e10, 9.9e10)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Turbulent viscoelastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',		args[0], 0., 9.9e10),
			('Parabolic visc:',	args[1], 0., 9.9e10),
			('Use threshold',	args[2]),
			('Threshold:',		args[3], 0., 9.9e10),
			('Use linear visc',	args[4]),
			('Linear visc:',	args[5], 0., 9.9e10, 'Linear viscous coefficient')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			self._args = [arg.val for arg in args]
		elif self.subtype == 'Linear viscoelastic bi-stop generic':
			tmp = ['select matrix']*2 + ['select drive']*2
			for i, link in enumerate(self.links):
				if link:
					tmp[i] = link.name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],			args[0], 'Change stiffness matrix'),
			(tmp[1],			args[1], 'Change viscosity matrix'),
			('Enable state',	args[2], 'Prescribe initial state'),
			('Active',			args[3], 'Initial state is inactive or active'),
			(tmp[2],			args[4], 'Change activitating condition'),
			(tmp[3], 			args[5], 'Change deactivating condition')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', self.eligible)
				if self._args[1] or not self.links[1]:
					self.persistantPupMenu('Select viscosity matrix:')
					self.select(1, 'Matrix', self.eligible)
				if self._args[4] or not self.links[2]:
					self.persistantPupMenu('Select activating condition:')
					self.select(2, 'Drive')
				if self._args[5] or not self.links[3]:
					self.persistantPupMenu('Select deactivating condition:')
					self.select(3, 'Drive')
				self._args[0:2] = [0]*2
				self._args[4:6] = [0]*2
				"""
		elif self.subtype == 'GRAALL damper':
		elif self.subtype == 'shock absorber':
		elif self.subtype == 'symbolic elastic':
		elif self.subtype == 'symbolic viscous':
		elif self.subtype == 'symbolic viscoelastic':
		elif self.subtype == 'ann elastic':
		elif self.subtype == 'ann viscoelastic':
				"""
		elif self.subtype == 'nlsf viscoelastic' or self.subtype == 'nlp viscoelastic':
			dim = self.dimension
			tmp = []
			for i, link in enumerate(self.links):
				if self.links[i]:
					tmp.append(self.links[i].name)
				else:
					tmp.append('select')
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],		args[0], 'Stiffness matrix')]
			for i in range(1, 1+dim):
				string += [(tmp[i], args[i], 'Stiffness function-'+str(i))]
			string += [
			(tmp[1+dim], 	args[1+dim], 'Viscosity matrix'),
			('Proportional',args[2+dim], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3+dim], 0., 9.9e10)]
			for i in range(4+dim, 4+2*dim):
				string += [(tmp[i-2], args[i], 'Viscosity function-'+str(i-3-dim))]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.subtype, string): return
			if not delete.val:
				self._args = [arg.val for arg in args]
				if self._args[0] or not self.links[0]:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', self.eligible)
				for i in range(1, 1+dim):
					if self._args[i] or not self.links[i]:
						self.persistantPupMenu('Select stiffness function-'+str(i))
						self.select(i, 'Function')
				if not self._args[2+dim] and (self._args[1+dim] or not self.links[1+dim]):
					self.persistantPupMenu('Select viscosity matrix')
					self.select(1+dim, 'Matrix', self.eligible)
				for i in range(4+dim, 4+2*dim):
						if self._args[i] or not self.links[i-2]:
							self.persistantPupMenu('Select viscosity function-'+str(i-3-dim))
							self.select(i-2, 'Function')
			self._args[0:2+dim] = [0]*(2+dim)
			self._args[4+dim:4+2*dim] = [0]*dim

#		elif self.subtype == 'invariant angular ':
		else:
			Draw.PupMenu('Error: '+self.subtype+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Constitutive)

	def string(self):
		if self.subtype in ['Linear elastic', 'Linear viscous']:
			if self.type == 'Const_1D':
				string = self.subtype+', '
			else:
				string = self.subtype+' isotropic, '
			string += str(self._args[0])
		elif self.subtype in ['Linear elastic generic', 'Linear viscous generic']:
			string = self.subtype+', '+self.links[0].string()
		elif self.subtype == 'Linear elastic generic axial torsion coupling':
			string = ('linear elastic generic axial torsion coupling, diag,'+
				self.links[0].string()+', '+str(self._args[1]))
		elif self.subtype == 'Cubic elastic generic':
			string = 'cubic elastic generic'
			for link in self.links:
				string += ','+link.string()
		elif self.subtype == 'Log elastic':
			string = 'log elastic, '+str(self._args[0])
		elif self.subtype == 'Linear elastic bi-stop generic':
			string = 'linear elastic bistop, '+self.links[0].string()+',\n\t\t'
			if self._args[1]:
				if self._args[2]:
					string += 'initial state, active, \n\t\t'
				else:
					string += 'initial state, inactive, \n\t\t'
			string += (self.links[1].string(self)+',\n\t\t'+self.links[2].string(self))
		elif self.subtype == 'Double linear elastic':
			string = 'double linear elastic,\n\t\t\t'+str(self._args[0:4]).strip('[]')
		elif self.subtype == 'Isotropic hardening elastic':
			string = 'isotropic hardening elastic, '+str(self._args[0:2]).strip('[]')
			if self._args[2]:
				string += ',\n\t\t\tlinear stiffness, '+str(self._args[3])
		elif self.subtype == 'Scalar function elastic':
			string = 'scalar function elastic isotropic, "'+self.links[0].name+'"'
		elif self.subtype == 'Scalar function elastic orthotropic':
			string = 'scalar function elastic orthotropic,\n\t\t\t'
			for link in self.links:
				string += '"'+link.name+'", '
			string = string[:-2]
		elif self.subtype == 'Linear viscoelastic':
			if self.type == 'Const_1D':
				string = ('linear viscoelastic, '+str(self._args[0])+', ')
			else:
				string = ('linear viscoelastic isotropic, '+str(self._args[0])+', ')
			if self._args[2]:
				string += 'proportional, '+str(self._args[3])
			else:
				string += str(self._args[1])
		elif self.subtype == 'Linear viscoelastic generic':
			string = 'linear viscoelastic generic, '+self.links[0].string()+', '
			if self._args[2]:
				string += '\n\t\tproportional, '+str(self._args[3])
			else:
				string += self.links[1].string()
		elif self.subtype == 'Linear viscoelastic generic axial torsion coupling':
			string = 'linear viscoelastic generic axial torsion coupling, \n\t\tdiag, '
			string += self.links[0].string()
			if self._args[2]:
				string += ',\n\t\tproportional, '+str(self._args[3])
			else:
				string += ',\n\t\tdiag, ' + self.links[1].string()
			string += ',\n\t\t'+str(self._args[4])
		elif self.subtype == 'Cubic viscoelastic generic':
			string = 'cubic viscoelastic generic'
			for i in range(3):
				string += ','+self.links[i].string()
			string += ','+self.links[3].string()
		elif self.subtype == 'Double linear viscoelastic':
			string = 'double linear viscoelastic,\n\t\t\t'
			for arg in self._args[:5]:
				string += str(arg)+', '
			string = string[:-2]
			if self._args[5]:
				string += ', \n\t\t\tsecond damping, '+str(self._args[6])
		elif self.subtype == 'Turbulent viscoelastic':
			string = 'turbulent viscoelastic,\n\t\t\t'
			for arg in self._args[:2]:
				string += str(arg)+', '
			string = string[:-2]
			if self._args[2]:
				string += ', '+str(self._args[3])
			if self._args[4]:
				string += ', '+str(self._args[5])
		elif self.subtype == 'Linear viscoelastic bi-stop generic':
			string = ('linear viscoelastic bistop, '+
				self.links[0].string()+', '+self.links[1].string())
			if self._args[2]:
				if self._args[3]:
					string += ',\n\t\tinitial state, active'
				else:
					string += ',\n\t\tinitial state, inactive'
			string += (',\n\t\t'+self.links[2].string(self)+
				',\n\t\t'+self.links[3].string(self))
#		elif self.subtype == 'GRAALL damper':
#		elif self.subtype == 'shock absorber':
#		elif self.subtype == 'symbolic elastic':
#		elif self.subtype == 'symbolic viscous':
#		elif self.subtype == 'symbolic viscoelastic':
#		elif self.subtype == 'ann elastic':
#		elif self.subtype == 'ann viscoelastic':
		elif (self.subtype == 'nlsf elastic' or self.subtype == 'nlp elastic' or
			self.subtype == 'nlsf viscous' or self.subtype == 'nlp viscous'):
			dim = self.dimension
			string = self.subtype+', '+self.links[0].string()
			for i in range(1, 1+dim):
				string += ',\n\t\t\t\t"'+self.links[i].name+'"'
		elif self.subtype == 'nlsf viscoelastic' or self.subtype == 'nlp viscoelastic':
			dim = self.dimension
			string = self.subtype+', '+self.links[0].string()+','
			for i in range(1, 1+dim):
				string += '\n\t\t\t\t"'+self.links[i].name+'",'
			if self._args[2+dim]:
				string += '\n\t\t\tproportional, '+str(self._args[3+dim])
			else:
				string += self.links[1+dim].string()
			for i in range(4+dim, 4+2*dim):
				string += ',\n\t\t\t\t"'+self.links[i-2].name+'"'

#		elif self.subtype == 'invariant angular ':
		return string

class Shape(Common):

	types = [
	'const',
	'linear',
	'piecewise linear',
	'parabolic']

	def __init__(self, shapeType, database, linking=False):
		self.type = shapeType
		self.name = shapeType
		self.database = database
		self.users = 0
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'piecewise linear':
			self._args = [0]
		self.name_check(self.database.Shape)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'const':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10)]
			if self in self.database.Shape:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'linear':
			string = [
			('Name: ', nval, 0, 30),
			('Value at -1:',	args[0], -9.9e10, 9.9e10),
			('Value at +1:',	args[1], -9.9e10, 9.9e10)]
			if self in self.database.Shape:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == 'piecewise linear':
			string = [
			('Name: ', nval, 0, 30),
			('Number of pts:', args[0], 2, 50)]
			if self in self.database.Shape:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				N = args[0].val
				if 2*N > len(self._series):
					self._series[len(self._series):2*N] = [0. for i in range(len(self._series),2*N)]
				series = [Draw.Create(self._series[i]) for i in range(2*N)]
				seq = ['']*(16*((N-1)/8+1))
				for col in range((N-1)/8+1):
					for i in range(min(8, N-col*8)):
						seq[16*col+i] = ('Point_'+str(8*col+i)+': ', series[16*col+2*i], -9.9e10, 9.9e10)
						seq[16*col+8+i] = ('Value_'+str(8*col+i)+': ', series[16*col+2*i+1], -9.9e10, 9.9e10)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [N]
		elif self.type == 'parabolic':
			string = [
			('Name: ', nval, 0, 30),
			('Value at -1:',	args[0], -9.9e10, 9.9e10),
			('Value at 0:',		args[1], -9.9e10, 9.9e10),
			('Value at +1:',	args[2], -9.9e10, 9.9e10)]
			if self in self.database.Shape:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Shape)

	def string(self, scale=1.):
		if self.type == 'const':
			return '\t\tconst, '+str(self._args[0]*scale)
		elif self.type == 'linear':
			string = '\t\tlinear'
			for arg in self._args[:2]:
				string += ', '+str(arg*scale)
			return string
		elif self.type == 'piecewise linear':
			string = '\t\tpiecewise linear, '+str(self._args[0])
			for i in range(self._args[0]):
				string += ',\n\t\t\t'+str(self._series[2*i]*scale)+', '
				string += str(self._series[2*i+1]*scale)
			return string
		elif self.type == 'parabolic':
			string = '\t\tparabolic'
			for arg in self._args[:3]:
				string += ', '+str(arg*scale)
			return string

class Driver(Common):

	types = [
	'File']

	def __init__(self, driverType, database, linking=False):
		self.type = driverType
		self.name = driverType
		self.database = database
		self.users = 0
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'File':
			self._args = [1, 0., 0., 0, '', 1]
		self.name_check(self.database.Driver)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'File':
			string = [
			('Name: ', nval, 0, 30),
			('Stream/Fixed',	args[0], 'Stream driver if pressed, else use a Fixed Step driver'),
			('Initial time:',	args[1], 0., 9.9e10, 'If zero, use simulation start time (Fixed step only)'),
			('Time step:',		args[2], 0., 1.e3, 'If zero, use simulation time step (Fixed step only)'),
			('Pad zeros',		args[3], 'Else use ends of series for out of time bound values '+
				'(Fixed step only)' ),
			('File name:',		args[4], 0, 50, 'If blank, use mbdyn "job name"_"driver name"'),
			('Save file',		args[5], 'If in Stream drive mode, save input values to the File name')]
			if self in self.database.Driver:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Driver)

	def write(self, text):
		if self.type == 'File':
			if self._args[0]:
				text.write('\tfile: '+str(self.database.Driver.index(self))+', stream,\n'+
				'\t\tstream drive name, "BX'+str(self.database.Driver.index(self)).zfill(4)+'",\n'+
				'\t\t\tport, '+str(self.database._port + self.database.Driver.index(self))+',\n'+
				'\t\t\thost, "'+self.database._hostname+'",\n'+
				'\t\tnon blocking,\n'+
				'\t\t'+str(len(self.columns))+';\n')
			else:
				filename, ext= Blender.sys.splitext(self.database.filename)
				if self._args[4]:
					filename += '.' + self._args[4]
				else:
					filename += '.' + self.name.replace(' ', '')
				try:
					tmp = open(filename, 'rt')
				except:
					Draw.PupMenu('Error: Cannot open File Driver file named '+filename)
					return
				lines = tmp.readlines()
				tmp.close()
				columns = lines[1].split()
				string = ('\tfile: '+str(self.database.Driver.index(self))+', fixed step, '+
				str(len(lines))+', '+str(len(columns))+', ')
				if self._args[1]:
					string += str(self._args[1])+', '
				else:
					string += str(self.database._t0)+', '
				if self._args[2]:
					string += str(self._args[2])+', '
				else:
					string += str(self.database._dt)+', '
				if self._args[3]:
					string += 'pad zeros, yes,\n'
				else:
					string += 'pad zeros, no,\n\t\t'
				text.write(string+'"'+filename+'";\n')

class Matrix(Common):

	types = [
	'1x1',
	'3x1',
	'6x1',
	'3x3',
	'6x6',
	'6xN']

	def __init__(self, matrixType, database, linking=False):
		self.type = matrixType
		self.name = matrixType
		self.database = database
		self.users = 0
		self._series = []
		self.links = []
		if self.type == '1x1':
			self._args = [0.]
		elif self.type == '3x1':
			self._args = [0, 0., 0., 0., 0, 1.]
		elif self.type == '6x1':
			self._args = [0]+[0.]*6
		elif self.type == '3x3':
			self._args = [1]+[0]*5
			self._series = [0.]*9
		elif self.type == '6x6':
			self._args = [1]+[0]*5
			self._series = [0.]*36
		elif self.type == '6xN':
			self._args = [1, 0, 2]
			self._series = [0.]*30
		self.name_check(self.database.Matrix)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == '1x1':
			string = [
			('Name: ', nval, 0, 30),
			('Value:',	args[0], -9.9e10, 9.9e10)]
			if self in self.database.Matrix:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == '3x1':
			string = [
			('Name: ', nval, 0, 30),
			('Default',	args[0], 'Override with either default or null values'),
			('x1:',		args[1], -9.9e10, 9.9e10),
			('x2:',		args[2], -9.9e10, 9.9e10),
			('x3:',		args[3], -9.9e10, 9.9e10),
			('Scale',	args[4], 'Scale the vector'),
			('Factor:',	args[5], -9.9e10, 9.9e10)]
			if self in self.database.Matrix:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:6]]
		elif self.type == '6x1':
			string = [
			('Name: ', nval, 0, 30),
			('Default',	args[0], 'Override with either default or null values'),
			('x1:',	args[1], -9.9e10, 9.9e10),
			('x2:',	args[2], -9.9e10, 9.9e10),
			('x3:',	args[3], -9.9e10, 9.9e10),
			('x4:',	args[4], -9.9e10, 9.9e10),
			('x5:',	args[5], -9.9e10, 9.9e10),
			('x6:',	args[6], -9.9e10, 9.9e10)]
			if self in self.database.Matrix:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:7]]
		elif self.type == '3x3':
			state = 0
			while state != 1 and not delete.val:
				string = [
			('Name: ', 	nval, 0, 30),
			('matr',	args[0], 'General matrix (select only one)'),
			('sym',		args[1], 'Symmetric (select only one)'),
			('skew',	args[2], 'Skew symmetric final time (select only one)'),
			('diag',	args[3], 'Diagonal matrix (select only one)'),
			('eye',		args[4], 'Identity matrix (select only one)'),
			('null',	args[5], 'Null matrix (select only one)')]
				if self in self.database.Matrix:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:6])
			self._args = [arg.val for arg in args[:6]]
			if not delete.val:
				N = 3			
				series = [Draw.Create(item) for item in self._series]
				if args[0].val == 1:
					string = []
					for i in range(N):
						string += ['']*2
						for j in range(N):
							ij = i+N*j
							string += [('', series[N*j+i], -9.9e10, 9.9e10)]
						string += ['']*3
					string[0] = 'matr'
				elif args[1].val == 1:
					string = []
					for i in range(N):
						string += ['']*2
						for j in range(i+1):
							ij = i+N*j
							string += [('x('+str(j+1)+','+str(i+1)+'):', series[ij], -9.9e10, 9.9e10)]
						string += ['']*(5-i)
					string[0] = 'sym'
				elif args[2].val == 1:
					string = [('skew')]
					for i,j in enumerate([7,2,3]):
						string += [('v('+str(i+1)+'):', series[j], -9.9e10, 9.9e10)]
				elif args[3].val == 1:
					string = [('diag')]
					for i in range(N):
						string += [('x('+str(i+1)+','+str(i+1)+'):', series[(N+1)*i], -9.9e10, 9.9e10)]
				Draw.PupBlock(self.type, string)
				self._series[:len(series)] = [item.val for item in series]
		elif self.type == '6x6':
			state = 0
			while state != 1 and not delete.val:
				string = [
			('Name: ', 	nval, 0, 30),
			('matr',	args[0], 'General matrix (select only one)'),
			('sym',		args[1], 'Symmetric (select only one)'),
			('diag',	args[2], 'Diagonal matrix (select only one)'),
			('eye',		args[3], 'Identity matrix (select only one)'),
			('null',	args[4], 'Null matrix (select only one)')]
				if self in self.database.Matrix:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:5])
			self._args = [arg.val for arg in args[:5]]
			if not delete.val:
				N = 6			
				series = [Draw.Create(item) for item in self._series]
				if args[0].val == 1:
					string = []
					for i in range(N):
						string += ['']
						for j in range(N):
							ij = i+N*j
							string += [('', series[N*j+i], -9.9e10, 9.9e10)]
						string += ['']
					string[0] = 'matr'
				elif args[1].val == 1:
					string = []
					for i in range(N):
						string += ['']*1
						for j in range(i+1):
							ij = i+N*j
							string += [('x('+str(j+1)+','+str(i+1)+'):', series[ij], -9.9e10, 9.9e10)]
						string += ['']*(6-i)
					string[0] = 'sym'
				elif args[2].val == 1:
					string = [('diag')]
					for i in range(N):
						string += [('x('+str(i+1)+','+str(i+1)+'):', series[(N+1)*i], -9.9e10, 9.9e10)]
				Draw.PupBlock(self.type, string)
				self._series[:len(series)] = [item.val for item in series]
		elif self.type == '6xN':
			state = 0
			while state != 1 and not delete.val:
				string = [
			('Name: ', 	nval, 0, 30),
			('matr',	args[0], 'General matrix (select only one)'),
			('null',	args[1], 'Null matrix (select only one)'),
			('N-cols',	args[2], 2, 5, 'Number of columns')]
				if self in self.database.Matrix:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:2])
			self._args = [arg.val for arg in args[:3]]
			if not delete.val:
				N = 6
				Ncols = args[2].val		
				series = [Draw.Create(item) for item in self._series]
				if args[0].val == 1:
					string = []
					for i in range(Ncols):
						string += ['']
						for j in range(N):
							ij = i+Ncols*j
							string += [('', series[Ncols*j+i], -9.9e10, 9.9e10)]
						string += ['']
					string[0] = 'matr'
				Draw.PupBlock(self.type, string)
				self._series[:len(series)] = [item.val for item in series]
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Matrix)

	def string(self):
		if self.type == '1x1':
			ret = '\n\t\t\t'+str(self._args[0])
		elif self.type == '3x1':
			if self._args[0]:
				ret = '\n\t\t\tdefault'
			else:
				ret = '\n\t\t\t'+str(self._args[1:4]).strip('[]')
				if self._args[4]:
					ret += ', scale, '+str(self._args[5])
		elif self.type == '6x1':
			if self._args[0]:
				ret = '\n\t\t\tdefault'
			else:
				ret = '\n\t\t\t'+str(self._args[1:7]).strip('[]')
		elif self.type == '3x3':
			if self._args[0] == 1:
				ret = ('\n\t\t\tmatr,\n'+
				'\t'*4+str(self._series[0:3]).strip('[]')+',\n'+
				'\t'*4+str(self._series[3:6]).strip('[]')+',\n'+
				'\t'*4+str(self._series[6:9]).strip('[]'))
			elif self._args[1] == 1:
				ret = ('\n\t\t\tsym,\n'+
				'\t'*4+str(self._series[0:3]).strip('[]')+',\n'+
				'\t'*5+str(self._series[4:6]).strip('[]')+',\n'+
				'\t'*6+str(self._series[8]))
			elif self._args[2] == 1:
				ret = '\n\t\t\tskew, '+str([self._series[i] for i in [7, 2, 3]]).strip('[]')
			elif self._args[3] == 1:
				ret = '\n\t\t\tdiag, '+str([self._series[i] for i in [0, 4, 8]]).strip('[]')
			elif self._args[4] == 1:
				ret = '\n\t\t\teye'
			elif self._args[5] == 1:
				ret = '\n\t\t\tnull'
		elif self.type == '6x6':
			if self._args[0] == 1:
				ret = ('\n\t\t\tmatr,\n'+
				'\t'*4+str(self._series[0:6]).strip('[]')+',\n'+
				'\t'*4+str(self._series[6:12]).strip('[]')+',\n'+
				'\t'*4+str(self._series[12:18]).strip('[]')+',\n'+
				'\t'*4+str(self._series[18:24]).strip('[]')+',\n'+
				'\t'*4+str(self._series[24:30]).strip('[]')+',\n'+
				'\t'*4+str(self._series[30:36]).strip('[]'))
			elif self._args[1] == 1:
				ret = ('\n\t\t\tsym,\n'+
				'\t'*4+str(self._series[0:6]).strip('[]')+',\n'+
				'\t'*5+str(self._series[7:12]).strip('[]')+',\n'+
				'\t'*6+str(self._series[14:18]).strip('[]')+',\n'+
				'\t'*7+str(self._series[21:24]).strip('[]')+',\n'+
				'\t'*8+str(self._series[28:30]).strip('[]')+',\n'+
				'\t'*9+str(self._series[35]))
			elif self._args[2] == 1:
				ret = '\n\t\t\tdiag, '+str([self._series[i] for i in [0, 7, 14, 21, 28, 35]]).strip('[]')
			elif self._args[3] == 1:
				ret = '\n\t\t\teye'
			elif self._args[4] == 1:
				ret = '\n\t\t\tnull'
		elif self.type == '6xN':
			if self._args[0] == 1:
				N = self._args[2]
				ret = ('\n\t\t\tmatr')
				for i in range(N):
					ret += ',\n\t'*4+str(self._series[i*N:(i+1)*N]).strip('[]')
			elif self._args[1] == 1:
				ret = '\n\t\t\tnull'
		return ret

class Friction(Common):

	types = [
	'Modlugre',
	'Discrete Coulomb']

	def __init__(self, frictType, database, linking=False):
		self.type = frictType
		self.name = frictType
		self.database = database
		self.users = 0
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'Modlugre':
			self._args = [0., 0., 0., 0., 1, 1, 0.]
			self.links = [None]
		elif self.type == 'Discrete Coulomb':
			self._args = [0, 0., 0, 0., 1, 1, 0.]
			self.links = [None]
		self.name_check(self.database.Friction)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Modlugre':
			tmp = 'select'
			if self.links[0] != None:
				tmp = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			('Sigma 0:',	args[0], -9.9e10, 9.9e10),
			('Sigma 1:',	args[1], -9.9e10, 9.9e10),
			('Sigma 2:',	args[2], -9.9e10, 9.9e10),
			('Kappa:',		args[3], -9.9e10, 9.9e10),
			('sf:'+tmp, 	args[4], 'Change friction function'),
			('Plane hinge',	args[5], 'Simple plane hinge (else just simple shape function)'),
			('Radius:',		args[6], 0., 9.9e10, 'Active only if plane hinge is selected')]
			if self in self.database.Friction:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[4] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select friction function')
					self.select(0, 'Function')
				self._args[4] = 0
		elif self.type == 'Discrete Coulomb':
			tmp = 'select'
			if self.links[0] != None:
				tmp = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			('Use Sigma 2',		args[0]),
			('Sigma 2:',		args[1], 0., 9.9e10),
			('Use vel ratio',	args[2], 'Use velocity ratio'),
			('Vel ratio:',		args[3], 0, 1., 'Velocity ratio'),
			('sf:'+tmp, 		args[4], 'Change friction function'),
			('Plane hinge',		args[5], 'Simple plane hinge (else just simple shape function)'),
			('Radius:',			args[6], 0., 9.9e10, 'Active only if plane hinge is selected')]
			if self in self.database.Friction:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args]
			if not delete.val:
				if self._args[4] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select friction function')
					self.select(0, 'Function')
				self._args[4] = 0
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Friction)

	def string(self):
		if self.frictType == 'Modlugre':
			string = ('modlugre, '+str(self._args[:4]).strip('[]')+',\n'+
			'\t\t\t"'+self._functions[0].name+'", ')
			if self._args[5]:
				string += 'simple plane hinge, '+str(self._args[6])
			else:
				string += 'simple'
			return string
		elif self.frictType == 'Discrete Coulomb':
			string = ('discrete coulomb')
			if self._args[0]:
				string += ',\n\t\t\tsigma2, '+str(self._args[1])
			if self._args[2]:
				string += ',\n\t\t\tvelocity ratio, '+str(self._args[3])
			string += ',\n\t\t\t"'+self._functions[0].name+'", '
			if self._args[5]:
				string += 'simple plane hinge, '+str(self._args[6])
			else:
				string += 'simple'
			return string

class Function(Common):

	types = [
	'Const',
	'Exp',
	'Log',
	'Pow',
	'Linear',
	'Cubic Natural Spline',
	'Multilinear',
	'Chebychev',
	'Sum',
	'Sub',
	'Mul',
	'Div']

	def __init__(self, functType, database, linking=False):
		self.type = functType
		self.name = functType
		self.database = database
		self.users = 0
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'Exp':
			self._args[0] = self._args[2] = 0
		elif self.type == 'Log':
			self._args[0] = self._args[2] = 0
		elif self.type in ['Cubic Natural Spline', 'Multilinear']:
			self._args[0:2] = [0]*2
		elif self.type == 'Chebychev':
			self._args[2:4] = [1, 0]
		elif self.type in ['Sum', 'Sub', 'Mul', 'Div']:
			self._args[0:2] = [1]*2
			self.links = [None]*2
		self.name_check(self.database.Function)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Const':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Exp':
			string = [
			('Name: ', nval, 0, 30),
			('Base',	args[0], 'Overrides default base (e)'),
			('Base:',	args[1], 0., 9.9e10),
			('Coef',	args[2], 'Overrides default coefficient (1)'),
			('Coef:',	args[3], -9.9e10, 9.9e10),
			('Mult:',	args[4], -9.9e10, 9.9e10)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Log':
			string = [
			('Name: ', nval, 0, 30),
			('Base',	args[0], 'Overrides default base (e)'),
			('Base:',	args[1], 0., 9.9e10),
			('Coef',	args[2], 'Overrides default coefficient (1)'),
			('Coef:',	args[3], 0., 9.9e10),
			('Mult:',	args[4], -9.9e10, 9.9e10)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Pow':
			string = [
			('Name: ', nval, 0, 30),
			('Exp:',	args[0], -9.9e10, 9.9e10, 'Exponent coefficient')]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Linear':
			string = [
			('Name: ', nval, 0, 30),
			('x1:',	args[0], -9.9e10, 9.9e10),
			('y1:',	args[1], -9.9e10, 9.9e10),
			('x2:',	args[2], -9.9e10, 9.9e10),
			('y2:',	args[3], -9.9e10, 9.9e10)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type in ['Cubic Natural Spline', 'Multilinear']:
			string = [
			('Name: ', nval, 0, 30),
			('Extrapolate',	args[0], 'Extrapolate'),
			('Number of pts:', args[1], 2, 50, 'Number of points')]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				N = args[1].val
				if 2*N > len(self._series):
					self._series[len(self._series):2*N] = [0. for i in range(len(self._series),2*N)]
				series = [Draw.Create(self._series[i]) for i in range(2*N)]
				seq = ['']*(16*((N-1)/8+1))
				for col in range((N-1)/8+1):
					for i in range(min(8, N-col*8)):
						seq[16*col+i] = ('x'+str(8*col+i)+': ', series[16*col+2*i], -9.9e10, 9.9e10)
						seq[16*col+8+i] = ('y'+str(8*col+i)+': ', series[16*col+2*i+1], -9.9e10, 9.9e10)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [arg.val for arg in args[:2]]
		elif self.type == 'Chebychev':
			string = [
			('Name: ', nval, 0, 30),
			('Lower bound:',	args[0], -9.9e10, 9.9e10),
			('Upper bound:',	args[1], -9.9e10, 9.9e10),
			('Extrapolate',		args[2]),
			('Number of pts:', args[3], 2, 50)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				N = args[3].val
				if N > len(self._series):
					self._series[len(self._series):N] = [0. for i in range(len(self._series),N)]
				series = [Draw.Create(self._series[i]) for i in range(N)]
				seq = []
				for i in range(N):
					seq.append(('coef_'+str(i)+':', series[i], -9.9e10, 9.9e10))
				Draw.PupBlock(self.type, seq)
				self._series[:N] = [item.val for item in series]
				self._args = [arg.val for arg in args[:4]]
		elif self.type in ['Sum', 'Sub', 'Mul', 'Div']:
			tmp = ['select']*2
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					args[i] = Draw.Create(0)
			string = [
			('Name: ', nval, 0, 30),
			('f1:'+tmp[0],	args[0]),
			('f2:'+tmp[1], 	args[1])]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				if args[0].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select f1')
					self.select(0, 'Function', head=head)
				if args[1].val == 1 or self.links[1] == None:
					self.persistantPupMenu('Select f2')
					self.select(1, 'Function', head=head)
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Function)

	def write(self, text):
		if self.written:
			return
		if self.links:
			for link in self.links:
				link.write(text)
		self.written = True
		if self.type == 'Const':
			text.write('scalar function: "'+self.name+'", const, '+str(self._args[0])+';\n')
		elif self.type == 'Exp':
			text.write('scalar function: "'+self.name+'", exp')
			text.write
			if self._args[0]:
				text.write(', base, '+str(self._args[1]))
			if self._args[2]:
				text.write(', coefficient, '+str(self._args[3]))
			text.write(', '+str(self._args[4])+';\n')
		elif self.type == 'Log':
			text.write('scalar function: "'+self.name+'", log')
			if self._args[0]:
				text.write(', base, '+str(self._args[1]))
			if self._args[2]:
				text.write(', coefficient, '+str(self._args[3]))
			text.write(', '+str(self._args[4])+';\n')
		if self.type == 'Pow':
			text.write('scalar function: "'+self.name+'", pow, '+str(self._args[0])+';\n')
		elif self.type == 'Linear':
			text.write('scalar function: "'+self.name+'", linear,\n\t'+
			str(self._args[:4]).strip('[]')+';\n')
		elif self.type == 'Cubic Natural Spline':
			text.write('scalar function: "'+self.name+'", cubicspline')
			if self._args[0]:
				text.write(', do not extrapolate')
			for i in range(self._args[1]):
				text.write(',\n\t'+str(self._series[2*i:2*(i+1)]).strip('[]'))
			text.write(';\n')
		elif self.type == 'Multilinear':
			text.write('scalar function: "'+self.name+'", multilinear')
			if self._args[0]:
				text.write(', do not extrapolate')
			for i in range(self._args[1]):
				text.write(',\n\t'+str(self._series[2*i:2*(i+1)]).strip('[]'))
			text.write(';\n')
		elif self.type == 'Chebychev':
			text.write('scalar function: "'+self.name+'", chebychev,\n\t'+
			str(self._args[:2]).strip('[]'))
			if self._args[2]:
				text.write(', do not extrapolate')
			N = self._args[3]/4
			for i in range(N):
				text.write(',\n\t'+str(self._series[4*i:4*(i+1)]).strip('[]'))
			if 4*N == self._args[3]:
				text.write(';\n')
			else:
				text.write(',\n\t'+str(self._series[4*N:self._args[3]]).strip('[]')+';\n')
		elif self.type == 'Sum':
			text.write('scalar function: "'+self.name+'", sum, "'+
			self.links[0].name+'", "'+self.links[1].name+'";\n')
		elif self.type == 'Sub':
			text.write('scalar function: "'+self.name+'", sub, "'+
			self.links[0].name+'", "'+self.links[1].name+'";\n')
		elif self.type == 'Mul':
			text.write('scalar function: "'+self.name+'", mul, "'+
			self.links[0].name+'", "'+self.links[1].name+'";\n')
		elif self.type == 'Div': 
			text.write('scalar function: "'+self.name+'", div, "'+
			self.links[0].name+'", "'+self.links[1].name+'";\n')

class Drive(Common):

	types = [
	'Null drive',
	'Unit drive',
	'Constant drive',
	'Time drive',
	'Linear drive',
	'Parabolic drive',
	'Cubic drive',
	'Step drive',
	'Double step drive',
	'Ramp drive',
	'Double ramp drive',
	'Piecewise linear drive',
	'Sine drive',
	'Cosine drive',
	'Tanh drive',
	'Fourier series drive',
	'Frequency sweep drive',
	'Exponential drive',
	'Random drive',
	'Meter drive',
	'File drive',
	'String drive',
	'Dof drive',
	'Node drive',
	'Element drive',
	'Drive drive',
	'Array drive',
	'Hints',
	'Template drive']

	def __init__(self, driveType, database, linking=False):
		self.type = driveType
		self.name = driveType
		self.database = database
		self.users = 0
		self.objects = []
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == 'Ramp drive':
			self._args[2] = 0
		elif self.type == 'Double ramp drive':
			self._args[5] = 0
		elif self.type == 'Piecewise linear drive':
			self._args = [0]
		elif self.type == 'Sine drive':
			self._args[3] = 0
			self._args[5:8] = [0]*3
		elif self.type == 'Cosine drive':
			self._args[3] = 0
			self._args[5:8] = [0]*3
		elif self.type == 'Fourier series drive':
			self._args[2] = 1
			self._series = [0.0]
			self._args[3:6]=[0]*3
		elif self.type == 'Frequency sweep drive':
			self._args[1:3] = [1]*2
			self._args[4] = 0
			self.links = [None]*2
		elif self.type == 'Random drive':
			self._args[3] = 0
			self._args[5:7] = [0]*2
		elif self.type == 'Meter drive':
			self._args[1] = self._args[3] = 0
		elif self.type == 'File drive':
			self._args = [1, 0., -1., 1., 0.1, 0, 0]
			self.links = [None]
		elif self.type == 'String drive':
			self._args = 'enter text'
		elif self.type == 'Dof drive':
				self._args = [1, 1]
				self.links = [None]*2
		elif self.type == 'Node drive':
			self._args = [1, '', 1]
			self.links = [None]
			self.objects = [None]
#		elif self.type == 'Element drive':
		elif self.type == 'Drive drive':
			self._args = [1]*2
			self.links = [None]*2
		elif self.type == 'Array drive':
			self._args = [0]
#		elif self.type == 'Hints':
#		elif self.type == 'Template drive':
		sel = Object.GetSelected()
		if sel and self.objects:
			if len(sel) == 2 and len(self.objects) >= 2:
				self.objects[:2] = sel
			else:
				self.objects[0] = sel[0]
		self.name_check(self.database.Drive)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Null drive':
			string = [
			('Name: ', nval, 0, 30)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = []
		elif self.type == 'Unit drive':
			string = [
			('Name: ', nval, 0, 30)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = []
		elif self.type == 'Constant drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Time drive':
			string = [
			('Name: ', nval, 0, 30)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = []
		elif self.type == 'Linear drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10),
			('Slope:',		args[1], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == 'Parabolic drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10),
			('Linear:',		args[1], -9.9e10, 9.9e10),
			('Parabolic:',	args[2], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		elif self.type == 'Cubic drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -9.9e10, 9.9e10),
			('Linear:',		args[1], -9.9e10, 9.9e10),
			('Parabolic:',	args[2], -9.9e10, 9.9e10),
			('Cubic:',		args[3], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Step drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Step:',			args[1], -9.9e10, 9.9e10),
			('Initial value:',	args[2], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		elif self.type == 'Double step drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Final time:',		args[1], 0., 9.9e10),
			('Step value:',		args[2], -9.9e10, 9.9e10),
			('Initial value:',	args[3], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Ramp drive':
			string = [
			('Name: ', nval, 0, 30),
			('Slope:',			args[0], -9.9e10, 9.9e10),
			('Initial time:',	args[1], 0., 9.9e10),
			('Forever',			args[2], 'Overrides final time'),
			('Final time:',		args[3], 0., 9.9e10),
			('Initial value:',	args[4], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Double ramp drive':
			string = [
			('Name: ', nval, 0, 30),
			('Slope-a:',		args[0], -9.9e10, 9.9e10),
			('Initial time:',	args[1], 0., 9.9e10),
			('Final time:',		args[2], 0., 9.9e10),
			('Slope-d:',		args[3], -9.9e10, 9.9e10),
			('Initial time:',	args[4], 0., 9.9e10),
			('Forever',			args[5], 'Overrides final time-d'),
			('Final time:',		args[6], 0., 9.9e10),
			('Initial value:',	args[7], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:8]]
		elif self.type == 'Piecewise linear drive':
			string = [
			('Name: ', nval, 0, 30),
			('Number of pts:', args[0], 2, 50)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				N = args[0].val
				if 2*N > len(self._series):
					self._series[len(self._series):2*N] = [0. for i in range(len(self._series),2*N)]
				series = [Draw.Create(self._series[i]) for i in range(2*N)]
				seq = ['']*(16*((N-1)/8+1))
				for col in range((N-1)/8+1):
					for i in range(min(8, N-col*8)):
						seq[16*col+i] = ('Point_'+str(8*col+i)+': ', series[16*col+2*i], -9.9e10, 9.9e10)
						seq[16*col+8+i] = ('Value_'+str(8*col+i)+': ', series[16*col+2*i+1], -9.9e10, 9.9e10)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [N]
		elif self.type == 'Sine drive':
			state = 2
			while not state in [0,1]:
				string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Omega:',			args[1], 0., 9.9e10, 'radians per second'),
			('Amplitude:',		args[2], -9.9e10, 9.9e10),
			('N-cycles:',		args[3], -1e6, 1e6),
			('Initial value:',	args[4], -9.9e10, 9.9e10),
			('Half',			args[5], 'Overrides final time (select only one override)'),
			('One',				args[6], 'Overrides final time (select only one override)'),
			('Forever',			args[7], 'Overrides final time (select only one override)')]
				if self in self.database.Drive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[5:8])
			self._args = [arg.val for arg in args[:8]]
		elif self.type == 'Cosine drive':
			state = 2
			while not state in [0,1]:
				string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Omega:',			args[1], 0., 9.9e10, 'radians per second'),
			('Amplitude:',		args[2], -9.9e10, 9.9e10),
			('N-cycles:',		args[3], -1e6, 1e6),
			('Initial value:',	args[4], -9.9e10, 9.9e10),
			('Half',			args[5], 'Overrides final time (select only one override)'),
			('One',				args[6], 'Overrides final time (select only one override)'),
			('Forever',			args[7], 'Overrides final time (select only one override)')]
				if self in self.database.Drive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[5:8])
			self._args = [arg.val for arg in args[:8]]
		elif self.type == 'Tanh drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Amplitude:',		args[1], -9.9e10, 9.9e10),
			('Slope:',			args[2], -9.9e10, 9.9e10),
			('Initial value:',	args[3], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Fourier series drive':
			state = 2
			while not state in [0,1]:
				string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Omega:',			args[1], -9.9e10, 9.9e10),
			('Number of terms:',args[2], 1, 50),
			('N-cycles:',		args[3], -1e6, 1e6),
			('One',				args[4], 'Overrides N-cycles (select only one override)'),
			('Forever',			args[5], 'Overrides N-cycles (select only one override)'),
			('Initial value:',	args[6], -9.9e10, 9.9e10)]
				if self in self.database.Drive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[4:6])
			if not delete.val:
				N = args[2].val + 1
				if N > (len(self._series)+1)/2:
					self._series.extend([0.]*2*(N-self._args[2]))
				self._args = [arg.val for arg in args]
				series = [Draw.Create(item) for item in self._series]
				seq = ['']*(16*((N-1)/8+1))
				for col in range((N-1)/8+1):
					for i in range(min(8, N-col*8)):
						seq[16*col+i] = ('a_'+str(8*col+i)+': ', series[16*col+2*i-1], -9.9e10, 9.9e10)
						seq[16*col+8+i] = ('b_'+str(8*col+i)+': ', series[16*col+2*i], -9.9e10, 9.9e10)
				seq[0] = ('Using only terms 0 to '+str(N-1))
				seq[8] = ('a_0: ', series[0], -9.9e10, 9.9e10)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
		elif self.type == 'Frequency sweep drive':
			tmp = ['select drive']*2
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			('Initial time: ',	args[0], 0., 9.9e10),
			(tmp[0],		 	args[1], 'Change angular velocity drive'),
			(tmp[1], 			args[2], 'Change amplitude drive.'),
			('Initial value:',	args[3], -9.9e10, 9.9e10),
			('Forever',			args[4], 'Overrides final time'),
			('Final time:',		args[5], 0., 9.9e10),
			('Final value:',	args[6], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:7]]
				if self._args[1] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select an angular velocity drive:')
					self.select(0, 'Drive', head=head)
				if self._args[2] == 1 or self.links[1] == None:
					self.persistantPupMenu('Select an amplitude drive:')
					self.select(1, 'Drive', head=head)
				self._args[1:3] = [0]*2
		elif self.type == 'Exponential drive':
			string = [
			('Name: ', nval, 0, 30),
			('Amplitude:',		args[0], -9.9e10, 9.9e10),
			('Time constant:',	args[1], -9.9e10, 9.9e10),
			('Initial time:',	args[2], 0., 9.9e10),
			('Initial value:',	args[3], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Random drive':
			string = [
			('Name: ', nval, 0, 30),
			('Amplitude:',		args[0], -9.9e10, 9.9e10),
			('Mean:',			args[1], -9.9e10, 9.9e10),
			('Initial time:',	args[2], 0., 9.9e10),
			('Forever',			args[3], 'Overrides final time'),
			('Final time:',		args[4], 0., 9.9e10),
			('Steps:',			args[5], 0, 1e6, 'Steps to hold value'),
			('Time seed',		args[6], 'Overrides Seed value with Time value'),
			('Seed:',			args[7], -9.9e10, 9.9e10)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:8]]
		elif self.type == 'Meter drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 9.9e10),
			('Forever',			args[1], 'Overrides final time'),
			('Final time:',		args[2], 0., 9.9e10),
			('Steps:',			args[3], 0, 1e6, 'Steps between spikes')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'File drive':
			tmp = 'select driver'
			if self.links[0] != None:
					tmp = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp,				args[0], 'Change driver'),
			('Initial value:',	args[1], -9.9e10, 9.9e10),
			('Minimum value:',	args[2], -9.9e10, 9.9e10),
			('Maximum value:',	args[3], -9.9e10, 9.9e10),
			('Increment:',		args[4], -9.9e10, 9.9e10),
			('Joystick',		args[5]),
			('Axis:',			args[6], 0, 5, 'Choose joystick axis')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				if args[0].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select driver:')
					self.select(0, 'Driver', clas_types=['File'])
			self._args = [0] + [arg.val for arg in args[1:]]
		elif self.type == 'String drive':
			string = [
			('Name: ', nval, 0, 30)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = Draw.PupStrInput('String: ', self._args, 100)
		elif self.type == 'Dof drive':
			link_name = ['select node', 'select drive']
			for i in range(2):
				if self.links[i] != None:
					link_name[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(link_name[0], 	args[0], 'Change abstract node'),
			(link_name[1], 	args[1], 'Change drive')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				if args[0].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select an abstract node')
					self.select(0, 'NS_Node', clas_types=['Abstract'])
				if args[1].val == 1 or self.links[1] == None:
					self.persistantPupMenu('Select a drive:')
					self.select(1, 'Drive', head=head)
				self._args = [0, 0]
		elif self.type == 'Node drive':
			obj_name = 'select node'
			if self.objects[0]:
				obj_name = self.objects[0].name
				args[0] = Draw.Create(0)
			tmp = 'select drive'
			if self.links[0] != None:
				tmp = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			(obj_name,		args[0], 'Change node'),
			('Reference:',	args[1], 0, 10, 'Private data of the node'),
			(tmp, 			args[2], 'Change drive.')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				title = 'Select node:'
				if args[0].val == 1 or self.objects[0] == None:
					self.persistantPupMenu(title)
					self.select(0, 'Object')
				if args[2].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select a drive:')
					self.select(0, 'Drive', head=head)
				self._args = [0, args[1].val, 0]
#		elif self.type == 'Element drive':
		elif self.type == 'Drive drive':
			tmp = ['select drive']*2
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change drive 1'),
			(tmp[1], 	args[1], 'Change drive 2')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				for i, arg in enumerate(args):
					if arg.val == 1 or self.links[i] == None:
						self.persistantPupMenu('Select drive '+str(i+1)+':')
						self.select(i, 'Drive', head=head)
				self._args[:2] = [0]*2
		elif self.type == 'Array drive': 
			toggles = [(drive.name, arg, 'Remove this drive from the list') for
				drive, arg in zip(self.links, args[1:])]
			string = [
			('Name: ', 				nval, 0, 30),
			('Number of drives:',	args[0], 1, 20, 'Number of drives in the array')]
			string += toggles
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				toBeDeleted = set([i for i in range(len(args)-1) if args[i+1].val == 1])
				if args[0].val < self._args[0]:
					toBeDeleted |= set(range(args[0].val, self._args[0]))
				toBeDeleted = list(toBeDeleted)
				toBeDeleted.sort(reverse=True)
				for i in toBeDeleted:
					self.links[i].users -= 1
					del self.links[i]
					del self._args[i+1]
				for i in range(len(self.links), args[0].val):
					self.persistantPupMenu('Select drive '+ str(i))
					self.links.append(None)
					self._args.append(1)
					self.select(i, 'Drive', head=head)
				self._args = [args[0].val] + [0]*args[0].val

#		elif self.type == 'Hints':
#		elif self.type == 'Template drive':
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Drive)

	def string(self, owner=None):
		if self.type == 'Null drive':
			return 'null'
		elif self.type == 'Unit drive':
			return 'unit'
		elif self.type == 'Constant drive':
			return str(self._args[0])
		elif self.type == 'Time drive':
			return 'time'
		elif self.type == 'Linear drive':
			return 'linear, '+str(self._args).strip('[]')
		elif self.type == 'Parabolic drive':
			return 'parabolic, '+str(self._args).strip('[]')
		elif self.type == 'Cubic drive':
			return 'cubic, '+str(self._args).strip('[]')
		elif self.type == 'Step drive':
			return 'step, '+str(self._args).strip('[]')
		elif self.type == 'Double step drive':
			return 'double step, '+str(self._args).strip('[]')
		elif self.type == 'Ramp drive':
			string = 'ramp, '+str(self._args[:2]).strip('[]')+', '
			if self._args[2]:
				string += 'forever, '
			else:
				string += str(self._args[3])+', '
			return string + str(self._args[4])
		elif self.type == 'Double ramp drive':
			string = 'double ramp, '+str(self._args[:5]).strip('[]')+', '
			if self._args[5]:
				string += 'forever, '
			else:
				string += str(self._args[6])+', '
			return string + str(self._args[7])
		elif self.type == 'Piecewise linear drive':
			N = self._args[0]
			string = 'piecewise linear, '+str(N)
			for i in range(N):
				string += ',\n\t\t'+str(self._series[2*i: 2*(i+1)]).strip('[]')
			return string
		elif self.type == 'Sine drive':
			string = 'sine, '+str(self._args[:3]).strip('[]')+', '
			if self._args[5]:
				string += 'half, '
			elif self._args[6]:
				string += 'one, '
			elif self._args[7]:
				string += 'forever, '
			else:
				string += str(self._args[3])+', '
			return string + str(self._args[4])			
		elif self.type == 'Cosine drive':
			string = 'cosine, '+str(self._args[:3]).strip('[]')+', '
			if self._args[5]:
				string += 'half, '
			elif self._args[6]:
				string += 'one, '
			elif self._args[7]:
				string += 'forever, '
			else:
				string += str(self._args[3])+', '
			return string + str(self._args[4])			
		elif self.type == 'Tanh drive':
			return 'tanh, '+str(self._args).strip('[]')
		elif self.type == 'Fourier series drive':
			string = ('fourier series, '+str(self._args[:3]).strip('[]')+
			'\n\t\t'+str(self._series[0]))
			N = self._args[2] + 1
			for i in range(1, N):
				string += ',\n\t\t'+str(self._series[2*i-1: 2*i+1]).strip('[]')
			if self._args[4]:
				string += ',\n\t\tone, '
			elif self._args[5]:
				string += ',\n\t\tforever, '
			else:
				string += ',\n\t\t'+str(self._args[3])+', '
			return string + str(self._args[6])			
		elif self.type == 'Frequency sweep drive':
			string = ('frequency sweep, '+str(self._args[0])+',\n'+
			'\t\t'+self.links[0].string(self)+',\n'+
			'\t\t'+self.links[1].string(self)+',\n'+
			'\t\t'+str(self._args[3])+', ')
			if self._args[4]:
				string += 'forever, '
			else:
				string += str(self._args[5])+', '
			return string + str(self._args[6])
		elif self.type == 'Exponential drive':
			return 'exponential, '+str(self._args).strip('[]')
		elif self.type == 'Random drive':
			string = 'random, '+str(self._args[:3]).strip('[]')+', '
			if self._args[3]:
				string += 'forever, '
			else:
				string += str(self._args[4])+', '
			string += 'steps, '+str(self._args[5])+', seed, '
			if self._args[6]:
				string += 'time'
			else:
				string += str(self._args[7])
			return string
		elif self.type == 'Meter drive':
			string = 'meter, '+str(self._args[0])+', '
			if self._args[1]:
				string += 'forever, '
			else:
				string += str(self._args[2])+', '
			return string+'steps, '+str(self._args[3])
		elif self.type == 'File drive':
			return ('file, '+str(self.database.Driver.index(self.links[0]))+', '+str(1+self.links[0].columns.index(self)))
		elif self.type == 'String drive':
			return 'string, "'+self._args+'"'
		elif self.type == 'Dof drive':
			string = 'dof, '+str(self.database.NS_Node.index(self.links[0]))+', abstract'
			if self.links[0]._args[1]:
				string += ', differential'
			else:
				string += ', algebraic'
			string += ',\n\t\t\t'+self.links[1].string(self)
			return string
		elif self.type == 'Node drive':
			if self.objects[0] in self.database.Node:
				ob = self.objects[0]
			elif self.database.rigid_dict.has_key(self.objects[0]):
				ob = self.database.rigid_dict[self.objects[0]]
			else:
				self.persistantPupMenu('Error: Node drive '+self.name+' is not assigned to a Node')
				print 'Error: Node drive '+self.name+' is not assigned to a Node'
				return
			string = 'node, '+str(self.database.Node.index(ob))
			if ob in self.database.structural_dynamic_nodes | self.database.structural_static_nodes:
				string += ', structural, string, "'+self._args[1]+'",\n'
			return (string+'\t\t'+self.links[0].string(self))
		elif self.type == 'Element drive': pass
		elif self.type == 'Drive drive':
			return ('drive,\n'+
			'\t\t\t'+self.links[0].string(self)+',\n'+
			'\t\t\t'+self.links[1].string(self))
		elif self.type == 'Array drive':
			string = ('array, '+str(self._args[0]))
			for i in range(self._args[0]):
				string += ',\n\t\t'+self.links[i].string(self)
			return string
		elif self.type == 'Hints': pass
		elif self.type == 'Template drive': pass

	def rigid_offset(self, i):
		if self.objects[i] in self.database.Node:
			ob = self.objects[i]
		elif self.database.rigid_dict.has_key(self.objects[i]):
			ob = self.database.rigid_dict[self.objects[i]]
		else:
			print 'Model Error: Object '+self.objects[i].name+' is not associated with a Node'
		rot = ob.getMatrix().toQuat().toMatrix()
		globalV = Vector(self.objects[i].loc) - Vector(ob.loc)
		return rot, globalV, self.database.Node.index(ob)


