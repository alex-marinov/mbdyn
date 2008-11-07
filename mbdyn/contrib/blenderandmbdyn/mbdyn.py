#!BPY

__author__ = "G. Douglas Baldwin, douglasbaldwin AT verizon.net"
__url__ = ["http://www.baldwintechnology.com"]
__version__ = "0.1.1"
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
# Copyright (C) 2008 G. Douglas Baldwin - http://www.baldwintechnology.com
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
			for link in new.links:
				link.users += 1
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

	entity_classes = ['Element', 'Constitutive', 'Drive', 'Friction', 'Shape', 'Function', 'Matrix']

	def __init__(self):
		self.Element = []
		self.Constitutive = []
		self.Drive = []
		self.Shape = []
		self.Function = []
		self.Friction = []
		self.Matrix = []
		self.Frame = []
		self.filename = ''
		self.defaults()

	def defaults(self):
		self._integrator = 'multistep'
		self._t0 = 0.
		self._tN = 0.
		self._dt = 1.e-3
		self._maxI = 10
		self._tol = 1.e-6
		self._dTol = 2.
		self._dC = 1.e-3
		self._fps = 30
		self._IPOs = 1
		self._reinit = 0

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
		reinit = Draw.Create(self._reinit)

		string = [
		('int', integrator, 0, 12, 'integrator'),
		('t0', t0, 0.0, 1.0e6, 'initial time'),
		('tN', tN, 0.0, 1.0e6, 'final time (tN==0. => Use Blender animation range)'),
		('dt*1.e3', dt_1e3, 1.0e-3, 1.0e3, 'time step'),
		('maxI', maxI, 1, 100, 'max iterations'),
		('tol*1.e6', tol_1e6, 1.0e-3, 1.0e3, 'tolerance'),
		('dTol', dTol, 1.0e-3, 1.0e3, 'deriatives tolerance'),
		('dC*1.e3', dC_1e3, 1.0e-3, 1.0e3, 'derivatives coefficient'),
		('fps', fps, 1, 1000, 'frames per second: used to format output'),
		('IPOs', IPOs, 'Create new IPs, else overwrite existing IPOs'),
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
		self._reinit = reinit.val

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
		structural_dynamic_nodes = set([])
		structural_static_nodes = set([])
		for clas in ['Element', 'Drive']:
			for entity in eval('self.'+clas):
				try:
					if entity.objects[0] not in rigids:
						nodes |= set([entity.objects[0]])
						if entity.type in Element.structural_dynamic:
							structural_dynamic_nodes |= set([entity.objects[0]])
						elif entity.type in Element.structural_static:
							structural_static_nodes |= set([entity.objects[0]])
					elif entity.objects[0] in rigids:
						ob = self.rigid_dict[entity.objects[0]]
						nodes |= set([ob])
						if entity.type in Element.structural_dynamic:
							structural_dynamic_nodes |= set([ob])
						elif entity.type in Element.structural_static:
							structural_static_nodes |= set([ob])
				except:
					pass
		nodes -= set([None])
		structural_static_nodes -= set([None])
		structural_dynamic_nodes -= set([None])
		structural_static_nodes -= structural_dynamic_nodes
		structural_node_count = len(structural_static_nodes | structural_dynamic_nodes)
		self.Node = list(nodes)
		self.Node.sort(key = lambda x: x.name)
		self.Frame.sort(key = lambda x: x.objects[0].name)
		self.Element = ([e for e in self.Element if e.name == 'Rotor']+
			[e for e in self.Element if e.name != 'Rotor'])
		frame_dict = {}
		for i, frame in enumerate(self.Frame):
			for ob in frame.objects[1:]:
				frame_dict[ob] = str(i)

		text.write(
		'\n/* Label Indexes\n')
		if self.Frame:
			text.write('\nReference:\n')
			for i, frame in enumerate(self.Frame):
				text.write('\t'+str(i)+'\t- '+frame.objects[0].name+'\n')
		for clas in ['Constitutive', 'Drive', 'Node', 'Element']:
			if eval('self.'+clas):
				text.write('\n'+clas+'s:\n')
				for i, entity in enumerate(eval('self.'+clas)):
					text.write('\t'+str(i)+'\t- '+entity.name+'\n')
		text.write('\n*/\n\n')

		text.write(
		'begin: data;\n'+
		'\tintegrator: '+self._integrator+';\n'+
		'end: data;\n\n')

		finalTime = self._tN
		if finalTime == 0.:
			finalTime = self._t0 + float(Blender.Get('endframe') - Blender.Get('staframe')
				)/float(self._fps)

		text.write(
		'begin: '+self._integrator+';\n'+
		'\tinitial time: '+str(self._t0)+';\n'+
		'\tfinal time: '+str(finalTime)+';\n'+
		'\ttime step: '+str(self._dt)+';\n'+
		'\tmax iterations: '+str(self._maxI)+';\n'+
		'\ttolerance: '+str(self._tol)+';\n'+
		'\tderivatives tolerance: '+str(self._dTol)+';\n'+
		'\tderivatives coefficient: '+str(self._dC)+';\n'+
		'end: '+str(self._integrator)+';\n\n'
		)

		text.write(
		'begin: control data;\n'+
		'\tdefault orientation: orientation matrix;\n'+
		'\toutput meter: meter, 0., forever, steps, '+
			str(int(round(1./(self._fps * self._dt))))+';\n')

		joint_count = 0
		rigid_body_count = 0
		air_properties = False
		gravity = False
		aerodynamic_element_count = 0
		rotor_count = 0
		for element in self.Element:
			cat = element.category()
			if cat == 'joint':  joint_count += 1
			elif cat == 'rigid body': rigid_body_count += 1
			elif cat == 'air properties': air_properties = True
			elif cat == 'aerodynamic element': aerodynamic_element_count += 1
			elif cat == 'gravity': gravity = True
			elif cat == 'rotor': rotor_count += 1

		if structural_node_count:
			text.write('\tstructural nodes: '+str(structural_node_count)+';\n')
		if joint_count:
			text.write('\tjoints: '+str(joint_count)+';\n')
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
			globalV = rot*localV
			if frame.links[1]._args[0]:
				localO = Mathutils.Vector(0., 0., 0.)
			else:
				localO = Mathutils.Vector([frame.links[1]._args[i] for i in range(1,4)])
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

		if self.Drive:
			text.write('\n')
		for i, drive in enumerate(self.Drive):
			text.write('drive caller: '+str(i)+',\n\t'+drive.string()+';\n')

#		if self.Constitutive:
#			text.write('\n')
#		for i, const in enumerate(self.Constitutive):
#			text.write('constitutive law: '+str(i)+',\n\t'+const.string()+';\n')

		if self.Function:
			text.write('\n')
		for function in self.Function:
			function.write(text)

		text.write('\nbegin: nodes;\n')

		for i, node in enumerate(self.Node):
			if node in structural_dynamic_nodes:
				text.write('\tstructural: '+str(i)+', dynamic,\n')
			elif node in structural_static_nodes:
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

		text.write('end: nodes;\n\n')
		text.write('begin: elements;\n')

		for element in self.Element:
			element.write(text)

		text.write('end: elements;\n')

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
				ret = Draw.PupName('Warning: Click here to rename object '+ob.name+' to be __'+ob.name)
				if ret:
					return
				ob.name = '__'+ob.name

		for ob in self._selected:
			if ob.type != 'Mesh':
				Draw.PupMenu('Error: Must select only Blender Mesh objects. |'+
					'Your selection included a '+ob.type+' object.')
				return

		if clas_types:
			for typ in clas_types:
				hold = [(entity.name, []) for entity in eval('self.database.'+self.classes[0]) if
					entity.type == typ]
				hold.sort()
				self.string += [(typ, [('New', [])]+hold)]
				self._assignevents(self.string)
				self._groups[typ] = self._event
			return

		for clas in self.classes:

			if clas == 'Matrix':
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

		if 'Object' in self.classes and sel < self._groups['Object']:
			for ob in scn.objects:
				if ob.name == self.events[sel]:
					return ob

		if 'All Elements' in self.classes and sel < self._groups['All Elements']:
			for e in self.database.Element:
				if e.name == self.events[sel]:
					return e

		for clas in [clas for clas in self.classes if clas not in ['Matrix']]:
			if self._groups.has_key('New '+clas) and sel < self._groups['New '+clas]:
				return eval(clas+'(self.events[sel], self.database, linking)')
			elif self._groups.has_key(clas) and sel < self._groups[clas]:
				for entity in eval('self.database.'+clas):
					if entity.name == self.events[sel]:
						return entity.modify(linking, self.head)

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

		coords = []
		scale = 1.
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

class Element(Common):

	types = [
	'Rotor',
	'Aerodynamic body',
	'Body',
	'Joint',
	'Gravity',
	'Air properties',
	'Driven']

	Joint = [
	'Axial rotation',
	'Clamp',
	'Distance',
	'Deformable hinge',
	'Revolute hinge',
	'Rod',
	'Spherical hinge']

	aerodynamic = [
	'Aerodynamic body']

	rigid_body = [
	'Body']

	structural_static = aerodynamic + ['Joint', 'Rotor']

	structural_dynamic = rigid_body

	def category(self):
		if self.type == 'Joint':
			return 'joint'
		elif self.type in Element.rigid_body:
			return 'rigid body'
		elif self.type == 'Air properties':
			return 'air properties'
		elif self.type in Element.aerodynamic:
			return 'aerodynamic element'
		elif self.type == 'Gravity':
			return 'gravity'
		elif self.type == 'Rotor':
			return 'rotor'

	def __init__(self, elemType, database, linking=False):
		if elemType in Element.Joint:
			self.type = 'Joint'
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
		elif self.type == 'Aerodynamic body':
			self._args = [1, 0, 0, 1., 1, 1, 1, 1, 1]
			self.objects = [None]
			self.links = [None]*5
		elif self.type == 'Body':
			self._args = [1, 1., 1]
			self.objects = [None]
			self.links = [None]
		elif self.type == 'Joint':
			if self.subtype == 'Axial rotation':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Clamp':
				self._args = [1]
				self.objects = [None]
			elif self.subtype == 'Deformable hinge':
				self._args = [1, 1, 1]
				self.objects = [None]*2
				self.links = [None]
			elif self.subtype == 'Distance':
				self._args = [1, 1, 1, 0]
				self.objects = [None]*2
				self.links = [None]
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
		elif self.type == 'Gravity':
			self._args = [1, 1]
			self.links = [None]*2
		elif self.type == 'Air properties':
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
		sel = Object.GetSelected()
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
			('Ref omega:',			args[7], 1., 1.e6),
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
		elif self.type == 'Aerodynamic body':
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
			('Mass:',				args[1], 0., 1.e6),
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
						self.select(i, 'Object')
					if self._args[1] or not self.objects[1]:
						self.persistantPupMenu('Select attached node:')
						self.select(i, 'Object', head=self.objects[0])
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
						self.select(0, 'Constitutive')
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
				('Initial theta:',		args[3], 0., 1.e6),
				('Enable friction',		args[4], 'Enable friction model'),
				('Avg radius:',			args[5], 0., 1.e6, 'Used in friction model'),
				('Enable preload',		args[6], 'Enable preload'),
				('Preload:',			args[7], 0., 1.e6, 'Used in friction model'),
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
						self.select(0, 'Constitutive')
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
				text.write(',\n\t\tdelay, reference, '+str(self.database.Drive.index(self.links[0])))
			text.write(',\n\t\tmax iterations, '+str(self._args[17])+
			',\n\t\ttolerance, '+str(self._args[18])+
			',\n\t\teta, '+str(self._args[19])+
			',\n\t\tcorrection, '+str(self._args[20:22]).strip('[]')+';\n')
		elif self.type == 'Aerodynamic body':
			rot, globalV, iNode = self.rigid_offset(0)
			localV = rot*globalV
			rot0 = self.objects[0].getMatrix().toQuat().toMatrix()
			text.write('\taerodynamic body: '+str(self.database.Element.index(self))+', '+
			str(iNode))
			if self._args[1]:
				text.write(', rotor, '+str(self.database.Element.index(self.links[0])))
			text.write(',\n\t\t'+str(localV[0])+', '+str(localV[1])+', '+str(localV[2])+
			',\n\t\tmatr,\n')
			self.rotationMatrix_write(rot0.transpose()*rot, text, '\t\t\t')
			text.write(',\n\t\t'+str(self._args[3])+',\n')
			for shape in self.links[1:4]:
				text.write(shape.string()+',\n')
			text.write(self.links[4].string(3.14159/180.)+',\n')
			text.write('\t\t'+str(self._args[8])+';\n')
		elif self.type == 'Body':
			rot, globalV, iNode = self.rigid_offset(0)
			localV = rot*globalV
			rot0 = self.objects[0].getMatrix().toQuat().toMatrix()
			text.write('\tbody: '+str(self.database.Element.index(self))+', '+
			str(iNode)+', '+str(self._args[1])+',\n'+
			'\t\t'+str(localV[0])+', '+str(localV[1])+', '+str(localV[2])+', '+
			self.links[0].string()+',\n\t\tinertial, matr,\n')
			self.rotationMatrix_write(rot0.transpose()*rot, text, '\t\t\t')
			text.write(';\n')
		elif self.type == 'Joint':
			if self.subtype == 'Axial rotation':
				self.write_hinge(text, 'axial rotation')
				text.write(',\n\t\treference, '+str(self.database.Drive.index(self.links[0]))+';\n')
			elif self.subtype == 'Clamp':
				text.write(
				'\tjoint: '+str(self.database.Element.index(self))+', clamp,\n'+
				'\t\t'+str(self.database.Node.index(self.objects[0]))+', node, node;\n')
			elif self.subtype == 'Deformable hinge':
				self.write_hinge(text, 'deformable hinge', V1=False, V2=False)
				text.write(',\n\t\t'+self.links[0].string()+';\n')
			elif self.subtype == 'Distance':
				self.write_offset(text, 'distance')
				if self._args[2]:
					text.write(', from nodes;\n')
				else:
					text.write(', reference, '+str(self.database.Drive.index(self.links[0]))+';\n')
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
				self.write_offset(text, 'rod')
				text.write(',\n\t\tfrom nodes, '+self.links[0].string()+';\n')
			elif self.subtype == 'Spherical hinge':
				self.write_hinge(text, 'spherical hinge')
				text.write(';\n')
		elif self.type == 'Gravity':
			text.write('\tgravity: '+self.links[0].string()+', reference, '+
				str(self.database.Drive.index(self.links[1]))+';\n')
		elif self.type == 'Air properties':
			if self._args[0]:
				text.write('\tair properties: std, SI,\n')
			elif self._args[0]:
				text.write('\tair properties: std, British,\n')
			text.write('\t\ttemperature deviation, '+str(self._args[2])+
				',\n\t\treference altitude, '+str(self._args[3])+
				',\n\t\t'+self.links[0].string()+', reference, '+
				str(self.database.Drive.index(self.links[1])))
			if self._args[5]:
				text.write(',\n\t\tgust, front 1D,\n\t\t'+
					self.links[2].string()+',\n\t\t\t'+
					self.links[3].string()+',\n\t\t\t'+
					str(self._args[9])+', reference, '+str(self.database.Drive.index(self.links[4])))
			text.write(';\n')
		elif self.type == 'Driven':
			text.write('\tdriven: '+str(self.database.Element.index(self.links[1]))+', '+
			'reference, '+str(self.database.Drive.index(self.links[0]))+',\n'+
			'\t\texisting: '+self.links[1].type+', '+str(self.database.Element.index(self.links[1]))+';\n')

	def write_hinge(self, text, name, V1=True, V2=True, M1=True, M2=True):
		rot0, globalV0, iNode0 = self.rigid_offset(0)
		localV0 = rot0*globalV0
		rot1, globalV1, iNode1 = self.rigid_offset(1)
		to_hinge = rot1*(globalV1 + Vector(self.objects[0].loc) - Vector(self.objects[1].loc))
		rot = self.objects[0].getMatrix().toQuat().toMatrix().transpose()
		text.write(
		'\tjoint: '+str(self.database.Element.index(self))+', '+name+',\n'+
		'\t\t'+str(iNode0))
		if V1:
			text.write(', '+str(localV0[0])+', '+str(localV0[1])+', '+str(localV0[2]))
		if M1:
			text.write(',\n\t\t\thinge, matr,\n')
			self.rotationMatrix_write(rot0*rot, text, '\t\t\t\t')
		text.write(', \n\t\t'+str(iNode1))
		if V2:
			text.write(', '+str(to_hinge[0])+', '+str(to_hinge[1])+', '+str(to_hinge[2]))
		if M2:
			text.write(',\n\t\t\thinge, matr,\n')
			self.rotationMatrix_write(rot1*rot, text, '\t\t\t\t')

	def write_offset(self, text, name):
		rot0, globalV0, iNode0 = self.rigid_offset(0)
		localV0 = rot0*globalV0
		rot1, globalV1, iNode1 = self.rigid_offset(1)
		localV1 = rot1*globalV1
		text.write(
		'\tjoint: '+str(self.database.Element.index(self))+', '+name+',\n'+
		'\t\t'+str(iNode0)+', position, '+
		str(localV0[0])+', '+str(localV0[1])+', '+str(localV0[2])+
		', \n\t\t'+str(iNode1)+', position, '+
		str(localV1[0])+', '+str(localV1[1])+', '+str(localV1[2]))

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

		if self.type == 'Aerodynamic body':
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
			for i in range(8):
				me.faces[i].smooth = 1
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
			for f in me.faces: f.smooth = 1
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
				for i in range(2,6): me.faces[i].smooth = 1
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
				for f in me.faces: f.smooth = 1
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
				for f in me.faces: f.smooth = 1
				for i in range(4,8): me.edges[i].crease = 255
			else:
				me = saveData
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
	'Linear elastic',
	'Linear elastic generic',
	'Linear elastic generic axial torsion coupling',
	'Cubic elastic generic',
	'Log elastic',
	'Linear elastic bi-stop generic',
	'Double linear elastic',
	'Isotropic hardening elastic',
	'Scalar function elastic',
	'Scalar function elastic orthotropic',
	'Linear viscous',
	'Linear viscous generic',
	'Linear viscoelastic',
	'Linear viscoelastic generic',
	'Linear viscoelastic generic axial torsion coupling',
	'Cubic viscoelastic generic',
	'Double linear viscoelastic',
	'Turbulent viscoelastic',
	'Linear viscoelastic bi-stop generic',
	'GRAALL damper',
	'shock absorber',
	'symbolic elastic',
	'symbolic viscous',
	'symbolic viscoelastic',
	'ann elastic',
	'ann viscoelastic',
	'nlsf viscoelastic',
	'nlp viscoelastic',
	'invariant angular ']

	def __init__(self, constType, database, linking=False):
		self.type = constType
		self.name = constType
		self.database = database
		self._args = [0.]*10
		self.users = 0
		self.links = []
		if self.type == 'Linear elastic':
			self._args[:3] = [1, 0, 0]
		elif self.type == 'Linear elastic generic':
			self._args[0] = 1
			self.links = [None]
		elif self.type == 'Linear elastic generic axial torsion coupling':
			self._args[0] = 1
			self.links = [None]
		elif self.type == 'Cubic elastic generic':
			self._args[0:3] = [1]*3
			self.links = [None]*3
		elif self.type == 'Log elastic': pass
		elif self.type == 'Linear elastic bi-stop generic':
			self._args = [1, 0, 0, 1, 1]
			self.links = [None]*3
		elif self.type == 'Double linear elastic':
			self._args[0:2] = [1, 0]
		elif self.type == 'Isotropic hardening elastic':
			self._args[0:3] = [1, 0, 0]
			self._args[5] = 0
		elif self.type == 'Scalar function elastic':
			self._args[0:3] = [1, 0, 0]
			self._args[3] = 1
			self._functions = [None]
		elif self.type == 'Scalar function elastic orthotropic':
			self._args = [1, 0, 0] + [1]
			self._functions = [None]
		elif self.type == 'Linear viscous':
			self._args[:3] = [1, 0, 0]
		elif self.type == 'Linear viscous generic':
			self._args[0] = 1
			self.links = [None]
		elif self.type == 'Linear viscoelastic':
			self._args[:3] = [1, 0, 0]
			self._args[5] = 0
		elif self.type == 'Linear viscoelastic generic':
			self._args[:3] = [1, 1, 0]
			self.links = [None]*2
		elif self.type == 'Linear viscoelastic generic axial torsion coupling':
			self._args[:3] = [1, 1, 0]
			self.links = [None]*2
		elif self.type == 'Cubic viscoelastic generic':
			self._args[0:4] = [1]*4
			self.links = [None]*4
		elif self.type == 'Double linear viscoelastic':
			self._args[0:2] = [1, 0]
			self._args[7] = 0
		elif self.type == 'Turbulent viscoelastic':
			self._args[2] = 0
			self._args[4] = 0
		elif self.type == 'Linear viscoelastic bi-stop generic':
			self._args = [1, 1, 0, 0, 1, 1]
			self.links = [None]*4
			"""
		elif self.type == 'GRAALL damper':
		elif self.type == 'shock absorber':
		elif self.type == 'symbolic elastic':
		elif self.type == 'symbolic viscous':
		elif self.type == 'symbolic viscoelastic':
		elif self.type == 'ann elastic':
		elif self.type == 'ann viscoelastic':
		elif self.type == 'nlsf viscoelastic':
		elif self.type == 'nlp viscoelastic':
		elif self.type == 'invariant angular ':
			"""
		self.name_check(self.database.Constitutive)
		self.modify(linking)

	def modify(self, linking, head=None):
		if not head:
			head = self
		args = [Draw.Create(self._args[i]) for i in range(len(self._args))]
		delete = Draw.Create(0)
		single = Draw.Create(0)
		nval = Draw.Create(self.name)
		if self.type == 'Linear elastic':
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',					args[0], 'Dimensionality (select only one)'),
			('3',					args[1], 'Dimensionality (select only one)'),
			('6',					args[2], 'Dimensionality (select only one)'),
			('Stiffness*1.e-6:',	args[3], 0., 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Linear elastic generic':
			tmp = ['select matrix']
			for i in range(1):
				if self.links[i] != None: tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change stiffness matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:1]]
				eligible = ['1x1', '3x3', '6x6']
				if self._args[0] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', eligible)
				self._args[0] = 0
		elif self.type == 'Linear elastic generic axial torsion coupling':
			tmp = ['select matrix']*1
			for i in range(1):
				if self.links[i] != None: tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change stiffness matrix'),
			('Coef:',	args[1], -1.e6, 1.e6, 'Coupling coefficient')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:2]]
				eligible = ['6x1']
				if self._args[0] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', eligible)
				self._args[0] = 0
		elif self.type == 'Cubic elastic generic':
			tmp = ['select']*3
			eligible = ['1x1', '3x1']
			for i in range(3):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					eligible = set(eligible) & set([self.links[0].type])
			string = [
			('Name: ', nval, 0, 30),
			('S1:'+tmp[0],	args[0], 'Stiffness matrix 1'),
			('S2:'+tmp[1],	args[1], 'Stiffness matrix 2'),
			('S3:'+tmp[2], 	args[2], 'Stiffness matrix 3')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:3]]
				for i in range(3):
					if self._args[i] == 1 or self.links[i] == None:
						self.persistantPupMenu('Select S'+str(i+1))
						self.select(i, 'Matrix', eligible)
						eligible = set(eligible) & set([self.links[i].type])
				self._args[:3] = [0]*3
		elif self.type == 'Log elastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',	args[0], 0., 1.e6)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Linear elastic bi-stop generic':
			tmp = ['select matrix']
			if self.links[0] != None:
				tmp[0] = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],				args[0], 'Change stiffness matrix'),
			('Initial state',		args[1], 'Enable inactive/active selection'),
			('Active',				args[2], 'Select initial state as inactive or active')]
			tmp = ['select drive']*2
			for i in range(2):
				if self.links[i+1] != None:
					tmp[i] = self.links[i+1].name
			string += [
			(tmp[0],	args[3], 'Change activitating condition'),
			(tmp[1], 	args[4], 'Change deactivating condition')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:5]]
				eligible = ['1x1', '3x3', '6x6']
				if self._args[0] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', eligible)
				if args[3].val == 1 or self.links[2] == None:
					self.persistantPupMenu('Select activating condition:')
					self.select(1, 'Drive')
				if args[4].val == 1 or self.links[3] == None:
					self.persistantPupMenu('Select deactivating condition:')
					self.select(2, 'Drive')
				self._args[0] = 0
				self._args[3:5] = [0]*2
		elif self.type == 'Double linear elastic':
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('Stiffness 1:',	args[2], -1.e6, 1.e6),
			('Upper strain:',	args[3], -1.e6, 1.e6),
			('Lower Strain:',	args[4], -1.e6, 1.e6),
			('Stiffness 2:',	args[5], -1.e6, 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:2])
			self._args = [arg.val for arg in args[:6]]
		elif self.type == 'Isotropic hardening elastic':
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('6',		args[2], 'Dimensionality (select only one)'),
			('Stiffness:',			args[3], 0., 1.e6),
			('Reference strain:',	args[4], 0., 1.e6),
			('Linear stiffness',	args[5], 'Use linear stiffness parameter'),
			('Linear stiffness:',	args[6], 0., 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			self._args = [arg.val for arg in args[:7]]
		elif self.type == 'Scalar function elastic':
			tmp = 'select'
			if self.links[0] != None:
				tmp = self.links[0].name
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('6',		args[2], 'Dimensionality (select only one)'),
			('sf:'+tmp, args[3], 'Scalar function')]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			self._args = [arg.val for arg in args[:4]]
			if not delete.val:
				if self._args[3] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select scalar function')
					self.select(0, 'Function')
				self._args[3] = 0
		elif self.type == 'Scalar function elastic orthotropic':
			state = 2
			while not state == 1:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('6',		args[2], 'Dimensionality (select only one)')]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			if not delete.val:
				qty = args[0].val + 3*args[1].val + 6*args[2].val
				_qty = self._args[0] + 3*self._args[1] + 6*self._args[2]
				for i in range(_qty, qty):
					self.links.append(None)
					args.append(Draw.Create(1))
				for i in range(_qty-1, qty-1, -1):
					self.links[i].users.remove(self)
					del	self.links[i]
				string = []
				for i in range(qty):
					if not self.links[i]:
						string += [('select  f-'+str(i+1)+':', args[i+3], 'Select function')]
					else:
						string += [(self.links[i].name, args[i+3],
							'Modify this function')]
				if not Draw.PupBlock(self.type, string): return
				for i, function in enumerate(self.links):
					self.persistantPupMenu('Select  f-'+str(i+1)+':')
					self.select(i, 'Function')
				self._args = [arg.val for arg in args[:3]]+[0]*qty
		elif self.type == 'Linear viscous':
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('6',		args[2], 'Dimensionality (select only one)'),
			('Viscosity:',	args[3], 0., 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Linear viscous generic':
			tmp = ['select matrix']
			if self.links[0] != None:
				tmp[0] = self.links[0].name
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],	args[0], 'Change viscosity matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:1]]
				eligible = ['1x1', '3x3', '6x6']
				if self._args[0] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select viscosity matrix:')
					self.select(0, 'Matrix', eligible)
				self._args[0] = 0
		elif self.type == 'Linear viscoelastic':
			state = 2
			while not state in [1]:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('6',		args[2], 'Dimensionality (select only one)'),
			('Stiffness:',	args[3], 0., 1.e6),
			('Viscosity:',	args[4], 0., 1.e6),
			('Proportional',args[5], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[6], 0., 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:3])
			self._args = [arg.val for arg in args[:7]]
		elif self.type == 'Linear viscoelastic generic':
			tmp = ['select']*2
			eligible = ['1x1', '3x3', '6x6']
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					eligible = set(eligible) & set([self.links[0].type])
			string = [
			('Name: ', nval, 0, 30),
			('S:'+tmp[0],	args[0], 'Stiffness matrix'),
			('V:'+tmp[1], 	args[1], 'Viscosity matrix'),
			('Proportional',args[2], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3], 0., 1.e6)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:3]]
				if args[0].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', eligible)
					eligible = set(eligible) & set([self.links[0].type])
				if args[1].val == 1 or self.links[1] == None:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(1, 'Matrix', eligible)
			self._args = [0, 0] + [arg.val for arg in args[2:4]]
		elif self.type == 'Linear viscoelastic generic axial torsion coupling':
			tmp = ['select']*2
			eligible = ['6x1']
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					eligible = set(eligible) & set([self.links[0].type])
			string = [
			('Name: ', nval, 0, 30),
			('S:'+tmp[0],	args[0], 'Stiffness matrix'),
			('V:'+tmp[1], 	args[1], 'Viscosity matrix'),
			('Proportional',args[2], 'Override viscosity with proportion of stiffness'),
			('Proportion:',	args[3], 0., 1.e6),
			('Coupling:',	args[4], 0., 1.e6)]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:3]]
				if args[0].val == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix')
					self.select(0, 'Matrix', eligible)
				if args[1].val == 1 or self.links[1] == None:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(1, 'Matrix', eligible)
			self._args = [0, 0] + [arg.val for arg in args[2:5]]
		elif self.type == 'Cubic viscoelastic generic':
			tmp = ['select']*4
			eligible = set(['1x1', '3x1'])
			for i in range(4):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					eligible = set(eligible) & set([self.links[0].type])
			string = [
			('Name: ', nval, 0, 30),
			('S1:'+tmp[0],	args[0], 'Stiffness matrix 1'),
			('S2:'+tmp[1],	args[1], 'Stiffness matrix 2'),
			('S3:'+tmp[2], 	args[2], 'Stiffness matrix 3'),
			('V:'+tmp[3], 	args[3], 'Viscosity matrix')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:4]]
				for i in range(3):
					if self._args[i] == 1 or self.links[i] == None:
						self.persistantPupMenu('Select S'+str(i+1))
						self.select(i, 'Matrix', eligible)
						eligible = set(eligible) & set([self.links[i].matrixType])
				if self._args[3] == 1 or self.links[3] == None:
					self.persistantPupMenu('Select viscosity matrix')
					self.select(3, 'Matrix', eligible)
				self._args[:4] = [0]*4
		elif self.type == 'Double linear viscoelastic':
			state = 2
			while not state == 1:
				string = [
			('Name: ', nval, 0, 30),
			('1',		args[0], 'Dimensionality (select only one)'),
			('3',		args[1], 'Dimensionality (select only one)'),
			('Stiffness-1:',	args[2], -1.e6, 1.e6),
			('Upper strain:',	args[3], -1.e6, 1.e6),
			('Lower Strain:',	args[4], -1.e6, 1.e6),
			('Stiffness-2:',	args[5], -1.e6, 1.e6),
			('Viscosity:',		args[6], -1.e6, 1.e6),
			('Second damping',	args[7], 'Use viscosity-2 for outside lower/upper strain range'),
			('Viscosity-2:',	args[8], -1.e6, 1.e6)]
				if self in self.database.Constitutive:
					self.permit_deletion(linking, string, delete, single)
				if not Draw.PupBlock(self.type, string): return
				state = sum(arg.val for arg in args[:2])
			self._args = [arg.val for arg in args[:9]]
		elif self.type == 'Turbulent viscoelastic':
			string = [
			('Name: ', nval, 0, 30),
			('Stiffness:',		args[0], 0., 1.e6),
			('Parabolic visc:',	args[1], 0., 1.e6),
			('Use threshold',	args[2]),
			('Threshold:',		args[3], 0., 1.e6),
			('Use linear visc',	args[4]),
			('Linear visc:',	args[5], 0., 1.e6, 'Linear viscous coefficient')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:6]]
		elif self.type == 'Linear viscoelastic bi-stop generic':
			tmp = ['select matrix']*2
			eligible = ['1x1', '3x3', '6x6']
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
					eligible = set(eligible) & set([self.links[0].type])
			string = [
			('Name: ', nval, 0, 30),
			(tmp[0],				args[0], 'Change stiffness matrix'),
			(tmp[1],				args[1], 'Change viscosity matrix'),
			('Enable state',		args[2], 'Enables inactive/active selection'),
			('Active',				args[3], 'Initial state is inactive or active')]
			tmp = ['select drive']*2
			for i in range(2):
				if self.links[i+2] != None: tmp[i] = self.links[i+2].name
			string += [
			(tmp[0],	args[4], 'Change activitating condition'),
			(tmp[1], 	args[5], 'Change deactivating condition')]
			if self in self.database.Constitutive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			if not delete.val:
				self._args = [arg.val for arg in args[:6]]
				if self._args[0] == 1 or self.links[0] == None:
					self.persistantPupMenu('Select stiffness matrix:')
					self.select(0, 'Matrix', eligible)
					eligible = set(eligible) & set([self.links[0].type])
				if self._args[1] == 1 or self.links[1] == None:
					self.persistantPupMenu('Select viscosity matrix:')
					self.select(1, 'Matrix', eligible)
				if args[4].val == 1 or self.links[2] == None:
					self.persistantPupMenu('Select activating condition:')
					self.select(2, 'Drive')
				if args[5].val == 1 or self.links[3] == None:
					self.persistantPupMenu('Select deactivating condition:')
					self.select(3, 'Drive')
				self._args[0:2] = [0]*2
				self._args[4:6] = [0]*2
				"""
		elif self.type == 'GRAALL damper':
		elif self.type == 'shock absorber':
		elif self.type == 'symbolic elastic':
		elif self.type == 'symbolic viscous':
		elif self.type == 'symbolic viscoelastic':
		elif self.type == 'ann elastic':
		elif self.type == 'ann viscoelastic':
		elif self.type == 'nlsf viscoelastic':
		elif self.type == 'nlp viscoelastic':
		elif self.type == 'invariant angular ':
				"""
		else:
			Draw.PupMenu('Error: '+self.type+' is not defined in this program.')
			return
		return self.finalize(nval, delete, single, self.database.Constitutive)

	def string(self):
		if self.type in ['Linear elastic', 'Linear viscous']:
			if self._args[0]:
				string = self.type+', '
			elif self._args[1]:
				string = self.type+' isotropic, '
			elif self._args[2]:
				string = self.type+' isotropic, '
			string += str(self._args[3])
		elif self.type in ['Linear elastic generic', 'Linear viscous generic']:
			string = self.type+', '+self.links[0].string()
		elif self.type == 'Linear elastic generic axial torsion coupling':
			string = ('linear elastic generic axial torsional coupling, '+
				self.links[0].string()+', '+self._args[1])
		elif self.type == 'Cubic elastic generic':
			string = 'cubic elastic generic, '+str(self._args[0:3]).strip('[]')
		elif self.type == 'Log elastic':
			string = 'log elastic, '+str(self._args[0])
		elif self.type == 'Linear elastic bi-stop generic':
			string = 'linear elastic bistop, '+self.links[0].string()+',\n\t\t'
			if self._args[1]:
				if self._args[2]:
					string += 'initial state, active, \n\t\t'
				else:
					string += 'initial state, inactive, \n\t\t'
			string += (',\n\t\treference, '+str(self.database.Drive.index(self.links[1]))+
				',\n\t\treference, '+str(self.database.Drive.index(self.links[2])))
		elif self.type == 'Double linear elastic':
			string = 'double linear elastic, '+str(self._args[2:6]).strip('[]')
		elif self.type == 'Isotropic hardening elastic':
			string = 'isotropic hardening elastic, '+str(self._args[3:5]).strip('[]')
			if self._args[5]:
				string += '\n\t\tlinear stiffness, '+str(self._args[6])
		elif self.type == 'Scalar function elastic':
			string = 'scalar function elastic isotropic, "'+self.links[0].name+'"'
		elif self.type == 'Scalar function elastic orthotropic':
			string = 'scalar function elastic orthotropic,\n\t\t'
			if self._args[0]:
				N = 1
			elif self._args[1]:
				N = 3
			elif self._args[2]:
				N = 6
			for i in range(N):
				string += '"'+self.links[i].name+'", '
			string = string[:-2]
		elif self.type == 'Linear viscoelastic':
			if self._args[0]:
				string = ('linear viscoelastic isotropic, '+str(self._args[1])+', ')
			else:
				string = ('linear viscoelastic, '+str(self._args[1])+', ')
			if args[5]:
				string += 'proportional, '+str(self._args[6])
			else:
				string += str(self._args[4])
		elif self.type == 'Linear viscoelastic generic':
			string = 'linear viscoelastic generic,\n\t\t'+self.links[0].string()+',\n\t\t'
			if self._args[2]:
				string += 'proportional, '+str(self._args[3])
			else:
				string += self.links[1].string()
		elif self.type == 'Linear viscoelastic generic axial torsion coupling':
			string = 'linear viscoelastic generic axial torsional coupling,\n\t\t'
			string += self.links[0].string()+',\n\t\t'
			if self._args[2]:
				string += 'proportional, '+str(self._args[3])
			else:
				string += self.links[1].string()
			string += ', '+str(self._args[4])
		elif self.type == 'Cubic viscoelastic generic':
			string = 'cubic viscoelastic generic'
			for link in self.links:
				string += ', \n\t\t'+link.string()
		elif self.type == 'Double linear viscoelastic':
			string = 'double linear viscoelastic'
			for arg in self._args[2:7]:
				string += ', '+str(arg)
			if self._args[7]:
				string += ', \n\t\tsecond damping, '+str(self._args[8])
		elif self.type == 'Turbulent viscoelastic':
			string = 'turbulent viscoelastic'
			for arg in self._args[:2]:
				string += ', '+str(arg)
			if self._args[2]:
				string += ', '+str(self._args[3])
			if self._args[4]:
				string += ', '+str(self._args[5])
		elif self.type == 'Linear viscoelastic bi-stop generic':
			string = ('linear viscoelastic bistop,\n\t\t'+
				self.links[0].string()+',\n\t\t'+ self.links[1].string()+',\n\t\t')
			if self._args[2]:
				if self._args[3]:
					string += 'initial state, active, \n\t\t'
				else:
					string += 'initial state, inactive, \n\t\t'
			string += (',\n\t\treference, '+str(self.database.Drive.index(self.links[2]))+
				',\n\t\treference, '+str(self.database.Drive.index(self.links[3])))
		"""
		elif self.type == 'GRAALL damper':
		elif self.type == 'shock absorber':
		elif self.type == 'symbolic elastic':
		elif self.type == 'symbolic viscous':
		elif self.type == 'symbolic viscoelastic':
		elif self.type == 'ann elastic':
		elif self.type == 'ann viscoelastic':
		elif self.type == 'nlsf viscoelastic':
		elif self.type == 'nlp viscoelastic':
		elif self.type == 'invariant angular ':
		"""
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
			('Constant:',	args[0], -1.e6, 1.e6)]
			if self in self.database.Shape:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'linear':
			string = [
			('Name: ', nval, 0, 30),
			('Value at -1:',	args[0], -1.e6, 1.e6),
			('Value at +1:',	args[1], -1.e6, 1.e6)]
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
						seq[16*col+i] = ('Point_'+str(8*col+i)+': ', series[16*col+2*i], -1.e6, 1.e6)
						seq[16*col+8+i] = ('Value_'+str(8*col+i)+': ', series[16*col+2*i+1], -1.e6, 1.e6)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [N]
		elif self.type == 'parabolic':
			string = [
			('Name: ', nval, 0, 30),
			('Value at -1:',	args[0], -1.e6, 1.e6),
			('Value at 0:',		args[1], -1.e6, 1.e6),
			('Value at +1:',	args[2], -1.e6, 1.e6)]
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
		self._args = [0.]*10
		self._series = []
		self.links = []
		if self.type == '1x1':
			self._args[0] = 0
		elif self.type == '3x1':
			self._args[0] = self._args[4] = 0
			self._args[5] = 1.
		elif self.type == '6x1':
			self._args[0] = 0
		elif self.type == '3x3':
			self._args[:6] = [1]+[0]*5
			self._series = [0.]*9
		elif self.type == '6x6':
			self._args[:6] = [1]+[0]*5
			self._series = [0.]*36
		elif self.type == '6xN':
			self._args[:2] = [1, 0]
			self._args[2] = 2
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
			('Value:',	args[0], -1.e6, 1.e6)]
			if self in self.database.Matrix:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == '3x1':
			string = [
			('Name: ', nval, 0, 30),
			('Default',	args[0], 'Override with either default or null values'),
			('x1:',		args[1], -1.e6, 1.e6),
			('x2:',		args[2], -1.e6, 1.e6),
			('x3:',		args[3], -1.e6, 1.e6),
			('Scale',	args[4], 'Scale the vector'),
			('Factor:',	args[5], -1.e6, 1.e6)]
			if self in self.database.Matrix:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:6]]
		elif self.type == '6x1':
			string = [
			('Name: ', nval, 0, 30),
			('Default',	args[0], 'Override with either default or null values'),
			('x1:',	args[1], -1.e6, 1.e6),
			('x2:',	args[2], -1.e6, 1.e6),
			('x3:',	args[3], -1.e6, 1.e6),
			('x4:',	args[4], -1.e6, 1.e6),
			('x5:',	args[5], -1.e6, 1.e6),
			('x6:',	args[6], -1.e6, 1.e6)]
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
							string += [('', series[N*j+i], -1.e6, 1.e6)]
						string += ['']*3
					string[0] = 'matr'
				elif args[1].val == 1:
					string = []
					for i in range(N):
						string += ['']*2
						for j in range(i+1):
							ij = i+N*j
							string += [('x('+str(j+1)+','+str(i+1)+'):', series[ij], -1.e6, 1.e6)]
						string += ['']*(5-i)
					string[0] = 'sym'
				elif args[2].val == 1:
					string = [('skew')]
					for i,j in enumerate([7,2,3]):
						string += [('v('+str(i+1)+'):', series[j], -1.e6, 1.e6)]
				elif args[3].val == 1:
					string = [('diag')]
					for i in range(N):
						string += [('x('+str(i+1)+','+str(i+1)+'):', series[(N+1)*i], -1.e6, 1.e6)]
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
							string += [('', series[N*j+i], -1.e6, 1.e6)]
						string += ['']
					string[0] = 'matr'
				elif args[1].val == 1:
					string = []
					for i in range(N):
						string += ['']*1
						for j in range(i+1):
							ij = i+N*j
							string += [('x('+str(j+1)+','+str(i+1)+'):', series[ij], -1.e6, 1.e6)]
						string += ['']*(6-i)
					string[0] = 'sym'
				elif args[2].val == 1:
					string = [('diag')]
					for i in range(N):
						string += [('x('+str(i+1)+','+str(i+1)+'):', series[(N+1)*i], -1.e6, 1.e6)]
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
							string += [('', series[Ncols*j+i], -1.e6, 1.e6)]
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
			ret = str(self._args[0])
		elif self.type == '3x1':
			if self._args[0]:
				ret = 'default'
			else:
				ret = str(self._args[1:4]).strip('[]')
				if self._args[4]:
					ret += ', scale, '+str(self._args[5])
		elif self.type == '6x1':
			if self._args[0]:
				ret = 'default'
			else:
				ret = '\n\t\t'+str(self._args[1:7]).strip('[]')
		elif self.type == '3x3':
			if self._args[0] == 1:
				ret = ('\n\t\tmatr,\n'+
				'\t'*3+str(self._series[0:3]).strip('[]')+',\n'+
				'\t'*3+str(self._series[3:6]).strip('[]')+',\n'+
				'\t'*3+str(self._series[6:9]).strip('[]'))
			elif self._args[1] == 1:
				ret = ('\n\t\tsym,\n'+
				'\t'*3+str(self._series[0:3]).strip('[]')+',\n'+
				'\t'*4+str(self._series[4:6]).strip('[]')+',\n'+
				'\t'*5+str(self._series[8]))
			elif self._args[2] == 1:
				ret = 'skew, '+str([self._series[i] for i in [7, 2, 3]]).strip('[]')
			elif self._args[3] == 1:
				ret = 'diag, '+str([self._series[i] for i in [0, 4, 8]]).strip('[]')
			elif self._args[4] == 1:
				ret = 'eye'
			elif self._args[5] == 1:
				ret = 'null'
		elif self.type == '6x6':
			if self._args[0] == 1:
				ret = ('\n\t\tmatr,\n'+
				'\t'*3+str(self._series[0:6]).strip('[]')+',\n'+
				'\t'*3+str(self._series[6:12]).strip('[]')+',\n'+
				'\t'*3+str(self._series[12:18]).strip('[]')+',\n'+
				'\t'*3+str(self._series[18:24]).strip('[]')+',\n'+
				'\t'*3+str(self._series[24:30]).strip('[]')+',\n'+
				'\t'*3+str(self._series[30:36]).strip('[]'))
			elif self._args[1] == 1:
				ret = ('\n\t\tsym,\n'+
				'\t'*3+str(self._series[0:6]).strip('[]')+',\n'+
				'\t'*4+str(self._series[7:12]).strip('[]')+',\n'+
				'\t'*5+str(self._series[14:18]).strip('[]')+',\n'+
				'\t'*6+str(self._series[21:24]).strip('[]')+',\n'+
				'\t'*7+str(self._series[28:30]).strip('[]')+',\n'+
				'\t'*8+str(self._series[35]))
			elif self._args[2] == 1:
				ret = '\n\t\tdiag, '+str([self._series[i] for i in [0, 7, 14, 21, 28, 35]]).strip('[]')
			elif self._args[3] == 1:
				ret = 'eye'
			elif self._args[4] == 1:
				ret = 'null'
		elif self.type == '6xN':
			if self._args[0] == 1:
				N = self._args[2]
				ret = ('\n\t\tmatr')
				for i in range(N):
					ret += ',\n\t'*3+str(self._series[i*N:(i+1)*N]).strip('[]')
			elif self._args[1] == 1:
				ret = 'null'
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
			('Sigma 0:',	args[0], -1.e6, 1.e6),
			('Sigma 1:',	args[1], -1.e6, 1.e6),
			('Sigma 2:',	args[2], -1.e6, 1.e6),
			('Kappa:',		args[3], -1.e6, 1.e6),
			('sf:'+tmp, 	args[4], 'Change friction function'),
			('Plane hinge',	args[5], 'Simple plane hinge (else just simple shape function)'),
			('Radius:',		args[6], 0., 1.e6, 'Active only if plane hinge is selected')]
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
			('Sigma 2:',		args[1], 0., 1.e6),
			('Use vel ratio',	args[2], 'Use velocity ratio'),
			('Vel ratio:',		args[3], 0, 1., 'Velocity ratio'),
			('sf:'+tmp, 		args[4], 'Change friction function'),
			('Plane hinge',		args[5], 'Simple plane hinge (else just simple shape function)'),
			('Radius:',			args[6], 0., 1.e6, 'Active only if plane hinge is selected')]
			if self in self.database.Friction:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.frictType, string): return
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
			('Constant:',	args[0], -1.e6, 1.e6)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Exp':
			string = [
			('Name: ', nval, 0, 30),
			('Base',	args[0], 'Overrides default base (e)'),
			('Base:',	args[1], 0., 1.e6),
			('Coef',	args[2], 'Overrides default coefficient (1)'),
			('Coef:',	args[3], -1.e6, 1.e6),
			('Mult:',	args[4], -1.e6, 1.e6)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Log':
			string = [
			('Name: ', nval, 0, 30),
			('Base',	args[0], 'Overrides default base (e)'),
			('Base:',	args[1], 0., 1.e6),
			('Coef',	args[2], 'Overrides default coefficient (1)'),
			('Coef:',	args[3], 0., 1.e6),
			('Mult:',	args[4], -1.e6, 1.e6)]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Pow':
			string = [
			('Name: ', nval, 0, 30),
			('Exp:',	args[0], -1.e6, 1.e6, 'Exponent coefficient')]
			if self in self.database.Function:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:1]]
		elif self.type == 'Linear':
			string = [
			('Name: ', nval, 0, 30),
			('x1:',	args[0], -1.e6, 1.e6),
			('y1:',	args[1], -1.e6, 1.e6),
			('x2:',	args[2], -1.e6, 1.e6),
			('y2:',	args[3], -1.e6, 1.e6)]
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
						seq[16*col+i] = ('x'+str(8*col+i)+': ', series[16*col+2*i], -1.e6, 1.e6)
						seq[16*col+8+i] = ('y'+str(8*col+i)+': ', series[16*col+2*i+1], -1.e6, 1.e6)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [arg.val for arg in args[:2]]
		elif self.type == 'Chebychev':
			string = [
			('Name: ', nval, 0, 30),
			('Lower bound:',	args[0], -1.e6, 1.e6),
			('Upper bound:',	args[1], -1.e6, 1.e6),
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
					seq.append(('coef_'+str(i)+':', series[i], -1.e6, 1.e6))
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
#		elif self.type == 'File drive':
		elif self.type == 'String drive':
			self._args = 'enter text'
#		elif self.type == 'Dof drive':
#		elif self.type == 'Node drive':
#		elif self.type == 'Element drive':
		elif self.type == 'Drive drive':
			self._args = [1]*2
			self.links = [None]*2
		elif self.type == 'Array drive':
			self._args = [0]
#		elif self.type == 'Hints':
#		elif self.type == 'Template drive':
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
			('Constant:',	args[0], -1.e6, 1.e6)]
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
			('Constant:',	args[0], -1.e6, 1.e6),
			('Slope:',		args[1], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:2]]
		elif self.type == 'Parabolic drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -1.e6, 1.e6),
			('Linear:',		args[1], -1.e6, 1.e6),
			('Parabolic:',	args[2], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		elif self.type == 'Cubic drive':
			string = [
			('Name: ', nval, 0, 30),
			('Constant:',	args[0], -1.e6, 1.e6),
			('Linear:',		args[1], -1.e6, 1.e6),
			('Parabolic:',	args[2], -1.e6, 1.e6),
			('Cubic:',		args[3], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Step drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 1.e6),
			('Step:',			args[1], -1.e6, 1.e6),
			('Initial value:',	args[2], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:3]]
		elif self.type == 'Double step drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 1.e6),
			('Final time:',		args[1], 0., 1.e6),
			('Step value:',		args[2], -1.e6, 1.e6),
			('Initial value:',	args[3], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Ramp drive':
			string = [
			('Name: ', nval, 0, 30),
			('Slope:',			args[0], -1.e6, 1.e6),
			('Initial time:',	args[1], 0., 1.e6),
			('Forever',			args[2], 'Overrides final time'),
			('Final time:',		args[3], 0., 1.e6),
			('Initial value:',	args[4], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:5]]
		elif self.type == 'Double ramp drive':
			string = [
			('Name: ', nval, 0, 30),
			('Slope-a:',		args[0], -1.e6, 1.e6),
			('Initial time:',	args[1], 0., 1.e6),
			('Final time:',		args[2], 0., 1.e6),
			('Slope-d:',		args[3], -1.e6, 1.e6),
			('Initial time:',	args[4], 0., 1.e6),
			('Forever',			args[5], 'Overrides final time-d'),
			('Final time:',		args[6], 0., 1.e6),
			('Initial value:',	args[7], -1.e6, 1.e6)]
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
						seq[16*col+i] = ('Point_'+str(8*col+i)+': ', series[16*col+2*i], -1.e6, 1.e6)
						seq[16*col+8+i] = ('Value_'+str(8*col+i)+': ', series[16*col+2*i+1], -1.e6, 1.e6)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
				self._args = [N]
		elif self.type == 'Sine drive':
			state = 2
			while not state in [0,1]:
				string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 1.e6),
			('Omega:',			args[1], 0., 1.e6, 'radians per second'),
			('Amplitude:',		args[2], -1.e6, 1.e6),
			('N-cycles:',		args[3], -1e6, 1e6),
			('Initial value:',	args[4], -1.e6, 1.e6),
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
			('Initial time:',	args[0], 0., 1.e6),
			('Omega:',			args[1], 0., 1.e6, 'radians per second'),
			('Amplitude:',		args[2], -1.e6, 1.e6),
			('N-cycles:',		args[3], -1e6, 1e6),
			('Initial value:',	args[4], -1.e6, 1.e6),
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
			('Initial time:',	args[0], 0., 1.e6),
			('Amplitude:',		args[1], -1.e6, 1.e6),
			('Slope:',			args[2], -1.e6, 1.e6),
			('Initial value:',	args[3], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Fourier series drive':
			state = 2
			while not state in [0,1]:
				string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 1.e6),
			('Omega:',			args[1], -1.e6, 1.e6),
			('Number of terms:',args[2], 1, 50),
			('N-cycles:',		args[3], -1e6, 1e6),
			('One',				args[4], 'Overrides N-cycles (select only one override)'),
			('Forever',			args[5], 'Overrides N-cycles (select only one override)'),
			('Initial value:',	args[6], -1.e6, 1.e6)]
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
						seq[16*col+i] = ('a_'+str(8*col+i)+': ', series[16*col+2*i-1], -1.e6, 1.e6)
						seq[16*col+8+i] = ('b_'+str(8*col+i)+': ', series[16*col+2*i], -1.e6, 1.e6)
				seq[0] = ('Using only terms 0 to '+str(N-1))
				seq[8] = ('a_0: ', series[0], -1.e6, 1.e6)
				Draw.PupBlock(self.type, seq)
				self._series[:len(series)] = [item.val for item in series]
		elif self.type == 'Frequency sweep drive':
			tmp = ['select drive']*2
			for i in range(2):
				if self.links[i] != None:
					tmp[i] = self.links[i].name
			string = [
			('Name: ', nval, 0, 30),
			('Initial time: ',	args[0], 0., 1.e6),
			(tmp[0],		 	args[1], 'Change angular velocity drive'),
			(tmp[1], 			args[2], 'Change amplitude drive.'),
			('Initial value:',	args[3], -1.e6, 1.e6),
			('Forever',			args[4], 'Overrides final time'),
			('Final time:',		args[5], 0., 1.e6),
			('Final value:',	args[6], -1.e6, 1.e6)]
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
			('Amplitude:',		args[0], -1.e6, 1.e6),
			('Time constant:',	args[1], -1.e6, 1.e6),
			('Initial time:',	args[2], 0., 1.e6),
			('Initial value:',	args[3], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
		elif self.type == 'Random drive':
			string = [
			('Name: ', nval, 0, 30),
			('Amplitude:',		args[0], -1.e6, 1.e6),
			('Mean:',			args[1], -1.e6, 1.e6),
			('Initial time:',	args[2], 0., 1.e6),
			('Forever',			args[3], 'Overrides final time'),
			('Final time:',		args[4], 0., 1.e6),
			('Steps:',			args[5], 0, 1e6, 'Steps to hold value'),
			('Time seed',		args[6], 'Overrides Seed value with Time value'),
			('Seed:',			args[7], -1.e6, 1.e6)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:8]]
		elif self.type == 'Meter drive':
			string = [
			('Name: ', nval, 0, 30),
			('Initial time:',	args[0], 0., 1.e6),
			('Forever',			args[1], 'Overrides final time'),
			('Final time:',		args[2], 0., 1.e6),
			('Steps:',			args[3], 0, 1e6, 'Steps between spikes')]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = [arg.val for arg in args[:4]]
#		elif self.type == 'File drive': (need a Driver Class)
		elif self.type == 'String drive':
			string = [
			('Name: ', nval, 0, 30)]
			if self in self.database.Drive:
				self.permit_deletion(linking, string, delete, single)
			if not Draw.PupBlock(self.type, string): return
			self._args = Draw.PupStrInput('String: ', self._args, 100)
#		elif self.type == 'Dof drive':
#		elif self.type == 'Node drive':
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

	def string(self):
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
			'\t\treference, '+str(self.database.Drive.index(self._drives[0]))+',\n'+
			'\t\treference, '+str(self.database.Drive.index(self._drives[1]))+',\n'+
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
		elif self.type == 'File drive': pass
		elif self.type == 'String drive':
			return 'string, "'+self._args+'"'
		elif self.type == 'Dof drive': pass
		elif self.type == 'Node drive': pass
		elif self.type == 'Element drive': pass
		elif self.type == 'Drive drive':
			return ('drive,\n'+
			'\t\treference, '+str(self.database.Drive.index(self._drives[0]))+',\n'+
			'\t\treference, '+str(self.database.Drive.index(self._drives[1])))
		elif self.type == 'Array drive':
			string = ('array, '+str(self._args[0]))
			for i in range(self._args[0]):
				string += ',\n\t\treference, '+str(self.database.Drive.index(self._drives[i]))
			return string
		elif self.type == 'Hints': pass
		elif self.type == 'Template drive': pass


