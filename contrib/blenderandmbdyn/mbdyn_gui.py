#!BPY
"""
Name: 'MBDyn'
Blender: 246
Group: 'Animation'
Tooltip: 'Create, Run, and Import the results of an MBDyn multibody dynamic model'
"""

__author__ = "G. Douglas Baldwin, douglasbaldwin AT verizon.net"
__url__ = ["http://www.baldwintechnology.com"]
__version__ = "0.3.0"
__bpydoc__ = """\
Description:

This script and its companion script, mbdyn.py, provide a design environment for
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

import mbdyn
reload(mbdyn)
from mbdyn import MBDyn, Frame, Element, Constitutive, Shape
from mbdyn import Matrix, Friction, Function, Drive, Menu
import Blender, pickle, bpy
from Blender import Draw, Object, Scene, Text, Window, BGL, Mathutils, Ipo, Mesh
from cStringIO import StringIO
from subprocess import call, Popen
from time import sleep, clock
from tempfile import TemporaryFile
from socket import socket, AF_INET, SOCK_STREAM, gethostbyname, getservbyname
from struct import unpack

mbdyn = MBDyn()

# Hooks for pickle to reference Blender Objects
def persistent_id(obj):
	if hasattr(obj, 'getName'):
		return 'the value %s' % obj.getName()
	else:
		return None

def persistent_load(persid):
	if persid.startswith('the value '):
		value = persid.split()[2]
		return Object.Get(value)
	else:
		raise pickle.UnpicklingError, 'Invalid persistent id'

def persistent_load_append(persid):
	if persid.startswith('the value '):
		value = persid.split()[2]
		if value in mbdyn.append_objects:
			ob = mbdyn.append_objects[value]
		else:
			scn= bpy.data.scenes.active                       # get current scene
			lib = bpy.libraries.load(mbdyn.append_filename)          # open file.blend
			pseudoOb = lib.objects.append(value)            # get an object wrapper
			ob = scn.objects.link(pseudoOb)                                   # link to scene
			mbdyn.append_objects[value] = ob
		return ob
	else:
		raise pickle.UnpicklingError, 'Invalid persistent id'

def append_callback(filename):
	if not Blender.sys.exists(filename):
		Draw.PupMenu(filename+' does not exist') 
		return
	if Draw.PupMenu("Append mbdyn model from%t|"+filename):
		mbdyn.append_filename = filename
		mbdyn.append_objects = {}
		lib = bpy.libraries.load(filename)
		try:
			pickle_file = lib.texts.link('_Pickle')
		except:
			Draw.PupMenu(filename+' does not contain an mbdyn model')
			return
		datastream = '\n'.join(pickle_file.asLines())
		dst = StringIO(datastream)
		up = pickle.Unpickler(dst)
		up.persistent_load = persistent_load_append
		mbdyn_append = MBDyn()
		mbdyn_append = up.load()
		dst.close()
		bpy.data.texts.unlink(pickle_file)
		for clas in MBDyn.entity_classes+['Frame']:
			for item in eval('mbdyn_append.'+clas):
				item.database = mbdyn
				item.name_check(eval('mbdyn.'+clas))
				eval('mbdyn.'+clas+'.append(item)')
		del mbdyn.append_filename
		del mbdyn.append_objects
		del mbdyn_append

def doPickle():
	src = StringIO()
	p = pickle.Pickler(src)
	p.persistent_id = persistent_id
	for text in Text.Get():
		if '_Pickle' == text.name:  Text.unlink(text)
	freeze = Text.New("_Pickle")          # create a new Text object
	p.dump(mbdyn)
	freeze.write(src.getvalue())
	src.close()

# Restore mbdyn database from _Pickle
for text in Text.Get():
	if '_Pickle' == text.name:
		datastream = '\n'.join(text.asLines())
		dst = StringIO(datastream)
		up = pickle.Unpickler(dst)
		up.persistent_load = persistent_load
		mbdyn = up.load()
		dst.close()

def builtin_event(evt, val):    # the function to handle input events
	actions = {
		Draw.ESCKEY : 'Draw.Exit()', 
		Draw.RIGHTMOUSE : 'start()',
		Draw.DKEY : 'mbdyn.duplicate()'}
	if actions.has_key(evt) and val==1:  exec(actions[evt])

def button_event(evt):  # the function to handle Draw Button events
	if evt == 1:
		mbdyn.modify()
		if mbdyn._reinit:
			mbdyn.defaults()
		if mbdyn._change_filename:
			filename, ext= Blender.sys.splitext(Blender.Get('filename'))
			if filename:
				filename += '.mbd'
			else:
				filename = 'untitled.mbd'
			Window.FileSelector(mbdyn_callback, 'mbdyn Job Name.', filename)
		mbdyn._change_filename = 0
	elif evt == 2:
		if not mbdyn.filename:
			Draw.PupMenu('Try again after selecting a file location.')
			button_event(1)
			return
		for text in Text.Get():
			if text.name == '_Input': text.setName('_Input.old')
		text = Text.New('_Input')
		if mbdyn.write(text):
			prepare()
		doPickle()
	elif evt == 3:
		run()
	elif evt == 4:
		display()
	elif evt == 5:
		rigids()
	elif evt == 6:
		Draw.PupMenu('IMPORTANT: Be certain .mov has "default orientation: orientation matrix"')
		Window.FileSelector(import_mov_callback, '.mov File Name.')
	elif evt == 7:
		Window.FileSelector(append_callback, '.blend File Name.')
	Draw.Redraw(0)

def gui():             # the function to draw the screen
	Draw.PushButton("Parameters", 1, 10, 10, 80, 20, "mbdyn execution environment")
	Draw.PushButton("Input", 2, 95, 10, 55, 20, "Generate mbdyn input file")
	Draw.PushButton("Run", 3, 155, 10, 55, 20, "Run mbdyn for input file")
	Draw.PushButton("Display", 4, 215, 10, 55, 20, "Import results into IPos starting at current frame")
	Draw.PushButton("Rigids", 5, 275, 10, 55, 20, "Parent Rigids")
	Draw.PushButton("Import", 6, 335, 10, 55, 20, "Import .mov from externally run MBDyn model")
	Draw.PushButton("Append", 7, 395, 10, 55, 20, "Append MBDyn model from .blend file")
	BGL.glRasterPos2i(10, 40)
	Draw.Text(mbdyn.filename)
 
Draw.Register(gui, builtin_event, button_event)  # registering the 3 callbacks

def start():
	menu = Menu(mbdyn)
	sel = Draw.PupTreeMenu(menu.string)
	if sel == -1: return
	menu.action(sel)
	doPickle()

def import_mov_callback(filename):                # callback for the FileSelector
	if Blender.sys.exists(filename):
		if Draw.PupMenu("Import%t|"+filename) != 1: return
	obname, ext= Blender.sys.splitext(filename)
	obname = Blender.sys.basename(obname)[:6]+'_'
	import_mov(filename, obname)

def mbdyn_callback(filename):                # callback for the FileSelector
	if Blender.sys.exists(filename):
		if Draw.PupMenu("Save Over%t|"+filename) != 1: return
	mbdyn.filename = filename
	doPickle()

def run():
	ready = False
	for text in Text.Get():
		if text.name == '_Input': ready = True
	if not ready:
		Draw.PupMenu('Try again after creating Input file.')
		button_event(2)
		return
	try:
		f = Text.Get('_Input').asLines()
		for line in f:
			if line.find('stream') != -1:
				Draw.PupMenu('Error: This model must be run in the Blender Game Engine')
				return
	except:
		Draw.PupMenu('No input file found.  Press "Input" button, then try again.')
		return
	tmp = open(mbdyn.filename, 'w')
	tmp.write('\n'.join(Text.Get('_Input').asLines()))
	tmp.close()
	sleep(1.)
#	Test for Linux or Windows
	if Blender.sys.sep == '/':
		runAsync()
	else:
#		Windows does not support the execution Progress Bar feature, 
#		and may require editing of the following lines
#		or may not work at all if Blender in Windows cannot find 
#		the mbdyn executable and input files
		Draw.PupMenu('MBDyn will now run asyncronously in a seperate Window.')
		command = 'start', 'cmd', '/K', 'C:\Program Files\MBDyn\mbdyn.exe', '-s', '-f', mbdyn.filename
		call([command])

def runAsync():
	f1 = TemporaryFile()
	command = 'mbdyn -s -f '+mbdyn.filename+' &'
	print command
	filename, ext= Blender.sys.splitext(mbdyn.filename)
	process = Popen(command, shell=True, stdout=f1)
	Window.DrawProgressBar(0., 'Running mbdyn...')
	command = "tail -n 1 "+filename+".out | awk '{print $3}'"
	tNow = mbdyn._t0
	tN = mbdyn._tN
	if tN == 0.:
		tN = tNow + float(Blender.Get('endframe') - Blender.Get('staframe')
			)/float(mbdyn._fps)
	Dt = tN - tNow
	tHold = -1.0
	sleep(0.1)
	call("touch "+filename+".out", shell=True)
	sleep(0.1)
	f2 = TemporaryFile()
	while tHold < tN and tHold != tNow:
		tHold = tNow
		call(command, shell=True, stdout=f2)
		try:
			f2.seek(0)
			tNow = float(f2.read().splitlines()[-1])
			fraction = 1.-(tN-tNow)/Dt
			Window.DrawProgressBar(fraction, 'mbdyn: '+str(int(100.*fraction))+'%')
		except:
			pass
		sleep(0.1)
	f2.close()
	error = 0
	if tHold < tN:
		Draw.PupMenu('mbdyn Error: check console for message')
	if f1:
		f1.seek(0)
		print f1.read()
	f1.close()

def display():
	filename, ext= Blender.sys.splitext(mbdyn.filename)
	filename += '.mov'
	Window.DrawProgressBar(0., 'Loading: '+'0%')
	if mbdyn._IPOs:
		for node in mbdyn.Node:
			if node.getIpo():
				node.clearIpo()
			temp = Ipo.New('Object', node.name)
			node.setIpo(temp)
		for key in mbdyn.rigid_dict.keys():
			if key.getIpo():
				key.clearIpo()
			temp = Ipo.New('Object', key.name)
			key.setIpo(temp)
	else:
		for node in mbdyn.Node:
			if not node.getIpo():
				temp = Ipo.New('Object', node.name)
				node.setIpo(temp)
		for key in mbdyn.rigid_dict.keys():
			if not key.getIpo():
				key.clearIpo()
				temp = Ipo.New('Object', key.name)
				key.setIpo(temp)
	key_Ipo = {}
	for key in mbdyn.rigid_dict.keys():
		key_Ipo[key] = key.getIpo()
	save = Blender.Get('curframe')
	frame = save
	for node in mbdyn.Node:
		node.insertIpoKey(Object.IpoKeyTypes.LOCROT)
	for key in mbdyn.rigid_dict.keys():
		key.insertIpoKey(Object.IpoKeyTypes.LOCROT)
		key.clearIpo()
		mbdyn.rigid_dict[key].makeParent([key])

	movFile = open(filename)
	lines = movFile.readlines()
	try:
		marker = int(lines[0].split()[0])
	except:
		Draw.PupMenu('Error: Display file, '+filename+
		', does not match with itemDict.')
		return
	nLines = len(lines)
	iLines = 0
	portion = 0.
	timeMark = clock()
	scn = Scene.GetCurrent()
	radConv = 18./3.14159265358979323846
	for line in lines:
		timeCheck = clock()
		if timeCheck-timeMark >0.1:
			portion = float(iLines)/float(nLines)
			Window.DrawProgressBar(portion, 'Loading: '+str(int(100.*portion))+'%')
			timeMark = timeCheck
		iLines +=1
		fields = line.split()
		fields = [int(fields[0])] + [float(field) for field in fields[1:13]]
		if fields[0] == marker:
			scn.update(1)
			for key in mbdyn.rigid_dict.keys():
				ipo = key_Ipo[key]
				matrix = key.mat
				euler = key.getEuler('worldspace')
				ipo[Ipo.OB_LOCX].append((float(frame), key.mat[3][0]))
				ipo[Ipo.OB_LOCY].append((float(frame), key.mat[3][1]))
				ipo[Ipo.OB_LOCZ].append((float(frame), key.mat[3][2]))
				ipo[Ipo.OB_ROTX].append((float(frame), euler[0]*radConv))
				ipo[Ipo.OB_ROTY].append((float(frame), euler[1]*radConv))
				ipo[Ipo.OB_ROTZ].append((float(frame), euler[2]*radConv))
			frame += 1
			Blender.Set('curframe',frame)
		euler = Mathutils.Matrix(
		fields[4:7], fields[7:10], fields[10:13]).transpose().toEuler()
		d2r = (3.14159265358979323846/180.)
		try:
			mbdyn.Node[fields[0]].setLocation(fields[1:4])
			mbdyn.Node[fields[0]].setEuler([d2r*euler.x, d2r*euler.y, d2r*euler.z])	
			mbdyn.Node[fields[0]].insertIpoKey(Object.IpoKeyTypes.LOCROT)
		except:
			Draw.PupMenu('Error: Display file, '+filename+
			', does not match with mbdyn database.')
			return
	for key in mbdyn.rigid_dict.keys():
		key.clrParent()
		key.setIpo(key_Ipo[key])
	Window.DrawProgressBar(1., 'Loading: '+'100%')
#	Blender.Set('curframe',1)
	Blender.Set('curframe',save)
	Blender.Redraw()

def prepare():
	filename, ext= Blender.sys.splitext(mbdyn.filename)
	tmp = open(filename+'.mb2', 'w')
	tmp.write('\n'.join(Text.Get('_Input').asLines()))
	tmp.close()
	try:
		text = Text.Get('Always')
		text.clear()
	except:
		text = Text.New('Always')
	if mbdyn._posixRT:
		sudo = 'sudo '
	else:
		sudo = ''
	text.write(
"""from struct import pack, unpack
import pygame
try:
	pygame.event.pump() 
	""")
	joystick = False
	for i, driver in enumerate(mbdyn.Driver):
		if driver.users and driver.type == 'File' and driver._args[0]:
			for j, col in enumerate(driver.columns):
				if col._args[5]:
					joystick = True
					half_range = 0.5*(col._args[3] - col._args[2])
					text.write('GameLogic.drive_'+str(i)+'['+str(j)+'] = '+
					str(col._args[2] + half_range)+' + '+str(half_range/.7)+
					' * pygame.joystick.Joystick(0).get_axis('+str(col._args[6])+')\n\t'+
					'GameLogic.ob.'+col.name.replace(' ', '')+' = GameLogic.drive_'+str(i)+'['+str(j)+']\n\t')
			text.write("""data = ''
	columns = ''
	for column in GameLogic.drive_"""+str(i)+""":
		data += pack('d', column)
		columns += str(column)+'\\t'
	GameLogic.sock_"""+str(i)+""".send(data)
	GameLogic.file_"""+str(i)+""".write(columns[:-1]+'\\n')
	""")
	text.write("""try:
		data = unpack('d'*12*GameLogic.qty, GameLogic.sock_in.recv(GameLogic.qty*96))
		for i, name in zip(range(0, GameLogic.qty*12, 12), GameLogic.names):
			GameLogic.obs['OB'+name].position = data[i:i+3]
			GameLogic.obs['OB'+name].orientation = [data[i+3:i+6], data[i+6:i+9], data[i+9:i+12]]
	except:
		pass
except:
	try:
		GameLogic.shell_pid
	except:
		cont = GameLogic.getCurrentController()
		GameLogic.ob = cont.getOwner()
		pygame.init()
		""")
	if joystick:
		text.write("""pygame.joystick.init() 
		pygame.joystick.Joystick(0).init()
		""")
	text.write("""GameLogic.names = """+str([n.name for n in mbdyn.Node])+"""
		GameLogic.qty = len(GameLogic.names)
		scene = GameLogic.getCurrentScene()
		GameLogic.obs = scene.getObjectList()
		from subprocess import Popen
		from tempfile import TemporaryFile
		from socket import socket, AF_INET, SOCK_STREAM
#		GameLogic.setLogicTicRate(100.)
		sock_in = socket(AF_INET, SOCK_STREAM)
		sock_in.bind(('', """+str(mbdyn._port + len(mbdyn.Driver))+"""))				
		sock_in.listen(5)
		""")
	for i, driver in enumerate(mbdyn.Driver):
		if driver.users and driver.type == 'File' and driver._args[0]:
			text.write('sock_'+str(i)+""" = socket(AF_INET, SOCK_STREAM)
		sock_"""+str(i)+""".bind(('', """+str(mbdyn._port+i)+"""))				
		sock_"""+str(i)+""".listen(5)
		""")
	text.write("""GameLogic.f1 = TemporaryFile()
		command = ['"""+sudo+"""mbdyn -s -f """+filename+'.mb2'+""" &']
		print command[0]
		process = Popen(command, shell=True, stdout=GameLogic.f1)
		GameLogic.shell_pid = process.pid
		GameLogic.sock_in, GameLogic.address_in = sock_in.accept()
		""")
	for i, driver in enumerate(mbdyn.Driver):
		if driver.users and driver.type == 'File' and driver._args[0]:
			if driver._args[4]:
				filename += '.' + driver._args[4]
			else:
				filename += '.' + driver.name.replace(' ', '')
			text.write('GameLogic.sock_'+str(i)+', GameLogic.address_'+str(i)+' = sock_'+str(i)+""".accept()
		GameLogic.drive_"""+str(i)+' = '+str([d._args[1] for d in driver.columns])+"""
		GameLogic.qty_"""+str(i)+""" = len(GameLogic.drive_"""+str(i)+""")
		GameLogic.file_"""+str(i)+" = open('"+filename+"'"+""", 'wt')
		GameLogic.file_"""+str(i)+""".truncate()
		sock_"""+str(i)+""".close()
		""")
	text.write('sock_in.close()\n')

	for i, driver in enumerate(mbdyn.Driver):
		if driver.users and driver.type == 'File' and driver._args[0]:
			for j, col in enumerate(driver.columns):
				if not col._args[5]:
					try:
						text = Text.Get(col.name.replace(' ', ''))
						text.clear()
					except:
						text = Text.New(col.name.replace(' ', ''))
					text.write("""cont = GameLogic.getCurrentController()
for sensor in cont.getSensors():
	GameLogic.drive_"""+str(i)+'['+str(j)+'] += '+str(col._args[4])+' * (float(sensor.isPositive() + '+ """sensor.getInvert()) - 2.)
if GameLogic.drive_"""+str(i)+'['+str(j)+'] > '+str(col._args[3])+""":
	GameLogic.drive_"""+str(i)+'['+str(j)+'] = '+str(col._args[3])+"""
if GameLogic.drive_"""+str(i)+'['+str(j)+'] < '+str(col._args[2])+""":
	GameLogic.drive_"""+str(i)+'['+str(j)+'] = '+str(col._args[2])+"""
ob = cont.getOwner()
ob."""+text.name.replace('.', '_')+' = GameLogic.drive_'+str(i)+'['+str(j)+']\n')
	try:
		text = Text.Get('Quit')
		text.clear()
	except:
		text = Text.New('Quit')
	text.write('import pygame\n')
	if joystick:
		text.write('pygame.joystick.quit()\n')
	text.write("""from subprocess import Popen
if GameLogic.f1:
	GameLogic.f1.seek(0)
	print GameLogic.f1.read()
	GameLogic.f1.seek(0)
command = '"""+sudo+"""kill -9 '+str(GameLogic.shell_pid + 1)
print command
Popen(command, shell=True, stdout=GameLogic.f1)
if GameLogic.f1:
	GameLogic.f1.seek(0)
	print GameLogic.f1.read()
	GameLogic.f1.seek(0)
GameLogic.sock_in.close()
""")
	for i, driver in enumerate(mbdyn.Driver):
		if driver.users and driver.type == 'File' and driver._args[0]:
			text.write('GameLogic.sock_'+str(i)+'.close()\n')
	text.write("""GameLogic.f1.close()
co = GameLogic.getCurrentController()
for act in co.getActuators():
	GameLogic.addActiveActuator(act, True)""")

def rigids():
	for ob in Object.GetSelected():
		ob.sel = 0
	for element in mbdyn.Element:
		if element.type == 'Rigid':
			element.objects[0].sel = 1
			element.objects[1].makeParent([element.objects[0]])

def import_mov(filename, obname):
	sce= bpy.data.scenes.active
	me = Mesh.Primitives.Cube(2.0)
	Window.DrawProgressBar(0., 'Loading: '+'0%')
#	frame = 0
#	Blender.Set('curframe',frame)
	save = Blender.Get('curframe')
	frame = save -1
	movFile = open(filename)
	lines = movFile.readlines()
	try:
		marker = int(lines[0].split()[0])
	except:
		Draw.PupMenu('Error: Can not read .mov file')
		return
	nLines = len(lines)
	iLines = 0
	portion = 0.
	timeMark = clock()
	radConv = 18./3.14159265358979323846
	d2r = 3.14159265358979323846/180.
	objects = []
	obj_dict = {}
	for line in lines:
		timeCheck = clock()
		if timeCheck-timeMark >0.1:
			portion = float(iLines)/float(nLines)
			Window.DrawProgressBar(portion, 'Loading: '+str(int(100.*portion))+'%')
			timeMark = timeCheck
		iLines +=1
		fields = line.split()
		fields = [int(fields[0])] + [float(field) for field in fields[1:13]]
		if fields[0] == marker:
			sce.update(1)
			frame += 1
			Blender.Set('curframe',frame)
		if frame == save:
			ob = sce.objects.new(me)
			ob.name = obname+str(fields[0])
			ob.drawMode = Object.DrawModes['NAME']
			objects.append(ob)
			obj_dict[fields[0]] = len(obj_dict)
		euler = Mathutils.Matrix(
		fields[4:7], fields[7:10], fields[10:13]).transpose().toEuler()
		try:
			objects[obj_dict[fields[0]]].setLocation(fields[1:4])
			objects[obj_dict[fields[0]]].setEuler([d2r*euler.x, d2r*euler.y, d2r*euler.z])	
			objects[obj_dict[fields[0]]].insertIpoKey(Object.IpoKeyTypes.LOCROT)
		except:
			Draw.PupMenu('Error: Problem loading data.')
			return
	Window.DrawProgressBar(1., 'Loading: '+'100%')
	Blender.Set('curframe',save)
	Blender.Redraw()


