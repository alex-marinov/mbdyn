#!BPY
"""
Name: 'MBDyn'
Blender: 246
Group: 'Animation'
Tooltip: 'Create, Run, and Import the results of an MBDyn multibody dynamic model'
"""

__author__ = "G. Douglas Baldwin, douglasbaldwin AT verizon.net"
__url__ = ["http://www.baldwintechnology.com"]
__version__ = "0.1.1"
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

import mbdyn
reload(mbdyn)
from mbdyn import MBDyn, Frame, Element, Constitutive, Shape
from mbdyn import Matrix, Friction, Function, Drive, Menu
import Blender, pickle
from Blender import Draw, Object, Scene, Text, Window, BGL, Mathutils, Ipo
from cStringIO import StringIO
from subprocess import call
from time import sleep, clock
from tempfile import TemporaryFile

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
		Window.FileSelector(mbdyn_callback, 'mbdyn Job Name.')
	elif evt == 2:
		for text in Text.Get():
			if text.name == '_Input': text.setName('_Input.old')
		text = Text.New('_Input')
		mbdyn.write(text)
		doPickle()
	elif evt == 3:
		run()
	elif evt == 4:
		display()
	Draw.Redraw(0)

def gui():             # the function to draw the screen
	Draw.PushButton("Parameters", 1, 10, 10, 80, 20, "mbdyn execution environment")
	Draw.PushButton("Input", 2, 95, 10, 55, 20, "Generate mbdyn input file")
	Draw.PushButton("Run", 3, 155, 10, 55, 20, "Run mbdyn for input file")
	Draw.PushButton("Display", 4, 215, 10, 55, 20, "Apply output position data to Blender objects")
	BGL.glRasterPos2i(10, 40)
	Draw.Text(mbdyn.filename)
 
Draw.Register(gui, builtin_event, button_event)  # registering the 3 callbacks

def start():
	menu = Menu(mbdyn)
	sel = Draw.PupTreeMenu(menu.string)
	if sel == -1: return
	menu.action(sel)
	doPickle()

def mbdyn_callback(filename):                # callback for the FileSelector
	if Blender.sys.exists(filename):
		if Draw.PupMenu("Save Over%t|"+filename) != 1: return
	mbdyn.filename = filename
	doPickle()

def run():
	if not mbdyn.filename:
		Draw.PupMenu('Try again after selecting a file location.')
		button_event(1)
		return
	try:
		Text.Get('_Input')
	except:
		Draw.PupMenu('No input file found.  Press "Input" button, then try again.')
		return
	tmp = open(mbdyn.filename, 'w')
	tmp.write('\n'.join(Text.Get('_Input').asLines()))
	tmp.close()
	sleep(1.)
#	Test for Linux or Windows
	if Blender.sys.sep == '/':
		error = runAsync(mbdyn.filename)
	else:
#		Windows does not support the execution Progress Bar feature, 
#		and may require editing C: in the following line
		error = call(['C:\Program Files\MBDyn\mbdyn.exe', '-s', '-f', mbdyn.filename])
	if error:
		Draw.PupMenu('mbdyn Error: check console for message')

def runAsync(filename):
	f1 = TemporaryFile()
	f2 = TemporaryFile()
	Window.DrawProgressBar(0., 'Running mbdyn...')
	command = 'mbdyn -s -f '+filename+' &'
	call(command, shell=True, stdout=f1)
	command = "tail -n 1 "+filename+".out | awk '{print $3}'"
	tNow = mbdyn._t0
	tN = mbdyn._tN
	if tN == 0.:
		tN = tNow + float(Blender.Get('endframe') - Blender.Get('staframe')
			)/float(mbdyn._fps)
	Dt = tN - tNow
	tHold = -1.0
	call("touch "+filename+".out", shell=True)
	sleep(0.1)
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
	error = 0
	if tHold < tN: error = 1
	if f1:
		f1.seek(0)
		print f1.read()
	f1.close()
	f2.close()
	return error

def display():
	filename = mbdyn.filename
	if Blender.sys.sep == '/':
		filename = filename.split('/')
		nameOnly = filename[-1].split('.')
		filename = '/'.join(filename[:-1])+'/'+nameOnly[0]+'.mov'
	else:
		filename = filename.split('\\')
		nameOnly = filename[-1].split('.')
		filename = '\\'.join(filename[:-1])+'\\'+nameOnly[0]+'.mov'
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
	key_Ipo = {}
	for key in mbdyn.rigid_dict.keys():
		key_Ipo[key] = key.getIpo()
	frame = Blender.Get('staframe')
	Blender.Set('curframe',frame)
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
		d2r = (3.14159/180.)
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
	Blender.Set('curframe',1)
	Blender.Redraw()

