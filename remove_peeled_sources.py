# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:33:08 2015
Wide-field VLBI Reduction
@author: Jack Radcliffe
### ONLY USED TO LOAD FOR J123642+621331
Notes to usage 
-- MUST BE USED ON AN EMPTY AIPS DRIVE
"""
import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize, isfile
from os import listdir
from datetime import date
from collections import deque
#import Utilities
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()
import numpy as np

def parse_inp(filename):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(filename, "r")
	control = dict()

	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')

	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')

			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)

			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(', ')
			control[param] = valuelist

	return control
	


def flagger(name, klass, seq, indisk):
	''' create directory, write control file and run SERPent '''
	mkdir = 'serpent_' + name + '_' + klass + '_' + str(seq)
	os.system('mkdir ' + mkdir)
	path2folder = os.path.join(os.getcwd(), mkdir)
	path2file = os.path.join(path2folder, 'SERPent_input.py')
	print "Writing ", path2file, "for file: ", name, klass, seq
	# write ctrl file
	inputfile = open(path2file, 'w')
	print >> inputfile, 'AIPS_user_number = ', AIPS.userno
	print >> inputfile, 'Name = \'' + name + '\''
	print >> inputfile, 'Klass = \'' + klass + '\''
	print >> inputfile, 'Disk = ' + str(indisk)
	print >> inputfile, 'Seq = ' + str(seq)
	print >> inputfile, 'path2folder = \'' + path2folder + '/\''
	print >> inputfile, 'which_baselines = \'all\''
	print >> inputfile, 'baselines = [\'1-2\']'
	print >> inputfile, 'flagging_options = \'choose\''
	print >> inputfile, 'aggressiveness_first_run = ', aggr1
	print >> inputfile, 'max_subset_first_run = ', max1
	print >> inputfile, 'aggressiveness_second_run = ', aggr2
	print >> inputfile, 'max_subset_second_run = ', max2
	print >> inputfile, 'rho = ', rho
	print >> inputfile, 'kickout_sigma_level = ', kickoutsigma
	print >> inputfile, 'NCPU = ', int(ncpu)
	print >> inputfile, 'zero_level = \'yes\'' 
	print >> inputfile, 'phasecal = \'no\'' 
	print >> inputfile, 'write2log = \'no\''
	print >> inputfile, 'do_lovell_cross = \'no\''
	print >> inputfile, 'coadd_polarization_flags = \'no\''
	print >> inputfile, 'coadd_zero_flags = \'no\'' 
	print >> inputfile, 'coadd_lovell_flags = \'no\'' 
	print >> inputfile, 'dobline_specific = \'no\''
	print >> inputfile, 'flag_coinc_chans = ', 0 
	print >> inputfile, 'flag_coinc_times = ', 0 
	print >> inputfile, 'cparm_1 = ', 3.00 
	print >> inputfile, 'cparm_2 = ', 32 
	print >> inputfile, 'cparm_3 = ', 1.00 
	print >> inputfile, 'cparm_4 = ', 1.00 
	print >> inputfile, 'cparm_5 = ', 1.00 
	print >> inputfile, 'cparm_6 = ', 1.00 
	print >> inputfile, 'cparm_7 = ', 0 
	print >> inputfile, 'baselines_0 = [\'2-6\',\'2-7\']'
	print >> inputfile, 'baselines_1 = [\'6-8\',\'5-6\']'
	print >> inputfile, 'parameters_0 = [25,16,25,128,1.5,3.0]'  
	print >> inputfile, 'parameters_1 = [25,16,2,128,1.5,3.0]'
	inputfile.close()

	# run the flagger
	os.chdir(path2folder)
	os.system('cp ../SERPent_240315.py .') 
	os.system('ParselTongue.old SERPent_240315.py')
	os.chdir('../')

def checkin(control):
	''' convert the control hash to a list of global variables '''

	global userno, indisk, fitsdir, fitsfil, toload, plotdir, fittpdir, msgkill
	global aggr1, max1, aggr2, max2, rho, ncpu, kickoutsigma
	global doflag, doimag, imsize, filenam, niter, sourceimsize
	AIPS.userno = int(control.get('userno',[0])[0])
	filenam = int(control.get('filenam',[10])[0])
	indisk = int(control.get('indisk', [0])[0])
	aggr1 = int(control.get('aggr1',[25])[0])
	max1 = int(control.get('max1',[32])[0])
	aggr2 = int(control.get('aggr2',[25])[0])
	max2 = int(control.get('max2',[256])[0])
	rho = float(control.get('rho',[1.5])[0])
	ncpu = float(control.get('ncpu',[8])[0])
	kickoutsigma = float(control.get('kickoutsigma',[1.5])[0])
	doflag = int(control.get('doflag',[1])[0])
	doimag = int(control.get('doimag',[1])[0])
	#Imaging things
	imsize = int(control.get('imsize', [1024])[0])
	sourceimsize = int(control.get('sourceimsize', [1024])[0])
	if imsize:
		imsize = imsize,imsize
	if imsize:
		sourceimsize = sourceimsize,sourceimsize
	niter = int(control.get('niter', [1024])[0])
	# Options relating to loading from disk
	fitsdir = control.get('fitsdir', [])
	toload = 0

	# plotdir defaults to the cwd if not specified
	plotdir = control.get('plotdir', [False])[0]
	if (not plotdir):
		plotdir = os.getcwd()
	if (not os.path.isdir(plotdir)):
		print "Error:", plotdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(plotdir, os.W_OK) ):
		print "Error:", plotdir, "is not writable by you. Check your inputs."

	# fittpdir defaults to the cwd if not specified
	fittpdir = control.get('fittpdir', [False])[0]
	if (not fittpdir):
		fittpdir = os.getcwd()
	if (not os.path.isdir(fittpdir)):
		print "Error:", fittpdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(fittpdir, os.W_OK) ):
		print "Error:", fittpdir, "is not writable by you. Check your inputs."

	
	msgkill = int(control.get('msgkill', [-5])[0])

	if (not AIPS.userno) or (not indisk):
		print 'Error: You must set both your AIPS userno and indisk.'
		sys.exit()

def runimagr(indata, sources, docalib, gainuse, flagver, doband, bpver, bchan, echan, nchav, chinc, cellsiz, imsiz, niter, dotv, outdisk):
	imagr = AIPSTask('IMAGR')
	imagr.indata = indata
	imagr.sources[1:] = sources
	source = str(sources[0])
	if len(source)>12 :
		source = source[0:12]
	imagr.outname = source
	imagr.outdisk = outdisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1:] = cellsiz
	imagr.imsize[1:] = imsiz
	imagr.niter = niter
	imagr.dotv = dotv
	imagr.inp()
	imagr.go()

# Check for an inputs file
if len(sys.argv)==2 :
	print "Using inputs specified in", sys.argv[1]
	afile = sys.argv[1]
	if os.path.isfile(afile) :
		control = parse_inp(afile)
	else :
		print "Error:" + afile + "does not exist, quitting."
		sys.exit()
else :
	print "Error: no parameter file specified, quitting."
	print "Usage: parseltongue pipeline.py inputs.txt"
	sys.exit()

checkin(control)
pca = AIPSCat(indisk)
AIPSTask.msgkill = msgkill

def findmaxb(uvdata):
	maxbaseline = 0
	antab = uvdata.table('AN',1)
	for row in antab :
		for row2 in antab :
			xsep = row.stabxyz[0] - row2.stabxyz[0]
			ysep = row.stabxyz[1] - row2.stabxyz[1]
			zsep = row.stabxyz[2] - row2.stabxyz[2]
			hypxy =  math.sqrt((xsep * xsep) + (ysep * ysep))
			hypxyz = math.sqrt((zsep * zsep) + (hypxy * hypxy))
			if hypxyz > maxbaseline :
				maxbaseline = hypxyz
	cellsize = (1./2.)*(1.22 * (300000000 / uvdata.header.crval[2]) / maxbaseline) / 3.141592 * 180 * 3600 / 5
	print "maxbaseline = ", maxbaseline, "cellsize = ", cellsize
	return cellsize,cellsize

directory = '/nas/radcliff/eMERGE/JVLA_A_Calibrated_UV/UV/JVLA'
for i in range(7):
	
	for file in listdir('../UV/JVLA'+str(i+2)):
		if file.endswith('peeled.fits'):
			uvdataname = str(file)
			uvdataname = uvdataname
			AIPSuvdataname = uvdataname[:5]
	


	AIPSCat().zap()
	fitld = AIPSTask('FITLD')
	fitld.datain = (directory+str(i+2)+'/' + uvdataname) 
	fitld.ncount = 1
	fitld.doconcat = 0
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = AIPSuvdataname
	fitld.outclass = 'UV'
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.go()

	fitld = AIPSTask('FITLD')
	fitld.datain = (directory+str(i+2)+'/S1_model.fits') 
	fitld.ncount = 1
	fitld.doconcat = 0
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = 'S1'
	fitld.outclass = 'IM'
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.go()

	fitld = AIPSTask('FITLD')
	fitld.datain = (directory+str(i+2)+'/S2_model.fits') 
	fitld.ncount = 1
	fitld.doconcat = 0
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = 'S2'
	fitld.outclass = 'IM'
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.go()

	uvdata = AIPSUVData(AIPSuvdataname,'UV',indisk,1)
	image = AIPSImage('S1','IM',indisk,1)
	uvsub = AIPSTask('UVSUB')
	uvsub.indata = uvdata
	uvsub.outdata = uvdata
	uvsub.in2data = image
	uvsub.ncomp[1] = -1000000
	uvsub.factor = 1.0
	uvsub.opcode = ''
	uvsub.go()

	uvdata = AIPSUVData(AIPSuvdataname,'UV',indisk,1)
	image = AIPSImage('S2','IM',indisk,1)
	uvsub = AIPSTask('UVSUB')
	uvsub.indata = uvdata
	uvsub.in2data = image
	uvsub.ncomp[1] = -1000000
	uvsub.outdata = uvdata
	uvsub.factor = 1.0
	uvsub.opcode = ''
	uvsub.go()

	fittp = AIPSTask('FITTP')
	fittp.indata = uvdata
	fittp.indisk = indisk
	fittp.dataout = directory+str(i+2)+'/'+AIPSuvdataname+'_peeled_uvsub.fits'
	fittp.go()
