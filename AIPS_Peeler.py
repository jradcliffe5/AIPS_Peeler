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


def headless(inputfile):
    ''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
    INPUTFILE = open(inputfile, "r")
    control = {}
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
            value = value.replace(' ','').strip()
            valuelist = value.split(',')
            if len(valuelist) == 1:
                if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
                    control[param] = int(valuelist[0])
                else:
                    control[param] = str(valuelist[0])
            else:
                control[param] = ','.join(valuelist)
    return control

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

def self_calibration(uvdata, source_no, RA_shift, Dec_shift, solints, sc_type, pols, imsize, niter,refant):
	nchan = uvdata.header.naxis[2]
	sutable = uvdata.table('SU',1)
	sourceID = str(sutable[0].source)
	for i in range(len(solints)):
		print('Self calibration iteration %d' % i)
		for j in range(len(pols)):
			print('Iteration %d: Imaging source number %d' % (i,source_no))
			imagr = AIPSTask('IMAGR')
			pol_data = AIPSUVData(uvdata.name,pols[j],indisk,uvdata.seq)
			print pol_data.klass
			print pols[j]
			imagr.indata = pol_data
			imagr.sources[1] = sourceID
			#imagr.cellsize = AIPSList([findmaxb(uvdata)[0],findmaxb(uvdata)[1]])
			imagr.cellsize = AIPSList([0.2,0.2])
			imagr.imsize = AIPSList(imsize)
			imagr.docalib = 2
			imagr.gainuse = i+1
			imagr.stokes = pols[j]
			imagr.rashift[1] = RA_shift
			imagr.decshift[1] = Dec_shift
			imagr.niter = niter
			imagr.nchav = nchan
			imagr.outname = 'S%s_SC%s' % (source_no,str(i+1))
			imagr.outdisk = indisk
			imagr.go()

			image = AIPSImage('S%s_SC%s' % (source_no,str(i+1)),'%sCL001'%pols[j][0],indisk,1)

			print('Iteration %d: Making calibrations soln w/ solint=%.1f min and %s type' % (i,solints[i],sc_type[i]))
			calib = AIPSTask('CALIB')
			calib.indata = pol_data
			calib.calsour[1] = sourceID
			calib.docalib = 2
			calib.in2data = image
			calib.refant = refant
			calib.solint = solints[i]
			calib.aparm[1:] = 3, 0, 0, 0, 0, 0, 3, 0
			calib.ncomp[1] = -1000000
			calib.soltype = 'L1'
			calib.solmode = sc_type[i]
			calib.snver = i+1
			calib.go()

		print('Iteration %d: Applying SN table %d' % (i,i))
		clcal = AIPSTask('CLCAL')
		clcal.indata = pol_data
		clcal.calsour[1] = sourceID
		clcal.sources[1] = sourceID
		clcal.opcode = 'CALP'
		clcal.interpol = 'AMBG'
		clcal.snver = i+1
		clcal.invers = i+1
		clcal.gainver = i+1
		clcal.gainuse = i+2
		clcal.cutoff = 60
		clcal.dobtween = 0
		clcal.refant = refant
		clcal.go()

def peel_source(uvdata, pols, source_no, RA_shift, Dec_shift):
	for i in range(len(pols)):
		pol_data = AIPSUVData(uvdata.name,pols[i],uvdata.indisk,uvdata.seq)

		tasav = AIPSTask('TASAV')
		tasav.indata = uvdata
		tasav.outname = 'SC_SN'
		tasav.outdisk = indisk
		tasav.go()

		sntables = AIPSUVData('SC_SN','TASAV',indisk,1)

		sncor = AIPSTask('SNCOR')
		sncor.indata = sntables
		sncor.opcode = 'PNEG'
		for i in range(5):
			sncor.snver = i+1
			sncor.go()


		split = AIPSTask('SPLIT')
		split.indata = uvdata
		split.sources[1] = sourceID
		split.docalib = 2
		split.outdisk = indisk
		split.go()

		uvdata.zap()


### Inputs
inputs = headless('inputs.txt')
do_load = int(inputs['do_load'])
do_self_cal = int(inputs['do_self_cal'])
do_peel = int(inputs['do_peel'])
indisk = int(inputs['indisk'])
uvdataname = str(inputs['uvdata'])
AIPS.userno = int(inputs['AIPS_user'])
split_pol = int(inputs['split_pol'])
refant = int(inputs['refant'])
####


AIPSuvdataname = uvdataname[0:8]

if do_load == 1:
	print "------------------------------------------------------------"
	print "              Loading your uv data file                     "
	print "------------------------------------------------------------"
	AIPSCat().zap()
	fitld = AIPSTask('FITLD')
	fitld.datain = ('PWD:' + uvdataname)
	fitld.ncount = 1
	fitld.doconcat = 0
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = AIPSuvdataname[0:8]
	fitld.outclass = 'UV'
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.go()

	print "------------------------------------------------------------"
	print "                       Complete                             "
	print "------------------------------------------------------------"
	print "------------------------------------------------------------"
	print "                   Turning into multi source                "
	print "------------------------------------------------------------"
	uvdata = AIPSUVData(AIPSuvdataname,'UV',indisk,1)

	multi = AIPSTask('MULTI')
	multi.indata = uvdata
	multi.outname = uvdata.name
	multi.outclass = 'MULTI'
	multi.outdisk = indisk
	multi.go()

	uvdata.zap()
	AIPSUVData(uvdata.name,'MULTI',indisk,1).rename(uvdata.name,uvdata.klass,0)

	print "------------------------------------------------------------"
	print "               TB sort of the UV data file                  "
	print "------------------------------------------------------------"

	uvsrt = AIPSTask('UVSRT')
	uvsrt.sort = 'TB'
	uvsrt.indata = uvdata
	uvsrt.outdata = uvdata
	uvsrt.outdisk = indisk
	uvsrt.go()


	print "------------------------------------------------------------"
	print "                       Complete                             "
	print "------------------------------------------------------------"
	print "------------------------------------------------------------"
	print "               Indexing of the UV file                      "
	print "------------------------------------------------------------"

	uvdata.zap_table('CL',1)
	uvdata.zap_table('NX',1)
	indxr = AIPSTask('INDXR')
	indxr.indata = uvdata
	#indxr.cparm[3] = 0.25
	indxr.go()

	if split_pol == 1:
		pols = ['RR','LL']
		for i in ['RR','LL']:
			splat = AIPSTask('SPLAT')
			splat.indata = uvdata
			splat.stokes = i
			splat.outname = uvdata.name
			splat.outdisk = indisk
			splat.outclass = i
			splat.go()
	else:
		pols = ['I']
if split_pol == 1:
	pols = ['RR','LL']
else:
	pols = ['I']
uvdata = AIPSUVData(AIPSuvdataname,'UV',indisk,1)

nchan = uvdata.header.naxis[2]
sutable = uvdata.table('SU',1)
sourceID = str(sutable[0].source)
if do_self_cal == 1:
	self_calibration(uvdata=uvdata, source_no=1, RA_shift=-322.555, Dec_shift=-1018.34, solints=[1], sc_type=['P'], pols=pols, imsize=[1024,1024], niter=100,refant=refant)
if do_peel == 1:
	peel_source(uvdata,source_no,RA_shift,Dec_shift)
