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

#### Inputs ####
for file in listdir('./'):
	if file.startswith('JVLA'): ## CHANGE WHEN DATA ARRIVES
		uvdataname = str(file)
		uvdataname = uvdataname
		AIPSuvdataname = uvdataname[:5]
for file in listdir('./'):
	if file.endswith('tasav'): ## CHANGE WHEN DATA ARRIVES
		tasavname = str(file)

AIPS.userno = 50 ## AIPS userno
indisk = 1  ## AIPS disk
sourceID = 'GOODSN3'
################

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

print "------------------------------------------------------------"
print "              Loading your uv data file                     "
print "------------------------------------------------------------"
print uvdataname
print AIPSuvdataname
#AIPSCat().zap()
print indisk

fitld = AIPSTask('FITLD')
fitld.datain = ('PWD:' + uvdataname)
fitld.ncount = 1
fitld.doconcat = 0
fitld.clint = 0
fitld.wtthresh = 0
fitld.outname = AIPSuvdataname
fitld.outclass = 'UV'
fitld.outdisk = indisk
fitld.digicor = -1
fitld.go()

print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"

uvdata = AIPSUVData(AIPSuvdataname,'UV',indisk,1)

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
indxr = AIPSTask('INDXR')
indxr.indata = uvdata
#indxr.cparm[3] = 0.25
indxr.go()

nxtable = uvdata.table('NX',1)
sourceID = str(nxtable[3].source_id)
print sourceID


split = AIPSTask('SPLIT')
split.indata = uvdata
split.sources[1] = sourceID
split.outdisk = indisk
split.outclass = 'SPLIT'
sourceID = source[i]
split.go()

# zap the original uv data to save space
uvdata.zap()
uvdata = AIPSUVData(sourceID,'SPLIT',indisk,1)
print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"
print "------------------------------------------------------------"
print "               UVFIX to first position                      "
print "------------------------------------------------------------"

uvfix = AIPSTask('UVFIX')
uvfix.indata = uvdata
uvfix.shift[1:] = -495.73, 394.87 ## in arcsec (need cos(dec) too)
uvfix.outclass = 'UVFIX'
uvfix.outdisk = indisk
uvfix.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'UVFIX',indisk,1)

multi = AIPSTask('MULTI')
multi.indata = uvdata
multi.outname = uvdata.name
multi.outclass = 'MULTI'
multi.outdisk = indisk
multi.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'MULTI',indisk,1)

uvdata.zap_table('CL',1)
indxr = AIPSTask('INDXR')
indxr.indata = uvdata
#indxr.cparm[3] = 0.25
indxr.go()

print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"
print "------------------------------------------------------------"
print "            Self calibration of source 1                    "
print "------------------------------------------------------------"

for i in range(5):

	imagr = AIPSTask('IMAGR')
	imagr.indata = uvdata
	imagr.sources[1] = sourceID
	imagr.cellsize[1:] = findmaxb(uvdata)
	imagr.imsize[1:] = sourceimsize
	imagr.docalib = 2
	imagr.gainuse = i+1
	imagr.niter = 1000
	imagr.nchav = 64
	imagr.outname = 'SC'+str(i+1)
	imagr.outdisk = indisk
	imagr.go()

	image = AIPSImage('SC'+str(i+1),'ICL001',indisk,1)

	calib = AIPSTask('CALIB')
	calib.indata = uvdata
	calib.calsour[1] = sourceID
	calib.docalib = 2
	calib.in2data = image
	calib.refant = 15
	if i+1 == 1:
		calib.solint = 5
	if i+1 == 2:
		calib.solint = 2
	if i+1 == 3:
		calib.solint = 2
	if i+1 == 4:
		calib.solint = 1
	if i+1 == 5:
		calib.solint = 1
	calib.aparm[1:] = 3, 0, 0, 0, 0, 0, 3, 0
	calib.ncomp[1] = -1000000
	calib.soltype = 'L1'
	calib.solmode = 'P'
	calib.snver = i+1
	calib.go()

	clcal = AIPSTask('CLCAL')
	clcal.indata = uvdata
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
	clcal.refant = 15
	clcal.go()

print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"
print "------------------------------------------------------------"
print "               Image and uvsub each IF                      "
print "------------------------------------------------------------"

#uvdata = AIPSUVData('1IF_TEST','SPLIT',2,1)
nchan1 = 64 #no of channels to average for source 1 ra 12:35:38.5 dec +62:19:32.87

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

uvdata = AIPSUVData(sourceID,'SPLIT',indisk,1)

imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.docalib = -1
imagr.gainuse = 0
imagr.outname = 'S1_B_UVSUB'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 1000
imagr.go()

for j in range(3):
	for i in range(16):

		print 'IF = '+str(i+1)
		imagr = AIPSTask('IMAGR')
		imagr.indata = uvdata
		imagr.cellsize[1:] = findmaxb(uvdata)
		imagr.imsize[1:] = sourceimsize
		imagr.sources[1] = sourceID
		imagr.outname = 'S1_IF'+str(i+1)
		imagr.outseq = j+1
		imagr.bif = i+1
		imagr.docalib = -1
		imagr.eif = i+1
		imagr.nchav = nchan1
		imagr.bchan = 1
		imagr.echan = 0
		imagr.chinc = 0
		imagr.outdisk = indisk
		imagr.niter = 600
		try:
			imagr.go()
		except RuntimeError:
			continue

		spectral_cube = AIPSImage('S1_IF'+str(i+1),'ICL001',indisk,j+1)


		print 'channel= '+str('all')
		print 'CC file= 1'#str(int(float(j)/float(nchan1))+1)
		print 'IF = ' +str(i+1)
		uvsub = AIPSTask('UVSUB')
		uvsub.indata = uvdata
		uvsub.in2data = spectral_cube
		uvsub.outdata = uvdata
		uvsub.outdisk = indisk
		uvsub.bif = i+1
		uvsub.eif = i+1
		#uvsub.channel = j+1
		uvsub.inver = 0 #int(float(j)/float(nchan1))+1
		uvsub.ncomp[1] = -1000000
		try:
			uvsub.go()
		except RuntimeError:
			continue


imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.docalib = -1
imagr.outname = 'S1_P_UVCL'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 1000
imagr.go()

multi = AIPSTask('MULTI')
multi.indata = uvdata
multi.outname = uvdata.name
multi.outclass = 'SRC1'
multi.outdisk = indisk
multi.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'SRC1',indisk,1)

uvdata.zap_table('CL',1)
indxr = AIPSTask('INDXR')
indxr.indata = uvdata
#indxr.cparm[3] = 0.25
indxr.go()

tacop = AIPSTask('TACOP')
tacop.indata = sntables
tacop.ncount = 1
tacop.inext = 'SN'
tacop.outdata = uvdata
for i in range(5):
	j = [5,4,3,2,1]
	tacop.invers = j[i]
	tacop.outvers = i+1
	tacop.go()

clcal = AIPSTask('CLCAL')
clcal.indata = uvdata
clcal.calsour[1] = sourceID
clcal.sources[1] = sourceID
clcal.opcode = 'CALP'
clcal.interpol = 'AMBG'
clcal.snver = 1
clcal.invers = 5
clcal.gainver = 1
clcal.gainuse = 2
clcal.cutoff = 60
clcal.dobtween = 0
clcal.refant = 15
clcal.go()

split = AIPSTask('SPLIT')
split.indata = uvdata
split.sources[1] = sourceID
split.docalib = 2
split.outdisk = indisk
split.go()

uvdata = AIPSUVData(sourceID,'SPLIT',indisk,1)
image = AIPSImage('SC5','ICL001',indisk,1)
uvsub = AIPSTask('UVSUB')
uvsub.indata = uvdata
uvsub.in2data = image
uvsub.outdata = uvdata
uvsub.ncomp[1] = -1000000
uvsub.factor = -1.0
uvsub.opcode = ''
uvsub.go()

imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.docalib = -1
imagr.outname = 'S1_MOD_B'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 1000
imagr.go()



print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"
print "------------------------------------------------------------"
print "            Recompute u,v,w for source 2 to remove          "
print "------------------------------------------------------------"

uvfix = AIPSTask('UVFIX')
uvfix.indata = uvdata
uvfix.shift[1:] = -322.555, -1018.34
uvfix.outclass = 'UVFIX'
uvfix.outdisk = indisk
uvfix.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'UVFIX',indisk,1)

multi = AIPSTask('MULTI')
multi.indata = uvdata
multi.outname = uvdata.name
multi.outclass = 'MULTI'
multi.outdisk = indisk
multi.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'MULTI',indisk,1)

print "------------------------------------------------------------"
print "                       Complete                             "
print "------------------------------------------------------------"
print "------------------------------------------------------------"
print "            Image spectral cube for source 2                "
print "------------------------------------------------------------"
nchan2 =64

for i in range(5):

	imagr = AIPSTask('IMAGR')
	imagr.indata = uvdata
	imagr.sources[1] = sourceID
	imagr.cellsize[1:] = findmaxb(uvdata)
	imagr.imsize[1:] = sourceimsize
	imagr.docalib = 2
	imagr.gainuse = i+1
	imagr.niter = 1000
	imagr.nchav = nchan2
	imagr.outname = '2SC'+str(i+1)
	imagr.outdisk = indisk
	imagr.go()

	image = AIPSImage('2SC'+str(i+1),'ICL001',indisk,1)

	calib = AIPSTask('CALIB')
	calib.indata = uvdata
	calib.calsour[1] = sourceID
	calib.docalib = 2
	calib.in2data = image
	calib.refant = 15
	if i+1 == 1:
		calib.solint = 5
	if i+1 == 2:
		calib.solint = 2
	if i+1 == 3:
		calib.solint = 2
	if i+1 == 4:
		calib.solint = 1
	if i+1 == 5:
		calib.solint = 1
	calib.aparm[1:] = 3, 0, 0, 0, 0, 0, 3, 0
	calib.ncomp[1] = -1000000
	calib.soltype = 'L1'
	calib.solmode = 'P'
	calib.snver = i+1
	calib.go()

	clcal = AIPSTask('CLCAL')
	clcal.indata = uvdata
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
	clcal.refant = 15
	clcal.go()

tasav = AIPSTask('TASAV')
tasav.indata = uvdata
tasav.outname = '2SC_SN'
tasav.outdisk = indisk
tasav.go()

sntables = AIPSUVData('2SC_SN','TASAV',indisk,1)

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

uvdata = AIPSUVData(sourceID,'SPLIT',indisk,1)

imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.outname = 'S2_B_UVSUB'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 500
imagr.go()


for j in range(5):
	for i in range(16):
		print 'IF = ' +str(i+1)
		imagr = AIPSTask('IMAGR')
		imagr.indata = uvdata
		imagr.cellsize[1:] = findmaxb(uvdata)
		imagr.imsize[1:] = sourceimsize
		imagr.outname = 'S2_IF'+str(i+1)
		imagr.bif = i+1
		imagr.eif = i+1
		imagr.outseq = j+1
		imagr.nchav = nchan2
		imagr.bchan = 1
		imagr.echan = 0
		imagr.chinc = 0
		imagr.outdisk = indisk
		imagr.niter = 2000
		try:
			imagr.go()
		except RuntimeError:
			continue

		spectral_cube = AIPSImage('S2_IF'+str(i+1),'ICL001',indisk,j+1)


		print 'channel= '+str('all')
		print 'CC file= 1'#str(int(float(j)/float(nchan1))+1)
		print 'IF = ' +str(i+1)
		uvsub = AIPSTask('UVSUB')
		uvsub.indata = uvdata
		uvsub.in2data = spectral_cube
		uvsub.outdata = uvdata
		uvsub.outdisk = indisk
		uvsub.bif = i+1
		uvsub.eif = i+1
		#uvsub.channel = j+1
		uvsub.inver = 0 #int(float(j)/float(nchan1))+1
		uvsub.ncomp[1] = -1000000
		try:
			uvsub.go()
		except RuntimeError:
			continue

imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.docalib = -1
imagr.outname = 'S2_P_UVCL'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 1000
imagr.go()

multi = AIPSTask('MULTI')
multi.indata = uvdata
multi.outname = uvdata.name
multi.outclass = 'SRC2'
multi.outdisk = indisk
multi.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'SRC2',indisk,1)

uvdata.zap_table('CL',1)
indxr = AIPSTask('INDXR')
indxr.indata = uvdata
#indxr.cparm[3] = 0.25
indxr.go()

tacop = AIPSTask('TACOP')
tacop.indata = sntables
tacop.ncount = 1
tacop.inext = 'SN'
tacop.outdata = uvdata
for i in range(5):
	j = [5,4,3,2,1]
	tacop.invers = j[i]
	tacop.outvers = i+1
	tacop.go()

clcal = AIPSTask('CLCAL')
clcal.indata = uvdata
clcal.calsour[1] = sourceID
clcal.sources[1] = sourceID
clcal.opcode = 'CALP'
clcal.interpol = 'AMBG'
clcal.snver = 1
clcal.invers = 5
clcal.gainver = 1
clcal.gainuse = 2
clcal.cutoff = 60
clcal.dobtween = 0
clcal.refant = 15
clcal.go()


split = AIPSTask('SPLIT')
split.indata = uvdata
split.sources[1] = sourceID
split.docalib = 2
split.outdisk = indisk
split.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'SPLIT',indisk,1)
image = AIPSImage('2SC5','ICL001',indisk,1)
uvsub = AIPSTask('UVSUB')
uvsub.indata = uvdata
uvsub.in2data = image
uvsub.ncomp[1] = -1000000
uvsub.outdata = uvdata
uvsub.factor = -1.0
uvsub.opcode = ''
uvsub.go()

imagr = AIPSTask('IMAGR')
imagr.indata = uvdata
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = sourceimsize
imagr.docalib = -1
imagr.outname = 'S2_MOD_B'
imagr.bif = 0
imagr.eif = 0
imagr.nchav = 64
imagr.outdisk = indisk
imagr.niter = 1000
imagr.go()

uvfix = AIPSTask('UVFIX')
uvfix.indata = uvdata
uvfix.shift[1:] = 819.47, 623.47
uvfix.outclass = 'UVFIX'
uvfix.outdisk = indisk
uvfix.outseq = 3
uvfix.go()

uvdata.zap()
uvdata = AIPSUVData(sourceID,'UVFIX',indisk,1)

fittp = AIPSTask('FITTP')
fittp.indata = uvdata
fittp.dataout = 'PWD:'+AIPSuvdataname+'_peeled.fits'
fittp.go()

#os.system('rsync -ar --progress '+AIPSuvdataname+'_split_uvsub_pass2.fits radcliff@marconi:/nas/radcliff/eMERGE/JVLA_A_Calibrated_UV/'+AIPSuvdataname)
