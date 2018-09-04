# Example inputs for the e-MERLIN end-to-end pipeline.
userno = 1061		# An AIPS user number with an empty catalogue!
filenam = 10            # Number of characters in first part e.g EG078C.data.1 = 6
indisk = 1		# AIPS disk to use
doflag = 0             # Use automated flagger?
doimag = 1              # Use imaging defaults?
msgkill = 0
#############################################
# AUTO FLAGGER OPTIONS                      #
#############################################

# Flagging options - see the SERPent documentation for their meaning and default values
aggr1 = 15
max1 = 16
aggr2 = 25
max2 = 256
rho = 1.5
ncpu = 32 # Number of CPUs available (TRAAL = 16, FLAGGER=8, KRIA =24, etc)
kickoutsigma = 2

#Image options
imsize = 512
sourceimsize = 300
niter = 1 # clean iterations


