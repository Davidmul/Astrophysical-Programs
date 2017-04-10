
from __future__ import division
import pyfits
import numpy as np
#import pylab as pyl
#import matplotlib.pyplot as plt
import optparse
#import os
import sys

#############################################PARSER OPTIONS#########################################################

parser = optparse.OptionParser()

required = optparse.OptionGroup(parser, "Required Attributes")

required.add_option('--P', help='Input P FITS cube', action='store', type='string', dest='inputp')
required.add_option('--Q', help='Input Q FITS cube', action='store', type='string', dest='inputq')
required.add_option('--U', help='Input U FITS cube', action='store', type='string', dest='inputu')
required.add_option('--freq', help='Input frequency list', action='store', type='string', dest='frequency')
required.add_option('--rms', help='The rms noise in a single Q- or U-channel map', action='store', type='float', dest='rms')
required.add_option('--fore', help='For subtract of the FD of Galactic Foreground. Default is 0 radm^-2', action='store', type='float', dest='fore',default = 0.0)
required.add_option('--cut', help='Cutoff in Jy. Normally five times the rms noise in the Faraday Q- or U-cubes', action='store', type='float', dest='cutoff')
parser.add_option("-o", help='Overwrite existing fits files', action="store_true", dest="overwrite",default=False)
required.add_option('--outfile', help='Output fits file',default = 'MAXFARADAYDEPTH.FITS', action='store', type='string', dest='outputfits')

parser.add_option_group(required)

options, arguments = parser.parse_args()

#######Define parabolic function for fitting##################

def parafunc(f, x):
	xv = 1/2. * (f[x-1] - f[x+1]) / (f[x-1] - 2 * f[x] + f[x+1]) + x
	yv = f[x] - 1/4. * (f[x-1] - f[x+1]) * (xv - x)
	return (xv, yv)
#This method comes from Smith and Serra (1987) ---peak detection from sampled spectra
#http://quod.lib.umich.edu/cgi/p/pod/dod-idx/parshl-an-analysissynthesis-program-for-non-harmonic-sounds.pdf?c=icmc;idno=bbp2372.1987.042
# This could simply be changed to find the maximum, this would be fine for LOFAR data.
###########################LOAD FITS FILE FOR RM cube####################


try:
    hdulist = pyfits.open(options.inputp)
    headerdata = hdulist[0].header
    scidata = hdulist[0].data
    x_size = scidata.shape[-1]
    y_size = scidata.shape[-2]
    z_size = scidata.shape[-3]
    print 'zsize is', z_size
except IOError:
    print('ERROR.Polarsed Intensity RM cube is not found! Please check your input.')
    sys.exit(0)
try:
    freq = np.loadtxt(options.frequency)
except IOError:
    print('ERROR.Frequency list is not found! Please check your input.')
    sys.exit(0)

c = 299792458.0 # Speed of light
lam2 = (c/freq)**2.0
lam02 = np.mean(lam2)

Faradaystart = headerdata['CRVAL3']
FDincrement = headerdata['CDELT3']

#######Define Output arrays##############

outputrmarray=np.nan*np.zeros((y_size,x_size), dtype=np.float32)#output array setup
outputfluxarray = np.nan*np.zeros((y_size,x_size), dtype=np.float32)
FDindexplane = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

###########Fitting Faraday Spectra#########

for x in range(x_size):
	for y in range(y_size):
            axis = scidata[:,y,x]
            #print np.argmax(axis)
            if np.argmax(axis)==0: #if the max peak is the first point the spectra, was having issues with large sidelobes at the last/first pixel 
                FDindex = Faradaystart
                PI = axis[0]
                outputfluxarray[y,x] = PI
                FDindexplane[y,x] = 0
                        
            elif np.argmax(axis)==z_size-1:#if the max peak is the last point the spectra
                FDindex = (z_size-1)*FDincrement+Faradaystart
                PI = axis[z_size-1]
                outputfluxarray[y,x] = PI
                FDindexplane[y,x] = z_size-1
                
            else:
                FDindex, PI = parafunc(axis, np.argmax(axis))
		outputfluxarray[y,x] = PI
			    
		if PI > options.cutoff:
                    print 'PI is',PI
                    print 'cutoff is:',options.cutoff
		    outputrmarray[y,x] = FDincrement*FDindex + Faradaystart - options.fore
		    FDindexplane[y,x] = FDindex #needed to find Q and U
		else:
                    print 'naning'
		    outputrmarray[y,x] = np.nan
		    FDindexplane[y,x] = FDindex #needed so there are no nans when computing Q and U
		    #print FDindex



########### Write out the RM and PI images   ############

headerdata['NAXIS'] = 2
del headerdata['NAXIS3']
del headerdata['CDELT3']
del headerdata['CRPIX3']
del headerdata['CRVAL3']
del headerdata['CTYPE3']
del headerdata['CUNIT3']
headerdata['BUNIT'] = 'rad/m2' 
hdulist[0].data = outputrmarray
hdulist.writeto(options.outputfits, clobber=options.overwrite)

headerdata['BUNIT'] = 'Jy/beam' 
hdulist[0].data = outputfluxarray
hdulist.writeto('outputflux.fits', clobber=options.overwrite)


#######Input Q & U cubes###########

try:
    hdulistq = pyfits.open(options.inputq)
    headerdataq = hdulistq[0].header
    qscidata = hdulistq[0].data
except IOError:
    print('ERROR.Stokes Q cube is not found! Please check your input.')
    sys.exit(0)

try:
    hdulistu = pyfits.open(options.inputu)
    uscidata = hdulistu[0].data
except IOError:
    print('ERROR.Stokes U cube is not found! Please check your input.')
    sys.exit(0)

##########Create Arrays for Output max Q,U,PA and error maps#######

output_u=np.nan*np.zeros((y_size,x_size), dtype=np.float32)#output array setup
output_q = np.nan*np.zeros((y_size,x_size), dtype=np.float32)
output_pa = np.nan*np.zeros((y_size,x_size), dtype=np.float32)
sigma_phi = np.nan*np.zeros((y_size,x_size), dtype=np.float32)
sigma_pa = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

#Finding Stokes Q and U from standard interpolation

for x in range(x_size):
	for y in range(y_size):
            if FDindexplane[y,x]==0 or FDindexplane[y,x]==z_size-1:
                output_u[y,x] = uscidata[int(FDindexplane[y,x]),y,x]
                output_q[y,x] = qscidata[int(FDindexplane[y,x]),y,x]
	    elif np.isnan(FDindexplane[y,x]):
		#Q_med, Q_high = qscidata[int(FDindexplane[y,x]),y,x],qscidata[int(FDindexplane[y,x])+1,y,x]
		#U_med, U_high = uscidata[int(FDindexplane[y,x]),y,x],uscidata[int(FDindexplane[y,x])+1,y,x]
		output_q[y,x] = np.nan
		output_u[y,x] = np.nan
            else:
                Q_med, Q_high = qscidata[int(FDindexplane[y,x]),y,x],qscidata[int(FDindexplane[y,x])+1,y,x]
                U_med, U_high = uscidata[int(FDindexplane[y,x]),y,x],uscidata[int(FDindexplane[y,x])+1,y,x]
                output_q[y,x] = Q_med + ((Q_high - Q_med)*(FDindexplane[y,x] - float(int(FDindexplane[y,x]))))
                output_u[y,x] = U_med + ((U_high - U_med)*(FDindexplane[y,x] - float(int(FDindexplane[y,x]))))
            if np.isnan(FDindexplane[y,x]):
                output_pa[y,x] = np.nan
            elif outputfluxarray[y,x] > options.cutoff:
                output_pa[y,x] = np.degrees(0.5*np.arctan2(output_u[y,x],output_q[y,x])) #computing the polarisation angle
            else:
                output_pa[y,x] = np.nan
                

print output_q[418,196]
print output_u[418,196]
print FDindexplane[418,196]

N=len(freq)

#######FIND ERRORS#########


#From equation A.17 from Brentjens and deBruyn (2005)
sig_lamb_square = (1.0/(N-1.0))*(np.sum(np.power(c/freq,4.0))-(1.0/N)*np.power(np.sum(np.power(c/freq,2.0)),2.0 ) )


#Finding error for Faraday Depth


for x in range(x_size):
	for y in range(y_size):
            if np.isnan(FDindexplane[y,x]):
                continue
            elif outputfluxarray[y,x] > options.cutoff:
               sigma_phi[y,x] = 0.5*(570/(outputfluxarray[y,x]/options.rms)) #Here 570 is the FWHM of the RMSF # Therese you will need to change this value to the rmsf of the your observation
            else:
                continue




#Finding error for polarization angle

for x in range(x_size):
	for y in range(y_size):
            if np.isnan(FDindexplane[y,x]):
                continue
            elif outputfluxarray[y,x] > options.cutoff:
                sigma_pa[y,x] = 0.5*(57.2958/(outputfluxarray[y,x]/options.rms)) #Here 57.2958 is 1 radian
            else:
                continue

#Output files here##

headerdataq['NAXIS'] = 2
del headerdataq['NAXIS3']
del headerdataq['CDELT3']
del headerdataq['CRPIX3']
del headerdataq['CRVAL3']
del headerdataq['CTYPE3']
del headerdataq['CUNIT3']
headerdataq['BUNIT'] = 'Jy/beam' 
hdulistq[0].data = output_q
hdulistq.writeto('outputqflux.fits', clobber=options.overwrite)

headerdataq['BUNIT'] = 'Jy/beam'
hdulistq[0].data = output_u
hdulistq.writeto('outputuflux.fits', clobber=options.overwrite)

headerdataq['BUNIT'] = 'degrees'
hdulistq[0].data = output_pa
hdulistq.writeto('outputpa.fits', clobber=options.overwrite)

headerdataq['BUNIT'] = 'rad/m2'
hdulistq[0].data = sigma_phi
hdulistq.writeto('outputsigmaphi.fits', clobber=options.overwrite)

headerdataq['BUNIT'] = 'degrees'
hdulistq[0].data = sigma_pa
hdulistq.writeto('outputsigmapa.fits', clobber=options.overwrite)
