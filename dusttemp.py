import numpy as np
import pyfits
import optparse
import sys
import pylab as pyl
import matplotlib.pyplot as plt
from scipy import stats

#define constants etc
beta = 2
spitzerlambda = 7e-5
spitzerfreq = 4.2827494e+12 #70 micros
hershellambda = 1.6e-4
hershelfreq = 1.8737e+12 #160 micros
c = 299792458 
hconstant = 6.626e-34
k = 1.38e-23
h = 0.1
initialT = 25
Te4 = 1 # electron temp in units of 10^4
fillfactor = 0.33 # filling factor
a = 1

#finding terms in planck function
freqfrac = (spitzerfreq**beta/hershelfreq**beta) #okay
outerfrachersc = (2*hconstant*(hershelfreq**3))/(c**2) #okay
innerfrachersc = (hconstant*hershelfreq)/k #defo okay
outerfracspitzer = (2*hconstant*(spitzerfreq**3))/(c**2) #okay
innerfracspitzer = (hconstant*spitzerfreq)/k


def derivative(f, x, h,fluxhers,fluxspitzer):
      return (f(x+h,fluxhers,fluxspitzer) - f(x-h,fluxhers,fluxspitzer)) / (2.0*h)

def planckfunc(T,fluxhers,fluxspitzer):
    return (freqfrac*(outerfracspitzer/outerfrachersc)*(np.exp(innerfrachersc/T)-1)/(np.exp(innerfracspitzer/T)-1)) - (fluxspitzer/fluxhers)

def newtonraph(fluxhers,fluxspitzer):
    lastT = initialT
    nextT = lastT + 10*h
    while (np.abs(lastT-nextT)>h):
        newY = planckfunc(nextT,fluxhers,fluxspitzer)
        #print "f(", nextT, ") = ", newY
        lastT = nextT
        nextT = lastT - newY/derivative(planckfunc, lastT, h ,fluxhers,fluxspitzer)
    return nextT

def opticaldepth(irflux,t):
    #
    return -1*np.log(1 - (irflux*((1e6*1e-26*c**2)/(2*hconstant*(hershelfreq**3))))*(np.exp(innerfrachersc/t)-1))

def findemissionmeasure(deredHA,hapixelsize):
    #finding the Emission measure 
    deredHAsr = deredHA*(206265**2/hapixelsize**2)
    EM = deredHAsr/((9.41e-8*(Te4**-1.017))*(10**(-0.029/Te4)))
    return EM

def continuumopt(EM,freq):
    #finding the optical depth in the contiuum 
    return 8.235e-2*a*(Te4**-1.35)*(freq**-2.1)*(1.08)*EM 

def findtherbright(contim):
    return Te4(1-np.exp(-contim))

def findthermalflux(tb,bmaj,bmin,freq):
    #finding the thermal flux
    return tb/((1.222e6)*(bmaj**-1)*(bmin**-1)*(freq**-2))

#parsers are located here

parser = optparse.OptionParser()

required = optparse.OptionGroup(parser, "Required Attributes")

required.add_option('--sp', help='Input SPITZER MAP (70 microns)', action='store', type='string', dest='input1')
required.add_option('--he', help='Input HERSCHEL MAP (240 microns)', action='store', type='string', dest='input2')
required.add_option('--ha', help='Input H-alpha MAP ', action='store', type='string', dest='inputha')
required.add_option('--co', help='Input continuum MAP ', action='store', type='string', dest='inputcon')
required.add_option('--cutoff', help='cutoff for finding dust temperature', action='store', type='float', dest='cutoff')
required.add_option('--o', help='Outputted color temp file', action='store', type='string', dest='outputfits')
parser.add_option("--ov", help='Overwrite existing fits files', action="store_true", dest="overwrite",default=False)

parser.add_option_group(required)

options, arguments = parser.parse_args()

#Need to read in frequency in GHZ of continuum  as well as beam size and pixel size of HA


try:
    hdulist1 = pyfits.open(options.input1)
    headerdata1 = hdulist1[0].header
    scidata1 = hdulist1[0].data
    x_size = scidata1.shape[-1]
    y_size = scidata1.shape[-2]
except IOError:
    print('ERROR.SPITZER MAP is not found! Please check your input.')
    sys.exit(0)
try:
    hdulist2 = pyfits.open(options.input2)
    headerdata2 = hdulist2[0].header
    scidata2 = hdulist2[0].data
    x_size2 = scidata2.shape[-1]
    y_size2 = scidata2.shape[-2]
except IOError:
    print('ERROR.HERSCHEL MAP is not found! Please check your input.')
    sys.exit(0)
try:
    hdulist3 = pyfits.open(options.inputha)
    headerdata3 = hdulist3[0].header
    scidata3 = hdulist3[0].data
    x_size3 = scidata3.shape[-1]
    y_size3 = scidata3.shape[-2]
except IOError:
    print('ERROR.H-alpha MAP is not found! Please check your input.')
    sys.exit(0)
try:
    hdulist4 = pyfits.open(options.inputcon)
    headerdata4 = hdulist4[0].header
    scidata4 = hdulist4[0].data
    x_size4 = scidata4.shape[-1]
    y_size4 = scidata4.shape[-2]
except IOError:
    print('ERROR.Contiuum MAP is not found! Please check your input.')
    sys.exit(0)

#checking size of maps and pixels are equal

if x_size!=x_size2!=x_size3!=x_size4:
    print ('Dimensions of maps are not equal! Program closing')
    sys.exit(0)

hapixelra = headerdata1['CDELT1']*3600
hapixeldec = headerdata1['CDELT2']*3600

if np.abs(hapixelra)!=hapixeldec:
    print ('Pixel of HII is not square! Program closing')
    sys.exit(0)

bmin = headerdata4['BMIN']*3600 #needs to be arcsec
bmaj = headerdata4['BMAJ']*3600 #needs to be arcsec
freq = headerdata4['CRVAL3']/1e9 #needs to be in GHz


###Define dust map output
outputdustarray = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

###Finding the dust temperature
for x in range(x_size):
	for y in range(y_size):
		if scidata1[0,y,x] > options.cutoff:
			outputdustarray[y,x] = newtonraph(scidata2[0,y,x],scidata1[0,y,x])
		else:
			outputdustarray[y,x] = np.nan

dustaverage = np.nanmean(outputdustarray[208:394,168:390])
dustmedian = stats.nanmedian(outputdustarray[208:394,168:390],axis=None)
duststd = np.nanstd(outputdustarray[208:394,168:390])
print (dustaverage)
print (dustmedian)
print (duststd)
plt.show()

#Finding the optical depth
opticaldeptharray = np.nan*np.zeros((y_size,x_size), dtype=np.float32)
for x in range(x_size):
	for y in range(y_size):
		if scidata1[0,y,x] > options.cutoff:
			opticaldeptharray[y,x] = opticaldepth(scidata2[0,y,x],outputdustarray[y,x])
		else:
			opticaldeptharray[y,x] = np.nan


haopticalarray = opticaldeptharray*2200*fillfactor

#Finding the dereddened Ha

deredHA = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

for x in range(x_size):
	for y in range(y_size):
		if np.isnan(haopticalarray[y,x]):
			deredHA[y,x] = scidata3[y,x]
                        print (scidata3[y,x])
		else:
                        deredHA[y,x] = scidata3[y,x]*np.exp(-1*haopticalarray[y,x])

EMarray = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

for x in range(x_size):
	for y in range(y_size):
            EMarray[y,x] = findemissionmeasure(deredHA[y,x],hapixeldec) #need to define hapixel


#FINDING coontinuum optical depth

coonopt = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

for x in range(x_size):
	for y in range(y_size):
            coonopt[y,x] = continuumopt(EMarray,freq) #need to define frequency


#Finding the brightness temperature

brighttemp = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

for x in range(x_size):
	for y in range(y_size):
            brighttemp[y,x] = findtherbright(coonopt) #need to define frequency

#Finding the thermal Flux

therflux = np.nan*np.zeros((y_size,x_size), dtype=np.float32)

for x in range(x_size):
	for y in range(y_size):
            therflux[y,x] = findthermalflux(brighttemp,bmaj,bmin,freq)


#Finding the non thermal image

nonthermal = scidata4 - therflux


#Outputting all Files

headerdata1['NAXIS'] = 2
del headerdata1['NAXIS3']
del headerdata1['CDELT3']
del headerdata1['CRPIX3']
del headerdata1['CRVAL3']
del headerdata1['CTYPE3']
del headerdata1['CUNIT3']
headerdata1['BUNIT'] = 'Kelvin' 
hdulist1[0].data = outputdustarray
hdulist1.writeto(options.outputfits, clobber=options.overwrite)

headerdata1['BUNIT'] = '' 
hdulist1[0].data = opticaldeptharray
hdulist1.writeto('dustopticaldepth160micron.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = '' 
hdulist1[0].data = haopticalarray
hdulist1.writeto('dustopticaldepthha.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = '' 
hdulist1[0].data = deredHA
hdulist1.writeto('dereddenHA.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = '' 
hdulist1[0].data = EMarray
hdulist1.writeto('EMISSONMEASURE.fits', clobber=options.overwrite)


headerdata1['BUNIT'] = '' 
hdulist1[0].data = coonopt
hdulist1.writeto('contiuumdepth.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = '' 
hdulist1[0].data = brighttemp
hdulist1.writeto('brighttemp.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = 'Jy/beam' 
hdulist1[0].data = therflux
hdulist1.writeto('thermalmap.fits', clobber=options.overwrite)

headerdata1['BUNIT'] = 'Jy/beam' 
hdulist1[0].data = nonthermal
hdulist1.writeto('nonthermalmap.fits', clobber=options.overwrite)
