# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 14:36:29 2016

@author: narsi
"""
import numpy as np
import matplotlib.pyplot as plt

import astropy
from astropy.io import fits

from scipy import signal
from scipy import misc


IM1=fits.open('template.fits')
template=IM1[0].data

IM2=fits.open('image.fits')
image=IM2[0].data

corr = signal.correlate2d(image, template, boundary='fill', mode='same')


hdu = fits.PrimaryHDU(corr)
hdulist = fits.HDUList([hdu])
hdulist.writeto('corr.fits', clobber=1)


plt.imshow(corr)