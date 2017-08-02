""" What does this thing do """

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_moon
from astropy.coordinates import get_sun
from astropy import units as u

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.pyplot import axes
import matplotlib.dates as mdates
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

import datetime

import ephem

from datetime import datetime, timedelta, date, time
import timeit

class MAGIC_visibility:
	
	def __init__(self, name, ra, dec):
		self.name = name
		self.ra = ra
		self.dec = dec 
		# define MAGIC Location
		self.MAGIC_loc = EarthLocation(lat=28.76194*u.deg, lon=-17.89*u.deg, height=2200*u.m)
		self.utcoffsetSummer = +1*u.hour
		self.utcoffsetWinter = 0*u.hour 

		# define the time array
		#midnight = Time('2017-12-31 00:00:00') + utcoffsetWinter
		#delta_midnight = np.linspace(0,525600,17520)*u.min
		#times = midnight + delta_midnight
		#days = np.split(times,365)




# Obtain the position of the moon, sun and object on a given day
	def day(self, date=None):		    		
		if date is None:
			date = Time.now()
		delta_midnight = np.linspace(0, 24, 100)*u.hour
		times = date + delta_midnight
		BIGB = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
		# MAGIC Horizontal Coordinate system
		altazframe = AltAz(obstime=times, location=self.MAGIC_loc)
		# Moon's position in it
		moonaltazs = get_moon(times).transform_to(altazframe)
		moon_alt = np.array(moonaltazs.alt)
		moon_az = np.array(moonaltazs.az)
		# Sun's position in it
		sunaltazs = get_sun(times).transform_to(altazframe)
		sun_alt = np.array(sunaltazs.alt)
		sun_az = np.array(sunaltazs.az)
		# Object's position in it
		BIGB_altazs = BIGB.transform_to(altazframe)  
		BIGB_alt = np.array(BIGB_altazs.alt)
		BIGB_az = np.array(BIGB_altazs.az)
		# Angular distance to the moon
		x = range(len(BIGB_altazs.alt))
		moon_dist = [np.arccos(np.sin(moon_alt[i]* np.pi / 180.)*np.sin(BIGB_alt[i]* np.pi / 180.)+ \
			     np.cos(moon_alt[i]* np.pi / 180.)*np.cos(BIGB_alt[i]* np.pi / 180.)*np.cos(moon_az[i]* np.pi/180.-BIGB_az[i]* np.pi / 180.)) for i in x]
		# Convert to degrees
		x = range(len(BIGB_alt))
		moon_dist=[(180/np.pi)*moon_dist[i] for i in x]
		moon_dist = np.array(moon_dist)
		return BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist

# Visibility plot of the object on a given day. Includes moon and sun altitude.
	def plot_day(self, date=None):
		delta_midnight = np.linspace(0, 24, 100)*u.hour
		BIGB = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
		BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist =  self.day(date)
		fig = plt.figure(figsize=(14, 8.4))
		plt.plot(delta_midnight, moon_alt, color='g', label='Moon')  
		plt.plot(delta_midnight, sun_alt, color='y', label='Sun')  
		plt.scatter(delta_midnight, BIGB_alt, c=BIGB_az, label=self.name, lw=0, s=8, zorder=1)  
		# plot a gray band corresponding to dark time i.e. negative altitude of the sun
		plt.fill_between(delta_midnight.to('hr').value, 0, 90, sun_alt < -18, color='0.2', zorder=0)  
		#plt.colorbar().set_label('Azimuth [deg]')  
		plt.title("Visibility of %s" %self.name)
		plt.legend(loc='upper left')  
		plt.xlim(0, 24)  
		plt.xticks(np.arange(0,24,1))  
		plt.ylim(0, 90)  
		plt.xlabel('Hours from start')  
		plt.ylabel('Altitude [deg]')  
		return plt.show()


	
# Obtain the position of the moon, sun and object a year from the starting date as well as the angular distance to the moon[deg]
	def year(self, date=None):
		if date is None:
			date = Time.now()
		delta_midnight = np.linspace(0,525600,17520)*u.min
		times = date + delta_midnight
		BIGB = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
		# MAGIC Horizontal Coordinate system
		altazframe = AltAz(obstime=times, location=self.MAGIC_loc)
		# Moon's position in it
		moonaltazs = get_moon(times).transform_to(altazframe)
		moon_alt = np.array(moonaltazs.alt)
		moon_az = np.array(moonaltazs.az)
		# Sun's position in it
		sunaltazs = get_sun(times).transform_to(altazframe)
		sun_alt = np.array(sunaltazs.alt)
		sun_az = np.array(sunaltazs.az)
		# Object's position in it
		BIGB_altazs = BIGB.transform_to(altazframe)  
		BIGB_alt = np.array(BIGB_altazs.alt)
		BIGB_az = np.array(BIGB_altazs.az)
		# Angular distance to the moon
		x = range(len(BIGB_altazs.alt))
		moon_dist = [np.arccos(np.sin(moon_alt[i]* np.pi / 180.)*np.sin(BIGB_alt[i]* np.pi / 180.)+ \
			     np.cos(moon_alt[i]* np.pi / 180.)*np.cos(BIGB_alt[i]* np.pi / 180.)*np.cos(moon_az[i]* np.pi/180.-BIGB_az[i]* np.pi / 180.)) for i in x]
		# Convert to degrees
		x = range(len(BIGB_alt))
		moon_dist=[(180/np.pi)*moon_dist[i] for i in x]
		moon_dist = np.array(moon_dist)
		return BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist


# Make the year data into a matrix and plot it as a visibility plot
	def plot_year(self, date=None):
		delta_midnight = np.linspace(0,525600,17520)*u.min
		# Obtain the data from the year function
		BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist =  self.year(date)
		# Split the data into a matrix of 365 columns and 48 rows (each entry represents half an hour)
		# First we split the days
		BIGB_spl = np.split(BIGB_alt,365)
		moon_spl = np.split(moon_alt,365)
		sun_spl = np.split(sun_alt,365)
		moondist_spl = np.split(moon_dist,365)
		# Create the empty matrices
		B_alt=np.zeros((48,365))
		M_alt=np.zeros((48,365))
		S_alt=np.zeros((48,365))
		M_dist=np.zeros((48,365))
		# Arrange the altitude data into the matrices
		for i in range(365):
	   		 B_alt[:,i]=BIGB_spl[i]
	    		 M_alt[:,i]=moon_spl[i]
			 S_alt[:,i]=sun_spl[i]
			 M_dist[:,i]=moondist_spl[i]
		# Mask the entries where the sun is above the horizon (>-18 deg)
		S1=ma.masked_where(S_alt <= -18, S_alt)
		# Make the rest zero for the plot (we could also make them zero directly on the object)
		np.putmask(S1, S1>-20, 0)
		# Make zero the entries where the moon is below the horizon (dark times)
		M_dist=ma.masked_where(M_alt<=0, M_dist)
		# Make zero when the moon is too close or too far
		np.putmask(M_dist, M_alt<30, 0)
		np.putmask(M_dist, M_alt>110, 0)
		fig = plt.figure(figsize=(14,8))
		plt.imshow(B_alt,origin='lower',cmap=cm.magma, aspect='auto',zorder=0)
		plt.colorbar().set_label('Altitude [deg]')
		plt.imshow(M_dist,origin='lower',cmap=cm.binary, aspect='auto',zorder=1)
		plt.colorbar().set_label('Angular distance to the moon [deg]')
		plt.imshow(S1,origin='lower', cmap=cm.binary, aspect='auto',zorder=2)
 		#plt.colorbar().set_label('Altitude [deg]')
		plt.title("Visibility of %s" %self.name)
		plt.xticks(np.array([0,32,60,91,121,152,182,213,244,274,305,335]),('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
		plt.yticks(np.arange(0, 48,2),('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'))   
		plt.xlabel('Days')  
		plt.ylabel('Hours (UTC)')
		return plt.show()



# Obtain the number of visible hours of a given object in a year starting from date

	def visible_hours(self, date=None):
		delta_midnight = np.linspace(0,525600,17520)*u.min
		# Obtain the data from the year function
		BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist =  self.year(date)
		# Split the data into a matrix of 365 columns and 48 rows (each entry represents half an hour)
		# First we split the days
		BIGB_spl = np.split(BIGB_alt,365)
		moon_spl = np.split(moon_alt,365)
		sun_spl = np.split(sun_alt,365)
		moondist_spl = np.split(moon_dist,365)
		# Create the empty matrices
		B_alt=np.zeros((48,365))
		M_alt=np.zeros((48,365))
		S_alt=np.zeros((48,365))
		M_dist=np.zeros((48,365))
		# Arrange the altitude data into the matrices
		for i in range(365):
	   		 B_alt[:,i]=BIGB_spl[i]
	    		 M_alt[:,i]=moon_spl[i]
			 S_alt[:,i]=sun_spl[i]
			 M_dist[:,i]=moondist_spl[i]
		# Make zero the entries where the sun is too high
		np.putmask(B_alt, S_alt>-18, 0)
		# And the moon too close or too far
		np.putmask(B_alt, M_dist<30, 0)
		np.putmask(B_alt, M_dist>110, 0)
		# Obtain the number of hours in each angle range
		h = np.zeros(4)
		for i in range(48):
			for j in range(365):
        			if B_alt[i,j]>20 and B_alt[i,j]<28:
           				 h[0]=h[0]+1
       				elif B_alt[i,j]>28 and B_alt[i,j]<40:
           				 h[1]=h[1]+1
        			elif B_alt[i,j]>40 and B_alt[i,j]<55:
           				 h[2]=h[2]+1
        			elif B_alt[i,j]>55 and B_alt[i,j]<85:
           				 h[3]=h[3]+1
            
            	h = h/2
		return h




