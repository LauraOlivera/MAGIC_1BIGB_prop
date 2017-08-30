################################################################################
# class for visibility studies for the MAGIC telescope
# 23 Aug 2017
# Authors: * Laura Olivera Nieto 
# 		   * Cosimo Nigro (cosimonigro2@gail.com)
################################################################################

from __future__ import division
import numpy as np
import numpy.ma as ma
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
from astropy import units as u
from scipy.interpolate import interp1d, interp2d

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.pyplot import axes
import matplotlib.dates as mdates
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'


class MAGIC_visibility:
	"""Class for visibility studies with the MAGIC telescopes contains:
	* functions for producing visibility plots
	* functions to make comparison with Fermi spectra and MAGIC sensitivity	
	"""
	
	def __init__(self, name, ra, dec, z, gamma, N0, E0, TS):
		"""The instance of the class will have name of the source in the catalogue, 
		RA and DEC, redshift, Fermi spectral info, array of times related to MAGIC 
		observatory.
		
		Parameters:
		-----------
		name : string
			source name in the catalogue
		ra : `~astropy.units.Quantity`
			right ascension of the source in deg
		dec : `~astropy.units.Quantity`
			declination of the source in deg
		z : float
			redshift of the source
		gamma : float
			Fermi spectral index (PowerLaw)
		N0 : `~astropy.units.Quantity`
			normalization of Fermi spectrum (PowerLaw)
		E0 : `~astropy.units.Quantity`
			pivot energy of Fermi spectrum (PowerLaw)
		TS : float
			TS in Fermi
		"""
		self.name = name
		self.ra = ra
		self.dec = dec 
		self.z = z
		self.gamma = gamma
		self.N0 = N0
		self.E0 = E0
		self.TS = TS

		# define MAGIC Location
		self.MAGIC_loc = EarthLocation(lat=28.76194*u.deg, lon=-17.89*u.deg, height=2200*u.m)
		# define the sky coordinates
		self.sky_coord = SkyCoord(ra = self.ra, dec = self.dec, frame = 'icrs')
		# define the list of times for 2018, we will do everything in UTC time
		# we will consider 20 minutes run from midnight of new years's eve spanning all the year
		start_year = Time('2017-12-31T00:00:00', format='isot', scale='utc')
		thirty_mins = np.linspace(0, 365*24*60, 365*24*2) * u.min # list of 30 minutes run until the end of the year
		self.time_year = start_year + thirty_mins 
	
	
	def zenith_year(self):
		"""function for obtaining a matrix of the zenith of the source during the year with values
		masked for Sun altitude above -18 deg and Moon altitude above 0
		Returns:
		--------
		`~numpy.ndarray`
			matrix of shape 48 x 365 with x axis day of the year, y axis hour of the day
		"""
		# yearly alt-azimuth frame
		altaz_frame = AltAz(obstime = self.time_year, location = self.MAGIC_loc)
		# Moon, Sun and source coordinates in the yearly alt-azimuth frame
		moon_altaz = get_moon( self.time_year ).transform_to( altaz_frame )
		sun_altaz = get_sun( self.time_year ).transform_to( altaz_frame )
		source_altaz = self.sky_coord.transform_to( altaz_frame )  
		# reshape the source_altaz array into a matrix (48,365), take the altitudes
		# first we select the values 48 to 48, obtaining 365 rows, then we transpose it
		source_alt = np.reshape(source_altaz.alt.value, (365, 48)).T
		moon_alt = np.reshape(moon_altaz.alt.value, (365, 48)).T
		sun_alt =  np.reshape(sun_altaz.alt.value, (365, 48)).T
		# turn the source altitude into zenith
		ninety = 90*np.ones((48,365))
		source_Zd = ninety - source_alt
		# mask where the sun is above -18 (astronomical sunrise) and moon above the horizon 
		source_sun_mask = ma.masked_where(sun_alt >= -18, source_Zd)
		source_moon_mask = ma.masked_where(moon_alt >= 0., source_Zd)
		source_mask = source_sun_mask + source_moon_mask
		
		return source_mask


	def plot_zenith_year(self, Zd_year_matrix):
		"""plotting the visibility for all the year in terms of zenith"""
		fig = plt.figure(figsize=(14,8))
		plt.imshow(Zd_year_matrix, origin='lower', cmap=cm.viridis_r, aspect='auto', interpolation='none',
                   vmin=np.min(Zd_year_matrix), vmax=np.max(Zd_year_matrix))
		cb = plt.colorbar(orientation='vertical')
		cb.set_label('Zenith [deg]', size=13)
		plt.xticks(np.array([0,32,60,91,121,152,182,213,244,274,305,335]),
                            ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
		plt.yticks(np.arange(0,48,2),
                            ('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
                             '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'))  
		
		anyArtist1 = plt.Line2D((0,1),(0,0), color='crimson', ls='-', linewidth=2)
		anyArtist2 = plt.Line2D((0,1),(0,0), color='crimson', ls='--', linewidth=2)
		# add the contours
		cs1 = plt.contour(Zd_year_matrix, [35], origin='lower', linewidths=2, linestyles='-', colors='crimson')
		cs2 = plt.contour(Zd_year_matrix, [50], origin='lower', linewidths=2, linestyles='--', colors='crimson')

		plt.legend([anyArtist1, anyArtist2],['35 deg', '50 deg'],loc='center left') 
		plt.title(self.name + ' all-year visibility')
		plt.xlabel('Months', size=13)  
		plt.ylabel('Hours (UTC)', size=13)
		plt.show()
		fig.savefig(self.name + '_yearly_visibility.png')
		

	def visible_hours(self, Zd_year_matrix):
		"""histogram with the number of hours visible in each zenith range
		Parameters:
		Zd_year_matrix : `~numpy.ndarray`
			yearly matrix with the zenith rangs output of zenith_year
		Returns:
		--------
		dictionary
			dictionary with the zenith range (string) and the number of hours the source spend inside
		"""	
		h = np.zeros(4)
		for i in range(48):
			for j in range(365):
					if 5 < Zd_year_matrix[i,j] < 35:
						h[0]+=1
					elif 35 < Zd_year_matrix[i,j] < 50:
						h[1]+=1
					elif 50 < Zd_year_matrix[i,j] < 62:
						h[2]+=1
					elif 62 < Zd_year_matrix[i,j] < 70:
						h[3]+=1
		h = h/2
		visible_hours = dict(Zd_05_35=h[0], Zd_35_50=h[1], Zd_50_62=h[2], Zd_62_70=h[3])
		return visible_hours


	def hist_visible_hours(self, visible_hours):	
		"""make an histogram of the visible hours"""
		fig = plt.figure()
		plt.fill_between([05,35], [0,0], [visible_hours['Zd_05_35'], visible_hours['Zd_05_35']], 
						 alpha=0.5, hatch='/', color = 'teal', label='5 < Zd < 35 [deg]')
		plt.fill_between([35,50], [0,0], [visible_hours['Zd_35_50'], visible_hours['Zd_35_50']], 
						 alpha=0.5, hatch='/', color = 'mediumspringgreen', label='35 < Zd < 50 [deg]')
		plt.fill_between([50,62], [0,0], [visible_hours['Zd_50_62'], visible_hours['Zd_50_62']], 
						 alpha=0.5, hatch='/', color='goldenrod', label='50 < Zd < 62 [deg]')
		plt.fill_between([62,70], [0,0], [visible_hours['Zd_62_70'], visible_hours['Zd_62_70']], 
						 alpha=0.5, hatch='/', color='darkorange', label='62 < Zd < 70 [deg]')
	
		plt.legend(loc='lower left')
		plt.xlim([0,80])
		plt.xlabel('Zenith [deg]', size=12)
		plt.ylabel('number of hours', size=12)
		plt.title(self.name + ', hours per zenith bin')
		plt.show()
		fig.savefig(self.name + '_hours_Zd_bin_histo.png')


	def MAGIC_sensitivity_checks(self):
		"""check if a simple extension of the source spectrum + EBL absorption overcome the sensitivity
		we are simply extending the Fermi spectrum to VHE and considering the absorption by ebl.
		"""
		# define the MAGIC differential sensitivity
		sens_from_string = """E0[GeV] E1[GeV] E2[GeV], sens[%Crab], sens_error[%Crab], sens[TeV cm^-2 s^-1], sens_error[TeV cm^-2 s^-1]
50	39.7	63.	26	2.7	1.79e-11	1.86e-12
79.4	63.	100.	7.1	0.3	4.87e-12	2.06e-13
125.8	100.	158.4	3.87	0.11	2.55e-12	7.25e-14
199.4	158.4	251.1	2.38	0.08	1.44e-12	4.86e-14
316.1	251.1	397.9	2.16	0.16	1.16e-12	8.64e-14
501.	397.9	630.7	1.61	0.18	7.43e-13	8.31e-14
794	630.7	1000.	1.64	0.09	6.22e-13	3.41e-14
1258.	1000.	1584.	1.8	0.09	5.40e-13	2.70e-14
1994.	1584.	2511.	2.6	0.19	5.93e-13	4.33e-14
3161.	2511.	3979.	4.3	0.4	7.18e-13	6.68e-14
5010.	3979.	6307.	6.8	0.7	7.99e-13	8.23e-14
7940	6307.	10000.	16	3	1.27e-12	2.38e-13
12500 	10000 	15736. 	37 	10 	1.89e-12 	4.92e-13"""

		# use a dictionary with astropy units for the sensitivity
		MAGIC_diff_sensitivity = dict()
		sens_rows = [row.split('\t') for row in sens_from_string.split('\n')[1:]]
		MAGIC_diff_sensitivity['energy'] = [float(info[0]) for info in sens_rows] * u.GeV
		MAGIC_diff_sensitivity['dfde_sens'] = [float(info[-2]) for info in sens_rows] * u.Unit('TeV cm-2 s-1')
		MAGIC_diff_sensitivity['dfde_sens_err'] = [float(info[-1]) for info in sens_rows] * u.Unit('TeV cm-2 s-1')

		def power_law(E, E0, N0, gamma):
			"""power-law SED"""
			return E**2 * N0 * (E/E0)**(-gamma)

		# interpolate Dominguez table
		tau_Dominguez_11 =  np.loadtxt('tau_dominguez11.out')

		z_tau = [ 0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 
				  0.10157895, 0.11684211, 0.13210526, 0.14736842, 0.16263158, 0.17789474,  
				  0.19315789,  0.20842105, 0.22368421, 0.23894737, 0.25421053, 0.26947368, 
				  0.28473684, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 
				  0.85, 0.9, 0.95, 1.,1.2,1.4,1.6,1.8,2.]

		E_tau = tau_Dominguez_11[:,0]
		tau_Dom_11_interp = interp2d(z_tau, E_tau, tau_Dominguez_11[:,1:], kind='cubic')

		# make the comparison plot
		sens_ene = MAGIC_diff_sensitivity['energy'].to('TeV') # all energies in TeV from now on
		MAGIC_pwl = power_law(sens_ene, self.E0, self.N0, self.gamma).to('TeV cm-2 s-1') 	
		absorp = np.array([np.exp(-tau_Dom_11_interp(self.z, _)[0]) for _ in sens_ene.value])
		MAGIC_extended_pwl = MAGIC_pwl * absorp
		# is the flux overcoming the sensitivity?
		greater = np.greater(MAGIC_extended_pwl.value, MAGIC_diff_sensitivity['dfde_sens'].value)
		ene = np.array(sens_ene)
		ene[greater]
		if ene[greater].any():
			print 'source overcoming sensitivity at'
			print ene[greater], 'GeV'
	
		# plot everything
		fig = plt.figure()

		Fermi_ene = np.logspace(np.log10(0.3), np.log10(500), 50) * u.GeV
		plt.plot(Fermi_ene.to('TeV'),  
				 power_law(Fermi_ene.to('TeV'), self.E0, self.N0, self.gamma).to('TeV cm-2 s-1'), 
				 ls='-', lw=1.5, marker='', color='teal', label='1BIGB catalog spectrum')

		plt.plot(sens_ene, MAGIC_extended_pwl, 
				 ls='--', lw=1.5, marker='', color='teal', label='extrapolation + EBL absorption')
		
		plt.plot(sens_ene, MAGIC_diff_sensitivity['dfde_sens'],  
				 ls='--', lw=1.5, marker='', color='k', label='MAGIC sensitivity')

		plt.title(self.name + ' sensitivity comparison \n z=' + str(self.z))
		plt.legend(loc='lower left')
		plt.xlabel('E [TeV]', size=12)
		plt.ylabel(r'E$^2$dF/dE [TeV cm$^{-2}$ s$^{-1}$]', size=12)	
		plt.xscale('log')
		plt.yscale('log')
		plt.ylim([1e-14, 1e-10])
		plt.show()
 		fig.savefig(self.name + '_MAGIC_sensitivity_comparison.png')
	
