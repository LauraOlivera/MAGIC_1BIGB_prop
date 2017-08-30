################################################################################
# macro for extending the Fermi spectrum to MAGIC energies, accounting for 
# EBL absorption
# 24 Aug 2017
# Author: Cosimo Nigro (cosimonigro2@gmail.com)
################################################################################
from __future__ import division
import numpy as np
import astropy.units as u
from astropy.table import Table
from scipy.interpolate import interp1d, interp2d

import matplotlib.pyplot as plt

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

print MAGIC_diff_sensitivity['energy']

# define the SED def power_law(E, E0, N0, Gamma, z):
def power_law(E, E0, N0, Gamma):
	return E**2 * N0 * (E/E0)**(-Gamma)

# interpolate Dominguez table
tau_Dominguez_11 =  np.loadtxt('tau_dominguez11.out')

z_tau = [ 0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842,    
		  0.16263158, 0.17789474,  0.19315789,  0.20842105, 0.22368421, 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3,  
          0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.,1.2,1.4,1.6,1.8,2.]

E_tau = tau_Dominguez_11[:,0]

tau_Dom_11_interp = interp2d(z_tau, E_tau, tau_Dominguez_11[:,1:], kind='cubic')


# read the info in the table
_1BIGB_table = Table.read('1BIGB_table.dat', format='ascii.basic')

for source in _1BIGB_table:
	N0 = source['N0'] * u.Unit('MeV-1 cm-2 s-1')
	E0 = 1 * u.GeV
	Gamma = source['$\Gamma$']
	z = source['z']
	sens_ene = MAGIC_diff_sensitivity['energy'].to('TeV') # all energies in TeV from now on
	MAGIC_pwl = power_law(sens_ene, E0, N0, Gamma).to('TeV cm-2 s-1') 	
	absorp = np.array([np.exp(-tau_Dom_11_interp(z, _)[0]) for _ in sens_ene.value])
	MAGIC_extended_pwl = MAGIC_pwl * absorp
	# is the flux overcoming the sensitivity?
	greater = np.greater(MAGIC_extended_pwl.value, MAGIC_diff_sensitivity['dfde_sens'].value)
	ene = np.array(sens_ene)
	ene[greater]
	BIGB_count = 0
	if ene[greater].any():
		print 'source overcoming sensitivity at'
		print ene[greater], 'GeV'
		print source['1BIGB source name']
	
		# plot everything
		fig = plt.figure()
		Fermi_ene = np.logspace(np.log10(0.3), np.log10(500), 50) * u.GeV
		plt.plot(Fermi_ene.to('TeV'),  power_law(Fermi_ene.to('TeV'), E0, N0, Gamma).to('TeV cm-2 s-1'), ls='-', lw=1.5, marker='', color='teal', label='1BIGB spectrum')
		plt.plot(sens_ene, MAGIC_extended_pwl, ls='--', lw=1.5, marker='', color='teal', label='extrapolation + EBL absorption')
		plt.plot(sens_ene, MAGIC_diff_sensitivity['dfde_sens'],  ls='--', lw=1.5, marker='', color='k', label='MAGIC sensitivity')
		plt.title(source['1BIGB source name'] + ', z=' + str(z) + '\n' + str(source['# of hours Zenith 5-50']) + ' dark hours')
		plt.legend(loc='lower left')
		plt.xlabel('E [TeV]', size=12)
		plt.ylabel(r'E$^2$dF/dE [TeV cm$^{-2}$ s$^{-1}$]')	
		plt.xscale('log')
		plt.yscale('log')
		plt.ylim([1e-14, 1e-10])
		plt.show()
 		fig.savefig('1BIGB_' + str(BIGB_count) + '.png')
		BIGB_count += 1
