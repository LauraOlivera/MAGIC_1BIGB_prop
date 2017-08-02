import MAGIC_visibility as MAGIC_vis
import timeit
from astropy.time import Time

BIGB = MAGIC_vis.MAGIC_visibility('BIGBwhatever', 70.66917,61.67750)
#print(BIGB.ra)
#print(BIGB.day(date=Time('2016-12-31T00:00:00.00', format='isot', scale='utc')))



start = timeit.default_timer()
#BIGB.plot_day(date=Time('2017-02-01T00:00:00.00', format='isot', scale='utc'))
#BIGB.plot_year(date=Time('2016-12-31T00:00:00.00', format='isot', scale='utc'))
#BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist   = BIGB.year(date=Time('2016-12-31T00:00:00.00', format='isot', scale='utc'))
#BIGB_alt, BIGB_az, moon_alt, sun_alt, moon_dist  = BIGB.day(date=Time('2016-12-31T00:00:00.00', format='isot', scale='utc'))
h=BIGB.visible_hours(date=Time('2017-12-31T00:00:00.00', format='isot', scale='utc'))
stop = timeit.default_timer()
print 'Computational time: '
print (stop - start)/60., ' mins'

#print moon_dist
#print BIGB_alt
print h




