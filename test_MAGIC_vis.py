import MAGIC_visibility as vis
import astropy.units as u

BIGB = vis.MAGIC_visibility('1BIGBJ194356.2+211821', 295.98417*u.deg, 21.30611*u.deg, 
                                          0.224, 1.43, 1.63*1e-13 * u.Unit('MeV-1 cm-2 s-1'), 1*u.GeV, 128)
Zd_matrix = BIGB.zenith_year()
BIGB.plot_zenith_year(Zd_matrix)

hours_dict = BIGB.visible_hours(Zd_matrix)
BIGB.hist_visible_hours(hours_dict)

BIGB.MAGIC_sensitivity_checks()
