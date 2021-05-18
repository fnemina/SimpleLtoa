"""
This code runs a simple simulation using OSOAA for a number of bands.
A single ocean configuration is used and the view, zenith and azimuth
angle for each band can be the same of different.

The simulation SIMULATIONS CONFIG section let's us introduce the surface, ocean
and atmospheric parameters to iterate over. Also this includes the angular 
parameters. The output name of the file with the configuration, results
and a FORTRAN readable output as well as the plot name is found there.

This code interpotales the view zenith angle variable to be faster,
however several tests have shown than the error introduced is not big.
"""

###################################################################
# SIMULATION CONFIG
###################################################################
# Bands wavelength - repeated numbers can be used
wl0 = [412.,  443.,  490.,  510.,  555.,  620.,  665.,  680.,
       710.,  750.,  765.,  865., 1044., 1240., 1640.]

# geometric
sun = 35    # Sun zenith angle
view = 0    # View zenith angle - can be array 
phi = 90    # Relative azimuth angle - can be array
            # 180 deg Sat and Sun in the same plane - backscatter
            # 0 deg   Sat and Sun opposite plane    - specular

# Atmosphere
pres=1013.2  # [hPa] set to 100 for no Rayleigh scattering
aot = 0.0    # aerosol optical thickness as 865
rh = 90      # Relative humidity
sfmodel = 3  # Shettle and Fenn
             # 1 : Tropospheric ,
             # 2 : Urban
             # 3 : Maritime,
             # 4 : Coastal

# Surface
wind = 0     # m/s - includes whitecaps over 6.33m/s

# Ocean
chl = 4      # [mg/m3]
sed = 0      # [mg/l]
ays = True   # 1/m - True to autocompute
sys = 0.017  # 1/nm
adet = 0     # 1/m
sdet = 0     # 1/nm

# Level
level = 1    # 1 - satellite level  
             # 3 - sea surface 0+

# Cleanup: 
# set to True to delete temporary files
# set to False to keep temporary files
cleanup = False

# Plot name - None for no plot
plotname = "fig:atm_cor_tmp"

# Output filename - None for no file
filename = "simulation.txt"
###################################################################
# END SIMULATION CONFIG
###################################################################


###################################################################
# IMPORTS AND LUTS
###################################################################
# General code
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Radiative transfer codes
import pyOSOAA

# Load ascii tables for simulations
models = {1 : r"Tropospheric S&F",
          2 : r"Urban S&F",
          3 : r"Maritime S&F",
          4 : r"Coastal S&F"}

# Whitecaps
def rho_wc(wl, wind):
    """ This function computes the whitecaps reflectance albedo for all
        SABIA-Mar bands for an specific wind speed.

        wl      Wavelength [nm]
        wind    Wind speed [m/s]
    """
    awc_wl = [412.0, 443.0, 490.0, 510.0, 555.0, 620.0, 665.0, 680.0, 
              710.0, 750.0, 765.0, 865.0, 1044.0, 1240.0, 1640.0]
    awc_rf = [1., 0.99868, 0.9975, 0.99759, 0.99079, 0.9574, 0.94296,
              0.92733, 0.86586, 0.7261, 0.73048, 0.64952, 0.50294,
              0.40446, 0.25752]
    wac_int = np.interp(wl, awc_wl, awc_rf, left=1, right=0)
    if wind <= 6.33:
        return 0*wl
    if wind >= 12:
        return wac_int*1.925*1e-5*(12.0-6.33)**3
    else:
        return wac_int*1.925*1e-5*(wind-6.33)**3


def taur(wl, pres):
    """ This function interpolates the Rayleigh optical thickness for a 
        certain wavelength using the Bodhaine LUT. It also corrects the
        optical thickness value according to the atmospheric pressure.

        wl      Wavelength [nm]
        pres    Atmospheric pressure [hPa]
    """
    # We read the Bodhaine LUT
    
    # Atmospheric pressure correction tau(p) = tau(p0) p/p0
    taur_wl = [412.0, 443.0, 490.0, 510.0, 555.0, 620.0, 665.0, 680.0, 
              710.0, 750.0, 765.0, 865.0, 1044.0, 1240.0, 1640.0]
    taur_aot = [0.317976, 0.23546, 0.15546, 0.13194, 0.0933786, 
                0.0594837, 0.0447581, 0.0408885, 0.0343288, 
                0.0275019,0.0253863, 0.0154609, 0.00724968, 
                0.00363123, 0.00118286]
    return np.round(np.interp(wl, taur_wl, taur_aot)*pres/1013.25, 8)
###################################################################
# END IMPORTS AND LUTS
###################################################################

###################################################################
# AUTHOR INFORMATION
###################################################################
__author__ = "Francisco Nemiña"
__copyright__ = "Copyright 2021, CONAE"
__credits__ = ["Francisco Nemiña"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Francisco Nemiña"
__email__ = "fnemina@conae.gov.ar"
__status__ = "Test"
###################################################################
# END AUTHOR INFORMATION
###################################################################


###################################################################
# SIMULATION
###################################################################
wl0 = np.array(wl0)

if not hasattr(view, "__len__"):
    view = np.zeros(wl0.size)+view

if not hasattr(phi, "__len__"):
    phi = np.zeros(wl0.size)+phi

# Auto calculate ays

if ays:
    ays = np.round(0.07*sed**0.8,2)

rho = wl0*0-1
for idx, wl in tqdm(np.ndenumerate(wl0), total=wl0.size, desc="Band "):

    s = pyOSOAA.OSOAA(cleanup=cleanup, logfile="/tmp/osoaa_log.tmp")

    # Wavelength
    s.wa = wl/1e3
    # We set the mot
    s.ap.SetMot(taur(wl, pres))
    # We configure the ocean
    s.phyto.chl = chl
    s.phyto.SetPrimaryMode()
    s.sea.depth=300
    s.sea.botalb=0
    s.sed.csed=sed
    s.sed.SetPrimaryMode()
    s.ys.abs440 = ays
    s.ys.swa = sys
    s.det.abs440 = adet
    s.det.swa = sdet
    # Ocean surface
    s.sea.wind = wind
    s.sea.surfalb = np.round(rho_wc(wl, wind), 8)
    # Maritime Shettle and Fenn model
    s.aer.SetModel(model=2, sfmodel=sfmodel, rh=rh)
    # Optical thickness
    s.aer.aotref = aot
    s.aer.waref = 865/1e3
    # We meassure the TOA reflectance
    s.view.level = level

    # Angular configuration 
    s.view.phi = phi[idx]
    s.ang.thetas = sun

    # run simulation and interpolate rho
    s.run()

    rho[idx] = np.interp(view[idx],s.outputs.vsvza.vza, s.outputs.vsvza.I)/np.cos(np.pi*sun/180)

rho = np.round(rho, 6)
###################################################################
# END SIMULATION
###################################################################

###################################################################
# REPORT
###################################################################
if not hasattr(sun, "__len__"):
    sun = np.zeros(wl0.size)+sun

header = "###################################################################\n"
header = header + "SABIA-mar simulations generator \n"
header = header + "###################################################################\n"
header = header + "\n"
header = header + f" Atmosphere configurations\n"
header = header + f" Atmosphere - Aerosol model                 - {models[sfmodel]}\n"
header = header + f" Atmosphere - AOT(865nm)                    - {aot}\n"
header = header + f" Atmosphere - Relative humidity [hPa]       - {rh}\n"
header = header + f" Atmosphere - Atmospheric pressure [%]      - {pres}\n"
header = header + "\n"
header = header + f" Surface - Wind speed [m/s]                 - {wind}\n"
header = header + "\n"
header = header + f" Ocean - Chlorophyll concentration [mg/m^3] - {chl}\n"
header = header + f" Ocean - Sediments [mg/l]                   - {sed}\n"
header = header + f" Ocean - CDOM absorption [1/m]              - {ays}\n"
header = header + f" Ocean - CDOM s [1/nm]                      - {sys}\n"
header = header + f" Ocean - Detritus absorption [1/m]          - {adet}\n"
header = header + f" Ocean - Detritus s [1/nm]                  - {sdet}\n"
header = header + "\n"
header = header + "###################################################################\n"
header = header + f"wl,      rho, view,   sun, phi"

array = [wl0, rho,  view,  sun, phi]

footer = "###################################################################\n"
footer = footer + "Copy here to paste in FORTRAN code\n"
footer = footer + f"wl:   {list(array[0])}\n"
footer = footer + f"rho:  {list(array[1])}\n"
footer = footer + f"view: {list(array[2])}\n"
footer = footer + f"sun:  {list(array[3])}\n"
footer = footer + f"phi:  {list(array[4])}\n"
footer = footer + "###################################################################"

array = np.transpose(array)
np.savetxt(filename, array, delimiter=",", fmt = "%4d, %0.6f, %0.2f, %0.2f, %0.2f",header=header, footer=footer)

###################################################################
# END REPORT
###################################################################

###################################################################
# PLOT
###################################################################
if plotname != None:
    plt.plot(wl0, rho, label=r"$\rho_{toa}$")

    # Plot configuration
    plt.xlabel(r"$\lambda$[nm]")
    plt.ylabel(r"$\rho$")

    title = r"$\theta_0=$"+f"{list(sun)}"+r"°"+"\n"+\
            r"$\theta_v=$"+f"{list(view)}"+r"°"+"\n"+\
            r"$\Delta\phi=$"+f"{list(phi)}"+r"°"+"\n"+\
            r"W$=$"+f"{wind}"+r"m/s"+"\t"+\
            r"[chla]$=$"+f"{chl}"+r"mg/m3 - "+"\t"+\
            r"[sed]$=$"+f"{sed}"+r"mg/l"+"\n"+\
            r"ays$=$"+f"{ays}"+r"1/n - "+"\t"+\
            r"sys$=$"+f"{sys}"+r"1/nm - "+"\t"+\
            r"adet$=$"+f"{adet}"+r"1/n - "+"\t"+\
            r"sdet$=$"+f"{sdet}"+r"1/nm"+"\n"+\
            r"P$=$"+f"{pres}"+r"hPa"+"\t"+\
            fr"{models[sfmodel]}"+"\t"+\
            r"RH$=$"+f"{rh}"+r"%"+"\t"+\
            r"$\tau_a=$"+f"{aot}"+r""

    plt.title(title, fontsize=8)
    plt.legend(loc='upper right')
    plt.ylim(bottom=0)
    plt.tight_layout()

    # Save figure
    plt.savefig(plotname+".png", dpi=300)
    plt.savefig(plotname+".pdf")
print("Finished")
###################################################################
# END PLOT
###################################################################