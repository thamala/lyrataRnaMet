#Script for inferring the DFE using dadi and fitdadi
#Based on the tutorials at https://dadi.readthedocs.io/en/latest/user-guide/dfe-inference/ and https://github.com/LohmuellerLab/fitdadi/blob/master/manual_examples/manual.pdf
#Tuomas Hämälä, January 2022

import dadi
import dadi.DFE as DFE
from dadi.DFE import *
from scipy.stats import gamma, norm
    
## this function needs to added to the "DemogSelModels.py" file in the dadi installation folder    
#def three_epoch(params, ns, pts):
#	nuB,nuF,TB,TF,gamma = params
#	xx = Numerics.default_grid(pts)
#	phi = PhiManip.phi_1D(xx, gamma=gamma)
#	phi = Integration.one_pop(phi, xx, TB, nuB, gamma=gamma)
#	phi = Integration.one_pop(phi, xx, TF, nuF, gamma=gamma)
#	fs = Spectrum.from_phi(phi, ns, (xx,))
#	return fs
	

## infer demographic parameters and theta0	
	
fs = dadi.Spectrum.from_file('J1_gw_4fold.sfs')
data = fs.fold()
ns = data.sample_sizes
pts_l = [40,50,60]
	
func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.three_epoch)
upper_bound = [1000, 1000, 3, 3]
lower_bound = [0.001, 0.001, 0, 0]

p0 = [0.5, 0.5, 0.05, 0.05]
p0 = dadi.Misc.perturb_params(p0, upper_bound=upper_bound,lower_bound=lower_bound)

demog_params = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=15)
                                   
model = func_ex(demog_params,ns,pts_l)
theta0 = dadi.Inference.optimal_sfs_scaling(model, data)

## infer DFE

fs = dadi.Spectrum.from_file('J1_gw_0fold.sfs')
data = fs.fold()
ns = data.sample_sizes
theta_ns = theta0 * 2.76

pts_l = [600,800,1000]
spectra = DFE.Cache1D(demog_params, ns, DFE.DemogSelModels.three_epoch, pts_l=pts_l, gamma_bounds=(1e-5, 500), gamma_pts=300, verbose=True, mp=True)

lower_bound = [0.001, 0]
upper_bound = [1, 1e6]

p0 = [0.1, 1000]
p0 = dadi.Misc.perturb_params(p0, lower_bound=lower_bound, upper_bound=upper_bound)

popt = dadi.Inference.optimize_log(p0, data, spectra.integrate, pts=None,
                                   func_args=[DFE.PDFs.gamma, theta_ns],
                                   lower_bound=lower_bound, upper_bound=upper_bound, 
                                   verbose=len(sel_params), maxiter=10, multinom=False)

# discretized DFE

cdf1 = gamma.cdf(1, popt[0], scale=popt[1])
cdf2 = gamma.cdf(10, popt[0], scale=popt[1]) - cdf1
cdf3 = 1 - gamma.cdf(10, popt[0], scale=popt[1])

# Note that this just an example. Both the demographic inference and fitting the DFE should be repeated multiple times, starting from different initial values and potentially using different ranges and grid points.
# More info here: https://github.com/LohmuellerLab/fitdadi/blob/master/manual_examples/manual.pdf
