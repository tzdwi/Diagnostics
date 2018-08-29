import numpy as np
from matplotlib import pyplot as plt, cm as cm, patheffects as pe, colors as colors, colorbar as cbar
import h5py as h5
import matplotlib.lines as mlines

ncounts = h5.File('../code/ncounts.hdf5','r')

modeltypes = ncounts['meta'].attrs['modeltypes']
zs = np.array([z.decode('ASCII') for z in ncounts['meta'].attrs['zs']])
z_float = [1.0e-5,1.0e-4,0.001,0.002,0.004,0.006,0.008,0.01,0.014,0.020,0.030,0.040]
subtypes = ncounts['meta'].attrs['subtypes']
Lcuts = ncounts['meta'].attrs['Lcuts']
f_bins = np.linspace(0,1,11)
logages = ncounts['logtime'].value
i_younger_than_100Myr = np.where(logages <= 8)
ts = logages[i_younger_than_100Myr]
dts = np.array([np.power(10.0,6.05)] + [np.power(10.0,6.15 + 0.1*i)-np.power(10.0,6.05 + 0.1*i) for i in range(1,51)])

bcmap = cm.get_cmap('plasma')
tcmap = cm.get_cmap('bone')
zcmap = cm.get_cmap('pink')

parname_dict = {'f_bin':f_bins,'z':zs,'logtime':ts}
cmap_dict = {'f_bin':bcmap,'z':zcmap,'logtime':tcmap}
scale_dict = {'WC/WN':'linear','WR/RSG':'log','BSG/RSG':'log','WR/O':'log','WR/YSG':'log'}

def get_arrs(subtype, z, Lcut = 0.0, SFH = 'burst'):
    """
    Given a subtype of star, gets the appropriate summed arrays from the data tables for both
    single and binary populations
    
    Parameters
    ----------
    subtype : str
        Subtype of star. Must be in subtypes
    z : str
        Metallicity, in format zXXX. Supports: zem5,zem4,z001,z002,z004,z006,z008,z010,
        z014,z020,z030,z040
    Lcut : float
        Minimum luminosity. Must be 0.0 or between 3.0 and 5.0, in steps of 0.1 dex
    SFH : str, or array-like
        Select a type of star forming history. Supported values: 'burst', 'const' or
        array-like with size 51, corresponding to SFR at each log time bin ago.
        Default: 'burst'
    
    Returns
    -------
    b_ncounts : `~numpy.ndarray`
        Sum of all of the subsubtypes that go into the desired subtype for binaries
    s_ncounts : `~numpy.ndarray`
        Sum of all of the subsubtypes that go into the desired subtype for singles
    """
    
    assert (type(SFH) == str)|(hasattr(SFH,'__len__')), "Please supply a string or array"
    if type(SFH) == str:
        assert SFH in ['burst','const'], "Only supported values for SFH are 'burst', 'const', or array of SFRs"
    else:
        assert len(SFH) == 51, "Please supply an array-like object of length 51"
        SFH = np.array(SFH)
    
    b_arr = ncounts['bin/{0}/{1}/{2}/ncounts'.format(z,subtype,str(Lcut))].value
    s_arr = ncounts['sin/{0}/{1}/{2}/ncounts'.format(z,subtype,str(Lcut))].value
    
    
    if type(SFH) ==  str:
        if SFH == 'burst':
            return b_arr,s_arr
    
        elif SFH == 'const':
            #multiply by width of time bin
            product_b = b_arr * dts
            product_s = s_arr * dts
        
    else:
        #multiply by width of time bin
        product_b = b_arr * dts * SFH
        product_s = s_arr * dts * SFH
        
    #sum up to current time bin (cumulative sum). If we're here then SFH != 'burst'
    sum_b = np.cumsum(product_b)
    sum_s = np.cumsum(product_s)

    return sum_b,sum_s

def get_ratio_at_parameter(ratio,z,logtime,f_bin,Lcut1 = 0.0,Lcut2 = 0.0,Lcut = None, SFH = 'burst'):
    """
    Given a ratio, gets the appropriate summed arrays from the data tables for both
    single and binary populations, interpolates to the binary fraction, and then the 
    time
    
    Parameters
    ----------
    ratio : str
        Name of ratio in the format 'X/Y'. You best hope 'X' and 'Y' are in subtypes.
    z : str
        Metallicity, in format zXXX. Supports: zem5,zem4,z001,z002,z004,z006,z008,z010,
        z014,z020,z030,z040
    logtime : float
        Log of time in years. Will interpolate.
    f_bin : float
        Binary fraction, must be in range [0,1]
    Lcut1 : float
        Minimum luminosity for species 1. Must be 0.0 or between 3.0 and 5.0, in steps 
        of 0.1 dex
    Lcut2 : float
        Minimum luminosity for species 2. Must be 0.0 or between 3.0 and 5.0, in steps 
        of 0.1 dex
    Lcut : float
        If given, sets a hard floor for both species
    SFH : str, or array-like
        Select a type of star forming history. Supported values: 'burst', 'const' or
        array-like with size 51, corresponding to SFR at each log time bin ago.
        Default: 'burst'
    
    Returns
    -------
    result : float
        The ratio at the given binary fraction, and metallicity, and time
    """
    if Lcut is not None:
        Lcut1 = np.clip(Lcut1,a_min=Lcut,a_max=None)
        Lcut2 = np.clip(Lcut2,a_min=Lcut,a_max=None)
    
    subtypes = ratio.split('/')
    
    subtype1_b,subtype1_s = get_arrs(subtypes[0],z,Lcut1,SFH=SFH)
    subtype2_b,subtype2_s = get_arrs(subtypes[1],z,Lcut2,SFH=SFH)
    
    subtype1 = f_bin*subtype1_b + (1.0-f_bin)*subtype1_s
    subtype2 = f_bin*subtype2_b + (1.0-f_bin)*subtype2_s
    
    subtype1_t = np.interp(logtime,logages,subtype1,left=0,right=0)
    subtype2_t = np.interp(logtime,logages,subtype2,left=0,right=0)
    
    return np.divide(subtype1_t,subtype2_t)

def z_to_col(z):
    """
    Given a BPASS metallicity string, returns a value to input into the Z colormap
    
    Parameter
    ---------
    z : str
        BPASS metallicity string
        
    Returns
    -------
    out : int
        Input to zcmap to get the right color out
    """
    pos_in_arr = np.where(np.array(zs) == z)[0]/len(zs)
    return pos_in_arr[0]

def z_to_val(z):
    """
    Given a BPASS metallicity string, returns a corresponding float value
    
    Parameter
    ---------
    z : str
        BPASS metallicity string
        
    Returns
    -------
    z_f : float
        Metallicity as a float
    """
    pos_in_arr = np.where(np.array(zs) == z)[0][0]
    z_f = z_float[pos_in_arr]
    return z_f

def t_to_col(t):
    """
    Given a BPASS time bin, returns a value to input into the time colormap
    
    Parameter
    ---------
    t : float
        BPASS log time bin between 6 and 11
        
    Returns
    -------
    out : int
        Input to tcmap to get the right color out
    """
    pos_in_arr = np.where(np.array(ts) == t)[0]/len(ts)
    return pos_in_arr[0]

def plot_ratios(ratio1,ratio2,par3,par3val,constraint_dict=None,SFH='burst'):
    """
    Plots ratio1 vs. ratio2. Between logtime, metallicity, and f_bin, choose one to freeze,
    and the ratios will be calculated on a grid of the other two options.
    
    Parameters
    ----------
    ratio1 : str
        Ratio you want to plot on the abscissa, in the form X/Y
    ratio2 : str
        Ratio you want to plot on the ordinate, in the form X/Y
    par3 : str
        Parameter you want frozen. Must be one of 'logtime', 'f_bin', or 'z'
    par3val : float or str
        Value of par3 to freeze at. If par3 = 'z', must be a BPASS metallicity string. Otherwise
        a float for the binary fraction (between 0 and 1), or log time (between 6 and 11)
    constraint_dict : dict
        Constraints to place on the grid. Keys are the same as values for par3, values are tuples of
        min/max parameter values. Example: constraint_dict = {'z':('z002','z014'),'logtime':(6,8)}
        will restrict the grid to being calculated for metallicities between 0.002 and 0.014, and 
        ages between 10^6 and 10^8 years (this assumes par3='f_bin'). Can also specify 'Lcut' 
        which specifies a lower luminosity bound for all four species, or 'Lcuts', which is a
        tuple of length 4.If ratio1=X/Y, ratio2=A/B, constraint_dict = {'Lcuts':(4.9,0.0,3.5,4.0)} 
        applies a minimum log luminosity of 4.9 to X, 0.0 to Y, 3.5 to A, and 4.0 to B.
    SFH : str, or array-like
        Select a type of star forming history. Supported values: 'burst', 'const' or
        array-like with size 51, corresponding to SFR at each log time bin ago.
        Default: 'burst'
        
    Returns
    -------
    fig : `~matplotlib.figure`
        Figure object containing the Axes in ax.
    ax : list
        Contains two `~matplotlib.axes.axes` objects. The first has the ratios, the second
        is a reference grid.
    """
    
    if constraint_dict is not None:
        
        if 'logtime' in constraint_dict:
            logtime_min,logtime_max = constraint_dict['logtime']
            ts_good = ts[(ts >= logtime_min) & (ts <= logtime_max)]
            
        else:
            ts_good = ts
            
        if 'f_bin' in constraint_dict:
            fbin_min,fbin_max = constraint_dict['f_bin']
            fbins_good = f_bins[(f_bins >= fbin_min) & (f_bins <= fbin_max)]
            
        else:
            fbins_good = f_bins
            
        if 'z' in constraint_dict:
            z_min,z_max = constraint_dict['z']
            zs_val = np.array([z_to_val(z_t) for z_t in zs])
            zs_good = zs[(zs_val >= z_to_val(z_min)) & (zs_val <= z_to_val(z_max))]
            
        else:
            zs_good = zs
            
        if 'Lcuts' in constraint_dict:
            Lcuts = constraint_dict['Lcuts']
            Lcut11 = Lcuts[0] #numerator of ratio1
            Lcut12 = Lcuts[1] #denominator of ratio1
            Lcut21 = Lcuts[2] #numerator of ratio2
            Lcut22 = Lcuts[3] #denominator of ratio2
        
        else:
            Lcut11 = 0.0
            Lcut12 = 0.0
            Lcut21 = 0.0
            Lcut22 = 0.0
        
        if 'Lcut' in constraint_dict:
            Lcut = constraint_dict['Lcut']
        
        else:
            Lcut = 0.0
            
    else:
        
        ts_good = ts
        fbins_good = f_bins
        zs_good = zs
        
        Lcut11 = 0.0
        Lcut12 = 0.0
        Lcut21 = 0.0
        Lcut22 = 0.0
        Lcut = 0.0
    
    
    fig,ax = plt.subplots(1,2,figsize=(12,6))
    
    if par3 == 'logtime':
        for f in fbins_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=f,logtime=par3val,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=f,logtime=par3val,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            ax[0].loglog(r1,r2,c=bcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axhline(y=f,c=bcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        for z_t in zs_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=f,logtime=par3val,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for f in fbins_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=f,logtime=par3val,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for f in fbins_good]
            ax[0].loglog(r1,r2,c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axvline(x=np.log10(z_to_val(z_t)),c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        ax[1].set(ylabel=r'$f_{bin}$',xlabel=r'Log $Z$') 
        
    elif par3 == 'f_bin':
        for t in ts_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=par3val,logtime=t,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=par3val,logtime=t,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            ax[0].loglog(r1,r2,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axvline(t,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        for z_t in zs_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=par3val,logtime=t,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for t in ts_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=par3val,logtime=t,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for t in ts_good]
            ax[0].loglog(r1,r2,c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axhline(y=np.log10(z_to_val(z_t)),c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])  
        ax[1].set(xlabel='Log Time [yr]',ylabel=r'Log $Z$') 
        
    elif par3 == 'z':
        for f in fbins_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=f,logtime=t,z=par3val,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for t in ts_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=f,logtime=t,z=par3val,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for t in ts_good]
            ax[0].loglog(r1,r2,c=bcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axhline(y=f,c=bcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        for t in ts_good:
            r1 = [get_ratio_at_parameter(ratio1,f_bin=f,logtime=t,z=par3val,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for f in fbins_good]
            r2 = [get_ratio_at_parameter(ratio2,f_bin=f,logtime=t,z=par3val,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for f in fbins_good]
            ax[0].loglog(r1,r2,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axvline(t,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        ax[1].set(xlabel='Log Time [yr]',ylabel=r'$f_{bin}$') 
    
    
    ax[0].set(xlabel=r'${}$'.format(ratio1),ylabel=r'${}$'.format(ratio2),title=par3+' = {0}'.format(par3val))
    
    return fig,ax

def calc_ratio_err(spec1,spec2):
    """
    Handy shorthand to calculates the observed ratio and error spec1/spec2, sigma_spec1/spec2, 
    assumes Poisson noise for both number counts. Note that this is an incorrect assumption.
    
    Parameters
    ----------
    spec1 : int or float
        Observed number of stars of the first species.
    spec2 : int or float
        Observed number of stars of the second species
        
    Returns
    -------
    ratio : float
        The quotient spec1/spec2
    error : float
        The standard error for ratio. Again: relies on some incorrect assumptions, as ratio is 
        not a Poisson variable.
    """
    
    result = spec1/spec2
    
    err = result*np.sqrt((1.0/spec1) + (1.0/spec2))
    
    return result,err
