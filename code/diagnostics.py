import numpy as np
from matplotlib import pyplot as plt, cm as cm, patheffects as pe, colors as colors, colorbar as cbar
import h5py as h5
import matplotlib.lines as mlines
import emcee

ncounts = h5.File('../code/ncounts.hdf5','r')

modeltypes = ncounts['meta'].attrs['modeltypes']
zs = np.array([z.decode('ASCII') for z in ncounts['meta'].attrs['zs']])
z_float = [1.0e-5,1.0e-4,4.0e-4,0.001,0.002,0.004,0.006,0.008,0.01,0.014,0.020,0.030,0.040]
subtypes = ncounts['meta'].attrs['subtypes']
Lcuts = ncounts['meta'].attrs['Lcuts']
models = ncounts['meta'].attrs['models']
f_bins = np.linspace(0,1,11)
f_rots = np.linspace(0,1,11)
logages = ncounts['logtime'].value
i_younger_than_100Myr = np.where(logages <= 8)
ts = logages[i_younger_than_100Myr]
dts = np.array([np.power(10.0,6.05)] + [np.power(10.0,6.15 + 0.1*i)-np.power(10.0,6.05 + 0.1*i) for i in range(1,51)])

bcmap = cm.get_cmap('plasma')
rcmap = cm.get_cmap('viridis')
tcmap = cm.get_cmap('bone')
zcmap = cm.get_cmap('pink')

parname_dict = {'f_bin':f_bins,'z':zs,'logtime':ts,'f_rot':f_rots}
cmap_dict = {'f_bin':bcmap,'z':zcmap,'logtime':tcmap,'f_rot':rcmap}
scale_dict = {'WC/WN':'linear','WR/RSG':'log','BSG/RSG':'log','WR/O':'log','WR/YSG':'log'}

def get_arrs(subtype, z, Lcut = 0.0, models = 'BPASS', SFH = 'burst'):
    """
    Given a subtype of star, gets the appropriate summed arrays from the data tables for both
    single and binary populations
    
    Parameters
    ----------
    subtype : str
        Subtype of star. Must be in subtypes
    z : str
        Metallicity, in format zXXX. Supports: zem5,zem4,z0004,z001,z002,z004,z006,z008,z010,
        z014,z020,z030,z040. Depends on value of the models parameter.
    Lcut : float
        Minimum luminosity. Must be 0.0 or between 3.0 and 5.0, in steps of 0.1 dex
    models : str
        Choice of evolutionary code. Supports 'BPASS' or 'Geneva'
    SFH : str, or array-like
        Select a type of star forming history. Supported values: 'burst', 'const' or
        array-like with size 51, corresponding to SFR at each log time bin ago.
        Default: 'burst'
    
    Returns
    -------
    b_ncounts : `~numpy.ndarray`
        Sum of all of the subsubtypes that go into the desired subtype for binaries/rotating stars
    s_ncounts : `~numpy.ndarray`
        Sum of all of the subsubtypes that go into the desired subtype for singles/nonrotating 
        stars
    """
    
    assert (type(SFH) == str)|(hasattr(SFH,'__len__')), "Please supply a string or array"
    if type(SFH) == str:
        assert SFH in ['burst','const'], "Only supported values for SFH are 'burst', 'const', or array of SFRs"
    else:
        assert len(SFH) == 51, "Please supply an array-like object of length 51"
        SFH = np.array(SFH)
        
    assert models in ['BPASS','Geneva'], "Only supported values for models are 'BPASS' or 'Geneva'"
    
    if models == 'BPASS':
        assert z != 'z0004', "That metallicity is only valid for models='Geneva'"
    else:
        assert z in ['z0004','z014','z002'], "That metallicity is only valid for models='BPASS'"
    
    if models == 'BPASS':
        b_arr = ncounts['{0}/{1}/{2}/{3}/{4}/ncounts'.format(models,'bin',z,subtype,str(Lcut))]
        s_arr = ncounts['{0}/{1}/{2}/{3}/{4}/ncounts'.format(models,'sin',z,subtype,str(Lcut))]
        
    else:
        b_arr = ncounts['{0}/{1}/{2}/{3}/{4}/ncounts'.format(models,'rot',z,subtype,str(Lcut))]
        s_arr = ncounts['{0}/{1}/{2}/{3}/{4}/ncounts'.format(models,'not',z,subtype,str(Lcut))]
    
    
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

def get_ratio_at_parameter(ratio, z, logtime, f_bin=None, f_rot=None, Lcut1=0.0, Lcut2=0.0, Lcut=None, SFH='burst'):
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
        Binary fraction, must be in range [0,1]. If not given, must provide f_rot.
    f_rot : float
        Fraction of rotating stars, must be in range [0,1]. If not given, must provide f_bin.
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
        The ratio at the given binary fraction/rotating fraction, and metallicity, and time
    """
    if Lcut is not None:
        Lcut1 = np.clip(Lcut1,a_min=Lcut,a_max=None)
        Lcut2 = np.clip(Lcut2,a_min=Lcut,a_max=None)
        
    assert not ((f_rot is not None)&(f_bin is not None)), "Only specify f_bin or f_rot, not both."
    
    assert (f_rot is not None)|(f_bin is not None), "Please specify f_bin or f_rot."
    
    if f_bin is not None:
        models = 'BPASS'
        f_mix = f_bin
    
    else:
        models = 'Geneva'
        f_mix = f_rot
    
    subtypes = ratio.split('/')
    
    subtype1_b,subtype1_s = get_arrs(subtypes[0],z,Lcut1,models=models,SFH=SFH)
    subtype2_b,subtype2_s = get_arrs(subtypes[1],z,Lcut2,models=models,SFH=SFH)
    
    subtype1 = f_mix*subtype1_b + (1.0-f_mix)*subtype1_s
    subtype2 = f_mix*subtype2_b + (1.0-f_mix)*subtype2_s
    
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

def plot_ratios(ratio1,ratio2,par3,par3val,models='BPASS',constraint_dict=None,SFH='burst',fig=None):
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
        Parameter you want frozen. Must be one of 'logtime', 'f_bin', 'f_rot', or 'z'
    par3val : float or str
        Value of par3 to freeze at. If par3 = 'z', must be a BPASS metallicity string. Otherwise
        a float for the binary/rotating fraction (between 0 and 1), or log time (between 6 and 11)
    models : str
        Model set to use if par3 is not 'f_bin' or 'f_rot'. Determines whether the other two 
        parameters to calculate the grid on is the binary fraction or the rotating fraction.
    constraint_dict : dict
        Constraints to place on the grid. Keys are the same as values for par3, values are tuples 
        of min/max parameter values. Example: constraint_dict = {'z':('z002','z014'),
        'logtime':(6,8)} will restrict the grid to being calculated for metallicities between 
        0.002 and 0.014, and ages between 10^6 and 10^8 years (this assumes par3='f_bin'). 
        Can also specify 'Lcut' which specifies a lower luminosity bound for all four species, 
        or 'Lcuts', which is a tuple of length 4. If ratio1=X/Y, ratio2=A/B, 
        constraint_dict = {'Lcuts':(4.9,0.0,3.5,4.0)} applies a minimum log luminosity of 4.9 to 
        X, 0.0 to Y, 3.5 to A, and 4.0 to B.
    SFH : str, or array-like
        Select a type of star forming history. Supported values: 'burst', 'const' or
        array-like with size 51, corresponding to SFR at each log time bin ago.
        Default: 'burst'
    fig : `matplotlib.Figure`
        If given, adds the axes into fig.
        
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
            
        if 'f_rot' in constraint_dict:
            frot_min,frot_max = constraint_dict['f_rot']
            frots_good = f_rots[(f_rots >= frot_min) & (f_rots <= frot_max)]
            
        else:
            frots_good = f_rots
            
        if 'z' in constraint_dict:
            z_min,z_max = constraint_dict['z']
            zs_val = np.array([z_to_val(z_t) for z_t in zs])
            zs_good = zs[(zs_val >= z_to_val(z_min)) & (zs_val <= z_to_val(z_max))]
            
        else:
            if par3 == 'f_bin':
                zs_good = [z for z in zs if z != 'z0004']
            elif par3 == 'f_rot':
                zs_good = ['z0004','z002','z014']
            elif models == 'BPASS':
                zs_good = [z for z in zs if z != 'z0004']
            elif models == 'Geneva':
                zs_good = ['z0004','z002','z014']
            else:
                assert False, 'You broke something, didnt you?'
            
            
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
        frots_good = f_rots
        if par3 == 'f_bin':
            zs_good = [z for z in zs if z != 'z0004']
        elif par3 == 'f_rot':
            zs_good = ['z0004','z002','z014']
        elif models == 'BPASS':
            zs_good = [z for z in zs if z != 'z0004']
        elif models == 'Geneva':
            zs_good = ['z0004','z002','z014']
        else:
            assert False, 'You broke something, didnt you?'
        
        Lcut11 = 0.0
        Lcut12 = 0.0
        Lcut21 = 0.0
        Lcut22 = 0.0
        Lcut = 0.0
    
    if fig is None:
        fig,ax = plt.subplots(1,2,figsize=(12,6))
    else:
        ax = fig.axes
    
    if par3 == 'logtime':
        if models == 'BPASS':
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
        elif models == 'Geneva':
            for f in frots_good:
                r1 = [get_ratio_at_parameter(ratio1,f_rot=f,logtime=par3val,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
                r2 = [get_ratio_at_parameter(ratio2,f_rot=f,logtime=par3val,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
                ax[0].loglog(r1,r2,c=rcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
                ax[1].axhline(y=f,c=rcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            for z_t in zs_good:
                r1 = [get_ratio_at_parameter(ratio1,f_rot=f,logtime=par3val,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for f in frots_good]
                r2 = [get_ratio_at_parameter(ratio2,f_rot=f,logtime=par3val,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for f in frots_good]
                ax[0].loglog(r1,r2,c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
                ax[1].axvline(x=np.log10(z_to_val(z_t)),c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].set(ylabel=r'$f_{rot}$',xlabel=r'Log $Z$')
        
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
        
    elif par3 == 'f_rot':
        for t in ts_good:
            r1 = [get_ratio_at_parameter(ratio1,f_rot=par3val,logtime=t,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            r2 = [get_ratio_at_parameter(ratio2,f_rot=par3val,logtime=t,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for z_t in zs_good]
            ax[0].loglog(r1,r2,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axvline(t,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
        for z_t in zs_good:
            r1 = [get_ratio_at_parameter(ratio1,f_rot=par3val,logtime=t,z=z_t,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for t in ts_good]
            r2 = [get_ratio_at_parameter(ratio2,f_rot=par3val,logtime=t,z=z_t,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for t in ts_good]
            ax[0].loglog(r1,r2,c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].axhline(y=np.log10(z_to_val(z_t)),c=zcmap(z_to_col(z_t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])  
        ax[1].set(xlabel='Log Time [yr]',ylabel=r'Log $Z$')
        
    elif par3 == 'z':
        if models == 'BPASS':
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
            
        elif models == 'Geneva':
            for f in frots_good:
                r1 = [get_ratio_at_parameter(ratio1,f_rot=f,logtime=t,z=par3val,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for t in ts_good]
                r2 = [get_ratio_at_parameter(ratio2,f_rot=f,logtime=t,z=par3val,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for t in ts_good]
                ax[0].loglog(r1,r2,c=rcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
                ax[1].axhline(y=f,c=rcmap(f),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            for t in ts_good:
                r1 = [get_ratio_at_parameter(ratio1,f_rot=f,logtime=t,z=par3val,Lcut1=Lcut11,Lcut2=Lcut12,Lcut=Lcut,SFH=SFH) for f in frots_good]
                r2 = [get_ratio_at_parameter(ratio2,f_rot=f,logtime=t,z=par3val,Lcut1=Lcut21,Lcut2=Lcut22,Lcut=Lcut,SFH=SFH) for f in frots_good]
                ax[0].loglog(r1,r2,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
                ax[1].axvline(t,c=tcmap(t_to_col(t)),lw=3,path_effects=[pe.Stroke(linewidth=4, foreground='0.75'), pe.Normal()])
            ax[1].set(xlabel='Log Time [yr]',ylabel=r'$f_{rot}$') 
    
    
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

def ratio_lnlikelihood(theta,X,Y):
    """
    Calculates the likelihood of the data assuming the underlying ratio, R, and the true
    value of the number of species Y, lambda_Y, given the data X and Y
    
    Parameters
    ----------
    theta : tuple
        model parameters, R and lambda_Y
    
    X : int
        Numerator of the observed ratio
    Y : int
        Denominator of the observed ratio
        
    Returns
    -------
    lnlike : float
        Natural log of the likelihood function.
    
    """
    
    R,lambda_Y = theta
    lognumerator = X*np.log(R) + (X+Y)*np.log(lambda_Y) - lambda_Y*(R+1.0)
    logdenominator = np.sum(np.log(np.arange(1,X+1))) + np.sum(np.log(np.arange(1,Y+1)))
    
    return lognumerator - logdenominator

def ratio_lnprior_phi_half(theta):
    """
    Calculates the prior probability of the underlying ratio, R, and the true
    value of the number of species Y, lambda_Y, given phi = 1/2. Ensures lambda_Y and R are 
    positive
    
    Parameter
    ---------
    theta : tuple
        model parameters, R and lambda_Y
        
    Returns
    -------
    lnprior : float
        Natural log of the prior probability.
    
    """
    
    R,lambda_Y = theta
    if (R <= 0) or (lambda_Y <= 0):
        return -np.inf
    phi = 0.5
    return (phi - 1.0)*np.log(R) + (2.0*phi - 1.0)*np.log(lambda_Y)

def ratio_lnprior_phi_zero(theta):
    """
    Calculates the prior probability of the underlying ratio, R, and the true
    value of the number of species Y, lambda_Y, given phi = 0. Ensures lambda_Y and R are positive
    
    Parameter
    ---------
    theta : tuple
        model parameters, R and lambda_Y
        
    Returns
    -------
    lnprior : float
        Natural log of the prior probability.
    
    """
    
    R,lambda_Y = theta
    if (R <= 0) or (lambda_Y <= 0):
        return -np.inf
    phi = 0.0
    return (phi - 1.0)*np.log(R) + (2.0*phi - 1.0)*np.log(lambda_Y)

def ratio_lnprior_phi_one(theta):
    """
    Calculates the prior probability of the underlying ratio, R, and the true
    value of the number of species Y, lambda_Y, given phi = 1. Ensures lambda_Y and R are positive
    
    Parameter
    ---------
    theta : tuple
        model parameters, R and lambda_Y
        
    Returns
    -------
    lnprior : float
        Natural log of the prior probability.
    
    """
    
    R,lambda_Y = theta
    if (R <= 0) or (lambda_Y <= 0):
        return -np.inf
    phi = 1.0
    return (phi - 1.0)*np.log(R) + (2.0*phi - 1.0)*np.log(lambda_Y)

lnprior_lookup = {0:ratio_lnprior_phi_zero,1/2:ratio_lnprior_phi_half,1:ratio_lnprior_phi_one}

def ratio_lnprob(theta, X, Y, phi):
    """
    Posterior distribution function
    
    Parameters
    ----------
    theta : tuple
        Contains R and lambda_Y
    
    X : int
        Observed number of species X
        
    Y : int
        Observed number of species Y
    
    phi : float
        Slope of prior probability function
        
    Returns
    -------
    lnprob : float
        The log posterior probability
    
    """
    lp = lnprior_lookup[phi](theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + ratio_lnlikelihood(theta, X, Y)


def MCMC_ratio_errors(X,Y,phi=0.5,nwalkers=100,nburnin=500,nsteps=3000,prob_width=68.0):
    """
    Performs a Markov-Chain Monte Carlo simulation to estimate the value and a confidence 
    interval for the ratio X/Y.
    
    Parameters
    ----------
    X : int
        Observed number of species X
    Y : int
        Observed number of species Y
    phi : float
        Slope of prior probability function
    nwalkers : int
        Number of MCMC walkers, default 100
    nburnin : int
        Number of burn-in steps to take that are then discarded, default 500
    nsteps : int
        Number of production steps to take, default 3000
    prob_width : float
        Desired size of the confidence interval, in percentage, default 68%.
        
    Returns
    -------
    R : `numpy.ndarray`
        Contains the median and upper/lower errors (68th percentile) for the true ratio
        
    lambda_Y : `numpy.ndarray`
        Contains the median and upper/lower errors (68th percentile) for the true lambda_Y
        
    sampler : `emcee.EnsembleSampler`
        The finished sampler object, which can be checked to ensure convergence.
    
    """
    
    R_test = np.clip(X/Y,1e-10,1e5) #if X or Y are zero, clips to a very large or small number
    
    pos = [np.array([R_test,Y]) + 1e-4*np.random.randn(2) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, 2, ratio_lnprob, args=(X,Y,phi))
    sampler.run_mcmc(pos, nburnin)
    p1 = sampler.chain[:, -1, :]
    sampler.clear_chain()
    sampler.run_mcmc(p1, nsteps)
    samples = sampler.flatchain
    
    R, lambda_Y = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, 
                                                [50 - prob_width/2, 50, 50 + prob_width/2],
                                                axis=0)))
    return np.array(R), np.array(lambda_Y), sampler