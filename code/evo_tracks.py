import numpy as np, os
from astropy.io import fits,ascii
from astropy.table import Table,vstack,Column
from scipy.ndimage.filters import gaussian_filter1d
from multiprocessing import Pool

#This is a version of evo_tracks.py for this project, just meant to grab some example tracks.
#It's been super messed with, so it only works for a few tracks from BPASS v2.2.
#The original is in Conti/code

bpass_dir = '/Volumes/shoobert/Research/UW/Diagnostics/data/BPASS/BPASSv2.2/example_tracks/'


#A dict to help us access the BPASS tracks
bpass_col_names = ['Step','Age','logR_1','logTe_1','logL_1','M_1','M_He_core_1','M_CO_core_1',
            'M_ONe_core_1','NULL','X','Y','C','N','O','Ne','M_H_star','M_He_star','M_C_star',
            'M_N_star','M_O_star','M_Ne_star','M_Mg_star','M_Si_star','M_Fe_star','E_bind_env',
            'E_bind_star','M_rem_weakSN','M_eject_weakSN','M_rem_SN','M_eject_SN',
            'M_rem_superSN','M_eject_superSN','L_ang_bin','P_bin','loga','M_1_alt','M_2',
            'M_tot','DMW_1','DMW_2','DMA_1','DMA_2','DMR_1','DMR_2','L_ang_dot','logR_2',
            'logTe_2','logL_2','RLOF_2','IMF_sys','IMF_rej','V-I','U','B','V','R','I','J','H',
            'K','u','g','r','i','z','f300w','f336w','f435w','f450w','f555w','f606w','f814w',
            'U_2','B_2','V_2','R_2','I_2','J_2','H_2','K_2','u_2','g_2','r_2','i_2','z_2',
            'f300w_2','f336w_2','f435w_2','f450w_2','f555w_2','f606w_2','f814w_2','logQ',
            'logFUV','logNUV']

bpass_col_dict = {col:i+1 for i,col in enumerate(bpass_col_names)}
bpass_index_dict = {i+1:col for i,col in enumerate(bpass_col_names)}

default_bpass_col_names = ['Age','logR_1','logTe_1','logL_1','M_1','X','Y','C','N','O','P_bin','M_2','DMW_1','DMW_2','DMA_1','DMA_2','DMR_1','DMR_2','logR_2','logTe_2','logL_2','U','B','V','R','I','J','H','K','u','g','r','i','z','f300w','f336w','f435w','f450w','f555w','f606w','f814w','U_2','B_2','V_2','R_2','I_2','J_2','H_2','K_2','u_2','g_2','r_2','i_2','z_2','f300w_2','f336w_2','f435w_2','f450w_2','f555w_2','f606w_2','f814w_2']

bpass_ms = np.array([300,200,150,120,100,80,70,60,50,40,35,30,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9.5,9,8.5,8,7.5,7,6.5,6,5.5,5,4.5,4,3.7,3.5,3.2,3,2.7,2.5,2.3,2.1,2,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1,0.9,0.8,0.6,0.5,0.4,0.3,0.2,0.1])
bpass_logPs = np.arange(0,4.2,0.2)
bpass_qs = np.arange(0.1,1,0.1)

class track:
    """#Object that contains, at minimum, an initial mass, an array of times, current mass
    luminosity, temperature, as well as any other info. Takes an astropy table or 'BPASS'
    as the track keyword. If BPASS tracks are enabled, the bpass parameters are used to 
    search for the right track and read it in. If you want parameters that aren't in the above
    list, give bpass_cols = [LIST OF COLUMN NAMES]. Make sure the names are from the above 
    list of available column names"""
    def __init__(self, mini = 9.0, mini_keyword = 'Mini', tracks = 'BPASS', Z = 'z014', bpass_M = 1.0, bpass_q = 0.5, bpass_logP = 2, bpass_cols = 'default', **kwargs):
        self.mini = mini
        self.kwargs = kwargs
        self.Z = Z
        assert tracks=='BPASS', "This version of evo_tracks is only for demo purposes"
        if bpass_M % 1 == 0:
            m_str = int(bpass_M)
        else:
            m_str = bpass_M
        q_str = bpass_q
        if not isinstance(bpass_logP, str):
            if bpass_logP % 1 == 0:
                P_str = int(bpass_logP)
        else:
            P_str = bpass_logP

        #We're only working with solar metallicity for now
        if not ((self.Z == 'z014')|(self.Z == 'z002')):
            raise NotImplementedError("Unsupported metallicity")

        #Check to see if user wants more columns than the default
        names = default_bpass_col_names
        if bpass_cols != 'default':
            if type(bpass_cols) != list:
                raise TypeError("This argument should be a list of additional columns to include")
            else:
                #For each name in the list, check if its already included, then include it
                for new_name in bpass_cols:
                    if new_name in names:
                        pass
                    else:
                        names.append(new_name)
            
        #Get the indices of the columns for reading in
        index = [bpass_col_dict[col] for col in names]
        if bpass_logP == 'inf':
            file_name = bpass_dir+'sneplot-{0}-{1}'.format(Z,str(m_str))
        else:
            file_name = bpass_dir+'/sneplot-{0}-{1}-{2}-{3}'.format(Z,str(m_str),str(q_str),str(P_str))

        #Check if the track exists in the directory...
        if os.path.isfile(file_name):

            #read in only given columns
            track = ascii.read(file_name,include_names = ['col{}'.format(i) for i in index])

            #rename columns for easy access later
            for i,name in zip(index,names):
                track['col{}'.format(i)].name = name

            #Create columns with info to remember the initial mass, mass ratio, period
            M1_i = Column([bpass_M for i in range(len(track))],name = 'M1_i')
            q_i = Column([bpass_q for i in range(len(track))],name = 'q_i')
            logP_i = Column([bpass_logP for i in range(len(track))],name = 'logP_i')

            track.add_columns([M1_i,q_i,logP_i])

            #Assign object variable with this track
            self.track = track

        else:
            print(file_name)
            raise FileNotFoundError('That combo of mass, q and logP does not exist')
        
    def get_time(self, time, time_keyword = 'Time'):
        #interpolate the values in the track to a given time. If the time is invalid, 
        #just return 0 (if too early) or nan (if it's exploded)
        #This is disgusting and hacky as all hell. 
        #I'm so sorry future me.
        out_arr = []
        out_types = []
        for col in self.track.colnames:
            if (self.track[col].dtype == '<U1')|(self.track[col].dtype == '<U13'): #if the column is a string and we can't interpolate values, just give me the value from the closest time step
                try:
                    close_idx = np.argmin(np.abs(self.track[time_keyword].data - time))
                    val = self.track[col].data[close_idx]
                except ValueError:
                    val = []
                    for t in time:
                        close_idx = np.argmin(np.abs(self.track[time_keyword].data - t))
                        val.append(self.track[col].data[close_idx])
                        
            else:
                val = np.interp(time, self.track[time_keyword].data, self.track[col].data, left = 0.0, right = np.nan)
                
            out_arr.append(np.array(val))
            
            if self.track[col].dtype == '<U13':
                out_types.append(float)
            else:
                out_types.append(self.track[col].dtype)
          
        return Table(np.array(out_arr).T, names = self.track.colnames, dtype = out_types)
                
