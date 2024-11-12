r'''Reader for FITS files (.fits)'''

import os

import numpy as np
from astropy.io import fits
from astropy.time import Time

from blimpy.io import sigproc
from blimpy.io.base_reader import Reader, logger, GIGA


class FITSReader(Reader):
    """ This class handles .fits files.
    """

    def __init__(self, filename,f_start=None, f_stop=None,t_start=None, t_stop=None, load_data=True, max_load=None):
        """ Constructor.

        Args:
            filename (str): filename of blimpy file.
            f_start (float): start frequency, in MHz
            f_stop (float): stop frequency, in MHz
            t_start (int): start time bin
            t_stop (int): stop time bin
            max_load (float): memory limit in gigabytes
        """
        super(FITSReader, self).__init__()

        self.header_keywords_types = sigproc.header_keyword_types

        if filename and os.path.isfile(filename):
            self.filename = filename
            self.load_data = load_data
            self.header = self.read_header()
            self.file_size_bytes = os.path.getsize(self.filename)
            self.n_channels_in_file  = self.header['nchans'] # FITS specific
            self.n_beams_in_file = self.header['nifs'] #Placeholder for future development.
            self.n_pols_in_file = 4 #Placeholder for future development.
            self._n_bytes = int(self.header['nbits'] / 8)  #number of bytes per digit
            self._d_type = self._setup_dtype()
            self.n_ints_in_file = self.header['NAXIS1'] # no need to calculate this for FITS
            self.file_shape = (self.n_ints_in_file,self.n_beams_in_file,self.n_channels_in_file)
            self.f_begin  = self.header['fch1']
            self.f_end  = self.f_begin + self.n_channels_in_file*self.header['foff']
            self.t_begin = 0
            self.t_end = self.n_ints_in_file # shouldn't this use the time increment???

            #Taking care all the frequencies are assigned correctly.
            self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop, init=True)
            #Convert input frequencies into what their corresponding channel number would be.
            self._setup_chans()
            #Update frequencies ranges from channel number.
            self._setup_freqs()

            self.freq_axis = 2
            self.time_axis = 0
            self.beam_axis = 1  # Place holder

#EE ie.
#           spec = np.squeeze(fil_file.data)
            # set start of data, at real length of header  (future development.)
#            self.datastart=self.hdrraw.find('HEADER_END')+len('HEADER_END')+self.startsample*self.channels

            #Applying data size limit to load.
            if max_load is not None and max_load > 0:
                self.max_data_array_size = max_load * GIGA

            if self.file_size_bytes > self.max_data_array_size:
                self.large_file = True
            else:
                self.large_file = False

            if self.load_data:
                if self.large_file:
                    if self.f_start or self.f_stop or self.t_start or self.t_stop:
                        if self.isheavy():
                            self.warn_memory("Selection", self._calc_selection_size())
                            self._init_empty_selection()
                        else:
                            self.read_data()
                    else:
                        self.warn_memory("File", self.file_size_bytes)
                        self._init_empty_selection()
                else:
                    self.read_data()
            else:
                logger.debug("Skipping loading data ...")
                self._init_empty_selection()
        else:
            raise IOError("Need a file to open, please give me one!")


    def read_header(self):
        """ Read FITS header and return a Python dictionary of key:value pairs

        Args:
            filename (str): name of file to open

        Returns:
            Python dict of key:value pairs.

        """
        hdr = fits.open(self.filename)[0].header 
        self.header = dict(hdr) # avoids warnings about FITS standards for blimpy keywords

        # add standard blimpy header keywords
        self.header['nchans'] = self.header['NAXIS2']
        self.header['foff'] = self.header['CDELT2']
        self.header['nifs'] = 1
        self.header['nbits'] = -1*self.header['BITPIX']
        self.header['source_name'] = self.header['NAME']
        self.header['tsamp'] = self.header['CDELT1']
        self.header['fch1'] = self.header['CRVAL2']
        self.header['tstart'] = Time(self.header['OBS-STAR'], format='isot', scale='utc').mjd # convert ISOT --> MJD

        return self.header

    def read_data(self, f_start=None, f_stop=None,t_start=None, t_stop=None):
        """ Read data.
        """

        self._setup_selection_range(f_start=f_start, f_stop=f_stop, t_start=t_start, t_stop=t_stop)

        #check if selection is small enough.
        if self.isheavy():
            self.warn_memory("Selection", self._calc_selection_size())
            self.data = np.array([0], dtype=self._d_type)
            return None

        #Convert input frequencies into what their corresponding channel number would be.
        self._setup_chans()

        #Update frequencies ranges from channel number.
        self._setup_freqs()

        n_chans = self.header['nchans']
        n_chans_selected = self.selection_shape[self.freq_axis]
        n_ifs   = 1

        # Load binary data
        f = fits.open(self.filename)
        fdata = f[0].data
        fdata = fdata[0,:,:].T # extract Stokes I; A --> A.T : (n_ints, n_chans)
        
        # now check to see how many integrations requested
        n_ints = self.t_stop - self.t_start

        #Loading data
        self.data = np.zeros((n_ints, n_ifs, n_chans_selected), dtype=self._d_type)
        
        # map selected:
        for ii in range(n_ints):
            for jj in range(n_ifs):
                dd = fdata[ii,self.chan_start_idx:self.chan_start_idx+n_chans_selected]
                self.data[ii, jj] = dd

        # Give the FD back to the O/S.
        f.close()

    



