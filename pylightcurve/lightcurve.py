from __future__ import absolute_import
# -*- coding: utf-8 -*-

__authors__ = ["Daniel Williams"]
__email__ = "mail@daniel-williams.co.uk"

#from __future__ import absolute_import

import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as pl
import os
import urllib2, shutil
import copy
from copy import deepcopy
import matplotlib.mlab as ml
import scipy.signal as signal
import pylightcurve
#from lightcurve import *
from .config import load_config
from .filters import *
from dateutil import parser
import datetime
config = load_config()

class Lightcurve():
    """
    This class is designed to handle light curves; time series of
    electromagnetic data, frequently encountered in astronomical
    data. The light curve can be supplemented with non-EM data,
    such as neutron flux measurements, to allow a 
    """

    data = pd.DataFrame()
    highlight = pd.DataFrame()
    meta = {}
    tags = {}
    title = "Light Curve"
    _dcoffsets = {}
    detrended = False
    detrend_method=None
    
    def __init__(self, title="Light Curve"):
        """
        Assemble the Lightcurve.
        
        """
        self.title=title

    def __getitem__(self, key):
        curve = copy.deepcopy(self)
        curve.data = curve.data[key]
        curve.cts = curve.time_seconds()
        return curve

    def time_seconds(self):
        """
        Returns the time stamps of the light curve in seconds since the
        beginning of the time series.
        """
        dt = (np.array(self.data.index.tolist()) - self.data.index[0].to_datetime())
        helper = np.vectorize(lambda x: x.total_seconds())
        ts = helper(dt)
        return ts

    def import_tags(self, tags, name):
        """
        This method takes a DataFrame with tags, for example flare locations, and
        adds them to the tag store for the lightcurve.

        Parameters
        ----------
        tags : Pandas DataFrame
           A Pandas DataFrame containing the tags.

        name : The name of the tag category.
        """

        self.tags[name] = tags
        return self

    
    def import_data(self, data, meta, **kwargs):
        self.data = self.data.join(data, how="outer")
        for column in data.columns.values.tolist():
            self.meta[column] = meta[column]
            self._dcoffsets[column] = np.median(data[column])
        if "cts" in kwargs:
            self.cts = kwargs["cts"]
        else:
            self.cts = self.time_seconds()
        return self

    def add_highlight(self, data):
        self.highlight = self.highlight.join(data, how="outer")
        
        return self

    def set_default(self, column):
        self.clc = self.data[column]
        self.default=column
        return self
    
    def header(self):
        print self.title
        print "---------------------------------------------"
        for column in self.meta:
            print column
            for field in self.meta[column]:
                print "\t"+str(field)+":\t"+str(self.meta[column][field])
            print "\n"
            
    def update_title(self, title):
        self.title = title
        
    def peek(self):
        """
        This method produces a plot of the data contained within the lightcurve
        object, with appropriate scaling. This method is suitable for glancing at
        the data, but is unsuitable for production or publication quality output.

        
        """
        #self.data = self.data.fillna(method='pad')
        
        # Create the underlying plot
        fig, ax = pl.subplots(figsize=(18,6))
        #fig.tight_layout()
        axes = {}
        position = 1
        number = 0
        for column in self.data.columns.values.tolist():
            if number==0:
                color = ax._get_lines.color_cycle
                ccolor=next(color)
                axes[column] = ax
                axes[column].set_ylabel(column, color=ccolor)
                axes[column].tick_params(axis='y', colors=ccolor)
                self.data[column].plot(ax=axes[column], color=ccolor)
                axes[column].set_xlabel('Time', color='black')
                axes[column].set_title(self.title, color='black')
                number=1
            else:
                axes[column] = ax.twinx()
                axes[column].yaxis.set_label_position("right")
                
                rspine = axes[column].spines['right']
                axes[column].yaxis.tick_right()
                rspine.set_position(('axes', position))
                axes[column].set_frame_on(True)
                axes[column].patch.set_visible(False)
                
                color = ax._get_lines.color_cycle
                ccolor=next(color)
                self.data[column].plot(ax=axes[column], color=ccolor)
                axes[column].set_ylabel(column, color=ccolor)
                axes[column].tick_params(axis='y', colors=ccolor)
                axes[column].yaxis.grid(False, which='major')
                position+=0.1

        pl.show()
        
    def _nan_buffer(self, data, size):
        """
        This function carries out the required NAN smoothing 
        which is controlled by `gap_smooth()`.
        """
        holes = np.where(np.isnan(data))[0]
        number_holes = len(holes)
        
        holes = np.array_split(holes,np.where(np.diff(holes)!=1)[0]+1)
        for hole in holes:
            start = hole[0]
            end = hole[-1]
            
            # Check we're not at the end or the beginning of the
            # array, and if we are, then adjust the size accordingly
            this_size = size
            if (start - size < 0 ):
                # If we're too near the start to complete the full
                # removal, then cut it down so it only goes to the 
                # beginning.
                this_size = start
            elif (start + size > len(data)):
                # Likewise, if we're going to hit the end of the 
                # array.
                this_size = len(data) - start
            
            if isinstance(data, pd.core.frame.DataFrame):
                data.iloc[start - this_size : start] = np.nan
                data.iloc[end: end + this_size] = np.nan
            else:
                data[start - this_size : start] = np.nan
                data[end : end+this_size] = np.nan
                    
        return data
        
    def gap_smooth(self, size, **kwargs):
        """
        This function finds np.nan values in a time series and 
        replaces the values close to the gap with np.nan values.
        This is designed as a way of coping with large spikes in
        the data which are caused by the instrument turning off
        and on, for eample.
        
        Parameters
        ----------
        size : int
           The size of the 'margin' around the hole which is 
           to be replaced by `NAN` values.
        column : str, optional
            The column which the process should be carried out on.
        inplace : bool, optional
           If set to True the original light curve will be modified; if set to False
           a new, modified light curve will be returned.
        
        Returns
        -------
        data : array-like
           The time series with the gaps smoothed around.
           
        """
        new_object = deepcopy(self)
        data = new_object.data
        if isinstance(data, pd.core.frame.DataFrame):
            # If the supplied data is a pandas DataFrame we'll need to
            # decide if we're working on just one column, or all of them.
            if "column" in kwargs:
                column = kwargs["column"]
                dataw = data[column]
                data[column] = self._nan_buffer(dataw, size)
            else:
                for column in data.columns.values.tolist():
                    dataw = data[column]
                    data[column] = self._nan_buffer(dataw, size)
        else:
            # If it's a numpy array, then we have an easier job
            data = self._nan_buffer(data, size)
            
        if "inplace" in kwargs:
            self.data = data
            return self
        else:
            return new_object

    def _download(self, uris, 
                 err='Unable to download data at specified URL',
                 **kwargs):
        """Attempts to download data at the specified URI"""

        # user specifies a download directory
        if "directory" in kwargs:
            download_dir = os.path.expanduser(kwargs["directory"])
        else:
            download_dir = config.get("downloads", "download_dir")

        # overwrite the existing file if the keyword is present
        if "overwrite" in kwargs:
            overwrite = kwargs["overwrite"]
        else:
            overwrite = False

        if not isinstance(uris, list):
            uris = [uris]

        filepaths = []

        for uri in uris:
            _filename = os.path.basename(uri).split("?")[0]

            # If the file is not already there, download it
            filepath = os.path.join(download_dir, _filename)

            if not(os.path.isfile(filepath)) or (overwrite and 
                                                 os.path.isfile(filepath)):
                try:
                    response = urllib2.urlopen(uri)
                except (urllib2.HTTPError, urllib2.URLError):
                    raise urllib2.URLError(err)
                with open(filepath, 'wb') as fp:
                    shutil.copyfileobj(response, fp)
            else:
                print "Using existing file rather than downloading, use overwrite=True to override."

            filepaths.append(filepath)

        return filepaths

    def psd(self, **kwargs):
        """
        Calculate the one-sided non-windowed power spectrum of the light curve. This uses the
        :func:`matplotlib.mlab.psd` function for computing the power spectrum, with a single
        non-overlapping FFT.

        Parameters
        ----------
        column : str, optional
           The column of the lightcurve which should be analysed.
        reduced: float, optional
           Removes sections of the array greater than its value times the mean
           of the array, in order to try and reduce the effect of very large events.

        Returns
        -------
        sk : array-like
           The Power spectral density of the light curve.
        f  : array-like
           An array of the frequencies.
        """

        data = deepcopy(self.data)
        
        if isinstance(data, pd.core.frame.DataFrame):
            # If the supplied data is a pandas DataFrame we'll need to
            # decide if we're working on just one column, or all of them.
            if "column" in kwargs:
                column = kwargs["column"]

                dataw = np.array(data[column])
                dataw = self.nan_interp(dataw)

                if "reduced" in kwargs and kwargs["reduced"] != None:
                    reduced = kwargs["reduced"]
                    mean = np.mean(dataw)
                    std = np.std(dataw)
                    np.delete(dataw,[dataw>mean+reduced*std])
                
                l = len(dataw)
                sk, f = ml.psd(x=dataw, window=signal.boxcar(l), noverlap=0, NFFT=l, Fs=self.fs(), sides='onesided')
            elif self.default:
                column = self.default
                dataw = np.array(data[column])
                dataw = self.nan_interp(dataw)
                if "reduced" in kwargs and kwargs["reduced"] != None:
                    reduced = kwargs["reduced"]
                    mean = np.mean(dataw)
                    std = np.std(dataw)
                    dataw[dataw>mean+reduced*std]=mean
                l = len(dataw)
                sk, f = ml.psd(x=dataw, window=signal.boxcar(l), noverlap=0, NFFT=l, Fs=self.fs(), sides='onesided')
            else:
                sk = {}
                f = {}
                for column in data.columns.values.tolist():
                    dataw = np.array(data[column])
                    dataw = self.nan_interp(dataw)
                    if "reduced" in kwargs and kwargs["reduced"] != None:
                        reduced = kwargs["reduced"]
                        mean = np.mean(dataw)
                        std = np.std(dataw)
                        dataw[dataw>mean+reduced*std]=mean
                    l = len(dataw)
                    sk[column], f[column] = ml.psd(x=dataw, window=signal.boxcar(l), noverlap=0, NFFT=l, Fs=self.fs(), sides='onesided')

        # return power spectral density and array of frequencies
        del data
        return sk, f

    def time_to_index(self, time):
      """
      Returns the index location which corresponds to a given time.
      """
      
      return self.data.index.searchsorted(parser.parse(time))

    def dt(self):
        """
        Calculate the sample separation for the light curve.

        Returns
        -------
        float
           The separation time for data in the light curve.
        """
        if isinstance(self.data.index[0], pd.tslib.Timestamp):
            return (self.data.index[1] - self.data.index[0] ).seconds
        else:
            return (self.data.index[1] - self.data.index[0])

    def fs(self):
        """
        Calculate the sample frequency of the light curve.

        Returns
        -------
        float
           The sample frequency of the light curve.
        """
        return 1.0 / self.dt()

    def nan_interp(self, data):
        nans, x= self._nan_helper(data)
        data = copy.deepcopy(data)
        data[nans]= np.interp(x(nans), x(~nans), data[~nans])
        return data

    def _nan_helper(self,y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]

    def _dcoffset(self, data=None, **kwargs):
        """
        Method to remove a DC offset from a light curve by subtracting the median value of the
        light curve from all values.
        """
        if data == None:
            data = self.clc
            
        self.dc  = np.median(data)
        data = data - self.dc
        return data

    def set_detrend(self, method='none', nbins=None, order=None, knee=None):
        """
        A method allowing the detrending parameters for the light curve to be changed.

        Parameters
        ----------
        method : string
           The detrending method. Can be 'savitzkygolay', 'runningmedian', or,
           'highpassfilter'.
        nbins : int
           The length of the detrending window, in bins.
        order : int
           The order of the detrending filter.

        See also
        --------
        detrend

        """
        self.detrend_method=method
        self.detrend_length=nbins
        self.detrend_nbins=nbins
        self.detrend_order=order
        self.detrend_knee=knee

    def detrend(self, method='none', nbins=None, order=None, knee=None, **kwargs):
        """
        A method to detrend the light curve using a Savitsky-Golay filter (:func:`.savitzky_golay`),
        a running median filter (:func:`.running_median`), or a high-pass filter
        (:func:`.highpass_filter_lightcurve`).

        Parameters
        ----------
        method : string, default: 'none'
           The detrending method. Either 'savitzkygolay', 'runningmedian', or
           'highpassfilter'.
        nbins : int, default: None
           The number of bins in the Savitsky-Golay, or running median detrend window
        order : int, default: None
           The polynomial order of the Savitsky-Golay detrending fit.
        knee : float, default: None
           The high-pass filter knee frequency (Hz).
        inplace : bool, optional
           Replaces the data in the object with the filtered data. By default this is
           False.

        """
        new_object = deepcopy(self)
        data = new_object.data

        self.set_detrend(method=method, nbins=nbins, order=order, knee=knee)
        self.detrended = True

        if "column" in kwargs:
            column = kwargs["column"]
            dataw = np.array(data[column])
        else:
            dataw = self.clc
        
        if method == 'savitzkygolay':
            if nbins is None or order is None:
                raise ValueError("Number of bins, or polynomial order, for Savitsky-Golay filter not set")
            nans, stuff = self._nan_helper(dataw)
            dataw = self.nan_interp(dataw)
            ffit = savitzky_golay(dataw, nbins, order)
            dataw = (dataw - ffit)
            dataw[nans] = np.nan
            
        elif method == 'runningmedian':
            if nbins is None:
                raise ValueError("Number of bins for running median filter not set")

            nans, stuff = self._nan_helper(dataw)
            dataw = self.nan_interp(dataw)
            ffit = running_median(dataw, nbins)
            dataw = (dataw - ffit)
            dataw[nans] = np.nan
            
        elif method == 'highpassfilter':
            if knee is None:
                raise ValueError("Knee frequency for high-pass filter not set.")

            nans, stuff = self._nan_helper(dataw)
            dataw = self.nan_interp(dataw)
            dataw = highpass_filter_lightcurve(self, knee=knee, **kwargs)
            dataw[nans] = np.nan
            #self.clc = np.copy(dlc.clc)
        else:
            raise ValueError("No detrend method set")

        if "column" in kwargs:
            column = kwargs['column']
            data[column]= dataw
        else:
            new_object.clc = dataw
            
        if "inplace" in kwargs:
            self.clc = dataw
            self.data = data
            return self
        else:
            return new_object.denorm()

    def norm(self, **kwargs):
        """
        Normalises all of the data in the light curve

        Parameters
        ----------
        inplace : bool, optional
           If set to True the original light curve will be modified; if set to False
           a new, modified light curve will be returned.
        """
        new_object = deepcopy(self)
        data = new_object.data
        for column in data.columns.values.tolist():
            self._dcoffsets[column] = np.median(data[column])
            data[column] = data[column] - self._dcoffsets[column]

        if self.default:
            self.dc = self._dcoffsets[self.default]
            
        if "inplace" in kwargs:
            self.data = data
            return self
        else:
            new_object.data = data
            return new_object

    def denorm(self, **kwargs):
        """
        De-normalises all of the data in the light curve

        Parameters
        ----------
        inplace : bool, optional
           If set to True the original light curve will be modified; if set to False
           a new, modified light curve will be returned.
        """
        new_object = deepcopy(self)
        data = new_object.data
        for column in data.columns.values.tolist():
            data[column] = data[column] + self._dcoffsets[column]
        if "inplace" in kwargs:
            self.data = data
            return self
        else:
            new_object.data = data
            return new_object

    def interpolate(self, **kwargs):
        """
        A method for interpolating the light curves, to compensate for NaN values.

        Examples
        --------

           >>> camelot = bf.Lightcurve(curves)
           >>> camelot.interpolate()

        """


        new_object = deepcopy(self)
        data = new_object.data
        if "column" in kwargs:
            column = kwargs["column"]
            z = data[column]
            nans, za= self._nan_helper(z)
            z[nans]= np.interp(za(nans), za(~nans), z[~nans]).astype('float32')
            data[column] = z
        else:
            for column in data.columns.values.tolist():
                z = data[column]
                nans, za= self._nan_helper(z)
                z[nans]= np.interp(za(nans), za(~nans), z[~nans]).astype('float32')
                data[column] = z
        if "inplace" in kwargs:
            self.data = data
            return self
        else:
            new_object.data = data
            return new_object