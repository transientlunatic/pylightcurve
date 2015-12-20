import numpy as np
from math import factorial
import scipy.signal as signal

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
      Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
      The Savitzky-Golay filter removes high frequency noise from data.
      It has the advantage of preserving the original shape and
      features of the signal better than other types of filtering
      approaches, such as moving averages techniques. This implementation is
      taken from [3]_.

      Parameters
      ----------
      y : array_like, shape (N,)
          The values of the time history of the signal.
      window_size : int
          The length of the window. Must be an odd integer number.
      order : int
          The order of the polynomial used in the filtering.
          Must be less then `window_size` - 1.
      deriv: int, default: 0
          the order of the derivative to compute (default = 0 means only smoothing)

      Returns
      -------
      ys : :class:`numpy.ndarray`, shape (N)
          the smoothed signal (or it's n-th derivative).

      Notes
      -----
      The Savitzky-Golay is a type of low-pass filter, particularly
      suited for smoothing noisy data. The main idea behind this
      approach is to make for each point a least-square fit with a
      polynomial of high order over a odd-sized window centered at
      the point.

      Examples
      --------
      >>> t = np.linspace(-4, 4, 500)
      >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
      >>> ysg = savitzky_golay(y, window_size=31, order=4)
      >>> import matplotlib.pyplot as plt
      >>> plt.plot(t, y, label='Noisy signal')
      >>> plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
      >>> plt.plot(t, ysg, 'r', label='Filtered signal')
      >>> plt.legend()
      >>> plt.show()

      References
      ----------
      .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
        Data by Simplified Least Squares Procedures. Analytical
        Chemistry, 1964, 36 (8), pp 1627-1639.
      .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Cambridge University Press ISBN-13: 9780521880688
      .. [3] http://wiki.scipy.org/Cookbook/SavitzkyGolay
    """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def running_median(y, window):
    """
    A method to subtract a running median for smoothing data.

    Parameters
    ----------
    y : :class:`numpy.ndarray`
       A 1D array containing the data time series.
    window : int
       The number of time bins to use for the running median. At edges
       the window will be shorter.
       
    Returns
    -------
    ffit : :class:`numpy.ndarray`
       A 1D array containing the running median of the data time series.
    """
    
    ffit = np.array([])
    idxs = range(len(y))
    halfwin = int(window/2)
    
    for i in range(len(y)):
        v = (idxs < (idxs[i]+halfwin)) & (idxs > (idxs[i]-halfwin))
        ffit = np.append(ffit, np.median(y[v]))

    return ffit

def median_filter(y, window):
    """A method to subtract a median-filtered signal from the original
    signal using the scipy implementaiton of a median filter.

    Parameters
    ----------
    y : :class:`numpy.ndarray`
       A 1-dimensional array containing the data time series.
    window : int
       An (odd) integer defining the length, in samples, of the filter window.

    Returns
    -------
    class:`numpy.ndarray`
       A 1-dimensional array of data with the median removed.
    """
    if not window % 2 == 1:
        window += 1

    # Calculate the median filter
    filt = signal.medfilt(y, window)
        
    return y - filt

def highpass_filter_lightcurve(lightcurve, knee=(1./(0.3*86400.)), **kwargs):
    """
    Detrends a light curve by high-pass filtering it using a third order Butterworth
    filter (:func:`scipy.signal.butter`).

    Parameters
    -----------
    x : :class:`numpy.ndarray`
       An array of time stamps
    z : :class:`numpy.ndarray`
       An array containing the time series data
    knee : float, optional, default: 3.858e-05
       The high-pass filter knee frequency in Hz (default is 3.858e-05 Hz or (1/0.3)/day).
    column : str, optional
       The column of data which the filter should be run over.
       Defaults to the default column if not set.

    Returns
    -------
    z : :class:`numpy.ndarray`
       An array which contains a time series which has been smoothed.

    """

    data = lightcurve.data
    if "column" in kwargs:
        column = kwargs["column"]
        z = np.array(data[column])
    else:
        z = lightcurve.clc
    
    x  = lightcurve.cts
    #z  = lightcurve.clc
    nans, stuff = lightcurve._nan_helper(z)
    z = lightcurve.nan_interp(z)
    dt = lightcurve.dt()
    if dt <= 0:
        raise(NameError("[ERROR] Sample time of 0 detected. Halting."))
    fs = lightcurve.fs()
    highcut = knee/(1./(2.*dt))
    zr = z[::-1]               # Reverse the timeseries to remove phase offset
    zd = np.concatenate((zr, z))
    b, a = signal.butter(3, highcut, btype='highpass')
    y = signal.lfilter(b, a, zd)
    z = y[np.floor(len(y)/2):]
    return z

def medfilt1(x=None,L=None):

    '''
    A simple median filter for 1d numpy arrays.
    Designed to be a replica of the medfilt1 function in MATLAB.

    performs a discrete one-dimensional median filter with window
    length L to input vector x. produces a vector the same size 
    as x. boundaries handled by shrinking L at edges; no data
    outside of x used in producing the median filtered output.
    (upon error or exception, returns None.)

    inputs:
        x, Python 1d list or tuple or Numpy array
        L, median filter window length
    output:
        xout, Numpy 1d array of median filtered result; same size as x
    
    bdj, 5-jun-2009
    '''

    # input checks and adjustments --------------------------------------------
    try:
        N = len(x)
        if N < 2:
            print 'Error: input sequence too short: length =',N
            return None
        elif L < 2:
            print 'Error: input filter window length too short: L =',L
            return None
        elif L > N:
            print 'Error: input filter window length too long: L = %d, len(x) = %d'%(L,N)
            return None
    except:
        print 'Exception: input data must be a sequence'
        return None

    xin = np.array(x)
    if xin.ndim != 1:
        print 'Error: input sequence has to be 1d: ndim =',xin.ndim
        return None
    
    xout = np.zeros(xin.size)

    # ensure L is odd integer so median requires no interpolation
    L = int(L)
    if L%2 == 0: # if even, make odd
        L += 1 
    else: # already odd
        pass 
    Lwing = (L-1)/2

    # body --------------------------------------------------------------------

    for i,xi in enumerate(xin):
  
        # left boundary (Lwing terms)
        if i < Lwing:
            xout[i] = np.median(xin[0:i+Lwing+1]) # (0 to i+Lwing)

        # right boundary (Lwing terms)
        elif i >= N - Lwing:
            xout[i] = np.median(xin[i-Lwing:N]) # (i-Lwing to N-1)
            
        # middle (N - 2*Lwing terms; input vector and filter window overlap completely)
        else:
            xout[i] = np.median(xin[i-Lwing:i+Lwing+1]) # (i-Lwing to i+Lwing)

    return xout
