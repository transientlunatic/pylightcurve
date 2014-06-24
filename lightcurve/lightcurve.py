# -*- coding: utf-8 -*-

__authors__ = ["Daniel Williams"]
__email__ = mail@daniel-williams.co.uk

import numpy as np
import pandas as import pd

class Lightcurve():
    """
    This class is designed to handle light curves; time series of
    electromagnetic data, frequently encountered in astronomical
    data. The light curve can be supplemented with non-EM data,
    such as neutron flux measurements, to allow a 
    """

    data = pd.DataFrame()
    title = "Light Curve"

    def __init__(self, title="Light Curve"):
        """
        Assemble the Lightcurve.
        
        """
        self.title = title
    
    def import_data(self, data):
        self.data = self.data.join(data, how="outer")
        return self
    
    def plot(self):
        self.data = self.data.fillna(method='pad')
        self.data.neutron_flux.plot(figsize=(18,6))
        self.data.sunspot.plot(secondary_y=True)
        
        
        