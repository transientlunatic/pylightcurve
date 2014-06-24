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
    meta = {}
    title = "Light Curve"

    def __init__(self, start, end):
        """
        Assemble the Lightcurve.
        
        """
        pass
    
    def import_data(self, data_list):
        data = data_list[0]
        meta = data_list[1]
        self.data = self.data.join(data, how="outer")
        for column in data.columns.values.tolist():
            self.meta[column] = meta[column]
        
        return self
    
    def header(self):
        for column in self.meta:
            print column
            for field in self.meta[column]:
                print "\t"+str(field)+":\t"+str(self.meta[column][field])
            print "\n"
    
    def plot(self):
        self.data = self.data.fillna(method='pad')
        self.data.neutron_flux.plot(figsize=(18,6))
        self.data.sunspot.plot(secondary_y=True)