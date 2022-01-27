#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import os
import datetime
import cartopy.crs as ccrs
import colorcet as cc
from utils import calculate_gridcell_areas, sum_total_emissions, read_tropomi, count_obs_in_mask, plot_field
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def imi_preview():
    '''
    Function to perform preview
    Requires preview simulation to have been run already (to generate HEMCO diags)
    Requires TROPOMI data to have been downloaded already
    '''

    # Open HEMCO diags file and calculate total emissions
    total_prior_emissions = 0

    # Read TROPOMI data and apply lat/lon/QA/albedo filters
    # Count number of observations in region of interest

    # Estimate DOFS from averaging kernel per grid cell

    # Write DOFS, averaging kernel, number of obs, etc. to text file

    # Make plot of prior emissions
    # Make plot of observations
    # Make plot of albedo
    # Save plots as .png


if __name__ == '__main__':
    import sys