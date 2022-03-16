# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 17:24:49 2022

@author: ethan
"""

#sample file to plot neuton flux on the left boundary in the highest energy group.

import numpy as np
import matplotlib.pyplot as plt

data = np.load('phi_Hist.npz')

t = data['t']
p = data['phi_hist']

plt.loglog(t,p[0,:])
plt.title('Flux History')