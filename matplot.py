# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:56:21 2020

@author: Breno
"""


import matplotlib.pyplot as plt

ax1 = plt.figure(1).add_subplot()
ax2 = plt.figure(2).add_subplot()
ax3 = plt.figure(3).add_subplot()

ax1.set_ylabel('y fig 1')
ax2.set_ylabel('y fig 2')
ax3.set_ylabel('y fig 3')

print(ax1)



