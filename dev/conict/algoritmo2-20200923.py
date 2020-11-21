# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 20:48:55 2020

@author: breno
"""
import numpy as np
i = np.array([-2.995, 2, -13, 2, 1, -4])
arr = np.array([-3, 2, -13, 2, 1, -4])
def condition(x): return abs(x - i) <= 0.005 #condicao
bool_arr = condition(arr)
output = np.where(bool_arr)[0] #indices que satisfazem condicao
print(output)