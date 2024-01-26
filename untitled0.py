#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 20:00:17 2023

@author: andrew
"""

model_name = 'assdf'
count = 0
while True:
    try:
        f = 5/count
        print('trying')
    except:
        count += 1
        print('that doesnt work')
    else:
        print(f'{count} {f} end')
        break
