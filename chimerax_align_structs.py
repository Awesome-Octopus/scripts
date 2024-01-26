#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 17:27:07 2023

@author: andrew
"""

from chimera import runCommand as rc
for s in range(13, 102):
    string = 'align #'+str(s)+':12-36 toAtom #1:12-36'
    runCommand('session', string)
