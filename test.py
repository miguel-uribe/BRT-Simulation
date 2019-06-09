# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 09:25:32 2017

@author: miguel
"""

import createSystem
import optimization
import numpy as np
import matplotlib.pylab  as pyl
import busC
import parameters


# The Line IDs
LineIDs=[0]

# Creating the system, the wagons are asigned
s=[1,2,3,3]    # wagon for lines [0 and 1, 2 and 3, 4 and 5, 6 and 7]

[lines,stations,limits]=createSystem.createsystem(s)
xmax=limits[1]
xmin=limits[0]

# The stopping time
stime=15

#The number of buses, real density and distance
[Nbuses,rdensity,dist]=optimization.getfromdensity(0.02,limits)

p=[Nbuses]

# Creting the distribution list
linelist=[]
for j in range(len(p)):
    for i in range(int(p[j])):
        linelist.append(LineIDs[j])
        
Time=0

#############################################
# Introducing the buses    
pos=parameters.Db  
busID=0
try:
    del buses
except:
    pass

#############################################
for line in linelist:
    if pos>xmax:
        print("There is a bus out of bounds at position %d"%pos)
    try:
        buses=np.append(buses,[busC.createbus(line,pos,busID%2,lines,stations,Time,stime,busID)], axis=0)
    except:
        buses=np.array([busC.createbus(line,pos,busID%2,lines,stations,Time,stime,busID)])
    busID=busID+1
    pos=pos+dist+parameters.Db  
    
#Gaps=[busC.gapsl(bus,buses) for bus in buses]

busC.calculategaps(buses,True,limits)

#buses=busC.updatebusorder(buses)

