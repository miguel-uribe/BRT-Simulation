    # -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 09:53:10 2017

@author: miguelurla
"""
from __future__ import division
import parameters
import stationC
import lineC


def createsystem (s):

    # The maximum size 
    xmin=0
    
    #Creating all stations
    # This system consisis of Nstations stations in a straight line separated by a distance DS in cell units
    stations=[]
    for i in range(parameters.NStations):
        name="S"+"%d"%i
        stations.append(stationC.stationC(name,parameters.gap+i*parameters.DS))
        
    #Creating all lines, three lines eastbound, three lines westbound
    lines=[]
    
    #E1-W1 lines stopping in every station
    stationIDs=[]
    stops=[]
    # This loop creates a line that stops in every station
    for i in range(parameters.NStations):
        stationIDs.append(stations[i].ID)
        stops.append(s[0])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E1",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W1",list(reversed(stationIDs)),stops,stations,-1))
    
    #E3-W3 lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops every three stations
    for i in range(0,parameters.NStations,3):
        stationIDs.append(stations[i].ID)
        stops.append(s[1])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E3",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W3",list(reversed(stationIDs)),stops,stations,-1))
    
    #E5-W5 lines stopping every five stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops every five station
    for i in range(0,parameters.NStations,5):
        stationIDs.append(stations[i].ID)
        stops.append(s[2])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E5",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W5",list(reversed(stationIDs)),stops,stations,-1))
    
    #E9-W9 lines stopping every five stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops every nine station
    for i in range(0,parameters.NStations,9):
        stationIDs.append(stations[i].ID)
        stops.append(s[3])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E9",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W9",list(reversed(stationIDs)),stops,stations,-1))
    
    # ET-WT lines with no stops
    
    stationIDs=[]
    stops=[]
        
    lines.append(lineC.lineC("ET",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("WT",list(reversed(stationIDs)),stops,stations,-1))
    
    #E3a-W3a lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops every three stations
    for i in range(1,parameters.NStations,3):
        stationIDs.append(stations[i].ID)
        stops.append(s[1])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E3",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W3",list(reversed(stationIDs)),stops,stations,-1))  
    
    
    #E3b-W3b lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops every three stations
    for i in range(2,parameters.NStations,3):
        stationIDs.append(stations[i].ID)
        stops.append(s[1])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E3",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W3",list(reversed(stationIDs)),stops,stations,-1))     
    
    #E9a-E9b lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops, in average, every nine station
    ranStations=[2,7,20,25,38]
    for i in ranStations:
        stationIDs.append(stations[i].ID)
        stops.append(s[3])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E9a",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W9a",list(reversed(stationIDs)),stops,stations,-1))

    #E9b-W9b lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops, in average, every nine station
    ranStations=[4,5,22,23,40]
    for i in ranStations:
        stationIDs.append(stations[i].ID)
        stops.append(s[3])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E9b",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W9b",list(reversed(stationIDs)),stops,stations,-1))    
    
     #E9c-W9c lines stopping every three stations
    stationIDs=[]
    stops=[]
    
    # This loop creates a line that stops, in average, every nine station
    ranStations=[3,6,21,24,39]
    for i in ranStations:
        stationIDs.append(stations[i].ID)
        stops.append(s[3])  # This number indicates the wagon where the line stops
        
    lines.append(lineC.lineC("E9c",stationIDs,stops,stations,1))
    lines.append(lineC.lineC("W9c",list(reversed(stationIDs)),stops,stations,-1))    
    
    
    
    # Establish the size of the system
    size=stations[-1].x+(stations[-1].wagons-1)*parameters.Ds+parameters.Dw+parameters.gap
    
    
    limits=[xmin,size]
    print(limits)
    return [lines,stations,limits]