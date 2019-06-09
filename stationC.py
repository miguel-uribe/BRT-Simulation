# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 06:50:22 2017

@author: miguelurla
"""
from __future__ import division
import parameters
# Importiung the anim parameter
anim=parameters.anim

# If anim, import the vpython module
if anim:
    from vpython import *

class stationC:
    "This class contains all information about the stations of the system"
    n=0     #number of stations
    # This is the constructor
    def __init__(self,name,x):
        self.name=name  # Station name
        self.x=x        # The x position of the station
        self.ID=stationC.n  # The station ID
        stationC.n=stationC.n+1  # The station counter is updated
        self.wagons=parameters.Nw    # By default all stations have three wagons
        self.boxes=None # Memory for the graphical representation
        self.lineIDs=[]   # The ID's of the lines
        self.lineTimes={}  # A dictionary storing all times a line reaches a station
        
    # This function overrides the representation    
    def __repr__(self):
        output="%d: %s; x=%d; %d wagons; lines:"%(self.ID,self.name,self.x, self.wagons)
        for line in self.lineIDs:
            output=output+" %d"%line
        return output
     
    def addline(self,lineID):
        self.lineIDs.append(lineID)
        print("Added line")

# This function updates the wagon number
def updatewagons(Stations,stationIDs,stops):
    for i in range(len(stationIDs)):
        for j in range(len(Stations)):
            if Stations[j].ID==stationIDs[i]:
                if stops[i]>Stations[j].wagons:
                    Stations[j].wagons=stops[i]
                    
                    
# defining the function 
# defining the function 
if anim:
    def setrepr(Stations):

        # Setting the wagon representation
          for station in Stations:   
              if station.boxes is None:
                station.boxes=[]
                station.streets=[]
                station.liasons=[]
                # The street between the wagons
                posst=vector(station.x+0.5*parameters.DS,-14,-1)
                station.streets.append(box(pos=posst, length=parameters.DS, width=1,height=6, color=color.gray(0.8) ))              
                
                for i in range(station.wagons):
                    posv=vector(station.x+0.5*parameters.Dw+i*(parameters.Ds),0,0)
                    station.boxes.append(box(pos=posv, length=parameters.Dw, width=1,height=10, color=color.blue))
                    posst=vector(station.x+parameters.Dw+i*(parameters.Ds),-11,-1)
                    station.streets.append(box(pos=posst, length=2*parameters.Dw, width=1,height=12, color=color.gray(0.8) ))
                    if i+1<station.wagons:
                        posv=vector(station.x+parameters.Dw+0.5*(parameters.Ds-parameters.Dw)+i*(parameters.Ds),0,0)
                        station.liasons.append(box(pos=posv, length=parameters.Dw, width=1,height=6, color=color.blue))
                    if i==0:
                        posst=vector(station.x-0.5*parameters.Dw,-11,-1)
                        station.streets.append(box(pos=posst, length=parameters.Dw, width=1,height=12, color=color.gray(0.8)))
                        
# This function retrieves the complete station length
def getstationlength(Stations,stationID):
    for station in Stations:
        if station.ID==stationID:
            size=(station.wagons-1)*parameters.Ds+parameters.Db
            break
    return size

# This function retrieves the station index from the ID
def getstationindexbyID(Stations,ID):
    for i in range(len(Stations)):
        if Stations[i].ID==ID:
            return i
    print("The given line ID has not been found in the Stations list")

# Get station from ID
def getstationbyID(Stations,ID):
    for station in Stations:
        if station.ID==ID:
            return station
    print("The given line ID has not been found in the Stations list")
    return None
    
# Remove a station from list
def removestation(Stations,ID):
    Stations.remove(getstationbyID(Stations,ID))
    
    

        