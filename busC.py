# -*- coding: utf-8 -*-
# This is 
"""
This is v2 of this module
Created on Wed Jun 07 06:32:11 2017

@author: miguelurla
"""
from __future__ import division


import lineC
import stationC
import parameters
import random
from operator import attrgetter
import numpy as np
import sys


# bus creation
def createbus(lineID,x,y,Lines,Stations,time,stime,busID):
    npars=22  # the number of parameters
    # A bus will be a numpy array with the following information
    A=-np.ones(npars, dtype=int)
    A[0]=x      # x, x position
    A[1]=0      # v, bus speed
    A[2]=y      # y, lane info, 0: main lane; 1: stoppin lane
    #A[3]=gapf   # gapf, Forward gap in the same lane
    #A[4]=gapb   # gapb, backward gap in the same lane
    #A[5]=gapfl   # gapfl, Forward gap in the opposite lane
    #A[6]=gapbl   # gapbl, backward gap in the opposite lane
    #A[7]=vbef    # vbef, Speed of the car behinf in the oppsite lane
    #A[8]=nextstop # nextstop, position of the next stop
    #A[9]=nextstationID # nextstation iD, The ID of the next station
    #A[10]=nextStationEnd # nextStationEnd, X position of the limit of the stoppimg lane in the next station
    #A[11]=nextStationindex # nextSDtationindex, The index of the next station in the line array
    acc=Lines[lineC.getlineindexbyID(Lines,lineID)].acc 
    A[12]=acc       # a, Acceleration
    A[14]=lineID  # lineID, Id of the line
    A[15]=busID      # ID, Bus ID
    A[16]=stime     # AVStopTime, Average stop time for the bus
    A[17]=1          # Active, 1: bus is Active, 0: Bus is not active
    A[18]=1      # Changing, 1: bus is able to change lanes, 0: bus is not able to change lane
    A[19]=1     # 1: Advancing, Bus is advancing, 0: Bus is stopped at a given station
    A[20]=0     # 1: Measuring, Bus is measuring, 0: Bus is not measuring
    A[21]=time          # time, bus creation time

    findnextstop(A,Lines,Stations)
    return A

def findnextstop(bus,Lines,Stations):
    line=lineC.getlinebyID(Lines,bus[14])
    # scanning all stops in the line
    for i in range(len(line.stops)):
        # Buses to the east
        if bus[12]>0 and line.stopx[i]-1.5*parameters.Db>bus[0]:
            bus[8]=line.stopx[i]
            bus[9]=line.stationIDs[i]
            nextstation=stationC.getstationbyID(Stations,bus[9])
            bus[10]=nextstation.x+(nextstation.wagons-1)*parameters.Ds+parameters.Dw+2*parameters.Db
            bus[11]=i
            break
        # Buses to the west
        elif bus[12]<0 and line.stopx[i]+1.5*parameters.Db<bus[12]:
            bus[8]=line.stopx[i]
            bus[9]=line.stationIDs[i]
            nextstation=stationC.getstationbyID(Stations,bus[9])
            bus[10]=nextstation.x-2*parameters.Db
            bus[11]=i
            break
        
# This function must be called after the bus stops in a given station
def updatestop(bus,Lines,Stations):
    line=lineC.getlinebyID(Lines,bus[14])
    i=bus[11]+1  # updating the index to the next station
    if i>=len(line.stops):     # in case the final station is reached
        bus[8]=-1
        bus[9]=-1
        bus[11]=-1 
    # Otherwise all paraneters are updated
    else:   
        bus[8]=line.stopx[i]
        bus[9]=line.stationIDs[i]
        bus[11]=i

# This function must be called after the bus stops in a given station AND returns to the main lane    
def updatestationend(bus,Stations):
    if bus[8]==-1:     # in case the final station is reached
        bus[10]=-1
    #Otherwise
    else:
        nextstation=stationC.getstationbyID(Stations,bus[9])
        if bus[12]>0:
            bus[10]=nextstation.x+(nextstation.wagons-1)*parameters.Ds+parameters.Dw+2*parameters.Db
        elif bus[12]<0:
            bus[10]=nextstation.x-2*parameters.Db
                

# This function calculates the forward gap on the opposite lane
# Here only the presence of other buses is taken into account, the physical boundaries are not expected to have an impact.
def gapsl(bus,Buses):
    
    #finding the bus in the Buses list
    i=np.where(np.all(Buses==bus,axis=1))[0][0]
    ##############################################
    # Buses moving to the east
    if bus[12]>0:
        [gapf,gapb,vbef]=[1000,-1000,0]  # Values by default
        ##############################################
        #Buses on the main lane
        if bus[2]==0:
            ###########################
            # gap forward
            # We create a view list of buses starting in i and in a range of 5 cars
            BusesLtd=Buses[i:min(i+5,len(Buses))]
            # We find the distance to the first bus in the list that is in the stopping lane
            try:
                gapf=BusesLtd[BusesLtd[:,2]==1][0][0]-bus[0]-parameters.Db
            # If no bus is found the gapf is not modified
            except:
                pass
            ###############################
            # gap backwards
            # We create a view list of buses ending in i-1 and in a range of 5 cars
            BusesLtd=Buses[max(0,i-5):i]
            # We find the distance to the last element of this list
            try:
                gapb=BusesLtd[BusesLtd[:,2]==1][-1][0]-bus[0]+parameters.Db
                vbef=BusesLtd[BusesLtd[:,2]==1][-1][1]
            # If no buses are found in the stopping lane, no changes are done to gapb or bvef
            except:
                pass

            
        ######################################
        #Buses on the stopping lane
        elif bus[2]==1:
            #################################
            # gap forward
            # We create a view list of buses starting in i
            BusesLtd=Buses[i:min(i+5,len(Buses))]
            # We find the distance to the first bus in the list that is in the stopping lane
            try:
                gapf=BusesLtd[BusesLtd[:,2]==0][0][0]-bus[0]-parameters.Db
            # If no bus is found the gapf is not modified
            except:
                pass
            ###############################
            # gap backwards
            # We create a view list of buses ending in i-1
            BusesLtd=Buses[max(0,i-5):i]
            # We find the distance to the last element of this list
            try:
                gapb=BusesLtd[BusesLtd[:,2]==0][-1][0]-bus[0]+parameters.Db
                vbef=BusesLtd[BusesLtd[:,2]==0][-1][1]
            # If no buses are found in the stopping lane, no changes are done to gapb or bvef
            except:
                pass
            
            
    #########################################       
    # Buses moving to the west
    elif bus[12]<0:
        [gapf,gapb,vbef]=[-1000,1000,0]  # Value by default
        #########################################
        #Buses on the main lane
        if bus[2]==0:
            #############################
            # gap forward
            # We create a view list of buses ending in i-1
            BusesLtd=Buses[max(0,i-5):i]
            # We find the distance to the last element of this list
            try:
                gapf=BusesLtd[BusesLtd[:,2]==1][-1][0]-bus[0]+parameters.Db
            # If no bus is found the gapf is not modified
            except:
                pass
            ##############################
            # gap backwards
            # We create a view list of buses starting in i
            BusesLtd=Buses[i:min(i+5,len(Buses))]
            # We find the distance to the first bus in the list that is in the stopping lane
            try:
                gapb=BusesLtd[BusesLtd[:,2]==1][0][0]-bus[0]-parameters.Db
                vbef=BusesLtd[BusesLtd[:,2]==1][0][1]
            # If no bus is found the gapf is not modified
            except:
                pass
        
        ###########################################
        #Buses on the stopping lane
        elif bus[2]==1:
            ################################
            # gap forward
            # We create a view list of buses ending in i-1
            BusesLtd=Buses[max(0,i-5):i]
            # We find the distance to the last element of this list
            try:
                gapf=BusesLtd[BusesLtd[:,2]==0][-1][0]-bus[0]+parameters.Db
            # If no bus is found the gapf is not modified
            except:
                pass
            ################################
            # gap backwards
            # gap backwards
            # We create a view list of buses starting in i
            BusesLtd=Buses[i:min(i+5,len(Buses))]
            # We find the distance to the first bus in the list that is in the stopping lane
            try:
                gapb=BusesLtd[BusesLtd[:,2]==1][0][0]-bus[0]-parameters.Db
                vbef=BusesLtd[BusesLtd[:,2]==1][0][1]
            # If no bus is found the gapf is not modified
            except:
                pass
    bus[5]=gapf
    bus[6]=gapb
    bus[7]=vbef

# This function calculates the gaps for vehicles on the same lane using only array operations
def calculategaps(Buses,periodic,limits):
    
    # East bound buses, main lane
    index0E=np.where(np.all(np.array([Buses[:,2]==0,Buses[:,12]>0]),axis=0))[0]
    # East bound buses, stopping lane
    index1E=np.where(np.all(np.array([Buses[:,2]==1,Buses[:,12]>0]),axis=0))[0]
    # West bound buses, main lane
    index0W=np.where(np.all(np.array([Buses[:,2]==0,Buses[:,12]<0]),axis=0))[0]
    # West bound buses, stopping lane
    index1W=np.where(np.all(np.array([Buses[:,2]==1,Buses[:,12]<0]),axis=0))[0]
    
    # The forward gap matrix only taking cars into account
    gapf0E=np.roll(Buses[index0E,0],-1)-Buses[index0E,0]-parameters.Db  # The last element is not well defined
    gapf1E=np.roll(Buses[index1E,0],-1)-Buses[index1E,0]-parameters.Db  # The last element is not well defined
    gapf0W=np.roll(Buses[index0W,0],1)-Buses[index0W,0]+parameters.Db  # The first element is not well defined
    gapf1W=np.roll(Buses[index1W,0],1)-Buses[index1W,0]+parameters.Db  # The first element is not well defined

    # Correcting for physical boundaries
    # For buses in the main lane, a new boundary appears at the end of the approximation zone
    gapf0E=np.minimum(gapf0E,Buses[index0E,8]-1.5*parameters.Db-Buses[index0E,0])
    gapf0W=np.maximum(gapf0W,Buses[index0W,8]+1.5*parameters.Db-Buses[index0W,0])
    
    # For buses in the stopping lane, a new boundary appears at the stop or at the end of the stopping lane
    gapf1E=np.minimum(gapf1E,np.minimum(Buses[index1E,8],Buses[index1E,10])-Buses[index1E,0])
    gapf1W=np.maximum(gapf1W,np.minimum(Buses[index1W,8],Buses[index1W,10])-Buses[index1W,0])
    
    # Correcting the negative or positive values in the gaps
    gapf0E[np.where(gapf0E<0)[0]]=0
    gapf1E[np.where(gapf1E<0)[0]]=0
    gapf0W[np.where(gapf0W>0)[0]]=0
    gapf1W[np.where(gapf1W>0)[0]]=0
    
#     Setting the buses in the main lane to surrender priority to buses in the stopping lane which are waiting
#     Buses changing and y==1, and  acc>0, and gapf=0
    indexPriorE=np.where(np.all([Buses[:,18],Buses[:,2]==1,Buses[:,12]>0,Buses[:,3]==0],axis=0))[0]
    for i in indexPriorE:
        try:
            # we find the gap to the car behind a>0 and y==0
            IndexLastE=np.where(np.all(Buses[:i,12]>0,Buses[:i,2]==0))[0][-1]
            gap=Buses[i,0]-Buses[IndexLastE,0]-parameters.Db
            if gap>0 and gap<parameters.Db and Buses[IndexLastE,1]<=3:
                Buses[IndexLastE,3]=0
        except:
            pass
    
    
    # The limits of the system
    xmin=limits[0]
    xmax=limits[1]
    Len=xmax-xmin
    
    # Correcting the ill-defined elements if periodic conditions are set
    if periodic:
        try:
            gapf0E[-1]=Len+Buses[index0E[0],0]-Buses[index0E[-1],0]-parameters.Db
        except:
            pass
        
        try:
            gapf1E[-1]=Len+Buses[index1E[0],0]-Buses[index1E[-1],0]-parameters.Db
        except:
            pass
        
        try:
            gapf0W[0]=-Len-Buses[index0W[0],0]+Buses[index0W[-1],0]+parameters.Db
        except:
            pass
        
        try:
            gapf1W[0]=-Len-Buses[index1W[0],0]+Buses[index1W[-1],0]+parameters.Db
        except:
            pass
        
    # if periodic conditions are not set:
    else:
        try:
            gapf0E[-1]=1000
        except:
            pass
        
        try:
            gapf1E[-1]=1000
        except:
            pass
        
        try:
            gapf0W[0]=-1000
        except:
            pass
        
        try:
            gapf1W[0]=-1000
        except:
            pass
    
    # Setting the forward gaps
    Buses[index0E,3]=gapf0E
    Buses[index1E,3]=gapf1E
    Buses[index0W,3]=gapf0W
    Buses[index1W,3]=gapf1W
    
#    # Setting the backward gaps
    Buses[index0E,4]=-np.roll(gapf0E,1)
    Buses[index1E,4]=-np.roll(gapf1E,1)
    Buses[index0W,4]=-np.roll(gapf0W,-1)
    Buses[index1W,4]=-np.roll(gapf1W,-1)



## This function updates the bus order in the system
def updatebusorder(Buses):
    BusesNew=Buses[Buses[:,0].argsort()]
    return BusesNew

# This function creates all the line changes in the system   
def buschangelane(Buses,Lines,Stations):
    newy=[]
    i=0
    
    index0E=np.where(np.all(np.array([Buses[:,17],Buses[:,12]>0,Buses[:,8]>0,Buses[:,2]==0,Buses[:,18],Buses[:,8]-Buses[:,0]-1.5*parameters.Db<parameters.Dc]),axis=0))[0]
    
    for bus in Buses:
        # default value if no change is performed
        newyaux=bus.y
        # only active buses
        if bus.active:
            ################################################
            if bus.acc>0:  # Eastbound buses
                ###################################################
                # Change from the main lane to the stopping lane
                if bus.nextstop and bus.y==0:       
                    # we check for the buses in the approximation zone
                    if bus.nextstop-bus.x-1.5*parameters.Db<parameters.Dc and bus.changing:
                        # we now calculate the bus gaps
                        bus.gapsl(Buses)
                        # And we check whether the bus is suitable to change
                        if -bus.gapbl>bus.vbef and bus.gapfl>bus.v:
                            # The change is proposed
                            newyaux=1
                            bus.changing=False
                ##################################################
                # Change from the stopping lane to the main lane
                # We check whether the bus already left its station and has some obstacles
                elif bus.nextstationend and bus.y==1 and bus.changing and bus.time-bus.laststop>5 and bus.gapf<=2*parameters.vmax:
                    # We calculate the bus gaps
                    bus.gapsl(Buses)
                    # we check its suitability for change
                    if  bus.gapfl>bus.v and -bus.gapbl>bus.vbef:
                        newyaux=0
                        bus.updatestationend(Stations)  # The next station information is updated
        #########################################################
            elif bus.acc<0:   #Westbound buses
                ########################################################
                #Change from the main lane to the stopping lane
                if bus.nextstop and bus.y==0:
                    # We first check whether the bus is in the approximation zone
                    if  -(bus.nextstop-bus.x+1.5*parameters.Db)<parameters.Dc and bus.changing:
                        # Next we calculate the gaps
                        bus.gapsl(Buses)
                        # Check if the bus is suitable for change
                        if bus.gapbl>-bus.vbef and bus.gapfl<bus.v:
                            newyaux=1
                            bus.changing=False
                ##################################################
                #Change from the stopping lane to the main lane
                # We check whether the bus already left uts station
                elif bus.nextstationend and bus.y==1 and bus.changing and bus.time-bus.laststop>5 and bus.gapf>=-2*parameters.vmax:
                    # We calculate the bus gaps
                    bus.gapsl(Buses)
                    #Check if it is suitable for change
                    if  bus.gapbl>-bus.vbef and bus.gapfl<bus.v:
                        newyaux=0
                        bus.updatestationend(Stations)
        newy.append(newyaux)
        i=i+1
    # All the lanes are updated in parallel
    i=0
    for bus in Buses:
        bus.y=newy[i]
        i=i+1
#
## This function introduces an advance in the system      
#def busadvance(Buses,Lines,Stations):
#    # First all the velocities are updated
#    for bus in Buses:
#        #Advance all buses' time
#        bus.time=bus.time+1
#        # Check whether the bus is advancing
#        if bus.advancing and bus.active:
#            #Check if the velocity is not vmax
#            if abs(bus.v)<abs(parameters.vmax):
#                bus.v=bus.v+bus.acc
#            #Check that velocity is not larger than the gapf
#            if abs(bus.v)>abs(bus.gapf):
#                bus.v=bus.gapf
#            #Apply random desceleration
#            if random.random()<parameters.p and abs(bus.v)>0:
#                bus.v=bus.v-bus.acc/abs(bus.acc)        
#            #if the velocity turns out to be opposite to the accelerations
#            if bus.v/(1.0*bus.acc)<0:
#                    bus.v=0
#                
#    # Now the position is updated
#    for bus in Buses:
#        #Only active buses
#        if bus.active:
#            # If the bus is advancing the position is updated
#            if bus.advancing:
#                #Updating the measured parameters
#                if bus.measuring:                
#                    bus.vsum=bus.vsum+bus.v
#                    bus.nmeas=bus.nmeas+1
#                #Updating the position
#                bus.x=bus.x+bus.v
#                # If after the movement the bus reaches its stop
#                if bus.x==bus.nextstop:
#                    bus.advancing=False
#                    bus.stoptime=0
#                    bus.expectedstoptime=np.random.poisson(bus.AVstoptime)
#                    bus.v=0
#                    bus.updatestop(Lines,Stations)
#            # If the bus is stopped the stopping time is increased
#            else:
#                #Updating the measured parameters
#                if bus.measuring:                
#                    bus.vsum=bus.vsum+0
#                    bus.nmeas=bus.nmeas+1
#                #Updating the stop time
#                bus.stoptime=bus.stoptime+1
#                #if the bus reaches the stoptime
#                if bus.stoptime>bus.expectedstoptime:
#                    bus.advancing=True
#                    bus.changing=True
#                    bus.laststop=bus.time
#    
#                
##This function removes the buses already out of the boundaries of the system
#def removebuses(Buses,limits):
#    # The limits of the system
#    xmin=limits[0]
#    xmax=limits[1]
#    
#    # Applying
#    for bus in Buses:
#        if bus.active:
#            if bus.x<xmin or bus.x>xmax:
#                bus.active=False
#                bus.advancing=False
#                bus.changing=False
#                bus.measuring=False
#                bus.x=parameters.Dh
#                bus.y=-1
#                if bus.box is not None:
#                    bus.box.visible=False
#                    del bus.box
#
#        
## This function needs to be called to introduce periodic boundary conditions
#def applyperiodic(Buses,limits,Lines,Stations):
#    # The limits of the system
#    xmin=limits[0]
#    xmax=limits[1]
#    Len=xmax-xmin
#    
#    for bus in Buses:
#        if bus.active:
#            # If the bus leaves the system in the west
#            if bus.x<xmin:
#                bus.x=bus.x+Len   # x=x+L
#                bus.findnextstop(Lines,Stations)
#            # If the bus leavyes the system on the east
#            elif bus.x>xmax:
#                bus.x=bus.x-Len   # x=x-L
#                bus.findnextstop(Lines,Stations)
#
## This function is called to print the data 
#def printtodata(Buses):
#    for bus in Buses:
#        if bus.active:
#            bus.data=np.append(bus.data,np.array([[bus.time,bus.x,bus.y,bus.v]]), axis=0)
#            
## With this function all active buses enter into measuring mode
#def startmeasuring(Buses):
#    for bus in Buses:
#        if bus.active:
#            bus.measuring=True
#
## With this function all the measurement information is erased
#def clearmeasurement(Buses):
#    for bus in Buses:
#        bus.vsum=0
#        bus.nmeas=0
#                        
#            
## This function retrieves a bus from its ID
#def getbusbyID(Buses,ID):
#    for bus in Buses:
#        if bus.ID==ID:
#            return bus
#    print("The given line ID has not been found in the Buses list")
#    
#    
## importing the anim parameters
#anim=parameters.anim
#
## if anim, import the vpython modules and functions
#if anim:
#    from visual import *
#    # This function needs to be called in case of an animation is wanted
#    def setrepr(Buses):
#        for bus in Buses:
#            if bus.active:
#                if bus.box is None:
#                    bus.box=box(pos=vector(bus.x-0.5*parameters.Db*bus.acc/abs(bus.acc),3*(2*bus.y-4)*bus.acc-0.5*parameters.Dy*bus.acc/abs(bus.acc)), length=parameters.Db, height=parameters.Dy, color=color.red)
#                else:
#                    bus.box.x=bus.x-0.5*parameters.Db*bus.acc/abs(bus.acc)
#                    bus.box.y=3*(2*bus.y-4)*bus.acc-0.5*parameters.Dy*bus.acc/abs(bus.acc)
#    
## This function returns the gap between two buses with ID1 and ID2
#def getgap(Buses,ID1,ID2):
#    if ID2<=ID1:
#        print("WARNING!! Bad usage of the getgap function %d %d"%(ID1,ID2))
#        sys.exit()
#    gap=Buses[ID2].x-Buses[ID1].x-parameters.Db
#    return gap
#    
## This function finds the gap to the next station in the main lane
#def getgapml(Buses,ID):
#    #setting the direction 1: eastbound, -1:westbound
#    direct=Buses[ID].acc/abs(Buses[ID].acc)
#    #calculating the gap
#    gap=Buses[ID].nextstop-1.5*direct*parameters.Db-Buses[ID].x
#    return gap
#    
## This function finds the gap to the next station of station end in the stopping lane
#def getgapsl(Buses,ID):
#    # If there is a next stop
#    if Buses[ID].nextstop:
#        # east bound buses
#        if Buses[ID].acc>0:
#            gap=min(Buses[ID].nextstop,Buses[ID].nextstationend)-Buses[ID].x
#        # west bound buses
#        elif Buses[ID].acc<0:
#            gap=max(Buses[ID].nextstop,Buses[ID].nextstationend)-Buses[ID].x 
#    # If there is not
#    else:   
#        # east bound buses
#        if Buses[ID].acc>0:
#            gap=Buses[ID].nextstationend-Buses[ID].x
#        # west bound buses
#        elif Buses[ID].acc<0:
#            gap=Buses[ID].nextstationend-Buses[ID].x 
#    return gap