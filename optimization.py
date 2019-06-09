# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 10:16:28 2017

@author: miguelurla
"""

from __future__ import division

import busC
import parameters
import random
import numpy as np
import multiprocessing


##$ Setting all the animation parameters
#if anim:
#    from visual import *
#    def canvas(scene):  # return canvas bounding box, excluding frames
#        bar, d = 30, 8  # title bar and frame thickness for Windows
#        return (int(scene.x+d), int(scene.y+bar), int(scene.width-d), int(scene.height-d))
#
#    scene.autoscale = False
#    scene.fullscreen = True
#
#    scene = display(title="Inertial reference frames",x=0, y=0, width=600, height=400, fov=pi/10,center=vector(1000,0,0), background=vector(1,1,1), range=1000, forward=vector(0,0,-1))
   
   
## Setting the evolution rule
def evolution(Buses,Lines,Stations,limits):
    periodic=True
#    busC.calculategapsl(Buses,periodic,limits)
    busC.buschangelane(Buses,Lines,Stations)
    busC.calculategaps(Buses,periodic,limits)
    busC.busadvance(Buses,Lines,Stations)    
    busC.applyperiodic(Buses,limits,Lines,Stations)



# With this function the system advances from time time to time+dt measuring and averging all the flows (overall, and lines)   
def advancetime(time,dt,buses,lines,stations,limits,LineIDs,Len,flowS,flowlinesS):
    #updating the time
    t1=time
    # Erasing the bus measurements
    busC.clearmeasurement(buses)
    # The system evolves for one cycle
    while(time<(t1+dt)):
        #Evolving the system
        buses=busC.updatebusorder(buses)
        evolution(buses,lines,stations,limits)
        #Evolving the time
        time=time+1
    
    # Calculating the average velocity and flux
    velav=0
    # The average velocity per line
    velavlines=[0 for line in LineIDs]
    fluxlines=[]
        
    for bus in buses:
        #getting the bus average speed         
        bus.velav=bus.vsum/bus.nmeas
        # The general flux
        velav=bus.velav+velav
        # Finding the flux for each line
        # retrieving the lineID
        tlineID=bus.lineID
        # finding the index of the id in lines and adding the speed
        index=LineIDs.index(tlineID)
        velavlines[index]=velavlines[index]+bus.velav
    # Normalizing the flow   
    flow=velav/Len  # 1/L*sum(vi)
    # Calculating the flow for each line
    for i in range(len(LineIDs)):
        fluxlines.append(velavlines[i]/Len)
    flowS.append(flow)
    flowlinesS.append(fluxlines)
    return time


# This functions calacultaes the average flow and standar deviation
def getaverage(flowS,flowlinesS,Len,prob,Nmeas):
    # First averagin over the overall flow
    flow=np.mean(flowS[-Nmeas:])
    flowSD=np.std(flowS[-Nmeas:],ddof=1)
    ratio=flowSD/flow
    # Now we fit the individual flow for each line
    # First we need to transpose the flowlinesS list
    flowlines=[]
    flowlinesSD=[]
    flowlinesT=map(list,zip(*flowlinesS))
    for flows in flowlinesT:
        # We perform the averages
        flowlines.append(np.mean(flows[-Nmeas:]))
        flowlinesSD.append(np.std(flows[-Nmeas:]))
    # Calculating the velocity per line
    vellines=[]
    vellinesSD=[]
    for i in range(len(flowlinesT)):
        if prob[i]!=0:
            vellines.append(flowlines[i]*Len/prob[i])
            vellinesSD.append(flowlinesSD[i]*Len/prob[i])
        else:
            vellines.append(None)
            vellinesSD.append(None)
    return [ratio,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]
    
    
## Fitting the flow values and obtaining the values and standard deviations   
#def getfits(timeS,flowS,flowlinesS,Len,prob):
#    #
#    AV=np.mean(flowS)
#    SD=np.std(flowS)
#    print(SD/AV)
#    # First, fitting the overall flow to exponential
#    try:
#        [p,pcov]=curve_fit(lambda t,a,b,c: c+a*np.exp(-t/b),  timeS,  flowS, p0=(flowS[0],10000,flowS[-1]))
#        plt.plot( timeS,  flowS)
#        plt.plot(timeS,p[2]+p[0]*np.exp(-timeS/p[1]) )
#    except:
#        plt.plot( timeS,  flowS)
#        print("ERROR! The is a problem with the flow fitting process in getfits")
#    # The results
#    decay=p[0]
#    flow=p[2]
#    flowSD=np.sqrt(pcov[2][2])
#    
#    # Now we fit the individual flow for each line
#    # First we need to transpose the flowlinesS list
#    flowlines=[]
#    flowlinesSD=[]
#    flowlinesT=map(list,zip(*flowlinesS))
#    for flows in flowlinesT:
#        # We perform the fit
#        try:
#            [p,pcov]=curve_fit(lambda t,a,b,c: c+a*np.exp(-t/b),  timeS,  flows, p0=(flows[0],10000,flows[-1]))
#            #plt.plot( timeS,  flows)
#            #plt.plot(timeS,p[2]+p[0]*np.exp(-timeS/p[1]) )
#        except:
#            #plt.plot( timeS,  flows)
#            print("ERROR! The is a problem with the flowlines fitting process in getfits")
#        # The results
#        flowlines.append(p[2])
#        flowlinesSD.append(np.sqrt(pcov[2][2]))
#    
#    # Calculating the velocity per line
#    vellines=[]
#    vellinesSD=[]
#    for i in range(len(flowlinesT)):
#        if prob[i]!=0:
#            vellines.append(flowlines[i]*Len/prob[i])
#            vellinesSD.append(flowlinesSD[i]*Len/prob[i])
#        else:
#            vellines.append(None)
#            vellinesSD.append(None)
#            
#    return [decay,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]
        
### singlefluxcalculation()
#def onesinglerun(S):
#    # decomposing the list S
#    lines=S[0]
#    stations=S[1]
#    dist=S[2]
#    limits=S[3]
#    stime=S[4]
#    linelist=S[5]
#    LineIDs=S[6]
#    p=S[7]
#    
#    if len(p)!=len(LineIDs):
#        print("ERROR!!! The length of p and LineIDs is different in onesinglerun")
#        print("p:")
#        print(p)
#        print("LineIDs")
#        print(LineIDs)
#        
#    # The simulation times and number of processes to average
#    dt=2000   # The measured interval
#    t1=5000  # The initial dead time
#    Nmeas=10 # The measured intervals of length dt after t1
#    dNmeas=5 # The additional intervals in case the system has not reached stability
#    Nav=10  # The number of cycles to average
#    tmax=200000
#
#    # The lengths
#    xmin=limits[0]
#    xmax=limits[1]
#    Len=xmax-xmin
#        
#   # Creating the bus matrix
#    buses=[]
#    # Shuffling the bus distribution
#    random.shuffle(linelist)
#    # Creating the time loop
#    time=0
#
#    # Introducing the buses    
#    pos=parameters.Db   
#    for line in linelist:
#        if pos>xmax:
#            print("There is a bus out of bounds at position %d"%pos)
#        buses.append(busC.busC(line,pos,0,lines,stations,time,stime))
#        pos=pos+dist+parameters.Db    
#    # The initial time loop
#    while time<t1:
#        #Evolving the system
#        buses=busC.updatebusorder(buses)
#        evolution(buses,lines,stations,limits)
#        #Evolving the time
#        time=time+1
#    # Activating the bus measurements
#    busC.startmeasuring(buses)
#        
#    timeS=[time+(n+0.5)*dt for n in range(Nmeas)] # The list of times, middle of the interval
#    flowS=[]  # The list of averaged flows
#    flowlinesS=[]  # The list of averaged flows per line
#    
#    # The system evolves until tmax
#    for i in range(Nmeas):
#       time=advancetime(time,dt,buses,lines,stations,limits,LineIDs,Len,flowS,flowlinesS)
#    
#    # Fitting the data
#    [ratio,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=getaverage(flowS,flowlinesS,Len,p,Nav)
#    
#    while ratio>0.01 and time<tmax:
#        timeS=timeS+[time+(n+0.5)*dt for n in range(dNmeas)] # The list of times, middle of the interval
#        # The system evolves dNmeas more cycles 
#        for i in range(dNmeas):
#            time=advancetime(time,dt,buses,lines,stations,limits,LineIDs,Len,flowS,flowlinesS)
#        # Fitting the data
#        [ratio,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=getaverage(flowS,flowlinesS,Len,p,Nav)
#    del buses
#    return [flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]

## singlefluxcalculation()
def onesinglerun(S):
    # decomposing the list S
    lines=S[0]
    stations=S[1]
    dist=S[2]
    limits=S[3]
    stime=S[4]
    linelist=S[5]
    LineIDs=S[6]
    p=S[7]
    
    if len(p)!=len(LineIDs):
        print("ERROR!!! The length of p and LineIDs is different in onesinglerun")
        print("p:")
        print(p)
        print("LineIDs")
        print(LineIDs)
        
    # The simulation times and number of processes to average
    dt=2000   # The measured interval
    t1=5000  # The initial dead time
    Nmeas=10 # The measured intervals of length dt after t1
    dNmeas=5 # The additional intervals in case the system has not reached stability
    Nav=10  # The number of cycles to average
    tmax=200000

    # The lengths
    xmin=limits[0]
    xmax=limits[1]
    Len=xmax-xmin
        
    # Shuffling the bus distribution
    random.shuffle(linelist)
    # Creating the time loop
    time=0

    # Introducing the buses    
    pos=parameters.Db   
    for line in linelist:
        if pos>xmax:
            print("There is a bus out of bounds at position %d"%pos)
        try:
            buses=np.append(buses,busC.busC(line,pos,0,lines,stations,time,stime))
        except:
            buses=np.array(busC.busC(line,pos,0,lines,stations,time,stime))
        pos=pos+dist+parameters.Db 
        
    # The initial time loop
    while time<t1:
        #Evolving the system
        buses=busC.updatebusorder(buses)
        evolution(buses,lines,stations,limits)
        #Evolving the time
        time=time+1
    # Activating the bus measurements
    busC.startmeasuring(buses)
        
    timeS=[time+(n+0.5)*dt for n in range(Nmeas)] # The list of times, middle of the interval
    flowS=[]  # The list of averaged flows
    flowlinesS=[]  # The list of averaged flows per line
    
    # The system evolves until tmax
    for i in range(Nmeas):
       time=advancetime(time,dt,buses,lines,stations,limits,LineIDs,Len,flowS,flowlinesS)
    
    # Fitting the data
    [ratio,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=getaverage(flowS,flowlinesS,Len,p,Nav)
    #print([flowlines,flowlinesSD])
    
    while ratio>0.01 and time<tmax:
        timeS=timeS+[time+(n+0.5)*dt for n in range(dNmeas)] # The list of times, middle of the interval
        # The system evolves dNmeas more cycles 
        for i in range(dNmeas):
            time=advancetime(time,dt,buses,lines,stations,limits,LineIDs,Len,flowS,flowlinesS)
        # Fitting the data
        [ratio,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=getaverage(flowS,flowlinesS,Len,p,Nav)
        #print([flowlines,flowlinesSD])
    del buses
    return [flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]

#
#
### singlefluxcalculation()
#def onesinglerun2(S,q,Ndts):
#    # decomposing the list S
#    lines=S[0]
#    stations=S[1]
#    dist=S[2]
#    limits=S[3]
#    stime=S[4]
#    linelist=S[5]
#    LineIDs=S[6]
#    p=S[7]
#    
#    # The simulation times and number of processes to average
#    t1=2000   
#    dt=500
#
#    
#    # The lengths
#    xmin=limits[0]
#    xmax=limits[1]
#    Len=xmax-xmin
#          
#
#   # Creating the bus matrix
#    buses=[]
#    # Shuffling the bus distribution
#    random.shuffle(linelist)
#    # Creating the time loop
#    time=0
#
#    # Introducing the buses    
#    pos=parameters.Db   
#    for line in linelist:
#        if pos>xmax:
#            print("There is a bus out of bounds at position %d"%pos)
#        buses.append(busC.busC(line,pos,0,lines,stations,time,stime))
#        pos=pos+dist+parameters.Db    
#    #The real density
##        densityreal=nbuses/(xmax-xmin)
##        print("Simulation started: rel_dens=%f, buses=%d"%(densityreal*parameters.Db,nbuses))
#    
#    # The initial time loop
#    while time<t1:
##            # If animation is set
##            if anim:
##                rate(20)
##                busC.setrepr(buses)  # This updates the animation
#        #Evolving the system
#        buses=busC.updatebusorder(buses)
#        evolution(buses,lines,stations,limits)
##            #The time check
##            if time%1000==0:
##                print("Elapsed time is %d"%time)
#        #Evolving the time
#        time=time+1
#    # Activating the bus measurements
#    busC.startmeasuring(buses)
#    # The system evolves for one cycle
#    while(time<(t1+dt)):
##            # If animation is set
##            if anim:
##                rate(20)
##                busC.setrepr(buses)  # This updates the animation
#        #Evolving the system
#        buses=busC.updatebusorder(buses)
#        evolution(buses,lines,stations,limits)
##            #The time check
##            if time%1000==0:
##                print("Elapsed time is %d"%time)
#        #Evolving the time
#        time=time+1
#        
#    # Calculating the average velocity and flux
#    velav=0
#    # The average velocity per line
#    velavlines=[0 for line in LineIDs]
#    fluxlines=[]
#    vellines=[]
#        
#    for bus in buses:
#        #getting the bus average speed         
#        bus.velav=bus.vsum/bus.nmeas
#        # The general flux
#        velav=bus.velav+velav
#        # Finding the flux for each line
#        # retrieving the lineID
#        tlineID=bus.lineID
#        # finding the index of the id in lines and adding the speed
#        index=LineIDs.index(tlineID)
#        velavlines[index]=velavlines[index]+bus.velav
#    # Normalizing the flow   
#    flux=velav/Len  # 1/L*sum(vi)
#    # Calculating the flow and speed for each line
#    for i in range(len(LineIDs)):
#        fluxlines.append(velavlines[i]/Len)
#        if p[i]!=0:
#            vellines.append(velavlines[i]/p[i])
#        else:
#            vellines.append(None)
#    
#    # Defining the list of parameters
#    fluxeS=[flux]
#    fluxlinesS=[fluxlines]
#    vellinesS=[vellines]
#    coinc=0  # The number of good coincidences
#    while(coinc<5):
#        #updating the time
#        t1=time
#        # Erasing the bus measurements
#        busC.clearmeasurement(buses)
#        # The system evolves for one cycle
#        while(time<(t1+dt)):
#    #            # If animation is set
#    #            if anim:
#    #                rate(20)
#    #                busC.setrepr(buses)  # This updates the animation
#            #Evolving the system
#            buses=busC.updatebusorder(buses)
#            evolution(buses,lines,stations,limits)
#    #            #The time check
#    #            if time%1000==0:
#    #                print("Elapsed time is %d"%time)
#            #Evolving the time
#            time=time+1
#            
#        # Calculating the average velocity and flux
#        velav=0
#        # The average velocity per line
#        velavlines=[0 for line in LineIDs]
#        fluxlines=[]
#        vellines=[]
#            
#        for bus in buses:
#            #getting the bus average speed         
#            bus.velav=bus.vsum/bus.nmeas
#            # The general flux
#            velav=bus.velav+velav
#            # Finding the flux for each line
#            # retrieving the lineID
#            tlineID=bus.lineID
#            # finding the index of the id in lines and adding the speed
#            index=LineIDs.index(tlineID)
#            velavlines[index]=velavlines[index]+bus.velav
#        # Normalizing the flow   
#        flux=velav/Len  # 1/L*sum(vi)
#        # Calculating the flow and speed for each line
#        for i in range(len(LineIDs)):
#            fluxlines.append(velavlines[i]/Len)
#            if p[i]!=0:
#                vellines.append(velavlines[i]/p[i])
#            else:
#                vellines.append(None)
#        # Checking the statistics:
#        fluxeS.append(flux)  # Adding the flux to the list
#        fluxlinesS.append(fluxlines)
#        vellinesS.append(vellines)
#        fluxav=np.average(fluxeS[-10:])  # Calcultaing the average flow for the last 5 scans
#        diff=np.abs((flux-fluxav)/fluxav) # Calculating the difference between the average and the lasty calculated flow
#        if (diff<0.007): # If the difference is a small proportion of the average
#            coinc=coinc+1
#        #print([flux,fluxSD/flux,diff,time])
#    
#    # Averaging over the last 10 scans
#    fluxlinesav=np.average(fluxlinesS[-10:],axis=0)
#    vellinesav=[]
#    vellinesS=map(list,zip(*vellinesS)) # We need the traspose of the speeds
#    # Calculating the flow and speed for each line
#    for i in range(len(LineIDs)):
#        if p[i]!=0:
#            vellinesav.append(np.average(vellinesS[i][-10:]))
#        else:
#            vellinesav.append(None)
#        
#    #velav=velav/len(buses)  # 1/N sum(vi)
#    del buses
#    #print([flux,fluxlines])
#    q.put([flux,fluxlinesav,vellinesav])

#
### One step cycle
#def findfluxav(lines,stations,dist,p,LineIDs,limits,stime,Ntimes,Nbuses):
#
#    if len(p)!=len(LineIDs):
#        print("The sizes of p and of LineIDs are different in findfluxav")
#        return
#    
#    # If we are out of the interest region, the flow is always 0
#    if sum(p)!=Nbuses:
#        return 0
#   
#    # Creating the distribution list
#    linelist=[]
#    for i in range(len(p)):
#        for j in range(int(p[i])):
#            linelist.append(LineIDs[i]) 
#    
#    # Creating the parameter vector
#    S=[lines,stations,dist,limits,stime,linelist,LineIDs,p]
#    
#    # We first calibrate the system
#    Ndist=calibratesystem(S)
#    
#    # Creating the queue
#    fluxes=Queue()
#    # The process list
#    procs=[]
#
#    # Creating the processes
#    for i in range(Ntimes):
#        #Creating one iteration process
#        ite=Process(target=onesinglerun,args=(S,fluxes))
#        ite.start()
#        procs.append(ite)
#        
#    #Waiting for all processes to stop
#    for process in procs:
#        process.join()
#    
#    #Finding the averages
#    fluxcount=[]
#    fluxAV=0         # The overall averaged flux
#    fluxlinesAV=[0 for i in range(len(LineIDs))]   # The flux per line
#    fluxlinescount=[[] for i in range(len(LineIDs))]
#    vellinesAV=[0 for i in range(len(LineIDs))]  # Average speed per line
#    for i in range(Ntimes):
#        [flux,fluxlines,vellines]=fluxes.get()
#        fluxcount.append(flux)
#        fluxAV=fluxAV+flux
#        for i in range(len(LineIDs)):
#            fluxlinesAV[i]=fluxlinesAV[i]+fluxlines[i]
#            fluxlinescount[i].append(fluxlines[i])
#            if p[i]!=0:
#                vellinesAV[i]=vellinesAV[i]+vellines[i]
#                
#    #print(fluxlinesAV)
#    #print(fluxlinescount)
#    #Calculating the standar deviation
#    fluxSD=0
#    fluxlinesSD=[0 for i in range(len(LineIDs))]
#    
#    for flux in fluxcount:
#        fluxSD=fluxSD+(flux-fluxAV/len(fluxcount))**2
#    
#    # The calculation for each line
#    for i in range(len(LineIDs)):
#        for flux in fluxlinescount[i]:
#            fluxlinesSD[i]=fluxlinesSD[i]+(flux-fluxlinesAV[i]/len(fluxlinescount[i]))**2       
#        
#    #Normalizing
#    fluxAV=fluxAV/(Ntimes*1.0)
#    fluxSD=np.sqrt(fluxSD/(Ntimes-1))  
#    for i in range(len(LineIDs)):
#      fluxlinesAV[i]=fluxlinesAV[i]/(Ntimes*1.0)
#      fluxlinesSD[i]=np.sqrt(fluxlinesSD[i]/(Ntimes-1)) 
#      if p[i]!=0:
#          vellinesAV[i]=vellinesAV[i]/(Ntimes*1.0)
#      else:
#          vellinesAV[i]=None
#    
##    print(p)ge
##    print("%f %f %d"%(fluxAV,fluxSD,len(fluxcount)))
#    return [fluxAV,fluxSD,fluxlinesAV,fluxlinesSD,vellinesAV]

# This function averages over mixed buses distributions with a given mixing probabibility
#def getmixedflux(prob,Nbuses,dist,lines,stations,LineIDs,limits,stime,Ntimes):    
#
#    if len(LineIDs)!=2:
#        print("Error! getmixedflux function called for different than two lines in LineIDs")
#        return
#    
#    # Creating the queue
#    fluxes=Queue()
#    # The process list
#    procs=[]
#
#    # Creating the processes
#    for i in range(Ntimes):
#        # The final bus distribution
#        p=[0,0]
#        # Creating the distribution list
#        linelist=[]
#        for j in range(Nbuses):
#            r=random.random()
#            if r<prob:
#                p[0]=p[0]+1
#                linelist.append(LineIDs[0])
#            else:
#                p[1]=p[1]+1
#                linelist.append(LineIDs[1])
#                
#        # Creating the parameter vector
#        S=[lines,stations,dist,limits,stime,linelist,LineIDs,p]
#                
#        #Creating one iteration process
#        ite=Process(target=onesinglerun,args=(S,fluxes,Ndist))
#        ite.start()
#        procs.append(ite)
#        
#    #Waiting for all processes to stop
#    for process in procs:
#        process.join()
#    
#       #Finding the averages
#    fluxcount=[]
#    fluxAV=0         # The overall averaged flux
#    fluxlinesAV=[0 for i in range(len(LineIDs))]   # The flux per line
#    fluxlinescount=[[] for i in range(len(LineIDs))] 
#    for i in range(Ntimes):
#        [flux,fluxlines,vellines]=fluxes.get()
#        fluxcount.append(flux)
#        fluxAV=fluxAV+flux
#        for i in range(len(LineIDs)):
#            fluxlinesAV[i]=fluxlinesAV[i]+fluxlines[i]
#            fluxlinescount[i].append(fluxlines[i])
#        
#    #Calculating the standar deviation
#    fluxSD=0
#    fluxlinesSD=[0 for i in range(len(LineIDs))]
#    
#    for flux in fluxcount:
#        fluxSD=fluxSD+(flux-fluxAV/len(fluxcount))**2
#    
#    # The calculation for each line
#    for i in range(len(LineIDs)):
#        for flux in fluxlinescount[i]:
#            fluxlinesSD[i]=fluxlinesSD[i]+(flux-fluxlinesAV[i]/len(fluxlinescount[i]))**2       
#        
#    #Normalizing
#    fluxAV=fluxAV/len(fluxcount)
#    fluxSD=np.sqrt(fluxSD/(len(fluxcount)-1))  
#    for i in range(len(LineIDs)):
#      fluxlinesAV[i]=fluxlinesAV[i]/len(fluxlinescount[i])
#      fluxlinesSD[i]=np.sqrt(fluxlinesSD[i]/(len(fluxlinescount[i])-1))    
#    
##    print(p)
##    print("%f %f %d"%(fluxAV,fluxSD,len(fluxcount)))
#    return [fluxAV,fluxSD,fluxlinesAV,fluxlinesSD]


#def genrandomp(Nbuses,dim):
#    # getting dim random numbers
#    pgen=[]
#    for i in range(dim):
#        pgen.append(random.random())
#    #normalizing to number of buses
#    pgen=[int(round(Nbuses*p/sum(pgen))) for p in pgen]
#    #getting the difference
#    diff=Nbuses-sum(pgen)
#    # Correcting the difference in any component
#    if diff!=0:
#        r=int(dim*random.random())
#        pgen[r]=pgen[r]+diff
#    #return the result
#    if sum(pgen)!=Nbuses:
#        print("The generated number of buses is different than requested")
#    return pgen

# This function calculates the number of buses and their separation from a given density
def getfromdensity(density,limits):
    # The limits
    xmin=limits[0]
    xmax=limits[1]    
    Len=xmax-xmin
    #Getting the real number of buses
    Nbuses=int(round(density*Len))
    #Getting the real density
    densityreal=Nbuses/(1.0*(Len))
    # Finding the distance between buses
    dist=int(round(1./densityreal-parameters.Db))
    # verifying buses are never generated out of bounds
    posmax=parameters.Db+(Nbuses-1)*(parameters.Db+dist)
    while posmax>xmax:
        dist=dist-1
        posmax=parameters.Db+(Nbuses-1)*(parameters.Db+dist)
    # returning
    return [Nbuses,densityreal,dist]

# This function optimizes the system for a given number of buses and lines
#def optimize(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes, *args):
#    if len(lineIDs)==1:
#        print("There is only one lineID, there is nothing to optimize")
#        return None
#    elif len(lineIDs)==2:
#        [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]=goldensection(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes)
#        print([Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt])
#        return [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]
#    else:
#        # if there is no seed
#        if len(args)<1:
#            [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]=dsa(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes)
#        # if there is a seed, we pass it
#        else:
#            [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]=dsa(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes,args[0])
#        print([Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt])
#        return [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]


# This function optimizes the system for a given number of buses and lines
def optimizeGA(lines,stations,Nbuses,lineIDs,dist,limits,stime,npopu,mprob,ntol,popguess=None):
    if len(lineIDs)==1:
        print("There is only one lineID, there is nothing to optimize")
        return None
    else:
        [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,vlinesSD,popt]=GAoptimize(lines,stations,Nbuses,lineIDs,dist,limits,stime,npopu,mprob,ntol,popguess)
        print([Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt])
        return [Neval,Fmax,FmaxSD,Flines,FlinesSD,vlines,popt]   
#        
## This function optimizes several times and averages the results
#def avoptimizations(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes,Nopt):
#    # The averaged flux
#    FmaxAV=0
#    FmaxSD=0
#    Fluxes=[]
#    # The averaged flux per lines
#    FlinesAV=[0 for i in range(len(lineIDs))]
#    FlinesSD=[0 for i in range(len(lineIDs))]
#    Flines=[]
#    # The averaged speed per lines
#    vlinesAV=[0 for i in range(len(lineIDs))]
#    vlinesSD=[0 for i in range(len(lineIDs))]
#    vNones=[0 for i in range(len(lineIDs))]
#    vlines=[]
#    # The averaged optimum distribution
#    pAV=[0 for i in range(len(lineIDs))]
#    pSD=[0 for i in range(len(lineIDs))]
#    popts=[]
#    # The averaged number of simulations
#    NsimulAV=0
#    NsimulSD=0
#    Nsimuls=[]
#
#    # Loop over several optimization processes
#    for i in range(Nopt):
#        
#        [nsimul,flux,fluxsd,fluxlines,fluxlinessd,vellines,popt]=optimize(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes)
#            
#        
#        # Adding to the averages
#        FmaxAV=FmaxAV+flux/(1.0*Nopt)
#        Fluxes.append(flux)
#
#        # Averaging per line
#        for j in range(len(lineIDs)):
#            FlinesAV[j]=FlinesAV[j]+fluxlines[j]/(1.0*Nopt)
#            if vellines[j]!=None:
#                vlinesAV[j]=vlinesAV[j]+vellines[j]
#            else:    
#                vNones[j]=vNones[j]+1            
#            pAV[j]=pAV[j]+popt[j]/(1.0*Nopt)
#        
#        # Saving all results for the standar deviation calculation
#        Flines.append(fluxlines)
#        vlines.append(vellines)
#        popts.append(popt)
#        
#   
#        # Averaging the number of simulations
#        NsimulAV=NsimulAV+nsimul/(1.0*Nopt)
#        Nsimuls.append(nsimul)
#    
#    # Normalizing the speed averages
#    for j in range(len(lineIDs)):        
#        if (Nopt-vNones[j])>0:
#            vlinesAV[j]=vlinesAV[j]/(1.0*Nopt-vNones[j])
#        else:
#            vlinesAV[j]=None    
#    
#    # Finding the standard deviations
#    # Average Flux
#    for flux in Fluxes:
#        FmaxSD=FmaxSD+(flux-FmaxAV)**2/(Nopt-1.0)
#    FmaxSD=np.sqrt(FmaxSD)
#    
#    # Average flux, speed and distribution
#    for i in range(len(lineIDs)):
#        for fluxline in Flines:
#            FlinesSD[i]=FlinesSD[i]+(fluxline[i]-FlinesAV[i])**2/(Nopt-1.0)
#        for velline in vlines:
#            if (Nopt-vNones[i]-1)>0:
#                if velline[i]!=None:
#                    vlinesSD[i]=vlinesSD[i]+(velline[i]-vlinesAV[i])**2/(Nopt-vNones[i]-1.0)
#            else:
#                vlinesSD[i]=None
#        for p in popts:
#            pSD[i]=pSD[i]+(p[i]-pAV[i])**2/(Nopt-1.0)
#        FlinesSD[i]=np.sqrt(FlinesSD[i])
#        if vlinesSD[i]!=None:
#            vlinesSD[i]=np.sqrt(vlinesSD[i])
#        pSD[i]=np.sqrt(pSD[i])
#        
#    # Average simulation
#    for N in Nsimuls:
#        NsimulSD=NsimulSD+(N-NsimulAV)**2/(Nopt-1.0)
#    NsimulSD=np.sqrt(NsimulSD)
#    
#    # Colecting the results
#    print([pAV,pSD,FmaxAV,FmaxSD,NsimulAV,NsimulSD,FlinesAV,FlinesSD,vlinesAV,vlinesSD])
#    return [pAV,pSD,FmaxAV,FmaxSD,NsimulAV,NsimulSD,FlinesAV,FlinesSD,vlinesAV,vlinesSD]
        
# This function optimizes several times and averages the results
def avoptimizationsGA(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes,Nopt,npopu=9,mprob=0.01,ntol=5):
    # The averaged flux
    FmaxAV=0
    FmaxSD=0
    Fluxes=[]
    # The averaged flux per lines
    FlinesAV=[0 for i in range(len(lineIDs))]
    FlinesSD=[0 for i in range(len(lineIDs))]
    Flines=[]
    # The averaged speed per lines
    vlinesAV=[0 for i in range(len(lineIDs))]
    vlinesSD=[0 for i in range(len(lineIDs))]
    vNones=[0 for i in range(len(lineIDs))]
    vlines=[]
    # The averaged optimum distribution
    pAV=[0 for i in range(len(lineIDs))]
    pSD=[0 for i in range(len(lineIDs))]
    popts=[]
    # The averaged number of simulations
    NsimulAV=0
    NsimulSD=0
    Nsimuls=[]

    # Loop over several optimization processes
    for i in range(Nopt):
        
        [nsimul,flux,fluxsd,fluxlines,fluxlinessd,vellines,popt]=optimizeGA(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes,npopu,mprob,ntol)
            
        
        # Adding to the averages
        FmaxAV=FmaxAV+flux/(1.0*Nopt)
        Fluxes.append(flux)

        # Averaging per line
        for j in range(len(lineIDs)):
            FlinesAV[j]=FlinesAV[j]+fluxlines[j]/(1.0*Nopt)
            if vellines[j]!=None:
                vlinesAV[j]=vlinesAV[j]+vellines[j]
            else:    
                vNones[j]=vNones[j]+1            
            pAV[j]=pAV[j]+popt[j]/(1.0*Nopt)
        
        # Saving all results for the standar deviation calculation
        Flines.append(fluxlines)
        vlines.append(vellines)
        popts.append(popt)
        
   
        # Averaging the number of simulations
        NsimulAV=NsimulAV+nsimul/(1.0*Nopt)
        Nsimuls.append(nsimul)
    
    # Normalizing the speed averages
    for j in range(len(lineIDs)):        
        if (Nopt-vNones[j])>0:
            vlinesAV[j]=vlinesAV[j]/(1.0*Nopt-vNones[j])
        else:
            vlinesAV[j]=None    
    
    # Finding the standard deviations
    # Average Flux
    for flux in Fluxes:
        FmaxSD=FmaxSD+(flux-FmaxAV)**2/(Nopt-1.0)
    FmaxSD=np.sqrt(FmaxSD)
    
    # Average flux, speed and distribution
    for i in range(len(lineIDs)):
        for fluxline in Flines:
            FlinesSD[i]=FlinesSD[i]+(fluxline[i]-FlinesAV[i])**2/(Nopt-1.0)
        for velline in vlines:
            if (Nopt-vNones[i]-1)>0:
                if velline[i]!=None:
                    vlinesSD[i]=vlinesSD[i]+(velline[i]-vlinesAV[i])**2/(Nopt-vNones[i]-1.0)
            else:
                vlinesSD[i]=None
        for p in popts:
            pSD[i]=pSD[i]+(p[i]-pAV[i])**2/(Nopt-1.0)
        FlinesSD[i]=np.sqrt(FlinesSD[i])
        if vlinesSD[i]!=None:
            vlinesSD[i]=np.sqrt(vlinesSD[i])
        pSD[i]=np.sqrt(pSD[i])
        
    # Average simulation
    for N in Nsimuls:
        NsimulSD=NsimulSD+(N-NsimulAV)**2/(Nopt-1.0)
    NsimulSD=np.sqrt(NsimulSD)
    
    # Colecting the results
    print([pAV,pSD,FmaxAV,FmaxSD,NsimulAV,NsimulSD,FlinesAV,FlinesSD,vlinesAV,vlinesSD])
    return [pAV,pSD,FmaxAV,FmaxSD,NsimulAV,NsimulSD,FlinesAV,FlinesSD,vlinesAV,vlinesSD]   
    
## This function implements the golden section algorithm
#def goldensection(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes):
#    # The golden number
#    c=(-1+np.sqrt(5))/2
#    # The tolerance
#    tol=int(max(1,0.01*Nbuses))
#    # If the number of lines is different than 2
#    if len(lineIDs)!=2:
#        print("Golden-Section called for different than two lines")
#        return None
#    # Otherwise
#    else:
#        # The initial interval limits
#        lim=[0,Nbuses]
#        # The step
#        d=c*(lim[1]-lim[0])
#        # The two internal numbers
#        x1=int(round(lim[0]+d))
#        x2=int(round(lim[1]-d))
#        # Creating the distribution vectors
#        p1=[x1,Nbuses-x1]
#        p2=[x2,Nbuses-x2]
#        # Calculating the flux for both
#        [F1,F1sd,F1lines,F1linesSD,v1lines]=findfluxav(lines,stations,dist,p1,lineIDs,limits,stime,Ntimes,Nbuses)
#        [F2,F2sd,F2lines,F2linesSD,v2lines]=findfluxav(lines,stations,dist,p2,lineIDs,limits,stime,Ntimes,Nbuses)
#        # The number of evaluations
#        Neval=2
#        # Print
#        print("%d %f %f"%(x1,F1,F1sd))
#        for i in range(len(lineIDs)):
#            print ("%f %f"%(F1lines[i],F1linesSD[i]))
#        print(v1lines)
#        print("%d %f %f"%(x2,F2,F2sd))
#        for i in range(len(lineIDs)):
#            print ("%f %f"%(F2lines[i],F2linesSD[i]))
#        print(v2lines)
#        # starting the iteration
#        while (lim[1]-lim[0])>tol:            
#            # Comparing the flux
#            if F1>F2:
#                # Printing all maximum information
#                print("%d %f %f"%(x1,F1,F1sd))
#                for i in range(len(lineIDs)):
#                    print("%f %f"%(F1lines[i],F1linesSD[i]))
#                print(v1lines)
#                # Updating tha maximum flux
#                [xmax,Fmax,FmaxSD,Fmaxlines,FmaxlinesSD,vmaxlines]=[x1,F1,F1sd,F1lines,F1linesSD,v1lines]
#                # updating the lowe limit
#                lim[0]=x2
#                # updating the flux for x2
#                [x2,F2,F2sd,F2lines,F2linesSD,v2lines]=[x1,F1,F1sd,F1lines,F1linesSD,v1lines]
#                # updating the upper inner point
#                x1=int(round(lim[0]+c*(lim[1]-lim[0])))
#                # updating the distribution
#                p1=[x1,Nbuses-x1]
#                # calculating the new flux at x1
#                [F1,F1sd,F1lines,F1linesSD,v1lines]=findfluxav(lines,stations,dist,p1,lineIDs,limits,stime,Ntimes,Nbuses)
#                # updating the evaluation counter
#                Neval=Neval+1
#            else:
#                # Printing the maximum information
#                print("%d %f %f"%(x2,F2,F2sd))
#                for i in range(len(lineIDs)):
#                    print ("%f %f"%(F2lines[i],F2linesSD[i]))
#                print(v2lines)
#                # Updating tha maximum flux
#                [xmax,Fmax,FmaxSD,Fmaxlines,FmaxlinesSD,vmaxlines]=[x2,F2,F2sd,F2lines,F2linesSD,v2lines]
#                # updating the upper limit
#                lim[1]=x1
#                # updating the flux for x1
#                [x1,F1,F1sd,F1lines,F1linesSD,v1lines]=[x2,F2,F2sd,F2lines,F2linesSD,v2lines]
#                # updating the upper lower point
#                x2=int(round(lim[1]-c*(lim[1]-lim[0])))
#                # updating the distribution
#                p2=[x2,Nbuses-x2]
#                # calculating the new flux at x2
#                [F2,F2sd,F2lines,F2linesSD,v2lines]=findfluxav(lines,stations,dist,p2,lineIDs,limits,stime,Ntimes,Nbuses)
#                # updating the evaluation counter
#                Neval=Neval+1
#                
#    # The final result
#    popt=[xmax,Nbuses-xmax]
#    # The output
#    return [Neval,Fmax,FmaxSD,Fmaxlines,FmaxlinesSD,vmaxlines,popt]
            
# This functions prints the optimization result to a file
def printresult(wfile,fluxav,fluxsd,p,nsimul):
    wfile.write("%d %f %f"%(nsimul,fluxav,fluxsd))
    for paux in p:
        wfile.write(" %d"%paux)
    wfile.write("\n")

#
#    
## This function finds a reflection on the bundary line for a given bus distrubution
#def dsa_boundaries(p,Nbuses):
#    # deleting the last element from p
#    del p[-1]
#    # in case any of the p is negative, it is reflected
#    for i in range(len(p)):
#        if p[i]<0:
#            p[i]=-p[i]
#  
#    # if p is inside
#    if sum(p)<=Nbuses:
#        pref=list(p)
#    
#    # if p is indeed outside, proceed with the reflection
#    else:
#        # Setting the sign negative to force the entrance to the loop
#        sign=-1
#        while(sign<0):
#            # in case any of the p is negative, it is reflected
#            for i in range(len(p)):
#                if p[i]<0:
#                    p[i]=-p[i]
#            # The normal vector
#            n=[1 for i in range(len(p))]
#            # The distance to the boundary
#            D=2*(sum(i*j for (i,j) in zip(p,n))-Nbuses)/(len(p)**2)
#            # The corrected normal vector
#            n=[D for i in range(len(p))]
#            # The reflected vector
#            pref=[a-b for (a,b) in zip(p,n)]
#            # Converting to integer numbers
#            pref=[int(round(a)) for a in pref]
#            # Correcting for possible outcomes out of the boundaries
#            while(sum(pref)>Nbuses):
#                # a random index is chosen
#                index=int(len(pref)*random.random())
#                # The value is decreased by one only when possible
#                if pref[index]>0:
#                    pref[index]=pref[index]-1
#            # Calculating the sign
#            sign=1
#            for prob in pref:    
#                sign=sign*prob
#            p=list(pref)
#            
#    # Assuring that the number of buses is Nbuses
#    pref.append(Nbuses-sum(pref))
#    return pref
#
## This function returns the initial simplex to start the DSA method
#def dsa_initiate(Nbuses,LineIDs):
#    # The number of vertices of the simplex n+1
#    dim=len(LineIDs)
#    # The central point
#    Nini=int(round(Nbuses/dim))
#    # Generating the initial p0
#    p0=[Nini for i in range(dim)]
#    # Normalizing to the number of buses
#    diff=Nbuses-sum(p0)
#    p0[-1]=p0[-1]+diff
#    # Now we proceed to generate dim different vortices around p0, see 2004-Wolf
#    # The size of the simplex
#    a=Nbuses/(2*dim*np.sqrt(dim-1))
#    # The parameters to generate the best possible initial distribution
#    r=(a/(dim*np.sqrt(2)))*(np.sqrt(dim-1)+dim-1)
#    t=(a/(dim*np.sqrt(2)))*(np.sqrt(dim-1)-1)
#    # The simplex is a list of bus distributions
#    simplex=[]
#    # We add the initial point
#    simplex.append(p0)
#    for i in range(dim-1):
#        # creating the new distribution, see 2004-wolf
#        p=list(p0)
#        p[i]=p[i]+r
#        for j in range(dim):
#            if j!=i:
#                p[j]=p[j]+t
#        p=dsa_integerp(p,Nbuses)
#        simplex.append(p)
#    # checking that everything is within the boundaries
#    for i in range(len(simplex)):
#        simplex[i]=dsa_boundaries(simplex[i],Nbuses)
#    return simplex

## This function returns the initial simplex to start the DSA method
#def dsa_initiate_random(Nbuses,LineIDs):
#    # The number of vertices of the simplex n+1
#    dim=len(LineIDs)
#    # creating the simplex
#    simplex=[]
#    # Generating dim points
#    for i in range(dim):
#        p=[]   # Creating the point
#        for j in range(dim-1):  # Scanning over all dimensions  but one
#            limit=0             # Creating a limit
#            for k in range(j):
#                limit=limit+p[k]  # Limiting the range so that trace(p)<Nbuses
#            if limit<Nbuses-2*dim: # If there are enough posibilities
#                if i>0: # To avoid duplicated copies we check other points
#                    equal=True
#                    while(equal):
#                        candidate=random.randint(0,Nbuses-limit) # Random
#                        equal=False  # We suppose it is different
#                        for m in range(i): # We check if its really different
#                            if simplex[m][j]==candidate:
#                                equal=True 
#                    p.append(candidate)  # The component is added to the point
#                else:
#                    p.append(random.randint(0,Nbuses-limit))
#            elif limit < Nbuses:
#                p.append(random.randint(0,Nbuses-limit)) # We assign and dont care about dupplicates
#            else: # If the limit has been reached, there is no longer room
#                p.append(0)
#        diff=Nbuses-sum(p)  # The last element is introduced to match trave(p)=Nbuses
#        p.append(diff)
#        # Inserting the point in the simplex
#        simplex.append(p)
#        
#    return simplex
        
    
## This function returns the initial simplex to start the DSA method using a seed
#def dsa_initiate_seed(LineIDs,pseed):
#    # The number of vertices of the simplex n+1
#    dim=len(LineIDs)
#    # The number of buses
#    Nbuses=sum(pseed)
#    # Generating the initial p0
#    p0=list(pseed)
#    # Normalizing to the number of buses
#    diff=Nbuses-sum(p0)
#    p0[-1]=p0[-1]+diff
#    # Now we proceed to generate dim different vortices around p0, see 2004-Wolf
#    # The size of the simplex
#    a=Nbuses/(2*dim*np.sqrt(dim-1))
#    # The parameters to generate the best possible initial distribution
#    r=(a/(dim*np.sqrt(2)))*(np.sqrt(dim-1)+dim-1)
#    t=(a/(dim*np.sqrt(2)))*(np.sqrt(dim-1)-1)
#
#    # The simplex is a list of bus distributions
#    simplex=[]
#    # We add the initial point
#    #simplex.append(p0)
#    for i in range(dim):
#        # creating the new distribution, see 2004-wolf
#        p=list(p0)
#        p[i]=p[i]-r
#        for j in range(dim):
#            if j!=i:
#                p[j]=p[j]-t
#        p=dsa_integerp(p,Nbuses)
#        simplex.append(p)
#    # checking that everything is within the boundaries
#    for i in range(len(simplex)):
#        simplex[i]=dsa_boundaries(simplex[i],Nbuses)
#    return simplex


## This function returns an integer round of a floating bus distribution
#def dsa_integerp(p,Nbuses):
#    pnew=[]
#    for prob in p:
#        pnew.append(int(round(prob)))
#    # checking that the total is indeed Nbuses, we correct with the last index
#    diff=Nbuses-sum(pnew)
#    pnew[-1]=pnew[-1]+diff
#    return pnew

## This function finds the centroid of a simplex
#def dsa_getcentroid(simplex,Nbuses):
#    # We create the list vector
#    pcen=[]
#    # We determine the dimension of the system
#    dim=len(simplex[0])
#    # We scan for each line
#    for i in range(dim):
#        pi=0
#        # We sum over every element of the simplex
#        for p in simplex:
#            pi=pi+p[i]/len(simplex)
#        # We add the component to the vector
#        pcen.append(pi)
#    # We normalize and make integer
#    pcen=dsa_integerp(pcen,Nbuses)
#    return pcen
    
## This function takes a simplex and orders it in terms of the values of the flux
#def dsa_sortsimplex(S,lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes):
#    fluxes=[]
#    fluxesSD=[]
#    fluxLines=[]
#    fluxLinesSD=[]
#    vLines=[]
#    for p in S:
#        [flux,fluxsd,fluxlines,fluxlinesSD,vlines]=findfluxav(lines,stations,dist,p,lineIDs,limits,stime,Ntimes,Nbuses)
#        fluxes.append(flux)
#        fluxesSD.append(fluxsd)
#        fluxLines.append(fluxlines)
#        fluxLinesSD.append(fluxlinesSD)
#        vLines.append(vlines)
#    # sorting the simplex and SDs according to the flux
#    S=[p for f,p in sorted(zip(fluxes,S))]
#    fluxLines=[p for f,p in sorted(zip(fluxes,fluxLines))]
#    fluxLinesSD=[p for f,p in sorted(zip(fluxes,fluxLinesSD))]
#    fluxesSD=[p for f,p in sorted(zip(fluxes,fluxesSD))]
#    vLines=[p for f,p in sorted(zip(fluxes,vLines))]
#    # sorting the flux
#    fluxes.sort()
#    return [S,fluxes,fluxesSD,fluxLines,fluxLinesSD,vLines]
    
## This function calculates the reflection of a point with p with respect to pcen
#def dsa_reflection(p,pcen,Nbuses):
#    print("Reflecting")
#    pr=[2*x-y for x,y in zip(pcen,p)]
#    pr=dsa_boundaries(pr,Nbuses)
#    return pr

## This function calculates the expansion of a simplex of point p with respect to pcen
#def dsa_expansion(p,pcen,Nbuses):
#    print("Expanding")
#    pe=[2*x-y for x,y in zip(p,pcen)]
#    pe=dsa_boundaries(pe,Nbuses)
#    return pe
    
## This function calculates the contraction of a simplex of point p with respect to pcen
#def dsa_contraction(p,pcen,Nbuses):
#    print("Contracting")
#    pc=[0.5*(x+y) for x,y in zip(p,pcen)]
#    pc=dsa_boundaries(pc,Nbuses)
#    pc=dsa_integerp(pc,Nbuses)
#    return pc
    
# This function produces a shrink of a simplex with respect to the last component, pmin is needed since it is not included in the simplex at this point
#def dsa_shrink(simplex,pmin,Nbuses):
#    print("Shrinking")
#    newsimplex=[]
#    newsimplex.append(dsa_contraction(pmin,simplex[-1],Nbuses))
#    for i in range(len(simplex)-1):
#       newsimplex.append(dsa_contraction(simplex[i],simplex[-1],Nbuses))
#    # The last point, which corresponds to simplex[-1] is not included
#    return newsimplex
    
#def dsa_getavdist(simplex,Nbuses):
#    pav=dsa_getcentroid(simplex,Nbuses)
#    distav=0
#    for p in simplex:
#        dist=0
#        for i in range(len(p)):
#            dist=dist+(pav[i]-p[i])**2
#        distav=distav+np.sqrt(dist)/len(pav)
#    return [distav,pav]

    
## This function executes the dsa algorithm
#def dsa(lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes, *args):
#    # The convergence constant
#    tol=0.01*Nbuses
#    equalpav=0
#    # Check that the number of lines is larger than two
#    if len(lineIDs)<2:
#        print("Downhill simplex method called for two or less lines!!!")
#        return 
#    # The method should work
#    else:
#        # if there is no seed
#        if len(args)<1:
#            # We inititiate the method
#            simplex=dsa_initiate_random(Nbuses,lineIDs)
#        else:
#            # We inititiate the method with a seed
#            # First we check the number of buses in Nbuses and in pseed is the same
#            pseed=list(args[0])  # The extra argument is the seed
#            if (Nbuses!=sum(pseed)):
#                print("WARNING sum(pseed) is not equal to Nbuses in dsa method")
#            else:
#                simplex=dsa_initiate_seed(lineIDs,pseed)
#        # We sort the simplex according to the flux, S[0] is minimum
#        [simplex,fluxes,fluxesSD,fluxlines,fluxlinesSD,vlines]=dsa_sortsimplex(simplex,lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes)
#        Neval=len(lineIDs)
#        # The convergence is determined by the variation in the flux of the simplex
#        [distav,pav]=dsa_getavdist(simplex,Nbuses)
#        #print(distav)
#        #print(simplex)
#        #while (distav>tol and equalpav<2):
#        counts=0
#        while(counts<30):
#            counts=counts+1
#            # The best point
#            [pmax,fmax,fmaxSD,fmaxlines,fmaxlinesSD,vmaxlines]=[simplex[-1],fluxes[-1],fluxesSD[-1],fluxlines[-1],fluxlinesSD[-1],vlines[-1]]
#            # printing best point information
#            print(distav)
#            print(simplex)
#            print(pmax)
#            print("%f %f"%(fmax,fmaxSD))
#            for i in range(len(lineIDs)):
#                print ("%f %f"%(fmaxlines[i],fmaxlinesSD[i]))
#            print(vmaxlines)
#
#            # The worst point
#            [pmin,fmin]=[simplex[0],fluxes[0]]
#            # The second worst point
#            [pv,fv]=[simplex[1],fluxes[1]]
#            # We remove the worst point and calculate the centroid
#            del simplex[0]
#            del fluxes[0]
#            del fluxlines[0]
#            del fluxesSD[0]
#            del fluxlinesSD[0]
#            del vlines[0]
#            
#            ps=dsa_getcentroid(simplex,Nbuses)
#            # We calculate the reflection point and evaluate the flux
#            pr=dsa_reflection(pmin,ps,Nbuses)
#            [fr,frSD,frlines,frlinesSD,vrlines]=findfluxav(lines,stations,dist,pr,lineIDs,limits,stime,Ntimes,Nbuses)
#            Neval=Neval+1
#            # if the new flux es better than the best
#            if fr>fmax:
#                # Attempt expansion
#                pe=dsa_expansion(pr,ps,Nbuses)
#                [fe,feSD,felines,felinesSD,velines]=findfluxav(lines,stations,dist,pe,lineIDs,limits,stime,Ntimes,Nbuses)
#                Neval=Neval+1
#                # In case the expansion is even better
#                if fe>fr:
#                    simplex.append(pe)
#                    fluxes.append(fe)
#                    fluxlines.append(felines)
#                    fluxesSD.append(feSD)
#                    fluxlinesSD.append(felinesSD)
#                    vlines.append(velines)
#                # In the opposite case we take pr
#                else:
#                    simplex.append(pr)
#                    fluxes.append(fr)
#                    fluxlines.append(frlines)
#                    fluxesSD.append(frSD)
#                    fluxlinesSD.append(frlinesSD)
#                    vlines.append(vrlines)
#            # If the new point is still better than the second worst
#            elif fr > fv:
#                # The reflection is taken
#                simplex.append(pr)
#                fluxes.append(fr)
#                fluxlines.append(frlines)
#                fluxesSD.append(frSD)
#                fluxlinesSD.append(frlinesSD)
#                vlines.append(vrlines)
#            # If fr<= fv, the new point is worst than the second worst
#            else:
#                # If the reflection is still better than the worst point
#                if fr>fmin:
#                    # Execute the outer contraction
#                    pc=dsa_contraction(pr,ps,Nbuses)
#                    [fc,fcSD,fclines,fclinesSD,vclines]=findfluxav(lines,stations,dist,pc,lineIDs,limits,stime,Ntimes,Nbuses)
#                    Neval=Neval+1
#                    # if there is an improvement
#                    if fc>fr:
#                        simplex.append(pc)
#                        fluxes.append(fc)
#                        fluxlines.append(fclines)
#                        fluxesSD.append(fcSD)
#                        fluxlinesSD.append(fclinesSD)
#                        vlines.append(vclines)
#                # If the reflection is even worst than the worst point
#                else:
#                    # Execute innter contraction
#                    pc=dsa_contraction(pmin,ps,Nbuses)
#                    [fc,fcSD,fclines,fclinesSD,vclines]=findfluxav(lines,stations,dist,pc,lineIDs,limits,stime,Ntimes,Nbuses)
#                    Neval=Neval+1
#                    # If there is an improvement
#                    if fc>fmin:
#                        simplex.append(pc)
#                        fluxes.append(fc)
#                        fluxlines.append(fclines)
#                        fluxesSD.append(fcSD)
#                        fluxlinesSD.append(fclinesSD)
#                        vlines.append(vclines)
#            # if so far nothing has worked to improve
#            # The shrink process has to be performed
#            if len(simplex)<len(lineIDs):
#                simplex=dsa_shrink(simplex,pmin,Nbuses)
#                # Finding all fluxes and sorting the simplex
#                [simplex,fluxes,fluxesSD,fluxlines,fluxlinesSD,vlines]=dsa_sortsimplex(simplex,lines,stations,Nbuses,lineIDs,dist,limits,stime,Ntimes)
#                Neval=Neval+len(lineIDs)-1
#                # Adding the pmax, fmax point at the end
#                simplex.append(pmax)
#                fluxes.append(fmax)
#                fluxesSD.append(fmaxSD)
#                fluxlines.append(fmaxlines)
#                fluxlinesSD.append(fmaxlinesSD)
#                vlines.append(vmaxlines)
#                
#            # Updating the average distance and average point
#                
#            pavold=pav
#            [distav,pav]=dsa_getavdist(simplex,Nbuses)
#            if pav==pavold:
#                equalpav=equalpav+1
#            
#            # sorting the simplex and others according to the flux
#            simplex=[p for f,p in sorted(zip(fluxes,simplex))]
#            fluxlines=[p for f,p in sorted(zip(fluxes,fluxlines))]
#            fluxlinesSD=[p for f,p in sorted(zip(fluxes,fluxlinesSD))]
#            fluxesSD=[p for f,p in sorted(zip(fluxes,fluxesSD))]
#            vlines=[p for f,p in sorted(zip(fluxes,vlines))]
#            # sorting the flux
#            fluxes.sort()            
#            
#        # After leaving the while loop  
#        # The best point
#        [pmax,fmax,fmaxSD,fmaxlines,fmaxlinesSD,vmaxlines]=[simplex[-1],fluxes[-1],fluxesSD[-1],fluxlines[-1],fluxlinesSD[-1],vlines[-1]]              
#        # We find the average distribution
#        fsd=np.std(fluxes,ddof=1)
#        fsdlines=np.std(fluxlines,axis=0,ddof=1)
#        # printing best point information
#        print(distav)
#        print(simplex)
#        print(pmax)
#        print("%f %f"%(fmax,fsd))
#        for i in range(len(lineIDs)):
#            print ("%f %f"%(fmaxlines[i],fsdlines[i]))
#        print(vmaxlines)
#        return [Neval,fmax,fsd,fmaxlines,fsdlines,vmaxlines,pmax]


# Introducing a genetic algorithm
# First, given a number of buses we define the working space:
def GAgetnbits(NBuses):
    s=bin(NBuses)
    return len(s)-2

# Getting the binary representaion of an integer
def GAinttobin(n,nbits):
    s=format(n,'#0%db'%(nbits+2))
    s=s[2:]
    if len(s)!=nbits:
        print("Warning, the length of the binary element is different than nbits in inttobin")
    return s

# given a point, translate it to a chromosome
def GAgetChromo(p,nbits):
    chrom=''
    for p0 in p[:-1]: # The chromosome representation of a point takes into account all points but the last one
        chrom=chrom+GAinttobin(p0,nbits)
    return chrom

# given a chromosome,  get p
def GAgetp(chrom,nbits,Nbuses):
    if len(chrom)%nbits!=0:
        print("Error!!! The length of the chromosome is not a multiple of nbits in getp")
        return None
    else:
        p=[]
        while(len(chrom)>0):
            pbin=chrom[:nbits]
            p0=int(pbin,2)
            p.append(p0)
            chrom=chrom[nbits:]
        diff=Nbuses-sum(p)
        if diff<0:
            p.append(0)
            ill=True
        else:
            p.append(diff)
            ill=False
        return [p,ill]
    
# This function checks whether two lists are equal
def areequal(X,Y):
    if len(X)!=len(Y):
        return False
    else:
        equal=True
        for i in range(len(X)):
            if X[i]==Y[i]:
                equal=equal & True
            else:
                equal=equal & False
    return equal
    

# This function creates a distribution with the same ratio of a given distribution
def GAgetpfromguess(Nbuses,popguess):
    Nbusesguess=sum(popguess) # The nbuimber of busses in the old distribution
    #the new distribution
    pnew=[]
    for p in popguess[:-1]:
        pnew.append(int(round(p*Nbuses/(Nbusesguess*1.0))))
    # Adding the last member
    pnew.append(Nbuses-sum(pnew))
    return pnew
    

# initialize a number N of random points inside the space phase
def GAinitialize(Nbuses,LineIDs,npopu, popguess):
    # First we establish the number of bits we are working with
    nbits=GAgetnbits(Nbuses)
    # Getting the dimension
    dim=len(LineIDs)
    # If there is an initial guess, we crete the initial population with that guess
    try:
        if len(popguess)!=len(LineIDs):
            print('Error!!! len(popguess) is different than len(LineIDs) in GAinitialize')
        else:
            population=[GAgetpfromguess(Nbuses,popguess)]
    # Otherwise, the population is empty
    except:
        population=[]
    # Generating npopu points (chromosomes)
    while len(population)<npopu:
        p=[]   # Creating the point
        repeated=True
        while(repeated):
            for j in range(dim-1):  # Scanning over all dimensions  but one
                # In order to keep sum(p)<Nbuses
                if sum(p) < Nbuses:
                    p.append(random.randint(0,Nbuses-sum(p))) # We assign and dont care about dupplicates
                else: # If the limit has been reached, there is no longer room
                    p.append(0)
            diff=Nbuses-sum(p)  # The last element is introduced to match trave(p)=Nbuses
            if diff<0:
                print("Warning, a point out of the phase space has been generated in GAinitialize")
            p.append(diff)
            # CHecking for duplicates
            repeated=False
            for point in population:
                if areequal(p,point):
                    repeated=True
                    p=[] # We start over with p
        # Inserting the point into the population
        population.append(p)
    
    
    return [population,nbits]

# One single run parallel process for the GA algorithm
def GAsinglerun(lines,stations,dist,limits,stime,lineIDs,p,q):
    # Creating the distribution list
    linelist=[]
    for i in range(len(p)):
        for j in range(int(p[i])):
            linelist.append(lineIDs[i])
            
    # Creating S
    S=[lines,stations,dist,limits,stime,linelist,lineIDs,p]
    # Running the simulation
    [flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=onesinglerun(S)
    #exporting the results
    q.put([p,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD])
    

        
# This function gets the fitness of a population and sorts it
def GAgetfitness(population,lines,stations,dist,lineIDs,limits,stime,Nbuses, *args):
    # getting the number of Cores available
    Ncpu=multiprocessing.cpu_count()
    
    fluxes=[]
    fluxSDs=[]
    fluxeslines=[]
    fluxSDslines=[]
    velslines=[]  
    vellinesSDs=[]
    pout=[]  # The population in order of process ending
    print("In GAgetfitness")
    # In case elitism is installed and there is already a best member
    if len(population)%2==1 and len(args)>0:       
        # importing the best values
        pout.append(population[-1])   # The best one is always the first one
        fluxes.append(args[0][0])
        fluxSDs.append(args[0][1])
        fluxeslines.append(args[0][2])
        fluxSDslines.append(args[0][3])
        velslines.append(args[0][4])
        vellinesSDs.append(args[0][5])
        # Now we scan over the remaining part of the population in parallel
        i=0  # The counter over the population
        while i<len(population)-1:
            # Creating  the sublist
            if i+Ncpu>len(population)-1:
                pi=population[i:-1]
            else:
                pi=population[i:i+Ncpu]
            i=i+Ncpu
            
            #Creating the queue
            results=multiprocessing.Queue()
            # The process list
            procs=[]
            # we scan over all points
            for point in pi:
                #Creating one iteration process
                simul=multiprocessing.Process(target=GAsinglerun,args=(lines,stations,dist,limits,stime,lineIDs,point,results))
                simul.start()
                procs.append(simul)
            #Waiting for all processes to stop
            for process in procs:
                process.join()
            
            # Retrieving the results
            for process in procs:
                [pres,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=results.get()
                pout.append(pres)
                fluxes.append(flow)
                fluxSDs.append(flowSD)
                fluxeslines.append(flowlines)
                fluxSDslines.append(flowlinesSD)
                velslines.append(vellines)
                vellinesSDs.append(vellinesSD)
                print([pres,flow,flowSD])
        # sorting the population
        population=[p for fs,p in sorted(zip(fluxes,pout))]
        fluxSDs=[p for fs,p in sorted(zip(fluxes,fluxSDs))]
        fluxeslines=[p for fs,p in sorted(zip(fluxes,fluxeslines))]
        fluxSDslines=[p for fs,p in sorted(zip(fluxes,fluxSDslines))]
        velslines=[p for fs,p in sorted(zip(fluxes,velslines))]
        vellinesSDs=[p for fs,p in sorted(zip(fluxes,vellinesSDs))]
        # sorting the flux
        fluxes.sort()
    # This one applies for the first fitness evaluation or the case without elitism
    else:
        # We scan over the population in parallel
        i=0  # The counter over the population
        while i<len(population):
            # Creating  the sublist
            if i+Ncpu>len(population):
                pi=population[i:]
            else:
                pi=population[i:i+Ncpu]
            i=i+Ncpu
            #Creating the queue
            results=multiprocessing.Queue()
            # The process list
            procs=[]
            # We scan over the entire population
            for point in pi:
                #Creating one iteration process
                simul=multiprocessing.Process(target=GAsinglerun,args=(lines,stations,dist,limits,stime,lineIDs,point,results))
                simul.start()
                procs.append(simul)
            #Waiting for all processes to stop
            for process in procs:
                process.join()
            # Retrieving the results
            for process in procs:
                [pres,flow,flowSD,flowlines,flowlinesSD,vellines,vellinesSD]=results.get()
                pout.append(pres)
                fluxes.append(flow)
                fluxSDs.append(flowSD)
                fluxeslines.append(flowlines)
                fluxSDslines.append(flowlinesSD)
                velslines.append(vellines)
                vellinesSDs.append(vellinesSD)
                print([pres,flow,flowSD])
        # sorting the population
        population=[p for fs,p in sorted(zip(fluxes,pout))]
        fluxSDs=[p for fs,p in sorted(zip(fluxes,fluxSDs))]
        fluxeslines=[p for fs,p in sorted(zip(fluxes,fluxeslines))]
        fluxSDslines=[p for fs,p in sorted(zip(fluxes,fluxSDslines))]
        velslines=[p for fs,p in sorted(zip(fluxes,velslines))]
        vellinesSDs=[p for fs,p in sorted(zip(fluxes,vellinesSDs))]
        # sorting the flux
        fluxes.sort()
    print("Out of GAgetfitness")
    return [population,fluxes,fluxSDs,fluxeslines,fluxSDslines,velslines,vellinesSDs]

# mutating a member of population
def GAmutate(pbin,mprob):
    pbinnew=''
    for b in pbin:
        dice=random.random()
        if dice<mprob:
            if b=='1':
                pbinnew=pbinnew+'0'
            elif b=='0':
                pbinnew=pbinnew+'1'
        else:
            pbinnew=pbinnew+b
    return pbinnew
            
# This function mates two members of the population
def GAmate(pbina,pbinb,mprob):
    if len(pbina)!=len(pbinb):
        print("Error!!! Length of two chromosomes is different in GAmating")
    # getting the length of the binaries
    nbin=len(pbina)
    # obtaining the crossover point
    cross=random.randint(0,nbin+1)
    # performing crossover and mutation
    pbinc=pbina[:cross]+pbinb[cross:]
    pbinc=GAmutate(pbinc,mprob)
    pbind=pbinb[:cross]+pbina[cross:]
    pbind=GAmutate(pbind,mprob)
    return [pbinc,pbind]

# Generating Random Parents
def GAgetRandomParents(population,popfitness):
    if len(population)!=len(popfitness):
        print("Error!! The length of the population and fitness is different in getRandomParents")
    # The probability to find each parent, normalized to the minimum fitness
    minfitness=min(popfitness)
    probs=[(p-minfitness)/(sum(popfitness)-len(popfitness)*min(popfitness)) for p in popfitness]
    # We generate a couple with both elements different
    equal=True
    while(equal):
        # The probability of choosing an element is proportional to its fitness
        couple=np.random.choice(range(len(population)),2,p=probs)
        if areequal(population[couple[0]],population[couple[1]]):
            equal=True
        else:
            equal=False

    return [population[couple[0]],population[couple[1]]]

# Creating a new generation
def GAnewgen(population,popfitness,mprob,nbits,Nbuses):
    newpopulation=[]
    # if not elitism is used, an even number of elements in population
    if len(population)%2==0:
        while(len(newpopulation)<len(population)):
            # To avoid duplicate sons
            aggr=0
            while(aggr==0):
                #The members in the new population
                lenold=len(newpopulation)
                # First, we get random parents
                [pA,pB]=GAgetRandomParents(population,popfitness)
                # We get the parents chromosomes
                pbinA=GAgetChromo(pA,nbits)
                pbinB=GAgetChromo(pB,nbits)
                # We mate the parents
                [sonbinA,sonbinB]=GAmate(pbinA,pbinB,mprob)
                # We translate the sons into points
                [sonA,illA]=GAgetp(sonbinA,nbits,Nbuses)
                [sonB,illB]=GAgetp(sonbinB,nbits,Nbuses)
                # Only one of the sons may result fit
                # Checking son A
                if not illA:
                    duplicate=False
                    # We check it is not a duplicate
                    for point in newpopulation:
                        if areequal(sonA,point):
                            duplicate=True
                    if not duplicate:
                        newpopulation.append(sonA)
                        # In case we reach the limits
                        if len(newpopulation)==len(population):
                            break
                # Checking son B
                if not illB:
                    duplicate=False
                    # We check it is not a duplicate
                    for point in newpopulation:
                        if areequal(sonB,point):
                            duplicate=True
                    if not duplicate:
                        newpopulation.append(sonB)
                aggr=len(newpopulation)-lenold     
                
    else:
        # We keep the most fitted element in the population
        newpopulation.append(population[-1])
        # And start mating
        while(len(newpopulation)<len(population)):
            # To avoid duplicate sons
            aggr=0
            while(aggr==0):
                lenold=len(newpopulation)
                # First, we get random parents
                [pA,pB]=GAgetRandomParents(population,popfitness)
                # We get the parents chromosomes
                pbinA=GAgetChromo(pA,nbits)
                pbinB=GAgetChromo(pB,nbits)
                # We mate the parents
                [sonbinA,sonbinB]=GAmate(pbinA,pbinB,mprob)
                # We translate the sons into points
                [sonA,illA]=GAgetp(sonbinA,nbits,Nbuses)
                [sonB,illB]=GAgetp(sonbinB,nbits,Nbuses)
                # Only one of the sons may result fit
                # Checking son A
                if not illA:
                    duplicate=False
                    # We check it is not a duplicate
                    for point in newpopulation:
                        if areequal(sonA,point):
                            duplicate=True
                    if not duplicate:
                        newpopulation.append(sonA)
                        # In case we reach the limits
                        if len(newpopulation)==len(population):
                            break
                # Checking son B
                if not illB:
                    duplicate=False
                    # We check it is not a duplicate
                    for point in newpopulation:
                        if areequal(sonB,point):
                            duplicate=True
                    if not duplicate:
                        newpopulation.append(sonB)
                aggr=len(newpopulation)-lenold

    return newpopulation      

# Running the optimization
def GAoptimize(lines,stations,Nbuses,lineIDs,dist,limits,stime,npopu,mprob,ntol, popguess):
    # Number of evaluations
    Neval=0
    # Check that the number of lines is larger than two
    if len(lineIDs)<2:
        print("Genetic algorithm method called for two or less lines!!!")
        return 
    # The method should work
    else:
        # We start generating a random population
        [population,nbits]=GAinitialize(Nbuses,lineIDs,npopu,popguess)
        # We calculate the fitness and sort the population
        print([population,nbits])
        [population,fluxes,fluxSDs,fluxeslines,fluxSDslines,velslines,velSDslines]=GAgetfitness(population,lines,stations,dist,lineIDs,limits,stime,Nbuses) # this might have to be changed
        Neval=Neval+npopu-npopu%2 # The number of evaluation is updated
        # We establish the best one
        bestp=population[-1]
        bestf=fluxes[-1]
        bestfSD=fluxSDs[-1]
        bestflines=fluxeslines[-1]
        bestfSDlines=fluxSDslines[-1]
        bestvellines=velslines[-1]
        bestvelSDlines=velSDslines[-1]
        print(population)
        print(bestp)
        print("%f %f"%(bestf,bestfSD))
        for i in range(len(lineIDs)):
            print ("%f %f"%(bestflines[i],bestfSDlines[i]))
        print(bestvellines)
        # We start the iteration
        notimproving=0
        while(notimproving<ntol):
            # we generate a new generation
            population=GAnewgen(population,fluxes,mprob,nbits,Nbuses)
            print("Exits the new generation")
            print(population)
            # We evaluate the fitness and sort the new generation
            # If there is elitism:
            if len(population)%2==1:
                # The bests correspond to the first element after GANewgen
                bests=[bestf,bestfSD,bestflines,bestfSDlines,bestvellines,bestvelSDlines]
                [population,fluxes,fluxSDs,fluxeslines,fluxSDslines,velslines,velSDslines]=GAgetfitness(population,lines,stations,dist,lineIDs,limits,stime,Nbuses,bests)
            # If there is not elitism
            else:
                [population,fluxes,fluxSDs,fluxeslines,fluxSDslines,velslines,velSDslines]=GAgetfitness(population,lines,stations,dist,lineIDs,limits,stime,Nbuses)
            Neval=Neval+npopu-npopu%2 # The number of evaluation is update
            # We check whether there has been an improvement
            if fluxes[-1]>bestf:
                bestf=fluxes[-1]
                bestp=population[-1]
                bestfSD=fluxSDs[-1]
                bestflines=fluxeslines[-1]
                bestfSDlines=fluxSDslines[-1]
                bestvellines=velslines[-1]
                bestvelSDlines=velSDslines[-1]
                notimproving=0
                print("There is an improvement")
            else:
                print("No improvement for this generation")
                notimproving=notimproving+1
            print("The population")
            print(population)
            print("The best specimen so far...")
            print(bestp)
            print("%f %f"%(bestf,bestfSD))
            for i in range(len(lineIDs)):
                print ("%f %f "%(bestflines[i],bestfSDlines[i]))            
        return [Neval,bestf,bestfSD,bestflines,bestfSDlines,bestvellines,bestvelSDlines,bestp]
          

# Finding the derivatives
def GetDerivatives(ps,step,fmax,lines,stations,dist,limits,stime,LineIDs):
    # ps has the structure [p,ph,p2h]
    # if there is no ph, The method does not work
    if ps[1]==None:
        return [None,None]
    # If there is no p2h, a linear approximation is used
    elif ps[2]==None:
        # We create the linelist
        linelist=[]
        for i in range(len(ps[1])):
            linelist=linelist+ps[1][i]*[LineIDs[i]]
        Sh=[lines,stations,dist,limits,stime,linelist,LineIDs,ps[1]]
        Rh=onesinglerun(Sh)
        fh=Rh[0]
#        print(ps[1])
#        print(fh)
#        print(fmax)
#        print((fh-fmax)/step)
        if fh>=fmax:  # If a new max is found
            return [None,None]
        else:
            return [(fh-fmax)/step,None]
    # Otherwise
    else:
        # The flows
        f=[]
        for p in ps[1:]:
            # We create the linelist
            linelist=[]
            for i in range(len(p)):
                linelist=linelist+p[i]*[LineIDs[i]]
            S=[lines,stations,dist,limits,stime,linelist,LineIDs,p]
            R=onesinglerun(S)
            f.append(R[0])
#            print(p)
#            print(R[0])
#            print(fmax)
        # The first derivative
        Dh=(-3*fmax-f[1]+4*f[0])/(2*step)
        # The second derivative
        D2h=(fmax+f[1]-2*f[0])/(step**2)
#        print([Dh,D2h])
        return [Dh,D2h]
            
    
# Creating two alternatives distributions after the step is given
def GAgetErrorij(i,j,p,tol,step,fmax,lines,stations,dist,limits,stime,LineIDs):   # Step can be positive or negative
    pnew1=list(p)
    pnew2=list(p)
    # Checking whether both lists are zero
    if p[i]==0 and p[j]==0:
        return 0
    # Checking whether the transformation goes OffBoundaries
    if pnew1[i]==0 and step<0:   # If there are no buses and we are reducing
#        print("The number of buses is 0, it is not possible to perform a backward error calculation")
        return 0 
    
    elif pnew1[i]<-step:# if there is not enough space to make the calculation
#        print("The avaiable space is smaller than the step, updating the step and performing a first order calculation.")
        step=-pnew1[i]
        pnew1[i]=0  # We change the step size to fit the available space
        pnew1[j]=pnew1[j]-step    #Pnew1[j] is reduced (really aumented because step is negative)  
        pnew2=None  # The second part is not performed
#        print([p,pnew1,pnew2])
        [dh,d2h]=GetDerivatives([p,pnew1,pnew2],step,fmax,lines,stations,dist,limits,stime,LineIDs)
        if dh==None:  # If a new maximum is found
            return step
        # Here comes the error calculation taking into account the tolerance
#        print("The error is %f"%(-tol*fmax/dh))
        else:
            return max(-tol*fmax/dh,step)
        
    elif pnew1[i]<-2*step: # There is not space for the second order calculation
#        print("The available space is not enough to perform the required second order calculation, updating the step size and performing a second order calculation.")
        step=int(-pnew2[i]/2.0) # We change the step size
        pnew1[i]=pnew1[i]+step    # Pnew[i] is increased by step
        pnew2[i]=pnew2[i]+2*step  # Pnew[i] is increased by 2*step
        pnew1[j]=pnew1[j]-step    #Pnew1[j] is reduced
        pnew2[j]=pnew2[j]-2*step  #Pnew1[j] is reduced
#        print([p,pnew1,pnew2])
        [dh,d2h]=GetDerivatives([p,pnew1,pnew2],step,fmax,lines,stations,dist,limits,stime,LineIDs)
        # We solve the quadratic equation
        coeff=[d2h/2,dh,tol*fmax]
        sols=np.roots(coeff)
        # if the solutions are not real, we set the error as the step
        if np.iscomplex(sols[0]):
            return 2*step
        # if the solutions are real
        else:
            # The solution depends on the step and the second derivative
            if step*d2h>0: # step >0 and d2h >0  or  step <0 and d2h<0
                dx=min(sols)
            else:         # step >0  and d2h <0  or step <0 and d2h>0
                dx=max(sols)
            return max(dx,2*step)   # step is negative in this case
        
    else:
        # The step should be adaptative
        dx=1e10
#        adapted=False
        reachlimit=False
        stepacc=step
        while(abs(dx)>=2*abs(stepacc) and reachlimit==False):
            # We apply a correction to the step in case we go off-bounds
            if pnew2[i]+2*step<0: # If the ith component becomes smaller than 0
                step=int(-pnew2[i]/2.0) # We change the step size
                reachlimit=True
                
            elif pnew2[j]-2*step<0: # if the jth component becomes smaller than 0
                step=int(pnew2[j]/2.0) # We change the step size
                reachlimit=True
            stepacc=stepacc+step # We update the accumulated step
            
            if step==0: # if the step has been set to zero
                dx=0
                break
                    
            pnew1[i]=pnew1[i]+step    # Pnew[i] is increased by step
            pnew2[i]=pnew2[i]+2*step  # Pnew[i] is increased by 2*step
            pnew1[j]=pnew1[j]-step    #Pnew1[maxind] is reduced
            pnew2[j]=pnew2[j]-2*step  #Pnew1[maxind] is reduced
            
#            print([p,pnew1,pnew2])
            [dh,d2h]=GetDerivatives([p,pnew1,pnew2],step,fmax,lines,stations,dist,limits,stime,LineIDs)
             # We solve the quadratic equation
            coeff=[d2h/2,dh,tol*fmax]
            sols=np.roots(coeff)
#            print("printing the results of the np.roots")
#            print(sols)
            # if the solutions are not real, we set the error as the step
            if np.iscomplex(sols[0]):
                dx=2.01*step  # The extra .01 is to ensure the while loop is repeated
#            print(step)
            # The solution depends on the step and the second derivative
            elif step*d2h>0: # step >0 and d2h >0  or  step <0 and d2h<0
                dx=min(sols)
            else:         # step >0  and d2h <0  or step <0 and d2h<0
                dx=max(sols)    
#            adapted=True
#            print("The results of the iteration")
#            print([step,stepacc,dx])
        return dx
        

    
    

# This method calculates the error bar based on a simple derivative process
def GAgetErrors(lines,stations,dist,p,LineIDs,limits,stime,Nbuses,flow,tol):
    # Finding the initial step guess
    step=max(5,int(round(Nbuses/20)))
    # List of errors
    pErrorF=[]
    pErrorB=[]

    # The list of covariances
    CovF=[[] for i in range(len(p))]
    CovB=[[] for i in range(len(p))]
    
    # scanning over all elements of p
    for i in range(len(p)):
        # we perform another scan over all elements in p such that j>i
        for j in range(i+1,len(p)):
            # Forward error
#            print("Calculating forward error for lines %d and %d"%(LineIDs[i],LineIDs[j]))
            if p[i]<p[j]:  # The GAgeterror function assumes p[i]<p[j]
                dp=GAgetErrorij(i,j,p,tol,step,flow,lines,stations,dist,limits,stime,LineIDs)
#                print("Forward error is %f"%dp)
                if dp*step<0: # If the sign is wrong
                    print("Warning!!! There is something wrong with the error calculation")
                CovF[i].append(dp)
                CovB[j].append(-dp)
            else:  # In this case we invert the function call
                dp=GAgetErrorij(j,i,p,tol,step,flow,lines,stations,dist,limits,stime,LineIDs)
#                print("Forward error is %f"%dp)
                if dp*step<0: # If the sign is wrong
                    print("Warning!!! There is something wrong with the error calculation")
                CovF[j].append(dp)
                CovB[i].append(-dp)
        
            
            # we find the backward error
#            print("Calculating backward error for lines %d and %d"%(LineIDs[i],LineIDs[j]))
            if p[i]<p[j]:  # The GAgeterror function assumes p[i]<p[j]
                dp=GAgetErrorij(i,j,p,tol,-step,flow,lines,stations,dist,limits,stime,LineIDs)
#                print("Backward error is %f"%dp)
                if dp*step>0: # If the sign is wrong
                    print("Warning!!! There is something wrong with the error calculation")
                pErrorB.append(abs(dp))
                CovF[j].append(-dp)
                CovB[i].append(dp)
            else:
                dp=GAgetErrorij(j,i,p,tol,-step,flow,lines,stations,dist,limits,stime,LineIDs)
#                print("Backward error is %f"%dp)
                if dp*step>0: # If the sign is wrong
                    print("Warning!!! There is something wrong with the error calculation")
                pErrorB.append(abs(dp))
                CovB[j].append(dp)
                CovF[i].append(-dp)
                
    # We calculate the error for the last line, the maximum
    pErrorF=[np.amax(np.abs(A)) for A in CovF]
    pErrorB=[np.amax(np.abs(A)) for A in CovB]
    print("Printing results in GA GetErrors")
    print(CovF)
    print(CovB)
    print(pErrorF)
    print(pErrorB)
    return [pErrorF,pErrorB]
    
