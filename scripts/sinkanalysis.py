#!/bin/python

#################################################################
#                                                               #
# Typing sinkanalysis.py in a phantom simulation folder will    #
# analyze the Sink001 and Sink002 files and plot eccentricity,  #
# semimajor axis and pericentre phase as a function of time.    #
#                                                               #
# Enrico Ragusa                                                 #
#                                                               #
#################################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import fnmatch

def loadSink():
        #find sink files in directory
        listFiles1=[]
        listFiles2=[]
        files=os.listdir("./") #can alaso use os.listdir("./")
        for sinkfiles in files:
                if fnmatch.fnmatch(sinkfiles, "*Sink0001*"):
                        listFiles1.append(sinkfiles)
                if fnmatch.fnmatch(sinkfiles, "*Sink0002*"):
                        listFiles2.append(sinkfiles)
        
#        filename= sys.argv[j]+"Sink"+
#       print "Prefix: ",filename

        ######################################################################
        #--------------------------------------------
        #    THIS IS JUST A HACK TO COMPUTE NUMBER OF COLUMNS
        #    Get column number and header lines 
        #    Print the number of columns for each line in a prov file.
        #    Number of columns given when ncols is equal for 4 lines.
        #    Assuming that the header lenght is less than 20 lines
        #----------------------------------------------

        listFiles1.sort()
        listFiles2.sort()
        for filename in listFiles1:
                print filename
                ncols = 0
                nrowh = 0
                n = 21

                command1="awk 'NR<"+ str(n)+" {print NF}' "+filename+">prov"
                os.system(command1)

                prov=np.loadtxt("prov")

                for i in range(n-4):
                    if (prov[i]==ncols and prov[i+1]==ncols and prov[i+2]==ncols and prov[i+3]==ncols): 
                        nrowh=i+1
                        break
                    ncols=prov[i]

                ncols=prov[i]
                os.system("rm prov")  
        #########################################################
                
                if "data1" in locals():
                        dataProv=np.loadtxt(filename,skiprows=int(nrowh))
                        data1=np.vstack((data1,dataProv))
                else:
                        data1=np.loadtxt(filename,skiprows=int(nrowh))

        for filename in listFiles2:
                print filename
                ncols = 0
                nrowh = 0
                n = 21

        ######################################################
                #Again computing number of columns for Sink2 files
                command1="awk 'NR<"+ str(n)+" {print NF}' "+filename+">prov"
                os.system(command1)

                prov=np.loadtxt("prov")

                for i in range(n-4):
                    if (prov[i]==ncols and prov[i+1]==ncols and prov[i+2]==ncols and prov[i+3]==ncols): 
                        nrowh=i+1
                        break
                    ncols=prov[i]

                ncols=prov[i]
                os.system("rm prov")  
        ###################################

                if "data2" in locals():
                        dataProv=np.loadtxt(filename,skiprows=int(nrowh))
                        data2=np.vstack((data2,dataProv))
                else:
                        data2=np.loadtxt(filename,skiprows=int(nrowh))


        return data1,data2

def FindOrbEvo(sink1,sink2):
        x1=sink1[:,1]
        y1=sink1[:,2]
        z1=sink1[:,3]
        vx1=sink1[:,5]
        vy1=sink1[:,6]
        vz1=sink1[:,7]
        M1=sink1[:,4]
        
        x2=sink2[:,1]
        y2=sink2[:,2]
        z2=sink2[:,3]
        vx2=sink2[:,5]
        vy2=sink2[:,6]
        vz2=sink2[:,7]
        M2=sink2[:,4]

        xx1=np.array([x1,y1,z1])
        xx2=np.array([x2,y2,z2])
        vv1=np.array([vx1,vy1,vz1])
        vv2=np.array([vx2,vy2,vz2])
        CMass=np.array([M1*x1+M2*x2,M1*y1+M2*y2,M1*z1+M2*z2])/(M1+M2)
        CMVel=np.array([M1*vx1+M2*vx2,M1*vy1+M2*vy2,M1*vz1+M2*vz2])/(M1+M2)
        
        #Coordinate system of the binary in CM frame
        x1=xx1[0]-CMass[0]
        y1=xx1[1]-CMass[1]
        z1=xx1[2]-CMass[2]

        x2=xx2[0]-CMass[0]
        y2=xx2[1]-CMass[1]
        z2=xx2[2]-CMass[2]

        vx1=vv1[0]-CMVel[0]
        vy1=vv1[1]-CMVel[1]
        vz1=vv1[2]-CMVel[2]

        vx2=vv2[0]-CMVel[0]
        vy2=vv2[1]-CMVel[1]
        vz2=vv2[2]-CMVel[2]

        #coordinate of the satellite in the star frame
        x2S=xx2[0]-xx1[0]
        y2S=xx2[1]-xx1[1]
        z2S=xx2[2]-xx1[2]

        vx2S=vv2[0]-vv1[0]
        vy2S=vv2[1]-vv1[1]
        vz2S=vv2[2]-vv1[2]

        j2S=x2S*vy2S-y2S*vx2S
        R2S=np.sqrt(x2S**2+y2S**2+z2S**2)


        j1=x1*vy1-y1*vx1
        R1=np.sqrt(x1**2+y1**2+z1**2)
        v1Tot=np.sqrt(vx1**2+vy1**2+vz1**2)


        j2=x2*vy2-y2*vx2
        R2=np.sqrt(x2**2+y2**2+z2**2)
        v2Tot=np.sqrt(vx2**2+vy2**2+vz2**2)


        Time=sink1[:,0]
        En=0.5*M1*v1Tot**2+0.5*M2*v2Tot**2-M1*M2/(R1+R2)#-1/2*GMm/a
        a=-0.5*M1*M2/En
        a1=a*M2/(M1+M2)
        a2=a*M1/(M1+M2)
        L=M1*j1+M2*j2
        j=(1./M1+1./M2)*L
        ecc=np.sqrt(1.-j**2/((M1+M2)*a))
        eccConf=np.sqrt(1-j1**2/((M1+M2)*a1**4)*a**3)
        #needs to be computed in the CM centred in the focus of the orbit (i.e. M_*)
        Phase=np.arctan2(-(j2S*vx2S/(M1+M2)+y2S/R2S),j2S*vy2S/(M1+M2)-x2S/R2S)+np.pi

        return Time,ecc,a,Phase

if __name__=="__main__":
        
        drawdirect=True
        sink1,sink2=loadSink()
        Time,ecc,a,Phase=FindOrbEvo(sink1,sink2)


        plt.figure(1)
        plt.plot(Time,ecc)
        plt.title("ecc vs time")
        plt.xlabel("$t$")
        plt.ylabel("$e$")

        plt.figure(2)
        plt.plot(Time,a)
        plt.title("a vs time")
        plt.xlabel("$t$")
        plt.ylabel("$a$")

        plt.figure(3)
        plt.plot(Time,Phase,marker=".",linestyle="")
        plt.title("Phase vs time")
        plt.xlabel("$t$")
        plt.ylabel("$\\Phi$")


        if(drawdirect):
            plt.draw()
            plt.pause(1)
            raw_input("<Hit enter to close the plots>")
            plt.close('all')



