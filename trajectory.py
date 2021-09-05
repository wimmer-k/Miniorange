#!/usr/bin/env python3
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Arc, Circle


ME = 0.510998   # MeV/c^2
QE = -1.602176E-19  # C
Mass2kg = 1.782662E-30 # kg
Energy2kg = 1.602176634E-19 # J = Nm/s^2

def main(argv):

    print( )
    
    g = -50
    b = +50

    BF = 0.150 # T
    EE = 1 # MeV

    eradius = np.sqrt(2*ME*Mass2kg * EE*1E6 *Energy2kg) / -QE / BF *1000
    
    # wedge
    wradius = 50
    wangle = 75
    #absorber
    absx = [-30,0,20]
    absy = [5,10,5]
    

    dashes = [8,8]
    
    fig = plt.figure(figsize=(8,4))
    ax = plt.subplot(1,1,1)
    plt.gca().axis([g-10,b+20,-wradius-10,wradius+10])
    ax.set_aspect('equal', 'box')
    # axis
    p, = plt.plot([g-10,b+20],[0,0],color='k',linewidth=1)
    p.set_dashes(dashes)

    # wedge
    wedge = Wedge([0,0],wradius,-wangle/2+90,+wangle/2+90,color="darkgray")
    plt.gca().add_patch(wedge)
    wedge = Wedge([0,0],wradius,-wangle/2-90,+wangle/2-90,color="darkgray")
    plt.gca().add_patch(wedge)
     
    # absorber
    coords = np.concatenate((absx,absx[::-1],absy, [-y for y in absy[::-1]]))
    coords = np.reshape(coords,(6,2), order='F')
    absorber = Polygon(coords,color="black")
    plt.gca().add_patch(absorber)

    
    # trajectories
    #for a in range(-90,95,10):
    #dx = 0.1
    for a in range(0,65,10):
        print("-----------------------------------")
        a = a*np.pi/180
        ## 
        #xx = np.linspace(g,b,(b-g)/dx)
        #yy = []
        #for x in xx:
        #    if len(yy) == 0:
        #        yy.append(0)
        #        continue
        #
        #    y = yy[-1]+np.tan(a)*dx
        #    #print(x,y,-x,-x/np.tan(wangle/2*np.pi/180))
        #    if y < -x/np.tan(wangle/2*np.pi/180):
        #        yy.append(y)
        #    else:
        #        yy.append(0)
        dt = 0.1
        tt = np.linspace(0,500,500/dt)
        xx = [g]
        yy = [0]
        t0 = -1
        Omega = 0.03
        vx1 = None
        vy1 = None
        for t in tt:
            x = xx[-1]+np.cos(a)*dt
            y = yy[-1]+np.sin(a)*dt
            #print(x,y,-x,-x/np.tan(wangle/2*np.pi/180))
            if t0<0 and y < -x/np.tan(wangle/2*np.pi/180):
                xx.append(x)
                yy.append(y)
            else:
                if t0<0:
                    t0=t
                x = xx[-1] + (np.cos(a)*np.cos(Omega*(t-t0)) + np.sin(a)*np.sin(Omega*(t-t0)))*dt
                y = yy[-1] + (np.sin(a)*np.cos(Omega*(t-t0)) - np.cos(a)*np.sin(Omega*(t-t0)))*dt
                if x < y*np.tan(wangle/2*np.pi/180):
                    xx.append(x)
                    yy.append(y)
                else:
                    if vx1 == None or vy1 == None:
                        vx1 = (np.cos(a)*np.cos(Omega*(t-t0)) + np.sin(a)*np.sin(Omega*(t-t0)))
                        vy1 = (np.sin(a)*np.cos(Omega*(t-t0)) - np.cos(a)*np.sin(Omega*(t-t0)))
                    
                    x = xx[-1] + vx1*dt
                    y = yy[-1] + vy1*dt
                    xx.append(x)
                    yy.append(y)
       
        
        plt.plot(xx,yy)
        ## analytical
        #
        #if abs(np.tan(a)*(abs(g-absx[0]))) < absy[0]:
        #    xintersect = absx[0]
        #else:
        #    if a >0:
        #        xintersect = g*np.tan(a)*np.tan(wangle/2*np.pi/180) / (1+np.tan(a)*np.tan(wangle/2*np.pi/180))
        #    else:
        #        xintersect = g*np.tan(a)*np.tan(-wangle/2*np.pi/180) / (1+np.tan(a)*np.tan(-wangle/2*np.pi/180))
        #
        #        
        #
        #if abs(np.tan(a)*(abs(g-xintersect))) > wradius*np.cos(wangle/2*np.pi/180):
        #    xintersect = 1000
        #
        #yintersect = (xintersect-g)*np.tan(a)
        #t = plt.plot([g,xintersect],[0, yintersect],lw=1)
        #if xintersect == 1000 or xintersect == absx[0]:
        #    continue
        #
        #print( xintersect , yintersect )
        #
        ### https://en.wikipedia.org/wiki/Tangent_lines_to_circles#With_analytic_geometry
        #xarccenter = xintersect + eradius*np.sin(a)
        #yarccenter = yintersect - eradius*np.cos(a)
        #print( xarccenter, yarccenter )
        #
        ##print( xarccenter-xintersect, yarccenter-yintersect )
        #
        ##print(np.arctan2(-eradius**2/g ,  eradius/g*np.sqrt(g*g-eradius**2)) *180/np.pi )
        #sta = np.arctan2(-eradius**2/g ,  eradius/g*np.sqrt(g*g-eradius**2)) *180/np.pi
        #
        ##circle = Arc([xarccenter , yarccenter ],eradius*2,eradius*2, +90, -50,180-sta, color=t[0].get_color(),fill=None,lw=1)
        #circle = Circle([xarccenter , yarccenter ],eradius, color=t[0].get_color(),fill=None,lw=1)
        #plt.gca().add_patch(circle)


        
        
    
    plt.axis('off')
    
    #plt.tight_layout()

    if len(sys.argv) == 2:
        plt.savefig('trajectory.pdf', format='pdf', bbox_inches='tight')
    else:
        plt.show()
               
if __name__=="__main__":
    sys.exit(main(sys.argv))
    
