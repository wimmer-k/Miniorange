#!/usr/bin/env python3
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Circle #, Arc, 

verbose = False # True # 
# constants
ME = 510.998   # keV/c^2
CC = 299792458 # m/s
# wedge
wradius = 50
wangle = 75
#absorber
absx = [-30,0,20]
absy = [5,10,5]
# detector
dradius = 10
dthick = 2

def main(argv):
    g = -70
    b = +40

    BF = 0.150 # T
    EE = 1000 # keV

    
    
    dashes = [8,8]
    
    fig = plt.figure(figsize=(8,4))
    ax = plt.subplot(1,1,1)
    plt.gca().axis([g-10,b+20,-wradius-5,wradius+5])
    ax.set_aspect('equal', 'box')

    # axis
    p, = plt.plot([g-10,b+20],[0,0],color='k',linewidth=1)
    p.set_dashes(dashes)
    for d in range(g,b+20,10):
        plt.plot([d,d],[-1,1],color='k',linewidth=1)
        plt.text(d,-25,d,ha="center",va="center",fontsize =8)
    
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

    # detector
    coords = [[b,dradius], [b+dthick,dradius], [b+dthick,-dradius], [b,-dradius]]
    print(coords)
    detector = Polygon(coords,color="black")
    plt.gca().add_patch(detector)
    
    
    # trajectories for different angles
    v = np.sqrt(1-(1/(EE/ME+1))**2)
    print("Eletron kinetic energy %f keV" % EE )
    print("Eletron velocity beta %f" % v )
    v = (v*CC*1000/1e12)
    print("Eletron velocity %f mm/ps" % v )
    
    Omega = BF/(1000*ME/CC/CC)/1e12
    print("Omega = qB/m %e 1/ps" % Omega )

    print("radius = v/Omega = %f mm" %  (v/Omega) )
    print("constant = Omega/v = %f mm" %  (Omega/v) )
    
    ## time is including velocity!
    angles = True
    energies = False
    if angles:
        for a in range(-50,55,5):
            trajectory(a,EE,BF,g,b)

    if energies:
        a = 30
        for e in range(500,2500,100):
            trajectory(a,e,BF,g,b)
        a = -30
        for e in range(500,2500,100):
            trajectory(a,e,BF,g,b)
        
        
        
    
    plt.axis('off')
    
    #plt.tight_layout()

    if len(sys.argv) == 2:
        plt.savefig('trajectory.pdf', format='pdf', bbox_inches='tight')
    else:
        plt.show()


def trajectory(a,E,BF,g,b):

    v = np.sqrt(1-(1/(E/ME+1))**2)
    if verbose:
        print("Eletron kinetic energy %f keV" % E )
        print("Eletron velocity beta %f" % v )
    v = (v*CC*1000/1e12)
    if verbose:
        print("Eletron velocity %f mm/ps" % v )
    
    Omega = BF/(1000*ME/CC/CC)/1e12
    if verbose:
        print("Omega = qB/m %e 1/ps" % Omega )
        print("radius = v/Omega = %f mm" %  (v/Omega) )
        print("constant = Omega/v = %f mm" %  (Omega/v) )
    
    a = a*np.pi/180
    dt = 1
    tt = np.linspace(0,1000,1000/dt)
    dt *= v
    xx = [g]
    yy = [0]
    t0 = -1
    vx1 = None
    vy1 = None
    for t in tt:
        x = xx[-1]+np.cos(a)*dt
        y = yy[-1]+np.sin(a)*dt
        #print(x,y,-x,-x/np.tan(wangle/2*np.pi/180))
        if x > absx[0] and x < 0 and abs(y) < absy[0]:
            continue
        if t0<0 and abs(y) < -x/np.tan(wangle/2*np.pi/180):
            #print("straight")
            xx.append(x)
            yy.append(y)
        else:
            if t0<0:
                t0=t
            sign = +1
            if y < 0:
                sign = -1
            x = xx[-1] + (np.cos(a)*np.cos(Omega/v*(t-t0)) + np.sin(a)*np.sin(Omega/v*(t-t0)) * sign)*dt/Omega*v
            y = yy[-1] + (np.sin(a)*np.cos(Omega/v*(t-t0)) - np.cos(a)*np.sin(Omega/v*(t-t0)) * sign)*dt/Omega*v 
            #print(x,y,y*np.tan(wangle/2*np.pi/180));
            if abs(x) < abs(y*np.tan(wangle/2*np.pi/180)) and x*x+y*y<wradius*wradius:
                xx.append(x)
                yy.append(y)
            else:
                if vx1 == None or vy1 == None:
                    vx1 = (np.cos(a)*np.cos(Omega/v*(t-t0)) + np.sin(a)*np.sin(Omega/v*(t-t0)) * sign)/Omega*v
                    vy1 = (np.sin(a)*np.cos(Omega/v*(t-t0)) - np.cos(a)*np.sin(Omega/v*(t-t0)) * sign)/Omega*v
                    
                x = xx[-1] + vx1*dt
                y = yy[-1] + vy1*dt
                xx.append(x)
                yy.append(y)
        #print(x,y,b)
        if x > b+20:
            break
        if x > b and abs(y) < dradius:
            break
        
    plt.plot(xx,yy,color="tab:green",lw=1)
        
if __name__=="__main__":
    sys.exit(main(sys.argv))
    
