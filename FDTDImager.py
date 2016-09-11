#!/usr/bin/python
import sys
import math
import timeit
import numpy as np
from PIL import Image
#todo - numpy arrays

class Solver():
    """ FDTD Solver for TM """
    def __init__(self,x_sz =360,y_sz = 360):
        """ Initialize 2d fields arrays """
        self.x_sz = x_sz
        self.y_sz = y_sz
        #allocate 2d arrays for fields
        self.Hx = np.zeros((x_sz+1,y_sz+1))
        self.Hy = np.zeros((x_sz+1,y_sz+1))
        self.Ez = np.zeros((x_sz+1,y_sz+1))

        #allocate 2d arrays for space characteristics


        #free space impedance
        self.imp0 = 377.0
        #courant number
        self.cn = (1.0 / math.sqrt(2))

        #solve step
        self.time = 0


    def _uH(self):
        """ Update H field """
        #space parameters TODO (need to be stored and calculated)
        c1 = 1.0
        c2 = self.cn/self.imp0
        c3 = 1.0
        c4 = self.cn/self.imp0
        # Hx update
        #for i in range(self.x_sz):
        #    for j in range(self.y_sz-1):
        #        self.Hx[i,j] = c1* self.Hx[i,j] - c2*(self.Ez[i,j+1]-self.Ez[i,j])
                # Hy update

        self.Hx[:,:-1] = c1 * self.Hx[:,:-1] - c2*( self.Ez[:,1:] - self.Ez[:,:-1] )
        #for i in range(self.x_sz-1):
        #    for j in range(self.y_sz):
        #        self.Hy[i,j] = c3* self.Hy[i,j] + c4*(self.Ez[i+1,j]-self.Ez[i,j])
        self.Hy[:-1,:] = c3 * self.Hy[:-1,:] + c4 * ( self.Ez[1:,:] - self.Ez[:-1,:] )
    def _uE(self):
        """ Update E field """
        #space parameters TODO
        c5 = 1.0
        c6 = 1.0
        c7 = self.cn*self.imp0

        # Ez update
        #for i in range (1,self.x_sz-1):
        #    for j in range (1,self.y_sz-1):
        #        self.Ez[i,j] = c5*self.Ez[i,j] + c7*((self.Hy[i,j] - self.Hy[i-1,j]) - (self.Hx[i,j] - self.Hx[i,j-1]))

        #numpy vectorization
        self.Ez[1:-1,1:-1] = c5*self.Ez[1:-1,1:-1] + c7*((self.Hy[1:-1,1:-1]-self.Hy[:-2,1:-1]) - (self.Hx[1:-1,1:-1] - self.Hx[1:-1,:-2]))
    def solveStep(self):
        """ Solve single 2D step within gived boundary conditions and source """
        self._uH()
        self._uE()
        self.time += 1
        #TODO Boundary conditions
        #TODO Source
        ppw = 70 #points per wavelength
        a = math.pi* ((self.cn*self.time)/ppw - 1.0)
        a = a*a
        self.Ez[math.floor(self.x_sz/2)][math.floor(self.y_sz/2)] = 10*( 1.0 - 2.0 * a)*math.exp(-a)
        #if self.time < 700:
        #    self.Ez[math.floor(self.x_sz/2)][math.floor(self.y_sz/2)] = math.sin(self.cn*self.time/50)*2

    def solveTime(self,times = 10):
        tmg = [0]*times
        for i in range(times):
            try:
                #preformance measurement for step
                i_t = timeit.default_timer()
                self.solveStep()
                tmg[i] = timeit.default_timer() - i_t
                print("PM:dt=%f,step=%d,xsz=%d,ysz=%d" % (tmg[i],self.time,self.x_sz,self.y_sz))
                #debug dump every 10 step
                #if self.time > 600:
                if self.time % 10 == 0:
                   self.dumpState()
                #except (ValueError):
                #    print ("Simulation reached incorrect values (NaN) at t=%d STOP" % self.time)
            except (OverflowError):
                print ("Simulation reached incorrent values (Too Big) at t=%dSTOP" % self.time)
                return 0
        dt_s = 0.0
        for i in range (times):
            dt_s += tmg[i]
        print ("SOLVED, dt_mean=%f" % (dt_s/i))


    def dumpState(self):
        fname = ("FDTD2D-T%04d.png" % self.time)
        img = Image.new("RGB",(self.x_sz,self.y_sz))
        rast = img.load()
        for i in range (self.x_sz):
            for j in range (self.y_sz):
                #TODO do correct remap
                Ep = int (abs(float(self.Ez[i,j])*600))
                #Hp = int (abs((self.Hx[i][j])*600*377 ))
                try:
                    rast[i, j] = (Ep,0,0)
                except (ValueError, OverflowError):
                    rast[i,j] = (255,255,255)
        img.save(fname)
        print ("Snapshot at time %d" % self.time)



#test code
print("START")
slv = Solver()
slv.solveTime(1000)
#slv.dumpState()
