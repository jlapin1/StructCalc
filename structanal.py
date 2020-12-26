#! /usr/bin/python3

import numpy as np
import os
import sys
import copy
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

pi = np.pi
d2r = pi/180
i = np.complex(0,1)
exp = np.exp
norm = np.linalg.norm
Bo = np.array([0,0,1])
sin = np.sin
cos = np.cos
asin = np.arcsin
acos = np.arccos
sys.path.append("/home/nevlab/Pyrosetta")

class env():
    def __init__(self):
        self.inner2 = np.matrix([
                [0.25, 0, -0.75],
                [0, -0.5, 0],
                [-0.75, 0, 0.25]],'complex')
        self.readang()
        self.readconf()

    def readang(self):
        with open("output/angles.csv", "r") as f:
            self.So = float(f.readline().split(",")[-1])
            self.HNCa = d2r*float(f.readline().split(",")[-1])
            self.NCCa = d2r*float(f.readline().split(",")[-1])
            self.HNCo = d2r*float(f.readline().split(",")[-1])
            self.a_nc = d2r*float(f.readline().split(",")[-1])
            self.tetra = d2r*float(f.readline().split(",")[-1])
            self.tetid = d2r*float(f.readline().split(",")[-1])
            self.omN = [d2r*float(m) for m in f.readline().split(",")[1:]]
            self.omH = [d2r*float(m) for m in f.readline().split(",")[1:]]
            self.omHN = [d2r*float(m) for m in f.readline().split(",")[1:]]
            self.chi0 = float(f.readline().split(",")[-1])
            self.chi1 = float(f.readline().split(",")[-1])
        self.gamma1 = self.omN[0]
        self.angle1 = (3*pi/2)-self.HNCa
        self.angle2 = (3*pi/2)-self.NCCa-self.HNCo
        self.sixty = 60*d2r
    
    def readconf(self):
        with open("config.conf", "r") as f:
            for line in f:
                l = line.split()
                if l[0]=='kb_max':
                    self.kb_max = int(l[1])
                elif l[0]=='CH':
                    self.CH = int(l[1])
                elif l[0]=='car':
                    self.car = float(l[1])
                elif l[0]=='wt':
                    self.wt = [float(m) for m in l[1].split(",")]
    
    def geninner(self, sig1, sig2, sig3):
        sig = [sig1,sig2,sig3]
        Dcs = wigner(self.omN[0], self.omN[1], self.omN[2])
        return self.car*Dcs*M_csa(sig)*dagger(Dcs)

    
class robj():
    def __init__(self, E):
        self.readoutput()
        
        tick=0
        for m in os.listdir("raw"):
            if m[-3:]=='iso':
                self.readiso(m)
                tick+=1
                break
            elif m[-3:]=='inp':
                self.readraw(m)
        if tick==0:
            self.iso = np.zeros((self.sz))
        
        self.readtarg()
        self.readsig()
        self.genspec(E)
    
    def readoutput(self):
        with open("output/output.csv", "r") as f:
            f.readline()
            [self.sz, self.amt] = [int(m) for m in f.readline().split(",")]
            self.res = f.readline()
            self.SA = np.zeros((self.amt, 2))
            self.pp = np.zeros((self.amt, 2, self.sz-1))
            for m in range(self.amt):
                self.SA[m] = [float(n) for n in f.readline().split(",")]
                for n in range(2):
                    self.pp[m][n] = [float(o) for o in f.readline().split(",")]
    
    def readiso(self, fn):
        with open("raw/"+fn,"r") as f:
            f.readline()
            sz = int(f.readline())
            self.iso = np.zeros((sz))
            for m in range(sz):
                self.iso[m] = float(f.readline())
    
    def readraw(self, fn):
        with open("raw/"+fn, "r") as f:
            f.readline()
            sz = int(f.readline())
            self.rawtarg = np.zeros((sz, 3))
            for m in range(sz):
                self.rawtarg[m] = [float(n) for n in f.readline().split(",")]
    
    def readtarg(self):
        with open("input/realtarg.csv", "r") as f:
            f.readline()
            sz = int(f.readline())
            self.targ = np.zeros((sz,3))
            for m in range(sz):
                self.targ[m] = [float(n) for n in f.readline().split(",")[1:]]
    
    def readsig(self):
        with open("output/sigout.csv", "r") as f:
            self.sigs = np.zeros((self.sz, 3))
            for m in range(self.sz):
                self.sigs[m] = [float(n) for n in f.readline().split(",")]
    
    def genspec(self, E):
        #self.y = np.zeros((self.amt, self.sz, 3), dtype='complex')
        self.rafreqs = np.zeros((self.amt, self.sz, 3))
        self.refreqs = np.zeros((self.amt, self.sz, 3))
        for m in range(self.amt):
            
            # First resonance
            y = Y(self.SA[m][0], self.SA[m][1])
            inner1 = E.geninner(self.sigs[0,0], self.sigs[0,1], self.sigs[0,2])
            self.rafreqs[m][0][0] = float(np.real(y*inner1*dagger(y))[0,0])
            self.rafreqs[m][0][1] = E.chi0*np.abs(float(np.real((y*E.inner2*dagger(y))[0,0]))) 
            
            for n in range(1,self.sz,1):
                PHI = self.pp[m][0][n-1]
                PSI = self.pp[m][1][n-1]
                
                #1H13C
                hp = 0
                if self.res[n-1]=='G':
                    ychb = y*wigner(E.angle1, PHI+E.sixty, (pi/2)-E.tetid)*wigner(0, -pi/2, 0)
                    hp = E.chi1*abs(float(np.real((3*(ychb[0,1]**2)-1)/2)))
                ych = y*wigner(E.angle1, PHI-E.sixty, (pi/2)-E.tetid)*wigner(0, -pi/2, 0)
                hn = E.chi1*abs(float(np.real((3*(ych[0,1]**2)-1)/2)))
                self.rafreqs[m][n-1][2] = (hn**2 + hp**2)**0.5
                
                # Propagate y
                y = y*wigner(E.angle1, PHI, E.tetra)*wigner(0, -PSI-pi, E.angle2)
                
                # 15N
                inner1 = E.geninner(self.sigs[n][0], self.sigs[n][1], self.sigs[n][2])
                self.rafreqs[m][n][0] = float(np.real((y*inner1*dagger(y))[0,0]))
                
                # 1H15N
                self.rafreqs[m][n][1] = E.chi0*np.abs(float(np.real((y*E.inner2*dagger(y))[0,0])))
        self.regenspec(E)
    
    def regenspec(self, E):
        self.refreqs = np.zeros((self.amt, self.sz, 3))
        for m in range(self.amt):
            for n in range(self.sz):
                self.refreqs[m][n][0] = (((self.rafreqs[m][n][0]/E.car)-self.iso[n])*(-E.So/2))+self.iso[n]
            self.refreqs[m][:,1] = self.rafreqs[m][:,1]*(E.So/2)
            self.refreqs[m][:,2] = self.rafreqs[m][:,2]*(E.So/2)
    
    def specplot(self, ind, CH=0, block=False):
        self.refreqs[ind][-1,2] = self.rawtarg[-1,2] # free variable
        fig,ax = plt.subplots()
        ax = plt.subplot(111)
        
        ax.set_xlabel('$^{15}N$ Chemical shift (Hz)')
        ax.set_ylabel('$^{1}H^{15}N$ Dipolar coupling (Hz)', color='red')
        ax.tick_params(axis='x', which='major', labelsize=7)
        ax.tick_params(axis='y', which='major', labelsize=7)
        
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.invert_xaxis()
        
        ax.plot(self.rawtarg[:,0], self.rawtarg[:,1], 'ro', markerfacecolor='none')
        ax.plot(self.refreqs[ind][:,0], self.refreqs[ind][:,1], 'rx')
        
        for m in range(self.sz):
            ax.text(self.rawtarg[m,0], self.rawtarg[m,1], '%c%d'%(self.res[m], m), fontsize=6)
            ax.text(self.refreqs[ind][m,0], self.refreqs[ind][m,1], '%c%d'%(self.res[m], m), fontsize=6)
        
        csrms = (sum((self.rawtarg[:,0]-self.refreqs[ind][:,0])**2)/self.sz)**0.5
        csrmshz = (sum((self.targ[:,0]-self.rafreqs[ind][:,0])**2)/self.sz)**0.5
        dcrms = (sum((self.rawtarg[:,1]-self.refreqs[ind][:,1])**2)/self.sz)**0.5
        chrms = 0
        
        if CH==1:
            ax2 = ax.twinx()
            ax2.set_ylabel(r'$^{1}H_{\alpha}^{13}C_{\alpha}$ Dipolar coupling (Hz)', color='blue')
            
            ax2.plot(self.rawtarg[:,0], self.rawtarg[:,2], 'bo', markerfacecolor='none')
            ax2.plot(self.refreqs[ind][:,0], self.refreqs[ind][:,2], 'bx')
            
            ax2.tick_params(axis='y', which='major', labelsize=7)
            ax2.yaxis.set_minor_locator(MultipleLocator(100))
            
            ax2.set_ylim((min(ax.get_ylim()+ax2.get_ylim()),max(ax.get_ylim()+ax2.get_ylim())))
            ax.set_ylim((min(ax.get_ylim()+ax2.get_ylim()),max(ax.get_ylim()+ax2.get_ylim())))
            
            for m in range(self.sz):
                ax.text(self.rawtarg[m,0], self.rawtarg[m,2], '%c%d'%(self.res[m], m), fontsize=6)
            
            chrms = (sum((self.rawtarg[:,2]-self.refreqs[ind][:,2])**2)/self.sz)**0.5
        
        fig.text(0.7,0.9, "$^{15}N_{rms}$ = %.3f ppm (%.0f Hz)\n$^{1}H^{15}N_{rms}$ = %.0f Hz\n$^{1}H^{13}C_{rms}$ = %.0f Hz"%(csrms, csrmshz, dcrms, chrms), size=6)
        
        plt.show(block=block)
    
    def ramachandran(self, ind, block=False):
        fig,ax = plt.subplots()
        ax.set_xlim([-180,180])
        ax.set_ylim([-180,180])
        ax.set_xlabel('$\phi$')
        ax.set_ylabel('$\psi$')
        ax.plot(self.pp[ind][0]/d2r, self.pp[ind][1]/d2r, 'rx')
        minors = [-150, -120, -90, -60, -30, 30, 60, 90, 120, 150] # No tick marks at 0
        ax.set_xticks(minors, minor=True)
        ax.set_xticklabels(['','-120','','-60','','','60','','120',''], minor=True)
        ax.set_yticks(minors, minor=True)
        ax.set_yticklabels(['','-120','','-60','','','60','','120',''], minor=True)
        ax.grid(which='minor',color='gray')
        ax.set_xticks([0]) # major tick marks
        ax.set_yticks([0])
        ax.grid(which='major',color='k')
        ax.set_aspect('equal', adjustable='box')
        for m in range(self.sz-1):
            ax.text(self.pp[ind][0,m]/d2r,self.pp[ind][1,m]/d2r,'%d'%(m),color='k',fontsize=6)
        
        plt.show(block=block)

class pose():
    def __init__(self, beta=0, alpha=0):
        self.sz = 1
        self.xyz = np.zeros((1, 7, 3))
        self.ba = np.array((beta, alpha))
        self.ram = np.zeros((1,2))
        
        y = Y(beta, alpha)
        self.y = y
        X = norm(y)*sin(beta)*cos(alpha)
        why = norm(y)*sin(beta)*sin(alpha)
        Z = norm(y)*cos(beta)
        pn = np.array((X, why, Z))
        self.pn = np.zeros((1,3))
        self.pn[0] = pn
        
        hold = np.cross(pn,Bo) # orthogonal to Bo and plane normal
        hold = hold/norm(hold)
        proj = np.cross(hold,pn) # projection of Bo in peptide plane
        proj = proj/norm(proj)
        self.proj = proj
        
        R = rotaxis(alpha,pn) # Rotate back to MF x-axis -> NH bond
        H = 0.99*np.dot(proj,R)
        R = rotaxis((360-151.8-90)*d2r,pn)
        Ca = 1.47*(np.dot(H,R)/norm(H))
        R = rotaxis(-119.5*d2r,pn)
        CO = 1.47*(np.dot(H,R)/norm(H))
        R = rotaxis(-28.07*d2r,pn)
        O = 2.5486*(np.dot(CO,R)/norm(CO))
        R = rotaxis(33.04*d2r,pn)
        Ca0 = 2.547*(np.dot(CO,R)/norm(CO))
        self.xyz[0][5,:] = H
        self.xyz[0][1,:] = [0,0,0]
        self.xyz[0][2,:] = Ca
        self.xyz[0][0,:] = CO
        self.xyz[0][3,:] = O
        self.xyz[0][4,:] = Ca0
        
        #freqs = measurables(y,y,0)
        #self.freqs = np.array(freqs)
        
        self.I = 0
        
    def addres(self, fi=0, sigh=0):
        self.sz += 1
        I = self.I
        self.I += 1
        add2 = np.zeros((1,7,3))
        self.xyz = np.vstack((self.xyz, add2))
        self.ram[I] = np.array((fi, sigh))
        self.ram = np.vstack((self.ram, np.zeros((1,2))))
        
        pn = self.pn[I]
        N   = self.xyz[I][1,:]
        CO  = self.xyz[I][0,:]-N
        O   = self.xyz[I][3,:]-N
        Ca  = self.xyz[I][2,:]-N
        Ca0 = self.xyz[I][4,:]-N
        H   = self.xyz[I][5,:]-N
        
        ra = 180-111-acos(np.dot(CO-Ca0,Ca)/(norm(CO-Ca0)*norm(Ca)))/d2r
        R = rotaxis(ra*d2r,pn) # rotate so that >N-Ca-CO2 = 111
        trans = Ca-np.dot(Ca0,R) # translate so that planes share Ca
        
        Ca0 = np.dot(Ca0,R)+trans
        CO1 = np.dot(CO, R)+trans
        O1  = np.dot(O, R)+trans
        N1  = trans
        Ca1 = np.dot(Ca,R)+trans
        H1  = np.dot(H, R)+trans
        # Create an Ha between two planes
        R = rotaxis(-27.4*d2r, pn)
        Ha = 2.028*(np.dot(Ca, R) / norm(Ca))
        R = rotaxis(60*d2r, -Ca)
        Ha = np.dot(Ha, R)
        
        # Rotate about bonds from origin, then translate back to their rightful place 
        # Apply phi -> 6 atoms
        R = rotaxis(fi, -Ca)
        Ha = np.dot(Ha, R)
        CO1 = np.dot(CO1, R)
        O1 = np.dot(O1, R)
        N1 = np.dot(N1, R)
        Ca1 = np.dot(Ca1, R)
        H1 = np.dot(H1, R)
        # Apply psi -> 4 atoms
        R = rotaxis(pi+sigh, Ca-CO1) # extra pi because original coordinates were in psi=180
        O1 = np.dot(O1-CO1, R) + CO1
        N1 = np.dot(N1-CO1, R) + CO1
        Ca1 = np.dot(Ca1-CO1, R) + CO1
        H1 = np.dot(H1-CO1, R) + CO1
        
        O1 += N
        N1 += N
        CO1 += N
        Ca1 += N
        Ca0 += N
        H1 += N
        Ha += N
        
        self.xyz[I+1][1] = N1
        self.xyz[I+1][0] = CO1
        self.xyz[I+1][3] = O1
        self.xyz[I+1][2] = Ca1
        self.xyz[I+1][4] = Ca0
        self.xyz[I+1][5] = H1
        self.xyz[I][6] = Ha
        NH = H1-N1
        
        pn1 = np.cross(Ca1-N1, NH)
        pn1 = pn1 / norm(pn1)
        self.pn = np.vstack((self.pn, pn1))
        hold = np.cross(pn1,Bo) # orthogonal to Bo and plane normal
        hold = hold/norm(hold)
        proj = np.cross(hold,pn1) # projection of Bo in peptide plane
        proj = proj/norm(proj)
        self.proj = np.vstack((self.proj, proj))
        
        beta = acos(np.dot(pn1, Bo)/(norm(pn1)*norm(Bo)))
        alpha = acos(np.dot(proj, NH) / (norm(proj)*norm(NH)))
        perp = np.cross(NH, proj)
        perp = perp / norm(perp)
        inorout = np.dot(perp, pn1)
        if round(inorout) == -1:
            alpha *= -1
        y = Y(beta, alpha)
        #freqs = measurables(self.y[I], y, fi)
        
        self.ba = np.vstack((self.ba, np.array((beta, alpha))))
        self.y = np.vstack((self.y, y))
        #self.freqs = np.vstack((self.freqs, np.array(freqs)))
    
    def delres(self):
        I = self.I
        self.I -= 1
        self.sz -= 1
        self.xyz = self.xyz[:I]
        self.ram = self.ram[:I]
        self.pn = self.pn[:I]
        self.proj = self.proj[:I]
        self.ba = self.ba[:I]
        self.y = self.y[:I]
        #self.freqs = self.freqs[:I]
    
    def plotfull(self, block=False):
        # https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
        Xmax = self.xyz[:,:,0].max()
        Xmin = self.xyz[:,:,0].min()
        Ymax = self.xyz[:,:,1].max()
        Ymin = self.xyz[:,:,1].min()
        Zmax = self.xyz[:,:,2].max()
        Zmin = self.xyz[:,:,2].min()
        max_range = np.array([Xmax-Xmin, Ymax-Ymin, Zmax-Zmin]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(Xmax+Xmin)
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Ymax+Ymin)
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Zmax+Zmin)
        plt.figure()
        ax = plt.subplot(111, projection='3d')
        #ax.set_aspect('equal')
        ax.set_xlabel('x-axis')
        ax.set_ylabel('y-axis')
        ax.set_zlabel('z-axis')
        for i in range(self.sz):
            H   = self.xyz[i][5]
            Ca  = self.xyz[i][2]
            N   = self.xyz[i][1]
            CO  = self.xyz[i][0]
            O   = self.xyz[i][3]
            Ca0 = self.xyz[i][4]
            if i>0:
                Ha  = self.xyz[i-1][6]
            NH  = H-N
            proj = copy.copy(self.proj[i])
            pn = copy.copy(self.pn[i])
            ax.plot([H[0]],[H[1]],[H[2]],color='gray',marker='.',markersize=4)
            ax.plot([Ca[0]],[Ca[1]],[Ca[2]],color='lightblue',marker='.',markersize=6)
            ax.plot([N[0]],[N[1]],[N[2]],color='blue',marker='.',markersize=6)
            ax.plot([CO[0]],[CO[1]],[CO[2]],color='lightblue',marker='.',markersize=6)
            ax.plot([O[0]],[O[1]],[O[2]],color='red',marker='.',markersize=6)
            ax.plot([Ca0[0]],[Ca0[1]],[Ca0[2]],color='lightblue',marker='.',markersize=6)
            if i>0:
                ax.plot([Ha[0]],[Ha[1]],[Ha[2]],color='gray',marker='.',markersize=2)
            
            ax.plot([N[0],H[0]],[N[1],H[1]],[N[2],H[2]],'gray', linewidth=0.5)
            ax.plot([N[0],Ca[0]],[N[1],Ca[1]],[N[2],Ca[2]],'gray', linewidth=0.5)
            ax.plot([N[0],CO[0]],[N[1],CO[1]],[N[2],CO[2]],'gray', linewidth=0.5)
            ax.plot([CO[0],O[0]],[CO[1],O[1]],[CO[2],O[2]],'gray', linewidth=0.5)
            ax.plot([CO[0],Ca0[0]],[CO[1],Ca0[1]],[CO[2],Ca0[2]],'gray', linewidth=0.5)
            if i>0:
                ax.plot([Ca0[0],Ha[0]],[Ca0[1],Ha[1]],[Ca0[2],Ha[2]],'gray')
            
            ax.plot([H[0],Ca[0]],[H[1],Ca[1]],[H[2],Ca[2]],'r:', linewidth=1)
            ax.plot([Ca[0],O[0]],[Ca[1],O[1]],[Ca[2],O[2]],'r:', linewidth=1)
            ax.plot([O[0],Ca0[0]],[O[1],Ca0[1]],[O[2],Ca0[2]],'r:', linewidth=1)
            ax.plot([Ca0[0],H[0]],[Ca0[1],H[1]],[Ca0[2],H[2]],'r:', linewidth=1)
            
            # Plane normal arrow (red)
            if i==0:
                mul = 1.5
            else:
                mul = 1
            z_arr = Arrow3D([N[0],N[0]+mul*pn[0]], [N[1],N[1]+mul*pn[1]], [N[2],N[2]+mul*pn[2]], mutation_scale=7, lw=1, arrowstyle="-|>", color="r")
            ax.add_artist(z_arr)
            # Bo arrow (black)
            z_arr2 = Arrow3D([N[0],N[0]], [N[1],N[1]], [N[2],N[2]+1.5], mutation_scale=7, lw=1, arrowstyle="-|>", color="k")
            ax.add_artist(z_arr2)
            
            # Plot Bo projection
            scale = 1.5
            ax.plot([N[0], N[0]+(scale*proj[0])],[N[1],N[1]+(scale*proj[1])],[N[2],N[2]+(scale*proj[2])],'k:')
            
            # phi/psi labels
#            ax.text(0.5*(Ca[0]-N[0])+N[0], 0.5*(Ca[1]-N[1])+N[1], 0.5*(Ca[2]-N[2])+N[2], r'$\phi_{%d}$'%(i+1), size=7)
#            ax.text(0.5*(CO[0]-Ca0[0])+Ca0[0], 0.5*(CO[1]-Ca0[1])+Ca0[1], 0.5*(CO[2]-Ca0[2])+Ca0[2], r'$\psi_{%d}$'%(i), size=7)
            
            # beta/alpha marker
            if i==0:
                ax.plot([0.55*Bo[0], 0.6*pn[0]], [0.55*Bo[1], 0.6*pn[1]], [0.55*Bo[2], 0.6*pn[2]],'g')
                ax.plot([N[0]+(0.5*(H[0]-N[0])),N[0]+(0.5*proj[0])],[N[1]+(0.5*(H[1]-N[1])),N[1]+(0.5*proj[1])],[N[2]+(0.5*(H[2]-N[2])),N[2]+(0.5*proj[2])],'g')
            
#            if i==0:
#                ax.text(-.2+0.4*(proj[0] + NH[0]) + N[0], -2 + 0.4*(proj[1] + NH[1]) + N[1], 3+ -.2 + 0.4*(proj[2] + NH[2]) + N[2],r'$\alpha_{1}$', size=10, color='black')
#                alpharr = Arrow3D([-.2+0.4*(proj[0] + NH[0]) + N[0],-.2+0.4*(proj[0] + NH[0]) + N[0]], [-1.5 + 0.4*(proj[1] + NH[1]) + N[1], 0.2+ 0.4*(proj[1] + NH[1]) + N[1]], [3+ -.2 + 0.4*(proj[2] + NH[2]) + N[2], -.2 + 0.4*(proj[2] + NH[2]) + N[2]], mutation_scale=7, lw=0.5, arrowstyle="-|>", color="gray")
#                ax.add_artist(alpharr)
#                ax.text( 0.8*(pn[0] + Bo[0]) + N[0], 1+ 0.8*(pn[1] + Bo[1]) + N[1], 2+ 0.8*(pn[2] + Bo[2]) + N[2], r'$\beta_{1}$', size=10, color='black')
#                boarr = Arrow3D([0.8*(pn[0] + Bo[0]) + N[0],0.8*(pn[0] + Bo[0]) + N[0]], [1+ 0.8*(pn[1] + Bo[1]) + N[1], -0.3+ 0.8*(pn[1] + Bo[1]) + N[1]], [2+ 0.8*(pn[2] + Bo[2]) + N[2],-0.8+ 0.8*(pn[2] + Bo[2]) + N[2]], mutation_scale=7, lw=0.5, arrowstyle="-|>", color="gray")
#                ax.add_artist(boarr)
                
        # Create bounding box to scale all distances correctly
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')
        
        plt.show(block=block)
    
    def writepdb(self, fn, rob):
        index = 0
        atomn = 1
        resn = 1
        fp = 'analysis/'
        ext = '.pdb'
        f = open(fp+fn+ext,'w')
        
        # Position NH2 terminus
        # Find normal to C-O and Ca-C
        CO = self.xyz[0][3]-self.xyz[0][0]
        CO = CO / norm(CO)
        CaC = self.xyz[0][4]-self.xyz[0][0]
        CaC = CaC / norm(CaC)
        n = np.cross(CaC,CO)
        # Rotate CCa bond to NCa
        R = rotaxis(116.2*d2r,n)
        NCa = 1.47*np.dot(-CaC,R)
        
        # Plot N(H2) by adding vector to first Ca atom
        f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                'ATOM',atomn,'N', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][4][0]+NCa[0], 
                self.xyz[index][4][1]+NCa[1], self.xyz[index][4][2]+NCa[2], 1, 0, 'N'))
        
        atomn += 1
        
        f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                'ATOM',atomn,'CA', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][4][0], 
                self.xyz[index][4][1], self.xyz[index][4][2], 1, 0, 'C'))
        
        atomn += 1
        
        for i in range(rob.sz-1):
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'C', AA[rob.res[i]], 'A', resn, self.xyz[index][0][0], 
                    self.xyz[index][0][1], self.xyz[index][0][2], 1, 0, 'C'))
            
            atomn += 1
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'O', AA[rob.res[i]], 'A', resn, self.xyz[index][3][0], 
                    self.xyz[index][3][1], self.xyz[index][3][2], 1, 0, 'O'))
            
            atomn += 1
            resn += 1
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'N', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][1][0], 
                    self.xyz[index][1][1], self.xyz[index][1][2], 1, 0, 'N'))
            
            atomn += 1
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'CA', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][2][0], 
                    self.xyz[index][2][1], self.xyz[index][2][2], 1, 0, 'C'))
            
            atomn += 1
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'H', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][5][0], 
                    self.xyz[index][5][1], self.xyz[index][5][2], 1, 0, 'H'))
            
            atomn += 1
            
            f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'HA', AA[rob.res[resn-1]], 'A', resn, self.xyz[index][6][0], 
                    self.xyz[index][6][1], self.xyz[index][6][2], 1, 0, 'H'))
            
            atomn += 1
            index += 1
        
        # Must add last C and O
        NH = self.xyz[index-1][5]-self.xyz[index-1][1]
        NH = NH / norm(NH)
        CaN = self.xyz[index-1][2] - self.xyz[index-1][1]
        CaN = CaN / norm(CaN)
        n = np.cross(CaN,NH)
        R = rotaxis(111*d2r,-n)
        CaC = 1.54*np.dot(-CaN,R)
        
        f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                    'ATOM',atomn,'C', AA[rob.res[resn-1]], 'A', resn, self.xyz[index-1][2][0]+CaC[0], 
                    self.xyz[index-1][2][1]+CaC[1], self.xyz[index-1][2][2]+CaC[2], 1, 0, 'C'))
            
        atomn += 1
        
        R = rotaxis(120.8*d2r,-n)
        CaC2 = CaC / norm(CaC)
        CO = 1.43*np.dot(CaC2,R)
        
        f.write("%4s %6d %3s %4s %c %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11c\n"%(
                'ATOM',atomn,'O', AA[rob.res[resn-1]], 'A', resn, self.xyz[index-1][2][0]+CaC[0]+CO[0], 
                self.xyz[index-1][2][1]+CaC[1]+CO[1], self.xyz[index-1][2][2]+CaC[2]+CO[2], 1, 0, 'O'))
        
        f.write('TER')
        f.close()
        
class Arrow3D(FancyArrowPatch):
    
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions( (xs[0], ys[0]), (xs[1], ys[1]) )
        FancyArrowPatch.draw(self, renderer)

def genrosterms(rob):
    import pyrosetta
    from pyrosetta.rosetta.protocols.membrane import AddMembraneMover,MembranePositionFromTopologyMover
    from pyrosetta import teaching as t
    pyrosetta.init()
    #XYZ = pyrosetta.rosetta.numeric.xyzVector_double_t
    
    # switching functions
    to_centroid = pyrosetta.SwitchResidueTypeSetMover('centroid')
    recover_sidechains = pyrosetta.SwitchResidueTypeSetMover('fa_standard')
    
    # Add Membrane Representation
    add_memb = AddMembraneMover()
    tick = 0
    for m in os.listdir("raw/"):
        if m[-4:]=='span':
            add_memb.spanfile("raw/"+m)
            tick+=1
            break
    if tick==0:
        print("\nError: No span file found in raw/ directory")
        return
    
    # Embed protein in the membrane
    init_mem_pos = MembranePositionFromTopologyMover()
    
    # Create score function
    # C:\Users\joell\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\usr\local\lib\python2.7\dist-packages\pyrosetta-2018.28+release.769b677-py2.7-linux-x86_64.egg\pyrosetta\database\scoring\weights
    #   - pyrosetta.rosetta.core.scoring
    customcen = pyrosetta.ScoreFunction()
    customcen.name("Customcen")
    customcen.set_weight(t.MPEnv, 1.0);customcen.set_weight(t.MPCbeta, 1.0);customcen.set_weight(t.MPPair, 1.0);customcen.set_weight(t.MPTermini, 1.0)
    customcen.set_weight(t.MPNonHelix, 1.0);customcen.set_weight(t.MPTMProj, 1.0);customcen.set_weight(t.MPHelicality, 1.0)
    customfa = pyrosetta.ScoreFunction()
    customfa.name("Customfa")
    customfa.set_weight(t.FaMPSolv, 1.0);customfa.set_weight(t.FaMPEnv, 1.0);customfa.set_weight(t.FaMPEnvSmooth, 1.0);customfa.set_weight(t.fa_atr, 1.0);customfa.set_weight(t.fa_rep, 1.0)
    customfa.set_weight(t.fa_sol, 1.0);customfa.set_weight(t.fa_intra_rep, 1.0);customfa.set_weight(t.fa_elec, 1.0);customfa.set_weight(t.pro_close, 1.0);customfa.set_weight(t.hbond_sr_bb, 1.0)
    customfa.set_weight(t.hbond_lr_bb, 1.0);customfa.set_weight(t.hbond_bb_sc, 1.0);customfa.set_weight(t.hbond_sc, 1.0);customfa.set_weight(t.dslf_fa13, 1.0);customfa.set_weight(t.rama, 1.0)
    customfa.set_weight(t.omega, 1.0);customfa.set_weight(t.fa_dun, 1.0);customfa.set_weight(t.p_aa_pp, 1.0);customfa.set_weight(t.ref, 1.0)
    
    f = open("analysis/termscen","w")
    g = open("analysis/termsfa", "w")
    namescen = "MPEnv MPCbeta MPPair MPTermini MPNonHelix MPTMProj MPHelicality"
    namesfa = "FaMPSolv FaMPEnv FaMPEnvSmooth fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbonb_sc dslf_fa13 rama omega fa_dun p_aa_pp ref"
    f.write(namescen+"\n")
    g.write(namesfa+"\n")
    termscen = np.zeros((rob.amt, len(namescen.split())))
    termsfa = np.zeros((rob.amt, len(namesfa.split())))
    for m in range(rob.amt):
        
        pose = pyrosetta.pose_from_sequence(rob.res)
        
        for n in range(rob.sz-1):
            pose.set_phi(n+2,rob.pp[m][0][n]/d2r)
            pose.set_psi(n+2,rob.pp[m][1][n]/d2r)

        # Apply membrane
        add_memb.apply(pose)
        init_mem_pos.apply(pose)
        
        # This code might be useful one day
        #mem = pose.membrane_info()
        #conf = pose.conformation()
        # pose.xyz(AtomID(x,y))[0]
        # AtomID(x,y) -> (atom#, residue#)
        # residue# of membrane: pose.size() + 1
        #center = mem.membrane_center(conf)
        #center = pose.xyz(AtomID(2,pose.size()+1)) 
        #norm = mem.membrane_normal(conf)
        #print("%f,%f,%f"%(center[0],center[1],center[2]))

        # Score full atom pose
        Efa = customfa.score(pose)
        fampsolv = customfa.score_by_scoretype(pose,t.FaMPSolv);fampenv = customfa.score_by_scoretype(pose,t.FaMPEnv);fampenvsmooth = customfa.score_by_scoretype(pose,t.FaMPEnvSmooth)
        faatr = customfa.score_by_scoretype(pose,t.fa_atr);farep = customfa.score_by_scoretype(pose,t.fa_rep);fasol = customfa.score_by_scoretype(pose, t.fa_sol)
        faintrarep = customfa.score_by_scoretype(pose,t.fa_intra_rep);faelec = customfa.score_by_scoretype(pose,t.fa_elec);proclose = customfa.score_by_scoretype(pose,t.pro_close);hbondsrbb = customfa.score_by_scoretype(pose,t.hbond_sr_bb)
        hbondlrbb = customfa.score_by_scoretype(pose,t.hbond_lr_bb);hbondbbsc = customfa.score_by_scoretype(pose,t.hbond_bb_sc);hbondsc = customfa.score_by_scoretype(pose,t.hbond_sc);dslffa13 = customfa.score_by_scoretype(pose,t.dslf_fa13)
        RAMA = customfa.score_by_scoretype(pose,t.rama);OMEGA = customfa.score_by_scoretype(pose,t.omega);fadun = customfa.score_by_scoretype(pose,t.fa_dun);paapp = customfa.score_by_scoretype(pose,t.p_aa_pp);REF = customfa.score_by_scoretype(pose,t.ref)
        # Switch to centroid and score
        to_centroid.apply(pose)
        Ecen = customcen.score(pose)
        mpenv = customcen.score_by_scoretype(pose, t.MPEnv);mpcbeta = customcen.score_by_scoretype(pose, t.MPCbeta);mppair = customcen.score_by_scoretype(pose, t.MPPair)
        mptermini = customcen.score_by_scoretype(pose, t.MPTermini);mpnonhelix = customcen.score_by_scoretype(pose, t.MPNonHelix);mptmproj = customcen.score_by_scoretype(pose, t.MPTMProj);mphelicality = customcen.score_by_scoretype(pose, t.MPHelicality)
        
        linecen = [mpenv, mpcbeta, mppair, mptermini, mpnonhelix, mptmproj, mphelicality]
        linefa = [fampsolv,fampenv,fampenvsmooth,faatr,farep,fasol,faintrarep,faelec,proclose,hbondsrbb,hbondlrbb,hbondbbsc,hbondsc,dslffa13,RAMA,OMEGA,fadun,paapp,REF]
        
        for n,o in enumerate(linecen):
            termscen[m][n] = o
            if n==(termscen.shape[1]-1):
                f.write("%f\n"%(o))
            else:
                f.write("%f "%(o))
        #f.write("%f %f %f %f %f %f %f\n"%(mpenv, mpcbeta, mppair, mptermini, mpnonhelix, mptmproj, mphelicality))
        for n,o in enumerate(linefa):
            termsfa[m][n] = o
            if n==(termsfa.shape[1]-1):
                g.write("%f\n"%(o))
            else:
                g.write("%f "%(o))
        #g.write("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n"%(fampsolv,fampenv,fampenvsmooth,faatr,farep,fasol,faintrarep,faelec,proclose,hbondsrbb,hbondlrbb,hbondbbsc,hbondsc,dslffa13,RAMA,OMEGA,fadun,paapp,REF))
        #recover_sidechains.apply(pose)
    f.close()
    g.close()
    
    return namescen.split(),termscen,namesfa.split(),termsfa
    
def Y(beta, alpha):
    return np.matrix([-1*(sin(beta)/2**0.5)*exp(i*alpha), cos(beta), (sin(beta)/2**0.5)*exp(-i*alpha)])

def wigner(one, two, three):
    cos_b = cos(two)
    sin_b_sq_2 = sin(two)/2**0.5
    ex_ia = exp(i*one)
    ex_ig = exp(i*three)
    D33 = ex_ia*((1+cos_b)/2)*ex_ig
    D31 = ex_ia*((1-cos_b)/2)*np.conj(ex_ig)
    D23 = sin_b_sq_2*ex_ig
    D32 = ex_ia*sin_b_sq_2
    return np.matrix([
            [np.conj(D33),  -np.conj(D32),  np.conj(D31)],
            [np.conj(D23),          cos_b,          -D23],
            [         D31,            D32,           D33]],'complex')

def M_csa(tsr):
    return np.matrix([
            [(tsr[0]+tsr[1])/2, 0     , (tsr[1]-tsr[0])/2],
            [0                , tsr[2], 0                ],
            [(tsr[1]-tsr[0])/2, 0     , (tsr[0]+tsr[1])/2]],'complex')
    
def dagger(A):
    return np.conj(A.transpose())

def rotaxis(ang, u):
    R = np.zeros((3,3))
    u = u/norm(u)
    cos = np.cos(ang)
    sin = np.sin(ang)
    # https://en.wikipedia.org/wiki/Rotation_matrix
    R[0,0] = cos + (u[0]**2)*(1-cos)
    R[1,1] = cos + (u[1]**2)*(1-cos)
    R[2,2] = cos + (u[2]**2)*(1-cos)
    R[0,1] = u[0]*u[1]*(1-cos) - u[2]*sin
    R[1,0] = u[0]*u[1]*(1-cos) + u[2]*sin
    R[0,2] = u[0]*u[2]*(1-cos) + u[1]*sin
    R[2,0] = u[0]*u[2]*(1-cos) - u[1]*sin
    R[1,2] = u[1]*u[2]*(1-cos) - u[0]*sin
    R[2,1] = u[1]*u[2]*(1-cos) + u[0]*sin
    return R

def printerms():
    sys.stdout.write("%5s%15s%15s%15s\n"%("index","term","mean","stdev"))
    for m in range(50):
        sys.stdout.write("-")
    sys.stdout.write("\n")
    for m,n in enumerate(namescen):
        print("%2d%18s%15.3f%15.3f"%(m,n+'(c)',np.mean(termscen[:,m],axis=0),np.std(termscen[:,m],axis=0)))
    for m,n in enumerate(namesfa):
        print("%2d%18s%15.3f%15.3f"%(m+len(namescen),n+'(f)',np.mean(termsfa[:,m],axis=0),np.std(termsfa[:,m],axis=0)))
    print("(c): centroid\n(f): full atom")

def writetcl(nums):
    with open("analysis/load.tcl", "w") as f:
        f.write("# run this script, in VMD, by typing \"source load.tcl\" in directory with PDB files\n")
        f.write("set arr {%d"%(nums[0]))
        if len(nums)>1:
            for m in nums[1:]:
                f.write(" %d"%(m))
        f.write("}\n")
        f.write("for {set i 0} {$i<%d} {incr i} {\n"%(len(nums)))
        f.write("        mol new PDB[lindex $arr $i].pdb\n")
        f.write("        mol modcolor top top colorID $i\n")
        f.write("        mol modstyle top top NewCartoon\n")
        f.write("}")
    print("\nWrote analysis/load.tcl for viewing structures in VMD\n")

# AA Used in pose.writepdb()
AA = {
      'G':'GLY',
      'A':'ALA',
      'L':'LEU',
      'M':'MET',
      'F':'PHE',
      'W':'TRP',
      'K':'LYS',
      'Q':'GLN',
      'E':'GLU',
      'S':'SER',
      'P':'PRO',
      'V':'VAL',
      'I':'ILE',
      'C':'CYS',
      'Y':'TYR',
      'H':'HIS',
      'R':'ARG',
      'N':'ASN',
      'D':'ASP',
      'T':'THR'
            }

E = env()
S = robj(E)

lev1 = 0
while lev1==0:
    ans1 = int(input("\nChoose an option\n0: Plot output\n1: Rosetta filtering\n2: Exit\n>>> "))
    if ans1==0:
        ans = int(input("\nType output index # (0-%d)\n>>> "%(S.amt-1)))
        
        S.specplot(ans, CH=E.CH)
        S.ramachandran(ans)
        
        P = pose(S.SA[ans][0], S.SA[ans][1])
        for m in range(S.sz-1):
            P.addres(S.pp[ans][0,m], S.pp[ans][1,m])
        P.addres(S.pp[ans][0,m], S.pp[ans][1,m]) # top it off
        P.plotfull(block=False)
        
    elif ans1==1:
        tick=0
        for m in os.listdir("analysis/"):
            if m=="termsfa":
                tick+=1
            elif m=="termscen":
                tick+=1
        print("\n%d Rosetta scoring term file(s) were found in analysis/"%(tick))
        if tick!=2:
            print("Must have both files to proceed with Rosetta filtering")
        ans = input("\nWould you like to generate Rosetta terms? (y/n)\n>>> ")
        if ans.upper()=='Y':
            namescen, termscen, namesfa, termsfa = genrosterms(S)
        else:
            with open("analysis/termscen", "r") as f:
                namescen = f.readline().split()
                termscen = []
                for line in f:
                    termscen.append([float(m) for m in line.split()])
            termscen = np.array(termscen)
            with open("analysis/termsfa", "r") as f:
                namesfa = f.readline().split()
                termsfa = []
                for line in f:
                    termsfa.append([float(m) for m in line.split()])
            termsfa = np.array(termsfa)
        
        inds = np.arange(S.amt)
        termsall = np.column_stack((termscen,termsfa))
        namesall = namescen+namesfa
        log = []
        lev2 = 0
        while lev2==0:
            
            printerms()
            
            ans = input("\nInput indices separated by correspending weights separated by # of top scores\ni1,i2,...in w1,w2,...wn t\n>>> ").split()
            if len(ans)!=3:
                print("\nInput error: Must have 3 space separated entries (counted %d)"%(len(ans)))
                input("Press enter to try again...\n")
                continue
            
            Is = np.array([int(m) for m in ans[0].split(",")])
            ws = np.array([float(m) for m in ans[1].split(",")])
            top = int(ans[2])
            
            if len(Is)!=len(ws):
                print("\nInput error: Index list did not match weight list")
                input("Press enter to try again...\n")
                continue
            else:
                
                Ws = np.zeros(( len(namescen) + len(namesfa) ))
                Ws[Is] = ws
                scores = termsall[inds].dot(Ws)
                sort = np.argsort(scores)[:top]
                inds = inds[sort]
                
                log.append([np.array(namesall)[Is],ws,top,inds])
                
                print(inds)
                
                ans = input("\nKeep filtering? (y/n)\n>>> ")
                if ans.upper()=='Y':
                    continue
                else:
                    ans = input("\nName the rosetta filtering logfile to write to analysis/")
                    if len(ans)==0:
                        ans = "logfile"
                    with open("analysis/%s"%(ans),"w") as f:
                        for m in log:
                            for n in m[0]:
                                f.write("%s "%(n))
                            f.write("\n")
                            for n in m[1]:
                                f.write("%.3f "%(n))
                            f.write("\n")
                            f.write("%d\n"%(m[2]))
                            for n in m[3]:
                                f.write("%d "%(n))
                            f.write("\n")
                    ans = input("\nWrite PDB files for final indices? (y/n)\n>>> ")
                    if ans.upper()=='Y':
                        for m in inds:
                            P = pose(S.SA[m][0], S.SA[m][1])
                            for n in range(S.sz-1):
                                P.addres(S.pp[m][0,n], S.pp[m][1,n])
                            P.addres(S.pp[m][0,n], S.pp[m][1,n]) # top it off
                            P.writepdb("PDB"+str(m), S)
                        writetcl(inds)
                    lev2=1
    else:
        lev1=1
