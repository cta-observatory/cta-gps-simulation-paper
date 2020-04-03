# Define the Supernova class and SNR class and connected procedures
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab as pl
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker
from itertools import zip_longest
import os
import csv
from random import choices
from matplotlib import rc

from astropy.io import ascii
from astropy.constants import c
import astropy.units as u
import naima
from naima.models import (ExponentialCutoffBrokenPowerLaw, Synchrotron,
                          InverseCompton,PionDecay, TableModel)




# CONSTANTS
definition_pevatron=500.  #TeV 
parsec=3.09*pow(10.,18.)
Msol=1.9*pow(10,33.)
E51=pow(10,51.)
masseproton=1.67*pow(10,-24.)
c=3.*pow(10,10) #  // cm/s
mc2_electron_TeV=511.*pow(10, -9.)
mc2_electron_erg=511.*1.602*pow(10, -9.)
k_boltzmann_erg=1.38*pow(10., -16)  #; // erg.
masseproton_TeV= 0.938*pow(10.,-3.)  # TeV
sigma_thomson=6.6524*pow(10., -25.) #; // en cm2
kyear_sec=np.pi*pow(10, 10.)
year=np.pi*10**7.

mu= 1.36
gamma_ad=4./3.
beta=6*(gamma_ad-1)/(gamma_ad+1)
sigma=4.
r00=8.5  # en kiloparsec
distance_SS_GC=8.5 # distance to galactic center in kiloparsec
erg_to_TeV=0.624151
TeV_to_erg=1./erg_to_TeV

age_sample=20. # kyear 20-40 enough to study Pevatrons
age_max_ST_pase=20. #kyears, typical end of ST phase

#--------------------------------------------------------
# GAS  Distribution functions (Shibata et al. 2010)
def rho_H1  ( r): # / r in kiloparsec , rho_H1 in 10^20 H atoms.cm^-3
   return np.exp(- ( 3.862+7.903*pow(10,-1.)*r -9.426*pow(10,-2.)*np.log(r) - 4.261*pow(r,1./2.)    ))

def rho_H2  ( r): #// r in kiloparsec , rho_H1 in 10^20 H atoms.cm^-3
    return np.exp(- (1.848+8.339*pow(10,-1.)*r -5.560*np.log(r) + 2.405*pow(10,-2.)*pow(r,2.)))


def XI_H1 ( r,  z ):
    return 1/(0.065*np.sqrt(np.pi)+0.160)*(0.4*np.exp(-pow(z/0.12,2.))+0.2*np.exp(-pow(z/0.35, 2.))+0.4*np.exp(-np.fabs(z)/0.40))

def z0 (  r):
    return 0.4*np.cosh(2*r/(3*r00))

def XI_H2 ( r, z ):
    return 1/(0.036*np.sqrt(np.pi)+0.2*z0(r))*( np.exp(-pow(z/0.071, 2.))+ 0.2*(np.fabs(z)/z0(r))*np.exp(-np.fabs(z)/z0(r)))

#nH1, nH2
def nH1 (  r,  z): # // r and z are in kilo parsecs
    return 0.57*XI_H1(r, z)/XI_H1(r00, 0)*rho_H1(r)/rho_H1(r00)

def nH2 (  r,  z): # // r and z are in kilo parsecs
    return 0.53*XI_H2(r, z)/XI_H2(r00, 0)*rho_H2(r)/rho_H2(r00)

def nH (  r,  z):
    return nH1(r,z)+2*nH2(r, z)


#--------------------------------------------------------
#preliminary functions for type 1 and type2 dynamics : all dynamics are from Ptuskin Zirahaskvili 2005
def alpha1_def ( Mdot, uw6):
    return masseproton*((Mdot*pow(10, -5)*Msol/(kyear_sec*pow(10, -3.)))/(4*np.pi*1.4*masseproton*uw6*pow(10, 6)))

def rho1_r (self, r ):
    return self.alpha1*pow(r*parsec, -2.)

def rho_r (self, r):
    if (r<self.r1):
        return self.rho1_r(r)
    else:
        if (r<self.r2):
            return self.rho2
        else :
            if (r<self.r2+0.5):
                return    (self.rho0-self.rho2)*(self.r2-self.r1)/0.5
            else :
                return float(self.rho0)

def INT1 ( self,  r ):
    return self.Mej*Msol*pow(r*parsec, self.beta)/beta+4*np.pi*self.alpha1*pow(r*parsec,self.beta+1.)/(self.beta+1.)

def INT2 (self, r ):
    return self.INT1(self.r1)+  (self.Mej*Msol+4*np.pi*self.alpha1*self.r1*parsec-4./3.*np.pi*self.rho2*pow(self.r1*parsec,3.))* (  pow(r*parsec,self.beta)- pow(self.r1*parsec, self.beta)   )/self.beta   +4./3.*np.pi*self.rho2*(    pow(r*parsec,self.beta+3)   -pow(self.r1*parsec, self.beta+3))/(self.beta+3)

def INT0 (self, r):
    return self.INT2(self.r2)+ (self.Mej*Msol+4*np.pi*self.alpha1*self.r1*parsec+4./3.*np.pi*self.rho2*(pow(self.r2*parsec,3.)-pow(self.r1*parsec,3.))-4./3.*np.pi*self.rho0*pow(self.r2*parsec,3))*(pow(r*parsec,self.beta)-pow(self.r2*parsec,self.beta))/self.beta + 4./3.*self.rho0*np.pi*(pow(r*parsec,self.beta+3.)-pow(self.r2*parsec,self.beta+3.))/(self.beta+3.)


def INT (self, r):
    if (r<self.r1):
        return self.INT1(r)
    else:
        if (r<self.r2):
            return self.INT2(r)
        else :
            return self.INT0(r)


def B ( self, r):
    return (3.*(self.gamma_ad-1)*(self.gamma_ad+1.)*self.E_SN*E51/(pow(self.M(r),2. )*pow(r*parsec,self.beta)))

def Ushock2_r(self,r) : # entrée en parsec , sortie en cm/s
    return pow(self.B(r)*self.INT(r), 1./2.)

def Ushock1_r(self ,r):
    tchange=260.*pow(self.Mej/1.4,5./6.)*pow(self.E_SN,-1./2.)*pow(self.n0,-1./3.)*pow(10,-3.)
    rchange= 5.3*pow((self.E_SN/(self.n0*self.Mej)),1./7. )*pow(tchange,4./7.)
    
    if ( r<rchange):
        return pow(10,5.)* 2.7*pow(10,3.)*pow((5.3/r),3./4.)*pow((self.E_SN/(self.n0*self.Mej)),1./4.)
    
    else:
        return pow(10,5.)*1.7*pow(10,3.)*pow((4.3/r),3./2.)*pow((self.E_SN/self.n0),1./2.)


def t_r (self, r ):  # entrée en parsec, sortie en kyear
    N=100
    RB=np.linspace(0.001,r,N)  # in parsec
    tb=0.
    for i in range (1,N):
        tb=tb+(1./self.Ushock2_r(RB[i])+1./self.Ushock2_r(RB[i-1]))*(RB[i]-RB[i-1])/2.*parsec

    return tb/(3.*pow(10.,10.))


def Rshock_type2 (self,  t ):
    step_integration1=200 # 200 ok value
    RB=np.logspace(-2.,2.,step_integration1)  #maximum distance 100 pc
    temp=0.
    p=1
    while(temp < t*3.*pow(10, 10.)and p<step_integration1-1 ):
        temp=temp+(1./self.Ushock2_r( RB[p]))*(RB[p]-RB[p-1])*parsec
        p=p+1
    return (RB[p-1]+RB[p])/2.

#--------------------------------------------------------
# PRELIMINARY FUNCTIONS :
def M1( self, r):
    return self.Mej*Msol+4*3.1415*self.alpha1*r*parsec
    
def M2( self, r):
    return self.M1(self.r1)+4./3.*3.1415*self.rho2*(pow(r*parsec,3.)-pow(self.r1*parsec,3.))
    
def M0( self, r):
    return self.M2(self.r2)+4./3.*3.1415*self.rho0*(pow(r*parsec,3.)-pow(self.r2*parsec,3.))
    
def M (self, r):
    if (r<self.r1):
        return self.M1(r)
    else :
        if (r<self.r2):
            return self.M2(r)
        else :
            return self.M0(r)

def Tequalitymass2(self):
    rr=0.0001
    rgrid=0.01 # 0.00005
    while ( self.M (rr)/Msol-self.Mej <self.Mej ):
        rr=rr+rgrid
    
    if (rr>= self.r1):
        rr=self.r1 # to make sure the knee is not in the bubble, otherwise Emax >10**5 at early stage
    return self.t_r(rr)

#--------------------------------------------------------
# IN SNR PROCEDURE, RSHOCK AND USHOCK DEFINITIONS
def Rshock_1(self,t):
    tchange=260.*pow(self.Mej/1.4,5./6.)*pow(self.E_SN,-1./2.)* pow(self.n0,-1./3.)*pow(10,-3.)
    if (t<tchange):
        return 5.3*pow(  pow(self.E_SN,2.)/(self.n0*self.Mej), 1./7.)*pow(t,4./7.)
        
    else :
        return 4.3*pow(self.E_SN/self.n0,1./5.)*pow(t,2./5.)*pow(1.- 0.06*pow(self.Mej,5./6.)/(pow(self.E_SN,1./2.)*pow(self.n0,1./3.)*t) ,2./5.)

def Ushock_1(self,t): # returned value in cm/s
    
    tchange=260.*pow(self.Mej/1.4,5./6.)*pow(self.E_SN,-1./2.)*pow(self.n0,-1./3.)*pow(10,-3.);
    
    if (t<tchange):
        return pow(10,5.)*2.7*pow(10,3.)*pow(  pow(self.E_SN,2.)/(self.n0*self.Mej), 1./7.)*pow(t,-3./7.)
        
    else :
        return pow(10,5.)*1.7*pow(10,3.)*pow(self.E_SN/self.n0,1./5.)*pow(t,-3./5.)*pow(1- 0.06*pow(self.Mej,5./6.)/(pow(self.E_SN,1./2.)*pow(self.n0,1./3.)*t) ,-3./5.)


def Rshock_2(self,t):
    r_trans=self.Mej*Msol/(4.*np.pi*self.alpha1*parsec)
    #    t_trans=0.1;//pow(r_trans/7.7,8./7.)*pow((pow(E_SN,7./2.)*uw6/(Mdot*pow(Mej,5./2.))),-1./7.);
    t_trans=self.Tequalitymass2()
    Norm_FE=self.Rshock_type2 (t_trans )/pow(t_trans,7./8.)
    if (t<t_trans):
        return Norm_FE*pow(t,7./8.)
    else :
        return self.Rshock_type2 ( t )

def Ushock_2(self,  t):
    r_trans=self.Mej*Msol/(4.*np.pi*self.alpha1*parsec)
    t_trans=self.Tequalitymass2()
    Norm_FE=self.Ushock2_r (self.Rshock_2(t_trans))*pow(t_trans,1./8.)
    if (t<t_trans and t>0.):
        return Norm_FE*pow(t,-1./8.)
    else :
        return self.Ushock2_r (  self.Rshock_2(t))



def Rshock_t(self,t):
    if (self.type==1):
        return self.Rshock_1(t)
    else :
        return self.Rshock_2(t)

def Ushock_t(self,t):
    if (self.type==1):
        return self.Ushock_1(t)
    else :
        return self.Ushock_2(t)



def Transitiontime1(self):
    return 260.*pow(self.Mej/1.4,5./6.)*pow(self.E_SN,-1./2.)*pow(self.n0,-1./3.)*pow(10,-3.)


def Transitiontime2(self):
    return self.Tequalitymass2()


def assign_type (self):
    a=np.random.uniform(0.,1.)
    if ( a<0.32 ) :
        self.type=1
    else :
        if (a<0.76) :
            self.type=2
        else :
            if(a<0.98):
                self.type=3
            else :
                self.type=4


def assign_age(self,age_sample):
    self.age= np.random.uniform(0.,1.)*age_sample
    self.TIME=np.linspace(self.age,self.age,1)


#--------------------------------------------------------
# DISTRIBUTION TYPE Ia  - Bad but OK, should be improved soon
DISTRIB_TYPE_IA=[2.21*pow(10, -2.),2.32*pow(10, -2.),2.38*pow(10, -2),2.27*pow(10, -2),2.12*pow(10, -2),1.89*pow(10, -2),1.68*pow(10, -2),1.49*pow(10, -2),1.33*pow(10, -2), 1.21*pow(10, -2),1.07*pow(10, -2),9.54*pow(10, -3), 8.48*pow(10, -3), 7.36*pow(10, -3), 6.55*pow(10, -3), 5.68*pow(10, -3), 4.09*pow(10, -3), 2.44*pow(10, -3), 1.17*pow(10, -3), 6.99*pow(10, -4)]

DISTANCE_DISTRIB_TYPE_IA=[2.02,2.73, 3.42,4.13,4.77, 5.39, 5.95, 6.42, 6.86, 7.30, 7.65, 8.04, 8.39, 8.75, 9.12, 9.52, 10.4, 11.7, 13.4,14.6]


def distribution_type_Ia ( r): #r in kiloparsec  5.55= to normalize to 1.
    if (r<DISTANCE_DISTRIB_TYPE_IA[0]):
        return 5.52*0.042/DISTANCE_DISTRIB_TYPE_IA[0] #  0.042= to have the ratio 0.12/0.4 btw bulge et whole galaxy
    else :
        if (r>DISTANCE_DISTRIB_TYPE_IA[19]):
            return 0.
        else :
            j=1
            while (r>DISTANCE_DISTRIB_TYPE_IA[j]):
                j=j+1 #when it stops r btw DIST[j-1]and DIST[j]
    
            return 5.52*DISTRIB_TYPE_IA[j-1]+(r-DISTANCE_DISTRIB_TYPE_IA[j-1])/(DISTANCE_DISTRIB_TYPE_IA[j]-DISTANCE_DISTRIB_TYPE_IA[j-1])*(DISTRIB_TYPE_IA[j]-DISTRIB_TYPE_IA[j-1])



def  integrate_distribution_type_Ia ( r):
    j=0
    step_int=100
    R_INTEGRATION=np.linspace(0,r,step_int)
    integrale=0.
    for j in range (1,len(R_INTEGRATION)):
        integrale=integrale+distribution_type_Ia((R_INTEGRATION[j]+R_INTEGRATION[j-1])/2.)*(R_INTEGRATION[j]-R_INTEGRATION[j-1])
    return integrale

step_integration_radial_distribution_typeIa=100
RADIAL_SNR_DISTRIBUTION_TYPEIA=np.linspace(0,17,step_integration_radial_distribution_typeIa)


def place_SN_typeIa (self):
    q=0
    a=np.random.uniform(0.,1.)

    while (integrate_distribution_type_Ia(RADIAL_SNR_DISTRIBUTION_TYPEIA[q])< a and q<step_integration_radial_distribution_typeIa-1):
        q=q+1
        self.pos_r=(RADIAL_SNR_DISTRIBUTION_TYPEIA[q]+RADIAL_SNR_DISTRIBUTION_TYPEIA[q-1])/2.


def place_SN_theta_typeIa (self):
    self.pos_theta=np.random.uniform(0.,2.*np.pi)

def place_SN_z_typeIa (self):
    size=1000
    Z=np.linspace(-1.5,1.5,size)
    I=np.zeros(size)
    for i in range (1,size):
        I[i]=I[i-1]+nH1(self.pos_r, (Z[i-1]+Z[i])/2.)*(Z[i]-Z[i-1])
    p=0
    a=np.random.uniform(0.,I[size-1])
    while (I[p]<a and p<len(I)-1):
        p=p+1
    self.pos_z=(Z[p]+Z[p+1])/2.



#--------------------------------------------------------
# DISTRIBUTION TYPE II  - Bad but OK, should be improved
rho0_SNR=1.96 # normalisation de rho_SNR;
r0=17.2 # //kiloparsec
beta0=0.13  # kiloparsec
theta0=0.08
rhomax_SNR=1.8506

def Rrho_SNR(r ): # // r entrée en kiloparsec
    if (r<=r0*(1-theta0/np.pi)):
        return r*rho0_SNR*np.sin(np.pi*(r/r0)+theta0)*np.exp(-r*beta0)
    else :
            return 0.

number_arm=4
K_galaxy = [4.25,4.25,4.89,4.89]
r0_galaxy=[3.48,3.48,4.90,4.90]
theta0_galaxy=[0.,np.pi,2.52,-0.62]

step_integration_radial_distribution=100
RADIAL_SNR_DISTRIBUTION=np.linspace(0.,17., step_integration_radial_distribution)
INTEGRALE_SNR_RADIAL_DISTRIBUTION=np.zeros(step_integration_radial_distribution)
for i in range (1,step_integration_radial_distribution):
    INTEGRALE_SNR_RADIAL_DISTRIBUTION[i]=INTEGRALE_SNR_RADIAL_DISTRIBUTION[i-1]+Rrho_SNR((RADIAL_SNR_DISTRIBUTION[i]+RADIAL_SNR_DISTRIBUTION[i-1])/2.)*(RADIAL_SNR_DISTRIBUTION[i]-RADIAL_SNR_DISTRIBUTION[i-1])



def place_SN_typeII (self):
    a=np.random.uniform(0., INTEGRALE_SNR_RADIAL_DISTRIBUTION[step_integration_radial_distribution-1])
    q=0
    while (INTEGRALE_SNR_RADIAL_DISTRIBUTION[q]<a):
        q=q+1
            
    self.pos_r=(RADIAL_SNR_DISTRIBUTION[q]+RADIAL_SNR_DISTRIBUTION[q-1])/2.



def place_SN_theta_typeII (self):
    p=np.random.randint(0, 3)
    theta_correction=np.random.uniform(0, 2*np.pi)
        
    theta=(K_galaxy[p]*np.log(self.pos_r/r0_galaxy[p])+theta0_galaxy[p]+theta_correction*np.exp(-0.35*self.pos_r))
    #       Correction postion with respect to R and theta
    self.pos_r=np.random.normal(self.pos_r, 0.07*self.pos_r)

    while (theta>2*np.pi):
        theta=theta-2*np.pi

    while (theta<0):
        theta=theta+2*np.pi
    self.pos_theta=theta


def place_SN_z_typeII (self):
    size=1000
    Z=np.linspace(-1.5,1.5,size)
    I=np.zeros(size)
    for i in range (1,size):
        I[i]=I[i-1]+nH2(self.pos_r, (Z[i-1]+Z[i])/2.)*(Z[i]-Z[i-1])
    p=0
    a=np.random.uniform(0.,I[size-1])
    while (I[p]<a and p<len(I)-2):
        p=p+1
    self.pos_z=(Z[p]+Z[p+1])/2.



def place_supernova(self):
    if(self.type==1):
        self.place_SN_typeIa()
        self.place_SN_theta_typeIa ()
        self.place_SN_z_typeIa ()
    
    else :
        self.place_SN_typeII()
        self.place_SN_theta_typeII ()
        self.place_SN_z_typeII()


def typical_associated_parameters(self):
    if (self.type==1) :
        self.Mej=1.4
        self.Mdot=1.
        self.E_SN=1.
    else:
        if (self.type==2):
            self.Mej=8.
            self.Mdot=1.
            self.E_SN=1.
        else:
            if(self.type==3):
                self.Mej=2.
                self.Mdot=1.

                self.E_SN=1.
            else:
                if(self.type==4):
                    self.Mej=1.
                    self.Mdot=10.
                    self.E_SN=3.



def GG (lambda1, x  ) :  # fonction réelle, la variable doit être entre 0 et 1 ;
    v= -lambda1
    s=0
    p=0 # variable compteur
    while (s<x):
        p=p+1
        v=v+np.log(lambda1)-np.log(p)
        s=s+np.exp(v)
    #    la boucle stop : donc variable est entre p-1 et p
    return p-1



def draw_number_SNRs_in_Galaxy(age_sample):
    SN_rate=3.
    mean=SN_rate*age_sample*10. #mean number of SN we want
    return GG(mean,np.random.uniform(0.,1.))

#end of Sedov Taylor phase:
def Time_cooling (self,t):# t in kyear , returns kyr

    if (self.type==1):
        temp=10**3.*(self.rho0/masseproton)**(-1.)*(self.Ushock_1(t)*10**-8.)**3.
    else:
        temp=10**3.*(self.rho_r(self.Rshock_2(t))/masseproton)**(-1.)*(self.Ushock_2(t)*10**-8.)**3.
    return temp


def find_t_end(self):
    t=10**(-3.)
    print('type =', self.type, '      Time_cooling(t)= ', self.Time_cooling(t), 'Ushock = ', self.Ushock_2(t), ' Rshock = ', self.Rshock_2(t)  )
    while(t<self.Time_cooling(t) and t< age_sample):
        t=t*10**0.05
    return t
    
    
    
    


def  density_around(self):
    self.tau=0.
    r=self.pos_r
    z=self.pos_z
    if (self.type==1):
        self.n0=nH(r,z)
        self.rho0=self.n0*(1.67*pow(10,-24.))
    
    else :
        self.alpha1=alpha1_def(self.Mdot,self.uw6)
        
        self.n0=nH(r, z)
        self.rho0=self.n0*(1.67*pow(10,-24.))
        
        self.n2=0.01*pow(pow(self.n0,19.)*pow(self.t6,-22.),1/35.)
        self.rho2= self.n2*(1.67*pow(10,-24.))
        
        self.Tb=1.6*pow(10., 6.)*pow(self.n0, 2/35.)*pow(self.V_2000*self.Mdot/10, 8/35.)*pow(self.t6,-6/35.)
        
       # self.r1=(self.alpha1*pow(self.uw6*pow(10, 6),2.))/(4*np.pi*self.n2*self.Tb*k_boltzmann_erg*pow(parsec, 2.))
        self.r1=(self.Mdot*10**(-5)*Msol/(year)*self.uw6*10**6./(4*np.pi*self.n2*k_boltzmann_erg*self.Tb))**0.5/parsec

        self.r2 = 28.*pow(self.L36/(mu*self.n0) ,1./5.)*pow(self.t6,3./5.) # 35 ?
        self.r2 = min(30.,self.r2) # 35 ?

    self.t_end= self.find_t_end()
        
        
        
        

def  set_density(self):

    if (self.type==1):
        self.rho0=self.n0*(1.67*pow(10,-24.))

    else :
        self.alpha1=alpha1_def(self.Mdot,self.uw6)

        self.rho0=self.n0*(1.67*pow(10,-24.))

        self.n2=0.01*pow(pow(self.n0,19.)*pow(self.t6,-22.),1/35.)
        self.rho2= self.n2*(1.67*pow(10,-24.))

        Mwolfrayet=1.
        V_2000=1.
        self.Tb=1.6*pow(10., 6.)*pow(self.n0, 2/35.)*pow(self.V_2000**2.*Mwolfrayet, 8/35.)*pow(self.t6,-6/35.)

        self.r1=(self.Mdot*10**(-5)*Msol/(year)*self.uw6*10**6./(4*np.pi*self.n2*k_boltzmann_erg*self.Tb))**0.5/parsec
        self.r2 = 28.*pow(self.L36/(mu*self.n0) ,1./5.)*pow(self.t6,3./5.) # 35 ?




def calculate_final_evolution (self):
    if(self.type==1):
        self.Rshock=self.Rshock_1(self.age)
        self.Ushock=self.Ushock_1(self.age)
    else :
        self.Rshock=self.Rshock_2(self.age)
        self.Ushock=self.Ushock_2(self.age)


def set_factor_Emax(self,KNEE): # modified 03.12.19
    if (self.type==1) :
        self.Emax_factor= 1. # KNEE/(self.Emax_t(self.Transitiontime1()))
    else :
        self.Emax_factor= 1. # KNEE/(self.Emax_t(self.Transitiontime2()))



#------ DISTANCE FUNCTIONS from NOVA
def distance_plan (self):
    X=np.sin(self.pos_theta)/(distance_SS_GC/self.pos_r-np.cos(self.pos_theta))
    d1=np.absolute(1/(np.sqrt(1+X**2.))*(-self.pos_r*np.cos(self.pos_theta) + self.pos_r*np.sin(self.pos_theta)*X + distance_SS_GC  ) )
    return d1

def b_calculated(self):
    return np.arctan(self.pos_z/self.distance_plan())*(180/np.pi)
    
def l_calculated(self):
    return np.sign(self.pos_theta)*np.arccos((-self.pos_r*np.cos(self.pos_theta)+distance_SS_GC)/self.distance_plan())*(180/np.pi)
    
    
def distance (self):
    X=np.sin(self.pos_theta)/(distance_SS_GC/self.pos_r-np.cos(self.pos_theta))
    d1=np.absolute(1/(np.sqrt(1+X**2.))*(-self.pos_r*np.cos(self.pos_theta) + self.pos_r*np.sin(self.pos_theta)*X + distance_SS_GC  ) )
    self.dist=np.sqrt(d1**2.+self.pos_z**2.)
    return self.dist
    


# ---------------------------------------------------------------------------#
#  EMAX FUNCTIONS to match the knee and corresponding amplified field

def Emax_t (self, t):
    self.Emax_mem=self.Emax_factor*self.Emax_simple(t)
    return  self.Emax_mem

def B_amplified_OK2015 (self, t):
    xi_B_correction_factor=1.
    self.M_A=23.
    if(self.type==1):
        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho0)**(0.5)
        return self.B0*self.sigma*np.sqrt(pow(self.Ushock_t(t),2.)/pow((self.M_A)*VA0*pow(10, -6.), 2.)+1.)
    else :
        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho_r(self.Rshock_t(t)))**(0.5)
        return self.B0*self.sigma*np.sqrt(pow(self.Ushock_t(t),2.)/pow((self.M_A)*VA0*pow(10, -6.), 2.)+1.)


def B_amplified_knee ( self, t):
    calc=self.Emax_factor*self.B_amplified_OK2015(t)
    return max(self.sigma*self.B0, calc)
    


#def Emax_OK2015 (self, t): # // in TeV
#    if (self.type==1):
#        self.delta_Emax=-np.log10(self.EMAX_B_dependant_t_OK2015(t+2.)/self.EMAX_B_dependant_t_OK2015(t))/np.log10(self.Rshock_t(t+2.)/self.Rshock_t(t))
#        A_test=self.EMAX_B_dependant_t_OK2015(t)/pow(self.Rshock_t(t),-self.delta_Emax)
#        return A_test*pow(self.Rshock_t(t),-self.delta_Emax)
#
#    else :
#        self.delta_Emax=2.;
#        return self.EMAX_B_dependant_t_OK2015(t)
#

#
#def EMAX_B_dependant_t_OK2015 ( self, t):
#    xi_B_correction_factor=1.
#
#    if(self.type==1):
#        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho0)**(0.5)
#        U1=self.Ushock_t(t)
#        R1=self.Rshock_t(t)
#
#        return pow(10, -12.)*3*pow(10, -6.)*pow(10, -8.)*self.chi*R1*parsec*U1*(self.B_amplified_OK2015(t)/self.sigma)*1./(pow(pow((VA0*self.M_A*pow(10, -6.))/U1, 2.)+1, 1.5))
#
#    else :
#        R2=self.Rshock_t(t)
#        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho_r(R2))**(0.5)
##     print ('JOJOOJJJO' )
##       print('self.rho_r(R2)=', self.rho_r(R2), 'self.rho0', self.rho0, '   R2 =', R2, '  t = ', t, '   r2 = ', self.r2)
#        U2=self.Ushock_t(t)
#        return pow(10, -12.)*3*pow(10, -6.)*pow(10, -8.)*self.chi*R2*parsec*U2*(self.B_amplified_OK2015(t)/self.sigma)*1./(pow(pow((VA0*self.M_A*pow(10, -6.))/U2, 2.)+1, 1.5))



def Emax_simple ( self, t):
    xi_B_correction_factor=1.
    
    if(self.type==1):
        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho0)**(0.5)
        U1=self.Ushock_t(t)
        R1=self.Rshock_t(t)

        return pow(10, -12.)*3*pow(10, -6.)*pow(10, -8.)*self.chi*R1*parsec*U1*(self.B_amplified_OK2015(t)/self.sigma)*1./(pow(pow((VA0*self.M_A*pow(10, -6.))/U1, 2.)+1, 1.5))

    else :
        R2=self.Rshock_t(t)
        VA0=xi_B_correction_factor*self.B0/(4*np.pi*self.rho_r(R2))**(0.5)
        U2=self.Ushock_t(t)
        return pow(10, -12.)*3*pow(10, -6.)*pow(10, -8.)*self.chi*R2*parsec*U2*(self.B_amplified_OK2015(t)/self.sigma)*1./(pow(pow((VA0*self.M_A*pow(10, -6.))/U2, 2.)+1, 1.5))





# ---------------------------------------------------------------------------#
#  EMAX electrons
    
def Emax_electron_vannoni_time (self, t): #  // in TeV  , B1 in microgauss
    rtot=6.  # 5.2
    B1=self.B0  #; // B0  // 3 ?
    beta_B=1./(pow(2*pow(rtot,2)/(3)+1/3,0.5)) # // compression ratio of the B field
    K_adjust=1.
    B2=self.B_amplified_knee(t) # //     double B2=etha_elec*B_volk(tage);

    Ush=self.Ushock_t(t)
    if (B1<5.): # { // B1< 10 formely
        E_temp= mc2_electron_TeV*pow(((1-1./rtot)*pow(Ush, 2.)*900./(K_adjust)  )/( ( ((B1*pow(10, -6.))/(8*np.pi)+self.Urad/(B1*pow(10, -6.)) +rtot*B2*pow(10, -6.))/(8*np.pi)+self.Urad/(B2*pow(10, -6.)))*4.*sigma_thomson*pow(c, 2.)*6.241*pow(10, 11.)),1./2.)
    else :
        E_temp= mc2_electron_TeV*pow(((1-1./rtot)*pow(Ush, 2.)*8*np.pi*900./(K_adjust)  )/( ( B1+rtot*B2)*pow(10, -6.)*4.*sigma_thomson*pow(c,2.)*6.241*pow(10,11.)),1./2.)

#  E_temp2=0.05*pow(10, -12.)*(900.)*Ush1*B1*pow(10, -6.)*self.Rshock_t(t)*parsec/c  # // si tacc=tage
    self.Emax_electron_vannoni_time_mem=min(self.Emax_t(t), E_temp)
    return self.Emax_electron_vannoni_time_mem



def Estar_electron_time ( self, t): # // E in TeV, B in MicroGauss
    B1=self.B0
    expression = 0.624*((pow(mc2_electron_erg,2.))/((4./3.)*sigma_thomson*c*(pow(self.B0*pow(10, -6.),2.)/(8*np.pi))))*(1./(t*kyear_sec)-self.Ushock_t(t)/(self.Rshock_t(t)*parsec))
    self.Estar_electron_time_mem=max(0.001, expression)
    return self.Estar_electron_time_mem
    

    

# ---------------------------------------------------------------------------#
#  NORMALIZATION FOR THE GAMMA RAYS FROM SNR

#def A( self, r,t): # return in units of p**-3 cm**-3 with p in TeV/c
#    a=2.-self.alpha
#    E0=1.
#    if (self.type==1):
#        RR=self.Rshock_t(t)
#        return ( 3.*self.eta*a*self.rho0*0.624*pow(self.Ushock1_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (1-self.sigma)*(-4+a) ) )/((pow(r/RR,-a*self.delta_Emax*self.sigma)*pow(self.Emax_t(t),a)-pow(masseproton_TeV,a)))
#
#    else :
#
#        RR=self.Rshock_t(t)
#        return ( 3.*self.eta*a*self.rho_r(pow(r/RR,self.sigma-1.)*r)*0.624*pow(self.Ushock2_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (1-self.sigma)*(-4+a)))/((pow(r/RR,-a*self.delta_Emax*self.sigma)*pow(self.Emax_t(t),a)-pow(masseproton_TeV,a )))


def v_E(E): # input E in TeV, output in cgs
    return c

def v_p (p): #p in mc units, out in cm.s**-1
    temp= (p*masseproton*c)*c**2./((p*masseproton*c)**2*c**2.+ masseproton**2.*c**4.)**0.5
    return temp


#def A( self, r,t): # return in units of p**-3 cm**-3 with p in TeV/c
#    a=self.alpha+2.
#    E0=1. #TeV
#    if (self.type==1):
#        inte=0.
#        Emin=10**-3.# TeV
#        E_GRID_LOCAL=np.logspace(np.log10(Emin),np.log10(self.Emax_t(t)),40.)
#        for i in range (1,len(E_GRID_LOCAL)):
#            inte=inte+(v_E(E_GRID_LOCAL[i])*(E_GRID_LOCAL[i]/E0)**(3-a) + v_E(E_GRID_LOCAL[i-1])*(E_GRID_LOCAL[i-1]/E0)**(3-a))*(E_GRID_LOCAL[i]-E_GRID_LOCAL[i-1])/2.
#        inte=inte*TeV_to_erg/c
#        RR=self.Rshock_t(t)
#        calc=3.*self.eta*a*self.rho0*pow(self.Ushock1_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (self.sigma-1)*(a/3.) )
#        calc=calc/(E0*TeV_to_erg*inte)
#        calc=calc/(self.sigma) #accounting for velocity downstream of the shock + heavier nuclei
#        return calc
#
#
#    else :
#
#        inte=0.
#        Emin=10**-3.# TeV
#        E_GRID_LOCAL=np.logspace(np.log10(Emin),np.log10(self.Emax_t(t)),20.)
#        for i in range (1,len(E_GRID_LOCAL)):
#            inte=inte+(v_E(E_GRID_LOCAL[i])*(E_GRID_LOCAL[i]/E0)**(3-a) + v_E(E_GRID_LOCAL[i-1])*(E_GRID_LOCAL[i-1]/E0)**(3-a))*(E_GRID_LOCAL[i]-E_GRID_LOCAL[i-1])/2.
#        inte=inte*TeV_to_erg/c
#
#        RR=self.Rshock_t(t)
#
#        calc=3.*self.eta*a*self.rho_r(pow(r/RR,self.sigma-1.)*r)*pow(self.Ushock2_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (self.sigma-1)*(a/3.) )
#        calc=calc/(E0*TeV_to_erg*inte)
#        calc=calc/(self.sigma) #accounting for velocity downstream of the shock + heavier nuclei
#        return calc


def A( self, r,t): # return in units of p**-3 cm**-3 with p in TeV/c
    a=self.alpha+2.
    E0=1. #TeV
    if (self.type==1):
        inte=0.
        pmin=10**-3.# mc units
        pmax=self.Emax_t(t)*10**3. # GeV/c
      #  print(' pmax = ', pmax)
        P_GRID_LOCAL=np.logspace(np.log10(pmin),np.log10(pmax),30)
     #   print('P_GRID_LOCAL = ', P_GRID_LOCAL)
     #   print('v_p(P_GRID_LOCAL)= ', v_p(P_GRID_LOCAL))
        P0=E0*10**3. # GeV/c
        for i in range (1,len(P_GRID_LOCAL)):
            inte=inte+(v_p(P_GRID_LOCAL[i])*(P_GRID_LOCAL[i]/P0)**(3-a) + v_p(P_GRID_LOCAL[i-1])*(P_GRID_LOCAL[i-1]/P0)**(3-a))*(P_GRID_LOCAL[i]-P_GRID_LOCAL[i-1])/2.
        inte=inte*10**-3.*TeV_to_erg/c
        RR=self.Rshock_t(t)
        calc=3.*self.eta*self.rho0*pow(self.Ushock1_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (self.sigma-1)*(a/3.) )
        calc=calc/(E0*TeV_to_erg*inte)
        calc=calc/(self.sigma) #accounting for velocity downstream of the shock + heavier nuclei
        return calc
    

    else :

        inte=0.
        pmin=10**-3.# mc units
        pmax=self.Emax_t(t)*10**3. # GeV/c
     #  print(' pmax = ', pmax)
        P_GRID_LOCAL=np.logspace(np.log10(pmin),np.log10(pmax),40.)
    #   print('P_GRID_LOCAL = ', P_GRID_LOCAL)
    #   print('v_p(P_GRID_LOCAL)= ', v_p(P_GRID_LOCAL))
        P0=E0*10**3. # GeV/c
        for i in range (1,len(P_GRID_LOCAL)):
            inte=inte+(v_p(P_GRID_LOCAL[i])*(P_GRID_LOCAL[i]/P0)**(3-a) + v_p(P_GRID_LOCAL[i-1])*(P_GRID_LOCAL[i-1]/P0)**(3-a))*(P_GRID_LOCAL[i]-P_GRID_LOCAL[i-1])/2.
        inte=inte*10**-3.*TeV_to_erg/c

        RR=self.Rshock_t(t)
        
        calc=3.*self.eta*self.rho_r(pow(r/RR,self.sigma-1.)*r)*pow(self.Ushock2_r(RR*pow(r/RR, self.sigma)),2.)*pow(r/RR , (self.sigma-1)*(a/3.) )
        calc=calc/(E0*TeV_to_erg*inte)
        calc=calc/(self.sigma) #accounting for velocity downstream of the shock + heavier nuclei
        return calc









def density_inside ( self, r,t): #// r en parsec, t en kyear
    if (self.type==1):
        return self.rho0*self.sigma*pow((r/self.Rshock_t(t)),3*(self.sigma-1))
    else :
        RR= self.Rshock_t(t)
        return self.sigma*self.rho_r(pow(r,self.sigma)/pow(RR,self.sigma-1.))*pow(r/RR,3*(self.sigma-1.))


def Norm_hadronic (self, t):
    N=50
    
    if (t>=self.t_end or self.Rshock_t(t)>=self.r2):
        self.Norm_hadronic_mem=0.
    else :
        Rsh=self.Rshock_t(t)
        R=np.linspace(0.,Rsh,N)
        norm=0.
        for i in range (0,len(R)-1):
            norm=norm+ pow(((R[i]+R[i+1])/2.)*parsec,2.)*self.density_inside(((R[i]+R[i+1])/2.), t)*self.A(((R[i]+R[i+1])/2.),t)*(R[i+1]-R[i])*parsec/masseproton
        self.Norm_hadronic_mem =4*np.pi*norm#/(4*np.pi*TeV_to_erg)

    #4*np.pi/c for conversion from f(p) to f(E)
    return self.Norm_hadronic_mem

# Norm leptonic, does not include denisty of target photons inside SNR:
def Norm_leptonic (self, t):
    N=50
    
    if (t>=self.t_end or self.Rshock_t(t)>=self.r2):
        self.Norm_leptonic_mem=0.
    else :
        Rsh=self.Rshock_t(t)
        R=np.linspace(0.,Rsh,N)
        norm=0.
        for i in range (0,len(R)-1):
            norm=norm+ pow(((R[i]+R[i+1])/2.)*parsec,2.)*self.A(((R[i]+R[i+1])/2.),t)*(R[i+1]-R[i])*parsec
        self.Norm_leptonic_mem =4*np.pi*norm*self.Kep#/(4*np.pi*TeV_to_erg)

    #4*np.pi/c for conversion from f(p) to f(E)
    return self.Norm_leptonic_mem


#--------------------------------------------------------
# Startting here we define functions outside the class
#calculate GAMMA - with the help of NAIMA
def spectrum_proton(self,time):
    return naima.models.ExponentialCutoffPowerLaw(self.Norm_hadronic(time)/u.erg,1*u.TeV,self.alpha,self.Emax_t(time) *u.TeV,beta=1)
    
def spectrum_proton_old_school (self, E,time): # this is needed for the secondaries
    return self.Norm_hadronic(time)/(0.624)*(E/1.)**(-self.alpha)*np.exp(-(E/self.Emax_t(time)))
    
    
def spectrum_electron(self,time):
    return naima.models.ExponentialCutoffBrokenPowerLaw(self.Norm_leptonic(time)/u.erg,1*u.TeV,self.Estar_electron_time(time) *u.TeV,self.alpha,self.alpha+1.,self.Emax_electron_vannoni_time(time) *u.TeV,beta=1)

    
def diff_spectrum_hadronic(self,time):
    self.dist=self.distance()
    PROTONS=self.spectrum_proton(time)
    PIONS_FROM_SHELL=PionDecay(PROTONS, nh=1.* u.cm** -3,Epmax=100* u.PeV)
    GAMMAS=PIONS_FROM_SHELL.sed(self.ENERGY,self.dist * u.kpc)
        #  GAMMAS.to(u.eV/(u.cm **2 * u.s))
    GAMMAS=GAMMAS.to(u.TeV/(u.cm**2 *u.s))
    return GAMMAS
    
    
def diff_spectrum_leptonic(self,time):
    self.dist=self.distance()
    ELECTRONS=self.spectrum_electron(time)
    IC = InverseCompton(ELECTRONS, seed_photon_fields=['CMB'])
    GAMMAS=IC.sed(self.ENERGY,self.dist * u.kpc)
    GAMMAS_TeV=GAMMAS.to(u.TeV/(u.cm**2 *u.s))
    return GAMMAS_TeV
    
def diff_spectrum_total (self,time):
    GAMMA_TOT=self.diff_spectrum_hadronic(time)+self.diff_spectrum_leptonic(time)
    return GAMMA_TOT
    
    
def calculate_alpha_gamma (self,E):
    p=1
    E_LOCAL=np.array(self.ENERGY)
    while (E>E_LOCAL[p] and p<len(E_LOCAL-1)):
        p=p+1
    p_up=1
    while (10*E>E_LOCAL[p_up] and p_up<len(E_LOCAL-1)):
        p_up=p_up+1
        
    GAMMA=np.array(self.diff_spectrum_total()/self.ENERGY[p]**2)
    self.alpha_gamma_memory=-(np.log10(GAMMA[p_up])-np.log10(GAMMA[p]))/(np.log10(E_LOCAL[p_up]-E_LOCAL[p]))
    return self.alpha_gamma_memory


def calculate_diff_spectrum_TIME (self):
    self.LGAMMA_DIFF_T=np.zeros((len(self.TIME), len(self.ENERGY)))
    for t in range (0,len(self.TIME)):
        SPEC=np.array(self.diff_spectrum_total(self.TIME[t]))

        for i in range (0,len(self.ENERGY)):
            self.LGAMMA_DIFF_T[t][i]=SPEC[i]

def calculate_diff_spectrum_PWNE_TIME (self):
    self.LGAMMA_DIFF_T=np.zeros((len(self.TIME), len(self.ENERGY)))


#--------------------------------------------------------#
#       DEFINING CLASSES

class SNR:
    #    SN=Supernova()
    t_end=10.
    age=1.
    pos_r=1.
    pos_z=1.
    pos_theta=1.
    type=1
    Rshock=1.
    ushock=1.
    alpha=1. #indice spectre proton
    Kep=1.
    dist=1.
    size=1.
    Mej=1.
    Mdot=1.  # Mass loss rate wind
    E_SN=1.
    r1=1.
    r2=1.
    rho0=1.
    rho2=1.
    alpha1=1. #coefficient of the wind
    t6=1.
    L36=1.   #énergie cinétique du vent
    uw6=1. # wind velocity
    
    Tb=1. # Temperature in the bubble in K.
    V_2000=1. # // speed in 2000 km/s units
    tau=0. # absorption
    Urad=0.25*1.602*pow(10, -12) # erg.cm-3
    M_A=23.

    
    n0=1.
    B0=3.
    gamma_ad=4./3.
    beta=6*(gamma_ad-1)/(gamma_ad+1)
    age=1.
    
    sigma=4.
    
    chi=0.1  # fraction of diffusion length before escape
    eta=0.1 #efficiency of particle acceleration at the shock
    Emax=1.
    delta_Emax=1.
    Emax_factor=1.  # to match the knee
    Emax_mem=1.
    Estar_electron_time_mem=1.
    Emax_electron_vannoni_time_mem=1.
    Norm_hadronic_mem=1.
    Norm_leptonic_mem=1.

    distance=distance
    distance_plan=distance_plan
    b_calculated=b_calculated
    l_calculated=l_calculated
    
    
    # Metthods used to parametrize the SNRs
    
    Rshock_1=Rshock_1
    Rshock_2=Rshock_2
    Ushock_1=Ushock_1
    Ushock_2=Ushock_2
    Rshock_type2=Rshock_type2
    Ushock2_r=Ushock2_r
    Ushock1_r=Ushock1_r
    Rshock_t=Rshock_t
    Ushock_t=Ushock_t
    B=B
    
    Tequalitymass2=Tequalitymass2
    rho0=rho0
    rho_r=rho_r
    rho1_r=rho1_r
    
    find_t_end=find_t_end
    Time_cooling=Time_cooling

    M=M
    M0=M0
    M1=M1
    M2=M2
    t_r=t_r
    Transitiontime1=Transitiontime1
    Transitiontime2=Transitiontime2
    Emax_t=Emax_t
   # Emax_OK2015=Emax_OK2015
    Emax_simple=Emax_simple
    B_amplified_OK2015=B_amplified_OK2015
    B_amplified_knee=B_amplified_knee
  #  EMAX_B_dependant_t_OK2015=EMAX_B_dependant_t_OK2015
    INT=INT
    INT0=INT0
    INT1=INT1
    INT2=INT2

    
    set_factor_Emax=set_factor_Emax
    assign_type=assign_type
    assign_age=assign_age
    place_supernova=place_supernova
    typical_associated_parameters=typical_associated_parameters
    place_SN_typeIa=place_SN_typeIa
    place_SN_theta_typeIa=place_SN_theta_typeIa
    place_SN_z_typeIa=place_SN_z_typeIa
    place_SN_typeII=place_SN_typeII
    place_SN_theta_typeII=place_SN_theta_typeII
    place_SN_z_typeII=place_SN_z_typeII
    calculate_final_evolution=calculate_final_evolution
    density_around=density_around
    set_density=set_density
    
    # Arrays for spectra
    ENERGY=np.logspace(-1,4,40)* u.TeV
    #for a given time
    GAMMAS_H=np.zeros(len(ENERGY)) * u.TeV/(u.cm**2 *u.s)
    GAMMAS_L=np.zeros(len(ENERGY)) * u.TeV/(u.cm**2 *u.s)
    integrated_hadronic_memory=1/(u.cm**2 * u.s )
    integrated_leptonic_memory=1/(u.cm**2 * u.s )
    integrated_memory_total=1.
    
    
    # ---------------------------------------------------------------------------#
    #  NORMALIZATION FOR THE GAMMA RAYS FROM SNR
    
    time_step=0.1 # kiloyear
    number=int((age)/time_step)
    TIME=np.linspace(1.,1.,1)
    
    LGAMMA_HADRONIC_T=np.zeros((len(TIME), len(ENERGY)))
    LGAMMA_LEPTONIC_T=np.zeros((len(TIME), len(ENERGY)))
    LGAMMA_DIFF_T=np.zeros((len(TIME), len(ENERGY)))
    
    Emax_electron_vannoni_time=Emax_electron_vannoni_time
    Estar_electron_time=Estar_electron_time
    A=A
    density_inside=density_inside
    Norm_hadronic=Norm_hadronic
    Norm_leptonic=Norm_leptonic

    spectrum_proton=spectrum_proton
    spectrum_proton_old_school=spectrum_proton_old_school
    spectrum_electron=spectrum_electron

    diff_spectrum_hadronic= diff_spectrum_hadronic
    diff_spectrum_leptonic=diff_spectrum_leptonic
    diff_spectrum_total=diff_spectrum_total
    calculate_alpha_gamma=calculate_alpha_gamma
    calculate_diff_spectrum_TIME=calculate_diff_spectrum_TIME
    calculate_diff_spectrum_PWNE_TIME=calculate_diff_spectrum_PWNE_TIME

#  calculate_integrated_spectrum_SGSO=calculate_integrated_spectrum_SGSO





def one_realization_only_pevatrons (slope, Kep, D,eta, KNEE):
    N=draw_number_SNRs_in_Galaxy(age_sample)
    print('Number of simulated objects= ', N )

    LIST_SNR=[]
    for i in range (0,N):
        print( 'realisation i/N=', i*1./N*100,'%')

        SNR_temp=SNR()
        SNR_temp.assign_type()
        SNR_temp.place_supernova()
        SNR_temp.assign_age(age_sample)
        SNR_temp.typical_associated_parameters()
        SNR_temp.density_around()
        
        if (SNR_temp.age < age_max_ST_pase):
            SNR_temp.calculate_final_evolution()
            SNR_temp.set_factor_Emax(KNEE)
            SNR_temp.alpha=slope
            SNR_temp.Kep=Kep
            SNR_temp.distance()
            SNR_temp.size=2*SNR_temp.Rshock/(SNR_temp.dist*pow(10., 3.))*3437.75
            #### Calculating the gammas from the SNR :
            SNR_temp.eta=eta
            SNR_temp.calculate_diff_spectrum_TIME()
          #  print(' KEP= ', SNR_temp.Kep)
          #  print(' SNR_temp.Norm_leptonic_mem= ', SNR_temp.Norm_leptonic_mem, ' SNR_temp.Norm_hadronic_mem  =', SNR_temp.Norm_hadronic_mem, ' Emax = ', SNR_temp.Emax_mem, 'time_end = ', SNR_temp.t_end)
        else:
            SNR_temp.calculate_diff_spectrum_PWNE_TIME()
            
            
        LIST_SNR.append(SNR_temp)
    return LIST_SNR



def many_realizations (slope,Kep, D, eta, KNEE, M):
    BIG_LIST=[]
    for i in range (0, M):
        LIST_SNR=one_realization(slope, Kep, D,eta, KNEE)
        BIG_LIST.append(LIST_SNR)

    return BIG_LIST

def print_one_SNR (SNR):
    print ('----------------')
    print ('            ')
    print ('pos_r =', SNR.pos_r)
    print ('pos_theta=', SNR.pos_theta)
    print ('pos_z =', SNR.pos_z)
    print ('type =', SNR.type)
    print ('age =', SNR.age)
    print ('size =', SNR.size)
    print ('Rshock =', SNR.Rshock_t(SNR.age))
    print ('alpha =', SNR.alpha)
    print ('Norm =', SNR.Norm_hadronic(SNR.age))
    print ('Emax_proton =', SNR.Emax_t(SNR.age))
    print ('Kep =', SNR.Kep)
    print ('Estar =', SNR.Estar_electron_time (SNR.age))
    print ('Emax_electron =', SNR.Emax_electron_vannoni_time (SNR.age))
    print ('diff spectrum :')
    for i in range (0,len(SNR.ENERGY)):
        for t in range (0,len(SNR.TIME)):
            print ( 'time= ', SNR.TIME[t], '    bin=', SNR.ENERGY[i], '     diff=', SNR.LGAMMA_DIFF_T[t][i])
           
    print ('            ')


def save_one_LIST_to_file (LIST,file):
    with open(file, 'w') as text_file:
        writer = csv.writer(text_file, delimiter='\t')
        writer.writerow(["Num_SNR","pos_r", \
                         "pos_theta","pos_z","n0", \
                         "type", "age" ,"size", \
                         "Rsh","a","Norm",\
                         "Emax_proton","Kep","Estar", \
                         "Emax_electron","Time[kyr]","E[TeV]", \
                         "diff_spectrum[TeVcm-2s-1]"])

        for i in range (0,len(LIST)):
            ENERGY=np.array(LIST[i].ENERGY)
            for j in range (0,len(LIST[i].ENERGY)):
                #          print('i/Len(LIST) =', i/len(LIST)*(100.) , '%')
                ##print( ' Norm = ',  LIST[i].Norm_hadronic_mem, ' dist = ',LIST[i].dist )

                writer.writerow((i,LIST[i].pos_r, \
                                        LIST[i].pos_theta, LIST[i].pos_z,LIST[i].n0, \
                                     LIST[i].type, LIST[i].age, LIST[i].size, \
                                     LIST[i].Rshock,LIST[i].alpha, LIST[i].Norm_hadronic_mem, \
                                     LIST[i].Emax_mem,LIST[i].Kep,LIST[i].Estar_electron_time_mem, \
                                     LIST[i].Emax_electron_vannoni_time_mem ,LIST[i].age, ENERGY[j], \
                                     LIST[i].LGAMMA_DIFF_T[len(LIST[i].TIME)-1][j]))

def save_one_LIST_to_file_cyril (LIST,file):
    with open(file, 'w') as text_file:
        writer = csv.writer(text_file, delimiter='\t')
        writer.writerow(["Num_SNR","POS_X", \
                         "POS_Y","POS_Z","n0", \
                         "type", "age[kyr]" ,"size", \
                         "SNR_radius[pc]","E[TeV]", \
                         "diff_spectrum"])
            
            
        for i in range (0,len(LIST)):
            ENERGY=np.array(LIST[i].ENERGY)
            for j in range (0,len(LIST[i].ENERGY)):
                                 
                writer.writerow((i,LIST[i].pos_r*np.cos(LIST[i].pos_theta), \
                                      LIST[i].pos_r*np.sin(LIST[i].pos_theta), LIST[i].pos_z,LIST[i].n0, \
                                      LIST[i].type, LIST[i].age, LIST[i].size, \
                                      LIST[i].Rshock, ENERGY[j], \
                                      LIST[i].LGAMMA_DIFF_T[len(LIST[i].TIME)-1][j]))

def get_first_nbr_from_str(input_str):
    if not input_str and not isinstance(input_str, str):
        return 0
    out_number = ''
    for ele in input_str:
        if (ele == '.' and '.' not in out_number) or ele.isdigit():
            out_number += ele
        elif out_number:
            break
    return int(float(out_number))

