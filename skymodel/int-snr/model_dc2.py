from numpy import *
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import units as u
import scipy.stats as stats

def prob(p0,M,a,m0=1e4):
    P=p0*(M/m0)**a
    return P
def m(f):
    mi=[]
    ma=[]
    for i in range(len(f)):
        if f[i]==[]:
            pass
        else:
            mi.append(min(f[i]))
            ma.append(max(f[i]))
    MI=min(mi)
    MA=max(ma)
    return MI,MA,mi,ma

database_real='gSNR4.fits'
database_pierre='../snrs/ctadc_skymodel_gps_sources_snr_1.ecsv'
database_nubi='Spectra.fits'


cloud=Table.read(database_nubi)
real=Table.read(database_real)
pierre=Table.read(database_pierre,format='ascii.ecsv')

alpha=0
p0=0.015

def isnr(alpha=[-1,0,1,2],p0=[0.05,0.015,0.0002,2e-6],database_synth='Spectra.fits',database_real='gSNR4.fits',n_rea=1000,rep=1,plot_rea='True',plot_mean='True',plot_real='True',save_rea='False'):
    
    if len(alpha) != len(p0):
        quit()
    
    cloud=Table.read(database_synth)
    real=Table.read(database_real)

    w=where((real['Type'] == 'INT' )* (real['SNR_Name']!= 'W28south'))
    fluxes_REAL=array(real['Flux100MeV'][w])
    ff_REAL=sort(fluxes_REAL)[::-1]
    nn_REAL=arange(len(ff_REAL))+1

    a1=[]
    x1=[]
    f_array=[]
    inter_reali=[]
    inter_x=[]
    label=[]
    prova_fit=[]


    color=['tab:blue','tab:orange','tab:green','violet']
    
    
    for c in range(rep):
        for j in range (len(alpha)):
            
            print (alpha[j],p0[j])
            
            pp=prob(p0[j],cloud['Mass'],alpha[j])
            #print ('M0, P0, pp : ',m0[i],p0[j],sum(pp))
            realization=[]
            for k in range(n_rea):
                
                rand=random.random(len(cloud))
                is_int = rand < pp
                true_index=is_int.nonzero()
                w=true_index[0]  # Indeces of iSNR random select from synthetic catalog
                a1.append(len(w))
                realization.append(k)
                rp=[]
                for l in range(len(w)):
                    rp.append(random.randint(1,len(pierre)))
                icloud=Table()
                icloud['Mass']=cloud['Mass'][w]*u.solMass
                icloud['Glat']=pierre['glat'][rp]
                icloud['Glon']=pierre['glon'][rp]
                icloud['Flux_int_100MeV']=cloud['Flux_int_100MeV'][w]
                icloud['Luminosity']=cloud['Luminosity'][w]
                icloud['index1']=cloud['index1'][w]
                icloud['index2']=cloud['index2'][w]
                icloud['Energy_break']=cloud['Energy_break'][w]
                icloud['Density']=cloud['Density'][w]

                if rep==1 and save_rea=='True':
    
                    icloud.write('LogNlogS/a0_p0015/a0_p0015'+ str(k) + '_' + str(len(w)) + '.fits',format='fits',overwrite=True)
        
        
                f=sort(icloud['Flux_int_100MeV'])[::-1]
                n=arange(len(f))+1
        #plt.plot(f,n)
        
        
                wv=where(f>=ff_REAL[6])
                pvalue=stats.ks_2samp(ff_REAL[0:7],f[wv]).pvalue
                if plot_rea=='True':
                    plt.plot(f,n,'-',alpha=.01,color=color[j])
                
                f_array.append(list(f.data))
            
            mi=m(f_array)[0]
            ma=m(f_array)[1]
            x=logspace(log10(mi),log10(ff_REAL.max()),100)
            inter_x.append(interp(x,f_array[k][::-1],n[::-1]))
            inter_reali.append(interp(ff_REAL[0:7][::-1],f_array[k][::-1],n[::-1]))
            
            prova_fit.append(x)
            
            mean=array(inter_reali).mean(axis=0)
            mean_sort=sort(mean)[::-1]
            if plot_mean=='True':
                p1=plt.plot(x,sort(array(inter_x).mean(axis=0))[::-1],color=color[j])
                plt.plot(ff_REAL[0:7],mean_sort[::-1],'-d',markeredgecolor='white',markerfacecolor=color[j],alpha=1)
            
            pvalue=stats.ks_2samp(nn_REAL[0:7],mean_sort).pvalue
            print(round(pvalue,7))
            x1=[]
            f_array=[]
            inter_reali=[]
            inter_x=[]
        
        
        flux=cloud['Flux_int_100MeV']
        ff=sort(flux)[::-1]
        nn=arange(0,len(ff),1)+1
        plt.plot(ff,nn,color='tab:red')
        
        w=where((real['Type']=='INT')*(real['SNR_Name']!='W28south'))
        flux_real=real['Flux100MeV'][w]
        ff_real=sort(flux_real)[::-1]
        nn_real=arange(0,len(ff_real),1)+1
        if plot_real=='True':
            plt.plot(ff_real,nn_real,'o-',color='black')
            plt.fill_between(ff_REAL,nn_REAL-sqrt(nn_REAL),nn_REAL+sqrt(nn_REAL),color='gray',alpha=.2)
        
        
        plt.xlabel(r'Flux $erg \, cm^2 s{^-1}$')
        plt.ylabel(r'N(>F)')
        plt.title('SNR-MC LogN-LogS',fontsize='large')
        plt.loglog()
    plt.show()


    

    return len(icloud)















