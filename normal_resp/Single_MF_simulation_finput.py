import numpy as np
'''
import sys
sys.path.append('../')
from synapses_and_connectivity.syn_and_connec_library import get_connectivity_and_synapses_matrix
'''
from syn_and_connec_library import get_connectivity_and_synapses_matrix
from cell_library import get_neuron_params
from theoretical_tools import get_fluct_regime_varsup, pseq_params,TF_my_templateup
from scipy.special import erf
from numpy import arange
from numpy import meshgrid
from scipy.stats import norm
import random
import math



def build_up_differential_operator_first_order(TF1, TF2, NRN1,NRN2,NTWK, T=5e-3):
    
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    params = get_neuron_params(NRN1, SI_units=True)
    reformat_syn_parameters(params, M)
    
    a, b, tauw = params['a'],\
        params['b'], params['tauw']
    Qe, Te, Ee = params['Qe'], params['Te'], params['Ee']
    Qi, Ti, Ei = params['Qi'], params['Ti'], params['Ei']
    Gl, Cm , El = params['Gl'], params['Cm'] , params['El']
    pconnec,Ntot,gei,ext_drive=params['pconnec'], params['Ntot'] , params['gei'],M[0,0]['ext_drive']




    def A0(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
        return 1./T*(TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[2])-V[0])
    
    def A1(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0):
        return 1./T*(TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[3])-V[1])
    
    def A2(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
        
        fe = (V[0]+exc_aff+pure_exc_aff)*(1.-gei)*pconnec*Ntot # default is 1 !!
        fi = V[1]*gei*pconnec*Ntot
        muGe, muGi = Qe*Te*fe, Qi*Ti*fi
        muG = Gl+muGe+muGi
        muV = (muGe*Ee+muGi*Ei+Gl*El-V[2])/muG
        

        
        return (-V[2]/tauw+(b)*V[0]+a*(muV-El)/tauw)
    

    
    def A3(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
        return  (-V[3]/1.0+0.*V[1])
    
    def Diff_OP(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0):
        return np.array([A0(V, exc_aff=exc_aff,inh_aff=inh_aff, pure_exc_aff=pure_exc_aff),\
                         A1(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract),\
                         A2(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff),\
                         A3(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff)])
    return Diff_OP


def build_up_differential_operator(TF1, TF2, NRN1,NRN2,NTWK,\
                                   Ne=8000, Ni=2000, T=5e-3):

    
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    params = get_neuron_params(NRN1, SI_units=True)
    reformat_syn_parameters(params, M)
    
    a, b, tauw = params['a'],\
        params['b'], params['tauw']
    Qe, Te, Ee = params['Qe'], params['Te'], params['Ee']
    Qi, Ti, Ei = params['Qi'], params['Ti'], params['Ei']
    Gl, Cm , El = params['Gl'], params['Cm'] , params['El']
    pconnec,Ntot,gei,ext_drive=params['pconnec'], params['Ntot'] , params['gei'],M[0,0]['ext_drive']



    def A0(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
    
        return 1./T*(\
                 .5*V[2]*diff2_fe_fe(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 .5*V[3]*diff2_fe_fi(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 .5*V[3]*diff2_fi_fe(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 .5*V[4]*diff2_fi_fi(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])-V[0])
        
    def A1(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0):
                     
        return 1./T*(\
                  .5*V[2]*diff2_fe_fe(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  .5*V[3]*diff2_fe_fi(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  .5*V[3]*diff2_fi_fe(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  .5*V[4]*diff2_fi_fi(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])-V[1])

    def A2(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
    
        return 1./T*(\
                 1./Ne*TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])*(1./T-TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5]))+\
                 (TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])-V[0])**2+\
                 2.*V[2]*diff_fe(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 2.*V[3]*diff_fi(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                 -2.*V[2])
        
    def A3(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0): # mu, nu = e,i, then lbd = e then i
                     
        return 1./T*(\
                  (TF1(V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])-V[0])*(TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])-V[1])+\
                  V[2]*diff_fe(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  V[3]*diff_fe(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                  V[3]*diff_fi(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                  V[4]*diff_fi(TF1, V[0]+exc_aff+pure_exc_aff, V[1]+inh_aff,V[5])+\
                  -2.*V[3])

    def A4(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0):
    
        return 1./T*(\
                 1./Ni*TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])*(1./T-TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6]))+\
                 (TF2(V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])-V[1])**2+\
                 2.*V[3]*diff_fe(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                 2.*V[4]*diff_fi(TF2, V[0]+exc_aff+inh_fract*pure_exc_aff, V[1]+inh_aff,V[6])+\
                 -2.*V[4])
        
    def A5(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
        fe = (V[0]+exc_aff+pure_exc_aff)*(1.-gei)*pconnec*Ntot # default is 1 !!
        fi = V[1]*gei*pconnec*Ntot
        muGe, muGi = Qe*Te*fe, Qi*Ti*fi
        muG = Gl+muGe+muGi
        muV = (muGe*Ee+muGi*Ei+Gl*El-V[5])/muG
                                         
  
        return (-V[5]/tauw+(b)*V[0]+a*(muV-El)/tauw)





    def A6(V, exc_aff=0, inh_aff=0, pure_exc_aff=0):
        return  (-V[6]/1.0+0.*V[1])
    
    def Diff_OP(V, exc_aff=0, inh_aff=0, pure_exc_aff=0,inh_fract=0):
        return np.array([A0(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff),\
                         A1(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract),\
                         A2(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff),\
                         A3(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract),\
                         A4(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract),\
                         A5(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff),\
                         A6(V, exc_aff=exc_aff, inh_aff=inh_aff, pure_exc_aff=pure_exc_aff)])
    return Diff_OP


















import matplotlib.pylab as plt


def heaviside(x):
    return 0.5*(1+np.sign(x))

def sinewave(t, t0, freq, amplitude, phase=0):
    return amplitude*(1-np.cos(2.*np.pi*freq*(t-t0)+phase))*heaviside(t-t0)/2




def double_gaussian(t, t0, T1, T2, amplitude):
    
    return amplitude*(\
                      np.exp(-(t-t0)**2/2./T1**2)*heaviside(-(t-t0))+\
                      np.exp(-(t-t0)**2/2./T2**2)*heaviside(t-t0))








def load_transfer_functions(NRN1, NRN2, NTWK):
    """
        returns the two transfer functions of the mean field model
        """
    
    # NTWK
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    
    # NRN1
    params1 = get_neuron_params(NRN1, SI_units=True)
    reformat_syn_parameters(params1, M)
    try:
        
        
        P1 = np.load('data/RS-cell_CONFIG1_fit.npy')
        
        
        params1['P'] = P1
        def TF1(fe, fi,XX):
            return TF_my_templateup(fe, fi,XX, *pseq_params(params1))
    
    
    
    
    
    
    except IOError:
        print('=======================================================')
        print('=====  fit for NRN1 not available  ====================')
        print('=======================================================')

    # NRN1
    params2 = get_neuron_params(NRN2, SI_units=True)
    reformat_syn_parameters(params2, M)
    try:
        
        P2 = np.load('data/FS-cell_CONFIG1_fit.npy')
        
        
        
        params2['P'] = P2
        def TF2(fe, fi,XX):
            return TF_my_templateup(fe, fi,XX, *pseq_params(params2))

    except IOError:
        print('=======================================================')
        print('=====  fit for NRN2 not available  ====================')
        print('=======================================================')
    
    return TF1, TF2








def reformat_syn_parameters(params, M):
    """
        valid only of no synaptic differences between excitation and inhibition
        """
    
    params['Qe'], params['Te'], params['Ee'] = M[0,0]['Q'], M[0,0]['Tsyn'], M[0,0]['Erev']
    params['Qi'], params['Ti'], params['Ei'] = M[1,1]['Q'], M[1,1]['Tsyn'], M[1,1]['Erev']
    params['pconnec'] = M[0,0]['p_conn']
    params['Ntot'], params['gei'] = M[0,0]['Ntot'], M[0,0]['gei']




##### Derivatives taken numerically,
## to be implemented analitically ! not hard...

def diff_fe(TF, fe, fi,XX, df=1e-5):
    return (TF(fe+df/2., fi,XX)-TF(fe-df/2.,fi,XX))/df

def diff_fi(TF, fe, fi,XX, df=1e-5):
    return (TF(fe, fi+df/2.,XX)-TF(fe, fi-df/2.,XX))/df

def diff2_fe_fe(TF, fe, fi,XX, df=1e-5):
    return (diff_fe(TF, fe+df/2., fi,XX)-diff_fe(TF,fe-df/2.,fi,XX))/df

def diff2_fi_fe(TF, fe, fi,XX, df=1e-5):
    return (diff_fi(TF, fe+df/2., fi,XX)-diff_fi(TF,fe-df/2.,fi,XX))/df

def diff2_fe_fi(TF, fe, fi,XX, df=1e-5):
    return (diff_fe(TF, fe, fi+df/2.,XX)-diff_fe(TF,fe, fi-df/2.,XX))/df

def diff2_fi_fi(TF, fe, fi,XX, df=1e-5):
    return (diff_fi(TF, fe, fi+df/2.,XX)-diff_fi(TF,fe, fi-df/2.,XX))/df




def run_mean_field_2order(NRN1, NRN2, NTWK,T=5e-3, dt=1e-4, tstop=2):
    
    freq=40.
    amp = 0
    t0 = 1.
    
    tstop= t0+10./freq
    tstop=2
    T=T
    dt=dt
    tstop=tstop
    
    
    #take network parameters
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    params = get_neuron_params(NRN1, SI_units=True)
    reformat_syn_parameters(params, M)
    
    a, b, tauw = params['a'],\
    params['b'], params['tauw']
    Qe, Te, Ee = params['Qe'], params['Te'], params['Ee']
    Qi, Ti, Ei = params['Qi'], params['Ti'], params['Ei']
    Gl, Cm , El = params['Gl'], params['Cm'] , params['El']
    pconnec,Ntot,gei,ext_drive=params['pconnec'], params['Ntot'] , params['gei'],M[0,0]['ext_drive']
    
    #IMPORTANT PARAMETERS
    
    
    #to change network size
    Ntot=10000
    
    
    extinp=2.5

    sigma_0=10.
    
    filesave="MF_UD_amplitude_1,9_2nd_order_ext_3_N_10000.npy"
    
    
    Ne=Ntot*(1-gei)
    Ni=Ntot*gei
    
    TF1, TF2 = load_transfer_functions(NRN1, NRN2, NTWK)

    t = np.arange(int(tstop/dt))*dt
    fe=0*t
    fi=0*t
    ww=0*t
    v2vec=0*t
    v3vec=0*t
    v4vec=0*t

    extinpnuovo=sinewave(t, t0, freq, amp)
   

    


    
    
    ornstein=0
    ornsteinin2=0
    extinpnoisex1=0
    extinpnoisex=0
    ornsteinA=0.0002 #time decay ornstein process
    
    
    
    
    
    
    #initial conditions 
    fecont=2.8
    ficont=50
    v2=0.5
    v3=0.5
    v4=0.5
    wcont=fecont*b*tauw
##    fecont= 1.4907357840560944
##    ficont=9.198083797700072
##    v2=0.15129606503776358
##    v3=0.3025868361216362
##    v4=0.7632335165277322
##    wcont=7.454102191637575e-11
    
    
    for i in range(len(t)):
        onlyfirstorder=1.
     
       
   
           
        #generate the wiener process sigma_0


                           
        #avoid external noise
        extinpnoisex=0.
        extinpnoisex1=0.
        
        fecontold=fecont
        ficontold=ficont
        

        azzera=1. #you can  avoid noise in covariance time evolution azzera=0
        '''
        fecont+=dt/T*(TF1(fecont+extinp, ficont,wcont)-fecont)
        ficont+=dt/T*(TF2(fecont+extinp, ficont,0)-ficont)
        '''
        
        fecont+=dt/T*(\
                      (.5*v2*diff2_fe_fe(TF1, fecont+extinp+azzera*extinpnoisex+extinpnuovo[i], ficont,wcont)+\
                      .5*v3*diff2_fe_fi(TF1, fecont+extinp+azzera*extinpnoisex+extinpnuovo[i], ficont,wcont)+\
                      .5*v3*diff2_fi_fe(TF1, fecont+extinp+azzera*extinpnoisex+extinpnuovo[i], ficont,wcont)+\
                      .5*v4*diff2_fi_fi(TF1, fecont+extinp+azzera*extinpnoisex+extinpnuovo[i], ficont,wcont))*onlyfirstorder+\
                      TF1(fecont+extinp+extinpnoisex+extinpnuovo[i], ficont,wcont)-fecont)
                      
        ficont+=dt/T*(\
                      (.5*v2*diff2_fe_fe(TF2, fecontold+extinp+azzera*extinpnoisex1, ficont,0.)+\
                      .5*v3*diff2_fe_fi(TF2, fecontold+extinp+azzera*extinpnoisex1, ficont,0.)+\
                      .5*v3*diff2_fi_fe(TF2,fecontold+extinp+azzera*extinpnoisex1, ficont,0.)+\
                      .5*v4*diff2_fi_fi(TF2, fecontold+extinp+azzera*extinpnoisex1, ficont,0.))*onlyfirstorder+\
                      TF2(fecontold+extinp+extinpnoisex1, ficont,0.)-ficont)
                      
          
        


        
        v2+=dt/T*(\
                 1./Ne*TF1(fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)*(1./T-TF1(fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont))+\
                 (TF1(fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)-fecontold)**2+\
                 2.*v2*diff_fe(TF1, fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)+\
                 2.*v3*diff_fi(TF1,fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)+\
                 -2.*v2)
           
        v3+=dt/T*(\
                (TF1(fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)-fecontold)*(TF2( fecontold+extinp+extinpnoisex1, ficontold,0.)-ficontold)+\
                v2*diff_fe(TF2,  fecontold+extinp+extinpnoisex1, ficontold,0.)+\
                v3*diff_fe(TF1, fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)+\
                v3*diff_fi(TF2,  fecontold+extinp+extinpnoisex1, ficontold,0.)+\
                v4*diff_fi(TF1, fecontold+extinp+extinpnoisex+extinpnuovo[i], ficontold,wcont)+\
                -2.*v3)
           
        v4+=dt/T*(\
                  1./Ni*TF2(fecontold+extinp+extinpnoisex1, ficontold,0.)*(1./T-TF2(fecontold+extinp+extinpnoisex1, ficontold,0.))+\
                  (TF2(fecontold+extinp+extinpnoisex1, ficontold,0.)-ficontold)**2+\
                  2.*v3*diff_fe(TF2, fecontold+extinp+extinpnoisex1, ficontold,0.)+\
                  2.*v4*diff_fi(TF2,fecontold+extinp+extinpnoisex1, ficontold,0.)+\
                  -2.*v4)


        
        

        #to ensure non-negative firing rate due to noise in the down state
        
        if(fecont<0.):
            fecont=1.e-8
        if(ficont<0.):
            ficont=1.e-2
        
        
        feA = (fecontold+extinp)*(1.-gei)*pconnec*Ntot # default is 1 !!
        fiA = ficontold*gei*pconnec*Ntot
        muGe, muGi = Qe*Te*feA, Qi*Ti*fiA
        muG = Gl+muGe+muGi
        muV = (muGe*Ee+muGi*Ei+Gl*El-wcont)/muG
        
        wcont+=dt*( -wcont/tauw+(b)*fecontold+a*(muV-El)/tauw)
    
        
       
        v2vec[i]=v2
        v3vec[i]=v3
        v4vec[i]=v4
        fe[i]=fecont
        fi[i]=ficont
        ww[i]=wcont
        print(t[i],fecont,ficont,v2,v3,v4,wcont)
 
 
 
 
    


    
    




    print("FIXED POINT",fecont,ficont,np.sqrt(v2),np.sqrt(v3),np.sqrt(v4),wcont)
    print("FIXED POINT",fecont,ficont,v2,v3,v4,wcont)


    #np.save(filesave,[t,fe,fi,np.sqrt(v4vec),np.sqrt(v2vec) ])
    plt.plot(t,fe,t,fi)
    plt.plot(t,v2vec,t,v3vec)
    plt.plot(t,v4vec,t,ww)
    plt.show()
    



def run_mean_field_2order_secondway(NRN1, NRN2, NTWK,T=5e-3, dt=1e-4, tstop=2):



    amp = 5
    t0 = 2.
    T1=0.01
    T2=0.1
    

    
    dt=dt
    tstop=tstop
    
    
    #take network parameters
    M = get_connectivity_and_synapses_matrix(NTWK, SI_units=True)
    params = get_neuron_params(NRN1, SI_units=True)
    reformat_syn_parameters(params, M)
    
    a, b, tauw = params['a'],\
    params['b'], params['tauw']
    Qe, Te, Ee = params['Qe'], params['Te'], params['Ee']
    Qi, Ti, Ei = params['Qi'], params['Ti'], params['Ei']
    Gl, Cm , El = params['Gl'], params['Cm'] , params['El']
    pconnec,Ntot,gei,ext_drive=params['pconnec'], params['Ntot'] , params['gei'],M[0,0]['ext_drive']
    
    
    
    
    #to change network size
    Ntot=10000
    
    
    extinp=2.5
    
    sigma_0=10.
    
    filesave='MF_Gaussian_response.npy'
    
    
    Ne=Ntot*(1-gei)
    Ni=Ntot*gei
    
    TF1, TF2 = load_transfer_functions(NRN1, NRN2, NTWK)
    
    
   
    
    t = np.arange(int(tstop/dt))*dt
    fe=0*t
    fi=0*t
    ww=0*t
    v2vec=0*t
    v3vec=0*t
    v4vec=0*t
    
    extinpnuovo=double_gaussian(t, t0, T1, T2, amp)
    
 
    
    
    
    
    
    ornstein=0
    ornsteinin2=0
    extinpnoisex1=0
    extinpnoisex=0
    ornsteinA=0.0002 #time decay ornstein process
    
    
    
    
    
    
    #initial conditions
    fecont=.8
    ficont=50
    v2=0.5
    v3=0.5
    v4=0.5
    wcont=fecont*b*tauw





    def dX_dt_scalar(X, t=0):
        exc_aff = extinp*t/t
        pure_exc_aff = extinpnuovo
        return build_up_differential_operator(TF1, TF2,NRN1, NRN2, NTWK,  Ne=Ne, Ni=Ni, T=T)(X,\
                                                                                             exc_aff=exc_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract)
    
    
    X0 = [1., 30, .5, .5, .5,1.e-10,0.]
    X = np.zeros((len(t), len(X0)))
    for i in range(0,6):
        X[0][i]=X0[i]



    inh_fract=1.
    for i in range(len(t)-1):


        stdvexc=np.sqrt(X[i,0] )
        stdvinh=np.sqrt(X[i,1] )
        
        


        
        
        
        exc_aff = extinp
        pure_exc_aff = (extinpnuovo[i])
        


        
        
        X[i+1,:] = X[i,:] + (t[1]-t[0])*build_up_differential_operator(TF1, TF2,NRN1, NRN2, NTWK, \
                                                                       Ne=Ne, Ni=Ni)(X[i,:], exc_aff=exc_aff, pure_exc_aff=pure_exc_aff,inh_fract=inh_fract)


        X[i+1,6] = 0
    
    fe, fi = X[:,0], X[:,1]
    sfe, sfei, sfi = [np.sqrt(X[:,i]) for i in range(2,5)]
    XXe,XXi= X[:,5], X[:,6]



    sig=0.8*fe+0.2*fi


    plt.plot(t,fe)
    plt.show()
    np.save(filesave,[t,fe, fi, sfe, sfei, sfi, XXe, XXi])






run_mean_field_2order('RS-cellbis', 'FS-cell', 'CONFIG1',T=5e-3, dt=5e-4, tstop=4)


