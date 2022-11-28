
import numpy as np
import pandas as pd
import pickle
import time
from scipy.stats import lognorm
from scipy.stats import gamma
from os import path

# VIRAL LOAD FUNCTIONS


def VL(t,T3,TP,YP,T6):
    '''
    VL(t,T3,TP,YP,T6)
        Computes the scalar value of a hinge function with control points:
        (T3, 3)
        (TP, YP)
        (T6, 6)
    Returns zero whenever (1) t<T3 or (2) hinge(t) is negative.
    '''
    if t < T3:
        return 0
    if t < TP:
        return (t-T3)*(YP-3)/(TP-T3)+3
    return np.max([(t-TP)*(6-YP)/(T6-TP)+YP,0])

    
def get_trajectory(is_symptomatic=False,full=False):
    '''
    get_trajectory(is_symptomatic=False,full=False)
        Stochastically draws a viral load trajectory by randomly drawing control points.
    
    is_symptomatic allows separate models for symptomatic and asymptomatic trajectories which
    the review by Cevik et al 2020 has shown to be different. (Asymptomatic clearance is faster)
    
    When full==False, returns a viral load trajectory evaluated on int days 0,1,2,...,21,
    where day 0 is the day of exposure. 

    When full==True, returns the integer-evaluated trajectory, as well as:
    v_fine the same trajectory evaluated at values between 0 and 28, with step size 0.01
    t_fine the steps at which the trajectory was evaluated
    t_3, t_peak, v_peak, t_6, the four stochastically drawn control values

    Use full==True for plotting; use full==False for more rapid simulation.
    '''
    if is_symptomatic==True:
        return get_symptomatic_trajectory(full=full)
    else:
        return get_asymptomatic_trajectory(full=full)


def get_symptomatic_trajectory(full=False):
    # Draw control points. 
    t_peak = gamma.rvs(1.5,loc=0.5)
    while t_peak>3 :
        t_peak = gamma.rvs(1.5,loc=0.5)
    v_peak = np.random.uniform(6,9)
    t_6 = np.random.uniform(3,4)
    t_3 = np.random.random()+2
    t_symptoms = np.random.uniform(1,2)
    # Compute t and v vectors and return 
    t = np.arange(21)
    v = np.array([VL(x,t_3,t_peak+t_symptoms+t_3,v_peak,t_6+t_3+t_peak+t_symptoms) for x in t])
    if full==True:
        tfine = np.arange(0,21,0.01)
        vfine = np.array([VL(x,t_3,t_peak+t_symptoms+t_3,v_peak,t_6+t_3+t_peak+t_symptoms) for x in tfine])
        return v,vfine,tfine,t_3,t_peak,v_peak,t_6,t_symptoms
    return v,int(np.round(t_symptoms+t_3))


def get_asymptomatic_trajectory(full=False):    
    # Draw control points. 
    t_peak = gamma.rvs(1.5,loc=0.5)
    while t_peak>5 :
        t_peak = gamma.rvs(1.5,loc=0.5)
    v_peak = np.random.uniform(6,9)
    t_6 = np.random.uniform(3,4)
    t_3 = np.random.random()+2
    # Compute t and v vectors and return 
    t = np.arange(21)
    v = np.array([VL(x,t_3,t_peak+t_3,v_peak,t_6+t_3+t_peak) for x in t])
    if full==True:
        tfine = np.arange(0,21,0.01)
        vfine = np.array([VL(x,t_3,t_peak+t_3,v_peak,t_6+t_3+t_peak) for x in tfine])
        return v,vfine,tfine,t_3,t_peak,v_peak,t_6
    return v


# INFECTIOUSNESS FUNCTIONS
def logproportional(v,cutoff=6):
    '''
    logproportional(v,cutoff=6)
        returns the infectiousness of a viral load v
        for the logproportional model
        (Manuscript default)
    '''
    if v < cutoff:
        return 0
    else:
        return v-cutoff


# 计算感染期间的平均感染性
def compute_factor_to_calibrate_R0_SQ(infectiousness,asymptomatic,cutoff):
    '''
    compute_factor_to_calibrate_R0_SQ(infectiousness,asymptomatic,cutoff)
        infectiousness: a function handle {proportional,logproportional,threshold}
        asymptomatic: the fraction of individuals who are asymptomatic [0,1]
        cutoff: the minimum value of log10 viral load for infectiousness. (Manuscript default: 6)

    Returns a constant mean_infectiousness which is used to scale absolute infectiousness in simulation
    '''
    total_infectiousness = 0
    n_draws = 10000
    for i in range(n_draws):
        if np.random.random() < asymptomatic:
            VL = get_trajectory(is_symptomatic=False)
            IN = [infectiousness(x,cutoff=cutoff) for x in VL]
            total_infectiousness += np.sum(IN)
        else:
            VL,t_symptoms = get_trajectory(is_symptomatic=True)
            IN = [infectiousness(x,cutoff=cutoff) for x in VL]
            total_infectiousness += np.sum(IN[:t_symptoms])
    mean_infectiousness = total_infectiousness/n_draws
    return mean_infectiousness


def SEIRsimulation(N,D,L,infectiousness_function,asymptomatic=0.77,results_delay=0,R0=9.5,cutoff=6,I_init=10,tmax=90):
    c = compute_factor_to_calibrate_R0_SQ(infectiousness=logproportional,asymptomatic=0.77,cutoff=6) 
    k = R0/(c*(N-1))  #
    was_infected_ext = np.zeros(N)
    was_infected_int = np.zeros(N)
    is_symptomatic = np.random.binomial(1,1-asymptomatic,size=N) 
    t_symptoms = np.zeros(N)
    viral_loads = {}
    infectious_loads = {}  
    for i in range(N):
        if is_symptomatic[i]:
            viral_loads[i],t_symptoms[i] = get_trajectory(is_symptomatic=True)
        else:
            viral_loads[i] = get_trajectory(is_symptomatic=False)
        infectious_loads[i] = [infectiousness_function(x,cutoff=cutoff)*k for x in viral_loads[i]]  


    test_schedule = np.random.choice(D,size=N)      
    ###赋初值
    S = set(np.arange(N))
    I = set()     
    R = set()
    Q = set()
    SQ = set()
    St = []
    It = []
    Rt = []
    Qt = []
    SQt = []
    infection_day = np.zeros(N,dtype=int)  
    days_till_results = -1*np.ones(N)

    for t in range(90):
        if t==0:
            initial_infections = np.random.choice(list(S),size=I_init,replace=False)  
            was_infected_int[initial_infections] = 1    
            S = S - set(initial_infections)  
            I = I.union(set(initial_infections)) 

            #确定初始10例感染的感染者的感染时间（假定10个体病毒载量已超过10^3）
            time=np.arange(0,21,1)
            for i in I:
            #print(i) 
                data=pd.DataFrame({'viral_load':viral_loads[i],'time':time})
                VL_location=np.where(data["viral_load"]>=6)
                #print(VL_location)
                VL_loc=VL_location[0].tolist()
                loc=np.random.choice(VL_loc,size=1,replace=False)
                infection_day[i]=data['time'][loc] 

         # Internal infections 内部感染
        infectiousnesses = [infectious_loads[i][infection_day[i]] for i in I]
        p_int_infection =1 - np.prod(1-np.array(infectiousnesses)) #内部感染率  
        #print(len(I),p_int_infection)
        #print(infectiousnesses)
        int_infections = np.random.choice(list(S),np.random.binomial(len(S),p_int_infection))
        #np.random.binomial(len(S),p_int_infection)
        was_infected_int[int_infections] = 1
        S = S - set(int_infections)
        I = I.union(set(int_infections))


        # Test
        tested = np.where((test_schedule + t) % D==0)[0] 
        for i in tested:
                if (i in Q) or (i in R) or (i in SQ):
                    continue
                if days_till_results[i] > -1:  
                    continue
                if viral_loads[i][infection_day[i]] > L: 
                    days_till_results[i] = results_delay  

        # Isolate
        for i in I:
            if days_till_results[i]==0:
                I = I - set([i])
                Q = Q.union(set([i]))

        # Self Quarantine
        for i in I:
            if is_symptomatic[i]==1:
                if infection_day[i]==t_symptoms[i]:
                        I = I-set([i])
                        SQ = SQ.union(set([i]))




         # Update all infection days  
        for i in I:
            if (infection_day[i] == 20) or ((viral_loads[i][infection_day[i]]<6) and (infection_day[i]>7)): 
                I = I - set([i])
                R = R.union(set([i]))
            infection_day[i] += 1
            #print(infection_day[i])
        for q in Q:
            if (infection_day[q] == 20) or ((viral_loads[q][infection_day[q]]<6) and (infection_day[q]>7)):
                Q = Q - set([q])
                R = R.union(set([q]))
            infection_day[q] += 1
            #print(infection_day[q])
        for q in SQ:
            if (infection_day[q] == 20) or ((viral_loads[q][infection_day[q]]<6) and (infection_day[q]>7)):
                SQ = SQ - set([q])
                R = R.union(set([q]))
            infection_day[q] += 1
           # print(infection_day)

        # Update the results delays:
        for i in range(N):
            if days_till_results[i] > -1:
                days_till_results[i] = days_till_results[i] - 1
            #print(days_till_results)

        St.append(len(S))
        It.append(len(I))
        Rt.append(len(R))
        Qt.append(len(Q))
        SQt.append(len(SQ))

    return St,It,Rt,Qt,SQt,sum(It)

