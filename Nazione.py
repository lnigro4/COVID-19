import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sc
import pandas as pd
import math
import json

from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.integrate import solve_ivp
from datetime import date, datetime


'''
dati presi da https://github.com/pcm-dpc/COVID-19

POI RIPETERE ANCHE PER LOGISTICHE/DIST CANONICA CON TAU AL POSTO DI K
'''

#-------------------------precompilatore----------------------------------------
def sir_model(I, R, N):
    S=N-np.asarray(R)-np.asarray(I)
    I=np.asarray(I)
    R=np.asarray(R)
    dS=np.zeros_like(S)
    dI=np.zeros_like(I)
    dR=np.zeros_like(R)
    for i in range(1, len(R)-1):
        dS[i]=S[i]-S[i-1]
        dI[i]=I[i]-I[i-1]
        dR[i]=R[i]-R[i-1]
    beta1 = -dS*N/(S*I)
    gamma1= dR/I
    beta2 = (dI+gamma1*I)*N/(S*I)
    gamma2= (beta1*S*I/N-dI)/(I)
    R01 = beta1[7:-1]/gamma1[7:-1]
    R02 = beta2[7:-1]/gamma1[7:-1]
    R03 = beta1[7:-1]/gamma2[7:-1]
    if R01.all()==R02.all() and R01.all()==R03.all():
        print('eq diff risolta')
        return R01
    else:
        print('eq diff non risolta in modo esatto')
        return R01, R02, R03

def seir_model(E, I, R, N):
    S=N-np.asarray(R)-np.asarray(I)-np.asarray(E)
    I=np.asarray(I)
    R=np.asarray(R)
    E=np.asarray(E)
    dS=np.zeros_like(S)
    dE=np.zeros_like(E)
    dI=np.zeros_like(I)
    dR=np.zeros_like(R)
    for i in range(1, len(R)-1):
        dS[i]=S[i]-S[i-1]
        dE[i]=E[i]-E[i-1]
        dI[i]=I[i]-I[i-1]
        dR[i]=R[i]-R[i-1]
    mu = -(dR+dI+dE+dS)/(I+R+S)
    eps = (dR+dI+mu*(I+R))/E
    gammu = (dR+mu*(I+R))/I
    beta = N*(dE+dR+dI+mu*(I+R))/(S*I)
    R0=beta[7:-1]*eps[7:-1]/(gammu[7:-1]*(eps[7:-1]+mu[7:-1]))
    return R0

def gaussian(x, a, k):
    return a*np.exp(k*x)

def linear(x, m, q):
    return m*x+q

inf_time = 5
exp_time = 17
N=60483973

#-------------------------lettura da file---------------------------------------

with open('COVID-19/dati-json/dpc-covid19-ita-andamento-nazionale.json') as nazione:
    dati_nazione = json.load(nazione)
    

#print(json.dumps(dati_nazione, indent=4, sort_keys=False))

#-------------------------analisi dati------------------------------------------

oggi = datetime.date(datetime.now())

n_giorni = (lambda delta: delta.days)(oggi - date(2020, 2, 24))

tot_casi = []
giorni = []
new_casi = []
x_giorni = []
deceduti = []
guariti = []
I = []
E = []
for d in range(0, n_giorni-1):
    tot_casi.append(dati_nazione[d]['totale_casi'])
    new_casi.append(dati_nazione[d]['nuovi_positivi'])
    deceduti.append(dati_nazione[d]['deceduti'])
    guariti.append(dati_nazione[d]['dimessi_guariti'])
    R = np.asarray(deceduti)
    I.append(dati_nazione[d]['totale_positivi'])
    E.append(dati_nazione[d]['isolamento_domiciliare'])
    giorni.append(dati_nazione[d]['data'][:10])
    x_giorni.append(d)

x=[]
y=[]
k_tot=[]
R0_tot = []
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y.append(tot_casi[j])
    popt, pcov = curve_fit(gaussian, x, y)
    k_tot.append(popt[1])
    R0_tot.append(1+popt[1]*(exp_time+inf_time)+(popt[1]**2)*(exp_time)*(inf_time))
    x.clear()
    y.clear()

k_new=[]
R0_new=[]
for i in range(0, len(tot_casi)):
    k_new.append(new_casi[i]/tot_casi[i])
    R0_new.append(1+(new_casi[i]/tot_casi[i])*(exp_time+inf_time)+((new_casi[i]/tot_casi[i])**2)*(exp_time)*(inf_time))

popt_1, pcov_1 = curve_fit(linear, tot_casi[-14:], k_new[-14:])
pred1 = -popt_1[1]/popt_1[0]
try:
    popt_2, pcov_2 = curve_fit(linear, tot_casi[n_giorni - int(inf_time/2)-14: n_giorni - int(inf_time/2)], k_tot[-14:])
    pred2 = -popt_2[1]/popt_2[0]
except:
    print("controlla il range di fit")

#-------------------------opzioni estetiche plot--------------------------------
    
fig = plt.figure(figsize=(15, 10))

ax1 = fig.add_subplot(221)
c1 = 'tab:blue'
ax1.grid()
ax1.xaxis.grid(True, which='minor', linestyle=':')
ax1.yaxis.grid(True, which='minor', linestyle=':')
ax1.set_title('Covid-19 Time Evolution - Italy')
ax1.set_ylabel('Total infected')
ax1.tick_params(axis='y')
ax1.set_yscale('log')
plt.xticks(rotation=45)
maj_loc = mpl.ticker.MultipleLocator(base=3.0)
min_loc = mpl.ticker.MultipleLocator(base=1.0)
ax1.xaxis.set_major_locator(maj_loc)
ax1.xaxis.set_minor_locator(min_loc)

ax2 = fig.add_subplot(222)
ax2.grid()
c2 = 'tab:red'
ax2.xaxis.grid(True, which='minor', linestyle=':')
ax2.yaxis.grid(True, which='minor', linestyle=':')
ax2.set_title('Growth k factor')
ax2.set_ylabel('k factor')
ax2.set_xlabel('Total infected')
ax2.tick_params(axis='y')
#ax2.set_xscale('log')
plt.xticks(rotation=45)
ax2.set_ylim(0,1)
ax2.set_xlim(10000, pred1)
ax2.set_ylim(0,.2)

ax3 = fig.add_subplot(223)
ax3.grid()
ax3.xaxis.grid(True, which='minor', linestyle=':')
ax3.yaxis.grid(True, which='minor', linestyle=':')
ax3.set_title('SEIR model (from k factor)')
ax3.set_ylabel('$R_0$ index')
ax3.set_xlabel('Total infected')
ax3.tick_params(axis='y')
#ax3.set_xscale('log')
plt.xticks(rotation=45)
ax3.set_xlim(10000, pred1)
ax3.set_ylim(0,9)

'''seir vs sir models
ax4 = fig.add_subplot(224)
ax4.grid()
ax4.xaxis.grid(True, which='minor', linestyle=':')
ax4.yaxis.grid(True, which='minor', linestyle=':')
ax4.set_title('SEIR vs SIR model (analitical)')
ax4.set_ylabel('$R_0$ index')
ax4.set_xlabel('Total infected')
ax4.tick_params(axis='y')
#ax4.set_xscale('log')
plt.xticks(rotation=45)
ax4.set_xlim(10000, pred1)
ax4.set_ylim(0,16)
'''

ax5 = fig.add_subplot(224)
ax5.grid()
ax5.xaxis.grid(True, which='minor', linestyle=':')
ax5.yaxis.grid(True, which='minor', linestyle=':')
ax5.set_title('Deaths vs Healed')
ax5.set_ylabel('number of people')
ax5.tick_params(axis='y')
#ax5.set_yscale('log')
plt.xticks(rotation=45)
ax5.xaxis.set_major_locator(maj_loc)
ax5.xaxis.set_minor_locator(min_loc)


plt.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.95, hspace=0.35)

#-------------------------plot dei dati-----------------------------------------

ax1.plot(giorni, tot_casi, 'ko--', label=r'Protezione Civile data - Totale infetti')
ax1.legend()

ax2.plot(tot_casi[int(inf_time/2)+1: n_giorni - int(inf_time/2)], k_tot[1:],
         's--', color=c2, label='based on total infected')
ax2.plot(tot_casi, k_new,
         's--', color=c1, label='based on new daily infected')
ax2.plot(tot_casi[-14:], linear(np.asarray(tot_casi[-14:]), popt_1[0], popt_1[1]),
         '--', color='tab:grey')
ax2.plot(tot_casi[-2-14:-2], linear(
    np.asarray(tot_casi[-2-14:-2]), popt_2[0], popt_2[1]),'--', color='tab:grey')
ax2.plot(pred1, 0, 'ko')
ax2.plot(pred2, 0, 'ko')
ax2.annotate('predicted\n %i\n infected' %pred1, xy=(pred1, 0),
             xytext=(pred1 -20000, 0.05), color='tab:blue',
             arrowprops=dict(arrowstyle='->'))
ax2.annotate('predicted\n %i\n infected' %pred2, xy=(pred2, 0),
             xytext=(pred2-50000, 0.0075), color='tab:red',
             arrowprops=dict(arrowstyle='->'))
ax2.legend(loc='best')

ax3.plot(tot_casi[int(inf_time/2)+12: n_giorni - int(inf_time/2)], R0_tot[12:],
         's--', color=c2, label='based on total infected')
ax3.plot(tot_casi[12:], R0_new[12:],
         's--', color=c1, label='based on new daily infected')
ax3.errorbar(tot_casi[-1], R0_new[-1], yerr=0.48, color='k', fmt = '.')
ax3.legend(loc='best')

'''seir vs sir models
ax4.plot(tot_casi[7:-1], seir_model(E, I, R, N),
         's--', label=r'$R_0=\beta\epsilon/(\gamma +\mu)(\epsilon +\mu)$')
ax4.plot(tot_casi[7:-1], sir_model(I, R, N), 's--', label=r'$R_0=\beta/\gamma$')
ax4.legend(loc='best')
'''
ax5.plot(giorni, deceduti, 'ro--', label=r'Protezione Civile data - Deceduti')
ax5.plot(giorni, guariti, 'ko--', label=r'Protezione Civile data - Guariti dimessi')
ax5.legend()

fig.savefig('/Users/nigresson/Desktop/COVID19/growth k-factor comparison')




