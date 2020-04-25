import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import scipy as sc
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import json
from datetime import date
from datetime import datetime


'''
dati presi da https://github.com/pcm-dpc/COVID-19
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

def k_factor(new, tot):
    k=[]
    for i in range(len(tot)):
        if tot[i]!= 0:
            k.append(new[i]/tot[i])
        else:
            k.append(None)
    return k

def R0(k, t_e, t_i):
    R0=[]
    for i in range(len(k)):
        if k[i] == None:
            R0.append(None)
        else:
            R0.append(1+k[i]*(t_e + t_i)+(k[i]**2)*t_e*t_i - 0.2209206990276071)
    return R0

def gaussian(x, a, k):
    return a*np.exp(k*x)

def linear(x, m, q):
    return m*x+q

inf_time = 5
exp_time = 12
N=10002615

regione='Lombardia'

asse_x = 'Data'
#asse_x = 'Infetti'

#-------------------------lettura da file---------------------------------------

with open('COVID-19/dati-json/dpc-covid19-ita-regioni.json') as regioni:
    dati_regioni = json.load(regioni)

#print(json.dumps(dati_regioni, indent=4, sort_keys=False))

#-------------------------analisi dati------------------------------------------

oggi = datetime.date(datetime.now())

n_giorni = (lambda delta: delta.days)(oggi - date(2020, 2, 24))

tot = []
giorni = []
new = []
x_giorni = []
deceduti = []
guariti = []
I = []
E = []
tamponi = []
counter=0
intensiva = []
new_intensiva = []
new_deceduti = []
new_guariti = []

for d in range(len(dati_regioni)):
    if dati_regioni[d]['denominazione_regione']==regione:
        counter+=1
        tot.append(dati_regioni[d]['totale_casi'])
        new.append(dati_regioni[d]['nuovi_positivi'])
        giorni.append(dati_regioni[d]['data'][:10])
        x_giorni.append(counter)
        deceduti.append(dati_regioni[d]['deceduti'])
        guariti.append(dati_regioni[d]['dimessi_guariti'])
        intensiva.append(dati_regioni[d]['terapia_intensiva'])
        R = np.asarray(deceduti)
        I.append(dati_regioni[d]['totale_positivi'])
        E.append(dati_regioni[d]['isolamento_domiciliare'])
        tamponi.append(dati_regioni[d]['tamponi'])

new_intensiva.append(intensiva[0])
new_deceduti.append(deceduti[0])
new_guariti.append(guariti[0])
for i in range(1, n_giorni):
    new_intensiva.append(intensiva[i]-intensiva[i-1])
    new_deceduti.append(deceduti[i]-deceduti[i-1])
    new_guariti.append(guariti[i]-deceduti[i-1])
    
k_inf = k_factor(new, tot)
k_int = k_factor(new_intensiva, intensiva)
k_dec = k_factor(new_deceduti, deceduti)
k_gua = k_factor(new_guariti, guariti)
R0_inf = R0(k_inf, exp_time, inf_time)
R0_int = R0(k_int, exp_time, inf_time)
R0_dec = R0(k_dec, exp_time, inf_time)
R0_gua = R0(k_gua, exp_time, inf_time)

        
#-------------------------opzioni estetiche plot--------------------------------

fig = plt.figure(figsize=(15, 10))

ax1 = fig.add_subplot(221)
c1 = 'tab:blue'
ax1.grid()
ax1.xaxis.grid(True, which='minor', linestyle=':')
ax1.yaxis.grid(True, which='minor', linestyle=':')
ax1.set_title('Covid-19 Time Evolution - '+regione)
ax1.set_ylabel(r'$n^°$ of people')
ax1.tick_params(axis='y')
ax1.set_yscale('log')
plt.xticks(rotation=45)
maj_loc = mpl.ticker.MultipleLocator(base=5.0)
min_loc = mpl.ticker.MultipleLocator(base=1.0)
ax1.xaxis.set_major_locator(maj_loc)
ax1.xaxis.set_minor_locator(min_loc)

ax2 = fig.add_subplot(222)
ax2.grid()
c2 = 'tab:red'
ax2.xaxis.grid(True, which='minor', linestyle=':')
ax2.yaxis.grid(True, which='minor', linestyle=':')
ax2.set_title('Growth rate')
ax2.set_ylabel('k factor')
plt.xticks(rotation=45)
ax2.xaxis.set_major_locator(maj_loc)
ax2.xaxis.set_minor_locator(min_loc)
ax6 = ax2.twinx()
ax6.tick_params(axis='y', labelcolor='tab:green')
ax6.xaxis.set_major_locator(maj_loc)
ax6.xaxis.set_minor_locator(min_loc)

ax3 = fig.add_subplot(223)
ax3.grid()
ax3.xaxis.grid(True, which='minor', linestyle=':')
ax3.yaxis.grid(True, which='minor', linestyle=':')
ax3.set_title('SEIR model (from k factor)')
ax3.set_ylabel('$R_0$ index')
ax3.set_xlabel('Total infected')
ax3.tick_params(axis='y')
ax3.xaxis.set_major_locator(maj_loc)
ax3.xaxis.set_minor_locator(min_loc)
ax3.yaxis.set_major_locator(maj_loc)
ax3.yaxis.set_minor_locator(min_loc)
ax3.set_ylim(0,deceduti[0]+0.8)
#ax3.set_xscale('log')
plt.xticks(rotation=45)
ax7 = ax3.twinx()
ax7.tick_params(axis='y', labelcolor='tab:green')
ax7.xaxis.set_major_locator(maj_loc)
ax7.xaxis.set_minor_locator(min_loc)

'''seir vs sir model
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
ax4.set_ylim(0,30)
'''

ax5 = fig.add_subplot(224)
ax5.grid()
ax5.xaxis.grid(True, which='minor', linestyle=':')
ax5.yaxis.grid(True, which='minor', linestyle=':')
ax5.set_title('Deaths vs Healed')
ax5.set_ylabel(r'$n^°$ of people')
ax5.tick_params(axis='y')
ax5.set_yscale('log')
plt.xticks(rotation=45)
ax5.xaxis.set_major_locator(maj_loc)
ax5.xaxis.set_minor_locator(min_loc)


plt.subplots_adjust(left=.08, right=.95, bottom=0.08, top=0.95, wspace=.15, hspace=0.35)

#-------------------------plot dei dati-----------------------------------------


ax1.plot(giorni, tot, 'ko--',
         label=r'Protezione Civile data - Totale infetti (Total Infected)')
ax1.plot(giorni, new, 'o--', color='tab:orange',
         label=r'Protezione Civile data - Nuovi infetti (New Infected)')
ax1.plot(giorni, intensiva, '.--', color='tab:red',
         label=r'Protezione Civile data - Terapia intensiva (Intensive Care)')
ax1.plot(giorni, tamponi, 'o--', color='tab:grey',
         label=r'Protezione Civile data - Tamponi (Swabs)')
ax1.legend()

ax2.plot(giorni[20:], k_inf[20:], 's--', label='Infections')
#ax2.plot(giorni[20:], k_int[20:], 's--', label='Intensive')
ax2.plot(giorni[20:], k_dec[20:], 's--', label='Deceased')
ax6.plot(giorni[20:], k_gua[20:], 's--', color='tab:green', label='Healed')
ax2.legend(loc='best')
ax6.legend(loc='upper center')

ax3.plot(giorni[20:], R0_inf[20:],
             's--', label='Infections')
#ax3.errorbar(giorni[-1], R0_inf[-1], yerr=0.48, color='k', fmt = '.')
ax3.plot(giorni[20:], R0_dec[20:], 's--', label='Deceased')
ax7.plot(giorni[20:], R0_gua[20:], 's--', color='tab:green', label='Healed')

ax3.legend(loc='best')
ax7.legend(loc='upper center')

'''
ax4.plot(tot[7:-1], seir_model(E, I, R, N),
         's--', label=r'$R_0=\beta\epsilon/(\gamma +\mu)(\epsilon +\mu)$')
ax4.plot(tot[7:-1], sir_model(I, R, N), 's--', label=r'$R_0=\beta/\gamma$')
ax4.legend(loc='best')
'''

ax5.plot(giorni, guariti, 'go--', label=r'Protezione Civile data - Guariti dimessi (Recovered)')
ax5.plot(giorni, deceduti, 'ko--', label=r'Protezione Civile data - Deceduti (Deceased)')
ax5.legend()

fig.savefig('/Users/nigresson/Desktop/COVID19/k factor - '+regione)



