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

def gaussian(x, a, k):
    return a*np.exp(k*x)

def linear(x, m, q):
    return m*x+q

inf_time = 5

regione='Lombardia'

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
tamponi = []
counter=0

for d in range(0, len(dati_regioni)-1):
    if dati_regioni[d]['denominazione_regione']==regione:
        counter+=1
        tot.append(dati_regioni[d]['totale_casi'])
        new.append(dati_regioni[d]['nuovi_positivi'])
        giorni.append(dati_regioni[d]['data'][:10])
        x_giorni.append(counter)
        tamponi.append(dati_regioni[d]['tamponi'])

x=[]
y=[]
k_tot=[]

for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y.append(tot[j])
    popt, pcov = curve_fit(gaussian, x, y)
    k_tot.append(popt[1])
    x.clear()
    y.clear()

k_new=[]

for i in range(0, len(tot)):
    if tot[i] != 0:
        k_new.append(new[i]/tot[i])
    else:
        k_new.append(0)

popt_1, pcov_1 = curve_fit(linear, tot[5:], k_new[5:])
pred1 = -popt_1[1]/popt_1[0]
try:
    popt_2, pcov_2 = curve_fit(linear, tot[2:-3], k_tot[1:])
    pred2 = -popt_2[1]/popt_2[0]
except:
    popt_2, pcov_2 = curve_fit(linear, tot[3:-3], k_tot[1:])
    pred2 = -popt_2[1]/popt_2[0]

        
#-------------------------opzioni estetiche plot--------------------------------

fig = plt.figure(figsize=(10, 9))

ax1 = fig.add_subplot(211)
c1 = 'tab:blue'
ax1.grid()
ax1.xaxis.grid(True, which='minor', linestyle=':')
ax1.yaxis.grid(True, which='minor', linestyle=':')
ax1.set_title('Covid-19 Time Evolution - '+regione)
ax1.set_ylabel('Total infected')
ax1.tick_params(axis='y')
#ax1.set_yscale('log')
plt.xticks(rotation=45)
maj_loc = mpl.ticker.MultipleLocator(base=3.0)
min_loc = mpl.ticker.MultipleLocator(base=1.0)
ax1.xaxis.set_major_locator(maj_loc)
ax1.xaxis.set_minor_locator(min_loc)

ax = fig.add_subplot(212)
ax.grid()
c2 = 'tab:red'
ax.xaxis.grid(True, which='minor', linestyle=':')
ax.yaxis.grid(True, which='minor', linestyle=':')
ax.set_title('Growth k factor')
ax.set_ylabel('k factor', color=c2)
ax.set_xlabel('Total infected')
ax.tick_params(axis='y', labelcolor=c2)
ax.set_xscale('log')
plt.xticks(rotation=45)
ax.set_ylim(0,1)

plt.subplots_adjust(left=.08, right=.95, bottom=0.08, top=0.95, wspace=.15, hspace=0.35)

#-------------------------plot dei dati-----------------------------------------


ax1.plot(giorni, tot, 'ko--', label=r'Protezione Civile data')
ax1.legend()

ax.axvline(x=tot[giorni.index('2020-02-25')], color='tab:grey')
ax.axvline(x=tot[giorni.index('2020-03-01')], color='tab:grey')
ax.axvline(x=tot[giorni.index('2020-03-04')], color='tab:grey')
ax.axvline(x=tot[giorni.index('2020-03-08')], color='tab:grey')
ax.axvline(x=tot[giorni.index('2020-03-09')], color='tab:grey')
ax.plot(tot, k_new,
         's--', color=c1, label='based on new daily infected')
ax.plot(tot[int(inf_time/2)+1: n_giorni - int(inf_time/2)], k_tot[1:],
         's--', color=c2, label='based on total infected')

ax.plot(tot[5:], linear(np.asarray(tot[5:]), popt_1[0], popt_1[1]),
         '--', color='tab:grey')
ax.plot(tot[3:-2], linear(
    np.asarray(tot[3:-2]), popt_2[0], popt_2[1]),'--', color='tab:grey')
ax.plot(pred1, 0, 'ko')
ax.plot(pred2, 0, 'ko')
ax.annotate('predicted\n %i\n infected' %pred1, xy=(pred1, 0),
             xytext=(pred1-10000, 0.25), color='tab:blue',
             arrowprops=dict(arrowstyle='->'))
ax.annotate('predicted\n %i\n infected' %pred2, xy=(pred2, 0),
             xytext=(pred2-15000, -0.2), color='tab:red',
             arrowprops=dict(arrowstyle='->'))
ax.legend(loc='best')

fig.savefig('/Users/nigresson/Desktop/COVID19/k factor - '+regione)



