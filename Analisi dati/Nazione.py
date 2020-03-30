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

FARE STESSA COSA ANCHE PER ANDAMENTO REGIONE(LOMBARDIA) E PROVINCIA(MILANO)

POI RIPETERE ANCHE PER LOGISTICHE/DIST CANONICA CON TAU AL POSTO DI K
'''

#-------------------------precompilatore----------------------------------------

def gaussian(x, a, k):
    return a*np.exp(k*x)

def linear(x, m, q):
    return m*x+q

inf_time = 5

#-------------------------lettura da file---------------------------------------

with open('/Users/nigresson/Desktop/COVID19/COVID-19/dati-json/dpc-covid19-ita-andamento-nazionale.json') as nazione:
    dati_nazione = json.load(nazione)
    

#print(json.dumps(dati_nazione, indent=4, sort_keys=False))

#-------------------------analisi dati------------------------------------------

oggi = datetime.date(datetime.now())

n_giorni = (lambda delta: delta.days)(oggi - date(2020, 2, 24))

tot_casi = []
giorni = []
new_casi = []
x_giorni = []
for d in range(0, n_giorni):
    tot_casi.append(dati_nazione[d]['totale_casi'])
    new_casi.append(dati_nazione[d]['nuovi_attualmente_positivi'])
    giorni.append(dati_nazione[d]['data'][:10])
    x_giorni.append(d)

x=[]
y=[]
k_tot=[]
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y.append(tot_casi[j])
    popt, pcov = curve_fit(gaussian, x, y)
    k_tot.append(popt[1])
    x.clear()
    y.clear()

k_new=[]
for i in range(0, len(tot_casi)):
    k_new.append(new_casi[i]/tot_casi[i])

popt_1, pcov_1 = curve_fit(linear, tot_casi[5:], k_new[5:])
pred1 = -popt_1[1]/popt_1[0]

popt_2, pcov_2 = curve_fit(linear, tot_casi[3:-2], k_tot[1:])
pred2 = -popt_2[1]/popt_2[0]
    
#-------------------------opzioni estetiche plot--------------------------------

fig = plt.figure(figsize=(10, 9))

ax1 = fig.add_subplot(211)
c1 = 'tab:blue'
ax1.grid()
ax1.xaxis.grid(True, which='minor', linestyle=':')
ax1.yaxis.grid(True, which='minor', linestyle=':')
ax1.set_title('Covid-19 Time Evolution - Italy')
ax1.set_ylabel('Total infected')
ax1.tick_params(axis='y')
#ax1.set_yscale('log')
plt.xticks(rotation=45)
maj_loc = mpl.ticker.MultipleLocator(base=3.0)
min_loc = mpl.ticker.MultipleLocator(base=1.0)
ax1.xaxis.set_major_locator(maj_loc)
ax1.xaxis.set_minor_locator(min_loc)

ax2 = fig.add_subplot(212)
ax2.grid()
c2 = 'tab:red'
ax2.xaxis.grid(True, which='minor', linestyle=':')
ax2.yaxis.grid(True, which='minor', linestyle=':')
ax2.set_title('Growth k factor')
ax2.set_ylabel('k factor', color=c2)
ax2.set_xlabel('Total infected')
ax2.tick_params(axis='y', labelcolor=c2)
ax2.set_xscale('log')
plt.xticks(rotation=45)
ax2.set_ylim(0,1)
#ax2.xaxis.set_major_locator(maj_loc)
#ax2.xaxis.set_minor_locator(min_loc)

plt.subplots_adjust(left=0.08, bottom=0.08, top=0.95, hspace=0.35)

#-------------------------plot dei dati-----------------------------------------

ax1.plot(giorni, tot_casi, 'ko--', label=r'Protezione Civile data')
ax1.legend()

ax2.plot(tot_casi[int(inf_time/2)+1: n_giorni - int(inf_time/2)], k_tot[1:],
         's--', color=c2, label='based on total infected')
ax2.plot(tot_casi, k_new,
         's--', color=c1, label='based on new daily infected')
ax2.plot(tot_casi[5:], linear(np.asarray(tot_casi[5:]), popt_1[0], popt_1[1]),
         '--', color='tab:grey')
ax2.plot(tot_casi[3:-2], linear(
    np.asarray(tot_casi[3:-2]), popt_2[0], popt_2[1]),'--', color='tab:grey')
ax2.plot(pred1, 0, 'ko')
ax2.plot(pred2, 0, 'ko')
ax2.annotate('predicted\n %i\n infected' %pred1, xy=(pred1, 0),
             xytext=(pred1 +40000, 0.1), color='tab:blue',
             arrowprops=dict(arrowstyle='->'))
ax2.annotate('predicted\n %i\n infected' %pred2, xy=(pred2, 0),
             xytext=(pred2-50000, -0.15), color='tab:red',
             arrowprops=dict(arrowstyle='->'))
ax2.legend(loc='best')

fig.savefig('/Users/nigresson/Desktop/COVID19/growth k-factor comparison')




