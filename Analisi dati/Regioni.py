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

pop_lomb = 10000000
pop_camp = 6000000
pop_cal = 2000000
pop_ven = 5000000


#-------------------------lettura da file---------------------------------------

with open('/Users/nigresson/Desktop/COVID19/COVID-19/dati-json/dpc-covid19-ita-regioni.json') as regioni:
    dati_regioni = json.load(regioni)

#print(json.dumps(dati_regioni, indent=4, sort_keys=False))

#-------------------------analisi dati------------------------------------------

oggi = datetime.date(datetime.now())

n_giorni = (lambda delta: delta.days)(oggi - date(2020, 2, 24))

tot_lombardia = []
tot_campania = []
tot_calabria = []
tot_veneto = []
giorni = []
new_lombardia = []
new_campania = []
new_calabria = []
new_veneto = []
x_giorni = []
new_tamponi_lomb = []
new_tamponi_camp = []
new_tamponi_cal = []
new_tamponi_ven = []
counter=0

for d in range(0, len(dati_regioni)):
    if dati_regioni[d]['denominazione_regione']=='Lombardia':
        counter+=1
        tot_lombardia.append(dati_regioni[d]['totale_casi'])
        new_lombardia.append(dati_regioni[d]['nuovi_attualmente_positivi'])
        giorni.append(dati_regioni[d]['data'][:10])
        x_giorni.append(counter)
        new_tamponi_lomb.append(dati_regioni[d]['tamponi'])
    elif dati_regioni[d]['denominazione_regione']=='Campania':
        tot_campania.append(dati_regioni[d]['totale_casi'])
        new_campania.append(dati_regioni[d]['nuovi_attualmente_positivi'])
        new_tamponi_camp.append(dati_regioni[d]['tamponi'])
    elif dati_regioni[d]['denominazione_regione']=='Calabria':
        tot_calabria.append(dati_regioni[d]['totale_casi'])
        new_calabria.append(dati_regioni[d]['nuovi_attualmente_positivi'])
        new_tamponi_cal.append(dati_regioni[d]['tamponi'])
    elif dati_regioni[d]['denominazione_regione']=='Veneto':
        tot_veneto.append(dati_regioni[d]['totale_casi'])
        new_veneto.append(dati_regioni[d]['nuovi_attualmente_positivi'])
        new_tamponi_ven.append(dati_regioni[d]['tamponi'])



tamponi_lomb = []
tamponi_camp = []
tamponi_cal = []
tamponi_ven = []
tamponi_lomb.append(new_tamponi_lomb[0])
tamponi_camp.append(new_tamponi_camp[0])
tamponi_cal.append(new_tamponi_cal[0])
tamponi_ven.append(new_tamponi_ven[0])
for i in range(1, len(new_tamponi_lomb)):
    tamponi_lomb.append(new_tamponi_lomb[i]+new_tamponi_lomb[i-1])
    tamponi_camp.append(new_tamponi_camp[i]+new_tamponi_camp[i-1])
    tamponi_cal.append(new_tamponi_cal[i]+new_tamponi_cal[i-1])
    tamponi_ven.append(new_tamponi_ven[i]+new_tamponi_ven[i-1])

x=[]
y_lomb=[]
y_camp=[]
y_cal=[]
y_ven=[]
k_tot_l=[]
k_tot_c=[]
k_tot_cal=[]
k_tot_ven=[]
err_rel_lomb = []
err_rel_camp = []
err_rel_cal = []
err_rel_ven = []

#lombardia
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y_lomb.append(tot_lombardia[j])
    popt_l, pcov_l = curve_fit(gaussian, x, y_lomb)
    k_tot_l.append(popt_l[1])
    x.clear()
    y_lomb.clear()
    err_rel_lomb.append(1-tamponi_lomb[i]/pop_lomb)

#campania
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y_camp.append(tot_campania[j])
        y_camp = list(filter(lambda num: num != 0, y_camp))
    if len(y_camp)==4:
        popt_c, pcov_c = curve_fit(gaussian, x, y_camp)
        k_tot_c.append(popt_c[1])
    x.clear()
    y_camp.clear()
    err_rel_camp.append(1-tamponi_camp[i]/pop_camp)

#calabria
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y_cal.append(tot_calabria[j])
        y_cal = list(filter(lambda num: num != 0, y_cal))
    if len(y_cal)==4:
        popt_cal, pcov_cal = curve_fit(gaussian, x, y_cal)
        k_tot_cal.append(popt_cal[1])
    x.clear()
    y_cal.clear()
    err_rel_cal.append(1-tamponi_cal[i]/pop_cal)

#veneto
for i in range(int(inf_time/2), n_giorni - int(inf_time/2)):
    for j in range(i-int(inf_time/2), i+int(inf_time/2)):
        x.append(x_giorni[j])
        y_ven.append(tot_veneto[j])
        y_ven = list(filter(lambda num: num != 0, y_ven))
    if len(y_ven)==4:
        popt_ven, pcov_cal = curve_fit(gaussian, x, y_ven)
        k_tot_ven.append(popt_ven[1])
    x.clear()
    y_ven.clear()
    err_rel_ven.append(1-tamponi_ven[i]/pop_ven)
    

k_new_l=[]
k_new_c=[]
k_new_cal=[]
k_new_ven=[]

for i in range(0, len(tot_lombardia)):
    if tot_lombardia[i] != 0:
        k_new_l.append(new_lombardia[i]/tot_lombardia[i])
    else:
        k_new_l.append(0)
    if tot_campania[i] != 0:
        k_new_c.append(new_campania[i]/tot_campania[i])
    else:
        k_new_c.append(0)
    if tot_calabria[i] != 0:
        k_new_cal.append(new_calabria[i]/tot_calabria[i])
    else:
        k_new_cal.append(0)
    if tot_veneto[i] != 0:
        k_new_ven.append(new_veneto[i]/tot_veneto[i])
    else:
        k_new_ven.append(0)


popt_lin_lomb, pcov_lin_lomb = curve_fit(linear, tot_lombardia[-25:], k_new_l[-25:])
print(-popt_lin_lomb[1]/popt_lin_lomb[0])
        
#-------------------------opzioni estetiche plot--------------------------------

fig = plt.figure(figsize=(18, 9))
c2 = 'tab:red'
c1 = 'tab:blue'

ax1 = fig.add_subplot(221)
ax1.grid()
ax1.xaxis.grid(True, which='minor', linestyle=':')
ax1.yaxis.grid(True, which='minor', linestyle=':')
ax1.set_title('Lombardia')
ax1.set_ylabel('k factor')
ax1.set_xlabel('Total infected')
ax1.set_xscale('log')
ax1.set_ylim(-.1,1)

ax2 = fig.add_subplot(222)
ax2.grid()
ax2.xaxis.grid(True, which='minor', linestyle=':')
ax2.yaxis.grid(True, which='minor', linestyle=':')
ax2.set_title('Campania')
ax2.set_ylabel('k factor')
ax2.set_xlabel('Total infected')
ax2.set_xscale('log')
ax2.set_ylim(-.1,1)

ax3 = fig.add_subplot(223)
ax3.grid()
ax3.xaxis.grid(True, which='minor', linestyle=':')
ax3.yaxis.grid(True, which='minor', linestyle=':')
ax3.set_title('Calabria')
ax3.set_ylabel('k factor')
ax3.set_xlabel('Total infected')
ax3.set_xscale('log')
ax3.set_ylim(-.1,1)

ax4 = fig.add_subplot(224)
ax4.grid()
ax4.xaxis.grid(True, which='minor', linestyle=':')
ax4.yaxis.grid(True, which='minor', linestyle=':')
ax4.set_title('Veneto')
ax4.set_ylabel('k factor')
ax4.set_xlabel('Total infected')
ax4.set_xscale('log')
ax4.set_ylim(-.1,1)

plt.subplots_adjust(left=.05, right=.95, bottom=0.08, top=0.95, wspace=.15, hspace=0.35)

#-------------------------plot dei dati-----------------------------------------

#lombardia
ax1.axvline(x=tot_lombardia[giorni.index('2020-02-25')], color='tab:grey')
ax1.axvline(x=tot_lombardia[giorni.index('2020-03-01')], color='tab:grey')
ax1.axvline(x=tot_lombardia[giorni.index('2020-03-04')], color='tab:grey')
ax1.axvline(x=tot_lombardia[giorni.index('2020-03-08')], color='tab:grey')
ax1.axvline(x=tot_lombardia[giorni.index('2020-03-09')], color='tab:grey')
ax1.plot(tot_lombardia, k_new_l,
         's--', color=c1, label='based on new daily infected')
ax1.plot(tot_lombardia[int(inf_time/2)+1: n_giorni - int(inf_time/2)], k_tot_l[1:],
         's--', color=c2, label='based on total infected')
ax1.legend()

#campania
ax2.axvline(x=tot_campania[giorni.index('2020-02-25')], color='tab:grey')
ax2.axvline(x=tot_campania[giorni.index('2020-03-01')], color='tab:grey')
ax2.axvline(x=tot_campania[giorni.index('2020-03-04')], color='tab:grey')
ax2.axvline(x=tot_campania[giorni.index('2020-03-08')], color='tab:grey')
ax2.axvline(x=tot_campania[giorni.index('2020-03-09')], color='tab:grey')
ax2.plot(tot_campania, k_new_c,
         's--', color=c1, label='based on new daily infected')
ax2.plot(tot_campania[int(inf_time/2)+4: n_giorni - int(inf_time/2)], k_tot_c[1:],
         's--', color=c2, label='based on total infected')
ax2.legend()

#calabria
ax3.axvline(x=tot_calabria[giorni.index('2020-02-25')], color='tab:grey')
ax3.axvline(x=tot_calabria[giorni.index('2020-03-01')], color='tab:grey')
ax3.axvline(x=tot_calabria[giorni.index('2020-03-04')], color='tab:grey')
ax3.axvline(x=tot_calabria[giorni.index('2020-03-08')], color='tab:grey')
ax3.axvline(x=tot_calabria[giorni.index('2020-03-09')], color='tab:grey')
ax3.plot(tot_calabria, k_new_cal,
         's--', color=c1, label='based on new daily infected')
ax3.plot(tot_calabria[int(inf_time/2)+5: n_giorni - int(inf_time/2)], k_tot_cal[1:],
         's--', color=c2, label='based on total infected')
ax3.legend()

#veneto
ax4.axvline(x=tot_veneto[giorni.index('2020-02-25')], color='tab:grey')
ax4.axvline(x=tot_veneto[giorni.index('2020-03-01')], color='tab:grey')
ax4.axvline(x=tot_veneto[giorni.index('2020-03-04')], color='tab:grey')
ax4.axvline(x=tot_veneto[giorni.index('2020-03-08')], color='tab:grey')
ax4.axvline(x=tot_veneto[giorni.index('2020-03-09')], color='tab:grey')
ax4.plot(tot_veneto, k_new_ven,
         's--', color=c1, label='based on new daily infected')
ax4.plot(tot_veneto[int(inf_time/2)+1: n_giorni - int(inf_time/2)], k_tot_ven[1:],
         's--', color=c2, label='based on total infected')
ax4.legend()

fig.savefig('/Users/nigresson/Desktop/COVID19/k factor')



