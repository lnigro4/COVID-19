# covid-19: a quick data analysis

In Nazione.py I analyze the protezione civile data regarding the total of infected people in Italy and the growth rate (or growth factor k). The blue line "based on new infected" is obtained by the ratio #daily_infected/#tot_infected_uptodate. The red line "based on total infected" is instead obtained by fitting in a range of 5 days across the whole period (i.e. a sort of mobile fit).

In Regione.py I analyze the protezione civile data regarding the total of infected people in Lombardia (it can be changed) and the growth rate (or growth factor k). The blue line "based on new infected" is obtained by the ratio #daily_infected/#tot_infected_uptodate. The red line "based on total infected" is instead obtained by fitting in a range of 5 days across the whole period (i.e. a sort of mobile fit).

In Provincia.py I analyze the protezione civile data regarding the total of infected people in Milan (it can be changed) and the growth rate (or growth factor k). The red line "based on total infected" is instead obtained by fitting in a range of 5 days across the whole period (i.e. a sort of mobile fit).

In Regioni.py I analyze 4 regions at the time.


In each case we also study the R0 reproduction index. It's obtained by the results for the k factor. A more analytical approach (i.e. with the ode's relative to the SEIR and SIR models) is put as a comment in Nazione.py and Regione.py. Couldn't do it for provinces as well because of the lack of data on the corresponding dataset. 