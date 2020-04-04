import os

os.system('rmdir COVID-19')
os.system('git clone https://github.com/pcm-dpc/COVID-19')
os.system('python3 Nazione.py')
os.system('python3 Regione.py')
os.system('python3 Provincia.py')
os.system('python3 Regioni.py')

os.system('git add growth\ k-factor\ comparison.png')
os.system('git add k\ factor\ -\ Lombardia.png')
os.system('git add k\ factor\ -\ Milano.png')
os.system('git add k\ factor.png')
os.system('git add README.md')
os.system('git commit -m "aggiornamento dati del giorno"')
os.system('git push -u origin master')
