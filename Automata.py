import os

os.system('git clone https://github.com/pcm-dpc/COVID-19')
os.system('echo Studio andamento nazionale...')
os.system('python3 Nazione.py')
os.system('echo Studio andamento regionale...')
os.system('python3 Regione.py')

os.system('echo Aggiorno i grafici...')
os.system('git add k-factor\ Italy.png')
os.system('git add k\ factor\ -\ Lombardia.png')
os.system('git add README.md')
date=os.popen('date').read()
os.system('git commit -m "Aggiornato al giorno '+date+' "')
os.system('git push -u origin master')
os.system('echo Tutto è andato a buon fine...spero!')
os.system('open k-factor\ Italy.png k\ factor\ -\ Lombardia.png')
