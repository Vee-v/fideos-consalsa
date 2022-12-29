import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord , Angle, get_sun, EarthLocation, AltAz, get_moon
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
from astropy.time import Time, TimeDelta
from datetime import datetime


print("Downloading the latest TOI list...")
url="https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
TOI_df=pd.read_csv(url, delimiter='|', index_col=1)
print("Download finished!")

print("Number of TOIs:", len(TOI_df))

# Creamos un dataframe con PCs y APCs con TESS mag <= 10.
bool_mask = TOI_df['TESS Mag'] <= 10
TOI_dfMag = TOI_df[bool_mask]
bool_mask = TOI_dfMag['TFOPWG Disposition'] == 'PC'
TOI_dfPC = TOI_dfMag[bool_mask]
bool_mask = TOI_dfMag['TFOPWG Disposition'] == 'APC'
TOI_dfAPC = TOI_dfMag[bool_mask]
TOI_df = pd.concat([TOI_dfPC, TOI_dfAPC])
TOI_df['Date TOI Alerted (UTC)'] = pd.to_datetime(TOI_df['Date TOI Alerted (UTC)'])
TOI_df = TOI_df.sort_values(by='Date TOI Alerted (UTC)', ascending=False)
print("Number of TOIs PCs and APCs under TESS mag 10:", len(TOI_df))
# Comienza el ciclo
days_range = input('Choose the number of days to simulate:\n')
days_range = int(days_range)
date_input = input('**ATTENTION, THIS PROGRAM IS MEANT TO BE RUN FOR THE NIGHT OF 12/30/2022**\nChoose \
a CLT friday for when to start the simulation with the format\n yyyy-mm-dd:\n')
date = date_input.split('-')
laSilla = EarthLocation.from_geodetic(-70.7375, -29.2575,[2347])
with open('nightStartTime.txt', 'w') as f:
    pass
    # Creates the text file or cleans it
with open('dontObserve.txt', 'r') as f:
    dontObserve = [line.rstrip() for line in f]
print('The following targets will not be considered for observations')
for i in dontObserve:
    print(i)
    # Reads the target names that won't be observed.
contador = 0
for day in range(days_range):
    if not((day % 7 == 0 or (day + 4) % 7 == 0 or (day + 5) % 7 == 0) or (day + 2) % 7 == 0):
        continue  # In this version this statement skips days that are not M,W,F,S.
    times = Time(datetime(int(date[0]), int(date[1]), int(date[2]), 23, 59, 0)) + day*u.day + np.linspace(-6, 9, 901)* u.hour
    # Ahora se busca el momento donde el sol cruza los -18° de altitud desde la posicion en La Silla
    
    sunPos = get_sun(times)
    aaFrame = AltAz(obstime=times, location=laSilla)
    aa_sunPos = sunPos.transform_to(aaFrame) 
    begin = np.argmax(aa_sunPos.alt.value < -12)
    end = np.argmax(np.flip(aa_sunPos.alt).value < -18)
    twiBegin = times[begin]
    twiEnd = np.flip(times)[end]
    # Ahora que tenemos la hora de inicio y de termino del twilight, tenemos el tiempo para calcular efemerides de luna y targets.
    obsTimeDelta = twiEnd + .5*u.hour - (twiBegin - .5*u.hour)
    obsTimeDelta = TimeDelta((obsTimeDelta.value * 24 // 1) / 24)
    steps = int(np.round(obsTimeDelta.sec,0) / 60)
    obsDateTime = twiBegin - .5*u.hour + obsTimeDelta * np.linspace(0, 1, steps+1)
    with open('nightStartTime.txt', 'a') as f:
        f.write(str(day)+'|'+
                str(obsDateTime[0].datetime.year)+'|'+
                str(obsDateTime[0].datetime.month)+'|'+
                str(obsDateTime[0].datetime.day)+'|'+
                str(obsDateTime[0].datetime.hour)+'|'+
                str(obsDateTime[0].datetime.minute)+'\n'
               )
    # Marco AltAz de La Silla
    aaFrame = AltAz(obstime=obsDateTime, location=laSilla)
    # Efemerides Luna
    moonPos = get_moon(obsDateTime, laSilla)
    aa_moonPos = moonPos.transform_to(aaFrame)  # TODO: output the Moon distance in the sky between targets and the moon.
    # Indices de muestreo para la altura del target cada una hora
    samplingIndexes=np.linspace(30, steps-30, int(steps/60), dtype=int)
    sampledOutput = np.zeros([len(TOI_df), int(steps/60)])
    y = 0
    for target in TOI_df.iloc:
        targetPos = SkyCoord(ra=target['RA'], dec=target['Dec'], distance=target['Stellar Distance (pc)'] ,
                             unit=(u.hourangle, u.deg, u.pc), obstime=obsDateTime)
        with erfa_astrom.set(ErfaAstromInterpolator(int(steps/60) * u.hour)): 
            temp = targetPos.transform_to(aaFrame)
            if np.array([temp.alt.value>=50]).sum() >= 60 and target['SG2'] <= 3 and not(target['TIC ID'] in dontObserve):
                # Si el target esta sobre 50 grados por 1 hora o mas, se muestrea cada 1 hora.
                x = 0
                for i in samplingIndexes:
                    # Si no esta sobre 50 grados por 1 hora minimo, no se rellena.
                    sampledOutput[y][x] = temp[i].alt.value
                    x += 1
        y += 1
    print(f'Simulating nights: {np.round(day+1 * 100 /91, 2)}%', end="\r", flush=True)
    pd.DataFrame(sampledOutput).to_csv(f"data{contador}.csv", index=False, header=False)
    contador += 1
with open('specObservations.txt', 'w') as f:
    pass
with open('specObservations.txt', 'a') as f:
    for target in TOI_df.iloc:
        f.write(str(target['TIC ID'])+'|'+str(target['Spectroscopy Observations'])+\
                '|'+str(target['SG2'])+'|'+str(target['Period (days)'])+'\n')
print('Simulation done! Now starting the optimization process...')