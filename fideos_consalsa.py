import pandas as pd
import numpy as np
import warnings
# import json
from astropy import units as u
from astropy.coordinates import SkyCoord , Angle, get_sun, EarthLocation, AltAz, get_moon
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
from astropy.time import Time, TimeDelta
from datetime import datetime
# from urllib.request import urlopen


print("Downloading the latest TOI list...")
TOI_url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
CTOI_url = 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=pipe'
TOI_df = pd.read_csv(TOI_url, delimiter='|', index_col=1)
CTOI_df = pd.read_csv(CTOI_url, delimiter='|', index_col=1)
print("Download finished!")

#  Filtering CTOIs promoted to TOIs
CTOI_df = CTOI_df[np.array(np.isnan(CTOI_df['Promoted to TOI']))]

# magMask = np.zeros(len(CTOI_df))
# i = 0
# for toi in CTOI_df.iloc:
#     print('looping :(')
#     url = f"https://exofop.ipac.caltech.edu/tess/target.php?id={toi['TIC ID']}&json"
#     response = urlopen(url)
#     data_json = json.loads(response.read())
#     for mag in data_json['magnitudes']:
#         if mag['band'] == 'V':
#             if float(mag['value']) < 9.5:
#                 magMask[i] = 1
#     i += 1

print("Number of TOIs:", len(TOI_df))
print("Number of CTOIs:", len(CTOI_df))

# Creamos un dataframe con PCs y APCs con TESS mag <= 10.
bool_mask = TOI_df['TESS Mag'] <= 9.5
TOI_dfMag = TOI_df[bool_mask]
bool_mask = CTOI_df['TESS Mag'] <= 9.5
CTOI_df = CTOI_df[bool_mask]
bool_mask = TOI_dfMag['TFOPWG Disposition'] == 'PC'
TOI_df = TOI_dfMag[bool_mask]
# bool_mask = TOI_dfMag['TFOPWG Disposition'] == 'APC'
# TOI_dfAPC = TOI_dfMag[bool_mask]
bool_mask = CTOI_df['User Disposition'] == 'PC'
CTOI_df = CTOI_df[bool_mask]
# TOI_df = pd.concat([TOI_dfPC, TOI_dfAPC])
TOI_df['Date TOI Alerted (UTC)'] = pd.to_datetime(TOI_df['Date TOI Alerted (UTC)'])
TOI_df = TOI_df.sort_values(by='Date TOI Alerted (UTC)', ascending=False)
singlePlanets = np.flatnonzero(np.core.defchararray.find(np.array(list(map(str, TOI_df.index))), '.01') != -1)
singlePlanets = TOI_df.index[singlePlanets]
TOI_df = TOI_df.loc[singlePlanets]
singlePlanets = np.flatnonzero(np.core.defchararray.find(np.array(list(map(str, CTOI_df.index))), '.01') != -1)
singlePlanets = CTOI_df.index[singlePlanets]
CTOI_df = CTOI_df.loc[singlePlanets]
with open('dontObserve.txt', 'r') as f:
    dontObserve = [line.rstrip() for line in f]
print('The following targets will not be considered for observations:')
for i in dontObserve:
    posTOI = np.flatnonzero(np.core.defchararray.find(np.array(list(map(str, TOI_df['TIC ID']))), i.strip('TIC')) != -1)
    posCTOI = np.flatnonzero(np.core.defchararray.find(np.array(list(map(str, CTOI_df['TIC ID']))), i.strip('TIC')) != -1)
    for j in posTOI:
        drop_index = TOI_df.index[j]
        TOI_df = TOI_df.drop(drop_index)
    for j in posCTOI:
        drop_index = CTOI_df.index[j]
        CTOI_df = CTOI_df.drop(drop_index)
    print(f'\t{i}')
    # Reads the target names that won't be observed.
print("Number of TOIs PCs under TESS mag 9.5:", len(TOI_df))
print("Number of CTOIs PCs under TESS mag 9.5:", len(CTOI_df),'\n')
# Comienza el ciclo
days_range = input('Choose the number of days to simulate:\n\t')
days_range = int(days_range)
date_input = input('Choose a CLT monday or wednesday for when to start the simulation with the format\n\tyyyy-mm-dd:\n\t')
date = date_input.split('-')
try:
    weekday = datetime(int(date[0]), int(date[1]), int(date[2])).weekday()
except:
    raise Exception('Wrong format for date. The format is:\n\tyyyy-mm-dd')
if weekday != 0 and weekday != 2:
    raise Exception('Please choose a valid date (Monday or Wednesday)')
if weekday == 0:
    print(f'Observing plan starts on\n\tMonday {date[0]}/{date[1]}/{date[2]}.')
else :
    print(f'Observing plan starts on\n\tWednesday {date[0]}, {date[1]}, {date[2]}.')
laSilla = EarthLocation.from_geodetic(-70.7375, -29.2575, [2347])
with open('nightStartTime.txt', 'w') as f:
    pass
    # Creates the text file or cleans it
contador = 0
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for day in range(days_range):
        if weekday == 0:  # Monday
            if not(day % 7 == 0 or (day + 5) % 7 == 0):
                continue  # In this version this statement skips days that are not M,W.
        elif weekday == 2:  # Wednesday
            if not(day % 7 == 0 or (day + 2) % 7 == 0):
                continue  # In this version this statement skips days that are not W,M.
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
        steps = int(np.round(obsTimeDelta.sec, 0) / 60)
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
        samplingIndexes = np.linspace(30, steps-30, int(steps/60), dtype=int)
        sampledOutput   = np.zeros([len(TOI_df) + len(CTOI_df), int(steps/60)])
        y = 0
        for target in TOI_df.iloc:
            targetPos = SkyCoord(ra=target['RA'], dec=target['Dec'], distance=target['Stellar Distance (pc)'] ,
                                 unit=(u.hourangle, u.deg, u.pc), obstime=obsDateTime)
            with erfa_astrom.set(ErfaAstromInterpolator(int(steps/60) * u.hour)): 
                temp = targetPos.transform_to(aaFrame)
                if np.array([temp.alt.value>=50]).sum() >= 60 and target['SG2'] <= 3:
                    # Si el target esta sobre 50 grados por 1 hora o mas, se muestrea cada 1 hora.
                    x = 0
                    for i in samplingIndexes:
                        # Si no esta sobre 50 grados por 1 hora minimo, no se rellena.
                        sampledOutput[y][x] = temp[i].alt.value
                        x += 1
            y += 1
        for target in CTOI_df.iloc:
            targetPos = SkyCoord(ra=target['RA'], dec=target['Dec'], distance=target['Stellar Distance (pc)'] ,
                                 unit=(u.hourangle, u.deg, u.pc), obstime=obsDateTime)
            with erfa_astrom.set(ErfaAstromInterpolator(int(steps/60) * u.hour)): 
                temp = targetPos.transform_to(aaFrame)
                if np.array([temp.alt.value>=50]).sum() >= 60:
                    # Si el target esta sobre 50 grados por 1 hora o mas, se muestrea cada 1 hora.
                    x = 0
                    for i in samplingIndexes:
                        # Si no esta sobre 50 grados por 1 hora minimo, no se rellena.
                        sampledOutput[y][x] = temp[i].alt.value
                        x += 1
            y += 1
        print(f'Simulating nights: {np.round((day+1) * 100 /days_range, 2)}%', end="\r", flush=True)
        pd.DataFrame(sampledOutput).to_csv(f"data{contador}.csv", index=False, header=False)
        contador += 1
with open('specObservations.txt', 'w') as f:
    pass
with open('specObservations.txt', 'a') as f:
    for target in TOI_df.iloc:
        if target['Period (days)'] > days_range / 2 + 1 or target['Period (days)'] == 0:
            prio = '0'
        else:
            prio = str(np.exp(-target['Spectroscopy Observations']) / target['TESS Mag'])
        f.write(str(target['TIC ID']) +'|'+ prio +'|'+ str(target['Period (days)'])+'\n')
    for target in CTOI_df.iloc:
        if target['Period (days)'] > days_range / 2 + 1 or target['Period (days)'] == 0:
            prio = '0'
        else:
            prio = str(np.exp(-(target['TESS Mag'] - 9)))
        f.write(str(target['TIC ID']) +'|'+ prio +'|'+ str(target['Period (days)'])+'\n')
print('Simulation done! Now starting the optimization process...')