######################################
# This code was made by Francisco Oyarzun in 2022
# Any questions send them to ftoyarzun at uc.cl
######################################


import numpy as np
from gurobipy import *


print("Reading data from files...")

specObs = np.genfromtxt('specObservations.txt', delimiter='|')
nightStart = np.genfromtxt('nightStartTime.txt', delimiter='|').astype(int)

#Accessing data from files
names = specObs[:,0].astype(int)
Priority = specObs[:,1]
Period = specObs[:,2]*24 // 1
horas_totales = (nightStart[-1,0]+1)*24

#Massage period for useful information
subPeriod = [int(min(5, x // 10 + 1)) for x in Period]

Priority[np.isnan(Period)] = 0
Priority[Period == 0] = 0
Period[Period == 0] = 1e9
Period[np.isnan(Period)] = 1e9
Priority[2*Period >= horas_totales] = 0

print("Reading data from files completed")

print("Building dictionaries...")

J_ = []
D_ = list(range(np.shape(nightStart)[0]))

a = {}
for d in D_:
    data = np.genfromtxt(f'data{d}.csv', delimiter=',')
    if d == 0:
        I_ = range(data.shape[0])
    J_.append(range(data.shape[1]))
    for j in J_[d]:
        for i in I_:
            a[i,j,d] = data[i,j]
            
H_ = range(horas_totales)
DU = list(nightStart[:,0])
DU.sort()

alt = {}
cnt = 0


for h in H_:
    d = h//24
    hd = h%24
    for i in I_:
        if (d in DU) and hd in J_[DU.index(d)]:
            alt[i,h] = a[i,hd,DU.index(d)]
        else:
            alt[i,h] = 0

V_ = []
for i in I_:
    temp = []
    if 2*Period[i] < horas_totales:
        for h in range(horas_totales - 2*int(Period[i])):
            temp.append(range(h, h+2*int(Period[i])-2))
    else:
        temp.append(range(horas_totales))
    V_.append(temp)
    
l = []
for i in I_:
    for h in H_:
        l.append((i,h))
        
vv = [(i,v) for i in I_ for v in range(len(V_[i]))]

print("Building dictionaries completed")

print("Building model...")

model = Model("Programa de observacion")

x = model.addVars(l,vtype = GRB.BINARY, name = "Observar")
y = model.addVars(vv,vtype = GRB.BINARY, name = "Observar")
z = model.addVars(I_,vtype = GRB.BINARY, name = "Activar target")


model.addConstrs((quicksum(x[i,h] for i in I_) <= 1 for h in H_), name = "Observar un solo target a la vez")
model.addConstrs(((x[i,h] * (alt[i,h] - 50)) >= 0 for i in I_ for h in H_), name = "respetar airmass")
model.addConstrs((quicksum(x[i,h] for h in H_) <= 10000 * z[i] for i in I_), name = "Activar target (1/2)")
model.addConstrs((z[i] <= quicksum(x[i,h] for h in H_) for i in I_), name = "Activar target (2/2)")
model.addConstrs((quicksum(x[i,h] for h in H_) <= 12 * z[i] for i in I_), name = "No observar mas de doce veces un mismo target")
model.addConstrs((quicksum(y[i,v] for v in range(len(V_[i]))) >= z[i] for i in I_), name = "Si se observa un target, tiene que haber alguna ventana activa")
model.addConstrs((quicksum(x[i,h] for h in V_[i][v]) >= 10*y[i,v] for i in I_ for v in range(len(V_[i]))), name = "Observar al menos cinco veces durante la ventana")
model.addConstrs((quicksum(x[i,k] for k in range(h,h+subPeriod[i])) <= 1 for i in I_ for h in range(horas_totales - subPeriod[i])), name = "Separacion entre observaciones proporcional al periodo")


print("Building model completed")

print("Begin optimization...")

obj = quicksum(quicksum(x[i,h] * Priority[i] * alt[i,h] for i in I_) for h in H_)
model.setObjective(obj, GRB.MAXIMIZE)
model.setParam('TimeLimit', 120)
                 

# Actualizar modelo
model.update()       

# Optimizar
model.optimize()

print("Optimization completed")

print("Generating text files ...")

f = open('Planer_por_fecha.txt', 'w')
f.write(f'El dia 0 corresponde al {nightStart[0,1]}-{nightStart[0,2]}-{nightStart[0,3]}')
f.close()
f = open('Planer_por_target.txt', 'w')
f.close()

previous_d = -1
cnt = 0
with open('Planer_por_fecha.txt', 'a') as file:
    for h in H_:
        for i in I_:
            if x[i,h].x>0:
                
                day = h//24
                hd = h%24
                if day != previous_d:
                    file.write('\n\n')
                    file.write(f'***** Dia {day}. La noche comienza a las {int(nightStart[cnt,-2]):02d}:{int(nightStart[cnt,-1]):02d} UTC*****\n')
                    cnt += 1
                file.write(f"Observar el target TIC{int(names[i])} el dia {day} a las {hd}:00 UTC\n")
                file.write(f'https://simbad.cds.unistra.fr/simbad/sim-id?Ident=TIC{int(names[i])}\n')
                previous_d = day
              
           
       
previous_i = -1
with open('Planer_por_target.txt', 'a') as file:
    for i in I_:
        for h in H_:
            if x[i,h].x>0:
                day = h//24
                hd = h%24
                if i != previous_i:
                    file.write('\n\n')
                    file.write(f'*****El periodo del target TIC{int(names[i])} es de {Period[i]} horas *****\n')
                file.write(f"Observar el target TIC{int(names[i])} el dia {day} a las {hd}:00 UTC\n")
                file.write(f'https://simbad.cds.unistra.fr/simbad/sim-id?Ident=TIC{int(names[i])}\n')
                previous_i = i

print("Generating text files completed")