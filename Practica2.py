from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

data = fits.getdata('Stokes_sunspot_HINODE.fits')
# Data es ahora un cubo de datos

# Primera dimension: parametros de Stokes (I, Q, U, V) (4)
# Segunda dimension: tiempo (399)
# Tercera dimension: puntos en los que esta dividida la rendija (401)
# Cuarta  dimension: diferentes longitudes de onda (96)

I = data[0,:,:,:]
Q = data[1,:,:,:]
U = data[2,:,:,:]
V = data[3,:,:,:]


# Para cada valor de tiempo existe un espectro 2D con 96 longitudes de onda y 401 puntos espaciales a lo largo de la rendija. Estos datos corresponden a una mancha solar.

# Segundo fichero que corresponde al atlas solar FTS que se usara para la calibracion

cal = fits.getdata('atlas_6301_6302.fits') # (2000, 2)


# Vamos a pintar una imagen del continuo del parametro I y seleccionar una zona del sol en calma. Claramente fuera de la mancha. Calcularemos el espectro medio de I en esa zona del sol en calma de los datos.




# Escogemos la longitud de onda en el continuo la cual es: I[:,:,0]
fig, ax = plt.subplots()
cmap = ax.imshow(I[:,:,0],origin='lower',cmap='binary_r')
cb = fig.colorbar(cmap)
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
plt.axvline(x=15, linestyle = '--',color='darkred')
plt.xlabel('x-px')
plt.ylabel('y-px')
plt.tight_layout()
plt.savefig('mapaI.png')

# En este caso mi y es el tiempo y mi x es el espacio de la rendija.

# Vamos a escoger una region del sol en calma perteneciente a (x,y) = (0:25,:)
Vprom = np.mean(V[:,15,:],axis=0)
Iprom = np.mean(I[:,15,:],axis=0) # calculamos el promedio
x = np.arange(Iprom.shape[0])
plt.figure()
for j in range(I.shape[0]):
    plt.plot(I[j,15,:],linewidth=0.3,color='gray')
plt.plot(x,Iprom,'darkred') # representamos el promedio
plt.xlabel(r'$\lambda$')
plt.ylabel(r'I (Cuentas)')
plt.tight_layout()
plt.savefig('Ipromd.png')
# YA ESTA NORMALIZADO :)



# VAMOS A LEER EL ATLAS SUMINISTRADO:

plt.figure()
plt.plot(cal[:,0],cal[:,1],'darkred')
plt.xlabel(r'$\lambda$ $(\AA)$')
plt.ylabel('Cuentas')


# Calculamos el minimo del atlas, si nos fijamos vemos que 
m1 = (cal[:,0]>6301.508)&(cal[:,0]<6301.515)
m2 = (cal[:,0]>6302.499)&(cal[:,0]<6302.505)


FeI_1 = cal[:,0][cal[:,1] == min(cal[:,1][m1])]
FeI_2 = cal[:,0][cal[:,1] == min(cal[:,1][m2])]

plt.plot(cal[:,0][cal[:,1] == FeI_1],cal[:,1][cal[:,1]==FeI_1],'x')
plt.plot(cal[:,0][cal[:,1] == FeI_2],cal[:,1][cal[:,1]==FeI_2],'x')


# Para calcular la precision del maximo con precision subpixel tenemos que ajustar a parabolas

# NORMALIZAMOS LOS DATOS DEL ATLAS

# -------
cal[:,1] = cal[:,1]/max(cal[:,1])

#plt.figure()
plt.plot(cal[:,0][m1],cal[:,1][m1])
x1 = np.linspace(cal[:,0][m1][0],cal[:,0][m1][-1],100)
fit1 = np.polyfit(cal[:,0][m1],cal[:,1][m1],2)
y1 = np.poly1d(fit1)
plt.plot(x1,y1(x1))


#plt.figure()
plt.plot(cal[:,0][m2],cal[:,1][m2])
x2 = np.linspace(cal[:,0][m2][0],cal[:,0][m2][-1],100)
fit2 = np.polyfit(cal[:,0][m2],cal[:,1][m2],2)
y2 = np.poly1d(fit2)
plt.plot(x2,y2(x2))

FeI1_teo = x1[y1(x1) == min(y1(x1))][0]
FeI2_teo = x2[y2(x2) == min(y2(x2))][0]

print('Los minimos se encuentran en:', FeI1_teo, FeI2_teo)

# Calculamos ahora los minimos de Hinode

plt.figure()
plt.plot(x,Iprom)

mh1 = (x>22.5)&(x<26)
fith1 = np.polyfit(x[mh1],Iprom[mh1],2)
yh1 = np.poly1d(fith1)
xh1 = np.linspace(x[mh1][0],x[mh1][-1],100)
plt.plot(xh1,yh1(xh1))



mh2 = (x>68.5)&(x<71.5)
fith2 = np.polyfit(x[mh2],Iprom[mh2],2)
yh2 = np.poly1d(fith2)
xh2 = np.linspace(x[mh2][0],x[mh2][-1],100)
plt.plot(xh2,yh2(xh2))



FeI1h_teo = xh1[yh1(xh1) == min(yh1(xh1))][0]
FeI2h_teo = xh2[yh2(xh2) == min(yh2(xh2))][0]

print('Los minimos se encuentran en:', FeI1h_teo, FeI2h_teo)

# Calculamos ahora la recta entre los dos minimos
linear = np.polyfit((FeI1h_teo, FeI2h_teo),(FeI1_teo,FeI2_teo),1)
lam = np.poly1d(linear)

# CAMBIAMOS LOS DATOS DE HINODE A LONGITUD DE ONDA Y NORMALIZAMOS
plt.figure()
lamI = lam(x)
Imean = np.mean(Iprom[lamI<6301.08])
Iprom = Iprom/Imean
Vprom = Vprom/max(Vprom)
# VAMOS A PINTAR EL PERFIL PROMEDIO FRENTE A DICHAS LONGITUDES DE ONDA.

plt.figure()
plt.plot(cal[:,0],cal[:,1],'black',linewidth=0.8,label='ALTAS Spectra')
plt.plot(lamI,Iprom,'darkred',label='Hinode Spectra')
plt.xlabel(r'$\lambda \AA$')
plt.ylabel('I')
plt.legend()
plt.tight_layout()
plt.savefig('twospectra.png')

# MEDIDA DEL CAMPO MAGNETICO.

# Vamos a calcular exactamente lo mismo en base a calibracion pero esta vez en la umbra. 


# Por otro lado, la separacion de los picos de V (lamdaB) se tiene que calcular evidentemente en la umbra porque es proporcional al campo y quermos medir el cmapo en la umbra.


# PARA UN SITIO EN LA UMBRA



fig,ax = plt.subplots()
cmap = plt.imshow(I[:,:,0],origin='lower',cmap='binary_r')
cb = fig.colorbar(cmap)
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
plt.plot(150,230,marker='o',fillstyle='none')
plt.xlabel('x-px')
plt.ylabel('y-px')
plt.tight_layout()
plt.savefig('puntoumbra.png')




V_u = V[150,230,:]
I_u = I[150,230,:]
Q_u = Q[150,230,:]
U_u = U[150,230,:]

x_u = np.arange(V_u.shape[0])

plt.figure()
for j in range(200,260):
    plt.plot(lam(x_u),V[j,150,:],linewidth=0.3,color='gray')
#    plt.plot(x_u,Vprom_u,'darkred') # representamos el promedio

plt.plot(lam(x_u),V_u,color='darkred')
plt.xlabel(r'$\lambda(\AA)$')
plt.ylabel('V (cuentas)')
plt.tight_layout()




cmap = plt.cm.coolwarm
rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, 4)))
fig,ax = plt.subplots()
ax.plot(lamI,I_u,label='I')
ax.plot(lamI,U_u,label='U')
ax.plot(lamI,Q_u,label='Q')
ax.plot(lamI,V_u,label='V')
plt.xlabel(r'$\lambda(\AA)$')
plt.ylabel('cuentas')
plt.legend()
plt.tight_layout()
plt.savefig('umbracomp.png')



plt.figure()
m1u = (lamI>6302.38)&(lamI<6302.44)
fit1u = np.polyfit(lamI[m1u],(I_u+V_u)[m1u],2)
y1u = np.poly1d(fit1u)
plt.plot(lamI,(I_u+V_u),'.',label='(I+V)')

x2u = np.linspace(6302.565,6302.655,100)
m2u = (lamI>6302.58)&(lamI<6302.63)
fit2u = np.polyfit(lamI[m2u],(I_u-V_u)[m2u],2)
y2u = np.poly1d(fit2u)
plt.plot(lamI,(I_u-V_u),'.',label='(I-V)')
plt.plot(x2u,y2u(x2u))

min1 = x1u[y1u(x1u) == min(y1u(x1u))][0]
min2 = x2u[y2u(x2u) == min(y2u(x2u))][0]

lambB = min2-min1

mfit = (lamI>6302.38)&(lamI<6302.44)
plt.plot(xfit,y1u(xfit))
s = np.std(xfit)
lambD = 2*np.sqrt(2*np.log(2))*s


h = 0.5*(max(yfit(xfit))+ min(yfit(xfit)))
xd = np.linspace(6302.35,6302.35+lambD,10)
yd = y = np.ones(10)*h

xb = np.linspace(min1,min2,100)
yb = np.linspace(min(y1u(x1u)),min(y2u(x2u)),100)


plt.plot(xd,yd,':k',label=r'$\Delta\lambda_D$')
plt.plot(xb,yb,'--k',label=r'$\Delta\lambda_B$')


plt.legend()
plt.xlabel(r'$\lambda(\AA)$')
plt.ylabel('Cuentas')
plt.tight_layout()
plt.savefig('anchurasumbra.png')


print(r'$\lambda_D$:' ,lambD)
print(r'$\lambda_B$:', lambB)

# Por tanto la intensidad de campo magnetico vendra dada por

C = 4.67e-13
lbd = 6302.51 # amstrong
geff = 2.5 # factor de lande efectivo

B_umbra = lambB/(2*C*geff*(lbd**2.))


# =========================================
# ==== CALCULAMOS LOS ANGULOS =============
# =========================================

inclinacion = (180/np.pi)*np.arctan2((Q_u**2. + U_u**2.)**0.5,V_u)[V_u == min(V_u)][0]
azimut = (180/np.pi)*(0.5*np.arctan2(U_u,Q_u))[U_u == min(U_u)][0]






fig, (ax1,ax2) = plt.subplots(1,2) 
cmap1 = ax1.imshow(abs(V[:,:,m])/I[:,:,m],origin='lower',cmap='binary_r')
cmap2 = ax2.imshow(P,origin='lower',cmap='binary_r') 
cmap2.axes.get_yaxis().set_visible(False) 
fig.colorbar(cmap1,ax=ax1,shrink=0.5)
fig.colorbar(cmap2,ax=ax2,shrink=0.5)

ax1.set_title('V')
ax2.set_title('P')

plt.savefig('polarizaciones.png')







# ===================================================

# APROXIMACION DE CAMPO DEBIL EN EL SOL EN CALMA

# ===================================================

fig,ax = plt.subplots()
cmap = plt.imshow(I[:,:,0],origin='lower',cmap='binary_r')
cb = fig.colorbar(cmap)
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
plt.plot(25,15,marker='o',fillstyle='none')
plt.savefig('puntocalma.png')


V_c = V[25,15,:]
I_c = I[25,15,:]
Q_c = Q[25,15,:]
U_c = U[25,15,:]

x_c = np.arange(V_c.shape[0])


plt.figure()
x1b = np.linspace(6301.37,6301.656,100)
m1b = (lamI>6301.485)&(lamI<6301.535)
fit1b = np.polyfit(lamI[m1b],(I_c+V_c)[m1b],2)
y1b = np.poly1d(fit1b)
plt.plot(lamI,(I_c+V_c),'.')
plt.plot(x1b,y1b(x1b))

x2b = np.linspace(6301.484,6301.54,100)
m2b = (lamI>6301.485)&(lamI<6301.535)
fit2b = np.polyfit(lamI[m2b],(I_c-V_c)[m2b],2)
y2b = np.poly1d(fit2b)
plt.plot(lamI,(I_c-V_c),'.')
plt.plot(x2b,y2b(x2b))

min1b = x1b[y1b(x1b) == min(y1b(x1b))][0]
min2b = x2b[y2b(x2b) == min(y2b(x2b))][0]

lambB_b = min2b-min1b






# CALCULAMOS LA ANCHURA DE DOPPLER AJUSTANDO A UNA GAUSSIANA

# Se hace con la gaussiana de y1b
x1d = np.linspace(6301.37,6301.656,100)
s = np.std(x1d)
lambD_b = 2*np.sqrt(2*np.log(2))*s


hb = 0.5*(max(y1b(x1d))+ min(y1b(x1d)))
xd = np.linspace(6301.41,6301.41+lambD_b,10)
yd = y = np.ones(10)*hb

xb = np.linspace(min1b,min2b,100)
yb = np.linspace(min(y1b(x1b)),min(y2b(x2b)),100)


plt.plot(xd,yd,':k',label=r'$\Delta\lambda_D$')
plt.plot(xb,yb,'--k',label=r'$\Delta\lambda_B$')

plt.savefig('anchurascalma.png')




fig,ax = plt.subplots()
ax.plot(lamI,I_c,label='I')
ax.plot(lamI,U_c,label='U')
ax.plot(lamI,Q_c,label='Q')
ax.plot(lamI,V_c,label='V')
plt.xlabel(r'$\lambda(\AA)$')
plt.ylabel('cuentas')
plt.legend()
plt.tight_layout()
plt.savefig('calmacomp.png')




lbd_b = 6301.52
geff_b = 1.667
mask = (lamI>6301.05)*(lamI<6301.98)                     
dI = (I_c[mask][2:]-I_c[mask][:-2])/(2*(lamI[1]-lamI[0]))


Blong = -np.sum(dI*V_c[mask][1:-1])/(np.sum(dI**2.)*(C*(lbd_b**2.)*geff_b))



mcont = (lamI>6301.75)&(lamI<6302.29)
sigma_c = np.std(V[25,25,:][mcont])

Berr = sigma_c/(C*(lbd_b**2.)*geff_b*(np.sum(dI**2.)**0.5))




