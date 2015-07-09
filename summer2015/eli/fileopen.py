import numpy as np
import matplotlib.pylab as plt

#LADO'S WAY
lado1=np.loadtxt('cmass_north_dr10.DD')

#OUR INTERPRETATION
lado2=np.loadtxt('DD_lado_tot.dat')

lado11=np.flipud(lado1)
lado111=lado1+lado11

lado22=np.fliplr(lado2)
lado222=lado2+lado22

lado333=lado222-lado111



plt.figure()
rangeval=300
extent=[-rangeval,rangeval,-rangeval,rangeval]

plt.subplot(1,3,1)
a=plt.imshow(lado1,extent=extent)
plt.title('Lado')
plt.ylabel(r'$r_\perp (h^{-1}$Mpc)')
plt.xlabel(r'$r_\parallel (h^{-1}$Mpc)')

plt.colorbar()
plt.subplot(1,3,2)
b=plt.imshow(lado2,extent=extent)
plt.title('Our Interpretation')
plt.ylabel(r'$r_\perp (h^{-1}$Mpc)')
plt.xlabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.colorbar()

plt.subplot(1,3,3)
b=plt.imshow(lado333,extent=extent)
plt.title('Difference')
plt.ylabel(r'$r_\perp (h^{-1}$Mpc)')
plt.xlabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.colorbar()


plt.show()
