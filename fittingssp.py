# Curve fitting to SSP:
from numpy import *
from matplotlib.pyplot import *

data = loadtxt('elba.ssp')

z = squeeze( data[:,0] ); Dmax = max( z )
c = squeeze( data[:,1] )

theorder = 6
pars = polyfit(z, c, theorder)
p    = poly1d(pars)
cp   = p(z)

figure(1)
plot(c,z)
plot(cp,z,'o')
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
title('Polynomial fitting to SSP')
ylim(Dmax,0)
grid(True)

show()

print('done.')
