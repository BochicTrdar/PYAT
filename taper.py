from numpy    import *
from scalecol import *

def taper( A=None, k=None, Nk=None, kleft=None, kright=None ):

# windowing to smooth out any discontinuities at the end of the spectrum
# adapted from taper in fieldsco.m 

    if ( kleft > k[0]  ):
       Nwinleft = 2 * round( ( kleft - k[0] ) / ( k[-1] - k[0] ) * Nk ) + 1  # odd number covering 10% of the spectrum
       winleft  = hanning( Nwinleft )
       Nwinleft = ( Nwinleft - 1 ) / 2
       winleft  = winleft[ 0 : Nwinleft+1 ] # taking just the left half of the symmetric Hanning window
    else:
       Nwinleft = 0
       winleft  = []

    if ( kright < k[-1]  ):
       Nwinright = 2 * round( ( k[-1] - kright ) / ( k[-1] - k[0] ) * Nk ) + 1 # odd number covering 10% of the spectrum
       winright  = hanning( Nwinright )
       Nwinright = ( Nwinright - 1 ) / 2
       winright  = winright[ end - Nwinright + 1 :] # taking just the right half of the symmetric Hanning window
    else:
       Nwinright = 0
       winright  = []

    if ( Nk < Nwinleft + Nwinright ):
       print( 'phase speed limits for windowing are invalid' )

    window = concatenate((winleft, ones(Nk - Nwinleft - Nwinright), winright), axis=None)
    
    G = scalecol( A, window )

    return G
