from numpy import *

def scalecol( A=None, mults=None ):

# Scales each column of A by the respective multiplier in mults
# Equivalent to B = A * diag( mults )
# except the full matrix doesn't need to be constructed
#
# mults can be either a row or a column vector
#
# usage: B = scalecol( A, mults )
# adapted from scalecol.m 

   s = shape( A )
   B = zeros( s ) + 1j*zeros( s )
   if size( s ) == 1:
      B = mults*A
   else:
      Ncols = s[1]
      for i in range(Ncols):
          B[:,i] = mults[i]*A[:,i]

   return B
