path = '../Output/'

;goto, plouf
file = path + 'cls.dat'
readcol, file, l, cl_hikap, cl_hihi, cl_kapkap

xrange = [50., 1000.]
xstyle = 1
;plot_oo, l, cl_kapkap, $
;         xrange = xrange, xstyle = xstyle 


yrange = [1.e-4, 1.e-2]
plot_oo, l, cl_hihi, $
         xrange = xrange, xstyle = xstyle, $
         yrange = yrange, ystyle = ystyle
stop
plot_oo, l, cl_hikap, $
         xrange = xrange, xstyle = xstyle
stop
;plouf:

xrange = [5.e-3,1.]
xstyle = 1
hub = 0.67
file = path + 'P_HI_z=0.1.dat'
readcol, file, k, PHI
plot_oo, k * hub , PHI / hub^3, $
         xrange = xrange, xstyle = xstyle 


end 
