
 a1 = readfits('outputsph.fits')
loadct,3
set_plot,'ps'
device,file='temp.eps',bits=24,/encapsulat
MAP_SET, /MOLLWEIDE, 0, 0, /ISOTROPIC,TITLE='North-south velocities (Rotation 1)', /NOBORDER
  contour,a1,/overplot,nlevels=20;,min_value=-10,max_value=10
device,/close
set_plot,'x'
read,pause
end
