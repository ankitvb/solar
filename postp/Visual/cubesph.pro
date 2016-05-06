PRO CUBESPH, filename
a = readfits(filename)
g = size(a)
nx = g[1]

f = dblarr(nx*4,nx*3)
f[0:nx-1,nx:2*nx-1] = a[*,*,4-1]
f[nx:2*nx-1,nx:2*nx-1] = a[*,*,1-1]
f[2*nx:3*nx-1,nx:2*nx-1] = a[*,*,2-1]
f[3*nx:4*nx-1,nx:2*nx-1] = a[*,*,3-1]

f[nx:2*nx-1,0:nx-1] = a[*,*,6-1]
f[nx:2*nx-1,2*nx:3*nx-1] = a[*,*,5-1]

tvim,f,/sc,noframe='no' 
oplot,[0, 4*nx-1],[nx,nx]
oplot,[0, 4*nx-1],[2*nx,2*nx]

oplot,[nx,nx],[0,3*nx-1]
oplot,[2*nx,2*nx],[0,3*nx-1]

oplot,[0,0],[nx,2*nx-1]

oplot,[3*nx,3*nx],[nx,2*nx-1]

oplot,[4*nx-1,4*nx-1],[nx,2*nx-1]

oplot,[nx,2*nx-1],[0,0]

oplot,[nx,2*nx-1],[3*nx-1,3*nx-1],thick=3

read,pause
end
