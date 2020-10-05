
da = dblarr(8, 61)

openr, 1, "output/energy.txt"
readf, 1, da
close, 1

ti =  da(0,*)   ; time
th =  da(1,*)   ; thermal energy
po =  da(2,*)   ; potential energy
ke =  da(3,*)   ; kinertic energy


tot = th + ke + po

window,xsize=1000,ysize=900

!p.multi=[0,1,2]

plot, ti, tot, charsize=2.0, yrange=[-2,1.5], ystyle=1

oplot,ti, ke, color=255*256L^2 + 100*256L + 100
oplot,ti, po, color=255*256L
oplot,ti, th, color=255

plot,ti,(tot-tot(0))/abs(tot(0)),charsize=2.0


end

