
da = dblarr(16, 61)

openr, 1, "output_old/energy.txt"
readf, 1, da
close, 1

ti =  da(0,*)   ; time
th =  da(1,*)   ; thermal energy
po =  da(2,*)   ; potential energy
ke =  da(3,*)   ; kinertic energy


tot = th + ke + po

window,xsize=1000,ysize=900

!p.multi=[0,1,2]

plot, ti,(po-po(0))/abs(po(0)),charsize=2.0, yrange=[-1,1] * max([ abs((po-po(0))/abs(po(0))), abs((ke-ke(0))/abs(ke(0)))])

oplot,ti,(ke-ke(0))/abs(ke(0)),linestyle=2, color=255
oplot,ti,(po-po(0))/abs(po(0)),linestyle=3, color=255*256L

plot,ti,(tot-tot(0))/abs(tot(0)),charsize=2.0


end

