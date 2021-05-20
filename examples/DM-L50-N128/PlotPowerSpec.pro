
fout  = "powerspec_z0.eps"
RunDir = "./output/"
Tag    = "DM-L50-N128"
inpsec_filename  = "./output/inputspec_snapshot.txt"

BoxSize = 50.0

Num = 2         ;  snapshot numbers at z=0

if num ge 1000 then begin
   exts='0000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
endif else begin
   exts='000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
endelse



mydevice=!d.name
set_plot,'PS'
!p.font=0
device,/times,/italic,font_index=20
device,/times
device,xsize=20.0,ysize=15.0
!x.margin=[10,3]
!p.thick=3.0
!p.ticklen=0.03

device, filename = fout, /encapsulated, /color

v1 = [255,   0,   0]
v2 = [  0, 185,   0]
v3 = [  0,   0, 255]
tvlct,v1,v2,v3,1

plot, [1], [1], /xlog, /ylog , $
      xtitle = "!20k!7 [ !20h!7 Mpc!U-1!N ]", ytitle = "!9D!7!U2!N(!20k!7)", thick=3.0, $
      xrange = 2*!PI/Boxsize*[0.7, 2000], xstyle = 1, yrange=[0.1, 8.0e3], ystyle=1,  xthick=3.5, ythick=3.5, charsize=1.2

oplot, 2.0*!PI/BoxSize * [1, 1], [1.0e-10,1.0e10], linestyle=2


Kall = [0]
Dall = [0]

; read in the three (folded) power spectrum measurements, and stick
; them together to reach higher k

FoldFac         = 16L
MinModeCount    = 8             ; rebin to get at least this number of modes per bin
TargetBinNummer = 50            ; aim for this number of bins 


fname = RunDir +  "/powerspecs/powerspec_" + exts + ".txt"
openr, 1, fname

for piece = 0, 2 do begin

   Time = 0.0D 
   Bins = 0L
   BoxSize = 0.0D
   Ngrid = 0.0D
   Dplus = 0.0D
   readf, 1, Time
   
   print, "piece=", piece, "                    TIME=", time
   
   readf, 1, Bins
   readf, 1, BoxSize
   readf, 1, Ngrid
   readf, 1, Dgrowth 
   da= dblarr(5, bins)
   readf, 1, da

   K         = da(0,*)
   Delta2    = da(1,*)
   ModePow   = da(2,*)
   ModeCount = da(3,*)
   Shot      = da(4,*)

   
   if (piece eq 0) then begin  ;; let's plot linear theory spectrum, scaled to present time

      openr,2, inpsec_filename 
      z = 0.0D
      dplus = 0.0D 
      readf,2,z,dplus
      dat = dblarr(4, 514)
      readf,2,dat
      close,2
      
      k_lin = dat(0,*)
      D2_lin = dat(1,*)
      D2_Lin /= Dgrowth^2

      oplot, k_lin, D2_lin

   endif

      
   
   oplot, K, Shot, color=13, thick=3.0
                                      
   Delta2 -= Shot  ; do a shot-noise subtraction
      
   SumPower =  Delta2 * ModeCount
   
   kmin = 2*!PI/BoxSize * double(FoldFac)^piece
   kmax = 2*!PI/BoxSize * Ngrid/2.0/4  * double(FoldFac)^piece

   if piece gt 0 then begin
      kmin = 2*!PI/BoxSize * Ngrid/2.0/4 * double(FoldFac)^(piece-1.0)
   endif

   print,  "Kmin=", Kmin, "  Kmax = ", Kmax


; we will do a band averaging of the finely binned points, 
; subject to two conditions:
; We want enough modes per bin in order to reduce the variance in a bin,
; and simultaneously, for large k, we don't want the bins become too narrow.
;
; The first condition is set by "MinModeCount",
; the second by "TargetBinNummer", which is used to compute a minimum 
; logarithmic bin-size.

      
   MinDlogK = (alog10(max(K)) - alog10(min(K)))/TargetbinNummer
   
   istart=0
   ind=[istart]
   k_list = [0]
   Delta2_list = [0]
   Delta2lin_list = [0]
   count_list = [0]
   repeat begin
      count = total(modecount(ind))
      deltak =  (alog10(max(K(ind))) - alog10(min(K(ind))))
      
      if (deltak ge mindlogk) or ((max(k(ind)) lt 2.0 * kmin) and (piece eq 0)) then begin
         
         kk = exp(total(alog(k(ind))*ModeCount(ind), /double)/total(ModeCount(ind),/double))
         d2 = total(Delta2(ind) * ModeCount(ind) , /double)/total(ModeCount(ind), /double)
         Delta2lin = exp(interpol(alog(D2_lin), alog(k_lin), alog(k(ind))))
         d2lin = total(Delta2lin * ModeCount(ind) , /double)/total(ModeCount(ind), /double)

         if d2 gt 0 then begin
            
            k_list = [k_list, kk]
            Delta2_list = [Delta2_list, d2]
            Delta2lin_list = [Delta2lin_list, d2lin]
            count_list = [count_list, total(ModeCount(ind))]
         endif
         
         istart = istart + 1
         ind = [istart]
      endif else begin
         istart = istart + 1
         ind = [ind, istart]
      endelse
   endrep until istart ge Bins

   
   K_list = k_list(1:*)
   Delta2_list = delta2_list(1:*)
   Count_list = count_list(1:*)
   Delta2lin_list = Delta2lin_list(1:*)

                                ; select the measurement of the range [kmin,kmax] and adopt it for our total measurement      
   ind = where((K_List ge kmin) and (K_list le kmax))
   
   Kall = [Kall, K_list(ind)]  
   Dall = [Dall, Delta2_list(ind)]

endfor

close,1
   
; get rid off the trailing zeros
Kall = Kall[1:*]
Dall = Dall[1:*]


;restrict the measurement to be at most a factor 1/2 below the shot noise

Sall = exp(interpol(alog(shot), alog(k), alog(Kall)))

ind = where(Dall gt Sall/2)
   
; plot it
oplot, Kall(ind), Dall(ind), color=1, thick=4

   
xyouts, 0.22, 0.76, /normal, tag, color = 3
xyouts, 0.22, 0.85, /normal, "!7power spectrum, !20z!7 = 0", charsize=1.2

device,/close
set_plot,"X"


end
