;*******************************
;massper_run.pro
;;06.23.16
;Chantanelle Nava
;
; Wrapper code to run a simple mass/period grid of planets.
;
; Deleted structure field per_found, as it is redundant and
; did not update when combining datasets - NM 11/1/2016
;
; If /usenewdata is not set then a suffix is needed (i.e. suffix = 20161114)
;*******************************


pro massper_run, subdir, ident, nit, fapnit, star, targetlist, weighting, eccmax = eccmax, permin=permin, $
                 suffix=suffix, usenewdata=usenewdata, msu=msu

;Set default values 
if not keyword_set(eccmax) then eccmax = 0.99
if not keyword_set(permin) then permin = 0.01
  
;;;;;Convert subdirectory name to string. 
subdir = str(subdir)

;;;;;Define file paths.
;If massper_run is called with /msu, change filepaths to work on MSU computers
if keyword_set(msu) then begin 
   rootdir = '/home/student/'
   savepath = rootdir+'masspergrid/massper_saves/' + subdir + '/'
   temppath = rootdir
;if /msu is not called, use filepaths for local machine
endif else begin
   rootdir = '~/Dropbox/UM_Minerva/'
   savepath = rootdir+'code/sim_saves/massper_saves/' + subdir + '/'
   temppath='~/Desktop/'
endelse
;stop
tempstr = 'simtemp'
;logfile = 'sim_info_' + str(ident) + '.txt'
datasave = savepath + 'datasave_'+ str(ident)+ '.sav'
schedule = str(star)+'.' + str(targetlist) + '_' + str(weighting) + '.daynum.txt'

if keyword_set(suffix) then logfile = 'sim_info_' + str(ident) + suffix + '.txt' $
else logfile = 'sim_info_' + str(ident) + '.txt'


;;;;;Begin sim timing.
spawn, 'rm ' + savepath + logfile
spawn, 'mkdir ' + savepath
spawn, 'echo sim_begin : >> ' + savepath + logfile
spawn, 'date >> ' + savepath + logfile


;;;;;Define constant star and planet values.
mstar = 0.87    ;;solar masses for HD185144
inc = 90.    ;;degrees
ecc = 0.
arg = 0.
tper = 0. 


;;;;;Define other constants used.
G = 6.6726d-11      ;;SI units
secprday = 3600. * 24.
mearth = 5.976d24  ;;kg 
msun = 1.989d30    ;;kg

;;;;;Set simulation parameters.
;;;;;Create planet population ditributed evenly over log space in a specified 
;;;;;mass/period grid. 
nmass = 7
nper = 7
;fapnit = 1
shortfap = 0
maxfalse = 1.  ; percent
massrange = [3,30]   ; MEarth
perrange = [50,400]  ; days  
nplan = nmass * nper

mass = 10 ^ (findgen(nmass) * (alog10(massrange[1]) - $
       alog10(massrange[0])) / float(nmass-1) + alog10(massrange[0]))

per = 10 ^ (findgen(nper) * (alog10(perrange[1]) - $
      alog10(perrange[0])) / float(nper-1) + alog10(perrange[0]))


;;;;;Save appropriate info to sim_info log file. 
spawn, 'echo nit = ' + str(nit) + ' >> ' + savepath + logfile
spawn, 'echo fapnit = ' + str(fapnit) + ' >> ' + savepath + logfile
spawn, 'echo shortfapt = ' + str(shortfap) + ' >> ' + savepath + logfile
spawn, 'echo nmass x nper = ' + str(nmass) + ' x ' + str(nper) + ' >> ' + $
       savepath + logfile
spawn, 'echo massrange = ' + str(massrange[0]) + ',' + str(massrange[1]) + $
       ' >> ' + savepath + logfile
spawn, 'echo perrange = ' + str(perrange[0]) + ',' + str(perrange[1]) + $
       ' >> ' + savepath + logfile
spawn, 'echo eccmax = ' + str(eccmax) + ' >> ' + savepath + logfile
spawn, 'echo permin = ' + str(permin) + ' >> ' + savepath + logfile

;Get observation times
if keyword_set(usenewdata) then begin

   readcol, rootdir + 'masspergrid/' + schedule, obs_ts, format= 'd'
   ;readcol, rootdir+'code/HD185144.daynum.txt', obs_ts, format='d' ; for local

   datablock = fltarr(n_elements(obs_ts), nplan, nit) 
endif else restore, datasave


;;;;;Build structures to be filled with planet info.
planets = replicate({mass: 0d, per:0d, K:0d, $
                     fitK:fltarr(nit), fit_per:fltarr(nit), $
                     fit_ecc:fltarr(nit), found:intarr(nit), $
                     fapt:intarr(nit)}, nplan)


;;;;;Loop over planets. Load OR create data set and fit with RVLIN.
for m = 0,nmass-1 do begin
   for p = 0,nper-1 do begin

      ind = nmass * m + p

      ;;;;;Calculate RV semi amp for planet.
      K = (2. * !pi * G / (per[p] * secprday))^(1./3) * mass[m] * mearth * $
          sin(inc*!dtor) / (mstar * msun + mass[m] * mearth)^(2./3) / sqrt(1 - ecc^2)

      ;print, 'pl = ' + str(ind+1), ', mass = '+str(mass[m]), ', period = '+str(per[p]), ', K = '+str(K)

      ;Create arrays to hold the RVLIN results
      found = intarr(nit)
      fitK = fltarr(nit)
      fit_per = fltarr(nit)
      fit_ecc = fltarr(nit)
      fap = intarr(nit)
      falserate = intarr(nit)
      
      ;;;;;Create RV curve for planet.
      if keyword_set(usenewdata) then $
      rv = KeplerEq(mstar, mass[m], per[p], ecc, inc, arg, tper, time = obs_ts) 
      
      for i = 0,nit-1 do begin
         ;print,'Starting iteration '+str(i+1)+' of '+str(nit)  ; XXX
         
         ;;;;;Add noise to curve for sim data. 
         
         if keyword_set(usenewdata) then begin
            noise_maker, obs_ts, noise
            datablock[*,ind,i] = rv + noise
         endif
         
         ;;;;;Fit the data with rv_lin.
         err = obs_ts*0 + 1.  ;;;m/s from MINERVA 
         fit = rv_fit_mp(obs_ts, datablock[*,ind,i], err, chi=chi, eccmax=eccmax, permin = permin, /quiet)
         fitK[i] = fit[4]
         fit_per[i] = fit[0]
         fit_ecc[i] = fit[2]

         ;;;;;Check the fit with the fap test.
         if shortfap eq 0 then begin
            falserate = fapt(fapnit, obs_ts, datablock[*,ind,i], chi, err)
            if falserate lt maxfalse then found[i] = 1
         endif
         
         if shortfap ne 0 then begin
            findfap = short_fapt(fapnit, obs_ts, datablock[*,ind,i], chi, err, k, fit[4])
            found[i] = findfap[0]
            fap[i] = findfap[1]
            falserate[i] = findfap[2]
         endif

      endfor

      ;;;;;Fill in structure info.
      planets[ind].mass = mass[m]
      planets[ind].per = per[p]
      planets[ind].K = K
      planets[ind].fitK = fitK
      planets[ind].fit_per = fit_per
      planets[ind].fit_ecc = fit_ecc
      planets[ind].found = found
      planets[ind].fapt = fap

      ; log progress during a run
      spawn, 'echo planet '+str(ind+1)+' completed at:  >> ' + savepath + logfile
      spawn, 'date >> ' + savepath + logfile

      ;STOP

   endfor  ; period loop

endfor  ; mass loop

;;;;;Save final structure.

if keyword_set(suffix) then save, planets, filename = savepath + 'massper_run_' + str(ident) + suffix + '.sav' $
else save, planets, filename = savepath + 'massper_run_' + str(ident)  + '.sav'

if keyword_set(usenewdata) then save, obs_ts, datablock, filename = datasave

;;;;;Remove temp file.
spawn, 'rm ' + temppath + '*' + tempstr + '*' 

;;;;;End sim timing.
spawn, 'echo sim_end : >> ' + savepath + logfile
spawn, 'date >> ' + savepath + logfile



END
