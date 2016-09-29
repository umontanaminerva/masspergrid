;*******************************
;massper_run.pro
;;06.23.16
;Chantanelle Nava
;
; Wrapper code to run a simple mass/period grid of planets.
;
;*******************************


pro massper_run, subdir, ident, nit, fapnit

old_sched = 0


;;;;;Convert subdirectory name to string. 
subdir = str(subdir)

;;;;;Define file paths.
;savepath = '~/Dropbox/UM_Minerva/code/sim_saves/massper_saves/' + subdir + '/'
;temppath='~/Desktop/'
tempstr = 'simtemp'
logfile = 'sim_info_' + str(ident) + '.txt'
savepath = '/home/student/masspergrid/massper_saves/' + subdir + '/'
temppath = '/home/student/'

;;;;;Begin sim timing.
spawn, 'rm ' + savepath + logfile
spawn, 'mkdir ' + savepath
spawn, 'echo sim_begin : >> ' + savepath + logfile
spawn, 'date >> ' + savepath + logfile


;;;;;Define constant star and planet values.
mstar = 0.87    ;;solar masses
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

if old_sched eq 1 then begin    ;;;XXX
;;;;;Sam's scheduled obs times.  
obsdir = '~/Dropbox/UM_Minerva/Scheduler/ThreeYearRuns/'
obsfile = 'exp_cad_6_3.txt'
;obsfile = 'three_obs.txt'
offset = -43. / 24.   ;;;;diff between Samson and Chani's 
feed1 =  "awk '($1=="
feed2 = '"HD185144"'
feed3 = "){print $6}' " + obsdir + obsfile + " > " + temppath + tempstr + $
        '_0.txt'
spawn, feed1 + feed2 + feed3
readcol, temppath + tempstr + '_0.txt', obs_ts
obs_ts = obs_ts + offset

endif


;;Brute force for now
;readcol, '~/Dropbox/UM_Minerva/eta_Earth/schedsample/HD185144.jd.txt', obs_ts
;readcol, '~/Dropbox/UM_Minerva/Scheduler/ThreeYearRuns/HD185144.exp_cad_6_3.txt', obs_ts
readcol, '/home/student/masspergrid/HD185144.exp_cad_6_3jd.txt', obs_ts
;obs_ts = obs_ts - min(obs_ts)  ; don't use this line! NM 9/27/16


;;;;;Build structures to be filled with planet info.
planets = replicate({mass: 0d, per:0d, K:0d, per_found:0d, $
                     fitK:fltarr(nit), fit_per:fltarr(nit), $
                     fit_ecc:fltarr(nit), found:intarr(nit), $
                     fapt:intarr(nit)}, nplan)

datablock = fltarr(n_elements(obs_ts), nplan, nit)


;;;;;Loop over planets.
for m = 0,nmass-1 do begin
   for p = 0,nper-1 do begin

      ind = nmass * m + p

      ;;;;;Calculate RV semi amp for planet.
      K = (2. * !pi * G / (per[p] * secprday))^(1./3) * mass[m] * mearth * $
          sin(inc*!dtor) / (mstar * msun + mass[m] * mearth)^(2./3) / sqrt(1 - ecc^2)

      print, 'pl = ' + str(ind+1), ', mass = '+str(mass[m]), ', period = '+str(per[p]), ', K = '+str(K)

      ;;;;;Create RV curve for planet.
      rv = KeplerEq(mstar, mass[m], per[p], ecc, inc, arg, tper, time = obs_ts)

      found = intarr(nit)
      fitK = fltarr(nit)
      fit_per = fltarr(nit)
      fit_ecc = fltarr(nit)
      fap = intarr(nit)
      falserate = intarr(nit)

      for i = 0,nit-1 do begin
         print,'Starting iteration '+str(i+1)+' of '+str(nit)  ; XXX
         
         ;;;;;Add noise to curve for sim data. 
                                ;noise_maker, obs_ts, noise
         noise = obs_ts*0. + randomn(seed,n_elements(obs_ts))
         data = rv + noise
         datablock[*,ind, i] = data

         ;;;;;Fit the data with rv_lin.
         err = 1.  ;;;m/s from MINERVA 
         fit = rv_fit_mp(obs_ts, data, err, chi=chi, /quiet)
         fitK[i] = fit[4]
         fit_per[i] = fit[0]
         fit_ecc[i] = fit[2]

         ;;;;;Check the fit with the fap test.
         if shortfap eq 0 then begin
            falserate = fapt(fapnit, obs_ts, data, chi, err)
            if falserate lt maxfalse then found[i] = 1
         endif
         
         if shortfap ne 0 then begin
            findfap = short_fapt(fapnit, obs_ts, data, chi, err, k, fit[4])
            found[i] = findfap[0]
            fap[i] = findfap[1]
            falserate[i] = findfap[2]
         endif

      endfor

      per_found = float(total(found)) / nit * 1d2


      ;;;;;Fill in structure info.
      planets[ind].mass = mass[m]
      planets[ind].per = per[p]
      planets[ind].K = K
      planets[ind].per_found = per_found
      planets[ind].fitK = fitK
      planets[ind].fit_per = fit_per
      planets[ind].fit_ecc = fit_ecc
      planets[ind].found = found
      planets[ind].fapt = fap

      ; log progress during a run
      spawn, 'echo planet '+str(ind+1)+' completed at:  >> ' + savepath + logfile
      spawn, 'date >> ' + savepath + logfile
      stop 
   endfor
endfor

;;;;;Save final structure.
save, planets, filename = savepath + 'massper_run_' + str(ident)+ '.sav'
save, obs_ts, datablock, filename = savepath + 'datasave_'+ str(ident)+ '.sav' 

;;;;;Remove temp file.
spawn, 'rm ' + temppath + '*' + tempstr + '*' 

;;;;;End sim timing.
spawn, 'echo sim_end : >> ' + savepath + logfile
spawn, 'date >> ' + savepath + logfile



END
