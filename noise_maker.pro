;******************************
;noise_maker.pro
;A procedure returning a set of noise values given the two parameters from
;an exponential squared gaussian fit and a set of observation times for
;which noise is desired. 
;Chantanelle Nava
;September 24, 2014 (at the CfA, pretty cool)


pro noise_maker, times, noise

;;a is the amplitude of the spread-> higher a gives larger vertical spread

;;g is the level of variance within a period-> higher g gives more
;;wiggles within period

;;P is the "period" of the quasi-periodic kernel-> pretty self-explanatory

;;l is the weight of correlation between points of a set distance->
;;higher l gives smoother, higher-correlated points


;;;eq form ---> k(r)=a^2 * exp(-g^2 * (sin(!pi * abs(r) / P)^2) -
;;;r^2/2./l^2)


;;;;;Define necessary arrays
noise=fltarr(n_elements(times))
kernel_vals=fltarr(n_elements(times), n_elements(times))
means=fltarr(n_elements(times))

;;;;;Read in amps and l_vals from 'noise_vars.py'.
readcol, 'noise_vars.txt', amps, l_vals, /silent

;;;;;Determine the covariance between a and l.
mat=transpose([[l_vals], [amps]])
cov=correlate(mat, /covariance)

;;;;;Determine means and standard deviations.
lmean=mean(l_vals)
ampmean=mean(amps)
lstdev=stdev(l_vals)
ampstdev=stdev(amps)

status=1 ;;to begin while loop

while status ne 0 do begin
;;;;;Return a and l values using the covariance matrix and stats.
   la=mrandomn(seed, cov, 1)
   l=la[0]*lstdev+lmean
   a=la[1]*ampstdev+ampmean

;;;;;Build kernel matrix.
   for i=0, n_elements(times)-1 do begin
      for j=0, n_elements(times)-1 do begin
         r=times[j]-times[i]
         ;kernel_vals[i,j]=a^2. * exp(-(g^2.) * (sin(!pi * abs(r) / P))^2. - r^2. / 2. / l^2.)
         kernel_vals[i,j]=a^2. * exp(-(r^2.) / 2. / l^2.)
      endfor
   endfor

   noise=mrandomn(seed, kernel_vals, status=status)

endwhile


END
