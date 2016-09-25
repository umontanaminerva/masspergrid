;**************************
;fapt.pro
;Function to run the false alarm probability test for given data set
;and number of iterations.
;
; INPUTS
;    nit     number of shuffle iterations 
;    times   time stamp array for RV data
;    data    RV array for target star
;    chi     chi squared statistic of best fit to actual (unscrambled) data
;    err     error array for RV data
;
; OUTPUT
;    falserate   returns percent of false positives (scramble chi > actual chi)
;
;Chantanelle Nava
;May 26, 2014
;**************************

function fapt, nit, times, data, chi, err

falsepos=intarr(nit)
for i=0, nit-1 do begin
   check=shuffle(data)
   testfit=rv_fit_mp(times, check, err, chi=testchi, /quiet)
   if testchi gt chi then falsepos[i]=0 else falsepos[i]=1   
endfor

falserate=float(total(falsepos)) / nit * 100.

return, falserate

END
