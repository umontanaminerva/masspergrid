;**************************                                                
;KeplerEq.pro
;February 2, 2013
;A program returning the radial velocity values of a theoretical
;planet as a function of star mass, planet mass, orbital parameters, and time.
;Written by: Chantanelle Nava
;**************************

;KeplerEq is a function that returns radial velocity values for a
;hypothetical planet/host star pair.


;Input:

;output=KeplerEq(mstar, mplan, perday, ecc, inc, arg, tper, time=time)

;mstar: mass of the host star in solar masses

;mplan: mass of the planet in earth masses

;perday: period of the planet's orbit in days

;ecc: eccentricity of the orbit

;inc: inclination of the orbit in degrees

;arg: argument of periastron of the orbit in degrees

;tper: time of periastron in days

;times: optional parameter designating a 1-D  array of  observation
;times in days for which radial velocities are to be solved. If
;undesignated, default times simulate daily observations over two orbital periods.    


;Output:

;output= array with same number of elements as 'times' containing
;corresponding radial velocity values.

; Modification History
;
; Added option to display syntax. NM 9-11-13

;**************************


function KeplerEq, mstarin, mplanin, perday, ecc, incin, argin, tper, time=time

;mstarin in solar masses
;mplanin in earth masses
;perday in days
;ecc
;inc in degrees
;arg in degrees
;tper in days
;times in days 

if(n_params() eq 0) then begin
   print,'syntax: result = KeplerEq(mstarin, mplanin, perday, ecc, incin, 
          $ argin, tper, time=time)'
   return,-1
endif

;Error checking
if(ecc ge 1. or ecc lt 0) then begin
   print, 'eccentricity out of range'
   stop 
endif

;;;;;List constants.                                                
G = 6.6730d-11         ;;kg*m^3/s^2
msun = 1.9891d30       ;;kg
mearth = 5.9722d24     ;;kg

;;;;;Convert planet and orbital parameters.                                    
mstar = mstarin * msun
mplan = mplanin * mearth
arg = argin * !dtor
inc = incin * !dtor
peryear = perday / 365.242
persec = perday * 24. * 3600.


;;;;;Calculate semi-major axis.                                        
aplan=((persec^2. * G * mstar) / (4. * !pi^2))^(1./3.)



;;;;;Solve Kepler's Equation to determine eccentric anomaly as
;;;;;a function of time.                                  
if (not keyword_set(time)) then begin
   numper = 2.
   time = findgen(perday * numper) + 1. + tper 
endif 


eccanom = fltarr(n_elements(time))

phase = (time-tper)/perday
manom = 2d*!dpi*(phase-floor(phase))
tol = 1d-5

for x=0L,n_elements(time)-1 do begin
   high = 2 * !dpi
   low = 0
   diff = 2 * tol        ;to start while loop
   E = !dpi             ;initial guess

   while abs(diff) gt tol do begin
      diff = E - manom[x] - ecc * sin(E)
      if abs(diff) gt tol then begin
         if diff lt 0 then low = E else high = E
         E = (low + high) / 2.
      endif
   endwhile

   eccanom[x] = E

endfor


;;;;;Calculate true anomaly from eccentric anomaly.
trueanom = acos((cos(eccanom) - ecc)/(1 - ecc * cos(eccanom)))
ind = where(eccanom gt !pi, n_ind)
if n_ind gt 0 then trueanom[ind] = 2 * !pi - trueanom[ind] 
;correction for trueanom bounds


;;;;; Calculate semi-major axis of star's orbit.                             
astar = aplan * mplan / mstar


;;;;;Calculate Radial Velocity as a function mstar, mplan, aplan, ecc,
;;;;;inc, arg, and time. 
n = 2 * !pi / persec  ;;;;for simplicity's sake

radvel= n * astar * sin(inc) / sqrt(1 - ecc^2) * (cos(trueanom + arg) + $
                                                 ecc * cos(arg)) ;in m/s  

return, radvel



END
