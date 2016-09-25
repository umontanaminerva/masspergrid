pro printpars, inpars, fixed, tied_in, message = message, quiet = quiet, errors = errors, tts = tt
  pars = inpars
  siz = size(pars)
  fix = n_elements(fixed) gt 0
  if ~fix then fixed = fix(inpars*0)
  if n_elements(tied_in) gt 0 then tied = tied_in else tied = fix(inpars*0)
  if siz[2] eq 8 then begin
    pars = pars(pars, sigpars)
    if n_elements(errors) eq 0 then errors = sigpars
    b = where(~finite(errors), nb)
    if nb gt 0 then errors[b] = 0
    b = where(~finite(pars), nb)
    if nb gt 0 then pars[b] = 0
  endif
  if keyword_set(quiet) then return
  names = ["Per", "Tp", "e", "omega", "K", "gamma", "dvdt"]
  if n_elements(message) gt 0 then print, message
  np = n_elements(pars)/7
  ntt = n_elements(tt)
;  if ntt gt 0 then tied[1+indgen(ntt)*7] = 1
  tie =  total(fixed) gt 0 or ntt gt 0
  errs = n_elements(errors)
  s = strarr(7)
  for i = 0, np-1 do begin
    for j = 0, 4 do begin
      ind = j+i*7
      snp = np gt 1 ? strtr(i+1)+' ' : ''
      if i gt 0 then snp = '  '+snp
      if fix then begin
        fixs = ' ('+(fixed[ind]?"Fixed)":"Free) ")
      endif else fixs = '' 
      if tie then begin
        if tied[ind] gt 0 then fixs = ' (Tied) '
        if i lt ntt and j eq 3 then fixs = ' (Tied) '
      endif
      if errs then err = sigfig(pars[ind], errors[ind], /par) else err = sigfig(pars[ind], pars[ind]/1d5)
      s[j] = s[j]+string(snp+names[ind mod 7], ': ', err, fixs, format = '(a-10,a,a15,a)')
    endfor
  endfor

  for i = 0, 4 do print, s[i]  

  for i = 0, n_elements(tt)-1 do begin
    if i eq 0 then s7 = ''
    snp = np gt 1 ? strtr(i+1)+' ' : ''
    if i gt 0 then snp = '  '+snp
    s7 = s7+string(snp+'t_transit', ': ', sigfig(tt[i], tt[i]/1d5), ' (Fixed)', format = '(a-10,a,a15,a)')
  endfor
  
  if n_elements(tt) gt 0 then print, s7

  for i = 5, 6 do begin
    if fix then begin
      fixs = ' ('+(fixed[i]?"Fixed)":"Free) ")
    endif else fixs = ''
    if errs then err = sigfig(pars[i], errors[i], /par) else err = sigfig(pars[i], pars[i]/1d5)
    s[i] = s[i]+string(names[i mod 7], ': ', err, fixs, format = '(a-10,a,a15,a)')
  endfor

  for i = 5, 6 do print, s[i]  
  
      
end

