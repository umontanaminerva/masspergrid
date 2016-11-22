
pro massper_start, ident, suffix

; unix calling sequence: nohup nice +15 idl -e "massper_start, 'ident', 'suffix'"
; ident is node dependent, suffix is process dependent
  
; wrapper code to call multiple massper_run instances on various
; nodes, etc.
  
; change ident, suffix for each MSU node

; should be same for each MSU node
  
usenewdata = 0
msu = 1
subdir = '100616'
nit = 10
fapnit = 1000
star = 'HD185144'
targetlist = 'random60'
weighting = 'ha'
eccmax = 0.4
permin = 1.2

massper_run, subdir, ident, nit, fapnit, star, targetlist, weighting, eccmax=eccmax, $
             permin=permin, suffix=suffix, msu=msu, usenewdata=usenewdata


end

