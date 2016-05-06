function readx,filename,ncol=ncol,nlin=nlin,status=status

status=0
nlin=0L
spawn,'wc -l '+filename,help
reads,help,nlin
if (nlin eq 0) then begin
  status=1
  print,'error: file '+filename+' is empty'
  return,0
endif
ntot=0L
spawn,'wc -w '+filename,help
reads,help,ntot
if ((ntot mod nlin) eq 0) then begin
  ncol=ntot/nlin
  q=dblarr(ncol,nlin)
  openr,37,filename
  readf,37,q
  close,37
endif else begin
  status=2
  print,'error: readx cannot read from file '+filename
endelse

return,q

end

