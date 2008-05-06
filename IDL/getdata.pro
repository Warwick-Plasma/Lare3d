FUNCTION getdata, snapshot, wkdir=wkdir, rho=rho, temp=temp, vx=vx,vy=vy,vz=vz, $
                  bx=bx,by=by,bz=bz,fields=fields,vel=vel, empty=empty

on_error, 2
IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

IF N_PARAMS() EQ 0 THEN BEGIN
    print, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
ENDIF
header = assoc(1,{nx: 0L, ny: 0L, nz: 0L, prec: 0L},0,/packed)
file = wkdir + string(snapshot,format='("/",I04,".llld")')
close,1
openr,1,file
h = header[0]
; set the array sizes
nx = h.nx+1
ny = h.ny+1
nz = h.nz+1
prec = h.prec

header_offset = 16 + prec * (nx+ny+nz + 1) ; the +1 is for time, 16 for the 4 x 4 Byte ints at start of file.
array_offset = prec * nx * ny * nz

IF (prec EQ 4) THEN  BEGIN
    header = assoc(1,{nx: 0L, ny: 0L, nz: 0L, prec: 0L, time: 0.0, x: fltarr(nx), y: fltarr(ny), z: fltarr(nz)$
                     },0,/packed)
    data = assoc(1,fltarr(nx,ny,nz),header_offset,/packed)
ENDIF ELSE BEGIN
    header = assoc(1,{nx: 0L, ny: 0L, nz: 0L, prec: 0L, time: 0.0D, x: dblarr(nx), y: dblarr(ny), z: dblarr(nz)$
                     },0,/packed)
    data = assoc(1,dblarr(nx,ny,nz),header_offset,/packed)
ENDELSE

f = {filename: file}
d = CREATE_STRUCT(f, header[0])

IF NOT KEYWORD_SET(empty) THEN BEGIN
    
    IF KEYWORD_SET(rho) THEN d = CREATE_STRUCT(d, {rho: data[0]})
    IF KEYWORD_SET(temp) THEN d = CREATE_STRUCT(d, {temp: data[1]})
    IF KEYWORD_SET(vx) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vx: data[2]})
    IF KEYWORD_SET(vy) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vy: data[3]})
    IF KEYWORD_SET(vz) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vz: data[4]})
    IF KEYWORD_SET(bx) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {bx: data[5]})
    IF KEYWORD_SET(by) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {by: data[6]})
    IF KEYWORD_SET(bz) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {bz: data[7]})

    IF N_TAGS(d) EQ 9 THEN BEGIN ; no data loaded so default to all
        IF (prec EQ 4) THEN BEGIN
            data = assoc(1,{nx: 0L, ny: 0L, nz: 0L, prec: 0L, time: 0.0, x: fltarr(nx), y: fltarr(ny), z: fltarr(nz),$
                            rho: fltarr(nx,ny,nz), temp: fltarr(nx,ny,nz), vx: fltarr(nx,ny,nz),$
                            vy: fltarr(nx,ny,nz), vz: fltarr(nx,ny,nz), bx: fltarr(nx,ny,nz),$
                            by: fltarr(nx,ny,nz), bz: fltarr(nx,ny,nz)},0,/packed)
            d = CREATE_STRUCT(f,data[0])
        ENDIF ELSE BEGIN
            data = assoc(1,{nx: 0L, ny: 0L, nz: 0L, prec: 0L, time: 0.0D, x: dblarr(nx), y: dblarr(ny), z: dblarr(nz),$
                            rho: dblarr(nx,ny,nz), temp: dblarr(nx,ny,nz), vx: dblarr(nx,ny,nz),$
                            vy: dblarr(nx,ny,nz), vz: dblarr(nx,ny,nz), bx: dblarr(nx,ny,nz),$
                            by: dblarr(nx,ny,nz), bz: dblarr(nx,ny,nz)},0,/packed)
            d = CREATE_STRUCT(f,data[0])
        ENDELSE
    ENDIF
 
ENDIF

close,1

RETURN, d

END
