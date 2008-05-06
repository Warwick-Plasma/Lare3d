; This version of getdata is for use with the old, many file output of
; Lare3d.

FUNCTION getdata, snapshot, wkdir=wkdir, rho=rho, temp=temp, vx=vx,vy=vy,vz=vz, $
                  bx=bx,by=by,bz=bz,fields=fields,vel=vel, empty=empty
COMMON control, nx, ny, nz, nx_local, ny_local, nz_local, nstep, nxp, nyp, nzp
COMMON grid, x, y, z
COMMON thetime, time

on_error, 2
IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

IF N_PARAMS() EQ 0 THEN BEGIN
    print, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
ENDIF

loadgrid,wkdir=wkdir,snapshot=snapshot

nx = nx+1
ny = ny+1
nz = nz+1

f = {snapshot: snapshot}
d = CREATE_STRUCT(f, {nx: nx, ny: ny, nz: nz, time: time, x: x, y:y, z:z})

IF NOT KEYWORD_SET(empty) THEN BEGIN
    IF KEYWORD_SET(rho) THEN d = CREATE_STRUCT(d, {rho: getrho(snapshot)})
    IF KEYWORD_SET(temp) THEN d = CREATE_STRUCT(d, {temp: gettemp(snapshot)})
    IF KEYWORD_SET(vx) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vx: getvx(snapshot)})
    IF KEYWORD_SET(vy) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vy: getvy(snapshot)})
    IF KEYWORD_SET(vz) OR KEYWORD_SET(vel) THEN d = CREATE_STRUCT(d, {vz: getvz(snapshot)})
    IF KEYWORD_SET(bx) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {bx: getbx(snapshot)})
    IF KEYWORD_SET(by) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {by: getby(snapshot)})
    IF KEYWORD_SET(bz) OR KEYWORD_SET(fields) THEN d = CREATE_STRUCT(d, {bz: getbz(snapshot)})

    IF N_TAGS(d) EQ 8 THEN BEGIN ; no data loaded so default to all
        d = CREATE_STRUCT(d,{$
                        rho: getrho(snapshot), temp: gettemp(snapshot), vx: getvx(snapshot),$
                        vy: getvy(snapshot), vz: getvz(snapshot), bx: getbx(snapshot),$
                        by: getby(snapshot), bz: getbz(snapshot)})
    ENDIF
 
ENDIF

close,1

RETURN, d

END


PRO addpip, data, perp=perp_current, par=par_current, eta=eta_perp

on_error,2
close,1
all = 1
IF (KEYWORD_SET(perp) OR KEYWORD_SET(par) OR KEYWORD_SET(eta)) THEN all = 0

IF N_ELEMENTS(data) NE 0 THEN BEGIN
    ; par_current, perp_current, eta_perp

    IF KEYWORD_SET(perp) OR (all EQ 1) THEN data = CREATE_STRUCT(data,{perp_current: getperp(data.snapshot)})
    IF KEYWORD_SET(par) OR (all EQ 1) THEN data = CREATE_STRUCT(data,{par_current: getpar(data.snapshot)})
    IF KEYWORD_SET(eta) OR (all EQ 1) THEN data = CREATE_STRUCT(data,{eta_perp: geteta(data.snapshot)})

ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addpip, <data sturcture>[, /perp_current, /par_current, /eta_perp]"
ENDELSE

END



;
; General function for loading the data, is never called directly by
; the user
;


FUNCTION lddata, variable, snapshot, wkdir

file = wkdir + string(0,snapshot,format='("/",I3.3,I4.4,".dat")')
openr,1,file

nx=0.0
ny=0.0
nz=0.0
nx_local=0.0
ny_local=0.0
nz_local=0.0

; Binary test
t = 1.0
readu,1,t
close,1

binary = (FIX(t) NE 0.0)
IF binary THEN BEGIN
    openr,1,file
    readu,1,nx,ny,nz
    readu,1,nx_local,ny_local,nz_local
    readu,1,time
    close,1
ENDIF ELSE BEGIN
    lw = 1L
    openr,1,file
    readu,1,lw,nx,ny,nz,lw
    readu,1,lw,nx_local,ny_local,nz_local,lw
    readu,1,lw,time,lw
    close,1
ENDELSE

nx = fix(nx)+1
ny = fix(ny)+1
nz = fix(nz)+1
nx_local = fix(nx_local)+1
ny_local = fix(ny_local)+1
nz_local = fix(nz_local)+1
nxp = nx / nx_local
nyp = ny / ny_local
nzp = nz / nz_local

var = fltarr(nx,ny,nz)

IF binary THEN BEGIN
    data = ASSOC(10,fltarr(nx_local,ny_local,nz_local),28,/PACKED)
    
    FOR k =  0, nzp - 1 DO BEGIN 
        FOR j =  0, nyp - 1 DO BEGIN
            FOR i =  0, nxp - 1 DO BEGIN
                rank = i + j*nxp + k*nxp*nyp
                file =  wkdir + string(rank,snapshot,format='("/",I3.3,I4.4,".dat")')
                openr,10,file
                x1 =  nx_local * i
                x2 =  nx_local * (i+1) - 1
                y1 =  ny_local * j
                y2 =  ny_local * (j+1) - 1
                z1 =  nz_local * k
                z2 =  nz_local * (k+1) - 1
                var(x1:x2, y1:y2, z1:z2) = data[variable]
                close,10
            ENDFOR
        ENDFOR
    ENDFOR
ENDIF ELSE BEGIN
    datasize = LONG(nx_local) * LONG(ny_local) * LONG(nz_local) * 4L
    coords = (2L * (nx_local + ny_local + nz_local) + 3L) * 4L
    CASE variable OF
        0: offset = 56L ; rho
        1: offset = 56L + datasize ; temp
        2: offset = 64L + 2L * datasize ; vx
        3: offset = 64L + 3L * datasize ; vy
        4: offset = 64L + 4L * datasize ; vz
        5: offset = 72L + 5L * datasize ; bx
        6: offset = 72L + 6L * datasize ; by
        7: offset = 72L + 7L * datasize ; bz
        8: offset = 96L + 8L * datasize + coords ; parallel current
        9: offset = 96L + 9L * datasize + coords ; perpendicular current
        10: offset = 96L + 10L * datasize + coords ; eta_perp
    ENDCASE

    data = ASSOC(10,fltarr(nx_local,ny_local,nz_local),offset,/PACKED)

    FOR k =  0, nzp - 1 DO BEGIN 
        FOR j =  0, nyp - 1 DO BEGIN
            FOR i =  0, nxp - 1 DO BEGIN
                rank = i + j*nxp + k*nxp*nyp
                file =  wkdir + string(rank,snapshot,format='("/",I3.3,I4.4,".dat")')
                openr,10,file
                x1 =  nx_local * i
                x2 =  nx_local * (i+1) - 1
                y1 =  ny_local * j
                y2 =  ny_local * (j+1) - 1
                z1 =  nz_local * k
                z2 =  nz_local * (k+1) - 1
                var(x1:x2, y1:y2, z1:z2) = data[0]
                close,10
            ENDFOR
        ENDFOR
    ENDFOR


ENDELSE

RETURN, var

END

;
; All the functions for loading the variables, they all take the form:
; get<variable> ( <snapshot> [, wkdir=<wkdir>] )
;

FUNCTION getrho, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(0,snapshot,wkdir)

  RETURN, rho

END


FUNCTION gettemp, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(1,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getvx, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(2,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getvy, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(3,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getvz, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(4,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getbx, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(5,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getby, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(6,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getbz, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(7,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getpar, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(8,snapshot,wkdir)

  RETURN, rho

END


FUNCTION getperp, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(9,snapshot,wkdir)

  RETURN, rho

END


FUNCTION geteta, snapshot, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

  rho = lddata(10,snapshot,wkdir)

  RETURN, rho

END

;
; fields loads in just the magnetic fields
;

PRO fields, snapshot
COMMON field, bx, by, bz, curlb
COMMON data_dir,  wkdir

on_error,2

IF N_ELEMENTS(snapshot) EQ 0 THEN BEGIN
    PRINT, 'Please use: fields, <snapshot>'
ENDIF ELSE BEGIN
    bx = getbx(snapshot,wkdir=wkdir)
    by = getby(snapshot,wkdir=wkdir)
    bz = getbz(snapshot,wkdir=wkdir)
ENDELSE

END

;
; velocities loads in just the velocities
;

PRO velocities, snapshot
COMMON fluid, vx, vy, vz, rho, temp

on_error,2


IF N_ELEMENTS(snapshot) EQ 0 THEN BEGIN
    PRINT, 'Please use: velocities, <snapshot>'
ENDIF ELSE BEGIN
    vx = getvx(snapshot,wkdir=wkdir)
    vy = getvy(snapshot,wkdir=wkdir)
    vz = getvz(snapshot,wkdir=wkdir)
ENDELSE

END

;
; loadgrid loads in the x,y,z arrays and sets nx,ny,nz etc.
;


PRO loadgrid, wkdir=theworkdir, snapshot=snapshot
COMMON control, nx, ny, nz, nx_local, ny_local, nz_local, nstep, nxp, nyp, nzp
COMMON grid, x, y, z
COMMON data_dir,  wkdir
COMMON thetime, time

on_error,2

IF NOT KEYWORD_SET(theworkdir) THEN theworkdir = wkdir
IF NOT KEYWORD_SET(snapshot) THEN snapshot = 0

file = theworkdir + string(0,snapshot,format='("/",I3.3,I4.4,".dat")')
openr,1,file

nx=0.0
ny=0.0
nz=0.0
nx_local=0.0
ny_local=0.0
nz_local=0.0

; Binary test
t = 1.0
readu,1,t
close,1

binary = (FIX(t) NE 0.0)
IF binary THEN BEGIN
    openr,1,file
    readu,1,nx,ny,nz
    readu,1,nx_local,ny_local,nz_local
    readu,1,time
    close,1
ENDIF ELSE BEGIN
    lw = 1L
    openr,1,file
    readu,1,lw,nx,ny,nz,lw
    readu,1,lw,nx_local,ny_local,nz_local,lw
    readu,1,lw,time,lw
    close,1
ENDELSE

nx = fix(nx)+1
ny = fix(ny)+1
nz = fix(nz)+1
nx_local = fix(nx_local)+1
ny_local = fix(ny_local)+1
nz_local = fix(nz_local)+1
nxp = nx / nx_local
nyp = ny / ny_local
nzp = nz / nz_local

x = fltarr(nx)
y = fltarr(ny)
z = fltarr(nz)
datasize = LONG(nx_local) * LONG(ny_local) * LONG(nz_local) * 4L

IF binary THEN BEGIN
    offset = 28L + 8L * datasize
ENDIF ELSE BEGIN
    offset = 80L + 8L * datasize
ENDELSE

coords = ASSOC(10,{x:fltarr(nx_local),y:fltarr(ny_local),z:fltarr(nz_local)},offset,/PACKED)

FOR k =  0, nzp - 1 DO BEGIN 
    FOR j =  0, nyp - 1 DO BEGIN
        FOR i =  0, nxp - 1 DO BEGIN
            rank = i + j*nxp + k*nxp*nyp
            file =  theworkdir + string(rank,snapshot,format='("/",I3.3,I4.4,".dat")')
            openr,10,file
            x1 =  nx_local * i
            x2 =  nx_local * (i+1) - 1
            y1 =  ny_local * j
            y2 =  ny_local * (j+1) - 1
            z1 =  nz_local * k
            z2 =  nz_local * (k+1) - 1
            thecoords = coords[0]
            x(x1:x2) = thecoords.x
            y(y1:y2) = thecoords.y
            z(z1:z2) = thecoords.z
            close,10
        ENDFOR
    ENDFOR
ENDFOR

END


