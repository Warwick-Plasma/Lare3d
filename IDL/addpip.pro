PRO addpip, data, perp=perp, par=par, eta_perp=eta_perp

on_error,2
close,1

IF N_ELEMENTS(data) NE 0 THEN BEGIN
    nx = data.nx+1
    ny = data.ny+1
    nz = data.nz+1
    prec = data.prec

    header_offset = 16 + prec * (nx+ny+nz+1)
    array_offset = prec * nx * ny * nz
    
    ; par_current, perp_current, eta_perp
    offset = header_offset + 8 * array_offset

    all = NOT(KEYWORD_SET(perp) OR KEYWORD_SET(par) OR KEYWORD_SET(eta))

    openr,1,data.filename
    pipvars = assoc(1,(prec EQ 4) ? fltarr(nx,ny,nz) : dblarr(nx,ny,nz),offset,/packed)
    IF KEYWORD_SET(perp) OR (all) THEN data = CREATE_STRUCT(data,{perp_current: pipvars[1]})
    IF KEYWORD_SET(par) OR (all) THEN data = CREATE_STRUCT(data,{par_current: pipvars[0]})
    IF KEYWORD_SET(eta) OR (all) THEN data = CREATE_STRUCT(data,{eta_perp: pipvars[2]})
    close,1
ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addpip, <data sturcture>[, /perp, /par, /eta_perp]"
ENDELSE

END

