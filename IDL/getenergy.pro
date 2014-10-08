FUNCTION getenergy, wkdir_in, _EXTRA=extra

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  retro = -1
  wkdir = ''

  IF (N_ELEMENTS(wkdir_in) NE 0) THEN BEGIN
    wkdir = wkdir_in
  ENDIF

  IF (KEYWORD_SET(extra)) THEN BEGIN
    extra_tags = TAG_NAMES(extra)
    gotextra = 0
    FOR i = 0, N_ELEMENTS(extra_tags)-1 DO BEGIN
      CASE extra_tags[i] OF
        'WKDIR': wkdir = extra.wkdir
        'RETRO': retro = extra.retro
        ELSE: BEGIN
          IF (gotextra EQ 0) THEN BEGIN
            new_extra = CREATE_STRUCT(extra_tags[i], extra.(i))
            gotextra = 1
          ENDIF ELSE BEGIN
            new_extra = CREATE_STRUCT(new_extra, extra_tags[i], extra.(i))
          ENDELSE
        END
      ENDCASE
    ENDFOR
  ENDIF

  IF (wkdir EQ '') THEN wkdir = wkdir_global
  IF (retro EQ -1) THEN retro = retro_global

  file = wkdir + '/en.dat'
  info = FILE_INFO(file)

  IF info.exists EQ 0 THEN BEGIN
    PRINT, "ERROR: Missing en.dat file."
    RETURN, "ERROR: Missing en.dat file."
  ENDIF

  OPENR, lun, file, /GET_LUN
  file_info = FSTAT(lun)
  fileint = ASSOC(lun, LONARR(1), 0, /packed)

  ; Set the size of the variables in bytes
  prec = REFORM(fileint[0])

  ; Set the number of variables (including time) in the energy data.
  ; If you change this you also need to alter the structure below
  nv = REFORM(fileint[1])

  ; Calculate the number of outputs (the -8 is due to the prec/columns output)
  outs = (file_info.size - 8) / nv / prec

  ; Again 8 offset due to prec/columns
  energy_mask = ASSOC(lun, $
      (prec EQ 4) ? FLTARR(nv,outs) : DBLARR(nv,outs) , 8, /packed)
  data = energy_mask[0]

  energy = {points:outs, time:REFORM(data[0,*]), en_b:REFORM(data[1,*]), $
            en_ke:REFORM(data[2,*]), en_int:REFORM(data[3,*]), $
            heating_visc:REFORM(data[4,*]), heating_ohmic:REFORM(data[5,*])}

  CLOSE, lun

  RETURN, energy
END
