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

  magic = 'mmm'
  READU, lun, magic

  IF magic EQ 'HIS' THEN BEGIN
    fileint = ASSOC(lun, LONARR(1), 3, /packed)

    endianness = REFORM(fileint[2])
    header_length = REFORM(fileint[3])

    ; Set the size of the variables in bytes
    prec = REFORM(fileint[4])

    ; Set the number of variables (including time) in the energy data.
    nv = REFORM(fileint[5])

    id_length = REFORM(fileint[6])

    ; Read in the variable id names
    var_id_start = 3 + 7 * 4
    varbytes = ASSOC(lun, BYTARR(id_length), var_id_start, /packed)
    varnames = STRARR(nv)
    FOR i = 0, nv[0]-1 DO BEGIN
      varnames[i] = STRTRIM(STRING(varbytes[i]), 2)
    ENDFOR

  ENDIF ELSE BEGIN
    fileint = ASSOC(lun, LONARR(1), 0, /packed)

    ; Set the size of the variables in bytes
    prec = REFORM(fileint[0])

    ; Set the number of variables (including time) in the energy data.
    ; If you change this you also need to alter the array below
    nv = REFORM(fileint[1])

    header_length = 8
    varnames = ['time', 'en_b', 'en_ke', 'en_int', $
                'heating_visc', 'heating_ohmic']
  ENDELSE

  ; Calculate the number of outputs
  outs = (file_info.size - header_length) / nv / prec

  energy_mask = ASSOC(lun, $
      (prec EQ 4) ? FLTARR(nv,outs) : DBLARR(nv,outs), header_length, /packed)
  data = energy_mask[0]

  energy = CREATE_STRUCT('points', outs)

  FOR i = 0, nv[0]-1 DO BEGIN
    energy = CREATE_STRUCT(energy, varnames[i], REFORM(data[i,*]))
  ENDFOR

  CLOSE, lun

  RETURN, energy
END
