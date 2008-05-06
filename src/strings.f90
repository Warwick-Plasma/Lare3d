MODULE strings

  USE shared_data

  IMPLICIT NONE

CONTAINS
  FUNCTION StrCmp(StrIn,StrTest)

    CHARACTER(*),INTENT(IN) ::  StrIn,StrTest
    CHARACTER(30) :: StrTrim
    LOGICAL :: StrCmp

    StrTrim=TRIM(ADJUSTL(StrIn))

    IF (LEN(StrTest) .GT. LEN(StrIn)) THEN
       StrCmp=.FALSE.
       return
    ENDIF

    IF (StrTrim(LEN(StrTest)+1:LEN(StrTest)+1) .NE. " ") THEN
       StrCmp=.FALSE.
       RETURN
    ENDIF

    StrCmp=StrTrim(1:Len(StrTest)) == StrTest

  END FUNCTION StrCmp

  FUNCTION AsInteger(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsInteger,Value,f
    READ(unit=StrIn,fmt=*,iostat=f) value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsInteger=value

  END FUNCTION AsInteger

  FUNCTION AsReal(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER::f
    REAL(num) :: AsReal
    REAL(num) :: Value

    READ(unit=StrIn,fmt=*,iostat=f) Value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsReal=Value
  END FUNCTION AsReal

  FUNCTION AsLogical(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    LOGICAL :: AsLogical

    AsLogical=.FALSE.
    IF(StrCmp(TRIM(ADJUSTL(StrIn)),"T")) THEN
       AsLogical=.TRUE.
       RETURN
    ENDIF
    IF(StrCmp(TRIM(ADJUSTL(StrIn)),"F")) THEN
       AsLogical=.FALSE.
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsLogical

  FUNCTION AsBC(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsBC

    AsBC=-1

    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"periodic")) THEN
       AsBC=periodic
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"open")) THEN
       AsBC=open
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"other")) THEN
       AsBC=other
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsBC
END MODULE strings
