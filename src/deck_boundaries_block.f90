MODULE deck_boundaries_block

  USE shared_data
  USE strings

SAVE

  INTEGER,PARAMETER :: BoundaryBlockElements=8
  LOGICAL, DIMENSION(BoundaryBlockElements)  :: BoundaryBlockDone
  CHARACTER(len=30),DIMENSION(BoundaryBlockElements) :: BoundaryBlockName=(/"xbc_left","xbc_right","ybc_up","ybc_down","zbc_front","zbc_back","farfield","damping"/)

CONTAINS



  FUNCTION HandleBoundaryDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleBoundaryDeck
    INTEGER :: loop,elementselected

    HandleBoundaryDeck=ERR_NONE

    elementselected=0

    DO loop=1,BoundaryBlockElements
       IF(StrCmp(Element,TRIM(ADJUSTL(BoundaryBlockName(loop))))) THEN
          elementselected=loop
          EXIT
       ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (BoundaryBlockDone(elementselected)) THEN
       HandleBoundaryDeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    BoundaryBlockDone(elementselected)=.TRUE.
    HandleBoundaryDeck=ERR_NONE

    SELECT CASE (elementselected)
    CASE(1)
       !nx
       xbc_left=AsBC(Value,HandleBoundaryDeck)
    CASE(2)
       !ny
       xbc_right=AsBC(Value,HandleBoundaryDeck)
    CASE(3)
       !nprocx
       ybc_up=AsBC(Value,HandleBoundaryDeck)
    CASE(4)
       !nprocy
       ybc_down=AsBC(Value,HandleBoundaryDeck)
    CASE(5)
       zbc_front=AsBC(Value,HandleBoundaryDeck)
    CASE(6)
       zbc_back=AsBC(Value,HandleBoundaryDeck)
    CASE(7)
       farfield=AsLogical(Value,HandleBoundaryDeck)
    CASE(8)
       damping=AsLogical(Value,HandleBoundaryDeck)
    END SELECT

  END FUNCTION HandleBoundaryDeck

END MODULE deck_boundaries_block
