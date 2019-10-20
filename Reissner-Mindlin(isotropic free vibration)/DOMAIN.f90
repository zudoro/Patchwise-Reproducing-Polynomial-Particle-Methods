MODULE DOMAIN

    USE NEWTYPE
    USE GLBVAR

    IMPLICIT NONE

CONTAINS

SUBROUTINE GETBGMESH(DOMAIN,BGMESH)

    TYPE(RECPATCH), INTENT(INOUT) :: BGMESH(NUMPATCHES)
    TYPE(RECPATCH), INTENT(IN) :: DOMAIN
    TYPE(POINT2D) :: XEND, YEND
    TYPE(POINT2D), ALLOCATABLE :: XCORD(:), YCORD(:)
    REAL*8 :: HX, HY    
	INTEGER :: I, J, K, NX, NY 
	
	PRINT*, ""
	PRINT*, "CONSTRUCT BACKGROUND MESHES"
	
	IF (NUMPATCHES .GT. INT(1)) THEN
		NX=INT(NUMPATCHES/DIVS(1)%B) ! NUMBER OF RECTANGULAR PATCHS FOR EACH Y-LAYER
		NY=INT(NUMPATCHES/DIVS(1)%A) ! NUMBER OF RECTANGULAR PATCHS FOR EACH X-LAYER
	
		ALLOCATE(XCORD(NX), YCORD(NY))
	
		! X COMPONENTS OF BOUNDARY
		XEND%X=DOMAIN%PT1%X
		XEND%Y=DOMAIN%PT2%X
	
		! Y COMPONENTS OF BOUNDARY
		YEND%X=DOMAIN%PT1%Y
		YEND%Y=DOMAIN%PT3%Y
	
		HX=(1.0D0/NX)*(XEND%Y-XEND%X)	! MESHSIZE ALONG X-AXIS
		HY=(1.0D0/NY)*(YEND%Y-YEND%X)	! MESHSIZE ALONG Y-AXIS
	
		! X AND Y COORDINATES OF PRE-BGMESHS
		DO I=1, NX
			XCORD(I)%X=XEND%X+(I-1)*HX
			XCORD(I)%Y=XCORD(I)%X+HX
		ENDDO
		DO I=1, NY
			YCORD(I)%X=YEND%X+(I-1)*HY
			YCORD(I)%Y=YCORD(I)%X+HY
		ENDDO
	
		DO I=1, NY	! CORRESPONDING TO Y COMPONENT
			DO J=1, NX	! CORRESPONDING TO X COMPONENT
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT1%X = XCORD(J)%X
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT4%X = XCORD(J)%X
			
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT2%X = XCORD(J)%Y
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT3%X = XCORD(J)%Y
			
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT1%Y = YCORD(I)%X
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT2%Y = YCORD(I)%X
			
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT3%Y = YCORD(I)%Y
				BGMESH(MOD(J,NX+1)+(I-1)*NX)%PT4%Y = YCORD(I)%Y
			ENDDO
		ENDDO
		
		DEALLOCATE(XCORD,YCORD)
			
	ELSE
		BGMESH(1) = DOMAIN
	ENDIF

	IF (NUMPATCHES .GT. INT(1)) THEN
		DO I=1, NUMPATCHES
			IF (MOD(I,NX) .EQ. 1) THEN
				BGMESH(I)%PT1%X = BGMESH(I)%PT1%X - DELTA
				BGMESH(I)%PT4%X = BGMESH(I)%PT4%X - DELTA
			ENDIF
			IF (MOD(I,NX) .EQ. 0) THEN
				BGMESH(I)%PT2%X = BGMESH(I)%PT2%X + DELTA
				BGMESH(I)%PT3%X = BGMESH(I)%PT3%X + DELTA
			ENDIF
			IF (I .LE. NX) THEN
				BGMESH(I)%PT1%Y = BGMESH(I)%PT1%Y - DELTA
				BGMESH(I)%PT2%Y = BGMESH(I)%PT2%Y - DELTA
			ENDIF
			IF (I .GT. NX*(NY-1)) THEN
				BGMESH(I)%PT3%Y = BGMESH(I)%PT3%Y + DELTA
				BGMESH(I)%PT4%Y = BGMESH(I)%PT4%Y + DELTA
			ENDIF
		ENDDO
	ELSE
		BGMESH(1) = BGMESH(1)+DELTA
	ENDIF
END SUBROUTINE GETBGMESH

SUBROUTINE GETDOMAIN(OMEGA)
	
TYPE(RECPATCH), INTENT(INOUT) :: OMEGA

    PRINT *,
    PRINT *, 'CONSTRUCT DOMAIN'
    IF (PROBLEM.EQ.0) THEN
		OMEGA%PT1=POINT2D(0.0D0, 0.0D0)
		OMEGA%PT2=POINT2D(L, 0.0D00)
		OMEGA%PT3=POINT2D(L, L)
		OMEGA%PT4=POINT2D(0.0D0, L)
	ELSEIF (PROBLEM.EQ.1) THEN
		OMEGA%PT1=POINT2D(0.0D0, 0.0D0)
		OMEGA%PT2=POINT2D(1.0D0, 0.0D00)
		OMEGA%PT3=POINT2D(1.0D0, 1.0D0)
		OMEGA%PT4=POINT2D(0.0D0, 1.0D0)
    ELSEIF (PROBLEM.EQ.2) THEN
		OMEGA%PT1=POINT2D(0.0D0, 0.0D0)
		OMEGA%PT2=POINT2D(2.0D0, 0.0D00)
		OMEGA%PT3=POINT2D(2.0D0, 2.0D0)
		OMEGA%PT4=POINT2D(0.0D0, 2.0D0)
    ELSEIF (PROBLEM.EQ.3) THEN
    	OMEGA%PT1=POINT2D(-1.0D0, 0.0D0)
		OMEGA%PT2=POINT2D(1.0D0, 0.0D00)
		OMEGA%PT3=POINT2D(1.0D0, 1.0D0)
		OMEGA%PT4=POINT2D(-1.0D0, 1.0D0)
    ENDIF
END SUBROUTINE GETDOMAIN

 END MODULE DOMAIN

