MODULE ERRORESTIMATE

	USE NEWTYPE
	USE GLBVAR
	USE EXACTSOL
	USE PACKAGE
	USE LOADFUNCTION

    IMPLICIT NONE
    
CONTAINS

!----MAX. NORM ESTIMATE----
SUBROUTINE MAXNORM(COEFF_SOL)

	REAL*8, INTENT(IN) :: COEFF_SOL(3*DOF)
	TYPE(POINT2D) :: MESH
	TYPE(INT2D) :: MAXGRID
	INTEGER :: II, I, J, NDX(3,DOF), INDICE(3), GRIDN
	TYPE(DISPLACEMENT), ALLOCATABLE :: SOL(:), DISP(:), ERROR_DISP(:)
	TYPE(MOMENT), ALLOCATABLE :: MOMENTS(:)
	REAL*8, ALLOCATABLE :: XCORDS(:), YCORDS(:)
	
    CALL STIFINDEX(NDX)
    
	MAXGRID = INT2D(100,100)
	
	ALLOCATE(XCORDS(MAXGRID%A+1), YCORDS(MAXGRID%B+1))
	GRIDN = UBOUND(XCORDS,1)*UBOUND(YCORDS,1)

	ALLOCATE(ERROR_DISP(GRIDN), MOMENTS(GRIDN))
	ALLOCATE(SOL(GRIDN), DISP(GRIDN))
		
	MESH=POINT2D((OMEGA%PT2%X-OMEGA%PT1%X)/MAXGRID%A, (OMEGA%PT4%Y-OMEGA%PT1%Y)/MAXGRID%B)
	
	DO II=1, MAXGRID%A+1
		XCORDS(II)=OMEGA%PT1%X+MESH%X*(II-1)
	ENDDO
	
	DO II=1, MAXGRID%B+1
		YCORDS(II)=OMEGA%PT1%Y+MESH%Y*(II-1)
	ENDDO
	!----- COMPUTE FOURIER SERIES OF ANALYTIC SOLUTION -----
	CALL ANALYTICSOL(SOL, XCORDS, YCORDS)
	!----- COMPUTE PUFEM APPROXIMATE SOLUTION -----
	CALL APPROXSOL(DISP, XCORDS, YCORDS, COEFF_SOL, NDX)
	!----- ESTIMATE MAXIMUM NORM ERROR -----
	DO I=1, GRIDN
		IF (DABS(SOL(I)%W0) .LE. EPS) THEN
			ERROR_DISP(I)%W0 = DABS(DISP(I)%W0 - SOL(I)%W0)
		ELSE
			ERROR_DISP(I)%W0 = DABS((DISP(I)%W0 - SOL(I)%W0)/SOL(I)%W0)
		ENDIF
		IF (DABS(SOL(I)%PIX) .LE. EPS) THEN
			ERROR_DISP(I)%PIX = DABS(DISP(I)%PIX - SOL(I)%PIX)
		ELSE
			ERROR_DISP(I)%PIX = DABS((DISP(I)%PIX - SOL(I)%PIX)/SOL(I)%PIX)
		ENDIF
		IF (DABS(SOL(I)%PIY) .LE. EPS) THEN
			ERROR_DISP(I)%PIY = DABS(DISP(I)%PIY - SOL(I)%PIY)
		ELSE
			ERROR_DISP(I)%PIY = DABS((DISP(I)%PIY - SOL(I)%PIY)/SOL(I)%PIY)
		ENDIF
	ENDDO
	
	!----- INDICE WHICH INDICATE MAXIMUM NORM ERROR AT EACH ERROR
	INDICE = (/ MAXLOC(ERROR_DISP(:)%W0), MAXLOC(ERROR_DISP(:)%PIX), MAXLOC(ERROR_DISP(:)%PIY) /)
	
	OPEN(1, FILE='exact_w0.dat')
	OPEN(2, FILE='exact_pix.dat')
	OPEN(3, FILE='exact_piy.dat')
	OPEN(4, FILE='exact_Mxx.dat')
	OPEN(5, FILE='exact_Myy.dat')
	OPEN(6, FILE='exact_Mxy.dat')
	OPEN(7, FILE='approx_w0.dat')
	OPEN(8, FILE='approx_pix.dat')
	OPEN(9, FILE='approx_piy.dat')
	OPEN(10, FILE='error_w0.dat')
	OPEN(11, FILE='error_pix.dat')
	OPEN(12, FILE='error_piy.dat')
	II=1;
	DO I=1, MAXGRID%A+1
		DO J=1, MAXGRID%B+1
			WRITE(1,*) XCORDS(I), YCORDS(J), SOL(II)%W0
			WRITE(2,*) XCORDS(I), YCORDS(J), SOL(II)%PIX
			WRITE(3,*) XCORDS(I), YCORDS(J), SOL(II)%PIY
			WRITE(4,*) XCORDS(I), YCORDS(J), MOMENTS(II)%MXX
			WRITE(5,*) XCORDS(I), YCORDS(J), MOMENTS(II)%MYY
			WRITE(6,*) XCORDS(I), YCORDS(J), MOMENTS(II)%MXY
			WRITE(7,*) XCORDS(I), YCORDS(J), DISP(II)%W0
			WRITE(8,*) XCORDS(I), YCORDS(J), DISP(II)%PIX
			WRITE(9,*) XCORDS(I), YCORDS(J), DISP(II)%PIY
			WRITE(10,*) XCORDS(I), YCORDS(J), ERROR_DISP(II)%W0
			WRITE(11,*) XCORDS(I), YCORDS(J), ERROR_DISP(II)%PIX
			WRITE(12,*) XCORDS(I), YCORDS(J), ERROR_DISP(II)%PIY
			II=II+1
		ENDDO
		WRITE(1,*) ""
		WRITE(2,*) ""
		WRITE(3,*) ""
		WRITE(4,*) ""
		WRITE(5,*) ""
		WRITE(6,*) ""
		WRITE(7,*) ""
		WRITE(8,*) ""
		WRITE(9,*) ""
		WRITE(10,*) ""
		WRITE(11,*) ""
		WRITE(12,*) ""
	ENDDO
	CLOSE(1);CLOSE(2);CLOSE(3);CLOSE(4);CLOSE(5);CLOSE(6)
	CLOSE(7);CLOSE(8);CLOSE(9);CLOSE(10);CLOSE(11);CLOSE(12)

	OPEN(13, FILE='error.dat')
	WRITE(13,*) "RELATIVE MAXIMUM ERROR FOR W0 : "
	WRITE(13,*) ERROR_DISP(INDICE(1))%W0
	WRITE(13,*) ""
	WRITE(13,*) "RELATIVE MAXIMUM ERROR FOR PI_X : "
	WRITE(13,*) ERROR_DISP(INDICE(2))%PIX
	WRITE(13,*) ""
	WRITE(13,*) "RELATIVE MAXIMUM ERROR FOR PI_Y : "
	WRITE(13,*) ERROR_DISP(INDICE(3))%PIY
	CLOSE(13)

	PRINT*, "RELATIVE MAXIMUM ERROR FOR W0 : "
	PRINT*, ERROR_DISP(INDICE(1))%W0
	PRINT*, ""
	PRINT*, "RELATIVE MAXIMUM ERROR FOR PI_X : "
	PRINT*, ERROR_DISP(INDICE(2))%PIX
	PRINT*, ""
	PRINT*, "RELATIVE MAXIMUM ERROR FOR PI_Y : "
	PRINT*, ERROR_DISP(INDICE(3))%PIY
		
	DEALLOCATE(XCORDS, YCORDS, SOL, MOMENTS, DISP, ERROR_DISP)

END SUBROUTINE MAXNORM

END MODULE ERRORESTIMATE
