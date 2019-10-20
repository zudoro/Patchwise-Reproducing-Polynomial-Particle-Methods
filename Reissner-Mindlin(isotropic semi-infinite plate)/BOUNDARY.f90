! FILE WRITTEN BY HYUNJU KIM
!
! DATE : TUE 30 2008

MODULE BOUNDARY

    USE GLBVAR
    USE NEWTYPE
	USE PACKAGE
    USE LOADFUNCTION
    
    IMPLICIT NONE

CONTAINS

SUBROUTINE IMPOSEBD(FULL_R, FULL_C, FULL_K, NZ, K, F)
	
	REAL*8, INTENT(OUT) :: FULL_K(:)
    INTEGER, INTENT(OUT) :: FULL_R(:), FULL_C(:)
    INTEGER, INTENT(OUT) :: NZ
	TYPE(SUBMAT), INTENT(INOUT) :: K
    REAL*8, INTENT(INOUT) :: F(3*DOF)
    TYPE(POINT2D) :: NODES(DOF)
    CHARACTER(LEN=2) :: EDGE(DOF)
    INTEGER :: I, J, II, JJ, KK
    LOGICAL :: ITIS
    REAL*8, ALLOCATABLE :: FAT_K(:,:)
    INTEGER, ALLOCATABLE :: INDICE(:)
    LOGICAL, ALLOCATABLE :: BD_MASK(:), FULL_MASK(:,:)

     PRINT *,
     PRINT *, 'IMPOSING BOUNDARY CONDITION'

	ALLOCATE(BD_MASK(3*DOF), FAT_K(3*DOF, 3*DOF), FULL_MASK(3*DOF, 3*DOF))

    NODES=NODE_ALL(NUMPATCHES)

    FAT_K(:,:) = 0.0D0

	DO I=1, DOF
		EDGE(I) = WHICH_BOUNDARY(NODES(I),OMEGA)
	ENDDO

!	OPEN(1, FILE='edge_node.dat')
!	DO I=1, DOF
!		WRITE(1,*) EDGE(I)
!	ENDDO
!	CLOSE(1)

	BD_MASK(:) = .FALSE.

!    DO I=1, DOF
!        ITIS=IS_BOUNDARY(NODES(I),OMEGA)
!         IF (ITIS.EQV..TRUE.) THEN
!            DO J=1, DOF
!                    B(J)=B(J)-U_BD(NODES(I))*FULLENTRY(J,I)
!                    IF (J.EQ.14) THEN
!                        WRITE(350,*) I, B(J), U_BD(NODES(I)), FULLA(J,I)
!                    ENDIF
!            ENDDO
!        ENDIF
!    ENDDO
    
	!---- W0 = 0 ON EVERY BOUNDARIES -----		
	WHERE (EDGE.NE.'NN')
!		F(1:DOF) = 0.0D0
		BD_MASK(1:DOF) = .TRUE.
	END WHERE
	!---- PI_X = 0 ON GAMMA 1 AND 3 -----		
	WHERE (EDGE.EQ.'FF' .OR. EDGE.EQ.'KK' .OR. & 
		   EDGE.EQ.'LF' .OR. EDGE.EQ.'RF' .OR. EDGE.EQ.'LK' .OR. EDGE.EQ.'RK')
!		F(DOF+1:2*DOF) = 0.0D0
		BD_MASK(DOF+1:2*DOF) = .TRUE.
	END WHERE
	!---- PI_Y = 0 ON GAMMA 2 AND 4 -----				
	WHERE (EDGE.EQ.'RR' .OR. EDGE.EQ.'LL' .OR. & 
		   EDGE.EQ.'LF' .OR. EDGE.EQ.'RF' .OR. EDGE.EQ.'LK' .OR. EDGE.EQ.'RK')
!		F(2*DOF+1:3*DOF) = 0.0D0
		BD_MASK(2*DOF+1:3*DOF) = .TRUE.
	END WHERE

	!----- CONSTRUCT FULL MATRIX CALLED FAT_K [K] -----
	FORALL (I=1:K%NZ(1))
		FAT_K(K%R11(I), K%C11(I)) = K%K11(I)
		FAT_K(K%C11(I), K%R11(I)) = K%K11(I)
	END FORALL

	FORALL (I=1:K%NZ(2))
		FAT_K(DOF+K%R22(I), DOF+K%C22(I)) = K%K22(I)
		FAT_K(DOF+K%C22(I), DOF+K%R22(I)) = K%K22(I)
	END FORALL

	FORALL (I=1:K%NZ(3))
		FAT_K(2*DOF+K%R33(I), 2*DOF+K%C33(I)) = K%K33(I)
		FAT_K(2*DOF+K%C33(I), 2*DOF+K%R33(I)) = K%K33(I)
	END FORALL

	FORALL (I=1:K%NZ(4))
		FAT_K(K%R12(I), DOF+K%C12(I)) = K%K12(I)
		FAT_K(DOF+K%C12(I), K%R12(I)) = K%K12(I)
	END FORALL

	FORALL (I=1:K%NZ(5))
		FAT_K(K%R13(I), 2*DOF+K%C13(I)) = K%K13(I)
		FAT_K(2*DOF+K%C13(I), K%R13(I)) = K%K13(I)
	END FORALL

	FORALL (I=1:K%NZ(6))
		FAT_K(DOF+K%R23(I), 2*DOF+K%C23(I)) = K%K23(I)
		FAT_K(2*DOF+K%C23(I), DOF+K%R23(I)) = K%K23(I)
	END FORALL

    DO I=1, DOF
        ITIS=IS_BOUNDARY(NODES(I),OMEGA)
         IF (ITIS.EQV..TRUE.) THEN
            DO J=1, DOF
            	F(J)=F(J)-U_BD(NODES(I))*FAT_K(J,I)
            ENDDO
        ENDIF
    ENDDO

	ALLOCATE(INDICE(COUNT(BD_MASK)))

	KK = 0

	DO I=1, 3*DOF
		WHERE (BD_MASK(:).EQV..TRUE.)
			FAT_K(I,:) = 0.0D0
			FAT_K(:,I) = 0.0D0
		END WHERE
		IF (BD_MASK(I).EQV..TRUE.) THEN
			KK = KK+1
			INDICE(KK) = I
		ENDIF
	ENDDO	
	
	FORALL (I=1:COUNT(BD_MASK)) FAT_K(INDICE(I), INDICE(I)) = 1.0D0
	
!	DO I=1, 3*DOF
!		IF (BD_MASK(I).EQV..TRUE.) THEN
!			DO J=1, 3*DOF
!				FAT_K(I,J) = 0.0D0
!				FAT_K(J,I) = 0.0D0
!			ENDDO
!			FAT_K(I,I) = 1.0D0
!		ENDIF
!	ENDDO

!	DO I=1, 3*DOF
!		DO J=1+(I-1), 3*DOF
!			FAT_K(J,I) = FAT_K(I,J)
!		ENDDO
!	ENDDO
	
	FULL_MASK(:,:) = .FALSE.
	
	WHERE (DABS(FAT_K(:,:)) .GT. EPS)
		FULL_MASK(:,:) = .TRUE.
	END WHERE
	NZ = UBOUND(PACK(FAT_K(:,:), FULL_MASK),1)
!	KK = 0
!	DO I=1, 3*DOF
!		DO J=1, 3*DOF
!			IF (DABS(FAT_K(I,J)) .GT. EPS) THEN
!				KK = KK+1
!			ENDIF
!		ENDDO
!	ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! CONSTRUCT FULL SPARSE VECTOR !!!!!
	KK = 0
	DO I=1, 3*DOF
		DO J=1, 3*DOF
			IF (DABS(FAT_K(I,J)) .GT. EPS) THEN
				KK = KK+1
				FULL_R(KK) = I
				FULL_C(KK) = J
				FULL_K(KK) = FAT_K(I,J)
			ENDIF
		ENDDO
	ENDDO

!	OPEN(1, FILE='k11.dat')
!	OPEN(2, FILE='k12.dat')
!	DO I=1, DOF
!		WRITE(1,*) (FAT_K(I,J), J=1,DOF)
!		WRITE(2,*) (FAT_K(I,DOF+J), J=1,DOF)
!	ENDDO
!	CLOSE(1);CLOSE(2)

!	OPEN(3, FILE='fat_k.dat')
!	DO I=1, 3*DOF
!		WRITE(3,*) (FAT_K(I,J), J=1, 3*DOF)
!	ENDDO
!	CLOSE(3)
	
	DEALLOCATE(FAT_K, BD_MASK, FULL_MASK, INDICE)
	
	IF (NZ .NE. KK) THEN
		PRINT*, "Oooops!!!!! You got error in BOUNDARY.f90"
		PRINT*, "Total number of nonzeros FAT_K isn't equal to FULL_k"
		STOP
	END IF
	
END SUBROUTINE IMPOSEBD

END MODULE BOUNDARY
