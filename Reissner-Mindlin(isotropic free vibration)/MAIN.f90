PROGRAM MAIN

     USE NEWTYPE
     USE GLBVAR
     USE DOMAIN
     USE PLOT
	 	 USE PACKAGE
     USE ASSEMBLE
     USE BOUNDARY
     USE EIGEN
     USE LUDECOMPOSITION
     USE ERRORESTIMATE

     IMPLICIT NONE

	TYPE(SUBMAT) :: K
    REAL*8 :: M(INT(0.5*DOF*(DOF+1)))
    REAL*8, ALLOCATABLE :: THIN_K(:,:), THIN_M(:,:), W(:), Z(:,:)
    INTEGER :: INDX(3*DOF), I, J, REDUCEDOF, IERR
    LOGICAL :: BD_MASK(3*DOF)
    REAL :: CPTS, CPTE
    
	 	 CALL CPU_TIME(CPTS)
!    Get domain information and plot the domain
     CALL GETDOMAIN(OMEGA)
!    CALL PLOTDOMAIN()

!    Get background mesh and plot the mesh
     CALL GETBGMESH(OMEGA,BGMESH)
!     CALL PLOTBGMESH()

!	 Check pu and special (interpolation) functions
!	 CALL CHECKBASIS()

!    Assemble stiffness matrix and mass matrix
     CALL STIF(K)
     CALL MASS(M)
	 
!    Find boundary nodes
	 CALL INFONODE(REDUCEDOF, BD_MASK)
	 
	 ALLOCATE(THIN_K(REDUCEDOF,REDUCEDOF), THIN_M(REDUCEDOF,REDUCEDOF), W(REDUCEDOF), Z(REDUCEDOF,REDUCEDOF))
	 
!	 Reduce matrices [K] and [M]
	 CALL REDUCEMATRIX(THIN_K, THIN_M, K, M, REDUCEDOF, BD_MASK)
!	 CALL TEST08(REDUCEDOF,THIN_K,THIN_M)
 	 OPEN(1, FILE='thin_k.dat')
 		DO I=1, REDUCEDOF
			WRITE(1,*) (THIN_K(I,J), J=1,REDUCEDOF)
		ENDDO
	 CLOSE(1)
	 
	 OPEN(1, FILE='thin_m.dat')
		DO I=1, REDUCEDOF
			WRITE(1,*) (THIN_M(I,J), J=1,REDUCEDOF)
		ENDDO
	 CLOSE(1)
! 
! !	 Find the smallest eigenvale (Fundamental frequency)
! !	 CALL PWMTD(THIN_K, THIN_M, W, IND)
! 	 CALL RSG ( REDUCEDOF, THIN_K, THIN_M, W, 1, Z, IERR )
! 	 PRINT*, ""
! 	 PRINT*, IERR
! 	 PRINT*, ''
! !	 PRINT*, "FUNDAMENTAL FREQUENCY : ", W*L*DSQRT(RHO/ G(1))
! 	 PRINT*, "THE FIRST SIX FREQUENCY PARAMETERS: "
! 	 DO I=1, 3
! 		PRINT*, W(I)*L*DSQRT(RHO/ G(1))
! 	 ENDDO
! 	 PRINT*, ""
! 
! 	 CALL CPU_TIME(CPTE)
! 	 PRINT*, ""
! 	 PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
! 	 PRINT*, ""
! 
! 	 OPEN(1, FILE = 'result.dat')
! 	 	WRITE(1,*) "-----------------------------------"
! 	 	WRITE(1,*) "PARAMETERS OF CONSIDERD MATERIAL : "
! 	 	IF (PROPERTY.EQ.'ISO') THEN
! 	 		WRITE(1,*) "THIS MATERIAL IS ISOTROPIC SQUARE PLATE."
! 	 	ELSEIF (PROPERTY.EQ.'ORT') THEN
! 	 		WRITE(1,*) "THIS MATERIAL IS ORTHOTROPIC SQUARE PLATE."
! 	 	ENDIF
! 		WRITE(1,*) "TYPE OF BOUNDARY CONDITION : ", BCTYPE
! 		WRITE(1,*) "SIDE OF THE PLATE : ", L
! 		WRITE(1,*) "THICKNESS OF THE PLATE : ", h
! 		WRITE(1,*) "DENSITY OF THE MATERIAL : ", RHO
! 		WRITE(1,*) "REACTION FORCE BY HOOK'S LAW : ", HOOKS_K
! 		WRITE(1,*) ""
! 		WRITE(1,*) "RATIO OF THICKNESS / SIDE : ", h / L
! 		WRITE(1,*) "YOUNG'S MODULUS E : ", E(1)
! 		WRITE(1,*) "SHEAR MODULUS(CONSTANT OF RIGIDITY) : ", G(1)
! 		WRITE(1,*) "POISSON RATIO : ", V(1)
! 		WRITE(1,*) "SHEAR CORRECTION FACTOR : ", KS
! 	 	WRITE(1,*) "-----------------------------------"
! 	 	WRITE(1,*) "PARAMETERS OF PUFEM : "
! 	 	WRITE(1,*) "PU ORDER : ", GORDER
! 	 	WRITE(1,*) "RPP ORDER : ", IORDER
! 	 	WRITE(1,*) "DELTA : ", DELTA
! 	 	WRITE(1,*) "NUMBER OF PATCHES : ", NUMPATCHES
! 	 	WRITE(1,*) "DEGREE OF FREEDOM : ", DOF
! 	 	WRITE(1,*) "NUMBER OF GAUSS POINTS ON EACH PATCH : ", NGAUSSPT
! 	 	WRITE(1,*) "-----------------------------------"
! 		WRITE(1,*) "THE FIRST SIX FREQUENCY PARAMETERS: "
! 		DO I=1, 6
! 			WRITE(1,*) W(I)*L*DSQRT(RHO/ G(1))
! 		ENDDO
! 		WRITE(1,*) ""
! 		WRITE(1,*) "ELAPSED CPU TIME : ", CPTE - CPTS
! 	 CLOSE(1)
! 	 DEALLOCATE(THIN_K,THIN_M, W, Z)
	 
END PROGRAM MAIN