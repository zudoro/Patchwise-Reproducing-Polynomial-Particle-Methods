! NAME : HYUNJU KIM
! DATE : 11 / 8 / 2010

MODULE GLBVAR

    USE NEWTYPE

    IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Global Parameter

	REAL(8), PARAMETER :: PI=3.141592653589793238462643d0
	INTEGER, PARAMETER :: GORDER=2
	INTEGER, PARAMETER :: IORDER=6
	REAL*8, PARAMETER :: DELTA=0.010D0
	INTEGER, PARAMETER :: PROBLEM=0

	INTEGER, PARAMETER :: NUMLAYERS=1
    TYPE(INT2D), DIMENSION(1:NUMLAYERS), PARAMETER :: DIVS= (/ INT2D(2,2) /)
    INTEGER, PARAMETER :: NUMPATCHES=DIVS(1)%A*DIVS(1)%B
    INTEGER, PARAMETER :: NUMINTSUBPATCH = 9

	REAL(8), PARAMETER:: EPS=5.0D-15

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Global Variables

    INTEGER, PARAMETER :: DOF=NUMPATCHES*(IORDER+1)**2
    INTEGER, PARAMETER :: NGAUSSPT=IORDER+GORDER+2

    TYPE(RECPATCH) :: OMEGA, BGMESH(NUMPATCHES)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Material Coefficients

	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER1=(/ 1, 2 /)
	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER2=(/ 2, 1 /)
												 
	CHARACTER(LEN=3), PARAMETER :: PROPERTY='ISO' ! ISO - ISOTROPIC,  ORT - ORTHOTROPIC
	CHARACTER(LEN=4), PARAMETER :: BCTYPE = 'SCSC' ! CLAM - CLAMPED(FIXED), HARD - HARD SIMPLY SUPPORTED
											! SOFT - SOFT SIMPLY SUPPORTED, FREE - JUST FREE

	REAL*8, PARAMETER :: L=1.0D0	! LENGTH OF THE PLATE											
	REAL*8, PARAMETER :: h=0.010D0	! THICKNESS OF THE PLATE
	REAL*8, PARAMETER :: RHO=1.0D0  ! MATERIAL DENSITY
	REAL*8, PARAMETER :: HOOKS_K=0.0D0	! FOUNDATION MODULUS
	REAL*8, PARAMETER :: LOAD_Q=1.0D0	! UNIFORMLY DISTRBUTED LOAD

	! ISOTROPIC : E1=2=E, V12=V21=V, G12=G13=G23=G = E/2(1+V)
		
	REAL*8, DIMENSION(1:2), PARAMETER :: E = (/1.0D0, 1.0D0/) 	! YOUNG'S MODULUS {E1, E2}
	REAL*8, DIMENSION(1:2), PARAMETER :: V = (/0.30D0, 0.30D0/)		! POISSON RATIO {V12, V21}
	
	REAL*8, DIMENSION(1:3), PARAMETER :: G = (E(1)/(2.0D0*(1+V(1))))*(/1.0D0, 1.0D0, 1.0D0/) ! SHEAR MODULUS {G12, G13, G23}

	REAL*8, DIMENSION(1:3), PARAMETER :: Q1 = (/E(1)/(1.0D0-V(1)*V(2)), V(1)*E(2)/(1.0D0-V(1)*V(2)), 0.0D0/)
	REAL*8, DIMENSION(1:3), PARAMETER :: Q2 = (/V(1)*E(2)/(1.0D0-V(1)*V(2)), E(2)/(1.0D0-V(1)*V(2)), 0.0D0/)
	REAL*8, DIMENSION(1:3), PARAMETER :: Q3 = (/0.0D0, 0.0D0, G(1)/)
	REAL*8, PARAMETER :: Q4=G(3)
	REAL*8, PARAMETER :: Q5=G(2)

!	REAL*8, PARAMETER :: KS=5.0D0/6.0D0		! SHEAR CORRECTION FACTOR
!	REAL*8, PARAMETER :: KS=0.8330D0		! SHEAR CORRECTION FACTOR		
!	REAL*8, PARAMETER :: KS=0.86010D0		! SHEAR CORRECTION FACTOR
	REAL*8, PARAMETER :: KS=0.8220D0		! SHEAR CORRECTION FACTOR
	! MORE EXACT EXPRESSION FOR THE CASE OF ISOTROPIC : KS = 20(1+V)/(24+25*V+V**2)
	
	! THE PLANE STRESS-REDUCED STIFFNESS COEFFICIENTS
	REAL*8, DIMENSION(1:3,1:3), PARAMETER :: STRESS_Q=RESHAPE((/ Q1, Q2, Q3 /), (/ 3, 3 /), ORDER=ORDER2)
	REAL*8, DIMENSION(1:2,1:2), PARAMETER :: SHEAR_Q=RESHAPE((/ Q4, 0.0D0, 0.0D0, Q5 /), (/ 2, 2 /), ORDER=ORDER2)
	! EXTENSIONAL STIFFNESS COEFFICIENTS
	REAL*8, DIMENSION(1:3,1:3), PARAMETER :: STRESS_A=h*STRESS_Q
	REAL*8, DIMENSION(1:2,1:2), PARAMETER :: SHEAR_A=h*SHEAR_Q
	! BENDING STIFFNESS COEFFICIENTS
	REAL*8, DIMENSION(1:3,1:3), PARAMETER :: STRESS_D=(2.0D0/3)*((0.50D0*h)**3)*STRESS_Q
	REAL*8, DIMENSION(1:2,1:2), PARAMETER :: SHEAR_D=(2.0D0/3)*((0.50D0*h)**3)*SHEAR_Q
	! HIGHER ORDER STIFFNESS COEFFICIENTS
	REAL*8, DIMENSION(1:3,1:3), PARAMETER :: STRESS_F=(2.0D0/5)*((0.50D0*h)**5)*STRESS_Q
	REAL*8, DIMENSION(1:3,1:3), PARAMETER :: STRESS_H=(2.0D0/7)*((0.50D0*h)**7)*STRESS_Q
	REAL*8, DIMENSION(1:2,1:2), PARAMETER :: SHEAR_F=(2.0D0/5)*((0.50D0*h)**5)*SHEAR_Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	ROW, COLUMN AND COMPONENTS OF SUB MATRICES

    
END MODULE GLBVAR
