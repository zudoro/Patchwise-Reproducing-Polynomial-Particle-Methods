MODULE BICGSTAB

!Assemble stiffness matrices and load vectors, and then compute
!solution vectors using BiCGSTAB method.

CONTAINS

SUBROUTINE SOLVER(X, DNR, IT, JC, C, B, N, NZ)

! OUTPUT --- X(N) : SOLUTION VECTOR
! INPUT    --- N : THE ORDER OF STIFFNESS MATRIX
! INPUT    --- NZ : THE NUMBER OF NON-ZEROS STORED IN STIFFNESS MATRIX

	  IMPLICIT NONE

    INTEGER, PARAMETER :: ITMAX=100
    REAL*8, PARAMETER :: TOL=1.0D-10
    INTEGER, INTENT(IN) :: N, NZ, IT(NZ), JC(NZ)
    REAL*8, INTENT(IN) :: C(NZ), B(N)
    REAL*8, INTENT(OUT) :: X(N), DNR
    REAL*8 :: A(NZ)
    INTEGER ::  IC(N+1)
    INTEGER :: IA(NZ), JA(N+1), JU(N)
	INTEGER :: I, J, K, NC

! A(*)  : Global stiffness matrix
! IA(*) : The row index of nonzero elements of A
! JA(I) : Pointer to the first element of the i-th column of A
! TOL   : Tolerance for stopping criterion of BiCGSTAB
! ITMAX : Maximum number of iteration to be allowed in the BiCGSTAB method
! B(*)  : Global load vector

!      OPEN(7,FILE='row.txt',STATUS='old',ACCESS='SEQUENTIAL', FORM='FORMATTED', ACTION='READWRITE')
!      OPEN(8,FILE='col.txt',STATUS='old',ACCESS='SEQUENTIAL', FORM='FORMATTED', ACTION='READWRITE')
!      OPEN(9,FILE='entry.txt',STATUS='old',ACCESS='SEQUENTIAL', FORM='FORMATTED', ACTION='READWRITE')
!      OPEN(10,FILE='load.txt',STATUS='old',ACCESS='SEQUENTIAL', FORM='FORMATTED', ACTION='READWRITE')

!      READ(7,*) (IT(J),J=1,NZ) ; READ(8,*) (JC(J),J=1,NZ) ; READ(9,*) (C(J),J=1,NZ) ; READ(10,*) (B(J),J=1,N)
!	  CLOSE(7, STATUS='DELETE') ; CLOSE(8, STATUS='DELETE') ; CLOSE(9, STATUS='DELETE') ; CLOSE(10, STATUS='DELETE')

	PRINT*, ""
	PRINT*, "SOLVE THE SYSTEM"
	
    K=1
	NC=1
	IC(K)=NC
	DO J=1,NZ-1
        NC=NC+1
        IF(IT(J).NE.IT(J+1)) THEN
           K=K+1
           IC(K)=NC
        END IF
	END DO
	IC(K+1)=NC+1


! Transform the Compressed row Storage Format into the Compressed column Storage Format

	NC=1
	DO J=1,N
	JA(J)=NC
        DO I=1,N
            DO K=IC(I),IC(I+1)-1
                IF(JC(K).EQ.J) THEN
                    IA(NC)=I
                    A(NC)=C(K)
                    NC=NC+1
                    GO TO 7
                ELSE IF(JC(K).GT.J) THEN
                    GO TO 7
                END IF
            END DO
7          CONTINUE
        END DO
	END DO
	JA(N+1)=NC

	CALL BCGSTAB(N, NC-1, A, IA, JA, JU, B, X, DNR, TOL, ITMAX)

END SUBROUTINE SOLVER

SUBROUTINE BCGSTAB(N, NZ, A, IA, JA, JU, B, X, DNR, TOL, ITMAX)

! Right-preconditioned BiCGSTAB Method for solving A*x=b.

      IMPLICIT REAL*8( A-H,O-Z)

      INTEGER, INTENT(IN) :: N, NZ, ITMAX
      INTEGER, PARAMETER :: ISYM=0
      INTEGER, INTENT(IN) ::  IA(NZ), JA(N+1)
      INTEGER, INTENT(INOUT) :: JU(N)
      REAL*8, INTENT(IN) ::    A(NZ), B(N), TOL
      REAL*8, INTENT(OUT) :: X(N), DNR
      REAL*8 :: M(NZ), P(N), AP(N), PM(N), RM(N), R(N), AR(N), HR(N)
      REAL*8 :: ALPHA, BETA, W1, TEMP1, TEMP2, BNORM, TEMP

!      OPEN(6,FILE='estimation_info.txt')
!      OPEN(11,FILE='out.txt')

! COPY A TO M

	DO I=1,NZ
	 M(I)=A(I)
	END DO

!      WRITE(6,*) 'The number of nonzero elements of A is ', NZ

! Generate a preconditioner M = LU such that A = LU - R using incomplete
! LU factorization, where L is a unit lower triangular matrix and
! U is an upper triangular matrix. Notice that both L and U are stored on
! array M.

      DO 480 I=1,N
      DO 470 J=JA(I),JA(I+1)-1
        IF(IA(J).EQ.I) THEN
          JU(I)=J
          GO TO 480
        END IF
470   CONTINUE
480   CONTINUE

      DO 460 I=1,N-1
       DO 450 J=JU(I)+1,JA(I+1)-1
         M(J)=M(J)/M(JU(I))
 450   CONTINUE

       DO 440 IK=I+1, N
        DO 435 J=JA(IK),JU(IK)
          IF(IA(J).EQ.I) THEN
            DO 425 L=J+1,JA(IK+1)-1
             DO 420 LL=JU(I)+1,JA(I+1)-1
               IF(IA(L).EQ.IA(LL)) THEN
                 M(L)=M(L)-M(LL)*M(J)
                 GO TO 425
               END IF
420         CONTINUE
425        CONTINUE
           GO TO 440
         END IF
435     CONTINUE
440    CONTINUE

460   CONTINUE

      BNORM=DSQRT(DDOT(N,B,B))

!   CHOOSE INITIAL GUESS FOR VECTOR X

      DO 5 I=1,N
        X(I)=0.D0
5     CONTINUE

      CALL DSMV(N,X,R,NZ,IA,JA,A,ISYM)
      DO 20 I=1,N
        R(I)=B(I) - R(I)
        HR(I)=R(I)
        P(I)=R(I)
   20 CONTINUE

      TEMP1=DDOT(N,HR,R)

      DO 400 I=1,ITMAX
        CALL DSOLVE(N,PM,P,NZ,IA,JA,JU,M)
        CALL DSMV(N,PM,AP,NZ,IA,JA,A,ISYM)
        ALPHA=TEMP1/DDOT(N,HR,AP)
        DO 50 J=1,N
         R(J)=R(J) - ALPHA * AP(J)
   50   CONTINUE
        CALL DSOLVE(N,RM,R,NZ,IA,JA,JU,M)
        CALL DSMV(N,RM,AR,NZ,IA,JA,A,ISYM)
        W1=DDOT(N,AR,R)/DDOT(N,AR,AR)
        DO 60 J=1,N
          X(J)=X(J) + ALPHA * PM(J) + W1 * RM(J)
          R(J)=R(J) - W1*AR(J)
   60   CONTINUE
        DNR=DSQRT(DDOT(N,R,R))
        IF(DNR/BNORM .LT. TOL) GO TO 210
        TEMP2=DDOT(N,HR,R)
        BETA=(TEMP2 / TEMP1) * (ALPHA/W1)
        DO 80 J=1,N
          P(J)=R(J) + BETA * (P(J) - W1 * AP(J))
   80  CONTINUE
       TEMP1=TEMP2
  400 CONTINUE

  210 CONTINUE
!      WRITE(6,99) I, DNR, BNORM
!99    FORMAT(1X,'Iter=',I4, '  Rnorm=',D22.14, '  Bnorm=',D22.14)

!      WRITE(11,98) (X(J),J=1,N)
!98    FORMAT(1(1X, D22.14))
!      CLOSE(6)
!      CLOSE(11)
END SUBROUTINE BCGSTAB

! SET OF SUBPROGRAMS ---------------------------------------.

REAL*8 FUNCTION DDOT(N,DX,DY)
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: N
!     .. Array Arguments ..
      REAL*8, INTENT(IN) :: DX(N),DY(N)
      REAL*8 :: DTEMP
      INTEGER :: I,IX,IY,M,MP1, INCX, INCY
      INCX=1 ; INCY=1 ; DDOT = 0.0d0 ; DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20

      IX = 1 ; IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        code for both increments equal to 1
!
!        clean-up loop
!
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
END FUNCTION DDOT

SUBROUTINE DSOLVE(N,X,Y,NELT,IA,JA,JU,M)

! This solves LUx = y, where L is a unit lower triangular matrix and U
! is an upper triangular matrix.

    IMPLICIT REAL*8( A-H,O-Z)

    INTEGER, INTENT(IN) :: N, NELT
    REAL*8, INTENT(OUT) :: X(N)
    REAL*8, INTENT(IN) ::   M(NELT),Y(N)
    INTEGER, INTENT(IN) ::  IA(NELT),JA(N+1),JU(N)

    DO 100 I=1,N
    X(I)=Y(I)
100 CONTINUE

! Solve Lw=x

      DO 90 I=1,N-1
      DO 90 J=JU(I)+1,JA(I+1)-1
        X(IA(J))=X(IA(J))-M(J)*X(I)
90    CONTINUE

! Solve Uw=x

      DO 80 I=N,2,-1
       X(I)=X(I)/M(JU(I))
       DO 70 J=JA(I),JU(I)-1
        X(IA(J))=X(IA(J))-M(J)*X(I)
70     CONTINUE
80    CONTINUE
      X(1)=X(1)/M(JU(1))

END SUBROUTINE DSOLVE


SUBROUTINE DSMV(N,X,Y,NZ,IA,JA,A,ISYM)

!***PURPOSE  Sparse Matrix Vector Product
!                Y = A*X.

! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! X      :IN       Double Precision X(N).
!         The vector that should be multiplied by the matrix.
! Y      :OUT      Double Precision Y(N).
!         The product of the matrix and the vector.
! NZ     :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NZ).
! JA     :IN       Integer JA(N+1).
! A      :IN       Integer A(NZ).
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all nonzero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle part of the matrix A is stored.


      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER, INTENT(IN) :: N, NZ, ISYM
      INTEGER, INTENT(IN) :: IA(NZ),JA(N+1)
      REAL*8, INTENT(IN) :: A(NZ),X(N)
      REAL*8, INTENT(OUT) :: Y(N)

!         Zero out the result vector.

      DO 10 I=1,N
         Y(I)=0.D0
 10   CONTINUE

!         Multiply by A.

      DO 30 ICOL=1,N
         IBGN=JA(ICOL)
         IEND=JA(ICOL+1)-1
         DO 20 I=IBGN,IEND
            Y(IA(I))=Y(IA(I))+A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE

      IF(ISYM.EQ.1) THEN

!         The matrix is symmetric.  Need to get the other half in...
!         This loops assumes that the diagonal is the first entry in
!         each column.

         DO 50 IROW=1,N
            JBGN=JA(IROW)+1
            JEND=JA(IROW+1)-1
            IF(JBGN.GT.JEND) GOTO 50
            DO 40 J=JBGN,JEND
               Y(IROW)=Y(IROW)+A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ENDIF

END SUBROUTINE DSMV

END MODULE BICGSTAB
