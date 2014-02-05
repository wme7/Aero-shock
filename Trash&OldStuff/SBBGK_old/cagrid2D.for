      program NACA
c
      PARAMETER (MAXWK=300000)
	implicit real*8 (a-h, o-z)
c
      INAME  = 70
	IGRID  = 71
      OPEN(INAME,FILE='B70.TEC' ,STATUS='UNKNOWN')
	OPEN(IGRID,FILE='B11.GRD' ,STATUS='UNKNOWN')
C
      PI=4.*ATAN(1.)
C
      NI=91
	NJ=91
	APHA=20.*PI/180.
C
      CALL GRIDIN(NI,NJ,APHA)

	STOP
	END
C
      SUBROUTINE GRIDIN(NI,NJ,APHA)

	implicit real*8 (a-h, o-z)
	DIMENSION X(NI,NJ),Y(NI,NJ)
	DIMENSION D(NI,NJ),XM(NI,NJ),P(NI,NJ),T(NI,NJ),U(NI,NJ),V(NI,NJ)
!	DIMENSION X1(NI,NJ),Y1(NI,NJ),U1(NI,NJ),V1(NI,NJ)
C
      WRITE(71,*) NI,NJ
C
      READ(70,*)    (( X(I,J),I=1,NI),J=1,NJ),
     1              (( Y(I,J),I=1,NI),J=1,NJ),
     1              (( D(I,J),I=1,NI),J=1,NJ),
     1              ((XM(I,J),I=1,NI),J=1,NJ),
     1              (( P(I,J),I=1,NI),J=1,NJ),
     1              (( T(I,J),I=1,NI),J=1,NJ),
     1              (( U(I,J),I=1,NI),J=1,NJ),
     1              (( V(I,J),I=1,NI),J=1,NJ)
C
      write(71,*)   ((X(I,J),I=1,NI),J=1,NJ),
     1              ((Y(I,J),I=1,NI),J=1,NJ)
!     1              ((D(I,J),I=1,NI),J=1,NJ),
!     1              ((XM(I,J),I=1,NI),J=1,NJ),
!     1              ((P(I,J),I=1,NI),J=1,NJ),
!     1              ((T(I,J),I=1,NI),J=1,NJ),
!     1              ((U(I,J),I=1,NI),J=1,NJ),
!     1              ((V(I,J),I=1,NI),J=1,NJ)
      write(71,*)" END OF GRID"
c
      close(71)
C
      close(70)
C
      RETURN
      END