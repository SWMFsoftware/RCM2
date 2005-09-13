!
   MODULE Hilmer_Voigt_subs_mod

     ! Added for SWMF
     use ModIoUnit, ONLY : io_unit_new

   IMPLICIT NONE
   INTEGER, PARAMETER :: LATDIM = 62, LTDIM = 51
   REAL :: vm_msm (latdim,ltdim), be_msm (latdim,ltdim), &
           xmin_msm (latdim,ltdim), ymin_msm (latdim,ltdim), &
           zmin_msm (latdim,ltdim)
      CHARACTER (LEN=37), PARAMETER :: prefix = '~/prj/rcm/inputs/bfield/hilmer/bin/bo'
      CHARACTER (LEN=*), PARAMETER :: format_prefix = "(A37,2I1,2I2.2,I1,'.bin')"
!     CHARACTER (LEN=34), PARAMETER :: &
!        prefix = 'c:\rcm\inputs\bfield\hilmer\bin\bo'
!     CHARACTER (LEN=*), PARAMETER :: format_prefix = "(A34, 2I1,2I2.2,I1,'.bin')"
!
!
!
!  BTRACE |
!         | <---- GETMAT | <---- FNDBRK
!         |              | <---- MEXSIT
!         |              | <---- LOADBM
!         |              | <---- ZEROBM
!         |
!         | <---- RMVBSH
!
!
   CONTAINS

      SUBROUTINE Getmat (iwant, latdim, ltdim, bfpar, work, bflim, bfexst)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: iwant, latdim, ltdim
      REAL,    INTENT (IN) :: bfpar (5)
      REAL,    INTENT (OUT):: bflim (2, 5)
      LOGICAL, INTENT (OUT) :: bfexst (2, 2, 2, 2, 2)
      REAL,    INTENT (OUT) :: work (latdim,ltdim,2,2,2,2,2)
!
!
!  copyright Rice University, 1996
!
!   THIS SUBROUTINE FINDS THE MAGNETIC FIELD MATRICES NEEDED TO
!   INTERPOLATE ONE MATRIX THAT IS REPRESENTATIVE OF CURRENT
!   CONDITIONS.
!
!   VERSION 1.00                           9/12/89
!           2.00 MSM DELIVERY VERSION      02.02.93
!           2.10                           06.11.93
!                ERROR OUTPUT ROUTED TO UNIT 9.
!                SOME INFORMATIONAL OUTPUT COMMENTED OUT.
!           2.20                           09.06.94
!                STOP STATEMENT REPLACED WITH CALL TO ERRSUB
!
!  PROGRAMMER: BRYAN BALES
!
!  PURPOSE:
!     FOR A GIVEN SET OF MAGNETIC FIELD INPUT PARAMETERS, FINDS THE
!     MAGNETIC FIELD MATRICES NEEDED FOR A RUN
!
!  INPUTS:
!    IWANT      DIFFERENT MAGNETIC FIELD PARAMETERS
!               =1 XMIN (GSM X LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =2 YMIN (GSM Y LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =3 ZMIN (GSM Z LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =4 BMIN (EQUATORIAL MAGNETIC FIELD STRENGTH (NT))
!               =5 VM (FLUX TUBE VOLUME**-(2/3) (1/WEBER))
!    LATDIM     NUMBER OF LATITUDINAL GRID LINES
!    LTDIM      NUMBER OF LOCAL TIME (LONGITUNDINAL) GRID LINES
!    BFPAR      ARRAY CONTAINING GEOPHYSICAL PARAMETERS TO USE TO FINE
!               MAGNETIC FIELD MATRICES
!               1 = STANDOFF DISTANCE (RE)
!               2 = TILT ANGLE (DEGREES)
!               3 = EQUATORWARD EDGE OF AURORAL OVAL
!               4 = DST (NT)
!               5 = COLLAPSE PARAMETER
!
!  OUTPUTS:
!    WORK       INTERNAL WORKING ARRAY
!    BFLIM      INDICES OF MAGNETIC FIELD MODELS READ
!    BFEXST     ARRAY OF LOCICAL VARIABLES OF WHICH MAGNETIC FIELD MODEL
!               MATRICES EXIST
!
      INTEGER, PARAMETER :: imstnd=6, imtilt=5, imined=18, imdst=10, imstch=2
      REAL, PARAMETER :: stndpr (imstnd) = (/ 4.0, 6.0, 8.0, 10.0, 12.0, 14.0/)
      REAL, PARAMETER :: tiltpr (imtilt) = (/ -35.0, -17.5, 0.0, 17.5, 35.0 /)
      REAL, PARAMETER :: finedp (imined) = &
                         (/ 34.38, 43.93, 49.47, 53.20, 55.89, 57.96, 59.60,&
                            60.95, 62.06, 63.02, 63.99, 65.01, &
                            65.91, 66.72, 67.45, 68.12, 68.73, 69.30 /)
      REAL, PARAMETER :: dstpr (imdst) = &
                         (/ -600.0, -500.0, -400.0, -300.0, -200.0, -150.0, &
                            -100.0,  -50.0,    0.0,   50.0 /)
      REAL, PARAMETER :: stchpr (imstch)  = (/ 0.0, 1.0 /)
      INTEGER :: bfndx (2,5), i, j, k, l, m, nummod
!
!
!
      CALL Fndbrk (bfpar(1), stndpr, imstnd, bfndx(1,1), bfndx(2,1))
      CALL Fndbrk (bfpar(2), tiltpr, imtilt, bfndx(1,2), bfndx(2,2))
      CALL Fndbrk (bfpar(3), finedp, imined, bfndx(1,3), bfndx(2,3))
      CALL Fndbrk (bfpar(4), dstpr,  imdst,  bfndx(1,4), bfndx(2,4))
      CALL Fndbrk (bfpar(5), stchpr, imstch, bfndx(1,5), bfndx(2,5))
!
      IF (iwant == 1) then
         WRITE (*,*) '***','stoff',bfpar(1),bfndx(1,1),bfndx(2,1)
         WRITE (*,*) 'tilt ',bfpar(2),bfndx(1,2),bfndx(2,2)
         WRITE (*,*) 'fmeb ',bfpar(3),bfndx(1,3),bfndx(2,3)
         WRITE (*,*) 'dst  ',bfpar(4),bfndx(1,4),bfndx(2,4)
         WRITE (*,*) 'clpse',bfpar(5),bfndx(1,5),bfndx(2,5)
      END IF
!
      DO i = 1, 2
         bflim(i,1) = stndpr (bfndx(i,1))
         bflim(i,2) = tiltpr (bfndx(i,2))
         bflim(i,3) = finedp (bfndx(i,3))
         bflim(i,4) = dstpr  (bfndx(i,4))
         bflim(i,5) = stchpr (bfndx(i,5))
      END DO
!
 6    nummod = 0
      DO I = 1, 2
      DO J = 1, 2
      DO K = 1, 2
      DO L = 1, 2
      DO M = 1, 2
         IF ( Mexist ( &
              bfndx(i,1),bfndx(j,2),bfndx(k,3), bfndx(l,4),bfndx(m,5))) THEN
              CALL Loadbm (iwant, latdim, ltdim, i, j, k, l, m, bfndx, work)
              bfexst(i,j,k,l,m) = .TRUE.
              nummod = nummod+1
         ELSE
              CALL Zerobm (latdim, ltdim, i, j, k, l, m, work)
              bfexst (i,j,k,l,m) = .FALSE.
         END IF
         IF (iwant == 1) THEN
            WRITE (*,1001) &
               bfndx(j,2), bfndx(i,1), bfndx(l,4), bfndx(k,3),(3-bfndx(m,5)),&
               bfexst(i,j,k,l,m), &
               bflim(i,1), bflim(j,2), bflim(k,3), bflim(l,4), bflim(m,5)
         END IF
      END DO
      END DO
      END DO
      END DO
      END DO
!
      IF (nummod == 0) THEN
         IF (bfndx(1,4) >= 2) THEN
            bfndx (1,4) = bfndx(1,4)-1
            bfndx (2,4) = bfndx(2,4)-1
            WRITE (*,'(A)',ADVANCE='NO') 'NO BF-MODELS. ADJUST DST. HIT ENTER..'
            READ (*,*)
            GO TO 6
         ELSE
           STOP 'SERIOUS PROBLEM FINDING CORRECT B-FIELD, SEE GETMAT CODE'
         END IF
      END IF
!
!
 1001 FORMAT (TR1, 'BF'        , 2I1, 2I2.2, I1, &
              TR1,  L1, &
              TR3, 'STANDOFF= ', F4.0, &
              TR1, 'TILT= ',     F4.0, &
              TR1, 'EQEDGE= ',   F5.2, &
              TR1, 'DST= ',      F5.0, &
              TR1, 'COLLAPSE= ', F3.0    )
      RETURN
      END SUBROUTINE Getmat
!
!
      SUBROUTINE Rmvbsh (work, latdim, ltdim, iwant)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: latdim, ltdim, iwant
      REAL,    INTENT (IN OUT) :: work (latdim,ltdim,2,2,2,2,2)
!
!
!  copyright Rice University, 1996
!
!  VERSION 1.00
!          2.00                                  02.03.93
!                                         MSM DELIVERY VERSION
!
!  PURPOSE:
!    FIX TO THE B-FIELD MATRICES TO REMOVE THE LARGE DATA
!    VALUES SIGNIFYING THAT THE FIELD LINES HAD 'GONE INTO THE BUSHES'.
!    IT LOOKS FOR THE CLOSEST VALUE IN J AND LINEARLY INTERPOLATES FOR
!    THE I VALUES.
!
!  INPUTS:
!    WORK       WORKING MAGNETIC FIELD MATRIX
!    LATDIM     NUMBER OF LATITUDINAL GRID LINES
!    LTDIM      NUMBER OF LOCAL TIME (LONGITUDINAL) GRID LINES
!               =1 XMIN (GSM X LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =2 YMIN (GSM Y LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =3 ZMIN (GSM Z LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =4 BMIN (EQUATORIAL MAGNETIC FIELD STRENGTH (NT))
!               =5 VM (FLUX TUBE VOLUME**-(2/3) (1/WEBER))
!
!  OUTPUTS:
!    WORK       REVISED WORKING MAGNETIC FIELD MATRIX
!
!
      REAL, PARAMETER :: bfbush = 1.E20
      INTEGER :: n1, n2, n3, n4, n5, ido, ichk, j
      REAL    :: coef1, coef2, temp, strtln
      DO N1=1,2
      DO N2=1,2
      DO N3=1,2
      DO N4=1,2
      DO N5=1,2
         DO J = 1, ltdim
             ido = 1
             IF ( ABS (work(ido,j,n1,n2,n3,n4,n5)) < bfbush) CYCLE
             DO
                ido=ido+1
                IF (ABS(work(ido,j,n1,n2,n3,n4,n5))<= bfbush .AND.&
                  ABS(work(ido+1,j,n1,n2,n3,n4,n5)) <= bfbush) EXIT
             END DO
             coef1 = work (ido,  j,n1,n2,n3,n4,n5)
             coef2 = work (ido+1,j,n1,n2,n3,n4,n5)
             strtln = coef2-coef1
             DO ichk = ido-1, 1, -1
                temp = work(ichk+1,j,n1,n2,n3,n4,n5)-strtln
                IF (iwant /= 5) THEN
                  work (ichk,j,n1,n2,n3,n4,n5) = temp
                ELSE
                  work (ichk,j,n1,n2,n3,n4,n5) = MAX (0.1, 0.5*coef1, temp)
                END IF
             END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END SUBROUTINE Rmvbsh
!
!
!
!
      SUBROUTINE Fndbrk (parval, pvals, ipdim, indmin, indmax)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: ipdim
      REAL,    INTENT (IN) :: parval, pvals (ipdim)
      INTEGER, INTENT (OUT):: indmin, indmax
!
!  copyright Rice University, 1996
!
!  PURPOSE:
!   THIS SUBROUTINE FINDS THE INDICES OF THE VALUES IN ARRAY
!   PVALS THAT BRACKET PARVAL.
!
!   VERSION 1.00                               9/12/89
!           2.00 MSM DELIVERY VERSION          02.02.93
!           2.10                               08.02.96
!                MINOR INTERNAL DOCUMENTATION ADDITIONS
!
!   PROGRAMMER: BRYAN BALES
!
!   INPUTS:
!     PARVAL       MAGNETIC FIELD GEOPHYSICAL PARAMETER
!     PVALS        ARRAY OF MAGNETIC FIELD PARAMETERS
!     IPDIM        DIMENSION OF PVALS
!
!   OUTPUTS:
!     MIN          LOWER ARRAY BOUND INDEX
!     MAX          UPPER ARRAY BOUND INDEX
!
!
      indmin = 0
      DO
         indmin = indmin + 1
         indmax = indmin + 1
         IF (parval <= pvals(indmax) .OR. indmax >= ipdim) EXIT
      END DO
      RETURN
      END SUBROUTINE Fndbrk
!
!
!
      SUBROUTINE Loadbm (iwant, latdim, ltdim, i, j, k, l, m, bfndx, work)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: iwant, latdim, ltdim, i, j, k, l, m, bfndx(2,5)
      REAL, INTENT (IN OUT) ::  work(latdim,ltdim,2,2,2,2,2)
!
!
!  copyright Rice University, 1996
!
!  CHANGE MADE 2/5/93 TO USE HILMER MARK 501 B-FIELD MODELS
!    LOCATED IN /usr/local/balestmp
!
!  VERSION 2.00 MSM DELIVERY VERSION              02.02.93
!          2.10                                   06.11.93
!               ERROR OUTPUT ROUTED TO UNIT 9
!          2.20                                   09.12.94
!               STOP STATEMENT REPLACED WITH CALL TO ERRSUB
!               ERR AND IOSTAT VARIABLES ADDED TO IO STATEMENTS
!
!  PURPOSE
!   THIS SUBROUTINE LOADS INDIVIDUAL B-MATRICES FROM THE OFFLINE
!   B-SUPERMATRIX INTO THE WORK B-MATRICES.  A CHECK IS MADE TO
!   VERIFY THAT THE CORRECT MATRIX HAS BEEN RETRIEVED.
!
!  INPUTS:
!    IWANT      DIFFERENT MAGNETIC FIELD PARAMETERS
!               =1 XMIN (GSM X LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =2 YMIN (GSM Y LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =3 ZMIN (GSM Z LOCATION OF WHERE FIELD LINE GOING THROUGH
!               GRID POINT CROSSES THE EQUATORIAL (B-FIELD MINIMUM)
!               PLANE (RE))
!               =4 BMIN (EQUATORIAL MAGNETIC FIELD STRENGTH (NT))
!               =5 VM (FLUX TUBE VOLUME**-(2/3) (1/WEBER))
!    LATDIM     NUMBER OF LATITUDINAL GRID LINES
!    LTDIM      NUMBER OF LOCAL TIME (LONGITUNDINAL) GRID LINES
!    I,J,K,L,M  INDICES OF WORKING MAGNETIC FIELD MATRICES
!    BFNDX      ARRAY CONTAINING THE INDICES OF THE MAGNETIC FIELD MATRICES
!               TO LOAD FOR GEOPHYSICAL CONDITIONS OF STANDOFF DISTANCE,
!               TILT ANGLE, EQUATORWARD EDGE OF THE AURORAL OVAL, DST, AND
!               TAIL COLLAPSE INDEX
!
!  OUTPUTS:
!    WORK       MAGNETIC FIELD MATRIX WORK ARRAY
!
      CHARACTER(LEN=80) :: filnam, form_char
      INTEGER :: LUN, ios, idf1, idf2, idf3, idf4, idf5, idf6, n, ii, jj
      LOGICAL :: Lexist, Lopen
!
!   THE FILE DEFINITIONS FOR THE B-MATRICES IN THE SUPERMATRIX:
!                    'BFvwxxyyz'.
!                         WHERE v IS THE TILT INDEX,
!                               w IS THE STANDOFF INDEX,
!                              xx IS THE DST INDEX,
!                              yy IS THE INNER EDGE INDEX,
!                               z IS THE COLLAPSE INDEX.

      filnam = Get_filnam (bfndx(i,1), bfndx (j,2), bfndx (k,3), &
                           bfndx (l,4), bfndx (m,5) )
      filnam = TRIM (filnam)
      IF (INDEX (filnam, '.dat') /= 0) THEN
!        dealing with formatted bfield files:
         form_char = 'FORMATTED'
         form_char = TRIM (form_char)
      ELSE IF (INDEX (filnam, '.bin') /= 0) THEN
!        dealing with unformatted bfield files:
          form_char = 'UNFORMATTED'
          form_char = TRIM (form_char)
      ELSE
          STOP 'UNABLE TO DETERMINE FORMAT OF BFIELD FILES'
      END IF
!
!
!     FIND AVAILABLE LOGICAL UNIT NUMBER FOR I/O:
!
      ! Added for SWMF
      LUN = io_unit_new()
      ! Commented out for SWMF
!      LUN = 99
!      DO LUN = 99, 11, -1
!         INQUIRE (UNIT = LUN, EXIST = Lexist, OPENED = Lopen )
!         IF (.NOT.Lexist) CYCLE
!         IF (Lopen)  CYCLE
!         EXIT
!      END DO
!      IF (LUN < 11) STOP 'no units available to open a file in LOADBM'
!
!
      OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//FILNAM, FORM = form_char, STATUS = 'OLD',&
            ACCESS = 'SEQUENTIAL', IOSTAT = ios)
!
            IF (ios /= 0) THEN
            WRITE (*,*) 'ERROR OPENING FILE IN LOADBM',trim(NameRcmDir),filnam
            STOP
         END IF
!
!
      IF (form_char == 'FORMATTED') THEN
!
!        READ THE ID INFO FOR THIS PARTICULAR B-FIELD MODEL:
!
         READ (LUN,800,IOSTAT=IOS) idf1, idf2, idf3, idf4, idf5, idf6
         IF (ios /= 0) THEN
            WRITE (*,*) 'ERROR READING ID IFO FROM FILE IN LOADBM,',filnam
            STOP
         END IF
!
!
!        CHECK IF THIS IS THE MODEL WE ACTUALLY REQUESTED:
!
         IF ( (bfndx(i,1) /= idf2) .or. (bfndx(j,2) /= idf1) .or.  &
              (bfndx(k,3) /= idf4) .or. (bfndx(l,4) /= idf3) .or.  &
              ((3-bfndx(m,5)) /= idf5) .or. (idf6 /= 501)) THEN
            WRITE (*,*) 'FILE   ',FILNAM,' IS INCORRECT.'
            WRITE (*,*) 'FILE ',FILNAM,' IS INCORRECT.'
            STOP
         END IF
!
!
!        READ DIFFERENT TYPES OF MATRICES FOR THE BFIELD MODEL:
!
         DO n = 1, iwant
            READ (LUN,801,IOSTAT=IOS) &
            ((work(ii,jj,i,j,k,l,m), ii=1,latdim),jj=1,ltdim)
            IF (ios /= 0) THEN
               WRITE (*,*) 'ERROR READING MATRICES IN LOADBM, FILE ', filnam
               STOP
            END IF
         END DO
!
      ELSE
!
!        READ THE ID INFO FOR THIS PARTICULAR B-FIELD MODEL:
!
         READ (LUN,IOSTAT=IOS) idf1, idf2, idf3, idf4, idf5, idf6
         IF (ios /= 0) THEN
            WRITE (*,*) 'ERROR READING ID IFO FROM FILE IN LOADBM,',filnam
            STOP
         END IF
!
!
!        CHECK IF THIS IS THE MODEL WE ACTUALLY REQUESTED:
!
         IF ( (bfndx(i,1) /= idf2) .or. (bfndx(j,2) /= idf1) .or.  &
              (bfndx(k,3) /= idf4) .or. (bfndx(l,4) /= idf3) .or.  &
              ((3-bfndx(m,5)) /= idf5) .or. (idf6 /= 501)) THEN
            WRITE (*,*) 'FILE   ',FILNAM,' IS INCORRECT.'
            WRITE (*,*) 'FILE ',FILNAM,' IS INCORRECT.'
            STOP
         END IF
!
!
!        READ DIFFERENT TYPES OF MATRICES FOR THE BFIELD MODEL:
!
         DO n = 1, iwant
            READ (LUN, IOSTAT = ios) work (:,:,i,j,k,l,m)
            IF (ios /= 0) THEN
               WRITE (*,*) 'ERROR READING MATRICES IN LOADBM, FILE ', filnam
               STOP
            END IF
         END DO
!
      END IF
      CLOSE (LUN)
      RETURN
  800 FORMAT(2I1,2I2.2,I1,I3)
  801 FORMAT(11E12.4)
      END SUBROUTINE Loadbm
!
!
!
      SUBROUTINE Zerobm (latdim, ltdim, i, j, k, l, m, work)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: latdim, ltdim, i, j, k, l, m
      REAL, INTENT (IN OUT) :: work (latdim, ltdim, 2, 2, 2, 2, 2)
!
!
!  copyright Rice University, 1996
!
!  VERSION 1.00 MSM DELIVERY VERSION            02.02.93
!
!  PURPOSE
!   THIS SUBROUTINE ERASES AN INDIVIDUAL B-MATRIX WITHIN THE
!   WORKING B-MATRICES.  THIS IS DONE WHEN AN INDIVIDUAL MATRIX
!   DOES NOT EXIST IN THE OFFLINE B-SUPERMATRIX.
!
!  INPUTS:
!   LATDIM      NUMBER OF LATITUDINAL GRID LINES
!   LTDIM       NUMBER OF LOCAL TIME (LONGITUDINAL) GRID LINES
!   I,J,K,L,M   INDICIES OF WORKING MAGNETIC FIELD MATRICES
!   WORK        WORKING MAGNETIC FIELD MATRIX
!
!  OUTPUTS:
!   WORK        ZEROED OUT WORKING MAGNETIC FIELD MATRIX
!
!
      work (1:latdim, 1:ltdim, i, j, k, l, m) = 0.0
      RETURN
      END SUBROUTINE Zerobm
!
!
!
!
!
      FUNCTION Mexist (mstnd, mtilt, mined, mdst, mstch)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: mstnd, mtilt, mined, mdst, mstch
      LOGICAL :: Mexist
!
!
!  copyright Rice University, 1996
!
!   VERSION 1.00 MSM DELIVERY VERSION            02.02.93
!           1.10                                 09.27.94
!                ERR AND IOSTAT VARIABLES ADDED TO IO STATEMENTS
!
! PURPOSE
!   THIS FUNCTION CHECKS TO SEE IF A PARTICULAR B-MATRIX EXISTS
!   IN THE B-SUPERMATRIX(TM).
!
! INPUTS:
!   MSTND       STANDOFF PARAMETER INDEX, MAX IS 6
!   MTILT       TILT ANGLE PARAMETER INDEX, MAX IS 1
!   MINED       EQUATORWARD EDGE OF THE AURORAL OVAL PARAMETER INDEX, MAX IS 18
!   MDST        DST PARAMETER INDEX, MAX IS 10
!   MSTCH       TAIL COLLAPSE PARAMETER INDEX, MAX IS 2
!
! OUTPUTS
!   MEXIST      LOGICAL VARIABLE ON THE STATUS OF A MAGNETIC FIELD MATRIX
!
!
      CHARACTER (LEN=80) :: filnam
!
      filnam = Get_filnam (mstnd, mtilt, mined, mdst, mstch)
      INQUIRE (FILE = TRIM (filnam), EXIST = Mexist)
      RETURN
      END FUNCTION Mexist
!
!
!
      FUNCTION Get_filnam (mstnd, mtilt, mined, mdst, mstch)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: mstnd, mtilt, mined, mdst, mstch
      CHARACTER (LEN=80) :: Get_filnam
!
      WRITE (Get_filnam, format_prefix) prefix, mtilt, mstnd, mdst, mined, (3-mstch)
      RETURN
      END FUNCTION Get_filnam
!
!
!
!
      SUBROUTINE Btrace (latdim, ltdim, itmdim, l, fstoff, feqedg, fdst, &
                         fclpse, ftilt, vm, bmin, xmin, ymin, zmin)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: latdim, ltdim, itmdim, l
      REAL,    INTENT (IN) :: fstoff, feqedg, fdst, fclpse, ftilt
      REAL, INTENT (IN OUT) :: vm (latdim,ltdim,itmdim), &
                            bmin (latdim,ltdim,itmdim), &
                            xmin (latdim,ltdim,itmdim), &
                            ymin (latdim,ltdim,itmdim), &
                            zmin (latdim,ltdim,itmdim)
!
!
!  copyright Rice University, 1996
!
!  VERSION 2.02                         DATE: 09.16.89
!          3.00  MSM DELIVERY VERSION         JANUARY 27, 1993
!          3.10                               APRIL 27, 1993
!               CHANGED TO GIVE TIME DEPTH TO XMIN, YMIN, ZMIN,
!          3.20                               JUNE 11, 1993
!               ERROR PRINTOUT ROUTED TO UNIT 9
!          3.30                               SEP 15, 1993
!               CODE CHANGED SO THAT B-FIELD MATRICES ARE NEVER
!               EXTRAPOLATATIONS, BUT INSTEAD CORRESPOND TO MATRICES
!               APPROPRIATE FOR THE LIMITS OF PARAMETER SPACE.
!          4.00                               JULY 9,1996
!               STOPS REMOVED WITH CALL TO ERRSUB SUBROUTINE
!
!  PURPOSE:  SUBROUTINE TO CALCULATE B FIELD MODELS APPROPRIATE FOR
!          GEOPHYSICAL CONDITIONS AS SPECIFIED IN BFPAR.  BTRACE
!          CALCULATES VM, BMIN, XMIN, YMIN, AND ZMIN ON THE MSM SPATIAL
!          GRID BY INTERPOLATING BETWEEN PRECOMPUTED MAGNETIC FIELD
!          MODELS AS PARAMETERIZED BY THE VARIABLES OF BFPAR.
!
!  INPUT:
!        LATDIM    NUMBER OF LATITUDINAL GRID POINTS
!        LTDIM     NUMBER OF LOCAL TIME GRID POINTS (INCL WRAPAROUND)
!        ITMDIM    MAX NUMBER OF TEMPORAL GRID POINTS
!        L         CURRENT TEMPORAL GRID INDEX
!        .....     GEOPHYSICAL PARAMETERS THAT DETERMINE WHICH BFIELD
!                  MODELS TO USE
!        FSTOFF    STANDOFF DIST = AUGPAR(STAND)
!        FEQEDG    EQUATORWARD EDGE OF AURORA = AUGPAR(EQEDGE)
!        FDST      DST = AUGPAR(DST)
!        FCLPSE    COLLAPSE PARAMETER = AUGPAR(CLAPSE)
!        FTILT     TILT PARAMETER = AUGPAR(TILTW)
!        ALOCT     LOCAL TIME HOUR ANGLE (MEASURED EAST FROM NOON)
!        TETA      COLATITUDE OF I GRID LINES (RADIANS)
!
!  OUTPUT:
!        VM        (FLUX TUBE VOLUME)**-2/3 ARRAY (RE/NT)**-2/3
!        BMIN      B VALUE AT MINIMUM VALUE ALONG FIELD LINE (NT)
!        XMIN      X (GSM) VALUE AT WHICH BMIN OCCURS
!        YMIN      Y (GSM) VALUE AT WHICH BMIN OCCURS
!        ZMIN      Z (GSM) VALUE AT WHICH BMIN OCCURS
!        R         R SQRT(X**2+Y**2+Z**2) AT WHICH BMIN OCCURS
!        P         HOUR ANGLE (RADIANS) MEASURED EASTWARD FROM NOON
!
!  WORK:
!        WORK      WORK ARRAY (NEEDED FOR BFIELD RETRIEVAL AND
!                    INTERPOLATION)
!
!  REFERENCE:
!        FIVE DIMENSIONAL INTERPOLATION IS BASED ON MULTILINEAR
!          INTERPOLATION SCHEME GIVEN BY PRESS ET AL. IN
!          NUMERICAL RECIPES, CAMBRIDGE UNIVERSITY PRESS, 1987,PP 95-97
!
!
!  PROGRAMMER:  R.W. SPIRO
!
!
!
      REAL, PARAMETER :: BFBUSH = 1.0E+20
      INTEGER :: nn, iwant, i, j, ii, n1, n2, n3, n4, n5, n (5), m
      LOGICAL :: bfexst (2, 2, 2, 2, 2)
      REAL :: work (latdim, ltdim, 2, 2, 2, 2, 2), SUMC (latdim, ltdim),&
              bfpar (5), bftemp, t (2,5), coeff, bflim (2,5)
!
!
!      pi = ATAN2 (0.0,-1.0)
      bfpar (1:5)= (/ fstoff, ftilt, feqedg, fdst, fclpse /)
      vm   (:,:,l) = 0.0
      bmin (:,:,l) = 0.0
      xmin (:,:,l) = 0.0
      ymin (:,:,l) = 0.0
      zmin (:,:,l) = 0.0
      sumc (:,:)   = 0.0
!
      WRITE (*,*) 'GETMAT IS REQUESTING THE FOLLOWING B-MATRICES:'
!
!     LOOP OVER FIVE FLAVORS OF B FIELD ARRAYS
      DO nn = 1, 5
         iwant = nn
!
!        1. CALL GETMAT TO READ APPROPRIATE ARRAYS FOR INTERPOLATION
!
         CALL Getmat (iwant, latdim, ltdim, bfpar, work, bflim, bfexst)
         CALL Rmvbsh (work, latdim, ltdim, iwant)
!
!        2. SET UP WEIGHTING PARAMETERS FOR INTERPOLATION
!
         DO ii = 1, 5
            bftemp = bfpar(ii)
            IF (bftemp > bflim(2,ii)) bftemp = bflim(2,ii)
            IF (bftemp < bflim(1,ii)) bftemp = bflim(1,ii)
            t(2,ii) = (bftemp-bflim(1,ii))/(bflim(2,ii)-bflim(1,ii))
            t(1,ii) = 1.0-t(2,ii)
         END DO
!
         IF (iwant == 1) WRITE (*,*) 'T ARRAY',T
!
!
!        3.  MAIN INTERPOLATION LOOPS
!
         DO N1=1,2
         N(1)=N1
         DO N2=1,2
         N(2)=N2
         DO N3=1,2
         N(3)=N3
         DO N4=1,2
         N(4)=N4
         DO N5=1,2
            N(5)=N5
!
            COEFF = 1.0
            DO M = 1,5
               COEFF=COEFF*T(N(M),M)
            END DO
!
!           3.1 CHECK TO SEE IF B FIELD MODEL EXISTS
!
            IF (iwant == 1) THEN
               WRITE(*,1221) 'N1=', N1, 'N2=', N2, 'N3=', N3, &
                             'N4=', N4, 'N5=', N5,            &
                             'COEFF=', COEFF,                 &
                             'BFEXST=', BFEXST(N1,N2,N3,N4,N5)
1221           FORMAT (5(TR3,A3,I2),TR3,A6,F9.3, TR3,A7,L4)
            END IF
            IF (.NOT.Bfexst(n1,n2,n3,n4,n5) .AND. coeff /= 0.0) THEN
                WRITE (*,'(T2,A)',ADVANCE='NO') &
                       & 'one B-field does not exist. press ENTER...'
                READ (*,*)
            END IF
!            WRITE(*,*) 'COEFF',  COEFF
!            WRITE(*,*) 'BFEXST', BFEXST(N1,N2,N3,N4,N5)
            IF (.NOT.bfexst(n1,n2,n3,n4,n5)) coeff = 0.0
!
!
!           3.2 LOOP OVER SPATIAL GRID
!
            IF (iwant == 1) THEN
               DO j = 1, ltdim
               DO i = 1, latdim
                  IF (ABS(work(i,j,n1,n2,n3,n4,n5)) < bfbush) THEN
                     xmin(i,j,l) = xmin(i,j,l)+coeff*work(i,j,n1,n2,n3,n4,n5)
                  ELSE
                     WRITE(*,*) 'DEBUG PRINTOUT IN BTRACE'
                     WRITE(*,*) 'IWANT,I,J,N1,N2,N3,N4,NT',IWANT,I,J,N1,N2,N3,N4,N5
                     WRITE(*,*) 'WORK(I,J,N1,N2,N3,N4,N5)',WORK(I,J,N1,N2,N3,N4,N5)
                  END IF
               END DO
               END DO
            ELSE IF (iwant == 2) THEN
               DO j = 1, ltdim
               DO i = 1, latdim
                  IF (work(i,j,n1,n2,n3,n4,n5) < bfbush) THEN
                     ymin(i,j,l)=ymin(i,j,l)+coeff*work(i,j,n1,n2,n3,n4,n5)
                  END IF
               END DO
               END DO
            ELSE IF (iwant == 3) THEN
               DO j = 1, ltdim
               DO i = 1, latdim
                  IF (work(i,j,n1,n2,n3,n4,n5) < bfbush) THEN
                     zmin(i,j,l)=zmin(i,j,l)+coeff*work(i,j,n1,n2,n3,n4,n5)
                  END IF
               END DO
               END DO
            ELSE IF (iwant == 4) THEN
               DO j = 1, ltdim
               DO i = 1, latdim
                  IF (work(i,j,n1,n2,n3,n4,n5) < bfbush) THEN
                     bmin(i,j,l)=bmin(i,j,l)+coeff*work(i,j,n1,n2,n3,n4,n5)
                  END IF
               END DO
               END DO
            ELSE
               DO j = 1, ltdim
               DO i = 1, latdim
                  IF (work(i,j,n1,n2,n3,n4,n5) < bfbush) THEN
                     vm(i,j,l)=vm(i,j,l)+coeff*work(i,j,n1,n2,n3,n4,n5)
                     sumc(i,j)=sumc(i,j)+coeff
                  END IF
               END DO
               END DO
            END IF
         END DO
         END DO
         END DO
         END DO
         END DO
!
      END DO
!
!
!     4. CORRECT FOR MISSING B FIELD MATRICES
!
      DO j = 1, ltdim
         DO i = 1, latdim
            IF (sumc(i,j) > 1.0E-10) THEN
               vm(i,j,l)   = vm(i,j,l)   / sumc(i,j)
               bmin(i,j,l) = bmin(i,j,l) / sumc(i,j)
               xmin(i,j,l) = xmin(i,j,l) / sumc(i,j)
               ymin(i,j,l) = ymin(i,j,l) / sumc(i,j)
               zmin(i,j,l) = zmin(i,j,l) / sumc(i,j)
            END IF
         END DO
      END DO
!
      RETURN
      END SUBROUTINE Btrace
!
   END MODULE Hilmer_Voigt_subs_mod
!
!
!
!
    PROGRAM Setup_bf
    USE Rcm_variables
    USE Rcm_io
    USE Hilmer_Voigt_subs_mod
    IMPLICIT NONE
!
    INTEGER (iprec) :: isize_msm, jsize_msm, jwrap_msm
    REAL (rprec) :: colat_msm (latdim,ltdim), aloct_msm (latdim,ltdim)
!
    INTEGER (iprec) :: n, nbf, i, j, i_msm, j_msm, istat, ncol=80
    REAL (rprec) :: bbi, bbj, &
                    Gntrp_2d_old, Dipole_bfield
!
    LOGICAL :: logic_flag_1, ST, use_dipole
!
!
!   READ MSM'S GRID
!
    OPEN (UNIT = LUN, FILE=trim(NameRcmDir)//'rcmcrd_msm', STATUS = 'OLD', &
          FORM = 'UNFORMATTED', IOSTAT = istat)
       IF (istat /= 0) THEN
          WRITE (*,*) 'ERROR OPENING ',trim(NameRcmDir),'rcmcrd_msm'
          STOP
       ELSE
          READ (LUN) isize_msm, jsize_msm, jwrap_msm
          IF (isize_msm /= latdim .OR. jsize_msm /= ltdim) THEN
              STOP 'DIMENSIONS OF MSM GRID IN rcmcrd_msm differ from THOSE HERE'
          END IF
          READ (LUN) colat_msm
          READ (LUN) aloct_msm
       END IF
    CLOSE (LUN)
!
!
      OPEN (UNIT = LUN_2, FILE = trim(NameRcmDir)//'setup_bf.list', STATUS = 'REPLACE',&
            ACTION = 'WRITE', FORM = 'FORMATTED', ACCESS = 'SEQUENTIAL')
!
    CALL Read_grid ()
!
!
    WRITE (*,'(A)',ADVANCE='NO') 'USE DIPOLAR FIELD? (T/F)_______'
    READ (*,*) use_dipole
! 
!
  IF (use_dipole) THEN
!  
    nbf = 1
!
  ELSE
!
!   Read number of bfields from 'definebf.dat':
!
    OPEN (LUN, FILE = trim(NameRcmDir)//'definebf.dat', STATUS = 'OLD')
    READ (LUN,*)
    READ (LUN,*)
    nbf = 0
    DO
        READ (LUN,'(A)', END = 20)
        nbf = nbf + 1
    END DO
 20 CLOSE (LUN)
!
  END IF
!
    ALLOCATE (ibtime (nbf))
!
!
!
!   Extract the BF-matrices from HILMER files:
!
  IF (.NOT.use_dipole) THEN
    OPEN (LUN_3, FILE = trim(NameRcmDir)//'definebf.dat', STATUS = 'OLD')
    READ (LUN_3,*)
    READ (LUN_3,*)
  END IF
!
    DO n = 1, nbf
!
       IF (.NOT. use_dipole) THEN
          READ (LUN_3,*) ibtime(n), fstoff, fmeb, fdst, fclps, ftilt
!
!         For each model, get hilmer-voight matrices on msm grid:
!
          CALL Btrace (latdim, ltdim, 1, 1,                &
                         fstoff, fmeb, fdst, fclps, ftilt, &
                         vm_msm , be_msm, xmin_msm, ymin_msm, zmin_msm )
!
!
!         NOW DEFINE BFIELD MATRICES ON OUR GRID:
!
          DO j = 1, jsize
          DO i = 1, isize
!
             IF (i >= latdim) THEN
!
!               Point is below MSM'S low-lat bndy, append dipole B:
!
                rmin (i,j)    = 1.0_rprec / SIN(colat(i,j))**2
                pmin (i,j)    = aloct (i,j)
                xmin (i,j) = rmin (i,j) * COS(pmin(i,j))
                ymin (i,j) = rmin (i,j) * SIN(pmin(i,j))
                bmin (i,j)   = besu * (1.0 / rmin(i,j))**3
                vm (i,j)   = (32./35.* rmin(i,j)**4 / Besu *  &
                              SQRT(1.-1./rmin(i,j))* &
                              (1.+0.5/rmin(i,j)+3./ 8./rmin(i,j)**2+ &
                              5./16./rmin(i,j)**3) &
                             ) ** (-2.0/3.0)
!
             ELSE
!
!               INTERPOLATE WITHIN THE MSM'S GRID:
!
                aloct_msm (:,jsize_msm) = aloct_msm(:,jwrap_msm) + 2.0*pi
                j_msm =  1
                logic_flag_1 = .FALSE.
                one_loop: DO i_msm = 1, isize_msm-1
                     IF ( colat(i,j) >= colat_msm (i_msm, j_msm+2) .AND. &
                          colat(i,j) <  colat_msm (i_msm+1,j_msm+2) ) THEN
                          logic_flag_1 = .TRUE.
                          EXIT one_loop
                     END IF
                END DO one_loop
                IF (.NOT. logic_flag_1) THEN
                   IF (ABS(colat(i,j)-colat_msm(isize_msm,j_msm+2)) < 1.0E-5) THEN
                       i_msm = isize_msm-1
                   ELSE
                       STOP 'PROBLEM 1'
                   END IF
                END IF
!
                logic_flag_1 = .FALSE.
                two_loop: DO j_msm = jwrap_msm, jsize_msm-1
                     IF ( aloct(i,j) >= aloct_msm (i_msm, j_msm) .AND. &
                          aloct(i,j) <  aloct_msm (i_msm,j_msm+1) ) THEN
                          logic_flag_1 = .TRUE.
                          EXIT two_loop
                     END IF
                END DO two_loop
                IF (.NOT. logic_flag_1) STOP 'PROBLEM 2'
!
!               We bracketed (i,j) point on the msm grid by (i_msm,j_msm),
!               (i_msm+1,j_msm), (i_msm,j_msm+1), and (i_msm+1,j_msm+1).
!
!
                bbi = ( i_msm    * colat_msm(i_msm+1,j_msm) - &
                       (i_msm+1) * colat_msm(i_msm  ,j_msm) + colat(i,j) ) / &
                     (colat_msm(i_msm+1,j_msm) - colat_msm(i_msm,j_msm))
                bbj = (j_msm * aloct_msm(i_msm,j_msm+1)- &
                     (j_msm+1)*aloct_msm(i_msm,j_msm) + aloct(i,j) ) / &
                     (aloct_msm(i_msm,j_msm+1) - aloct_msm (i_msm,j_msm))
                vm (i,j)   = Gntrp_2d_old (vm_msm, isize_msm, jsize_msm, &
                                      jwrap_msm, bbi, bbj)
                bmin (i,j) = Gntrp_2d_old (be_msm, isize_msm, jsize_msm, &
                                       jwrap_msm, bbi, bbj)
                xmin (i,j) = Gntrp_2d_old (xmin_msm, isize_msm, jsize_msm, &
                                       jwrap_msm, bbi, bbj)
                ymin (i,j) = Gntrp_2d_old (ymin_msm, isize_msm, jsize_msm, &
                                       jwrap_msm, bbi, bbj)
!
                aloct_msm (:,jsize_msm) = aloct_msm(:,jwrap_msm) - 2.0*pi
!
             END IF
!
          END DO
          END DO
!
       ELSE
!
          ibtime (1) = 99999
!
          DO i = 1, isize
          DO j = 1, jsize
             xmin (i,j) = Dipole_Bfield (colat(i,j), aloct(i,j), 1)
             ymin (i,j) = Dipole_Bfield (colat(i,j), aloct(i,j), 2)
             bmin (i,j) = Dipole_Bfield (colat(i,j), aloct(i,j), 3)
             vm   (i,j) = Dipole_Bfield (colat(i,j), aloct(i,j), 4)
          END DO
          END DO
!
       END IF
!
       CALL Wrap_around_ghostcells (vm, isize, jsize, n_gc)
       CALL Wrap_around_ghostcells (xmin, isize, jsize, n_gc)
       CALL Wrap_around_ghostcells (ymin, isize, jsize, n_gc)
       CALL Wrap_around_ghostcells (bmin, isize, jsize, n_gc)
!
       label%intg(6)  = ibtime(n)
       label%real(12) = fmeb
       label%real(13) = fstoff
       label%real(14) = fdst
       label%real(15) = fclps
       label%real(16) = ftilt
!
       label%intg (20)= nbf
!
       label%char = 'OCTOBER_1998 RUN'
!
       CALL Outp (vm,   1,isize,1,1,jsize,1, 0.0, label%intg,'VM  ', LUN_2, ncol)
       CALL Outp (bmin, 1,isize,1,1,jsize,1, 0.0, label%intg,'BMIN', LUN_2, ncol)
       CALL Outp (xmin, 1,isize,1,1,jsize,1, 0.0, label%intg,'XMIN', LUN_2, ncol)
       CALL Outp (ymin, 1,isize,1,1,jsize,1, 0.0, label%intg,'YMIN', LUN_2, ncol)
!
       IF (n == 1) THEN
          ST = .TRUE.
       ELSE
          ST = .FALSE.
       END IF
       CALL Write_array ('rcmbmin_inp', n, label, bmin, SETUP = ST, ASCI = asci_flag)
       CALL Write_array ('rcmvm_inp',   n, label, vm,   SETUP = ST, ASCI = asci_flag)
       CALL Write_array ('rcmxmin_inp', n, label, xmin, SETUP = ST, ASCI = asci_flag)
       CALL Write_array ('rcmymin_inp', n, label, ymin, SETUP = ST, ASCI = asci_flag)
!
    END DO
    CLOSE (UNIT = LUN_3)
    CLOSE (UNIT = LUN_2 )
!
    STOP
    END PROGRAM Setup_bf

    FUNCTION Gntrp_2d_old (array, i_max, j_max, jwrap, bi, bbj)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i_max, j_max, jwrap
    REAL (rprec), INTENT (IN)    :: array (i_max,j_max), bi, bbj
    REAL (rprec)                 :: Gntrp_2d_old
!
    INTEGER (iprec) :: ii, jn, jj, jp1
    REAL (rprec)    :: fi,fj,a1,a2,bj, &
                       Bjmod_real_old
!
!
    ii = MAX (1, MIN (INT (bi), i_max-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
    bj    = Bjmod_real_old ( bbj, jwrap, j_max)
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp.
!
      Gntrp_2d_old = (1.0_rprec-(bi-fi)) * array(ii,jn) + &
                     (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interp.
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (1.0_rprec-(bi-fi))*array(ii,jj) + (bi-fi)*array(ii+1,jj)
       a2 = (1.0_rprec-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
       Gntrp_2d_old = (1.0_rprec - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Gntrp_2d_old
!

    FUNCTION Bjmod_real_old (bj, jwrap, jsize)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    REAL (rprec),    INTENT (IN) :: bj
    INTEGER(iprec),  INTENT (IN) :: jwrap, jsize
    REAL (rprec)                 :: Bjmod_real_old
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_real_old = bj
!                                                                       
    do_1: DO
       IF (Bjmod_real_old > REAL (jsize - 1,rprec)) THEN
          Bjmod_real_old = Bjmod_real_old - REAL (jsize - jwrap,rprec)
       ELSE
          EXIT do_1
       END IF
    END DO do_1
!
    do_2: DO
       IF (Bjmod_real_old < REAL (jwrap,rprec)) THEN
           Bjmod_real_old = Bjmod_real_old + REAL(jsize - jwrap,rprec)
       ELSE
           EXIT do_2
       END IF
    END DO do_2
!
    RETURN
    END FUNCTION Bjmod_real_old
!
