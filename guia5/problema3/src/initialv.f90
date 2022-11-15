! Initial condition for the velocity.
! This file contains the expression used for the initial 
! velocity field. You can use temporary real arrays R1-R3 
! of size (1:n,1:n,ksta:kend) and temporary complex arrays 
! C1-C8 of size (n,n,ista:iend) to do intermediate 
! computations. The variable u0 should control the global 
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Superposition of harmonic modes with random phases 
! correlated to give a tunable total relative helicity 
! [after Pouquet & Patterson, JFM 1978] with a ~k^(-2) slope
!     kdn    : minimum wave number
!     kup    : maximum wave number
!     vparam0: sin(2*vparam0) is the total relative helicity
!     vparam8, vparam9 are used for temporal storage

      vparam0 = 0.0_GP
  
      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))**3
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,1,1))**3
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,j,1))**3
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,C4,1)
      CALL rotor3(C1,C3,C5,2)
      CALL rotor3(C1,C2,C6,3)
      CALL normalize(C4,C5,C6,u0,1,MPI_COMM_WORLD)

!C4,C5,C6 are the components of the first random vector (say v1)
!Repeating this procedure we create the second random vector (v2)

      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))**3
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,1,1))**3
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,j,1))**3
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,C7,1)
      CALL rotor3(C1,C3,C8,2)
      CALL rotor3(C1,C2,C3,3)
      CALL normalize(C7,C8,C3,u0,1,MPI_COMM_WORLD)

! So far, v1 = (C4,C5,C6) and v2 = (C7,C8,C3) are two random 
! normalized vectors in Fourier space. Correlating vectors 
! v1 and v2 we create the velocities vectors vx, vy, and vz.

      vparam8=SIN(vparam0)
      vparam9=COS(vparam0)

! vx: components y and z of u2 are calculated in order to get 
! the x component of curl(u2)

      DO i = ista,iend
          DO j = 1,ny
             DO k = 1,nz
                IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = vparam8*C5(k,j,i)+vparam9*C8(k,j,i)
                  C2(k,j,i) = vparam8*C6(k,j,i)+vparam9*C3(k,j,i)
                ENDIF
             END DO
          END DO
      END DO

      CALL rotor3(C1,C2,C1,1)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  vx(k,j,i) = vparam9*C4(k,j,i)+vparam8*C7(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

! vy: components x and z of u2 are calculated in order to get 
! the y component of curl(u2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = vparam8*C4(k,j,i)+vparam9*C7(k,j,i)
                  C2(k,j,i) = vparam8*C6(k,j,i)+vparam9*C3(k,j,i)
               ENDIF
            END DO
         END DO
      END DO

      CALL rotor3(C1,C2,C1,2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  vy(k,j,i) = vparam9*C5(k,j,i)+vparam8*C8(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

! vz: components x and y of u2 are calculated in order to get
! the z component of curl(u2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = vparam8*C4(k,j,i)+vparam9*C7(k,j,i)
                  C2(k,j,i) = vparam8*C5(k,j,i)+vparam9*C8(k,j,i)
               ENDIF
            END DO
         END DO
      END DO

      CALL rotor3(C1,C2,C1,3)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  vz(k,j,i) = vparam9*C6(k,j,i)+vparam8*C3(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

      CALL normalize(vx,vy,vz,u0,1,MPI_COMM_WORLD)
