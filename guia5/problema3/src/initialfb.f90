! External electromotive forcing.
! This file contains the expression used for the external 
! electromotive forcing. You can use temporary real arrays 
! R1-R3 of size (1:nx,1:ny,ksta:kend) and temporary complex 
! arrays C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable m0 should control the global 
! amplitude of the forcing, and variables mparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays mx, my, and mz.

! Null electromotive forcing

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               mx(k,j,i) = 0.0_GP
               my(k,j,i) = 0.0_GP
               mz(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
