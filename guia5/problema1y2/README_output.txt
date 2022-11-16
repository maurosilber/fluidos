Output for this solver includes the following files:
(note some output files may not be created, depending
on specific solver and options being used).

From subroutines in pseudospec3D_mod.f90 (see file for more details):

From subroutines in module_dns.f90 (see file for more details):

From subroutines in pseudospec3D_hd.f90 (see file for more details):
! Output files contain:
! 'balance.txt':  time, <v^2>, <omega^2>, mechanic injection rate
! 'helicity.txt': time, kinetic helicity
! 'divergence.txt' [OPTIONAL]: time, <(div.v)^2>
!

! Output files contain:
! 'kspectrum.XXX.txt': k, Ev(k)
! 'mspectrum.XXX.txt': k, Eb(k)
! 'khelicity.XXX.txt': k, Hv(k) (kinetic helicity spectrum)
! 'mhelicity.XXX.txt': k, Hb(k) (magnetic helicity spectrum)
! 'ghelicity.XXX.txt': k, G(k)  (generalized helicity in Hall-MHD)
!

! Output files contain:
! 'Pspec1dIJ.XXX.txt' [I=x,y,z (field component)][J=x,y,z (dir)]
!     [P=prefix: 'k' for kinetic energy, 'm' for magnetic, etc.]:
!     k_dir, E_comp(k_dir)
!

! Output files contain:
! 'ktransfer.XXX.txt': k, Tv(k) (kinetic energy transfer function)
! 'mtransfer.XXX.txt': k, Tb(k) (magnetic energy transfer)
! 'jtransfer.XXX.txt': k, Tj(k) (Lorentz force work)
! 'kcrostran.XXX.txt': k, TC_v(k) (kinetic cross-helicity transfer)
! 'mcrostran.XXX.txt': k, TC_b(k) (magnetic cross-helicity transfer)
!

! Output files contain:
! 'hktransfer.XXX.txt': k, TH_v(k) (transfer of kinetic helicity)
! 'hmtransfer.XXX.txt': k, TH_b(k) (transfer of magnetic helicity)
!


From subroutines in pseudospec3D_rot.f90 (see file for more details):
! Output files contain:
! 'kspecpara.XXX.txt': kz, Ev(kz), Ev_perp(kz), Ev_z(kz)
!   [Ev: kinetic energy; Ev_perp: energy in vx,vy; Ev_z: same in vz]        
! 'mspecpara.XXX.txt': kz, Eb(kz), Eb_perp(kz), Eb_z(kz)
! 'khelipara.XXX.txt': kz, Hv(kz), Hv_perp(kz), Hv_z(kz)
!   [Hv: kinetic helicity; Hv_perp: v_perp.w_perp, Hv_z: vz.wz]
! 'mhelipara.XXX.txt': kz, Hb(kz), Hb_perp(kz), Hb_z(kz)
! 'ghelipara.XXX.txt': kz, G(kz)  ,G_perp(kz) , G_z(kz) 
!   [Generalized helicity in Hall-MHD]
!

! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'kspecperp.XXX.txt': kp, E_v(kp), ev_x,y(kp,kz=0), ev_z(kp,kz=0)
!   [E_v: kin. ener.; ev_x,y: 2D spec. for vx,vy; ev_z: same for vz]
! 'mspecperp.XXX.txt': kp, E_b(kp), eb_x,y(kp,kz=0), eb_z(kp,kz=0)
! 'kheliperp.XXX.txt': kp, H_v(kp), hv_x,y(kp,kz=0), hv_z(kp,kz=0)
! 'mheliperp.XXX.txt': kp, H_b(kp), hb_x,y(kp,kz=0), hb_z(kp,kz=0)
! 'gheliperp.XXX.txt': kp, G(kz),   g(kp,kz=0),      g(kp,kz=0)
!

! Output files contain:
! 'ktranpara.XXX.txt': kz, Tv(kz) (kinetic energy transfer function)
! 'mtranpara.XXX.txt': kz, Tb(kz) (magnetic energy transfer)
! 'jtranpara.XXX.txt': kz, Hj(kz) (Lorentz force work)
!

! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'ktranperp.XXX.txt': kp, Tv(kp) (kinetic energy transfer function)
! 'mtranperp.XXX.txt': kp, Tb(kp) (magnetic energy transfer)
! 'jtranperp.XXX.txt': kp, Hj(kp) (Lorentz force work)
!

! Output files contain:
! 'hktranpara.XXX.txt': kz, TH_v(kz) (kinetic helicity transfer)
! 'hmtranpara.XXX.txt': kz, TH_b(kz) (magnetic helicity transfer)
!

! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'hktranperp.XXX.txt': kp, TH_v(kp) (kinetic helicity transfer)
! 'hmtranperp.XXX.txt': kp, TH_b(kp) (magnetic helicity transfer)
!

! Output files contain:
! 'odir/kspec2D.XXX.out': kinetic energy 2D spectrum ev(kperp,kpara)
! 'odir/mspec2D.XXX.out': magnetic energy spectrum   eb(kperp,kpara)
! 'odir/kheli2D.XXX.out': kinetic helicity spectrum  hv(kperp,kpara)
! 'odir/mheli2D.XXX.out': magnetic helicity spectrum hv(kperp,kpara)
! 'odir/gheli2D.XXX.out': generalized helicity spec. g(kperp,kpara)
!        

! Output files contain [field = vx,vy,...]:
! 'odir/field_kx': field(kz,ky,kx=0) with size (nz,ny)
! 'odir/field_ky': field(kz,ky=0,kx) with size (nz,nx/2+1)
! 'odir/field_kz': field(kz=0,ky,kx) with size (ny,nx/2+1)
!   [Wavenumbers are k_i=(0,1,...)/L_i (i=x,y,z)]
!


