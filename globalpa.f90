! -*- Fortran -*-
real omega0, omegab, compressmax, hubbleh0, redshiftzi, boxsize, &
     sigma8, omegal, rcompress1, rcompress2,rcut0
integer  magic, iseed(4), icrecut(3)
logical lrecut
common /globalpa/ omega0, omegal, omegab, compressmax, &
         rcompress1, rcompress2, rcut0, &
         hubbleh0, redshiftzi, boxsize, sigma8, iseed, magic, &
         icrecut, lrecut
! \Omega in baryons (i.e. gas)
! \Omega_0 = \Omega_{DM} + \Omega_b
! hubbleh0 is H_0 in units of 100 km/s/Mpc
! boxsize is in units of h^{-1} Mpc.      
