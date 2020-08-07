!
!======================================================================
! This file defines variables related to Symmetric instability
! The variables here is global and used for CROCO (www.croco-ocean.org)
! 
! Other related files:
!       si_scheme.F: main code file calculating mixing due to symmetric instability
!                    
!    lmd_skpp1994.F: kpp scheme based on Large et al. (1994) 
!
!    lmd_skpp2005.F: kpp scheme based on Shchepetkin (2005)
!
!        step3d_t.F: update temperature and salt due to convection related to surface buoyancy loss
!
!      t3dmix_ISO.F: update isopycnal diffusivity, if you want to turn on SI_SCHEME, you need to define TS_MIX_ISO 
!
!======================================================================
! written by Jihai Dong, jihai_dong@nuist.edu.cn
!======================================================================
/* This is module file. 
  --------------------------------------------
*/

/*# if defined SI_SCHEME*/
       real Kiso_si(GLOBAL_2D_ARRAY, 1:N)                 /*isopycnal diffusivity*/
       real Ftemp_si(GLOBAL_2D_ARRAY, 1:N)                /*vertical convection of temperature*/
       real Fsalt_si(GLOBAL_2D_ARRAY, 1:N)                /*vertical convection of salt*/

       common / si_temp / Ftemp_si
       common / si_salt / Fsalt_si
       common / si_kiso / Kiso_si

/*# endif*/
