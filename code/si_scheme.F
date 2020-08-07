!
!======================================================================
! This file is the main F-file related to Symmetric instability,
! computing viscocity, diffusivities
! 
! Other related files:
!                    
!    lmd_skpp1994.F: kpp scheme based on Large et al. (1994)
!      
!    lmd_skpp2005.F: kpp scheme based on Shchepetkin (2005)
!
!              si.h: define global variables related to symmetric instability
!
!        step3d_t.F: update temperature and salt due to convection related to surface buoyancy loss
!
!      t3dmix_ISO.F: update isopycnal diffusivity, if you want to turn on SI_SCHEME, you need to define TS_MIX_ISO 
!
!======================================================================
! written by Jihai Dong, jihai_dong@nuist.edu.cn
!====================================================================== 
!================================================================================
!
#include "cppdefs.h"

      subroutine si_scheme (tile)

!

!      
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
 
      trd=omp_get_thread_num()
      call si_scheme_tile(Istr,Iend,Jstr,Jend)
      
      return
      end subroutine si_scheme
      
      
      subroutine si_scheme_tile(Istr,Iend,Jstr,Jend)  
      
!
# include "param.h"
# include "private_scratch.h"     
# include "grid.h"
# include "ocean3d.h"
# include "mixing.h"
# include "scalars.h"
# include "forces.h"   
# include "si.h"      
      
!      
      integer  ::Istr,Iend,Jstr,Jend                                  !index for each tile， PRIVATE_2D_SCRATCH_ARRAY Istr-2:Iend+2,Jstr-2:Jend+2
!      real     ::Kv      (PRIVATE_2D_SCRATCH_ARRAY, 0:N),             !vertical viscosity, on RHO grid horizontally and W grid vertically
!     &           Kt      (PRIVATE_2D_SCRATCH_ARRAY, 0:N),             !vertical diffusion for temperature, same with Kv
!     &           Ks      (PRIVATE_2D_SCRATCH_ARRAY, 0:N)              !vertical diffusion for salt, same with Kv
# define tind nstp      
      
      !local index only used in this subroutine
      
      integer  ::i,j,k,k_si                                                !
!     
      !parameters directly relatied to  quartic polynomial solving
      
      real     :: alpha_si,                                           !parameter for solving quartic polynomial
     &            x_si    ,                                           !ratio h_si/H_si, real root of quartic polynomial
     &            wstar_si,                                           !convective velocity 
     &            ustar_si,                                           !friction velocity
     &            U_si    ,                                           !along-front (i.e., geostrophic) velocity difference within H_si
     &            theta_si,                                           !angle between along-front direction and wind direction
     &            H_si    ,                                           !total SI layer thickness
     &            hc_si,                                              !upper convective layer thickness 
     &            c_si=14.0                                           !constant in quartic polynomial
     
      
!      parameter ( c_si=14, pi=3.1415927 )                             
     
      
      !some related variables 

      real     :: dUde_si(Istr:Iend+1,Jstr:Jend+1),                   !y-direction gradient of u, on PSI grid 
     &            dVdx_si(Istr:Iend+1,Jstr:Jend+1),                   !x-direction gradient of v, on PSI grid
     &            rvor_si(Istr:Iend+1,Jstr:Jend+1),                   !vertical vorticity on PSI grid
     &            rvor_r (Istr:Iend  ,Jstr:Jend  ),                   !vertical vorticity on RHO grid
     &            dRdx_si(Istr:Iend  ,Jstr:Jend  ),                   !x-direction gradient of buoyancy, on RHO grid
     &            dRde_si(Istr:Iend  ,Jstr:Jend  ),                   !y-direction gradient of buoyancy, on RHO grid
     &            dRdz_si(Istr:Iend  ,Jstr:Jend  ),                   !vertical difference of buoyance at surface, on RHO grid
     &            dUdz_si(Istr:Iend  ,Jstr:Jend  ),                   !vertical difference of u at surface, on RHO grid                    
     &            dVdz_si(Istr:Iend  ,Jstr:Jend  ),                   !vertical difference of v at surface, on RHO grid
     &            alpha  (Istr:Iend  ,Jstr:Jend  ),                   !vertical difference of u at surface, on RHO grid                    
     &            beta   (Istr:Iend  ,Jstr:Jend  ),                   !vertical difference of v at surface, on RHO grid
     &            GSP_si(0:N),                                        !conversion rate from background kinetic energy to SI
     &            wtemp_si(0:N),                                      !vertical convection of temperature
     &            wsalt_si(0:N)                                       !vertical convection of salt
      
      real     :: B0,ustr0,vstr0,orho0,pvor,EBF0,dRdx0,dRde0,pvor0,  
     &            dRdx_r,dRde_r,dUdz0,dVdz0,dRdz0,thetaB,dRdb_si,    
     &            thetaU,thetaW,zw,dRdx_Ri,dRde_Ri,                  
     &            Akv_si,Rib,Akt_si,cff,cff1
      

!================================================================================
!    
      

      call alfabeta_si (Istr,Iend,Jstr,Jend, alpha,beta)
      
      !---------firstly, relative vertical vortiticy-----------------

      !calculate surface velocity gradient, u_y and v_x, on PSI

      DO j=Jstr,Jend+1   
          DO i=Istr,Iend+1  
              dUde_si(i,j)=om_u(i,j  )*u(i,j  ,N,tind)-              
     &                  om_u(i,j-1)*u(i,j-1,N,tind)   
#  ifdef MASKING
              dUde_si(i,j)=dUde_si(i,j)*pmask(i,j)
#  endif
              dVdx_si(i,j)=on_v(i  ,j)*v(i  ,j,N,tind)-              
     &                  on_v(i-1,j)*v(i-1,j,N,tind)
#  ifdef MASKING
              dVdx_si(i,j)=dVdx_si(i,j)*pmask(i,j)
#  endif
          END DO
      END DO
!
      !calculate vertical vorticity
!
      DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
               cff=pm(i,j)*pn(i,j)
               rvor_si(i,j)=cff*(dVdx_si(i,j)-dUde_si(i,j))   ! vertical vorticity
            END DO
      END DO
!
      !move vorticity onto rho grid
!
      DO j=Jstr,Jend        
          DO i=Istr,Iend
              rvor_r(i,j)=0.25*(rvor_si(i,j  )+rvor_si(i+1,j  )+     
     &                          rvor_si(i,j+1)+rvor_si(i+1,j+1))
          END DO
      END DO
 
      
      !---------secondly, buoyancy gradient----------------------- 
      
      orho0=1.0_r8/rho0
         
      DO j=Jstr,Jend
          DO i=Istr,Iend
              ! x direction
              cff =0.5_r8*(pm(i-1,j)+pm(i  ,j))
              cff1=0.5_r8*(pm(i  ,j)+pm(i+1,j))
              !
              dRdx_si(i,j)=-0.5*g*orho0*
     &                     (cff*(rho(i  ,j,N)-rho(i-1,j,N))+       !rho1 is the potential density anomaly (-1000), rho is in-situ density, refer to rho_eos.F         
     &                     cff1*(rho(i+1,j,N)-rho(i  ,j,N)))                              
              ! y direction
              cff =0.5_r8*(pn(i,j-1)+pn(i,j  ))
              cff1=0.5_r8*(pn(i,j  )+pn(i,j+1))
              dRde_si(i,j)=-0.5*g*orho0*
     &                     (cff*(rho(i,j  ,N)-rho(i,j-1,N))+              
     &                     cff1*(rho(i,j+1,N)-rho(i,j  ,N)))
          END DO
      END DO

      !---------------thirdly, vertical difference-------------------
      !--------------------Here, we use bulk PV----------------------
      !---------its difference in vertical, not gradient-------------
!            
      DO j=Jstr,Jend        ! buoyancy
          DO i=Istr,Iend
              cff=-g*orho0
              dRdz_si(i,j)=cff*(rho(i  ,j  ,N)-rho(i  ,j  ,N-1))
          END DO
      END DO
         
      DO j=Jstr,Jend       ! u and v
          DO i=Istr,Iend
              cff=0.5_r8
              dUdz_si(i,j)=cff*(u(i  ,j,N,tind)-u(i  ,j,N-1,tind)+   
     &                          u(i+1,j,N,tind)-u(i+1,j,N-1,tind))
              dVdz_si(i,j)=cff*(v(i,j  ,N,tind)-v(i,j  ,N-1,tind)+   
     &                          v(i,j+1,N,tind)-v(i,j+1,N-1,tind))
          END DO
      END DO
          
!================================================================================
!-------------------let initial value as zero------------------------------------  
     ! print *,AKv(1:30,30,18)     
          
      DO j=Jstr,Jend
         DO i=Istr,Iend
             DO k=0,N
                 Ftemp_si(i,j,k)=0.0
                 Fsalt_si(i,j,k)=0.0
                  Kiso_si(i,j,k)=0.0
             END DO
         END DO
      END DO      
      
!================================================================================ 
!=======================main code to calculate SI mixing=========================
!--------------------------------------------------------------------------------
      DO j=Jstr,Jend                                                  !loop each grid point, to determine
        DO i=Istr,Iend                                                !whether SI is applied
   !  
            !----------------criterior: surface bulk PV > 0 ---------
            
            pvor=f(i,j)*((f(i,j)+rvor_r(i,j))*dRdz_si(i,j)+         
     &                   dUdz_si(i,j)*dRde_si(i,j)-
     &                   dVdz_si(i,j)*dRdx_si(i,j))
            B0= -g*alpha(i,j)*stflx(i,j,itemp)
            
# ifdef SALINITY
            B0= -g*( alpha(i,j)*stflx(i,j,itemp)
     &        -beta(i,j)*stflx(i,j,isalt))                           !surface buoyancy flux
            
# endif
            ustr0=0.5_r8*(sustr(i,j)+sustr(i+1,j))                    !x-component wind stress
            vstr0=0.5_r8*(svstr(i,j)+svstr(i,j+1))                    !y-component wind stress
            ustar_si=SQRT(SQRT(ustr0**2.0+vstr0**2.0))                    !friction velocity
            cff=1.0_r8/(f(i,j))                                                                                                                 !!!============!+rvor_r(i,j))       
            EBF0=cff*(vstr0*dRdx_si(i,j)-ustr0*dRde_si(i,j))          !Ekman buoyancy flux
            
            IF (pvor.lt.0.0 .and. rmask(i,j).eq.1 .and.                
     &          B0.gt.0.0 .and. EBF0.gt.0.0 .and.
     &          rvor_r(i,j)/f(i,j).gt.-1.0)             THEN            !pvor is negative and water grid, grid with centrifugal instability is excluded 
               dRdx0 = dRdx_si(i,j)
               dRde0 = dRde_si(i,j)
               pvor0 = pvor   
               !print*,ustr0,vstr0,dRdx_si(i,j),dRde_si(i,j)
               k=N-1
               DO WHILE (pvor0.lt.0.0 .and. k.gt.1)                    !loop the depth from surface to get H_si, N→surface, 0→bottom
                  k=k-1
                  
                  !--------------buoyancy gradient------------------- 
                  cff =0.5_r8*(pm(i-1,j)+pm(i  ,j))
                  cff1=0.5_r8*(pm(i  ,j)+pm(i+1,j))
                  !
                  dRdx_r=-0.5*g*orho0*
     &                     (cff*(rho(i  ,j,k)-rho(i-1,j,k))+       !!x-direction gradient on RHO at k          
     &                     cff1*(rho(i+1,j,k)-rho(i  ,j,k)))
                  cff =0.5_r8*(pn(i,j-1)+pn(i,j  ))
                  cff1=0.5_r8*(pn(i,j  )+pn(i,j+1))
                  dRde_r=-0.5*g*orho0*
     &                     (cff*(rho(i,j  ,N)-rho(i,j-1,N))+         !y-direction gradient on RHO at k     
     &                     cff1*(rho(i,j+1,N)-rho(i,j  ,N)))
                  !--------------------------------------------------
                  dRdx0 = (dRdx0*(z_w(i,j,N)-z_w(i,j,k+1))+          
     &                     dRdx_r*(z_w(i,j,k+1)-z_w(i,j,k)))
                  cff=1.0_r8/(z_w(i,j,N)-z_w(i,j,k))
                  dRdx0=cff*dRdx0                                     !it is a depth averaged buoyancy gradient
                  dRde0 = (dRde0*(z_w(i,j,N)-z_w(i,j,k+1))+          
     &                     dRde_r*(z_r(i,j,k+1)-z_w(i,j,k)))
                  cff=1.0_r8/(z_w(i,j,N)-z_w(i,j,k))
                  dRde0=cff*dRde0                                     !it is a depth averaged buoyancy gradient     
                  
                  !---------------velocity difference----------------
                   dUdz0=0.5*(u(i  ,j,N,tind)-u(i  ,j,k,tind)+            !u difference bwteen surface and k-layer at RHO
     &                    u(i+1,j,N,tind)-u(i+1,j,k,tind))
                   dVdz0=0.5*(v(i,j  ,N,tind)-v(i,j  ,k,tind)+            !v difference bwteen surface and k-layer at RHO
     &                    v(i,j+1,N,tind)-v(i,j+1,k,tind))
                  
                  !---------------buoyancy difference---------------
                  cff=-g*orho0
                  dRdz0=cff*(rho(i  ,j  ,N)-rho(i  ,j  ,k))
                  
                  !----------once again, calculate PV----------------
                  pvor0=f(i,j)*((f(i,j)+rvor_r(i,j))*dRdz0            !surface vertical vorticity is used, according to Dong et al. (2020), vorticity varies slightly in the SI layer
     &                   +dUdz0*dRde0-dVdz0*dRdx0)
               END DO
              !--------------now, we get H_si-----------------------
               k_si=k+1;
               H_si=ABS(z_r(i,j,k)-z_w(i,j,N))                        !strictly, Hsi is thickness between surface and k-layer
               wstar_si=(B0*H_si)**(1.0_r8/3.0_r8)                    !convective velocity
               
               !------------get front direction----------------------
               thetaB=ATAN(dRde0/dRdx0)
               IF (dRdx0.lt.0.0 .and. dRde0.ge.0.0)       THEN            !second quadrant 
                   thetaB=pi+thetaB
               ELSE IF (dRdx0.lt.0.0 .and. dRde0.gt.0.0)  THEN            !third quadrant
                   thetaB=pi+thetaB
               ELSE IF (dRdx0.gt.0.0 .and. dRde0.lt.0.0)  THEN            !forth quadrant
                   thetaB=2.0_r8*pi+thetaB
               END IF
              !
                !----------get along front shear and its direction----
               dRdb_si=dRdx0*COS(thetaB)+dRde0*SIN(thetaB)
               U_si = dRdb_si*H_si/ABS(f(i,j))
               IF (f(i,j).le.0.0)                       THEN            !along-front flow direction
                   thetaU=thetaB-pi/2.0                               !in southern hemisphere, 90 clockwise
               ELSE
                   thetaU=thetaB+pi/2.0                              !in the northern, 90 anticlockwise
               END IF
               IF (thetaU.lt.0.0)                       THEN
                   thetaU=thetaU+2.0*pi                             !make sure thetaU is in [0 360]
               ELSE IF (thetaU.gt.2.0*pi) THEN
                   thetaU=thetaU-2.0*pi
               END IF
               
               !-------------get wind direction----------------------
               thetaW=ATAN(vstr0/ustr0)
               !print *, thetaU,thetaW
               IF (ustr0.lt.0.0 .and. vstr0.ge.0.0)      THEN           !second quadrant 
                   thetaW=pi+thetaW
               ELSE IF (ustr0.lt.0.0 .and. vstr0.gt.0.0) THEN           !third quadrant
                   thetaW=pi+thetaW
               ELSE IF (ustr0.gt.0.0 .and. vstr0.lt.0.0) THEN           !forth quadrant
                   thetaW=2.0_r8*pi+thetaW
               END IF
               IF (thetaU.lt.0.0)                       THEN
                   thetaU=thetaU+2.0*pi                             !make sure thetaU is in [0 360]
               ELSE IF (thetaU.gt.2.0*pi) THEN
                   thetaU=thetaU-2.0*pi
               END IF
               theta_si=thetaW-thetaU                                 !theta in the quartic equation
               IF (abs(theta_si).ge.pi/2.0) THEN
                   x_si=1.0
               ELSE
                  alpha_si=c_si*c_si*c_si*(wstar_si**3.0/ABS(U_si**3.0)+    
     &                  ustar_si**2.0/u_si**2.0*COS(theta_si))**2.0
      !            print *, alpha_si
                  IF (alpha_si.gt.1.0E6)                 THEN            !if alpha is too large
                      x_si=1.0_r8-alpha_si**(-1.0_r8/3.0_r8)
                  ELSE IF (alpha_si.lt.1.0E-6)           THEN            !if alpha is too small
                      x_si=alpha_si**(1.0_r8/4.0_r8)
                  ELSE
                     CALL quartic_poly(1.0,alpha_si,-3.0_r8*alpha_si,       
     &                               3.0_r8*alpha_si,-alpha_si,x_si)   ! solve quartic polynomail
                  END IF
               END IF
               !
               !-------To make simulation stable---------------------
               !-------we cannot let x_si too big--------------------
               !-------which implies buoyancy gradient---------------
               !-------is close to zero------------------------------
               !--------HERE, WE CHOOSE X_SI<0.95-------------------- !in Bachman et al. (2017), 0.9 is suggested, but according to Dong et al. (2020), an average percentage is 0.2 globally, so we can use a smaller threshold
               hc_si=H_si*x_si                                        !upper convective layer
               !print *,EBF0,B0
               !=====================================================
               !----------------now, we can get GSP------------------
               !------------Note that GSP is on W grid---------------
               IF (x_si.le.0.95)                    THEN               ! &             0.06*U_si/f(i,j).ge.55.0
               DO k=k_si,N
                   zw=-ABS(z_w(i,j,k)-z_w(i,j,N))                      !to make the surface layer as 0
                   IF (zw.gt.-H_si .and. zw.le.-hc_si)  THEN
                       GSP_si(k)=(EBF0+B0)*(zw+H_si)/H_si
                   ELSE IF (zw.gt.-hc_si .and. zw.lt.0.0) THEN
                       GSP_si(k)=(EBF0+B0)*(zw+H_si)/H_si-           
     &                           B0*(zw+hc_si)/hc_si
                   ELSE 
                       GSP_si(k)=0.0
                   END IF
                   
                   !-----------------tracer diffusivity--------------
                   IF (k.gt.k_si .and. k.lt.N)           THEN
                       
                       !balanced Ri
                       cff =0.5_r8*(pm(i-1,j)+pm(i  ,j))
                       cff1=0.5_r8*(pm(i  ,j)+pm(i+1,j))
                       !
                       dRdx_Ri=-0.25*g*orho0*
     &                     (cff*(rho(i  ,j,k  )-rho(i-1,j,k  ))+       !!x-direction gradient on RHO at k          
     &                     cff1*(rho(i+1,j,k  )-rho(i  ,j,k  ))+
     &                      cff*(rho(i  ,j,k+1)-rho(i-1,j,k+1))+       !!x-direction gradient on RHO at k+1          
     &                     cff1*(rho(i+1,j,k+1)-rho(i  ,j,k+1)))
                       cff =0.5_r8*(pn(i,j-1)+pn(i,j  ))
                       cff1=0.5_r8*(pn(i,j  )+pn(i,j+1))
                       dRde_Ri=-0.25*g*orho0*
     &                     (cff*(rho(i,j  ,k  )-rho(i,j-1,k  ))+         !y-direction gradient on RHO at k     
     &                     cff1*(rho(i,j+1,k  )-rho(i,j  ,k  ))+
     &                     cff*(rho(i,j  ,k+1)-rho(i,j-1,k+1))+         !y-direction gradient on RHO at k+1     
     &                     cff1*(rho(i,j+1,k+1)-rho(i,j  ,k+1)))
                       cff=f(i,j)*f(i,j)/(dRdx_Ri**2.0+dRde_Ri**2.0)
                       !print *,dRdx_Ri
                       Akv_si=cff*GSP_si(k)                           !Akv due to SI
                       !cff=f(i,j)*f(i,j)/(dRdx_Ri**2.0+dRde_Ri**2.0)
                       Rib=cff*bvf(i,j,k)                            !balanced Ri
                       !print *, k,Rib
                       Akt_si=2.0_r8*Akv_si/(1.0_r8+(10.0_r8*        
     &                                       MAX(0.0,Rib))**0.8)
                       !print *, GSP_si(k)
                       !----finally, we get vertical diffusivity-----
                       !-------and, replace Akv,Akt, and Aks---------
                       !IF (AKv_si.lt.0.1) THEN
                       AKv(i,j,k)=Akv_si
                       AKt(i,j,k,itemp)=Akt_si
#  ifdef SALINITY
                       AKt(i,j,k,isalt)=Akt_si
# endif                      
                       
                  END IF
               END DO
               
              !===================================================
               !===============isopycnal tracer diffusivity========
               DO k=k_si,N
                    cff =0.5_r8*(pm(i-1,j)+pm(i  ,j))
                    cff1=0.5_r8*(pm(i  ,j)+pm(i+1,j))
                    !
                    dRdx_Ri=-0.5*g*orho0*
     &                  (cff*(rho(i  ,j,k  )-rho(i-1,j,k  ))+       !!x-direction gradient on RHO at k          
     &                  cff1*(rho(i+1,j,k  )-rho(i  ,j,k  )))
                    cff =0.5_r8*(pn(i,j-1)+pn(i,j  ))
                    cff1=0.5_r8*(pn(i,j  )+pn(i,j+1))
                    dRde_Ri=-0.5*g*orho0*
     &                  (cff*(rho(i,j  ,k  )-rho(i,j-1,k  ))+         !y-direction gradient on RHO at k     
     &                  cff1*(rho(i,j+1,k  )-rho(i,j  ,k  )))
                    cff=f(i,j)*f(i,j)/(dRdx_Ri**2.0+dRde_Ri**2.0)
                    !
                    Rib=0.5*cff*(bvf(i,j,k)+bvf(i,j,k-1))              !balanced Ri               
                    !------isopycnal tracer diffusivity-----------
                    Kiso_si(i,j,k)=0.5_r8*(GSP_si(k)+GSP_si(k-1))* 
     &                                MIN(1.0_r8,Rib**2.0)/f(i,j)/f(i,j) 
              END DO
               
              !---------------vertical convection--------------------
              !
              DO k=k_si,N
                  zw=-ABS(z_w(i,j,k)-z_w(i,j,N)) 
                  IF (zw.gt.-hc_si .and. zw.lt.0.0)      THEN
                     wtemp_si(k)=-stflx(i,j,itemp)*(zw+hc_si)/hc_si  !surface heat flux
                     wsalt_si(k)=-stflx(i,j,isalt)*(zw+hc_si)/hc_si  !surface freshwater flux
                  ELSE
		     wtemp_si(k)=0.0
	             wsalt_si(k)=0.0
		  END IF
                  
              END DO
               
              !------------vertical temperature and salt fluxes-------
              DO k=k_si,N
                  Ftemp_si(i,j,k)=wtemp_si(k)-wtemp_si(k-1)       !flux difference at RHO
                  Fsalt_si(i,j,k)=wsalt_si(k)-wsalt_si(k-1)       !determines the change of tracer concentration, deltaTemp=Ftemp_si*dt/Hz 
              END DO
              END IF
              
              
          END IF
        END DO
      END DO
      !
       !=====================================================
      !====================================================
   

# define k0 0 
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Akv(istr-1,j,k)=Akv(istr,j,k)
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
#  ifdef SALINITY
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Akv(iend+1,j,k)=Akv(iend,j,k)
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
#  ifdef SALINITY
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
#  endif
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Akv(i,jstr-1,k)=Akv(i,jstr,k)
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
#  ifdef SALINITY
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
#  endif
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Akv(i,jend+1,k)=Akv(i,jend,k)
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
#  ifdef SALINITY
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
#  endif
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Akv(istr-1,jstr-1,k)=Akv(istr,jstr,k)
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
#  ifdef SALINITY
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
#  endif
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Akv(istr-1,jend+1,k)=Akv(istr,jend,k)
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
#  ifdef SALINITY
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Akv(iend+1,jstr-1,k)=Akv(iend,jstr,k)
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
#  ifdef SALINITY
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
#  ifdef SALINITY
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
#  endif
        enddo
      endif
#  endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (istr,iend,jstr,jend, Akv)
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                        Akt(START_2D_ARRAY,0,itemp))
#  ifdef SALINITY
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                        Akt(START_2D_ARRAY,0,isalt))
#  endif
# endif
!!=============================================================================

!!!==============================================================================
# define k0 1 
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Kiso_si(istr-1,j,k)=Kiso_si(istr,j,k)        
            Ftemp_si(istr-1,j,k)=Ftemp_si(istr,j,k)
#  ifdef SALINITY
            Fsalt_si(istr-1,j,k)=Fsalt_si(istr,j,k)
#  endif
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Kiso_si(iend+1,j,k)=Kiso_si(iend,j,k)
            Ftemp_si(iend+1,j,k)=Ftemp_si(iend,j,k)
#  ifdef SALINITY
            Fsalt_si(iend+1,j,k)=Fsalt_si(iend,j,k)
#  endif
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Kiso_si(i,jstr-1,k)=Kiso_si(i,jstr,k)
            Ftemp_si(i,jstr-1,k)=Ftemp_si(i,jstr,k)
#  ifdef SALINITY
            Fsalt_si(i,jstr-1,k)=Fsalt_si(i,jstr,k)
#  endif
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Kiso_si(i,jend+1,k)=Kiso_si(i,jend,k)
            Ftemp_si(i,jend+1,k)=Ftemp_si(i,jend,k)
#  ifdef SALINITY
            Fsalt_si(i,jend+1,k)=Fsalt_si(i,jend,k)
#  endif
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Kiso_si(istr-1,jstr-1,k)=Kiso_si(istr,jstr,k)
          Ftemp_si(istr-1,jstr-1,k)=Ftemp_si(istr,jstr,k)
#  ifdef SALINITY
          Fsalt_si(istr-1,jstr-1,k)=Fsalt_si(istr,jstr,k)
#  endif
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Kiso_si(istr-1,jend+1,k)=Kiso_si(istr,jend,k)
          Ftemp_si(istr-1,jend+1,k)=Ftemp_si(istr,jend,k)
#  ifdef SALINITY
          Fsalt_si(istr-1,jend+1,k)=Fsalt_si(istr,jend,k)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Kiso_si(iend+1,jstr-1,k)=Kiso_si(iend,jstr,k)
          Ftemp_si(iend+1,jstr-1,k)=Ftemp_si(iend,jstr,k)
#  ifdef SALINITY
          Fsalt_si(iend+1,jstr-1,k)=Fsalt_si(iend,jstr,k)
#  endif
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Kiso_si(iend+1,jend+1,k)=Kiso_si(iend,jend,k)
          Ftemp_si(iend+1,jend+1,k)=Ftemp_si(iend,jend,k)
#  ifdef SALINITY
          Fsalt_si(iend+1,jend+1,k)=Fsalt_si(iend,jend,k)
#  endif
        enddo
      endif
#  endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (istr,iend,jstr,jend, Kiso_si)
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                        Ftemp_si(START_2D_ARRAY,0))
#  ifdef SALINITY
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                        Fsalt_si(START_2D_ARRAY,0))
#  endif
# endif




       
           
      END subroutine si_scheme_tile
      !===========================END OF MAIN CODE===============================
      !----------------subroutine for solving quartic polynomial-----------------           
      !==========================================================================            
      SUBROUTINE quartic_poly(a_si, b_si, c_si, d_si, e_si,x_si)
    
         implicit none

         real, intent(in)   :: a_si, b_si, c_si, d_si, e_si
         real, intent(out)  :: x_si
         real               :: delta1, delta2, delta3, m1, m2, m3, m4
         real               :: m5, m6a
   
         delta1 = c_si*c_si - 3.0*b_si*d_si + 12.0*a_si*e_si
         delta2 = 2.0*c_si*c_si*c_si - 9.0*b_si*c_si*d_si +           
     &            27.0*a_si*d_si*d_si + 27.0*b_si*b_si*e_si -        
     &            72.0*a_si*c_si*e_si
         m1 = (delta2 + SQRT(-4.0*delta1*delta1*delta1 +             
     &                       delta2*delta2))**(1/3.0)
         delta3 = 2.0**(1/3.0)*delta1/3.0/a_si/m1 +                  
     &            m1/3.0/a_si/2.0**(1/3.0)
         m2 = b_si*b_si/4.0/a_si/a_si - 2.0*c_si/3.0/a_si
         m3 = -b_si*b_si*b_si/a_si/a_si/a_si +                       
     &        4.0*b_si*c_si/a_si/a_si - 8.0*d_si/a_si
         m4 = -b_si/4.0/a_si
         m5 = 0.5*SQRT(m2 + delta3)
         m6a = 0.5*SQRT(2.0*m2 - delta3 - m3/8.0/m5)
   
         x_si = m4 - m5 + m6a ! only this root is real and less than 1
   
   
   !===== other roots not used here====================
   !x_si = m4 - m5 - m6a
   !x_si = m4 + m5 - m6a
   !x_si = m4 + m5 + m6a
   !===========================
   
   
      END SUBROUTINE quartic_poly                  
