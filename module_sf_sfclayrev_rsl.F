!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_sfclayrev_rsl

 REAL    , PARAMETER ::  VCONVC=1.
 REAL    , PARAMETER ::  CZO=0.0185
 REAL    , PARAMETER ::  OZO=1.59E-5

 REAL,   DIMENSION(0:1000 ),SAVE :: psim_stab,psim_unstab,psih_stab,psih_unstab

CONTAINS

!-------------------------------------------------------------------
   SUBROUTINE SFCLAYREV_rsl(U3D,V3D,T3D,QV3D,P3D,dz8w,                    &
                     CP,G,ROVCP,R,XLV,PSFC,CHS,CHS2,CQS2,CPM,      &
                     ZNT,UST,PBLH,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                     FM,FH,                                        &
                     XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                     U10,V10,TH2,T2,Q2,                            &
                     GZ1OZ0,WSPD,BR,ISFFLX,DX,                     &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,EOMEG,STBOLT,                          &
                     P1000mb,                                      &
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
                     itimestep, &
                     LAI, VEGTYP, MMINLU,                          &
!                     psim_hat_jun, psih_hat_jun,                   &
!                     psix_jun, psiq_jun, &
                     wspd_jun, ust_jun,        &
                     hc_jun, lai_jun, lc_jun, m10_wind, m2_temp,   &
                     beta_jun, zpd_jun,                            &
                     CHS_ori, &
!                     jpsim1,jpsim2,jpsim3,jpsim4,jpsiq1,jpsiq2,jpsiq3,jpsiq4, &
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte,                    &
                     ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,scm_force_flux           )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!   Changes in V3.7 over water surfaces:
!          1. for ZNT/Cd, replacing constant OZO with 0.11*1.5E-5/UST(I)
!             the COARE 3.5 (Edson et al. 2013) formulation is also available
!          2. for VCONV, reducing magnitude by half
!          3. for Ck, replacing Carlson-Boland with COARE 3
!-------------------------------------------------------------------
!-- U3D         3D u-velocity interpolated to theta points (m/s)
!-- V3D         3D v-velocity interpolated to theta points (m/s)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- dz8w        dz between full levels (m)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G           acceleration due to gravity (m/s^2)
!-- ROVCP       R/CP
!-- R           gas constant for dry air (J/kg/K)
!-- XLV         latent heat of vaporization for water (J/kg)
!-- PSFC        surface pressure (Pa)
!-- ZNT         roughness length (m)
!-- UST         u* in similarity theory (m/s)
!-- USTM        u* in similarity theory (m/s) without vconv correction
!               used to couple with TKE scheme
!-- PBLH        PBL height from previous time (m)
!-- MAVAIL      surface moisture availability (between 0 and 1)
!-- ZOL         z/L height over Monin-Obukhov length
!-- MOL         T* (similarity theory) (K)
!-- REGIME      flag indicating PBL regime (stable, unstable, etc.)
!-- PSIM        similarity stability function for momentum
!-- PSIH        similarity stability function for heat
!-- FM          integrated stability function for momentum
!-- FH          integrated stability function for heat
!-- XLAND       land mask (1 for land, 2 for water)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          net upward latent heat flux at surface (W/m^2)
!-- TSK         surface temperature (K)
!-- FLHC        exchange coefficient for heat (W/m^2/K)
!-- FLQC        exchange coefficient for moisture (kg/m^2/s)
!-- CHS         heat/moisture exchange coefficient for LSM (m/s)
!-- QGH         lowest-level saturated mixing ratio
!-- QSFC        ground saturated mixing ratio
!-- U10         diagnostic 10m u wind
!-- V10         diagnostic 10m v wind
!-- TH2         diagnostic 2m theta (K)
!-- T2          diagnostic 2m temperature (K)
!-- Q2          diagnostic 2m mixing ratio (kg/kg)
!-- GZ1OZ0      log(z/z0) where z0 is roughness length
!-- WSPD        wind speed at lowest model level (m/s)
!-- BR          bulk Richardson number in surface layer
!-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
!-- DX          horizontal grid size (m)
!-- SVP1        constant for saturation vapor pressure (kPa)
!-- SVP2        constant for saturation vapor pressure (dimensionless)
!-- SVP3        constant for saturation vapor pressure (K)
!-- SVPT0       constant for saturation vapor pressure (K)
!-- EP1         constant for virtual temperature (R_v/R_d - 1) (dimensionless)
!-- EP2         constant for specific humidity calculation 
!               (R_d/R_v) (dimensionless)
!-- KARMAN      Von Karman constant
!-- EOMEG       angular velocity of earth's rotation (rad/s)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- isftcflx    =0, (Charnock and Carlson-Boland); =1, AHW Ck, Cd, =2 Garratt
!-- iz0tlnd     =0 Carlson-Boland, =1 Czil_new
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
      INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
!                                                               
      INTEGER,  INTENT(IN )   ::        ISFFLX
      REAL,     INTENT(IN )   ::        SVP1,SVP2,SVP3,SVPT0
      REAL,     INTENT(IN )   ::        EP1,EP2,KARMAN,EOMEG,STBOLT
      REAL,     INTENT(IN )   ::        P1000mb
!
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           dz8w
                                        
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           QV3D, &
                                                              P3D, &
                                                              T3D

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::             MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                              TSK
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )               ::                U10, &
                                                              V10, &
                                                              TH2, &
                                                               T2, &
                                                               Q2, &
                                                             QSFC

!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)               ::             REGIME, &
                                                              HFX, &
                                                              QFX, &
                                                               LH, &
                                                          MOL,RMOL
!m the following 5 are change to memory size
!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                 GZ1OZ0,WSPD,BR, &
                                                  PSIM,PSIH,FM,FH

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                            U3D, &
                                                              V3D
                                        
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::               PSFC

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                            ZNT, &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                             CHS2, &
                                                             CQS2, &
                                                              CHS

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                      FLHC,FLQC

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                                 &
                                                              QGH
                                    
      REAL,     INTENT(IN   )               ::   CP,G,ROVCP,R,XLV,DX
 
      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(OUT)     ::                  ck,cka,cd,cda

      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(INOUT)   ::                           USTM

      INTEGER,  OPTIONAL,  INTENT(IN )   ::     ISFTCFLX, IZ0TLND
      INTEGER,  OPTIONAL,  INTENT(IN )   ::     SCM_FORCE_FLUX
! LOCAL VARS

      REAL,     DIMENSION( its:ite ) ::                       U1D, &
                                                              V1D, &
                                                             QV1D, &
                                                              P1D, &
                                                              T1D

      REAL,     DIMENSION( its:ite ) ::                    dz8w1d

      INTEGER ::  I,J
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      REAL,DIMENSION(ims:ime, jms:jme),INTENT(INOUT) :: psim_hat_jun, psih_hat_jun
      REAL,DIMENSION(ims:ime, jms:jme),INTENT(INOUT) :: wspd_jun, ust_jun
      REAL,DIMENSION(ims:ime, jms:jme),INTENT(INOUT) :: hc_jun, lai_jun, lc_jun, m10_wind, m2_temp, beta_jun, zpd_jun, CHS_ori
!      REAL,DIMENSION(ims:ime, jms:jme),INTENT(INOUT) ::  jpsim1,jpsim2,jpsim3,jpsim4,jpsiq1,jpsiq2,jpsiq3,jpsiq4
!      REAL,DIMENSION(its:ite) ::  psim_hat_jun1d, psih_hat_jun1d
!      REAL,DIMENSION(ims:ime) ::  psim_hat_jun1d, psih_hat_jun1d
      REAL,DIMENSION(ims:ime, jms:jme),INTENT(IN) :: LAI
      INTEGER,DIMENSION(ims:ime, jms:jme),INTENT(IN) :: VEGTYP
      CHARACTER(LEN=*),INTENT(IN) :: MMINLU
!      REAL,DIMENSION(its:ite) :: LAI1D
      REAL,DIMENSION(ims:ime) :: LAI1D
!      INTEGER,DIMENSION(its:ite) :: VEGTYP1D
      INTEGER,DIMENSION(ims:ime) :: VEGTYP1D
      INTEGER,INTENT(IN) :: itimestep
!      REAL,DIMENSION(ims:ime, jms:jme),INTENT(INOUT) :: psix_jun, psiq_jun
!      REAL,DIMENSION(its:ite) :: psix_jun1d, psiq_jun1d, wspd_jun1d, ust_jun1d
!      REAL,DIMENSION(its:ite) :: hc_jun1d, lai_jun1d, lc_jun1d, m10_wind1d, m2_temp1d, beta_jun1d, zpd_jun1d, CHS_ori1d
!      REAL,DIMENSION(its:ite) :: jpsim11d,jpsim21d,jpsim31d,jpsim41d,jpsiq11d,jpsiq21d,jpsiq31d,jpsiq41d
!      REAL,DIMENSION(its:ite) :: psix_jun1d, psiq_jun1d, 
      REAL,DIMENSION(its:ite) :: wspd_jun1d, ust_jun1d
      REAL,DIMENSION(its:ite) :: lc_jun1d, lai_jun1d, m10_wind1d, m2_temp1d, beta_jun1d, zpd_jun1d, CHS_ori1d
      REAL,DIMENSION(ims:ime) :: hc_jun1d
!      REAL,DIMENSION(its:ite) :: jpsim11d,jpsim21d,jpsim31d,jpsim41d,jpsiq11d,jpsiq21d,jpsiq31d,jpsiq41d
!      REAL :: time_totalS, time_totalE
!call cpu_time(time_totalS)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      DO J=jts,jte
        DO i=its,ite
          dz8w1d(I) = dz8w(i,1,j)
        ENDDO
   
        DO i=its,ite
           U1D(i) =U3D(i,1,j)
           V1D(i) =V3D(i,1,j)
           QV1D(i)=QV3D(i,1,j)
           P1D(i) =P3D(i,1,j)
           T1D(i) =T3D(i,1,j)
        ENDDO

        !  Sending array starting locations of optional variables may cause
        !  troubles, so we explicitly change the call.

        CALL SFCLAYREV1D(J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,               &
                CP,G,ROVCP,R,XLV,PSFC(ims,j),CHS(ims,j),CHS2(ims,j),&
                CQS2(ims,j),CPM(ims,j),PBLH(ims,j), RMOL(ims,j),   &
                ZNT(ims,j),UST(ims,j),MAVAIL(ims,j),ZOL(ims,j),    &
                MOL(ims,j),REGIME(ims,j),PSIM(ims,j),PSIH(ims,j),  &
                FM(ims,j),FH(ims,j),                               &
                XLAND(ims,j),HFX(ims,j),QFX(ims,j),TSK(ims,j),     &
                U10(ims,j),V10(ims,j),TH2(ims,j),T2(ims,j),        &
                Q2(ims,j),FLHC(ims,j),FLQC(ims,j),QGH(ims,j),      &
                QSFC(ims,j),LH(ims,j),                             &
                GZ1OZ0(ims,j),WSPD(ims,j),BR(ims,j),ISFFLX,DX,     &
                SVP1,SVP2,SVP3,SVPT0,EP1,EP2,KARMAN,EOMEG,STBOLT,  &
                P1000mb,                                           &
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!@                LAI1D, VEGTYP1D, MMINLU,                           &
                LAI(ims,j), VEGTYP(ims,j), MMINLU,                           &
!                psim_hat_jun1d, psih_hat_jun1d,                    &
!                psix_jun1d, psiq_jun1d,& 
                wspd_jun1d, ust_jun1d,      &
!@                hc_jun1d, lai_jun1d, lc_jun1d, m10_wind1d, m2_temp1d, &
                hc_jun(ims,j), lai_jun1d, lc_jun1d, m10_wind1d, m2_temp1d, &
                beta_jun1d, zpd_jun1d,                             &
                CHS_ori1d, &
                itimestep,                                         &
!                jpsim11d,jpsim21d,jpsim31d,jpsim41d,jpsiq11d,jpsiq21d,jpsiq31d,jpsiq41d,&
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                ids,ide, jds,jde, kds,kde,                         &
                ims,ime, jms,jme, kms,kme,                         &
                its,ite, jts,jte, kts,kte                          &
#if ( EM_CORE == 1 )
                ,isftcflx,iz0tlnd,scm_force_flux,                               &
                USTM(ims,j),CK(ims,j),CKA(ims,j),                  &
                CD(ims,j),CDA(ims,j)                               &
#endif
                                                                   )
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        DO i=its,ite
!           psim_hat_jun(i,j) = psim_hat_jun1d(i)
!           psih_hat_jun(i,j) = psih_hat_jun1d(i)
!           psix_jun(i,j)     = psix_jun1d(i)
!           psiq_jun(i,j)     = psiq_jun1d(i)
           wspd_jun(i,j)     = wspd_jun1d(i)
           ust_jun(i,j)      = ust_jun1d(i)
!           hc_jun(i,j)       = hc_jun1d(i)
           lai_jun(i,j)      = lai_jun1d(i)
           lc_jun(i,j)       = lc_jun1d(i)
           m10_wind(i,j)       = m10_wind1d(i)
           m2_temp(i,j)       = m2_temp1d(i)
           beta_jun(i,j)     = beta_jun1d(i)
           zpd_jun(i,j)      = zpd_jun1d(i)
           CHS_ori(i,j)      = CHS_ori1d(i)
!           jpsim1(i,j)        = jpsim11d(i)
!           jpsim2(i,j)        = jpsim21d(i)
!           jpsim3(i,j)        = jpsim31d(i)
!           jpsim4(i,j)        = jpsim41d(i)
!           jpsiq1(i,j)        = jpsiq11d(i)
!           jpsiq2(i,j)        = jpsiq21d(i)
!           jpsiq3(i,j)        = jpsiq31d(i)
!           jpsiq4(i,j)        = jpsiq41d(i)
        ENDDO
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      ENDDO

   END SUBROUTINE SFCLAYREV_rsl


!-------------------------------------------------------------------
   SUBROUTINE SFCLAYREV1D(J,UX,VX,T1D,QV1D,P1D,dz8w1d,                &
                     CP,G,ROVCP,R,XLV,PSFCPA,CHS,CHS2,CQS2,CPM,PBLH,RMOL, &
                     ZNT,UST,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH,FM,FH,&
                     XLAND,HFX,QFX,TSK,                            &
                     U10,V10,TH2,T2,Q2,FLHC,FLQC,QGH,              &
                     QSFC,LH,GZ1OZ0,WSPD,BR,ISFFLX,DX,             &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,EOMEG,STBOLT,                          &
                     P1000mb,                                      &
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
                     LAI, VEGTYP, MMINLU,                          &
!                     psim_hat_jun, psih_hat_jun,                   &
!                     psix_jun, psiq_jun,& 
                     wspd_jun, ust_jun,        &
                     hc_jun, lai_jun, lc_jun, m10_wind, m2_temp,   &
                     beta_jun, zpd_jun,                            &
                     CHS_ori,         &
                     Itimestep,                                    &
!                     jpsim1,jpsim2,jpsim3,jpsim4,jpsiq1,jpsiq2,jpsiq3,jpsiq4,&
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte,                    &
                     isftcflx, iz0tlnd,scm_force_flux,                            &
                     ustm,ck,cka,cd,cda                            )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      REAL,     PARAMETER     ::        XKA=2.4E-5
      REAL,     PARAMETER     ::        PRT=1.

      INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte, &
                                        J
!                                                               
      INTEGER,  INTENT(IN )   ::        ISFFLX
      REAL,     INTENT(IN )   ::        SVP1,SVP2,SVP3,SVPT0
      REAL,     INTENT(IN )   ::        EP1,EP2,KARMAN,EOMEG,STBOLT
      REAL,     INTENT(IN )   ::        P1000mb

!
      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(IN   )               ::             MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                              TSK
!
      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(IN   )               ::             PSFCPA

      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(INOUT)               ::             REGIME, &
                                                              HFX, &
                                                              QFX, &
                                                         MOL,RMOL
!m the following 5 are changed to memory size---
!
      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(INOUT)   ::                 GZ1OZ0,WSPD,BR, &
                                                  PSIM,PSIH,FM,FH

      REAL,     DIMENSION( ims:ime )                             , &
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
                INTENT(INOUT)   ::                             &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                             CHS2, &
                                                             CQS2, &
                                                              CHS
!      REAL,     DIMENSION( ims:ime ), INTENT(OUT)   ::        ZNT
      REAL,     DIMENSION( ims:ime ), INTENT(INOUT)   ::        ZNT
!      REAL,     DIMENSION( its:ite ), INTENT(INOUT)   ::       CHS_ori
      REAL,     DIMENSION( its:ite ), INTENT(INOUT)   ::       CHS_ori
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(INOUT)   ::                      FLHC,FLQC

      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(INOUT)   ::                                 &
                                                              QGH

      REAL,     DIMENSION( ims:ime )                             , &
                INTENT(OUT)     ::                        U10,V10, &
                                                TH2,T2,Q2,QSFC,LH

                                    
      REAL,     INTENT(IN   )               ::   CP,G,ROVCP,R,XLV,DX

! MODULE-LOCAL VARIABLES, DEFINED IN SUBROUTINE SFCLAY
      REAL,     DIMENSION( its:ite ),  INTENT(IN   )   ::  dz8w1d

      REAL,     DIMENSION( its:ite ),  INTENT(IN   )   ::      UX, &
                                                               VX, &
                                                             QV1D, &
                                                              P1D, &
                                                              T1D
 
      REAL, OPTIONAL, DIMENSION( ims:ime )                       , &
                INTENT(OUT)     ::                  ck,cka,cd,cda
      REAL, OPTIONAL, DIMENSION( ims:ime )                       , &
                INTENT(INOUT)   ::                           USTM

      INTEGER,  OPTIONAL,  INTENT(IN )   ::     ISFTCFLX, IZ0TLND
      INTEGER,  OPTIONAL,  INTENT(IN )   ::     SCM_FORCE_FLUX

! LOCAL VARS

      REAL,     DIMENSION( its:ite )        ::                 ZA, &
                                                        THVX,ZQKL, &
                                                           ZQKLP1, &
                                                           THX,QX, &
                                                            PSIH2, &
                                                            PSIM2, &
                                                           PSIH10, &
                                                           PSIM10, &
                                                           DENOMQ, &
                                                          DENOMQ2, &
                                                          DENOMT2, &
                                                            WSPDI, &
                                                           GZ2OZ0, &
                                                           GZ10OZ0
!
      REAL,     DIMENSION( its:ite )        ::                     &
                                                      RHOX,GOVRTH, &
                                                            TGDSA
!
      REAL,     DIMENSION( its:ite)         ::          SCR3,SCR4
      REAL,     DIMENSION( its:ite )        ::         THGB, PSFC
!
      INTEGER                               ::                 KL

      INTEGER ::  N,I,K,KK,L,NZOL,NK,NZOL2,NZOL10

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      REAL    ::  PL,THCON,TVCON,E1
      REAL    ::  THCON,E1
      REAL,DIMENSION( its:ite ) ::  PL(i), TVCON(i)
!!      REAL    ::  ZL,TSKV,DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
!      REAL    ::  TSKV,DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
      REAL    ::  DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
      REAL,DIMENSION( its:ite ) :: TSKV,DTHVDZ
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      REAL    ::  DTG,PSIX,DTTHX,PSIX10,PSIT,PSIT2,PSIQ,PSIQ2,PSIQ10
      REAL    ::  FLUXC,VSGD,Z0Q,VISC,RESTAR,CZIL,GZ0OZQ,GZ0OZT
      REAL    ::  ZW, ZN1, ZN2
!
! .... paj ...
!
!      REAL    :: zolzz,zol0
!     REAL    :: zolri,zolri2
!     REAL    :: psih_stable,psim_stable,psih_unstable,psim_unstable
!     REAL    :: psih_stable_full,psim_stable_full,psih_unstable_full,psim_unstable_full
      REAL    :: zl2,zl10,z0t
      REAL,     DIMENSION( its:ite )        ::         pq,pq2,pq10

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      REAL,DIMENSION( its:ite ), INTENT(INOUT) :: psim_hat_jun, psih_hat_jun
!@      REAL,DIMENSION( its:ite ), INTENT(IN) :: LAI
      REAL,DIMENSION( ims:ime ), INTENT(IN) :: LAI
!@      INTEGER,DIMENSION( its:ite ), INTENT(IN) :: VEGTYP
      INTEGER,DIMENSION( ims:ime ), INTENT(IN) :: VEGTYP
      CHARACTER(LEN=*),INTENT(IN) :: MMINLU
      REAL, DIMENSION( its:ite ) :: GZZ0
!      REAL,DIMENSION( its:ite ), INTENT(OUT)   :: psix_jun, psiq_jun,
      REAL,DIMENSION( its:ite ), INTENT(OUT)   :: wspd_jun, ust_jun
      REAL,DIMENSION( its:ite ), INTENT(OUT)   :: lai_jun, lc_jun, m10_wind, m2_temp, beta_jun, zpd_jun 
!@      REAL,DIMENSION( its:ite ), INTENT(INOUT) :: hc_jun
      REAL,DIMENSION( ims:ime ), INTENT(IN) :: hc_jun
!      REAL,DIMENSION( its:ite ), INTENT(OUT)   :: jpsim1,jpsim2,jpsim3,jpsim4,jpsiq1,jpsiq2,jpsiq3,jpsiq4
      integer,intent(in) :: itimestep
      REAL,dimension( its:ite ) :: PSIM_hat, PSIH_hat
      REAL,dimension( its:ite ) :: zl, zol0, zzzol, zolzz, zzzolzz

      REAL,dimension( its:ite )  ::  Hc, VAI, D1, Lc, Prc, ff, Uhc, Thc, &
!dwind      & znt_new, znt_old, zpd_dt, dtemp, PSIH_hat_ori, pq_chs, chs_i1, chs_i2!, chs_i3, chs_i4
      & znt_new, znt_old, zpd_dt, dtemp, dwind, PSIH_hat_ori, pq_chs, chs_i1, chs_i2!, chs_i3, chs_i4
!chs_temp
      real,dimension( its:ite) :: chs_temp
!chs_temp
      REAL,dimension( its:ite )  :: beta, beta_old, betaHF, betaNO
      REAL,dimension( its:ite, -50:50 ) :: Uwind!, Ustar
!      REAL,dimension( its:ite, -2000:2000 ) :: Uwind1, Ustar1
!      REAL,dimension( its:ite, -2000:2000 ) :: Uwind2, Ustar2
      REAL,dimension( its:ite, -50:50 ) :: Twind!, Tstar
!      REAL,dimension( its:ite, -2000:2000 ) :: Twind1, Tstar1
!      REAL,dimension( its:ite, -2000:2000 ) :: Twind2, Tstar2

      integer :: zzz
      integer :: ITER, ITER_CONFIG
      REAL :: factorm1, factorh1
      REAL :: factorm11, factorm12, factorh11, factorh12
      REAL :: zll
      real,dimension( its:ite ) :: Psim_temp, Psih_temp
!chs_temp      real,dimension( its:ite ) :: Psim_hat_temp, Psih_hat_temp
!chs_temp
      real,dimension( its:ite ) :: Psim_hat_temp, Psih_hat_temp, psih_hat_chs
!chs_temp

integer :: ierr
real, dimension(1:50  ) ::  SNUPTBL, RSTBL, RGLTBL, HSTBL,                &
                            SHDTBL, MAXALB,                               &
                            EMISSMINTBL, EMISSMAXTBL,                     &
                            LAIMINTBL, LAIMAXTBL,                         &
                            Z0MINTBL, Z0MAXTBL,                           &
                            ALBEDOMINTBL, ALBEDOMAXTBL,                   &
                            ZTOPVTBL,ZBOTVTBL
INTEGER, dimension(1:50):: NROTBL
INTEGER :: LUMATCH, IINDEX, LUCATS, LCC
CHARACTER(LEN=256) :: LUTYPE
real :: G1, G2, startz, endz, startz1, endz1, zhat
real,dimension(its:ite) :: cm1, ch1, cm2, ch2
LOGICAL,parameter :: ONOFF = .T.  ! profile on? off?


 REAL :: P,pp
 INTEGER :: pr_VEGTYP
 REAL    :: pr_znt, pr_zl, pr_zol, pr_br, pr_lai,  pr_ust1, pr_ust2, pr_l
 REAL    :: pr_psim_hat, pr_psih_hat, pr_psim, pr_psim1, pr_psim2, pr_psih, pr_psih1, pr_psih2
 REAL    :: pr_pq, pr_pq1, pr_pq2, pr_gz1Oz0, pr_gzz0
 REAL    :: pr_psix, pr_psit, pr_psiq, pr_lnpq, pr_wspd, pr_lnpq2, pr_psiq2
 REAL    :: pr_HC, pr_VAI, pr_PRC, pr_FF, pr_LC, pr_BETA, pr_CHS, Pr_zpd, pr_phihc, pr_zpd2, pr_mol2, pr_ITER
 REAL    :: pr_term1, pr_term2, pr_term3, pr_term4

 REAL,dimension( its:ite )  :: beta1, beta2

 REAL, PARAMETER :: betaN = 0.374

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MPI191107
REAL :: TIMES, TIMEE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



     ITER_CONFIG = 1
     Uwind = 0!; Ustar = 0
!     Uwind1 = 0; Ustar1 = 0
!     Uwind2 = 0; Ustar2 = 0
     Twind = 0!; Tstar = 0
!     Twind1 = 0; Tstar1 = 0
!     Twind2 = 0; Tstar2 = 0
     zpd_jun = 0
     beta_jun = 0

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!-------------------------------------------------------------------
      KL=kte

!IF ( J == 93 ) CALL cpu_time(TIMES)

      DO i=its,ite
! PSFC cb
         PSFC(I)=PSFCPA(I)/1000.
      ENDDO
!                                                      
!----CONVERT GROUND TEMPERATURE TO POTENTIAL TEMPERATURE:  
!                                                            
      DO 5 I=its,ite                                   
        TGDSA(I)=TSK(I)                                    
! PSFC cb
!        THGB(I)=TSK(I)*(100./PSFC(I))**ROVCP                
        THGB(I)=TSK(I)*(P1000mb/PSFCPA(I))**ROVCP   
    5 CONTINUE                                               
!                                                            
!-----DECOUPLE FLUX-FORM VARIABLES TO GIVE U,V,T,THETA,THETA-VIR.,
!     T-VIR., QV, AND QC AT CROSS POINTS AND AT KTAU-1.  
!                                                                 
!     *** NOTE ***                                           
!         THE BOUNDARY WINDS MAY NOT BE ADEQUATELY AFFECTED BY FRICTION,         
!         SO USE ONLY INTERIOR VALUES OF UX AND VX TO CALCULATE 
!         TENDENCIES.                             
!                                                           
   10 CONTINUE                                                     

!     DO 24 I=its,ite
!        UX(I)=U1D(I)
!        VX(I)=V1D(I)
!  24 CONTINUE                                             
                                                             
   26 CONTINUE                                               
                                                   
!.....SCR3(I,K) STORE TEMPERATURE,                           
!     SCR4(I,K) STORE VIRTUAL TEMPERATURE.                                       
                                                                                 
      DO 30 I=its,ite
! PL cb
         PL(i)=P1D(I)/1000.
         SCR3(I)=T1D(I)                                                   
!         THCON=(100./PL)**ROVCP                                                 
         THCON=(P1000mb*0.001/PL(i))**ROVCP
         THX(I)=SCR3(I)*THCON                                               
         SCR4(I)=SCR3(I)                                                    
         THVX(I)=THX(I)                                                     
         QX(I)=0.                                                             
   30 CONTINUE                                                                 
!                                                                                
      DO I=its,ite
         QGH(I)=0.                                                                
         FLHC(I)=0.                                                               
         FLQC(I)=0.                                                               
         CPM(I)=CP                                                                
      ENDDO
!                                                                                
!     IF(IDRY.EQ.1)GOTO 80                                                   
      DO 50 I=its,ite
         QX(I)=QV1D(I)                                                    
         TVCON(i)=(1.+EP1*QX(I))                                      
         THVX(I)=THX(I)*TVCON(i)
         SCR4(I)=SCR3(I)*TVCON(i) 
!         print*,"THETA", itimestep, i,j, THVX(I)*100, THX(I)*100, TVCON(i), EP1, QX(I)*100, QV1D(I)*100
   50 CONTINUE                                                                 
!                                                                                
      DO 60 I=its,ite
        E1=SVP1*EXP(SVP2*(TGDSA(I)-SVPT0)/(TGDSA(I)-SVP3))                       
!  for land points QSFC can come from previous time step
        if(xland(i).gt.1.5.or.qsfc(i).le.0.0)QSFC(I)=EP2*E1/(PSFC(I)-E1)                                                 
! QGH CHANGED TO USE LOWEST-LEVEL AIR TEMP CONSISTENT WITH MYJSFC CHANGE
! Q2SAT = QGH IN LSM
        E1=SVP1*EXP(SVP2*(T1D(I)-SVPT0)/(T1D(I)-SVP3))                       
        PL(i)=P1D(I)/1000.
        QGH(I)=EP2*E1/(PL(i)-E1)                                                 
        CPM(I)=CP*(1.+0.8*QX(I))                                   
   60 CONTINUE                                                                   
   80 CONTINUE
                                                                                 
!-----COMPUTE THE HEIGHT OF FULL- AND HALF-SIGMA LEVELS ABOVE GROUND             
!     LEVEL, AND THE LAYER THICKNESSES.                                          
                                                                                 
      DO 90 I=its,ite
        ZQKLP1(I)=0.
        RHOX(I)=PSFC(I)*1000./(R*SCR4(I))                                       
   90 CONTINUE                                                                   
!                                                                                
      DO 110 I=its,ite                                                   
           ZQKL(I)=dz8w1d(I)+ZQKLP1(I)
  110 CONTINUE                                                                 
!                                                                                
      DO 120 I=its,ite
         ZA(I)=0.5*(ZQKL(I)+ZQKLP1(I))                                        
! print*, "za test", itimestep, i,j, za(i), ZQKL(I), ZQKLP1(I), SCR4(I), PSFC(I), R, dz8w1d(I) 
  120 CONTINUE                                                                 
!                                                                                
      DO 160 I=its,ite
        GOVRTH(I)=G/THX(I)                                                    
  160 CONTINUE                                                                   
                                                                                 
!-----CALCULATE BULK RICHARDSON NO. OF SURFACE LAYER, ACCORDING TO               
!     AKB(1976), EQ(12).                                                         
                   
      DO 260 I=its,ite
        GZ1OZ0(I)=ALOG((ZA(I)+ZNT(I))/ZNT(I))   ! log((z+z0)/z0)                                     
        GZ2OZ0(I)=ALOG((2.+ZNT(I))/ZNT(I))      ! log((2+z0)/z0)                           
        GZ10OZ0(I)=ALOG((10.+ZNT(I))/ZNT(I))    ! log((10+z0)z0)                    
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        GZZ0  (I)=ALOG((ZA(I))/ZNT(I))
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        IF((XLAND(I)-1.5).GE.0)THEN                                            
          zl(i)=ZNT(I)                                                            
        ELSE                                                                     
          zl(i)=0.01                                                                
        ENDIF                                                                    
        WSPD(I)=SQRT(UX(I)*UX(I)+VX(I)*VX(I))                        

        TSKV(i)=THGB(I)*(1.+EP1*QSFC(I))                     
        DTHVDZ(i)=(THVX(I)-TSKV(i))
! write(*,'(A4,2x,2(I3,2x),5(f15.5,2x))') "TEMP", itimestep, i, j, DTHVDZ, (THVX(I)-TSKV), THVX(I), TSKV, hfx(i), qfx(i)
! write(*,'(4A,2x,3(I3,2x),6(f15.5,2x))') "TEMP", itimestep, i, j, DTHVDZ, (THVX(I)-TSKV), THVX(I), TSKV, hfx(i), qfx(i)
! write(*,*) "TEMP", itimestep, i, j, DTHVDZ(i), (THVX(I)-TSKV(i)), THVX(I), TSKV(i), hfx(i), qfx(i)
! print*,"TEMP1", EP1, QSFC(i), THGB(i), P1000mb, PL(i), ROVCP, SCR3(I), THCON, THX(I) 

!  Convective velocity scale Vc and subgrid-scale velocity Vsg
!  following Beljaars (1994, QJRMS) and Mahrt and Sun (1995, MWR)
!                                ... HONG Aug. 2001
!
!       VCONV = 0.25*sqrt(g/tskv*pblh(i)*dthvm)
!      Use Beljaars over land, old MM5 (Wyngaard) formula over water
        if (xland(i).lt.1.5) then
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        fluxc = max(hfx(i)/rhox(i)/cp                    &
              + ep1*tskv(i)*qfx(i)/rhox(i),0.)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        VCONV = vconvc*(g/tgdsa(i)*pblh(i)*fluxc)**.33
        else
        IF(-DTHVDZ(i).GE.0)THEN
          DTHVM=-DTHVDZ(i)
        ELSE
          DTHVM=0.
        ENDIF
!       VCONV = 2.*SQRT(DTHVM)
! V3.7: reducing contribution in calm conditions
        VCONV = SQRT(DTHVM)
        endif
! Mahrt and Sun low-res correction
        VSGD = 0.32 * (max(dx/5000.-1.,0.))**.33
        WSPD(I)=SQRT(WSPD(I)*WSPD(I)+VCONV*VCONV+vsgd*vsgd)
        WSPD(I)=AMAX1(WSPD(I),0.1)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        WSPD_JUN(I) = WSPD(I)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        BR(I)=GOVRTH(I)*ZA(I)*DTHVDZ(i)/(WSPD(I)*WSPD(I))                        

!  IF PREVIOUSLY UNSTABLE, DO NOT LET INTO REGIMES 1 AND 2
        IF(MOL(I).LT.0.)BR(I)=AMIN1(BR(I),0.0)
!jdf
        RMOL(I)=-GOVRTH(I)*DTHVDZ(i)*ZA(I)*KARMAN
!jdf

  260 CONTINUE                                                                   

!                                                                                
!-----DIAGNOSE BASIC PARAMETERS FOR THE APPROPRIATED STABILITY CLASS:            
!                                                                                
!                                                                                
!     THE STABILITY CLASSES ARE DETERMINED BY BR (BULK RICHARDSON NO.)           
!     AND HOL (HEIGHT OF PBL/MONIN-OBUKHOV LENGTH).                              
!                                                                                
!     CRITERIA FOR THE CLASSES ARE AS FOLLOWS:                                   
!                                                                                
!        1. BR .GE. 0.0;                                                         
!               REPRESENTS NIGHTTIME STABLE CONDITIONS (REGIME=1),               
!                                                                                
!        3. BR .EQ. 0.0                                                          
!               REPRESENTS FORCED CONVECTION CONDITIONS (REGIME=3),              
!                                                                                
!        4. BR .LT. 0.0                                                          
!               REPRESENTS FREE CONVECTION CONDITIONS (REGIME=4).                
!                                                                                
!CCCCC                                                                           

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      DO 320 I=its,ite
DO I=its,ite ! LOOP1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                                                                           

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
! define constant variables

  VAI(i) = LAI(i)
  IF ( VAI(I) <= 0 ) then
    VAI(i) = 0.1
  ENDIF
  Lc(i) = 4.0 * Hc_jun(i) / VAI(i)

  IF ( Lc(i) >= 500. ) Lc(i) = 500.
  lai_jun(i) = VAI(i)
  lc_jun(i) = Lc(i)
  Hc(i) = Hc_jun(i)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      if (br(I).gt.0) then
!        if (br(I).gt.250.0) then
!        zol(I)=zolri(250.0,ZA(I),ZNT(I),)
!        else
!        zol(I)=zolri(br(I),ZA(I),ZNT(I))
!        endif
!      endif
!!
!      if (br(I).lt.0) then
!       IF(UST(I).LT.0.001)THEN
!          ZOL(I)=BR(I)*GZ1OZ0(I)
!        ELSE
!        if (br(I).lt.-250.0) then
!        zol(I)=zolri(-250.0,ZA(I),ZNT(I))
!        else
!        zol(I)=zolri(br(I),ZA(I),ZNT(I))
!        endif
!       ENDIF
!      endif

    IF ( (XLAND(I)-1.5) .GT. 0. .or.(VEGTYP(I)==1.or.VEGTYP(I)>=16) ) then
      if (br(I).gt.0) then
        if (br(I).gt.250.0) then
        zol(I)=zolri(250.0,ZA(I),ZNT(I))
        else
        zol(I)=zolri(br(I),ZA(I),ZNT(I))
        endif
      endif
!
      if (br(I).lt.0) then
       IF(UST(I).LT.0.001)THEN
          ZOL(I)=BR(I)*GZ1OZ0(I)
        ELSE
        if (br(I).lt.-250.0) then
        zol(I)=zolri(-250.0,ZA(I),ZNT(I))
        else
        zol(I)=zolri(br(I),ZA(I),ZNT(I))
        endif
       ENDIF
      endif
    ELSE
      if (br(I).gt.0) then
        if (br(I).gt.250.0) then
        zol(I)=zolri3(250.0,ZA(I),ZNT(I),0.01,Lc(i),Hc(i),ust(i),xka)
        else
        zol(I)=zolri3(br(I),ZA(I),ZNT(I),0.01,Lc(i),Hc(i),ust(i),xka)
        endif
      endif
!
      if (br(I).lt.0) then
       IF(UST(I).LT.0.001)THEN
          ZOL(I)=BR(I)*GZZ0(I)
!          ZOL(I)=BR(I)*GZ1OZ0(I)
        ELSE
        if (br(I).lt.-250.0) then
        zol(I)=zolri3(-250.0,ZA(I),ZNT(I),0.01,Lc(i),Hc(i),ust(i),xka)
        else
        zol(I)=zolri3(br(I),ZA(I),ZNT(I),0.01,Lc(i),Hc(i),ust(i),xka)
        endif
       ENDIF
      endif
!      IF (zol(i) > 500. ) zol(i)=500.
    ENDIF

    IF (BR(I) == 0) zol(i) = 0.

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
ENDDO !LOOP1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
znt_new(:)=znt(:)
znt_old(:)=0
DO I = its,ite !LOOP2
zpd_dt(i) = Hc(i) * 0.35

IF (ITER_CONFIG .EQ. 1 .and.  XLAND(I) .LT. 1.5 .and.(VEGTYP(I)>=2.and.VEGTYP(I)<=15)) THEN !LANDUSE_only_for_RSL

ITER = 1


IF (zol(i)==0.) THEN
 MOL(i) = 999999.
ELSE
 MOL(i) = za(i) / zol(i)
ENDIF
beta(i) = betaN /PHIM(zpd_dt(i),MOL(i))

betaHF(i) = betaN /PHIM(zpd_dt(i),MOL(i))
betaNO(i) = 0.4/2./PHIM(zpd_dt(i),MOL(i))
beta(i)   = betaNO(i) + (BetaHF(i)-BetaNO(i))/(1+2* ((abs(Lc(i)/MOL(i)+0.15))**1.5)  )
IF ( Lc(i)/MOL(i) > -0.15 ) beta(i) = betaHF(i)

!IF ( beta(i) <= 0.01 ) beta(i) = 0.01
beta_old(i) = beta(i)

ENDIF !LANDUSE_only_for_RSL


IF (ITER_CONFIG .EQ. 1 .and.  XLAND(I) .LT. 1.5 .and.(VEGTYP(I)>=2.and.VEGTYP(I)<=15)) THEN !LANDUSE_only_for_RSL

DO WHILE  ( (abs (beta1(i) - beta(i)) > 0.0001 .and. ITER < 5 ).or. ITER == 1 )

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!update znt

  zpd_dt(i) = beta_old(i)*beta_old(i)*Lc(i) ! dt in HF07

 IF ( zpd_dt(i) .ge. 0.90*za(i) ) THEN
   zpd_dt(i) = 0.90*za(i)
 ENDIF
 IF ( zpd_dt(i) .ge. 0.90*Hc(i) ) THEN
   zpd_dt(i) = 0.90*Hc(i)
 ENDIF

IF ( ITER >= 2 ) beta2(i) = beta1(i)
beta1(i) = beta(i)

beta(i) = betaN /PHIM(zpd_dt(i),MOL(i))

betaHF(i) = betaN /PHIM(zpd_dt(i),MOL(i))
betaNO(i) = 0.4/2./PHIM(zpd_dt(i),MOL(i))
beta(i)   = betaNO(i) + (BetaHF(i)-BetaNO(i))/(1+2* ((abs(Lc(i)/MOL(i)+0.15))**1.5)  )
IF ( Lc(i)/MOL(i) > -0.15 ) beta(i) = betaHF(i)

beta_old(i) = beta(i)

IF ( abs(beta2(i)-beta(i)) < 0.001 .and. abs(beta1(I)-beta(i)) > 0.001 .and. ITER >= 2 ) then
beta(i) = ( beta(i)+beta1(i) )/ 2.
beta_old(i) = beta(i)
go to 987
ENDIF
IF ( beta(i) < beta1(i) .and. beta(i) < beta2(i) .and. ITER >= 2 ) then
beta(i) = ( beta(i)+beta1(i) )/ 2.
beta_old(i) = beta(i)
go to 987
ENDIF
IF ( beta(i) > beta1(i) .and. beta(i) > beta2(i) .and. ITER >= 2) then
beta(i) = ( beta(i)+beta1(i) )/ 2.
beta_old(i) = beta(i)
go to 987
ENDIF

987 continue

beta_old(i) = beta(i)


ITER = ITER + 1
ENDDO  ! WHILE LOOP

beta_jun(I) = beta_old(i)
zpd_jun(I) = zpd_dt(i)

Prc(I) = 0.5 + 0.3 * tanh ( 2. * Lc(I) * zol(I) / za(I) )
ff(I)     = 0.5 * ( sqrt(1. + 4.*0.1*Prc(I) ) -1. )
IF ( BR(i) > 0 ) THEN
  cm2(i) = KARMAN        * ( 2.       -2.*beta_old(i)*beta_old(i)*Lc(i)*5./MOL(i)/PHIM(zpd_dt(i),MOL(i)) ) / (2.*beta_old(i)*PHIM(zpd_dt(i),MOL(i)) - KARMAN         )
  ch2(i) = KARMAN *Prc(i)* ( 2. +ff(i)-2.*beta_old(i)*beta_old(i)*Lc(i)*5./MOL(i)/PHIM(zpd_dt(i),MOL(i)) ) / (2.*beta_old(i)*PHIM(zpd_dt(i),MOL(i)) - KARMAN *Prc(i) )

ELSE IF ( BR(i) < 0 ) THEN
  cm2(i) = KARMAN        * (2.       - 8.*beta_old(i)*beta_old(i)*Lc(i)/MOL(i)*(PHIM(zpd_dt(i),MOL(i))**4))/(2.*beta_old(i)*PHIM(zpd_dt(i),MOL(i)) - KARMAN         )
  ch2(i) = KARMAN *Prc(i)* (2. +ff(i)- 8.*beta_old(i)*beta_old(i)*Lc(i)/MOL(i)*(PHIM(zpd_dt(i),MOL(i))**4))/(2.*beta_old(i)*PHIM(zpd_dt(i),MOL(i)) - KARMAN *Prc(i) )
ELSE
  cm2(i) = KARMAN        * (2.        )/(2*beta_old(i) - KARMAN         )
  ch2(i) = KARMAN *Prc(i)* (2. +ff(i) )/(2*beta_old(i) - KARMAN *Prc(i) )
ENDIF
IF ( cm2(i) >= 5. ) cm2(i) = 5.
IF ( cm2(i) <= 0.  ) cm2(i) = 0.01
IF ( ch2(i) >= 5. ) ch2(i) = 5.
IF ( ch2(i) <= 0.  ) ch2(i) = 0.01


ITER = 1
DO WHILE  ( (abs (znt_old(i) - znt(i)) > 0.001 .and. ITER < 5 ).or. ITER == 1 )

  zzzol(I) = zol(i) / za(i) * zpd_dt(i)
  zol0(i)=zol(I)*znt(I)/za(I)          ! z0/L
  zzzolzz(I) = zzzol(I)+zol0(I)


IF (zzzol(i)==0.) THEN
 MOL(i) = 999999.
ELSE
 MOL(i) = zpd_dt(i) / zzzol(i)
ENDIF

cm1(i) = (1.-(0.4/(betaN*2.)))*exp(cm2(i)/2.)

startz = zpd_dt(i)
endz   = za(i)*2.

PSIM_hat_temp(i) = 0
startz1 = startz
endz1   = za(i)
G1  = (endz1-startz1)/2.*(-1/sqrt(3.)) + (endz1+startz1)/2.
G2  = (endz1-startz1)/2.*( 1/sqrt(3.)) + (endz1+startz1)/2.
PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * (     (PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )

DO zhat = za(i), endz-za(i)/2., za(i)/2.
startz1 = zhat
endz1 = zhat+za(i)/2.
G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.

PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * (     (PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
ENDDO


  IF ( BR(i) > 0 ) THEN
    psim_temp(i) = psim_stable(zzzol(I))-psim_stable(zol0(i))
  ELSE IF ( BR(i) < 0 ) THEN
    psim_temp(i) = psim_unstable(zzzol(I))-psim_unstable(zol0(i))
  ELSE
    psim_temp(i) = 0
  ENDIF

  ZNT_OLD(i) = ZNT(i)
  ZNT(i) = zpd_dt(i) * exp ( -1.*0.4 / beta_old(i) ) * exp ( -1.*psim_temp(i) ) * exp( PSIM_hat_temp(i) )

 IF ( ZNT(i) > 3./7. * zpd_dt(i) ) THEN
   ZNT(i) = 3./7. * zpd_dt(i)
 ELSE IF ( ZNT(i) < 0.0001 ) THEN
   ZNT(i) = 0.0001
 ENDIF


ITER = ITER + 1
ENDDO  ! WHILE LOOP
ENDIF !LANDUSE_only_for_RSL
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ... paj: compute integrated similarity functions.
!
        zolzz(i)=zol(I)*(za(I)+znt(I))/za(I) ! (z+z0/L
        zol10=zol(I)*(10.+znt(I))/za(I)   ! (10+z0)/L
        zol2=zol(I)*(2.+znt(I))/za(I)     ! (2+z0)/L
        zol0(i)=zol(I)*znt(I)/za(I)          ! z0/L
        ZL2=(2.)/ZA(I)*ZOL(I)             ! 2/L      
        ZL10=(10.)/ZA(I)*ZOL(I)           ! 10/L
        GZ1OZ0(I)=ALOG((ZA(I)+ZNT(I))/ZNT(I))   ! log((z+z0)/z0)                                     

        IF((XLAND(I)-1.5).LT.0.)THEN
        ZLL=(0.01)/ZA(I)*ZOL(I)   ! (0.01)/L     
        ELSE
        ZLL=zol0(i)                     ! z0/L
        ENDIF


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        IF((XLAND(I)-1.5).LT.0  .and.(VEGTYP(I)>=2.and.VEGTYP(I)<=15))THEN

          call RPSI_hat (zl(i), zol0(i), ZOL(I), ZA(i), ZA(I), Lc(I), zpd_dt(i), cm2(i), ch2(i),&  !IN
                       &D1(i), beta_old(i), PSIM_hat(i), PSIH_hat(i) )          !OUT
        ELSE
          PSIM_hat(i) = 0
          PSIH_hat(i) = 0
        ENDIF

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong

IF ( BR(I) .GT. 0. ) THEN
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!                                                                                
!-----CLASS 1; STABLE (NIGHTTIME) CONDITIONS:                                    
!
        REGIME(I)=1.
!
! ... paj: psim and psih. Follows Cheng and Brutsaert 2005 (CB05).
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        psim(I)=psim_stable(zol(I))-psim_stable(zol0(i))-PSIM_hat(i)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        psih(I)=psih_stable(zolzz(i))-psih_stable(zol0(i))
        psih(I)=psih_stable(zol(i))-psih_stable(zol0(i))
!
        psim10(I)=psim_stable(zol10)-psim_stable(zol0(i))
        psih10(I)=psih_stable(zol10)-psih_stable(zol0(i))
!
        psim2(I)=psim_stable(zol2)-psim_stable(zol0(i))
        psih2(I)=psih_stable(zol2)-psih_stable(zol0(i))
!
! ... paj: preparations to compute PSIQ. Follows CB05+Carlson Boland JAM 1978.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        pq(I)=psih_stable(zol(I))-psih_stable(zlL)-PSIH_hat(i)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        pq2(I)=psih_stable(zl2)-psih_stable(zol0(i))
        pq10(I)=psih_stable(zl10)-psih_stable(zol0(i))
!
!       1.0 over Monin-Obukhov length
        RMOL(I)=ZOL(I)/ZA(I) 
!                                                                                

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        GOTO 320                                                                 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
ELSE IF ( BR(i) .EQ. 0. ) THEN
!                                                                                
!-----CLASS 3; FORCED CONVECTION:                                                
!                                                                                
!  280   REGIME(I)=3.                                                           
        REGIME(I)=3.                                                           
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!xx        PSIM(I)=0.0                                                              
!xx        PSIH(I)=PSIM(I)                                                          
        PSIM(I)=-PSIM_hat(i)
        PSIH(I)=0.
        PSIM10(I)=0.                                                   
        PSIH10(I)=PSIM10(I)                                           
        PSIM2(I)=0.                                                  
        PSIH2(I)=PSIM2(I)                                           
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!
! paj: preparations to compute PSIQ.
!
!xx        pq(I)=PSIH(I)
        pq(I)=-PSIH_hat(i)
        pq2(I)=PSIH2(I)
        pq10(I)=0.
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!
        ZOL(I)=0.                                             
        RMOL(I) = ZOL(I)/ZA(I)  

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        GOTO 320                                                                 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
ELSE
!                                                                                
!-----CLASS 4; FREE CONVECTION:                                                  
!                                                                                
!  310   CONTINUE                                                                 
        REGIME(I)=4.                                                           
!
! ... paj: PSIM and PSIH ...
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!print*,"befo unstab add psim", I, J
!        psim(I)=psim_unstable(zolzz(i))-psim_unstable(zol0(i))
!        psim(I)=psim_unstable(zolzz(I))-psim_unstable(zol0(i))-PSIM_hat(i)
        psim(I)=psim_unstable(zol(I))-psim_unstable(zol0(i))-PSIM_hat(i)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        psih(I)=psih_unstable(zolzz(i))-psih_unstable(zol0(i))
        psih(I)=psih_unstable(zol(i))-psih_unstable(zol0(i))
!
        psim10(I)=psim_unstable(zol10)-psim_unstable(zol0(i))
        psih10(I)=psih_unstable(zol10)-psih_unstable(zol0(i))
!
        psim2(I)=psim_unstable(zol2)-psim_unstable(zol0(i))
        psih2(I)=psih_unstable(zol2)-psih_unstable(zol0(i))
!
! ... paj: preparations to compute PSIQ 
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!print*,"befo unstab add psih", I, J
!        pq(I)=psih_unstable(zol(I))-psih_unstable(zol0(i))
!        pq(I)=psih_unstable(zol(I))-psih_unstable(zol0(i))-PSIH_hat(i)
        pq(I)=psih_unstable(zol(I))-psih_unstable(zlL)-PSIH_hat(i)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        pq2(I)=psih_unstable(zl2)-psih_unstable(zol0(i))
        pq10(I)=psih_unstable(zl10)-psih_unstable(zol0(i))
!
!---LIMIOT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND HIGH ROUGHNESS            
!---  THIS PREVENTS DENOMINATOR IN FLUXES FROM GETTING TOO SMALL                 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        pq(I)=AMIN1(pq(I),0.9*GZZ0(I))
        PSIM(I)=AMIN1(PSIM(I),0.9*GZZ0(I))
!        pq(I)=AMIN1(pq(I),0.9*ALOG(KARMAN*UST(I)*ZA(I)/XKA+ZA(I)/ZL(I)))
!        PSIM(I)=AMIN1(PSIM(I),0.9*GZ1OZ0(I))

        PSIH(I)=AMIN1(PSIH(I),0.9*GZ1OZ0(I))
!on_gzz        PSIM(I)=AMIN1(PSIM(I),0.9*GZ1OZ0(I))
        PSIH2(I)=AMIN1(PSIH2(I),0.9*GZ2OZ0(I))
        PSIM10(I)=AMIN1(PSIM10(I),0.9*GZ10OZ0(I))
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! AHW: mods to compute ck, cd
        PSIH10(I)=AMIN1(PSIH10(I),0.9*GZ10OZ0(I))

        RMOL(I) = ZOL(I)/ZA(I)  

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!  320 CONTINUE                                                                   
ENDIF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                                                                                
!-----COMPUTE THE FRICTIONAL VELOCITY:                                           
!     ZA(1982) EQS(2.60),(2.61).                                                 
!                                                                                

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!      DO 330 I=its,ite
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        DTG=THX(I)-THGB(I)                                                   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        PSIX=GZ1OZ0(I)-PSIM(I)                                                   
        GZZ0(I) = ALOG((ZA(I))/ZNT(I))
        PSIX=GZZ0(I)-PSIM(I)
!        PSIX=GZ1OZ0(I)-PSIM(I)
!        psix_jun(i) = PSIX

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        PSIX10=GZ10OZ0(I)-PSIM10(I)

!     LOWER LIMIT ADDED TO PREVENT LARGE FLHC IN SOIL MODEL
!     ACTIVATES IN UNSTABLE CONDITIONS WITH THIN LAYERS OR HIGH Z0
!       PSIT=AMAX1(GZ1OZ0(I)-PSIH(I),2.)
       PSIT=GZ1OZ0(I)-PSIH(I)
       PSIT2=GZ2OZ0(I)-PSIH2(I)
!
        IF((XLAND(I)-1.5).GE.0)THEN                                            
          ZL(I)=ZNT(I)                                                            
        ELSE                                                                     
          ZL(I)=0.01                                                                
        ENDIF                                                                    
!
        PSIQ=ALOG(KARMAN*UST(I)*ZA(I)/XKA+ZA(I)/ZL(I))-pq(I)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        psiq_jun(i) = PSIQ

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        PSIQ2=ALOG(KARMAN*UST(I)*2./XKA+2./ZL(I))-pq2(I)

! AHW: mods to compute ck, cd
        PSIQ10=ALOG(KARMAN*UST(I)*10./XKA+10./ZL(I))-pq10(I)


!print*, "PSIQ!?!", i, j
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
IF ( ONOFF == .TRUE. ) THEN ! Profile on? / off? 
IF ( XLAND(i) .LT. 1.5  .and.(VEGTYP(I)>=2.and.VEGTYP(I)<=15)) THEN
! MOST profile - fixed znt
!IF ( ITER .EQ. 1 ) THEN
  zzzol(I) = zol(i)
  zol0(i)  = zol(I)*znt(I)/za(I)          ! z0/L
  zzzolzz(i) = zzzol(i)+zol0(i)

Prc(I) = 0.5 + 0.3 * tanh ( 2. * Lc(I) * zol(I) / za(I) )
ff(I)     = 0.5 * ( sqrt(1. + 4.*0.1*Prc(I) ) -1. )

! RSL profile
  zzz   = za(i) * 1.
  zzzol(I) = zol(i)
  zol0(i)  = zol(I)*znt(I)/za(I)          ! z0/L
  zzzolzz(i) = zzzol(i)+zol0(i)

!##############################################################################_171201

IF (zol(i)==0.) THEN
 MOL(i) = 999999.
ELSE
 MOL(i) = za(i) / zol(i)
ENDIF

!  IF ( MOL(i) < 0 ) THEN
!    IF ( abs(MOL(i)) < 0.1 ) MOL(i) = -0.1
!  ELSE
!    IF ( abs(MOL(i)) < 0.1 ) MOL(i) =  0.1
!  ENDIF
!  cm1(i) = (1-(0.4/(betaN*2)))*exp(cm2(i)/2.) 
  cm1(i) = (1. - ( 0.4        / ( 2.*beta_old(i)*PHIM(zpd_dt(i),MOL(i)) ) ))*exp(cm2(i)/2.)
  ch1(i) = (1. - ( 0.4*PRc(I) / ( 2.*beta_old(i)*PHIH(zpd_dt(i),MOL(i)) ) ))*exp(cm2(i)/2.)
  cm1(i) = (1. - ( 0.4        / ( 2.*betaN )))*exp(cm2(i)/2.)
  ch1(i) = (1. - ( 0.4*PRc(I) / ( 2.*betaN )))*exp(cm2(i)/2.)


  startz = za(i)
  endz   = za(i)*2.
!  endz   = za(i)*10.

  PSIM_hat_temp(i) = 0
  PSIH_hat_temp(i) = 0
  DO zhat = startz, endz-za(i)/2., za(i)/2.
  startz1 = zhat 
  endz1   = zhat + za(i)/2.

  G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
  G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.
  PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * (     (PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
  PSIH_hat_temp(i) = PSIH_hat_temp(i) + (endz1-startz1)/2. * (     (PHIH(G1,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )

  ENDDO


!chs_temp
  PSIH_hat_chs(i) = 0
  startz1 = zpd_dt(i) 
  endz1   = za(i)
  G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
  G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.
  PSIH_hat_chs(i) = (endz1-startz1)/2. * (     (PHIH(G1,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
!chs_temp

!  call PSI_hat (ZL(I), zol0(i), zzzol(I), real(zzz)/10., ZA(I), Lc(I), zpd_dt(i),&    !IN
!               &D1(i), beta(i), PSIM_hat_temp(i), PSIH_hat_temp(i) )         !OUT
!##############################################################################

  IF ( BR(i) > 0 ) THEN
!    psim_temp(i) = psim_stable(zzzolzz(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_stable(zzzol(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_stable(zzzol(i))-psih_stable(zol0(i))-PSIH_hat_temp(i)
  ELSE IF ( BR(i) < 0 ) THEN
!    psim_temp(i) = psim_unstable(zzzolzz(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_unstable(zzzol(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_unstable(zzzol(i))-psih_unstable(zol0(i))-PSIH_hat_temp(i)
  ELSE
!xx    psim_temp(i) = 0
!xx    psih_temp(i) = 0
    psim_temp(i) = -PSIM_hat_temp(i)
    psih_temp(i) = -PSIH_hat_temp(i)
  ENDIF

  factorm1 = WSPD(I) / ( alog( ZA(i) / znt(i) ) - psim_temp(i)  )
  factorh1 = ( SCR3(I) - TSK(i) ) / ( alog( 0.4*0.4*factorm1*za(I)/xka+za(i)/zl(i)  ) - psih_temp(i)  )
  chs_i1(i) =  alog( za(i)/znt(i) ) - psih_temp(i)
!chs_temp
  chs_i1(i) =  - psih_temp(i) - PSIM_hat_temp(i)
!chs_temp


  DO zzz = int(za(I)*1.), int(zpd_dt(i)*1.)+1, -1
  IF ( zzz == 10 .or. zzz == 2 .or. zzz == int(zpd_dt(i)*1)+1 ) THEN

!IF ( I == 16 .and. J == 48 ) print*, " TEST end1",  itimestep, I,J,zzz
    zzzol(I) = zol(i) / ZA(i) * zzz/1.
    zol0(i)=zol(I)*znt(I)/za(I)          ! z0/L
    zzzolzz(i) = zzzol(i)+zol0(i)

!##############################################################################_171201
    startz = real(zzz)/1. 
    endz   = za(i)*2.
!    endz   = za(i)*10.

PSIM_hat_temp(i) = 0
PSIH_hat_temp(i) = 0

IF ( startz <= za(i) ) then
startz1 = startz
endz1 = za(i)

G1  = (endz1-startz1)/2*(-1/sqrt(3.)) + (endz1+startz1)/2
G2  = (endz1-startz1)/2*( 1/sqrt(3.)) + (endz1+startz1)/2
PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2 * (     (PHIM(G1,MOL(i)) * (1-(1.-cm1(i)*exp( cm2(i)/(2*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1-(1.-cm1(i)*exp( cm2(i)/(2*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
PSIH_hat_temp(i) = PSIH_hat_temp(i) + (endz1-startz1)/2 * (     (PHIH(G1,MOL(i)) * (1-(1.-ch1(i)*exp( ch2(i)/(2*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1-(1.-ch1(i)*exp( ch2(i)/(2*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
ENDIF

DO zhat = za(i), endz-za(i)/2., za(i)/2.
startz1 = zhat
endz1 = zhat+za(i)/2.

    G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
    G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.
    PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * (     (PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
    PSIH_hat_temp(i) = PSIH_hat_temp(i) + (endz1-startz1)/2. * (     (PHIH(G1,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )

ENDDO
!    call PSI_hat (ZL(I), zol0(i), zzzol(I), real(zzz)/10., ZA(I), Lc(I), zpd_dt(i),&    !IN
!                 &D1(i), beta(i), PSIM_hat_temp(i), PSIH_hat_temp(i) )         !OUT
!##############################################################################

! print*, "windp1", itimestep, i, j
  IF ( BR(i) > 0 ) THEN
!    psim_temp(i) = psim_stable(zzzolzz(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_stable(zzzol(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_stable(zzzol(i))-psih_stable(zol0(i))-PSIH_hat_temp(i)
  ELSE IF ( BR(i) < 0 ) THEN
!    psim_temp(i) = psim_unstable(zzzolzz(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_unstable(zzzol(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_unstable(zzzol(i))-psih_unstable(zol0(i))-PSIH_hat_temp(i)
  ELSE
!xx    psim_temp(i) = 0
!xx    psih_temp(i) = 0
    psim_temp(i) = -PSIM_hat_temp(i)
    psih_temp(i) = -PSIH_hat_temp(i)
  ENDIF

    Uwind(i,zzz) = factorm1 * ( alog( zzz/1. / znt(i) ) - psim_temp(i) )
    Twind(i,zzz) = factorh1 * ( alog( 0.4*0.4*factorm1*(zzz/1.)/xka+(zzz/1.)/zl(i) ) - psih_temp(i) ) + TSK(i)

  ENDIF
  ENDDO

  zzzol(I) = zol(i) / ZA(i) * zpd_dt(i) 
  zol0(i)  = zol(I)*znt(I)/za(I)          ! z0/L
  zzzolzz(i) = zzzol(i)+zol0(i)

!##############################################################################_171201
  startz = zpd_dt(i) 
  endz   = za(i)*2.
!  endz   = za(i)*10.

PSIM_hat_temp(i) = 0
PSIH_hat_temp(i) = 0

startz1 = startz
endz1 = za(i)
G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.
PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * ((PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
PSIH_hat_temp(i) = PSIH_hat_temp(i) + (endz1-startz1)/2. * ((PHIH(G1,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )

DO zhat = za(i), endz-za(i)/2., za(i)/2.
startz1 = zhat
endz1 = zhat+za(i)/2.

  G1  = (endz1-startz1)/2.*(-1./sqrt(3.)) + (endz1+startz1)/2.
  G2  = (endz1-startz1)/2.*( 1./sqrt(3.)) + (endz1+startz1)/2.
  PSIM_hat_temp(i) = PSIM_hat_temp(i) + (endz1-startz1)/2. * ((PHIM(G1,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIM(G2,MOL(i)) * (1.-(1.-cm1(i)*exp( cm2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )
  PSIH_hat_temp(i) = PSIH_hat_temp(i) + (endz1-startz1)/2. * ((PHIH(G1,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G1*(-1.) ))) / G1)&
& + (PHIH(G2,MOL(i)) * (1.-(1.-ch1(i)*exp( ch2(i)/(2.*beta_old(i)*beta_old(i)*Lc(i))*G2*(-1.) ))) / G2) )

ENDDO

!  call PSI_hat (ZL(I), zol0(i), zzzol(I), zpd_dt(i), ZA(I), Lc(I), zpd_dt(i),&     !IN
!               &D1(i), beta(i), PSIM_hat_temp(i), PSIH_hat_temp(i) )      !OUT
!##############################################################################

  IF ( BR(i) > 0 ) THEN
!    psim_temp(i) = psim_stable(zzzolzz(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_stable(zzzol(I))-psim_stable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_stable(zzzol(i))-psih_stable(zol0(i))-PSIH_hat_temp(i)
  ELSE IF ( BR(i) < 0 ) THEN
!    psim_temp(i) = psim_unstable(zzzolzz(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psim_temp(i) = psim_unstable(zzzol(I))-psim_unstable(zol0(i))-PSIM_hat_temp(i)
    psih_temp(i) = psih_unstable(zzzol(i))-psih_unstable(zol0(i))-PSIH_hat_temp(i)
  ELSE
!xx    psim_temp(i) = 0
!xx    psih_temp(i) = 0
    psim_temp(i) = -PSIM_hat_temp(i)
    psih_temp(i) = -PSIH_hat_temp(i)
  ENDIF

 Uhc(I) = factorm1 * ( alog( zpd_dt(i) / znt(i) ) - psim_temp(i) ) 
 Thc(I) = factorh1 * ( alog( 0.4*0.4*factorm1*zpd_dt(i)/xka+zpd_dt(i)/zl(i) ) - psih_temp(i) ) + TSK(i)
!chs_temp
 chs_i2(i) =  - psih_temp(i) -PSIM_hat_temp(i)
!chs_temp
!dwind
 dwind(i) = (Uwind(i,int(zpd_dt(i)*1.)+1)-Uhc(I))/1
!dwind


 UST(I)=KARMAN*WSPD(I)/PSIX                                             
!dwind
  zzz = int(2)
  Uwind(i,zzz) = Uhc(i) *       exp( (real(zzz)/1.-zpd_dt(i)) * dwind(i) / Uhc(i) )
!dwind
  Twind(i,zzz) = (Thc(i)-TSK(I) ) * exp( dtemp(i)/(Thc(i)-TSK(I))*(real(zzz)/1.-zpd_dt(i)) ) + TSK(I)
  zzz = int(10)
  Uwind(i,zzz) = Uhc(i) *       exp( (real(zzz)/1.-zpd_dt(i)) * dwind(i) / Uhc(i) )
  Twind(i,zzz) = (Thc(i)-TSK(I) ) * exp( dtemp(i)/(Thc(i)-TSK(I))*(real(zzz)/1.-zpd_dt(i)) ) + TSK(I)


ENDIF
ENDIF ! profile on? / off? ve
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

m10_wind(i) = uwind(i,10)/WSPD(I) * sqrt(UX(i)*UX(i)+VX(i)*VX(i))
m2_temp(i) = Twind(i,2)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


! V3.7: using Fairall 2003 to compute z0q and z0t over water:
!       adapted from module_sf_mynn.F
        IF ( (XLAND(I)-1.5).GE.0. ) THEN
              VISC=(1.32+0.009*(SCR3(I)-273.15))*1.E-5
              RESTAR=UST(I)*ZNT(I)/VISC
              Z0T = (5.5e-5)*(RESTAR**(-0.60))
              Z0T = MIN(Z0T,1.0e-4)
              Z0T = MAX(Z0T,2.0e-9)
              Z0Q = Z0T

              PSIQ=max(ALOG((ZA(I)+Z0Q)/Z0Q)-PSIH(I), 2.)
              PSIT=max(ALOG((ZA(I)+Z0T)/Z0T)-PSIH(I), 2.)
              PSIQ2=max(ALOG((2.+Z0Q)/Z0Q)-PSIH2(I), 2.)
              PSIT2=max(ALOG((2.+Z0T)/Z0T)-PSIH2(I), 2.)
              PSIQ10=max(ALOG((10.+Z0Q)/Z0Q)-PSIH10(I), 2.)
        ENDIF

        IF ( PRESENT(ISFTCFLX) ) THEN
           IF ( ISFTCFLX.EQ.1 .AND. (XLAND(I)-1.5).GE.0. ) THEN
! v3.1
!             Z0Q = 1.e-4 + 1.e-3*(MAX(0.,UST(I)-1.))**2
! hfip1
!             Z0Q = 0.62*2.0E-5/UST(I) + 1.E-3*(MAX(0.,UST(I)-1.5))**2
! v3.2
              Z0Q = 1.e-4
!
! ... paj: recompute psih for z0q
!
           zolzz(i)=zol(I)*(za(I)+z0q)/za(I)    ! (z+z0q)/L
           zol10=zol(I)*(10.+z0q)/za(I)   ! (10+z0q)/L
           zol2=zol(I)*(2.+z0q)/za(I)     ! (2+z0q)/L
           zol0(i)=zol(I)*z0q/za(I)          ! z0q/L
!
              if (zol(I).gt.0.) then
              psih(I)=psih_stable(zolzz(i))-psih_stable(zol0(i))
              psih10(I)=psih_stable(zol10)-psih_stable(zol0(i))
              psih2(I)=psih_stable(zol2)-psih_stable(zol0(i))
              else
                if (zol(I).eq.0) then
                psih(I)=0.
                psih10(I)=0.
                psih2(I)=0.
                else
                psih(I)=psih_unstable(zolzz(i))-psih_unstable(zol0(i))
                psih10(I)=psih_unstable(zol10)-psih_unstable(zol0(i))
                psih2(I)=psih_unstable(zol2)-psih_unstable(zol0(i))
                endif
              endif
!
              PSIQ=ALOG((ZA(I)+z0q)/Z0Q)-PSIH(I)
              PSIT=PSIQ
              PSIQ2=ALOG((2.+z0q)/Z0Q)-PSIH2(I)
              PSIQ10=ALOG((10.+z0q)/Z0Q)-PSIH10(I)
              PSIT2=PSIQ2
           ENDIF
          IF ( ISFTCFLX.EQ.2 .AND. (XLAND(I)-1.5).GE.0. ) THEN
! AHW: Garratt formula: Calculate roughness Reynolds number
!        Kinematic viscosity of air (linear approc to
!                 temp dependence at sea level)
! GZ0OZT and GZ0OZQ are based off formulas from Brutsaert (1975), which
! Garratt (1992) used with values of k = 0.40, Pr = 0.71, and Sc = 0.60
              VISC=(1.32+0.009*(SCR3(I)-273.15))*1.E-5
!!            VISC=1.5E-5
              RESTAR=UST(I)*ZNT(I)/VISC
              GZ0OZT=0.40*(7.3*SQRT(SQRT(RESTAR))*SQRT(0.71)-5.)
!
! ... paj: compute psih for z0t for temperature ...
!
              z0t=znt(I)/exp(GZ0OZT)
!
           zolzz(i)=zol(I)*(za(I)+z0t)/za(I)    ! (z+z0t)/L
           zol10=zol(I)*(10.+z0t)/za(I)   ! (10+z0t)/L
           zol2=zol(I)*(2.+z0t)/za(I)     ! (2+z0t)/L
           zol0(i)=zol(I)*z0t/za(I)          ! z0t/L
!
              if (zol(I).gt.0.) then
              psih(I)=psih_stable(zolzz(i))-psih_stable(zol0(i))
              psih10(I)=psih_stable(zol10)-psih_stable(zol0(i))
              psih2(I)=psih_stable(zol2)-psih_stable(zol0(i))
              else
                if (zol(I).eq.0) then
                psih(I)=0.
                psih10(I)=0.
                psih2(I)=0.
                else
                psih(I)=psih_unstable(zolzz(i))-psih_unstable(zol0(i))
                psih10(I)=psih_unstable(zol10)-psih_unstable(zol0(i))
                psih2(I)=psih_unstable(zol2)-psih_unstable(zol0(i))
                endif
              endif
!
!              PSIT=GZ1OZ0(I)-PSIH(I)+RESTAR2
!              PSIT2=GZ2OZ0(I)-PSIH2(I)+RESTAR2
              PSIT=ALOG((ZA(I)+z0t)/Z0t)-PSIH(I)
              PSIT2=ALOG((2.+z0t)/Z0t)-PSIH2(I)
!
              GZ0OZQ=0.40*(7.3*SQRT(SQRT(RESTAR))*SQRT(0.60)-5.)
              z0q=znt(I)/exp(GZ0OZQ)
!
           zolzz(i)=zol(I)*(za(I)+z0q)/za(I)    ! (z+z0q)/L
           zol10=zol(I)*(10.+z0q)/za(I)   ! (10+z0q)/L
           zol2=zol(I)*(2.+z0q)/za(I)     ! (2+z0q)/L
           zol0(i)=zol(I)*z0q/za(I)          ! z0q/L
!
              if (zol(I).gt.0.) then
              psih(I)=psih_stable(zolzz(i))-psih_stable(zol0(i))
              psih10(I)=psih_stable(zol10)-psih_stable(zol0(i))
              psih2(I)=psih_stable(zol2)-psih_stable(zol0(i))
              else
                if (zol(I).eq.0) then
                psih(I)=0.
                psih10(I)=0.
                psih2(I)=0.
                else
                psih(I)=psih_unstable(zolzz(i))-psih_unstable(zol0(i))
                psih10(I)=psih_unstable(zol10)-psih_unstable(zol0(i))
                psih2(I)=psih_unstable(zol2)-psih_unstable(zol0(i))
                endif
              endif
!
              PSIQ=ALOG((ZA(I)+z0q)/Z0q)-PSIH(I)
              PSIQ2=ALOG((2.+z0q)/Z0q)-PSIH2(I)
              PSIQ10=ALOG((10.+z0q)/Z0q)-PSIH10(I)
!              PSIQ=GZ1OZ0(I)-PSIH(I)+2.28*SQRT(SQRT(RESTAR))-2.
!              PSIQ2=GZ2OZ0(I)-PSIH2(I)+2.28*SQRT(SQRT(RESTAR))-2.
!              PSIQ10=GZ10OZ0(I)-PSIH(I)+2.28*SQRT(SQRT(RESTAR))-2.
           ENDIF
        ENDIF
        IF(PRESENT(ck) .and. PRESENT(cd) .and. PRESENT(cka) .and. PRESENT(cda)) THEN
           Ck(I)=(karman/psix10)*(karman/psiq10)
           Cd(I)=(karman/psix10)*(karman/psix10)
           Cka(I)=(karman/psix)*(karman/psiq)
           Cda(I)=(karman/psix)*(karman/psix)
        ENDIF
        IF ( PRESENT(IZ0TLND) ) THEN
           IF ( IZ0TLND.EQ.1 .AND. (XLAND(I)-1.5).LE.0. ) THEN
              ZL(I)=ZNT(I)
!             CZIL RELATED CHANGES FOR LAND
              VISC=(1.32+0.009*(SCR3(I)-273.15))*1.E-5
              RESTAR=UST(I)*ZL(I)/VISC
!             Modify CZIL according to Chen & Zhang, 2009

              CZIL = 10.0 ** ( -0.40 * ( ZL(I) / 0.07 ) )
!
! ... paj: compute phish for z0t over land
!
              z0t=znt(I)/exp(CZIL*KARMAN*SQRT(RESTAR))
!
           zolzz(i)=zol(I)*(za(I)+z0t)/za(I)    ! (z+z0t)/L
           zol10=zol(I)*(10.+z0t)/za(I)   ! (10+z0t)/L
           zol2=zol(I)*(2.+z0t)/za(I)     ! (2+z0t)/L
           zol0(i)=zol(I)*z0t/za(I)          ! z0t/L
!
              if (zol(I).gt.0.) then
              psih(I)=psih_stable(zolzz(i))-psih_stable(zol0(i))
              psih10(I)=psih_stable(zol10)-psih_stable(zol0(i))
              psih2(I)=psih_stable(zol2)-psih_stable(zol0(i))
              else
                if (zol(I).eq.0) then
                psih(I)=0.
                psih10(I)=0.
                psih2(I)=0.
                else
                psih(I)=psih_unstable(zolzz(i))-psih_unstable(zol0(i))
                psih10(I)=psih_unstable(zol10)-psih_unstable(zol0(i))
                psih2(I)=psih_unstable(zol2)-psih_unstable(zol0(i))
                endif
              endif
!
              PSIQ=ALOG((ZA(I)+z0t)/Z0t)-PSIH(I)
              PSIQ2=ALOG((2.+z0t)/Z0t)-PSIH2(I)
              PSIT=PSIQ
              PSIT2=PSIQ2
!
!              PSIT=GZ1OZ0(I)-PSIH(I)+CZIL*KARMAN*SQRT(RESTAR)
!              PSIQ=GZ1OZ0(I)-PSIH(I)+CZIL*KARMAN*SQRT(RESTAR)
!              PSIT2=GZ2OZ0(I)-PSIH2(I)+CZIL*KARMAN*SQRT(RESTAR)
!              PSIQ2=GZ2OZ0(I)-PSIH2(I)+CZIL*KARMAN*SQRT(RESTAR)

           ENDIF
        ENDIF
! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!        UST(I)=0.5*UST(I)+0.5*KARMAN*WSPD(I)/PSIX                                             
        UST(I)=KARMAN*WSPD(I)/PSIX                                             
        UST_JUN(I)=UST(I)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TKE coupling: compute ust without vconv for use in tke scheme
        WSPDI(I)=SQRT(UX(I)*UX(I)+VX(I)*VX(I))
        IF ( PRESENT(USTM) ) THEN
        USTM(I)=0.5*USTM(I)+0.5*KARMAN*WSPDI(I)/PSIX
        ENDIF

        U10(I)=UX(I)*PSIX10/PSIX                                    
        V10(I)=VX(I)*PSIX10/PSIX                                   
        TH2(I)=THGB(I)+DTG*PSIT2/PSIT                                
        Q2(I)=QSFC(I)+(QX(I)-QSFC(I))*PSIQ2/PSIQ                   
        T2(I) = TH2(I)*(PSFCPA(I)/P1000mb)**ROVCP                     
!                                                                                
        IF((XLAND(I)-1.5).LT.0.)THEN                                            
          UST(I)=AMAX1(UST(I),0.001)
        ENDIF                                                                    
        MOL(I)=KARMAN*DTG/PSIT/PRT                              
        DENOMQ(I)=PSIQ
        DENOMQ2(I)=PSIQ2
        DENOMT2(I)=PSIT2
        FM(I)=PSIX
 !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
        IF ((XLAND(I)-1.5).LT.0.) THEN
          FH(I)=PSIQ
        ELSE
          FH(I)=PSIT
        ENDIF

!        FH(I)=PSIT

IF( (XLAND(I)-1.5).LT.0  .and. (VEGTYP(I)>=2.and.VEGTYP(I)<=15) )THEN
  !m10_wind(i) = uwind(i,100)/WSPD(I) * sqrt(UX(i)*UX(i)+VX(i)*VX(i))
ELSE
  m10_wind(i) = sqrt(U10(I)*U10(I)+V10(I)*V10(I))
  m2_temp(i) = T2(I)
ENDIF

 !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!  330 CONTINUE                                                                   

888 continue
ENDDO  ! LOOP2


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                                                                                
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
!  335 CONTINUE                                                                   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                                                                                  
!-----COMPUTE THE SURFACE SENSIBLE AND LATENT HEAT FLUXES:                       
      IF ( PRESENT(SCM_FORCE_FLUX) ) THEN
         IF (SCM_FORCE_FLUX.EQ.1) GOTO 350
      ENDIF
      DO i=its,ite
        QFX(i)=0.                                                              
        HFX(i)=0.                                                              
      ENDDO
  350 CONTINUE                                                                   

      IF (ISFFLX.EQ.0) GOTO 410                                                
                                                                                 
!-----OVER WATER, ALTER ROUGHNESS LENGTH (ZNT) ACCORDING TO WIND (UST).          
                                                                                 
      DO 360 I=its,ite
        IF((XLAND(I)-1.5).GE.0)THEN                                            
!         ZNT(I)=CZO*UST(I)*UST(I)/G+OZO                                   
! Since V3.7 (ref: EC Physics document for Cy36r1)
          ZNT(I)=CZO*UST(I)*UST(I)/G+0.11*1.5E-5/UST(I)
! COARE 3.5 (Edson et al. 2013)
!         CZC = 0.0017*WSPD(I)-0.005
!         CZC = min(CZC,0.028)
!         ZNT(I)=CZC*UST(I)*UST(I)/G+0.11*1.5E-5/UST(I)
! AHW: change roughness length, and hence the drag coefficients Ck and Cd
          IF ( PRESENT(ISFTCFLX) ) THEN
             IF ( ISFTCFLX.NE.0 ) THEN
!               ZNT(I)=10.*exp(-9.*UST(I)**(-.3333))
!               ZNT(I)=10.*exp(-9.5*UST(I)**(-.3333))
!               ZNT(I)=ZNT(I) + 0.11*1.5E-5/AMAX1(UST(I),0.01)
!               ZNT(I)=0.011*UST(I)*UST(I)/G+OZO
!               ZNT(I)=MAX(ZNT(I),3.50e-5)
! AHW 2012:
                ZW  = MIN((UST(I)/1.06)**(0.3),1.0)
                ZN1 = 0.011*UST(I)*UST(I)/G + OZO
                ZN2 = 10.*exp(-9.5*UST(I)**(-.3333)) + &
                       0.11*1.5E-5/AMAX1(UST(I),0.01)
                ZNT(I)=(1.0-ZW) * ZN1 + ZW * ZN2
                ZNT(I)=MIN(ZNT(I),2.85e-3)
                ZNT(I)=MAX(ZNT(I),1.27e-7)
             ENDIF
          ENDIF
          ZL(I) = ZNT(I)
        ELSE
          ZL(I) = 0.01
        ENDIF                                                                    
        FLQC(I)=RHOX(I)*MAVAIL(I)*UST(I)*KARMAN/DENOMQ(I)
!       FLQC(I)=RHOX(I)*MAVAIL(I)*UST(I)*KARMAN/(   &
!               ALOG(KARMAN*UST(I)*ZA(I)/XKA+ZA(I)/ZL(I))-PSIH(I))
        DTTHX=ABS(THX(I)-THGB(I))                                            
        IF(DTTHX.GT.1.E-5)THEN                                                   
          FLHC(I)=CPM(I)*RHOX(I)*UST(I)*MOL(I)/(THX(I)-THGB(I))          
!         write(*,1001)FLHC(I),CPM(I),RHOX(I),UST(I),MOL(I),THX(I),THGB(I),I
 1001   format(f8.5,2x,f12.7,2x,f12.10,2x,f12.10,2x,f13.10,2x,f12.8,f12.8,2x,i3)
        ELSE                                                                     
          FLHC(I)=0.                                                             
        ENDIF                                                                    
  360 CONTINUE                                                                   

!                                                                                
!-----COMPUTE SURFACE MOIST FLUX:                                                
!                                                                                
!     IF(IDRY.EQ.1)GOTO 390                                                
!                                                                                
     IF ( PRESENT(SCM_FORCE_FLUX) ) THEN
        IF (SCM_FORCE_FLUX.EQ.1) GOTO 405
     ENDIF

      DO 370 I=its,ite
        QFX(I)=FLQC(I)*(QSFC(I)-QX(I))                                     
        QFX(I)=AMAX1(QFX(I),0.)                                            
        LH(I)=XLV*QFX(I)
  370 CONTINUE                                                                 
                                                                                
!-----COMPUTE SURFACE HEAT FLUX:                                                 
!                                                                                
  390 CONTINUE                                                                 
      DO 400 I=its,ite
        IF(XLAND(I)-1.5.GT.0.)THEN                                           
          HFX(I)=FLHC(I)*(THGB(I)-THX(I)) 
!         IF ( PRESENT(ISFTCFLX) ) THEN
!            IF ( ISFTCFLX.NE.0 ) THEN
! AHW: add dissipative heating term (commented out in 3.6.1)
!               HFX(I)=HFX(I)+RHOX(I)*USTM(I)*USTM(I)*WSPDI(I)
!            ENDIF
!         ENDIF 
        ELSEIF(XLAND(I)-1.5.LT.0.)THEN                                       
          HFX(I)=FLHC(I)*(THGB(I)-THX(I))                                
          HFX(I)=AMAX1(HFX(I),-250.)                                       
        ENDIF                                                                  
  400 CONTINUE                                                                 

  405 CONTINUE                                                                 
         
      DO I=its,ite
         IF((XLAND(I)-1.5).GE.0)THEN
           ZL(I)=ZNT(I)
         ELSE
           ZL(I)=0.01
         ENDIF
!v3.1.1
!         CHS(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*ZA(I) &
!                /XKA+ZA(I)/ZL(I))-PSIH(I))
         CHS(I)=UST(I)*KARMAN/DENOMQ(I)
         CHS_ori(i) = CHS(i)

IF( (XLAND(I)-1.5).LT.0  .and. (VEGTYP(I)>=2.and.VEGTYP(I)<=15) )THEN

Prc(I) = 0.5 + 0.3 * tanh ( 2. * Lc(I) * zol(I) / za(I) )
!chs_temp
  chs_temp(i) = ( alog(za(i)/zpd_dt(i)) + chs_i1(i) - chs_i2(i) - PSIH_hat_chs(i))
  IF ( alog(za(i)/zpd_dt(i)) *0.9 < - chs_i1(i) + chs_i2(i) + psih_hat_chs(i) ) chs_temp(i) = alog(za(i)/zpd_dt(i)) * 0.1 
  CHS(I) = 1./ ( 1./ (KARMAN * UST(i) / ( chs_temp(i) )) + &
!chs_temp
!zl_add
 &               1./ (KARMAN * UST(i) / ( ALOG(KARMAN*UST(I)*ZL(I)/XKA+1) )) + &
 &  Prc(i) / (beta_old(i)*beta_old(i)*Uhc(i)) * ( exp( (zpd_dt(i)-zl(i))/(2.*beta_old(i)*beta_old(i)*Lc(i)) ) -1. ) )
!zl_add


IF (zpd_dt(i) <= znt(i) )  print*, "RSH WARN !!", itimestep, i, J,zpd_dt(i), znt(i)

 IF ( CHS(i) <= 0.00001 ) then

print*,"PRC2", i,j, PRC(i), 2. * Lc(I) * zol(I) / za(I), tanh(2. * Lc(I) * zol(I) / za(I))
 CHS(i) = CHS_ori(i)
 endif
 IF ( CHS(i) >= 1     ) then

print*,"PRC1", i,j, PRC(i), 2. * Lc(I) * zol(I) / za(I), tanh(2. * Lc(I) * zol(I) / za(I))
ENDIF


ENDIF
!        GZ2OZ0(I)=ALOG(2./ZNT(I))
!        PSIM2(I)=-10.*GZ2OZ0(I)
!        PSIM2(I)=AMAX1(PSIM2(I),-10.)
!        PSIH2(I)=PSIM2(I)
! v3.1.1
!         CQS2(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*2.0  &
!               /XKA+2.0/ZL(I))-PSIH2(I))
!         CHS2(I)=UST(I)*KARMAN/(GZ2OZ0(I)-PSIH2(I))
         CQS2(I)=UST(I)*KARMAN/DENOMQ2(I)
         CHS2(I)=UST(I)*KARMAN/DENOMT2(I)
      ENDDO
                                                                        
  410 CONTINUE                                                                   
!jdf
!     DO I=its,ite
!       IF(UST(I).GE.0.1) THEN
!         RMOL(I)=RMOL(I)*(-FLHC(I))/(UST(I)*UST(I)*UST(I))
!       ELSE
!         RMOL(I)=RMOL(I)*(-FLHC(I))/(0.1*0.1*0.1)
!       ENDIF
!     ENDDO
!jdf

!                                                                                


   END SUBROUTINE SFCLAYREV1D

!====================================================================
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
   SUBROUTINE sfclayrevinit_rsl

    INTEGER                   ::      N
    REAL                      ::      zolf

    DO N=0,1000
! stable function tables
       zolf = float(n)*0.01
       psim_stab(n)=psim_stable_full(zolf)
       psih_stab(n)=psih_stable_full(zolf)
 
! unstable function tables
       zolf = -float(n)*0.01
       psim_unstab(n)=psim_unstable_full(zolf)
       psih_unstab(n)=psih_unstable_full(zolf)

    ENDDO

   END SUBROUTINE sfclayrevinit_rsl
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
subroutine RPSI_hat(ZL, ZOL0, ZOL, ZZ, ZLVL, Lc, zpd_dt, cm2, ch2, & !IN
                  &D1, beta,PSIM_hat, PSIH_hat )     !OUT
IMPLICIT NONE

!IN
REAL,intent(IN) :: ZL,ZOL0,ZOL, ZZ, ZLVL, Lc, zpd_dt  ! zl=0.01/L, ZOL0=z0/L, ZOL=z/L, ZZ=z, ZLVL=zr
REAL,intent(IN) :: beta
REAL,intent(IN) :: cm2, ch2

!OUT
REAL,intent(INOUT) :: PSIM_hat, PSIH_hat
REAL,intent(INOUT) :: D1

!Local
REAL :: startz, endz,ii, startz1,  endz1, zhat
REAL :: MOL
REAL :: cm1, ch1
REAL :: Prc, ff, R, VKC, phim_hat, phih_hat
REAL :: betaNO, betaHF
REAL :: G1, G2

REAL, PARAMETER :: betaN = 0.374

PSIM_hat = 0
PSIH_hat = 0

IF ( ZL == 0 .or. ZOL0 == 0 .or. ZOL == 0 .or. ZLVL == 0 ) then
ENDIF

IF (zol==0.) THEN
 MOL = 999999.
ELSE
 MOL    = ZZ/zol      !  z / ( z/L )
ENDIF

R      = 0.1
VKC    = 0.4

Prc    = 0.5 + 0.3 * tanh ( 2 * Lc / MOL )
ff     = 0.5 * ( sqrt(1 + 4*R*Prc ) -1 )
Cm1 = (1 - ( VKC       / ( 2*beta*PHIM(zpd_dt,MOL) ) ))*exp(Cm2/2) ! v2.2-1 change
Ch1 = (1 - ( VKC * Prc / ( 2*beta*PHIH(zpd_dt,MOL) ) ))*exp(Ch2/2)
Cm1 = (1 - ( VKC       / ( betaN*2 ) ))*exp(Cm2/2) ! v2.2-1 change
Ch1 = (1 - ( VKC * Prc / ( betaN*2 ) ))*exp(Ch2/2)
D1  = 2*beta*beta*Lc    ! 1/Zrsl ??

startz = zz
endz   = zlvl*2

PSIM_hat = 0; PSIH_hat = 0
IF ( zz <= zlvl ) THEN
startz1 = startz 
endz1 = zlvl
G1 = (endz1-startz1)/2 * (-1/sqrt(3.)) + (startz1+endz1)/2
G2 = (endz1-startz1)/2 * ( 1/sqrt(3.)) + (startz1+endz1)/2
PSIM_hat = PSIM_hat + (endz1-startz1)/2 * ( (PHIM(G1,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIM(G2,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G2  * (-1.) ))) / G2) )
PSIH_hat = PSIH_hat + (endz1-startz1)/2 * ( (PHIH(G1,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIH(G2,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G2  * (-1.) ))) / G2) )
ENDIF

DO zhat = zlvl, endz-zlvl/2, zlvl/2
startz1 = zhat
endz1 = zhat+zlvl/2

G1 = (endz1-startz1)/2 * (-1/sqrt(3.)) + (startz1+endz1)/2
G2 = (endz1-startz1)/2 * ( 1/sqrt(3.)) + (startz1+endz1)/2
PSIM_hat = PSIM_hat + (endz1-startz1)/2 * ( (PHIM(G1,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIM(G2,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G2  * (-1.) ))) / G2) )
PSIH_hat = PSIH_hat + (endz1-startz1)/2 * ( (PHIH(G1,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIH(G2,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G2  * (-1.) ))) / G2) )

ENDDO
999 continue

end subroutine


subroutine PSI_hat(ZL, ZOL0, ZOL, ZZ, ZLVL, Lc, hc,& !IN
                  &D1, beta,PSIM_hat, PSIH_hat )     !OUT
IMPLICIT NONE

!IN
!REAL,intent(IN) :: ZL,ZOL0,ZOL, ZZ, ZLVL, Lc, zpd_dt  ! zl=0.01/L, ZOL0=z0/L, ZOL=z/L, ZZ=z, ZLVL=zr
REAL,intent(IN) :: ZL,ZOL0,ZOL, ZZ, ZLVL, Lc, Hc

!OUT
REAL,intent(INOUT) :: PSIM_hat, PSIH_hat
REAL,intent(INOUT) :: D1, beta

!Local
REAL :: startz, endz,ii, zhat, startz1, endz1 
REAL :: MOL
REAL :: cm1, cm2, ch1,ch2
REAL :: Prc, ff, R, VKC, phim_hat, phih_hat
REAL :: betaNO, betaHF, beta1, beta2
REAL :: G1, G2
REAL :: zpd_dt
INTEGER :: ITER

REAL, PARAMETER :: betaN = 0.374

PSIM_hat = 0
PSIH_hat = 0

IF ( ZL == 0 .or. ZOL0 == 0 .or. ZOL == 0 .or. ZLVL == 0 ) then
ENDIF

IF (zol==0.) THEN
 MOL = 999999.
ELSE
 MOL    = ZZ/zol      !  z / ( z/L )
ENDIF

R      = 0.1
VKC    = 0.4


ITER = 1
zpd_dt = Hc*0.35

beta = betaN   /PHIM(zpd_dt,MOL)

betaHF = betaN   /PHIM(zpd_dt,MOL)
betaNO = 0.4/2./PHIM(zpd_dt,MOL)
beta   = betaNO + (BetaHF-BetaNO)/(1+2* ((abs(Lc/MOL+0.15))**1.5)  )
IF ( Lc/MOL > -0.15 ) beta = betaHF

IF ( beta <= 0.01 ) beta = 0.01

Prc    = 0.5 + 0.3 * tanh ( 2 * Lc / MOL )
ff     = 0.5 * ( sqrt(1 + 4*R*Prc ) -1 )
Cm2 = 0.5; Ch2 = 0.5

Cm1 = (1 - ( VKC       / ( 2*beta*PHIM(zpd_dt,MOL) ) ))*exp(Cm2/2) ! v2.2-1 change
Ch1 = (1 - ( VKC * Prc / ( 2*beta*PHIH(zpd_dt,MOL) ) ))*exp(Ch2/2)
Cm1 = (1 - ( VKC       / ( betaN*2 ) ))*exp(Cm2/2) ! v2.2-1 change
Ch1 = (1 - ( VKC * Prc / ( betaN*2 ) ))*exp(Ch2/2)
D1  = 2*beta*beta*Lc    ! 1/Zrsl ?? 

startz = zz
endz   = zlvl*2

PSIM_hat = 0; PSIH_hat = 0
IF ( zz <= zlvl ) THEN
startz1 = startz
endz1 = zlvl
G1 = (endz1-startz1)/2 * (-1/sqrt(3.)) + (startz1+endz1)/2
G2 = (endz1-startz1)/2 * ( 1/sqrt(3.)) + (startz1+endz1)/2
PSIM_hat = PSIM_hat + (endz1-startz1)/2 * ( (PHIM(G1,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIM(G2,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G2  * (-1.) ))) / G2) )
PSIH_hat = PSIH_hat + (endz1-startz1)/2 * ( (PHIH(G1,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIH(G2,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G2  * (-1.) ))) / G2) )
ENDIF

DO zhat = zlvl, endz-zlvl/2, zlvl/2
startz1 = zhat
endz1 = zhat+zlvl/2

G1 = (endz1-startz1)/2 * (-1/sqrt(3.)) + (startz1+endz1)/2
G2 = (endz1-startz1)/2 * ( 1/sqrt(3.)) + (startz1+endz1)/2
PSIM_hat = PSIM_hat + (endz1-startz1)/2 * ( (PHIM(G1,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIM(G2,MOL) * (1-(1. - cm1 * exp( cm2/d1 * G2  * (-1.) ))) / G2) )
PSIH_hat = PSIH_hat + (endz1-startz1)/2 * ( (PHIH(G1,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G1  * (-1.) ))) / G1)&
& + (PHIH(G2,MOL) * (1-(1. - ch1 * exp( ch2/d1 * G2  * (-1.) ))) / G2) )

ENDDO

999 continue

end subroutine
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_junhong
function PHIM(zzz,MOLe)
real :: PHIM
real :: zzz, MOLe
IF (MOLe >= 999999. ) then
  PHIM=1.
ELSE IF (MOLe > 0 .and. MOLe < 999999.) then
  PHIM=1+5.*(zzz/MOLe)
ELSE IF (MOLe < 0 .and. MOLe > -999999.) then
  PHIM=(1-16.*(zzz/MOLe))**(-0.25)
ELSE
  PHIM=1.
ENDIF
return
end function
function PHIH(zzz,MOLe)
real :: PHIH
real :: zzz, MOLe
IF (MOLe >= 999999. ) then
  PHIH=1.
ELSE IF (MOLe > 0 .and. MOLe < 999999. ) then
  PHIH=1+5.*(zzz/MOLe)
ELSE IF (MOLe < 0 .and. MOLe > -999999.) then
  PHIH=(1-16.*(zzz/MOLe))**(-0.50)
ELSE
  PHIH=1.
ENDIF
return
end function
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      function zolri(ri,z,z0)
!
      if (ri.lt.0.)then
        x1=-5.
        x2=0.
      else
        x1=0.
        x2=5.
      endif
!
      fx1=zolri2(x1,ri,z,z0)
      fx2=zolri2(x2,ri,z,z0)
      Do While (abs(x1 - x2) > 0.01)
      if(abs(fx2).lt.abs(fx1))then
        x1=x1-fx1/(fx2-fx1)*(x2-x1)
        fx1=zolri2(x1,ri,z,z0)
        zolri=x1
      else
        x2=x2-fx2/(fx2-fx1)*(x2-x1)
        fx2=zolri2(x2,ri,z,z0)
        zolri=x2
      endif
!
      enddo
!

      return
      end function

!
! -----------------------------------------------------------------------
!
      function zolri2(zol2,ri2,z,z0)
!
      if(zol2*ri2 .lt. 0.)zol2=0.  ! limit zol2 - must be same sign as ri2
!
      zol20=zol2*z0/z ! z0/L
      zol3=zol2+zol20 ! (z+z0)/L
!
      if (ri2.lt.0) then
      psix2=log((z+z0)/z0)-(psim_unstable(zol3)-psim_unstable(zol20))
      psih2=log((z+z0)/z0)-(psih_unstable(zol3)-psih_unstable(zol20))
      else
      psix2=log((z+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
      psih2=log((z+z0)/z0)-(psih_stable(zol3)-psih_stable(zol20))
      endif
!
      zolri2=zol2*psih2/psix2**2-ri2
!
      return
      end function

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_Junhong
      function zolri3(jri,jz,jz0,jzl,Lc,Hc,ust,xka)
real :: zolri3
real :: jri,jz,jz0,jzl, Lc, Hc, ust, xka
real :: jx1, jx2, fjx1, fjx2
real :: D1, beta, prc, ff
Integer :: jjun ! iteration countor

      jjun =0
      if (jri.lt.0.)then
        jx1=-5.
        jx2=0.
      CALL PSI_hat(jzl*jx1/jz,jz0*jx1/jz,jx1,jz,jz, Lc, Hc, & !IN
                  &D1, beta, PSIM_hat, PSIH_hat)           !OUT
      fjx1=zolri4(jx1,jri,jz,jz0,jzl,PSIM_hat,PSIH_hat,ust,xka)
      fjx2=zolri4(jx2,jri,jz,jz0,jzl,0.,0.,ust,xka)
      else
        jx1=0.
        jx2=5.
      CALL PSI_hat(jzl*jx2/jz,jz0*jx2/jz,jx2,jz,jz, Lc, Hc,&  !IN
                  &D1, beta, PSIM_hat, PSIH_hat)           !OUT
      fjx1=zolri4(jx1,jri,jz,jz0,jzl,0.,0.,ust,xka)
      fjx2=zolri4(jx2,jri,jz,jz0,jzl,PSIM_hat,PSIH_hat,ust,xka)
      endif


      Do While (abs(jx1 - jx2) > 0.01)
      if(abs(fjx2).lt.abs(fjx1))then
        jx1=jx1-fjx1/(fjx2-fjx1)*(jx2-jx1)
        CALL PSI_hat(jzl*jx1/jz,jz0*jx1/jz,jx1,jz,jz, Lc, Hc,&  !IN
                    &D1, beta, PSIM_hat, PSIH_hat)           !OUT
        fjx1=zolri4(jx1,jri,jz,jz0,jzl,PSIM_hat,PSIH_hat,ust,xka)
        zolri3=jx1

      else

        jx2=jx2-fjx2/(fjx2-fjx1)*(jx2-jx1)
        CALL PSI_hat(jzl*jx2/jz,jz0*jx2/jz,jx2,jz,jz, Lc, Hc, & !IN
                    &D1, beta, PSIM_hat, PSIH_hat)           !OUT
        fjx2=zolri4(jx2,jri,jz,jz0,jzl,PSIM_hat,PSIH_hat,ust,xka)
        zolri3=jx2
      endif
!
      jjun = jjun + 1
!
      return
      end function
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$_Junhong
      function zolri4(kzol2,kri2,kz,kz0,kzl,PSIM_hat,PSIH_hat,ust,xka)
real :: zolri4
real :: kzol20,kzol2,kri2,kz,kz0,kzl,PSIM_hat,PSIH_hat,ust,xka
real :: kpsix2, kpsih2
real :: temp, kzol3, kzol2l
!
      if(kzol2*kri2 .lt. 0.)kzol2=0.  ! limit zol2 - must be same sign as ri2
!
      kzol20=kzol2*kz0/kz ! z0/L
      kzol3=kzol2+kzol20 ! (z+z0)/L
      kzol2l=kzol2*kzl/kz ! zl/L
!
      if (kri2.lt.0) then
      kpsix2=log(kz/kz0)-(psim_unstable(kzol2)-psim_unstable(kzol20)-PSIM_hat)
      kpsih2=log(0.4*ust*kz/xka  +kz/kzl)-(psih_unstable(kzol2)-psih_unstable(kzol2l)-PSIH_hat)
      else

      kpsix2=log(kz/kz0)-(psim_stable(kzol2)-psim_stable(kzol20)-PSIM_hat)
      kpsih2=log(0.4*ust*kz/xka  +kz/kzl)-(psih_stable(kzol2)-psih_stable(kzol2l)-PSIH_hat)
      endif
!
      zolri4=kzol2*kpsih2/kpsix2**2-kri2
!
      return
      end function
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ... integrated similarity functions ...
!
      function psim_stable_full(zolf)
        psim_stable_full=-6.1*log(zolf+(1+zolf**2.5)**(1./2.5))
      return
      end function

      function psih_stable_full(zolf)
        psih_stable_full=-5.3*log(zolf+(1+zolf**1.1)**(1./1.1))
      return
      end function
      
      function psim_unstable_full(zolf)
        x=(1.-16.*zolf)**.25
        psimk=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*ATAN(1.)
!
        ym=(1.-10.*zolf)**0.33
        psimc=(3./2.)*log((ym**2.+ym+1.)/3.)-sqrt(3.)*ATAN((2.*ym+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
!
        psim_unstable_full=(psimk+zolf**2*(psimc))/(1+zolf**2.)

      return
      end function

      function psih_unstable_full(zolf)
        y=(1.-16.*zolf)**.5
        psihk=2.*log((1+y)/2.)
!
        yh=(1.-34.*zolf)**0.33
        psihc=(3./2.)*log((yh**2.+yh+1.)/3.)-sqrt(3.)*ATAN((2.*yh+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
!
        psih_unstable_full=(psihk+zolf**2*(psihc))/(1+zolf**2.)

      return
      end function

! look-up table functions
      function psim_stable(zolf)
      integer :: nzol
      real    :: rzol
        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .le. 1000)then
           psim_stable = psim_stab(nzol) + rzol*(psim_stab(nzol+1)-psim_stab(nzol))
        else
           psim_stable = psim_stable_full(zolf)
        endif
      return
      end function

      function psih_stable(zolf)
      integer :: nzol
      real    :: rzol
        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .le. 1000)then
           psih_stable = psih_stab(nzol) + rzol*(psih_stab(nzol+1)-psih_stab(nzol))
        else
           psih_stable = psih_stable_full(zolf)
        endif
      return
      end function
      
      function psim_unstable(zolf)
      integer :: nzol
      real    :: rzol
        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .le. 1000)then
           psim_unstable = psim_unstab(nzol) + rzol*(psim_unstab(nzol+1)-psim_unstab(nzol))
        else
           psim_unstable = psim_unstable_full(zolf)
        endif
      return
      end function

      function psih_unstable(zolf)
      integer :: nzol
      real    :: rzol
        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .le. 1000)then
           psih_unstable = psih_unstab(nzol) + rzol*(psih_unstab(nzol+1)-psih_unstab(nzol))
        else
           psih_unstable = psih_unstable_full(zolf)
        endif
      return
      end function

!-------------------------------------------------------------------          

END MODULE module_sf_sfclayrev_rsl

!
! ----------------------------------------------------------
!


