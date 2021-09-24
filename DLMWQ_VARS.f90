	   MODULE DLMWQ_VARIABLES
	   IMPLICIT NONE
!**********************************************************************

      INTEGER*4  MAXPAR
      INTEGER*4  MAXOUT
      INTEGER*4  MAXINF
      INTEGER*4  MAXSTO
      INTEGER*4  MAXNS
      INTEGER*4  MAXDIF
      INTEGER*4  MAXCHL
	   INTEGER*4  MAXNLAYERS
	   INTEGER*4  MAXSEG
      
      REAL*8  VISC

! changed all layer array dimensions 100=>500, 1000=>2000 LAKE TAHOE
      PARAMETER (MAXOUT=10,MAXINF=100,MAXSTO=20000,MAXDIF=45)
      PARAMETER (MAXNS=20000,MAXPAR=500,MAXCHL=7)
	   PARAMETER (MAXNLAYERS=20000)   	  
      PARAMETER (VISC=0.00000114)
	   PARAMETER (maxseg=10)

      REAL*8 DOLEN(MAXNS)
      REAL*8 AREA(MAXNS)
      REAL*8 DEN(MAXNS)
	   REAL*8 DEPF(MAXNS)
      REAL*8 DEPTH(MAXNS)
      REAL*8 DEPTHM(MAXNS)
      REAL*8 VOL(MAXNS)
      REAL*8 VOL1(MAXNS)
      REAL*8 EP(MAXNS)
      REAL*8 TEMP(MAXNS)
      REAL*8 SAL(MAXNS)
      REAL*8 WQUAL(MAXNS,100)
      REAL*8 CF(MAXNS,7)
	   REAL*8 FISH(MAXNS)
	   REAL*8 idepth (0:maxnlayers)
      REAL*8 density(0:maxnlayers)
      REAL*8 tempy  (0:maxnlayers)
	   REAL*8 salty  (0:maxnlayers)
	   REAL*8 vector_thermdepth(200000),vector_day(200000)
	   REAL*8 vector_thermtemp(200000),n_low(1000),n_up(1000)
	   REAL*8 n_low_D(1000),n_up_D(1000), n_up_U(1000)
      INTEGER*4 NS

 !     COMMON /ALL/ AREA,DEN,DEPF,DEPTH,DEPTHM,VOL,VOL1,EP,TEMP,   &
 !        SAL,WQUAL,CF,NS,FISH,DOLEN,vector_thermdepth,vector_day

      REAL*8 BASE
	   REAL*8 CRL
      REAL*8 HLE
      REAL*8 LC
      REAL*8 OLEV(MAXOUT)
      REAL*8 OLEN(MAXOUT)
      REAL*8 OWID(MAXOUT)
      REAL*8 WC
      REAL*8 VCRL
      INTEGER*4 NUMOUT

! wef 24mar05 added CRLNGTH, the spillway length
	   REAL*8 CRLNGTH

!!wef 24mar05      COMMON /BASIN/ CRL,HLE,LC,OLEV,OLEN,OWID,WC,NUMOUT,VCRL,BASE
!      COMMON /BASIN/ CRL,HLE,LC,OLEV,OLEN,OWID,WC,NUMOUT,VCRL,BASE, CRLNGTH

      REAL*8 DIFF(MAXDIF)
      REAL*8 DISS
      REAL*8 EINFF
      REAL*8 H1
      REAL*8 HSIG
      REAL*8 VEL
      REAL*8 WNSQ
      REAL*8 XMOM

      INTEGER*4 NUMDIF
      INTEGER*4 ZOOSW
      INTEGER*4 PARTSW
      INTEGER*4 NMETSW

!      COMMON /ENRDIF/ DIFF,DISS,EINFF,H1,HSIG,VEL,WNSQ,XMOM,NUMDIF,ZOOSW,PARTSW,NMETSW

      REAL*8 HTSAVE

      INTEGER*4 NOSECS,ICLOCK,ITIMES

 !     COMMON /TMSTEP/ HTSAVE,NOSECS,ICLOCK,ITIMES

      REAL*8 DRW(MAXOUT),DELTA_VOL,PREV_DRW(MAXOUT)
      REAL*8 ALPHA(MAXINF,maxseg)
      REAL*8 CDRAG(MAXINF,maxseg)
!wef      REAL*8 PHI(MAXINF)
      REAL*8 seglngth(MAXINF,maxseg)
      REAL*8 bgnwdth(MAXINF,maxseg)
      REAL*8 bgnele(MAXINF,maxseg)
	   REAL*8 plngdpth(MAXINF)
	   REAL*8 FLOINF(MAXINF),PRE_FLOINF(MAXINF)
      REAL*8 TEMINF(MAXINF),PRE_TEMINF(MAXINF)
      REAL*8 SALINF(MAXINF),PRE_SALINF(MAXINF)
      REAL*8 WQINF(MAXINF,100),PRE_WQINF(MAXINF,100)
      REAL*8 CFINF(MAXINF,7),PRE_CFINF(MAXINF,7)
      REAL*8 DDOWN(MAXINF,MAXPAR)
      REAL*8 DIINS(MAXINF,MAXPAR), INDIINS(MAXPAR)
      REAL*8 DLWST(MAXINF)
      REAL*8 DOLD(MAXINF,MAXPAR)
      REAL*8 HFLOW(MAXINF)
      REAL*8 QDOWN(MAXINF,MAXPAR)
      REAL*8 QINS(MAXINF,MAXPAR), INQINS(MAXPAR)
      REAL*8 SDOWN(MAXINF,MAXPAR)
      REAL*8 WQDOWN(MAXINF,100,MAXPAR), INWQ3(MAXPAR), INTIME (MAXPAR)
      REAL*8 CFDOWN(MAXINF,7,MAXPAR)
      REAL*8 SINS(MAXINF,MAXPAR), INSINS(MAXPAR),INCNT(MAXPAR)
      REAL*8 WQINS(MAXINF,100,MAXPAR)
      REAL*8 CFINS(MAXINF,7,MAXPAR)
      REAL*8 TDOWN(MAXINF,MAXPAR),TIMDOWN(MAXINF,MAXPAR)
      REAL*8 TINS(MAXINF,MAXPAR),INTINS(MAXPAR), INPINS(MAXPAR)
      REAL*8 TOTIN(MAXINF),TIMINS(MAXINF,MAXPAR)

      INTEGER*4 ICNT(MAXINF),PARCNT,TOT_CNT
      INTEGER*4 INPAR(MAXINF,MAXPAR),PDOWN(MAXINF,MAXPAR) 
      INTEGER*4 NOINS(MAXINF),PINS(MAXINF,MAXPAR)
      INTEGER*4 NUMINF
	   INTEGER*4 numseg(MAXINF)
      
!	   COMMON /INFLW/ ALPHA,CDRAG,FLOINF,TEMINF,SALINF,WQINF,       &
!     CFINF,DDOWN,QDOWN,SDOWN,WQDOWN,CFDOWN,TDOWN,DOLD,            &
!     ICNT,INPAR,NOINS,NUMINF,HFLOW,DLWST,TOTIN,DIINS,QINS,SINS,   &
!     WQINS,CFINS,TINS,DRW,numseg,seglngth,bgnwdth,bgnele,plngdpth

! Meteorological variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!*
!* RAIN		Rain
!* SRAT		LW incoming (could be 
!* SVPD		Vapor pressure of atmosphere
!* SW		Short wave radiation
!* T4		Air Temperature
!* U6		Wind Speed
!* relhum	relative humidity {decimal}
!* XWIND	Wind multiplicative factor
!* Sunrise  time since midnight that the sun rises {s}
!* Sunset   time the sun sets {s}
!* daylength length fo daylight {s} 
!* RH		Relative humidity {%}
!* HUMIDITY	Humidity switch (0 or 1)
!* zu		Height of anemometer
!* zq		Height of humidity sensor
!* ztair	Height of temperature sensor
!***********************************************************************
      REAL*8 RAIN
      REAL*8 SRAT
      REAL*8 SVPD
!      REAL*8 SVPW			
      REAL*8 SW
      REAL*8 T4
      REAL*8 U6
	   REAL*8 relhum		    
	   REAL*8 XWIND
	   REAL*8 RH
	   REAL*8 U6X
      REAL*8 sunrise	
      REAL*8 sunset
      REAL*8 daylength
	   INTEGER*4 BULK
	   INTEGER*4 HUMIDITY
	   REAL*8 zu, zq, ztair 

!      COMMON /METEO/ RAIN,SRAT,SVPD,SW,T4,U6,relhum,XWIND,RH,U6X,    &
!         sunrise, sunset, daylength, BULK,HUMIDITY,zu,zq,ztair

! Water Quality variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!*
!*
!***********************************************************************
 	   REAL*8 faco(maxns)
      REAL*8 K_SOD
!      new rated definition ....
	   REAL*8 cn1,cn2,cn3,cn4,cn5,biolc
	   REAL*8 cp1,cp2,cp3
      REAL*8 COAG
      REAL*8 CONNH3
      REAL*8 CONSTBO
	   REAL*8 constnh
      REAL*8 CONSTM(MAXCHL)
      REAL*8 CONSTR(MAXCHL)
	   REAL*8 COPDRP
!     zoo feeding rate ....
	   REAL*8 czoo(MAXCHL)
!     Fixed stoichiometry...
	   REAL*8 death(maxns,maxchl)
!     Densy of BOD (not very useful...)
      REAL*8 DENSY(1)
      REAL*8 DENSCF(7)
      REAL*8 ETWAT
	   REAL*8 ETCA(MAXCHL)
  	   REAL*8 ETPART(7)
      REAL*8 ET1(MAXNS)
	   REAL*8 FERED
      REAL*8 FEOXY
!     Fixed stoichiometry
	   REAL*8 growth(maxns,maxchl)
	   REAL*8 grazing(maxns,maxchl)
      REAL*8 GROMAX(MAXCHL)
!     Half Saturation constants
       REAL*8 halfc(maxchl)
	   REAL*8 halfn(maxchl)
	   REAL*8 halfp(maxchl)
      REAL*8 halfsi(maxchl) 
      REAL*8 HSCN(MAXCHL)
      REAL*8 HSCP(MAXCHL)
      REAL*8 HSCSI(MAXCHL)
      REAL*8 HSCZZ
      REAL*8 LIGHT(MAXCHL)	! converted to PAR in fendwq.f
	   REAL*8 MNRED
      REAL*8 MNOXY
	   REAL*8 NINMIN(MAXCHL)
      REAL*8 PINMIN(MAXCHL)
      REAL*8 NINMAX(MAXCHL)
!     for algae conversion...
      REAL*8 phyto_part_fac(MAXCHL)
      REAL*8 org_part_fac
	   REAL*8 paro(maxns)
	   REAL*8 PINMAX(MAXCHL)
      REAL*8 RDSAT
      REAL*8 SATFZ1
      REAL*8 SATFZ2
      REAL*8 sedthp
	   REAL*8 SEDPO4
      REAL*8 SEDNH3
	   REAL*8 SEDNO3
	   REAL*8 SEDMN
      REAL*8 SEDFE
      REAL*8 SEDTEM
      REAL*8 SEGMIN(MAXNS,MAXCHL)
!     for algae settling...
      REAL*8 phyto_setl_vel(MAXCHL)
      REAL*8 org_setl_vel
      REAL*8 STARV1
      REAL*8 STARV2
      REAL*8 THETABO
      REAL*8 THETANH
      REAL*8 THETAT(MAXCHL)
      REAL*8 UNMAX(MAXCHL)
      REAL*8 UPMAX(MAXCHL)
      REAL*8 C_CHLA(MAXNS,MAXCHL)
	   REAL*8 ZRESP
!     Nutrient cholophyll ratios (constant)...
       REAL*8 ycchla(maxchl)
       REAL*8 ypchla(maxchl)
	   REAL*8 ynchla(maxchl)
	   REAL*8 ysichla(maxchl)
!     Oxygen parameters ...
	   REAL*8 KH_OSOD
	   REAL*8 hscdo
	   REAL*8 hscnit
  	   REAL*8 thetabs
!      COMMON /WQUAL/K_SOD,COAG,CONNH3,CONSTBO,constnh,CONSTM,           &
!      CONSTR,COPDRP,DENSY,DENSCF,ETWAT,ETCA,ETPART,ET1,FERED,FEOXY,     &
!      GROMAX,HSCN,HSCP,HSCSI,HSCZZ,LIGHT,MNRED,MNOXY,NINMIN,PINMIN,     &
!      NINMAX,PINMAX,RDSAT,SATFZ1,SATFZ2,sedthp,SEDPO4,SEDNH3,SEDNO3,    &
!      SEDFE,SEDTEM,SEGMIN,STARV1,STARV2,THETABO,THETANH,THETAT,UNMAX,   &
!      UPMAX,C_CHLA,ZRESP,part_fac,setl_vel,faco,paro, cn1,cn2,cn3,cn4,  &
!      cn5,cp1,cp2,cp3,growth,grazing,death, czoo,halfn,halfp,ypchla,    &
!      ynchla,ysichla, KH_OSOD,hscdo,hscnit,thetabs


      REAL*8 B1
      REAL*8 AEM
      REAL*8 AKH
      REAL*8 CK
      REAL*8 CS
      REAL*8 CT
      REAL*8 DB
      REAL*8 DEPMX
	   REAL*8 DEPMX2
      REAL*8 DF
      REAL*8 ETA
      REAL*8 FO
      REAL*8 FSUM
      REAL*8 FTIME
      REAL*8 GPEFF
      REAL*8 H
      REAL*8 OLDSL
      REAL*8 SM
      REAL*8 WQUALM(100)
      REAL*8 CFM(7)
      REAL*8 SPE
      REAL*8 THR
      REAL*8 TI
      REAL*8 TIMEI
      REAL*8 TIMEFI
      REAL*8 TM
      REAL*8 UAV
      REAL*8 UF
      REAL*8 UI
      REAL*8 VF
      REAL*8 VM

      INTEGER*4 J1
      INTEGER*4 K1
      INTEGER*4 MSTEP

      LOGICAL*4 PMIXER

!      COMMON /MIXING/ AEM,AKH,CK,CS,CT,DB,DEPMX,DEPMX2,DF,ETA,FO,FSUM,     &
!         FTIME,GPEFF,H,OLDSL,SM,WQUALM,CFM,SPE,THR,TI,TIMEI,TIMEFI,TM,     &
!         UAV,UF,UI,VF,VM,J1,K1,MSTEP,PMIXER,B1
      REAL*8 A(MAXSTO)
      REAL*8 DADZ(MAXSTO)
      REAL*8 DVV(MAXSTO)
      REAL*8 DVVDA(MAXSTO)
      REAL*8 VV(MAXSTO)
      REAL*8 VVDASH(MAXSTO)

      INTEGER*4 NUMSTO

!      COMMON /TABLE/ A,DVV,VV,DVVDA,VVDASH,DADZ,NUMSTO

      REAL*8 DMAX
      REAL*8 DMIN
      REAL*8 VMAX
      REAL*8 VMIN

!      COMMON /THCK/ DMAX,DMIN,VMAX,VMIN

! Bubbler variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!*
!*KLG BUBFLAG and DIM added March 1992
!*KLG DIM removed June 1992
!*    BUBFLAG	Flag used to determine whether to write insertion data
!*		to the inflow debug file.  If the bubbler is on i.e.
!*		BUBFLAG is true, then the insertion data is not written
!*		to the inflow debug file
!*    DIM	Variable used in the interacting plumes model written
!*		by Dale Robertson
!***********************************************************************

!     COMMON /BUBB/ AFLOW, NPORTS, BDEPTH, IDAYST, IDAYFN,     &
!                  TON, TOFF, BUBLEN, BLEVEL, cdifbeg, BUBFLAG

      REAL*8 AFLOW
      REAL*8 bdepth
      REAL*8 bublen
      REAL*8 ton
      REAL*8 toff
      REAL*8 cdifbeg

      INTEGER*4 nports
      INTEGER*4 idayst
      INTEGER*4 idayfn
      INTEGER*4 BLEVEL

      LOGICAL*4 BUBFLAG


!***********************************************************************
!mech these are the additions for the mechanical mixer
! started 23 Sept 98 by Stephen McCord (SAM)

!	common /mechvars/ mday,xflow,xangle,xdepth,xdiam,xnumber

	   REAL*8 xflow, xdepth, xangle, xdiam		
	   INTEGER*4 mday, xnumber

! Simulation details variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!*
!* JYEAR1	year to start simulation
!* JDAY1	julian day to start simulation
!* NDAYS	number of days to simulate
!* PRODAT	date of inital profile
!* SIMDAY	day of simulation
!*
!***********************************************************************

!	COMMON /SIMDET/ JYEAR1,JDAY1,NDAYS,PRODAT,ITMPR,SIMDAY

	   INTEGER*4 JYEAR1
	   INTEGER*4 JDAY1
	   INTEGER*4 NDAYS
	   INTEGER*4 PRODAT
	   INTEGER*4 ITMPR
	   INTEGER*4 SIMDAY

! Charater variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!*
!* RESNAM	reservoir name
!* RIVNAM	inflowing stream names
!*
!************************************************************************

!	COMMON /NAMES/ RESNAM, RIVNAM

	   CHARACTER*20 RESNAM
	   CHARACTER*10 RIVNAM(MAXINF)
!     DEBUGGING VARIABLES
!	   COMMON /DEBUG/ PBUBB, PENT,z_light
       LOGICAL*4 PBUBB
       LOGICAL*4 PENT
	    REAL*8 z_light !csam for light penetration depth

!     sam for postpro_sum...
!	   COMMON /newcalcs/ light_sum, v_hypo, H_tilda, DO_input
	   REAL*8 light_sum(maxns), v_hypo, H_tilda, DO_input

!   for Error estimates and for calibration of WQ
!     COMMON /GAcalibration/ calibration

	   CHARACTER*3 calibration
!   Zoo variables
      REAL*8 Shit_P, Shit_N, Caca_P, Caca_N
!      COMMON /ZooMysys/ Shit_P, Shit_N, Caca_P, Caca_N

!   Secchi Depth Tahoe Lake Project
      REAL*8 SD,msecchi_dep,A_Tahoe,B_Tahoe,Bsed_Tahoe,Chloro_Tahoe,Kd_Tahoe
	   REAL*8 Parts_Tahoe(7),POP_Tahoe,PON_Tahoe,alg_min(MAXCHL)
	   REAL*8 Nitrate_Tahoe,Ammonia_Tahoe,THP_Tahoe
      REAL*8 fraction,fraction_mean    
!	   COMMON /Secchi/ SD,A_Tahoe,B_Tahoe,Bsed_Tahoe,Chloro_Tahoe,    &
!        Kd_Tahoe,Parts_Tahoe,POP_Tahoe,PON_Tahoe,alg_min,            &
!        Nitrate_Tahoe,Ammonia_Tahoe,THP_Tahoe,fraction,fraction_mean

!     Bulk or Turbulent heat budget formulae switch


!     Output files control variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!* PINSRT	Insertion
!* PHEATR	Heat calculations
!* POUT	    Outflow
!* PLAKE	Lake Number
!* PNUTRIT	Nutrit 
!*
!**********************************************************************

!	COMMON /Switch/ PINSRT,PHEATR,POUT,PLAKE,PWORK,PNUTRIT

	   LOGICAL*4 PINSRT
	   LOGICAL*4 PHEATR
	   LOGICAL*4 POUT
	   LOGICAL*4 PLAKE
	   LOGICAL*4 PWORK
	   LOGICAL*4 PNUTRIT
	   LOGICAL*4 GROW_LIMIT

	   PARAMETER (GROW_LIMIT = .false.)

! Schmidt Work Variables

!***********************************************************************
!*
!* VARIABLE DICTIONARY
!* Z_red	Reduced Depth
!* DepthMax Maximum Depth
!*
!***********************************************************************

	   REAL*8 Z_red
	   REAL*8 DepthMax
	   REAL*8 Shape
!	COMMON /SchmidtWork/ Z_red,DepthMax, Shape
       REAL*8 Bedslope(maxinf)
       REAL*8 TTEMP(20000)
       INTEGER*4 JJDAY
       INTEGER*4 NCHL
       INTEGER*4 JDAY
       
       INTEGER*4 RIVFLO_TIMESTEP,OUTFLO_TIMESTEP,MET_TIMESTEP
       INTEGER*4 RIV_SUBTIMESTEP,OUT_SUBTIMESTEP,MET_SUBTIMESTEP
!GOLOKA
 
!******************************************************************
      END MODULE DLMWQ_VARIABLES
!****************Descriptions ***************************************************	

		 






