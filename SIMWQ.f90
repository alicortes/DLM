!**********************************************************************************
      SUBROUTINE SIMWQ(workfile,HeatDay,HeatFile,LakeFile,						&
      SchmidtFile,SchmidtFile2,NutFile,MixFile,InFile,OutFile,					&
      sename,i_nexp,limname,SDname,red_part,red_TP,red_TN)
!***********************************************************************************
      USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

      REAL*8 CFOLD(maxinf,7)
      REAL*8 CFNEW(maxinf,7)
      REAL*8 DELTNHO(maxns)
      REAL*8 DRWNEW(maxout),DRWOLD(maxout)
      REAL*8 EXTEMP,EXSALT
      REAL*8 EXIFT,EXIFS
      REAL*8 EXCF(7),EXIFCF(7)
      REAL*8 EXWQ(28),EXIFWQ(28)
      REAL*8 HH
      REAL*8 OFLOW1
	  REAL*8 latit
      REAL*8 par
      REAL*8 VOLSUM
      REAL*8 WDL(2)
      REAL*8 xinf(maxinf),xout(maxout)
      REAL*8 sum
      REAL*8 TAECO, TAEWS, TAKETIL, TAKE, SPE_AL
      REAL*8 ZERO,THSND
		REAL*8 AECO, AEWS, AKETIL, AKE
	                   
! 
!     SWITCHES FOR OPTIONAL OUTPUT FILES
!  
      INTEGER*4 ULKN, UFLO, UWOR,UNUT, uwor2, yeari, i_nexp      
		LOGICAL*4 DUNOF,DUNLF
		REAL(8) :: Z_min,solar1, solary,qac,offtake,overflow,conc,D_ANGLE,EQ_TIME,DECL
		CHARACTER*1 lkn, oflo
! 
!     File name declarations
! 
		CHARACTER*8  BITNAME
		CHARACTER*20 WorkFile
		CHARACTER*20 LakeFile,NutFile,MixFile,InFile,OutFile
		CHARACTER*20 SchmidtFile,SchmidtFile2
      CHARACTER*20 Nitrogen,  Phosphoros 
      CHARACTER*20 HeatFile,  HeatDay
      CHARACTER*20 sename(16),SDname
		CHARACTER*20 limname(1)

!  VARIABLE TO CONTROL CREATION OF SIM file
      LOGICAL*4 SWITCH
! 
!     Julain day
! 
		INTEGER*4 MDAYS
! 
!     Cloudiness index
! 
		REAL*8    cloud
! 
!     Heat fluxes variables 
! 
      REAL*8 qsn,XLW,qw,xev,xco,QSNM,QACM,QWM,XEVM,XCOM
		REAL*8 sumaflux,fluxradi,fluxnoradi
      REAL*8 Tsup

! 
!     Night/Day variability
! 
		REAL(8) T_day
		REAL(8) T_night
		REAL(8), PARAMETER	:: Tair_Ratio  = 0.0
! 
!     Sensitivity Analisis parameters
! 
      REAL*8 Relative_U
! 
!     Nutrient Budget Control
! 
		REAL*8 Balance_N,Balance_P
		REAL*8 Balance_N_Stack,Balance_P_Stack
		REAL*8 Balance_Part(7),Balance_Part_Stack(7)
		REAL*8 Part_Stack_F(7),Part_Stack_I(7)
		REAL*8 Nitrogen_T, Phosphorus_T, Part_T(7)
		REAL*8 Nitrogen_T_1, Phosphorus_T_1,Part_T_1(7)
		REAL*8 InLake_Delta_N,InLake_Delta_P,InLake_Delta_Part(7)
		REAL*8 N_Ground,P_Ground
		REAL*8 N_Atm,P_Atm, Part_Atm(7)
		REAL*8 Part_shore(7)
		REAL*8 N_Ins,P_Ins, Part_Ins(7)
		REAL*8 N_Out,P_Out, Part_Out(7)
		REAL*8 N_Ent,P_Ent, Part_Ent(7) 
		REAL*8 N_Inf,P_Inf, Part_Inf(7) 
		REAL*8 Lost_1(9),Lost_2(11)
		REAL*8 Set_Lost_N,Set_Lost_P
		REAL*8 Part_Lost_N,Part_Lost_P, Part_Lost(7)
		REAL*8 Al_THP,Al_NH4,Al_NO3
		REAL*8 Sed_Al_P,Sed_Al_N, Sed_Al_NO3
		REAL*8 Anox_Atm_N
		REAL*8 N_anoxic
!		REAL*8 overflow, offtake

!		Include these parameters IF wanting to know the offtake properties
      REAL*8 EXPROP(3,maxout)
		REAL*8 exprop_WQ(28,maxout)
		REAL*8 tolay(maxout), BOLAY(maxout)
!
      PARAMETER (zero=0.0d0,thsnd=1000.0d0)
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0

      INTEGER*4 i,j,k,kk
      INTEGER*4 ICODE,ISEQ(4),ITMFN
!IF      INTEGER*4 nchl
      INTEGER*4 jstart,JYEAR
      INTEGER*4 LAYNO,lk
      INTEGER*4 ntot
      INTEGER*4 LWIND
      INTEGER*4 zz
! 
!     Switch for addicional output files
! 
      LOGICAL*4 NEWRUN
      LOGICAL*4 NEWFLS
      LOGICAL*4 INIFILEYN

      CHARACTER*1  MOPT
      CHARACTER*1  TKOPT
      INTEGER*4    QUALOPT
      CHARACTER*20 FIN1,FRUN,FMET,FRIV,FOUT,FPAR,FWQD,FLKN,FFLO
      CHARACTER*20 PARTNAME,FAIR,FEFF
! 
!     Time step variable 
! 
      CHARACTER*10 thestamp 
		REAL*8       caca_Part_I(7),caca_Part_F(7)
      INTEGER*4    i_p,j_p
		REAL(8), PARAMETER :: arfac  = 1.0d+6,volfac = 1000.0d0
		REAL(8) AVERAGE_CHLA,SUM_VOL_CHLA, SUM_VOL 

!gbs********River temperature variation************************
      INTEGER*4 totaltimestep, lateflow,delayedtimestep,t3
      REAL*8 increment,decrease,t0,t1,t2
		REAL*8 mintemp(100),maxtemp(100),Avetemp(100),temperature(100,100)
!gbs***********************************************************
		REAL*8 reduction, total_reduction,total_years,outflow 
		REAL*8 red_part,red_TP,red_TN
		REAL*8 satDO(maxns)
!gbs--------------water balance-------------------------
		REAL*8 LAKE_EVA,LAKE_PREC,daily_EVA,daily_PREC,daily_gw
		REAL*8 vol_gw,allflows,ts_flows
		REAL*8 alloutflows,ts_outflow,ts_overflow,alloverflows
      REAL*8 res_storchkB,res_storchkH, res_storchkI, res_storchkO,res_storchkE
!gbs--------------MET DATA-------------------------------
      REAL*8 prev_sw,prev_srat,prev_t4,prev_svpd,prev_rh,prev_U6,prev_rain, prev_msecchi_dep
!GBS------------INFLOW-----------------------------------------------------------
!      REAL*8 PRE_floinf(100), PRE_teminf(100)
!GBS---------------------OUTFLOW----------------
!      REAL*8 PREV_DRW(100)
!gbs*********************************************************
		PBUBB = .FALSE.
10    CONTINUE

		SWITCH = .TRUE.
!  
!     Initialise MENU items
! 
      INIFILEYN = .FALSE.
      NEWRUN    = .FALSE.
      NEWFLS    = .FALSE.
      EINFF     = 0.0	
      ntot      = 0		
      iclock    = 0
      nchl      = 1
      ITMFN     = 86400.0       ! GLOBAL TIME STEP 
      XWIND     = 1.0D0      !Alternate value from subdaily equations
      !XWIND     = 0.87D0 
      ZOOSW     = 0
      PARTSW    = 0
      NMETSW    = 0
! 
!     Outflow flow factor
! 
      DO 20 j=1,maxout
			xout(j)= 1.0D0	!xout = outflow factor.     ALAN PROJECT 
20    CONTINUE
! 
!     Inflow flow factor
!  
      DO 21 j=1,maxinf       
			xinf(j)= 1.0D0	!xinf = inflow factor.      ALAN PROJECT
21    CONTINUE      
!     Initialise bubbler details

      sum=0.0
! 
!     Initialise water quality parameters
! 

! 	Number of difusibles species
      numdif=37	 
! 
!     Get the date part of output files
! 
      CALL timestamp(thestamp) 
 30   CONTINUE
! 
!     SELECT .INI file
! 
		CALL SELFILE(FIN1,PARTNAME,INIFILEYN,WORKFILE)          
! 
!     CREATE the FRUN,FIN1,FPAR file names
! 
		IF (INIFILEYN) THEN
			FRUN = trim(adjustl(partname)) // '.run'
            FMET= trim(adjustl(partname)) // '.CLI'
            FRIV= trim(adjustl(partname)) // '.QIN'
            FOUT= trim(adjustl(partname)) // '.QOT'
			fwqd = trim(adjustl(partname)) // '.wqd'
			fpar = trim(adjustl(partname)) // '.par'

			WRITE(*,*)FPAR,FIN1,FRUN,FMET,FRIV,FOUT,FWQD 
		ENDIF
! 
!     Choose optional units for additional output files
!       
		ulkn  = 30
		uflo  = 31
		unut  = 32
		uwor  = 33
		uwor2 = 34
! 
!     Set to zero the clock for sub-daily loop
! 			
      iclock = 0
      lk     = 0
! 
!  READ DATA from module 4 file
! 
      WRITE(*,50)
50    FORMAT(/1X,'Reading Input DATA',//)

		Nitrogen_T_1   = 0.0D0
		Phosphorus_T_1 = 0.0D0
		
		DO i_p =1,7
			Part_T_1(i_p) = 0.0D0
		ENDDO
! xxxxxxxxxxxxxxx INITIAL PROFILE xxxxxxxxxxxxxxxxxxxxxxxxxxxx
! gbs READ the initial file of the Lake/Reservoir
!  gbs
		OPEN(7,file=FIN1,form='UNFORMATTED',status='OLD')
		READ(7)numout,numinf,numdif,nchl,NUMSTO,ns,ZOOSW,PARTSW,NMETSW,			&
			(TEMP(i),sal(i),DEPTH(i),CF(i,1),CF(i,2),CF(i,3),CF(i,4),				&
			CF(i,5),CF(i,6),CF(i,7),wqual(i,1),wqual(i,2),wqual(i,3),				&
			wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),wqual(i,8),					&
			wqual(i,9),wqual(i,10),wqual(i,11),wqual(i,12),wqual(i,13),				&
			wqual(i,14),wqual(i,15),wqual(i,16),wqual(i,17),wqual(i,18),			&	
			wqual(i,19),wqual(i,20),wqual(i,21),wqual(i,22),wqual(i,23),			&
			wqual(i,24),wqual(i,25),wqual(i,26),wqual(i,27),wqual(i,28),i=1,ns)
      CLOSE(7) 
      DO 60 i=1,ns
			den(i)=densty(temp(i),sal(i))
!	 WRITE(*, fmt='(i, 2f8.1,f15.0)') i, depth(i), temp(i), cf(i,1)
60    CONTINUE
! 
!   READ PARAMETER file
! 
      OPEN(4,file=FPAR,form='UNFORMATTED',status='OLD')

!gbs	READ(4)ndays,LWIND,(DIFF(i),
!gbs     & i=1,numdif),crl,LC,WC,latit,(olev(i),OLEN(i),OWID(i),
!gbs     & i=1,numout),(alpha(i),phi(i),CDRAG(i),i=1,numinf),
!gbs     & CK,ETA,CT,CS,AKH,(VV(i),DVV(i),A(i),DVVDA(i),VVDASH(i),
!gbs     & DADZ(i),i=1,NUMSTO),VMAX,VMIN,VCRL,DMAX,DMIN,
!gbs     & (((WQDOWN(i,KK,j),WQINS(i,KK,j),i=1,numinf),KK=1,28),j=1,maxpar),
!gbs     & (((CFDOWN(i,KK,j),CFINS(i,KK,j),i=1,numinf),KK=1,7),j=1,maxpar),
!gbs     & ((QDOWN(i,j),QINS(i,j),SDOWN(i,j),
!gbs     & SINS(i,j),TDOWN(i,j),TINS(i,j),DDOWN(i,j),DOLD(i,j),DIINS(i,j),
!gbs     & INPAR(i,j),j=1,maxpar),HFLOW(i),TOTIN(i),DLWST(i),ICNT(i),
!gbs     & NOINS(i),i=1,numinf),MSTEP,FSUM,UI,UAV,TIMEI,THR,TI,AEM,FO,
!gbs     & OLDSL,UF,HTSAVE,jstart,RESNAM

		READ(4)NDAYS,LWIND,(DIFF(I),I=1,NUMDIF),CRL,CRLNGTH,LC,WC,LATIT,			&
				(OLEV(I),OLEN(I),OWID(I),I=1,NUMOUT),(numseg(I),I=1,NUMINF),		&
     			((ALPHA(I,k),seglngth(I,k),bgnwdth(I,k),bgnele(I,k),					&
				CDRAG(I,k),I=1,NUMINF),k=1,10),CK,ETA,CT,CS,AKH,						&
				(VV(I),DVV(I),A(I),DVVDA(I),VVDASH(I),DADZ(I),I=1,NUMSTO),			&
				VMAX,VMIN,VCRL,DMAX,DMIN,(((WQDOWN(I,KK,J),								&
				WQINS(I,KK,J),i=1,NUMINF),KK=1,28),j=1,5),(((CFDOWN(I,KK,J),		&
				CFINS(I,KK,J),i=1,NUMINF),KK=1,7),j=1,5),((QDOWN(I,J),QINS(I,J),	&
				SDOWN(I,J),SINS(I,J),TDOWN(I,J),TINS(I,J),DDOWN(I,J),DOLD(I,J),	&
				DIINS(I,J),INPAR(I,J),J=1,MAXPAR),HFLOW(I),TOTIN(I),DLWST(I),		&
				ICNT(I),NOINS(I),I=1,NUMINF),MSTEP,FSUM,UI,UAV,TIMEI,THR,TI,		&
				AEM,FO,OLDSL,UF,HTSAVE,JSTART,RESNAM
      CLOSE(4)

102   CONTINUE

! 
!gbs New WQ definitions
! 
		OPEN(6,file=FWQD,form='UNFORMATTED',status='OLD')

      READ(6)(gromax(i),constr(i),constm(i),thetat(i),light(i),					    &
			ycchla(i),ynchla(i),ypchla(i),ysichla(i),halfc(i),halfn(i),halfp(i),	&
			halfsi(i),hscn(i),czoo(i),etca(i),i=1,nchl),densy(1),etwat,				&
			K_SOD,KH_OSOD,constbo,constnh,cn1,cn2,cn3,cn4,cn5,cp1,cp2,cp3,			&
			thetabo,thetanh,thetabs,sedthp,sednh3,sedno3,sedtem,					&
			coag,(denscf(i),etpart(i),i=1,7),(phyto_setl_vel(i),phyto_part_fac(i),i=1,nchl), &
			org_setl_vel,org_part_fac 
      CLOSE(6)
!-----------------------------------------------------------------------
!--------------FILES OPENED FOR RESULTS AND DEBUGGINH-------------------
!-----------------------------------------------------------------------	
		IF(PMIXER)OPEN(13,file = MixFile,form='FORMATTED',status='UNKNOWN')
		IF(PINSRT)OPEN(18,file = InFile,form='FORMATTED',status='UNKNOWN')
		IF(PHEATR) THEN
			OPEN(19,file = HeatFile,form='FORMATTED',status='UNKNOWN',access='append')
			OPEN(22,file = HeatDay,status='unknown',access='append')      
		ENDIF  
		IF(POUT) OPEN(UFLO,file = OutFile,form='FORMATTED',status='UNKNOWN')
      IF(PLAKE)OPEN(ulkn,file = LakeFile,form='formatted',status='unknown')
      IF(PWORK) THEN
			OPEN(uwor,file = SchmidtFile,form='formatted',status='unknown')
			OPEN(uwor2,file = SchmidtFile2,form='formatted',status='unknown')
		ENDIF 
        	       
      IF(PNUTRIT) OPEN(unut,file = NutFile,form='formatted',status='unknown')	       
! 
!   WRITE title pages of extra output files IF requested
! 
      IF (PMIXER) WRITE (13,540)
540   FORMAT (2X,'Year',2X,'jday',5X,'AECOOLING',3X,'AEWINDSPEED',8X,'AKETIL',2X,'AVAILABLE KE')    
      IF(PHEATR) WRITE(19,550) !RESNAM,jday
!550   FORMAT(' HEATR file FOR ',A20,' RUN ON ',10X,' AT ',5X,
!     &' FOR DAY ',I10,/,
!     &'iclock',8X,'LWnet',8X,'EVA',4X,'CON',8X,'LWout',4X,
!     &'LWin',4X,'SW',3X,'LWnet+sw',5X,'sum',5x,'Tsup',3x,
!     &'HeatConten (Mj)')
550   FORMAT('jday',1x,'Clock',5x,'NetLWR'6X,'EVAPO',5X,'CONDUC',6X,'LWout',7X,	&
			'LWin',8X,'SWR',2X,'LWnet+SWR',4X,'NETHEAT',1X,'Airtemp',1X,'SurTemp',1X,'BotTemp',  &
			6X,'Ea',6X,'Es',2X,'eva_loss')
!      FORMAT(I4,1X,I5,1X,8(i6,2X),2f8.2)
!-------written in INSERT SUBROUTINE line 2879 see DCALCS.for
      IF(PINSRT) THEN 
			WRITE(18,560) RESNAM
560		FORMAT(/,'(INSERT SUBROUTINE, INSERTION file for ',A20,5X,/)
			WRITE(18,790)
790		FORMAT('Jday',3x,'timestep',4x,'River #',6x,'Layer',8X,'Depth',5X,'Volume 10^3m3')
		ENDIF
!----------------------------------------------------------------------------------------------OUTPUT FILE (UFLO)
		IF(POUT) WRITE(UFLO,570) RESNAM
570   FORMAT (1X,A35,7X,'---------------- OUTLET #1 ---------------',			&
			/,'jday',8X,'OUT VOLUME   TEMP    sal','   TOP  BOTTOM',					&
		' CHLA1	CHLA2	CHLA3	DOtot	BOD	THP	Pin1	Pin2	POP	RP ',				&
		' NO3	NH3   Nin1	Nin2	PON	  DON   Silica'									&
		'Part1   Part2    Part3    Part4   Part5    Part6    Part7')

		IF(PLAKE) WRITE(ulkn,580) RESNAM
580   FORMAT (1x,A35,/,4x,'jday',3x,'LAKE NUMBER',2x, 'Mixing_Eleva',			&
            'ThermEleva',2x,'Dev_Low_Eleva',4x'Dev_Up_Eleva',4x,					&
            'ThermTemp',4x,'Dev_Low_Temp',x,'Dev_Up_Temp',x)

		IF(PWORK) WRITE(uwor,590)
!590	FORMAT('     Jday', x,'  Mix_Elev',x,
!     &'Therm_El_F',x,'Therm_El_C', 5x,
!     &'Schmidt_W',5x,'Total_W')! ,'Heat_Content')
590	FORMAT(4X,'Jday', 2x,'Mix_Elev',x,'T_average',x,'T_ave_vol',4x,'Schmidt_W',	&
				7x,'Total_W',x,'T_ave_vol_LT100',x,'T_ave_vol_GT100')      ! ,'Heat_Content')

     
      IF(PNUTRIT) WRITE(unut,595) 
595   FORMAT('Units Kg per day  ',/,														&            
			2X,'Year',2X,'jday',07X,'Balance_N',01X,'Balance_Stack_N',04X,			&
			'Acumulated_N',06X,'Nitrogen_T',11X,'N_Out',11X,'N_Atm',11X,			&
			'N_Ins',11X,'N_Inf',08X,'N_Ground',07X,'N_Settled',06X,					&
			'N_Sed_Rele',06X,'N_Anox_Atm',07X,'Balance_P',01X,							&
			'Balance_Stack_P',04X,'Acumulated_P',04X,'Phosphorus_T',11X,			&
			'P_Out',11X,'P_Atm',11X,'P_Ins',11X,'P_Inf',08X,'P_Ground',07X,		&
			'P_Settled',06X,'P_Released')
        
	  
!	  FORMAT(2i6,1X,23(F15.2,1X))
!---------------SECCHI DEPTH -------------------------------------------	
		OPEN(27,file= 'secchi.txt',status='unknown',access='append')
		WRITE(27,596)
!	596	FORMAT(2X,'Year',2X,'Jday', 5x,'Secchi(m)')
596	FORMAT(2X,'Jday', 5x,'Secchi(m)')
		CLOSE(27)
!---------------SEDIMENT RELEASE FILE-----------------------------------
		OPEN(28,file= 'sediments_NP.txt',status='unknown',access='append')
		WRITE(28,597)
597	FORMAT(2X,'Year',2x,'Jday', 13x,'SRP(mg)',13x,'NH4(mg)',13x,'NO3(mg)')
		CLOSE(28)
!--------------WATER BALANCE DAILY----------------------------------------------
		OPEN(29,file='water_balance.txt',status='unknown',access='append')
		WRITE(29,598)
598	FORMAT(2x,'Jday',10x, 'Stream(m3)',14x, 'GroundWater(m3)',12x,			&
       'RAIN(m3)',10x,'EVAPOR(m3)',8X,'Outflows(m3)',7X,'Overflows(m3)',7x,'Lake_surface_area(m2)')
		CLOSE(29)
!---------------WATER SURFACE ELEVATION---------------------------------
		OPEN(73,FILE="surface.txt",STATUS='UNKNOWN',access='append')
		WRITE(73,599)
599	FORMAT(2X,'Jday',4x, 'WSE(m)')
		CLOSE(73)
!------------------------------------------------------------------------
!--------------WATER BALANCE- TIME STEP---------------------------------------------
		OPEN(21,file='water_balance_timestep.txt',status='unknown',access='append')
		WRITE(21,601)
		
601	FORMAT (4X,'JDAY',3X,'ICLOCK',1X,'res_storchkB',3X,'-LAKE_EVA',3X,'depth(ns)', &
        1X,'res_storchkH',4X,'ts_flows',1X, 'res_storchkI', 3X,'depth(ns)',&
       3X,'-offtake',3X,'-overflow',1X,'res_storchkO',3X,'depth(ns)',1X,'res_storchkE')       

		CLOSE(21)
!------------------------------------------------------------------------------------		
!  Calculate cumulative layer volumes from layer depths

		ICODE=1
		LAYNO=1

		CALL RESINT(ICODE,LAYNO)
!       
!   Check layers for vmax,vmin
!  
      CALL THICK

      DEPMX=DEPTH(ns-1)
! 
!     Check Nutrient Balance. Fill with initial conditions
! 
		DO i = 1,ns
			Nitrogen_T_1 = Nitrogen_T_1 + (wqual(i,21) + wqual(i,22) + wqual(i,23) +	&
     				 wqual(i,24) + wqual(i,25) + wqual(i,26))*vol(i)

			Phosphorus_T_1 = Phosphorus_T_1 +(wqual(i,15)+ wqual(i,16) + wqual(i,17) + &
			wqual(i,18) + wqual(i,19) + wqual(i,20))*vol(i)
	 
			DO i_p = 1,7
				Part_T_1(i_p) = Part_T_1(i_p) + cf(i,i_p)*vol(i)*1000.0
			ENDDO  ! i_p
      ENDDO  ! i 

      caca_part_I = 0.0

      DO j_p = 1,ns
			DO i_p =1,7
				caca_part_I(i_p) = caca_part_I(i_p) + cf(j_p,i_p)*vol(j_p)
			ENDDO
		ENDDO

      OPEN(8,file=FRUN,form='UNFORMATTED',status='OLD')
      OPEN(9011,file=FRIV,form='UNFORMATTED',status='OLD')
      OPEN(9012,file=FOUT,form='UNFORMATTED',status='OLD')
      OPEN(9013,file=FMET,form='UNFORMATTED',status='OLD')

! Convert latitude from degrees to radians

!		latit=2*PI+latit*PI/180   !DAVIS, SEE DUFFIE and BECKMAN,1990
		latit=latit*PI/180	    

!***************************************************************************************
!************************DAILY LOOP COMMENCES*******************************************
!***************************************************************************************
460   lk		= lk+1
      ntot	= ntot+1
      AECO   = 0.0d0
      AEWS   = 0.0d0
      AKETIL = 0.0d0
		AKE    = 0.0d0
! 
!     Nutrient Budget inicialize
! 
		Balance_N_Stack	= 0.0D0
		Balance_P_Stack	= 0.0D0
		Nitrogen_T			= 0.0D0
		Phosphorus_T		= 0.0D0
		N_Ground			= 0.0D0
		P_Ground			= 0.0D0
		N_Atm				= 0.0D0 
		P_Atm				= 0.0D0
		N_Ins				= 0.0D0
		P_Ins				= 0.0D0
		N_Out				= 0.0D0
		P_Out				= 0.0D0
		N_Ent				= 0.0D0
		P_Ent				= 0.0D0
		N_Inf				= 0.0D0
		P_Inf				= 0.0D0

		DO i_p = 1,7
			Part_Inf(i_p)           = 0.0D0
			Part_Ent(i_p)           = 0.0D0
			Balance_Part(i_p)       = 0.0D0
			Balance_Part_Stack(i_p) = 0.0D0
			Part_T(i_p)             = 0.0D0
			Part_Ins(i_p)           = 0.0D0
			Part_Out(i_p)           = 0.0D0
			Part_Lost(i_p)          = 0.0D0
			Part_Atm(i_p)			 = 0.0D0
			Part_shore(i_p)         = 0.0D0			 
		ENDDO 

		Lost_1 = 0.0d0
		Lost_2 = 0.0d0
      
!		caca_part_I = 0.0
		caca_part_F = 0.0d0
		Set_Lost_N  = 0.0D0
		Part_Lost_N = 0.0D0
		Set_Lost_P  = 0.0D0
		Part_Lost_P = 0.0D0
		Sed_Al_P    = 0.0D0
		Sed_Al_N    = 0.0D0
		Sed_Al_NO3	= 0.0D0
		Anox_Atm_N  = 0.0D0
! 
!     Set to zero daily heat fluxes
! 
		QSNM       = 0.0d0
   	    QACM       = 0.0d0
		QWM        = 0.0d0
		XEVM       = 0.0d0
		XCOM       = 0.0d0
		fluxradi   = 0.0d0
     	fluxnoradi = 0.0d0
		sumaflux   = 0.0d0
!gbs  water balance terms
		daily_EVA  = 0.0d0
		daily_PREC = 0.0d0
		daily_gw   = 0.0d0
		allflows   = 0.0d0
		alloutflows= 0.0d0
		alloverflows=0.0d0
! 
!     Set to zero daily heat fluxes
! 
		QSNM       = 0.0d0
    	QACM       = 0.0d0
		QWM        = 0.0d0
		XEVM       = 0.0d0
		XCOM       = 0.0d0
		fluxradi   = 0.0d0
     	fluxnoradi = 0.0d0
		sumaflux   = 0.0d0

! 
!gbs READ run-time DATA all STREAM INPUT DATA ON DAILY, outflow and other from file *.INI. 
!
!--------------OLD----------------------------------------------------------------
!2015/9      IF(humidity.eq.1) THEN
!2015/9			READ(8)(ISEQ(i),i=1,2),(floinf(i),teminf(i),salinf(i),				&
!2015/9			wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),wqinf(i,5),				&
!2015/9			wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),wqinf(i,10),			&
!2015/9			wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),					&
!2015/9			wqinf(i,15),wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),		&
!2015/9			wqinf(i,20),wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),		&
!2015/9			wqinf(i,25),wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),			&
!2015/9			cfinf(i,2),cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),				&
!2015/9			cfinf(i,7),i=1,numinf),(drw(i),i=1,numout),sw,srat,t4,rh,U6,rain
!2015/9      ELSE
!2015/9			READ(8)(ISEQ(i),i=1,2),(floinf(i),teminf(i),salinf(i),wqinf(i,1),	&
!2015/9			wqinf(i,2),wqinf(i,3),wqinf(i,4),wqinf(i,5),wqinf(i,6),				&
!2015/9			wqinf(i,7),wqinf(i,8),wqinf(i,9),wqinf(i,10),wqinf(i,11),			&
!2015/9			wqinf(i,12),wqinf(i,13),wqinf(i,14),wqinf(i,15),wqinf(i,16),		&
!2015/9			wqinf(i,17),wqinf(i,18),wqinf(i,19),wqinf(i,20),wqinf(i,21),		&
!2015/9			wqinf(i,22),wqinf(i,23),wqinf(i,24),wqinf(i,25),wqinf(i,26),		&
!2015/9			wqinf(i,27),wqinf(i,28),cfinf(i,1),cfinf(i,2),cfinf(i,3),			&
!2015/9			cfinf(i,4),cfinf(i,5),cfinf(i,6),cfinf(i,7),i=1,numinf),				&
!2015/9			(drw(i),i=1,numout),sw,srat,t4,svpd,U6,rain 
!2015/9      ENDIF
!-------------------------END OLD------------------------------------------------------------      
!------------READ DAILY RIVERS' INFLOWS----------------------------------------------
      IF(RIVFLO_TIMESTEP.EQ.86400) THEN
            READ(9011)(ISEQ(1),ISEQ(2),ISEQ(3),ISEQ(4),floinf(i),teminf(i),		&
			salinf(i),wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),				&
			wqinf(i,5),wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),				&
			wqinf(i,10),wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),		&
			wqinf(i,15),wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),		&
			wqinf(i,20),wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),		&
			wqinf(i,25),wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),			&
			cfinf(i,2),cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),				&
			cfinf(i,7),i=1,numinf)		
      ENDIF      
!------------READ DAILY RESERVOIR OUTFLOWS----------------------------------------------      
      IF(OUTFLO_TIMESTEP.EQ.86400) THEN
            READ(9012)(ISEQ(i),i=1,4),(drw(i),i=1,numout)             
      ENDIF
!      PAUSE
!------------READ DAILY METEOROLOGICAL DATA------------------------------------------------------------    
    IF(MET_TIMESTEP.EQ.86400) THEN
      IF(humidity.eq.1) THEN
			READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,rh,U6,rain,msecchi_dep
      ELSE
			READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,svpd,U6,rain,msecchi_dep 
      ENDIF    
      !     Put in wind factor common for all routines	
!    	     
		U6X  = U6*XWIND

!	Relative_U = 0.15*sw
!	sw = sw + Relative_U
! 
!     iseq(1) = year, iseq(2) = julian day		
! 
      jday=iseq(1)*1000+iseq(2)
		JYEAR=ISEQ(1) 
! wef 29jan05 added a running count of the simulation day to make it easier to WRITE debug output files
		SIMDAY = ISEQ(2)
! wef 29jan05	

! 
!     Astronomical Calculations
!
!  
!		sct 2Aug18		replaced with 1972 TVA equations 
!		solar1 = 23.45*sin((284+iseq(2))*2.0d0*pi/365.0d0)
!		solar1 = solar1*pi/180
!		solary = 2.0*acos(-tan(latit)*tan(solar1))
!		daylength =  3600.0d0 * solary/(15.0D0*pi/180.0d0)
!		sunrise   = 43200.0d0 - daylength/2.0d0
!		sunset    = 43200.0d0 + daylength/2.0d0

	D_ANGLE= 2.0*pi* ((ISEQ(2))-1.0)/365.242		
	EQ_TIME= -(0.12357*sin(D_ANGLE)-0.004289*cos(D_ANGLE) +0.153809*sin(2.0*D_ANGLE)+0.060783*cos(2.0*D_ANGLE))
    DECL= asin(sin((23.445*pi/180)) * sin((2.0*pi/360.0) * ( 279.9348 + (D_ANGLE*360.0/2.0/pi) + (1.914827*sin(D_ANGLE)) - (0.079525*cos(D_ANGLE)) + (0.019938*sin(2*D_ANGLE)) - (0.001620*cos(2*D_ANGLE)))))
    sunset=(((12)/pi) * acos((sin(0.0)-sin(latit)*sin(DECL))/(cos(latit)*cos(DECL)))+0.0-EQ_TIME+(12))*3600 
    sunrise= (sunset - (acos((sin(0.0)-sin(latit)*sin(DECL))/(cos(latit)*cos(DECL)))*(24)/pi)*3600)     
	daylength=sunset-sunrise 

				
!		IF(t4.lt.0.0d0) Tair_Ratio = - Tair_Ratio
!		T_night   = (86400/(86400 +daylength*(Tair_Ratio-1)))*t4
!		T_day     = Tair_Ratio*T_night

		T_night   = (86400.0d0*t4 - daylength*Tair_Ratio)/86400.0d0
		T_day     = Tair_Ratio + T_night
		IF(itimes .ge. sunrise .and. itimes.le. sunset) THEN
			T4 = T_day
		ELSE
			T4 = T_night
		ENDIF 
! 
!  Convert inputs to rates - leave inflow and withdrawal
!  IF cloud cover is used to calculate LONG WAVE IN THEN 
!  don't convert to a rate
!     Change units from kj/m2day to W/m2 

      IF(LWIND.eq.2.OR.LWIND.eq.3) srat = srat*1000.0d0/86400.0D0

    ENDIF    
!----------------END READING DAILY TIME SERIES DATA-----------------------


!  Start forcing-mixing-diffusion loop. Index is iclock      

		t0=pi   !SUB-DAILY INFLOW VARIATION
		t1=pi   !SUB-DAILY INFLOW VARIATION
		t3=0    !SUB-DAILY INFLOW VARIATION

600   CONTINUE         ! SUB-DAILY TIMESTEP START
!	print*,wqual(ns-30,10),wqual(ns-30,7), 'step 1'
! 
!      Thermal transfers are done by heatr
! 
!*****************************************************************************************
!*************************START SUB-DAILY TIMESTEP****************************************
!*****************************************************************************************
! 
!      Set initial time step 
! 
! 		IF(itimes.eq.1440) THEN   !variable time step
! 			nosecs=900
!     ELSE
			nosecs=itimes*60         !fixed time step
!     ENDIF
        res_storchkB=0.0d0
        res_storchkH=0.0d0
        res_storchkI=0.0d0
        res_storchkO=0.0d0
        res_storchkE=0.0d0
        do i=1,ns
            res_storchkB=res_storchkB+vol(i)            
        enddo
       
!------------READ SUB_DAILY RIVERS' INFLOWS----------------------------------------------
      IF(RIVFLO_TIMESTEP.LT.86400) THEN
         IF(MOD(ICLOCK,RIVFLO_TIMESTEP).EQ.0)THEN
            READ(9011)(ISEQ(1),ISEQ(2),ISEQ(3),ISEQ(4),floinf(i),teminf(i),		&
			salinf(i),wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),				&
			wqinf(i,5),wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),				&
			wqinf(i,10),wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),		&
			wqinf(i,15),wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),		&
			wqinf(i,20),wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),		&
			wqinf(i,25),wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),			&
			cfinf(i,2),cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),				&
			cfinf(i,7),i=1,numinf)
               DO I=1,numinf
                   PRE_floinf(i)=floinf(i)
                   PRE_teminf(i)=teminf(i)
                   PRE_salinf(i)=salinf(i)
                   DO j=1,28
                       PRE_wqinf(i,j)=wqinf(i,j)
                   ENDDO
                   DO k=1,7
                       PRE_cfinf(i,k)=cfinf(i,k)
                   ENDDO                   
               ENDDO

         ELSEIF(MOD(ICLOCK,RIVFLO_TIMESTEP).NE.0)THEN
            READ(9011)(ISEQ(1),ISEQ(2),ISEQ(3),ISEQ(4),floinf(i),teminf(i),		&
			salinf(i),wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),				&
			wqinf(i,5),wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),				&
			wqinf(i,10),wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),		&
			wqinf(i,15),wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),		&
			wqinf(i,20),wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),		&
			wqinf(i,25),wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),			&
			cfinf(i,2),cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),				&
			cfinf(i,7),i=1,numinf)
               BACKSPACE 9011
               DO I=1,numinf
                   floinf(i)=0.5*(PRE_floinf(i)+floinf(i))
                   teminf(i)=0.5*(PRE_teminf(i)+teminf(i))
                   salinf(i)=0.5*(PRE_salinf(i)+salinf(i))
                   DO j=1,28
                       wqinf(i,j)=0.5*(PRE_wqinf(i,j)+wqinf(i,j))
                   ENDDO
                   DO k=1,7
                      cfinf(i,k)=0.5*(PRE_cfinf(i,k)+cfinf(i,k))
                   ENDDO                   
               ENDDO               
			ENDIF			
      ENDIF
!------------READ SUB_DAILY RESERVOIR OUTFLOWS----------------------------------------------      
      IF(OUTFLO_TIMESTEP.LT.86400) THEN
         IF(MOD(ICLOCK,OUTFLO_TIMESTEP).EQ.0)THEN
            READ(9012)(ISEQ(i),i=1,4),(drw(i),i=1,numout)
            DO I=1,NUMOUT
                PREV_DRW(I)=DRW(I)
            ENDDO
         ELSEIF(MOD(ICLOCK,OUTFLO_TIMESTEP).NE.0)THEN
             READ(9012)(ISEQ(i),i=1,4),(drw(i),i=1,numout)
             BACKSPACE 9012
             DO I=1,NUMOUT
                DRW(I)=0.5*(PREV_DRW(I)+DRW(I))
             ENDDO
         ENDIF            
      ENDIF
!------------READ SUB_DAILY METEOROLOGICAL DATA------------------------------------------------------------    
    IF(MET_TIMESTEP.LT.86400) THEN 
!    PRINT*,'GBS',ICLOCK     
      IF(MOD(ICLOCK,MET_TIMESTEP).EQ.0)THEN
!      PRINT*,ICLOCK,MOD(ICLOCK,MET_TIMESTEP)
         IF(humidity.eq.1) THEN
			   READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,rh,U6,rain,msecchi_dep
         ELSE
			   READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,svpd,U6,rain,msecchi_dep
         ENDIF		 
         prev_sw = sw
         prev_srat = srat
         prev_t4 = t4
         prev_svpd = svpd
         prev_rh = rh
         prev_U6 = U6
         prev_rain = rain
         prev_msecchi_dep = msecchi_dep

!         print*,'met',(ISEQ(i),i=1,4),sw,srat,t4
       ELSEIF(MOD(ICLOCK,MET_TIMESTEP).NE.0)THEN
         IF(humidity.eq.1) THEN
			   READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,rh,U6,rain,msecchi_dep
         ELSE
			   READ(9013)(ISEQ(i),i=1,4),sw,srat,t4,svpd,U6,rain,msecchi_dep
         ENDIF
!         print*,'met1',(ISEQ(i),i=1,4),sw,srat,t4
         backspace 9013
          sw   = (prev_sw   + sw)  *0.5
          srat = (prev_srat + srat)*0.5
          t4   = (prev_t4   + t4)  *0.5
          svpd = (prev_svpd + svpd)*0.5
          rh   = (prev_rh   + rh)  *0.5
          U6   = (prev_U6   + U6)  *0.5
          rain = (prev_rain + rain)*0.5
          msecchi_dep =(prev_msecchi_dep+msecchi_dep)*0.5
 !        print*,'met2',(ISEQ(i),i=1,4),sw,srat,t4
       ENDIF 
       IF(ICLOCK.EQ.84600) ISEQ(2)=ISEQ(2)-1 
!      PAUSE  
!     Put in wind factor common for all routines	
!    	
!      IF(U6.LT.0.00001) U6=0.001
!       IF(U6.gt.2.5) U6=2.5
!		U6X  = U6*XWIND
      U6X  = U6*XWIND
!	Relative_U = 0.15*sw
!	sw = sw + Relative_U
! 
!     iseq(1) = year, iseq(2) = julian day		
! 
      jday=iseq(1)*1000+iseq(2)
		JYEAR=ISEQ(1) 
! wef 29jan05 added a running count of the simulation day to make it easier to WRITE debug output files
		SIMDAY = ISEQ(2)
! wef 29jan05	
!print*,jday, jyear, iseq(1), iseq(2)
! 
!     Astronomical Calculations
! 

		
!		sct 2Aug18		replaced with 1972 TVA equations  flag for recheck
!		solar1 = 23.45*sin((284+iseq(2))*2.0d0*pi/365.0d0)
!		solar1 = solar1*pi/180
!		solary = 2.0*acos(-tan(latit)*tan(solar1))
!		daylength =  3600.0d0 * solary/(15.0D0*pi/180.0d0)
!		sunrise   = 43200.0d0 - daylength/2.0d0
!		sunset    = 43200.0d0 + daylength/2.0d0

		! sunrise = (acos((sin(0.05)-sin(40.0629*pi/180.0)*sin(asin(sin(23.445*pi/180.0)*sin(2.0*pi/360.0* &
		! (279.9348+360.0*(2.0*pi/365.242*(ISEQ(2)-1.0))/2.0/pi+1.914827*sin((2.0*pi/365.242*(ISEQ(2)-1.0)))- &
		! 0.079525*cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.019938*sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))- &
		! 0.001620*cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0))))))))/cos(40.0629*pi/180.0)* &
		! cos(asin(sin(23.445*pi/180.0)*sin(2.0*pi/360.0*(279.9348+360.0*(2.0*pi/365.242*(ISEQ(2)-1.0))/2.0/pi &
		! +1.914827*sin((2.0*pi/365.242*(ISEQ(2)-1.0)))-0.079525*cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.019938* &
		! sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))-0.001620*cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0))))))))*12.0/pi- &
		! 12.0+(120.0-(119.0+34.0/60.0+16.48/3600.0))/15.0+(0.123570*sin((2.0*pi/365.242*(ISEQ(2)-1.0)))-0.004289* &
		! cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.153809*sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))+0.060783* &
		! cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))))*60.0*60.0

		! sunset = (acos((sin(0.075)-sin(40.0629*pi/180.0)*sin(asin(sin(23.445*pi/180.0)*sin(2.0*pi/360.0* &
		! (279.9348+360.0*(2.0*pi/365.242*(ISEQ(2)-1.0))/2.0/pi+1.914827*sin((2.0*pi/365.242*(ISEQ(2)-1.0)))- &
		! 0.079525*cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.019938*sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))-0.001620* &
		! cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0))))))))/cos(40.0629*pi/180.0)*cos(asin(sin(23.445*pi/180.0)* &
		! sin(2.0*pi/360.0*(279.9348+360.0*(2.0*pi/365.242*(ISEQ(2)-1.0))/2.0/pi+1.914827* &
		! sin((2.0*pi/365.242*(ISEQ(2)-1.0)))-0.079525*cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.019938* &
		! sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))-0.001620*cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0))))))))*12.0/pi+12.0+ &
		! (120.0-(119.0+34.0/60.0+16.48/3600.0))/15.0+(0.123570*sin((2.0*pi/365.242*(ISEQ(2)-1.0)))-0.004289* &
		! cos((2.0*pi/365.242*(ISEQ(2)-1.0)))+0.153809*sin(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))+0.060783* &
		! cos(2.0*(2.0*pi/365.242*(ISEQ(2)-1.0)))))*60.0*60.0
				
		! daylength =  sunset-sunrise		
				
	D_ANGLE= 2.0*pi* ((ISEQ(2))-1.0)/365.242		
	EQ_TIME= -(0.12357*sin(D_ANGLE)-0.004289*cos(D_ANGLE) +0.153809*sin(2.0*D_ANGLE)+0.060783*cos(2.0*D_ANGLE))
    DECL= asin(sin((23.445*pi/180)) * sin((2.0*pi/360.0) * ( 279.9348 + (D_ANGLE*360.0/2.0/pi) + (1.914827*sin(D_ANGLE)) - (0.079525*cos(D_ANGLE)) + (0.019938*sin(2*D_ANGLE)) - (0.001620*cos(2*D_ANGLE)))))
    sunset=(((12)/pi) * acos((sin(0.0)-sin(latit)*sin(DECL))/(cos(latit)*cos(DECL)))+0.0-EQ_TIME+(12))*3600 
    sunrise= (sunset - (acos((sin(0.0)-sin(latit)*sin(DECL))/(cos(latit)*cos(DECL)))*(24)/pi)*3600)     
	daylength=sunset-sunrise

		
!		IF(t4.lt.0.0d0) Tair_Ratio = - Tair_Ratio ! flag for recheck sct
!		T_night   = (86400/(86400 +daylength*(Tair_Ratio-1)))*t4
!		T_day     = Tair_Ratio*T_night

		T_night   = (86400.0d0*t4 - daylength*Tair_Ratio)/86400.0d0
		T_day     = Tair_Ratio + T_night
		IF(itimes .ge. sunrise .and. itimes.le. sunset) THEN
			T4 = T_day
		ELSE
			T4 = T_night
		ENDIF 
! 
!  Convert inputs to rates - leave inflow and withdrawal
!  IF cloud cover is used to calculate LONG WAVE IN THEN !  don't convert to a rate
      WRITE(996,FMT='(2I8,8F15.6)') JDAY,ICLOCK, FLOAT(ICLOCK)/3600.0,sw,srat,t4,rh,U6,rain,msecchi_dep
      IF(LWIND.eq.2.OR.LWIND.eq.3) srat = srat*1000.0d0/REAL(MET_TIMESTEP)           !     Change units from kj/m2day to W/m2 
    ENDIF
!----------------END READING SUB_DAILY TIME SERIES DATA-----------------------
! 	print*,iclock,nosecs
!	pause
      CALL HEATR(lwind,sum,latit,par,cloud,qsn,qac,qw,xev,xco,tsup,LAKE_EVA,LAKE_PREC, jyear)
!     Sumation of heat flux (W/m2)
!   	        	
		QSNM = QSNM + qsn
		QACM = QACM + qac
		QWM  = QWM  + qw
		XEVM = XEVM + xev 
		XCOM = XCOM + xco
		daily_EVA =daily_EVA +LAKE_EVA
		daily_PREC=daily_PREC+LAKE_PREC
		
        do i=1,ns
            res_storchkH=res_storchkH+vol(i)
        enddo
! 
!   Mixing is done by mixer
! 
        
		CALL MIXER (AECO, AEWS, AKETIL, AKE)	
!  Mix out instabilities
!  IF reservoir is completely mixed (ns=1) THEN skip diffusion
!CMB   EXTRA CALL TO THICK AS IN DYRESM
      CALL THICK       		
      CALL STABLE
      IF (ns .eq. 1) GO TO 602
!  DO calculations on energy inputs, buoyancy frequency, and mixing time
!wef	IF (ntot .EQ. 1 .AND. iclock .EQ. 0.0) THEN
!wef		OPEN(61,FILE="ener.txt",FORM='FORMATTED',STATUS='UNKNOWN')
!wef		WRITE(61,89)"iclock","einff","einfw","eww","diss"
!wef89		FORMAT(A6,T10,A6,T20,A6,T30,A6,T40,A6)
!wef	END IF

      CALL ENER 

      IF (ntot .eq. ndays) THEN				! CLOSE plunge.txt for output
			CLOSE(61)
		ENDIF

602   CONTINUE
!gbs*****************ADDED STUFF FROM BILL
		IF (ntot .eq. 1) THEN					! OPEN plunge.txt for output
			OPEN(63,FILE="plunge.txt",FORM='FORMATTED',STATUS='UNKNOWN')
		ENDIF      

!wef 11dec03  prevent ns = 1 when calling inflow
		IF (ns .eq. 1) THEN
			CALL thick
		ENDIF   
   
!gbs	IF (jday .eq. 2002365) THEN
!		CALL INFLOW(XINF,N_Inf,P_Inf)		! sub-daily inflow 
!gbs         CALL INFLOW(xinf,N_Ins,P_Ins,N_Ent,P_Ent,N_Inf,P_Inf,
!gbs     &      Part_Ent,Part_Inf,Part_Ins,Balance_N_Stack,Balance_P_Stack,
!gbs     &      Balance_Part_Stack,jday,Part_Stack_F,Part_Stack_I)		! sub-daily inflow   In Dcalcs.for
!gbs	ELSE
!gbs*******************Sub-daily river temperature variation******************
 
!gbs***********************END of river temperature variation***************** 
!	print*,'b_inflow',wqual(ns,28)
!    write(*, fmt='(2i8,f12.6)')jday,iclock,wqual(ns,28)
        t3=t3+1         
		CALL INFLOW(xinf,N_Ins,P_Ins,N_Ent,P_Ent,N_Inf,P_Inf,						&
           Part_Ent,Part_Inf,Part_Ins,Balance_N_Stack,Balance_P_Stack,		&
           Balance_Part_Stack,Part_Stack_F,Part_Stack_I,t1,t2,t3,ts_flows)		! sub-daily inflow   In Dcalcs.for    
            t1=t2
		allflows=allflows+ts_flows
        do i=1,ns
            res_storchkI=res_storchkI+vol(i)
        enddo
!		print*,'e_inflow',wqual(ns,28)
!    write(*, fmt='(2i8,f12.6)')jday,iclock,wqual(ns,28)
!	print*,allflows,ts_flows, iclock,nosecs
!	pause
!         CALL INFLOW(XINF,N_Inf,P_Inf)		! sub-daily inflow   
!gbs	ENDIF
!wef 26feb04 ************ catch here when jday = 2002365 ******************************* 
		IF (ntot .eq. ndays) THEN				! CLOSE plunge.txt for output
			CLOSE(63)
		ENDIF
! DO withdrawal for each offtake
! Include water quality and particles
!wef	IF (nosecs .EQ. 0) THEN					! wef for each day only 
      DO 680 I=1,NUMOUT

!-------commented out for other cases--------------------------
!	 DRW(I) = DRW(I)*XOUT(I)
			IF(DRW(I).GE.0.0D0) THEN
				IF(OLEV(I).gt.depth(ns)) THEN       
					HH=depth(ns-1)
				ELSEIF(OLEV(I).lt.depth(ns)) THEN
					HH = OLEV(I)
				ENDIF
          ENDIF	
!            PRINT*,HH, DEPTH(NS), OLEV(I)
!----***-Only for lake Tahoe special climate change CASE-****** comment for other CASE------        
         
!			CALL OUTFLOW_GENERATION (jyear,outflow)
!			DRW(I)=outflow 
!----******------------------------------------------------------------------------	
!2015/9			offtake = drw(i)*REAL(nosecs)/REAL(86400.0d0) 
            offtake = 1.0*drw(i)*REAL(nosecs)/REAL(OUTFLO_TIMESTEP)
!            print*,drw(i), OUTFLO_TIMESTEP
! 
!     BUGG FOUND AND MODIFY BY QUIM PEREZ (09/08/98) 
!     DRW(i) (TRUE) INSTEAD OF OUTFLOW1 (BUG)
! 	 
!!!wef eliminate CALL IF there is no flow
       IF ((depth(ns)-OLEV(I)).gt.0.0) THEN	! greater than the crest or maximum legal limit	       
			IF (drw(i)>0.0) THEN
				CALL OUTFLO(HH,offtake,EXTEMP,EXSALT,EXWQ,EXIFT,EXIFS,EXIFWQ,WDL,EXCF,EXIFCF)
				alloutflows=alloutflows+offtake
			ENDIF	
       ENDIF
! 
!gbs Nutrient Budget							
! 
			N_Out = N_Out+(exwq(15)+exwq(16)+exwq(17)+exwq(18)+exwq(19)+  &
     				exwq(20))*drw(i)*REAL(nosecs)/REAL(86400.0D0) 

			P_Out = P_Out+(exwq(11)+exwq(12)+ exwq(13)+				& 
					exwq(14))*drw(i)*REAL(nosecs)/REAL(86400.0D0)
! 
!     Determine the properties of the offtake
!     Include water quality at the moment						       
! 
			EXPROP(1,I) = EXTEMP					! needs to move into sub-daily timestep
			EXPROP(2,I) = EXSALT
			EXPROP(3,I) = DENSTY(EXTEMP,EXSALT)+1000.0d0
			IF (EXTEMP.LT.1d-5 .AND. EXSALT.LT.1d-5) EXPROP(3,I)=0.0d0
			TOLAY(I) = WDL(1)
			BOLAY(I) = WDL(2)
		  DO j = 1,28
				exprop_WQ(j,i)= exwq(j)
		  ENDDO
!!!wef add up overflow and CALL outflo 
			oflow1=0.0				! calculate the overflow
!gbs Changes of Bill
!!wef 24mar05	DO k=1,numout - added CRLNGTH to the *.fiz file
!	oflow1=oflow1+(1.382d0*CRLNGTH*(depth(ns)-(CRL-BASE))
!     &		**3.577d0)*86.4d0
!	oflow1=MIN(oflow1,vol1(ns)-vcrl)
!wef 24mar05	END DO
!hanges of Bill 
			DO k=1,numout	
!gbs2010_02_01----bug----olev(k), base is already deducted in DLM1D_WQ SUBROUTINE------
!		oflow1=oflow1+(1.382d0*owid(k)*(depth(ns)-(olev(k)-BASE))
!          print*,depth(ns), olev(k)
!----------------Bill's formual-----------------
!gbs		IF ((depth(ns)-olev(k)).gt.0.0) THEN		          
!gbs             oflow1=oflow1+(1.382d0*owid(k)*(depth(ns)-olev(k))
!gbs     &		  **3.577d0)*86.4d0	 !86.4=86400/1000.0 (daily/10^3)
!gbs          ELSE
!gbs             oflow1 = 0.0
!gbs	    ENDIF
!gbs		oflow1=MIN(oflow1,vol1(ns)-vcrl)  
!-----------Broad crested weir formula-----------------------------------------------------	
				IF ((depth(ns)-CRL).gt.0.0) THEN	! greater than the crest or maximum legal limit	          
					oflow1=oflow1+(1.705d0*owid(k)*(depth(ns)-CRL)**1.5)*(REAL(nosecs))/1000.0d0                       ! convert to m3 to 10^3 m3 		  	
				ELSE
					oflow1 = 0.0
				ENDIF
!		oflow1=MIN(oflow1,vol1(ns)-vcrl)  !Bill's formula
			ENDDO
			overflow=oflow1  !/(REAL(86400.0d0)/REAL(nosecs))		!limit outflow to volume over the crest   
			IF (OFLOW1 .GT. zero) THEN
				CALL OUTFLO(HH,overflow,EXTEMP,EXSALT,EXWQ,EXIFT,EXIFS,EXIFWQ,WDL,EXCF,EXIFCF)
			ENDIF
			ts_overflow=overflow
			alloverflows=alloverflows+ts_overflow
! 
!gbs Nutrient Budget							
! 
			N_Out = N_Out+(exwq(18)+exwq(19)+exwq(20)+exwq(21)+exwq(22))*overflow 
			P_Out = P_Out+(exwq(13)+exwq(14)+ exwq(15)+exwq(16)+exwq(17))*overflow
			conc=exwq(3)
            if (k.eq.1) WRITE(UFLO, fmt='(2i8,12f12.3)')jday,iclock,float(iclock)/3600.0,extemp,exwq(8)
680   CONTINUE 
       do i=1,ns
            res_storchkO=res_storchkO+vol(i)
       enddo
    

      !**********************************************************
		CALL mass_bal(3,oflow1, conc)		
!	CALL h2o_bal(jday,oflow1)				
!**********************************************************
!wef	IF ((ntot .EQ. ndays) .AND. (nosecs .EQ. 86400)) THEN
!wef		CLOSE(83)
!wef	END IF
       do i=1,ns
            res_storchkE=res_storchkE+vol(i)            
       enddo
       !--------------Daily WATER BALANCE --------------------------------------------
      OPEN(21,file='water_balance_timestep.txt',status='unknown',access='append')
		 WRITE(21,fmt='(2i8,12f12.2)') jday,iclock, res_storchkB,-LAKE_EVA,depth(ns),res_storchkH,ts_flows,res_storchkI, depth(ns),&
       -offtake,-overflow,res_storchkO,depth(ns),res_storchkE
		CLOSE (21)
       
      
!gbs*************************END OF ADDED STUFF******************************************    
! 
!     Determine the change in chlorophyll concentration
! 

1        CALL ALGAE_TAHOE(par,limname,SDname,latit) !line 651 to 688
!		print*,'algae',wqual(ns,1),cf(ns,1)
!		CALL zooplankton
!       write(*, fmt='(2i8,f12.6)')jday,iclock,wqual(ns,28)
!       pause
! 
!     Determine the change in nutrient concentrations
! 
	
		CALL NUTRIT_TAHOE(deltnho,Al_THP,Al_NH4,Al_NO3,N_anoxic)
		Sed_Al_P   = Sed_Al_P + Al_THP
		Sed_Al_N   = Sed_Al_N + Al_NH4
      Sed_Al_NO3 = Sed_Al_NO3 + Al_NO3
		Anox_Atm_N = Anox_Atm_N + N_anoxic
! 
!     Determine change in oxygen
!    
		CALL OXYGEN(deltnho) 
!***********dissolved oxygen should not be more than it gas holding capacity***
!		DO i=1,ns
!			satDO(i) = 14.5532-0.38217*temp(i)+0.0054258*temp(i)*				&
!             temp(i)-(sal(i)/1.80655)*(1.665d-4 - 5.866d-6*temp(i)+9.796d-8*temp(i)*temp(i))
!			IF(wqual(i,8).gt.satDO(i)) THEN
!				wqual(i,8) =satDO(i)
!	      ELSE
!            wqual(i,8) =wqual(i,8)
!			ENDIF	
!        ENDDO
!******************************************************************************
! 
!     Settling of algae, phyto N & P, and Detritus N & P (new)
!     Modified from McCord routine
! 

!		CALL ALG_SETL(Lost_1, par)
		Set_Lost_N = Set_Lost_N + Lost_1(7) + Lost_1(8)  
		Set_Lost_P = Set_Lost_P + Lost_1(4) + Lost_1(5) 
!		print*,'algaesetl',wqual(ns,1),cf(ns,1)
! 
!     Settling inorganic particles and BOD (new)... 
! 
!		CALL PARTICLES(Lost_2)
!		print*,'particle',wqual(ns,1),cf(ns,1)		
!	   CALL COAGULATION (Lost_2)
!gbs09.10	ENDIF

      DO i_p = 1,7
			part_Lost(i_p) = part_lost(i_p) + Lost_2(i_p)
		ENDDO

		Part_Lost_N = Part_Lost_N + Lost_2(10) 
		Part_Lost_P = Part_Lost_P + Lost_2(9)  
619   CONTINUE

      IF (ns.eq.1) GO TO 620

! DO diffusion integrations

      CALL DIFUSE

620   CONTINUE

!  Check mixed layers for volume

      CALL THICK
		CALL STABLE
!		print*,'thick',wqual(ns,1),cf(ns,1), jday
	

!***********dissolved oxygen should be more than it gas holding capacity***
!		DO i=1,ns
!			satDO(i) = 14.5532-0.38217*temp(i)+0.0054258*temp(i)*				&
!             temp(i)-(sal(i)/1.80655)*(1.665d-4 - 5.866d-6*temp(i)+9.796d-8*temp(i)*temp(i))
!			IF(wqual(i,8).gt.satDO(i)) THEN
!				wqual(i,8) =satDO(i)
!	      ELSE
!            wqual(i,8) =wqual(i,8)
!			ENDIF
!			IF(wqual(i,8).lt.0.0d0) wqual(i,8) = 0.001d0
!		ENDDO
!***********************************************************************
      iclock = iclock+nosecs

		IF(SD.le.0.0d0) SD=0.0d0
		IF(iclock.eq.itmpr.or.itmpr.eq.86400) THEN
			OPEN(27,file= 'secchi.txt',status='unknown',access='append')
			WRITE(27,1732)jyear, jday,SD
			OPEN(203,file= SDname,status='unknown',access='append')
			WRITE(203,1731) jyear,jday,SD,A_Tahoe,B_Tahoe,Kd_Tahoe,Bsed_Tahoe,		&
     				Chloro_Tahoe,POP_Tahoe,PON_Tahoe,Nitrate_Tahoe,Ammonia_Tahoe,		&
     				THP_Tahoe,(Parts_Tahoe(k),k=1,7)
1731		FORMAT(2I6,x,f8.2,x,10(f10.2,x),7(f16.1,x)) 
			CLOSE(203)
1732		FORMAT(2I6,3x,f8.2)
			CLOSE(27)
      ENDIF
	
! 
!     WRITE on the simulation file.
!     Add all the new output files (LakeNumber, Heat, Nutrit,...)
! 
!2015/09      ITMPR=sunrise
        ITMPR=54000
!		IF (iclock .eq. ITMPR .and. ITMPR .ne. 86400) THEN
		IF (iclock .eq. 54000.and.ITMPR .ne. 86400) THEN
			CALL POSTPRO(switch,thestamp) !writing in the simulation file	simfile
			IF(PLAKE) THEN
				CALL LAKENUM (ULKN,lk)
			ENDIF
	 
			IF(PWORK) THEN
				CALL SCHMIDT (uwor,uwor2,lk,Z_min,yeari)
			ENDIF			
!     WRITE on the mixing file
! 
!			IF (PMIXER) WRITE (13,670) jday,AECO, AEWS, AKETIL, AKE
			IF(PMIXER) WRITE (13,670) jyear,AECO,AEWS,AKETIL,AKE
670		FORMAT (I6, 4(2X, F12.6))	
		ENDIF
      
! 
!     END of simulation, go to the beginning
! 
      
		IF(ntot .eq. ndays .and. iclock .eq. ITMFN.and. ITMFN .ne. 86400)GOTO 30	
! 
!     Same day,CONTINUE simulation for another time step 
! 
      IF(iclock .lt. 86400) GOTO 600	
	
!********************************************************************************************
!**********************************END SUB-DAILY TIMESTEP***********************************
!*******************************************************************************************
! 
!     END of forcing-mixing-diffusion loop
!     DO inflow, Ground Water, and Outflow on a daily basis 
! 	

      sum=0.0D0
!-----------------------------------------------
		 red_part=1.0 - 0.0
		 red_TP  =1.0 - 0.0
		 red_TN  =1.0 - 0.0	 
!		red_part=reduction
!		red_TP  =reduction
!		red_TN  =reduction
!!		CALL ATMOS_DEPOSITION(jday,jyear,reduction,N_Atm,P_Atm,Part_Atm)
!2015/09		CALL ATMOS_DEPOSITION(jyear,red_part,red_TP, red_TN,N_Atm,P_Atm,Part_Atm)! 
!2015/09      CALL shore_erosion(jyear,red_part,Part_shore)! 
		CALL GROUND_WATER(jyear,red_TP, red_TN,N_Ground,P_Ground,vol_gw)
           daily_gw=vol_gw
		IF (ntot .EQ. 1) THEN
			OPEN(73,FILE="surface.txt",STATUS='UNKNOWN',access='append')
		ENDIF

		WRITE(73,84)jday, depth(ns)
84		FORMAT(I10,F16.4)

		IF (ntot .EQ. ndays) THEN
			CLOSE(73)
		ENDIF
!--------------Daily SEDIMENT RELEASED NUTRIENTS----------------------------
		OPEN(28,file= 'sediments_NP.txt',status='unknown',access='append')
		WRITE(28,fmt='(2i6, 3f20.2)') jyear, jday,Sed_Al_P, Sed_Al_N,Sed_Al_NO3
		CLOSE (28)
!--------------Daily HEAT BUDGET--------------------------------------------
!		QSNM =  QSNM/float(iclock)
! 		QACM =  QACM/float(iclock) 
!		QWM  =  QWM /float(iclock)  
!		XEVM =  XEVM/float(iclock) 
!		XCOM =  XCOM/float(iclock)
	 
		fluxradi   = QSNM + QACM + QWM
     	fluxnoradi = XEVM + XCOM 
		sumaflux   = QSNM + QACM + QWM + XEVM + XCOM 
! 
!		WRITE the surface energy fluxes	
		IF(PHEATR) WRITE(22,4444) iseq(2),QSNM,QACM,QWM,(QACM+QWM),				&
     	XEVM,XCOM,fluxradi,fluxnoradi, sumaflux,temp(ns),temp(1)	     		
4444	FORMAT(I5,9F12.1,2F6.1)
!--------------Daily WATER BALANCE --------------------------------------------
      OPEN(29,file='water_balance.txt',status='unknown',access='append')
		WRITE(29,fmt='(1i8, 7f20.2)') jday, allflows*1000.0,				&
          daily_gw*1000.0,daily_PREC,daily_EVA,alloutflows*1000.0d0,			&
          alloverflows*1000.0d0, area(ns)*arfac 
		CLOSE (29)
!------------------------------------------------------------------------------           
!     Nutrient Budget
!---------------------------------------------------------------------------- 
		N_Out = N_Out + (exwq(21) + exwq(22)+ exwq(23)+ exwq(24) + exwq(25))*oflow1

		P_Out = P_Out + (exwq(15)+ exwq(16)+ exwq(17) + exwq(18)+ exwq(19))*oflow1
! 
!     Check thick of layers
! 
      CALL THICK
! 
!   END of single day loop, RETURN to start of daily loop
! 
! 
!     Check Nutrient Balance
! 
		DO i = 1,ns
			Nitrogen_T = Nitrogen_T + (wqual(i,21) + wqual(i,22) +	wqual(i,23) &
			+ wqual(i,24) + wqual(i,25) + wqual(i,26))*vol(i)
			Phosphorus_T = Phosphorus_T + (wqual(i,15) + wqual(i,16) + wqual(i,17) &
			+ wqual(i,18) + wqual(i,19)+ wqual(i,20))*vol(i)
      ENDDO
	      
      DO i_p = 1,7
			DO i = 1,ns
				Part_T(i_p) = Part_T(i_p) + cf(i,i_p)*vol(i)*1000.0D0
			ENDDO
      ENDDO
! 
!     Perform Nutrient Budget calculations. Units are g_day
! 
		InLake_Delta_N = Nitrogen_T   - Nitrogen_T_1
		InLake_Delta_P = Phosphorus_T - Phosphorus_T_1

      Balance_N = InLake_Delta_N + N_Out+ N_Ent + Anox_Atm_N							&  
     	      + Set_Lost_N + Part_Lost_N - (N_Atm + N_Ground + N_Ins + Sed_Al_N)
      
		Balance_P = InLake_Delta_P + P_Out + P_Ent											& 
     	      + Set_Lost_P + Part_Lost_P- (P_Atm + P_Ground + P_Ins + Sed_Al_P)
    
      DO i_p =1,7
 			InLake_Delta_Part(i_p) = Part_T(i_p) - Part_T_1(i_p)
     
			Balance_Part(i_p) = InLake_Delta_Part(i_p) + Part_Out(i_p)+					&
     		 Part_Ent(i_p) +Part_Lost(i_p) - Part_Ins(i_p) - Part_Atm(i_p)- Part_shore(i_p)
      ENDDO	 

!MB    WRITE the screen message describing the progress of the simulation

       IF(jday.gt.0)THEN
			WRITE(*,462)jday,REAL(ntot*100)/REAL(ndays)
       ELSE
			WRITE(*,464)jstart,ntot*100/ndays
       ENDIF	       
	 	
462	FORMAT(' Running day ',i8,',',2x,f7.3,'% of days complete')
464   FORMAT(' Start date  ',I8,',',2X,I5,' % complete')


      IF (ITMPR .eq. 86400) THEN
!			CALL POSTPRO(jday,nchl,switch,thestamp)
			SWITCH = .FALSE.  
! 
!     Calculate and WRITE Lake Number
! 
		IF(PLAKE) THEN
			CALL LAKENUM (ULKN,lk)
		ENDIF

		IF(PWORK) THEN
			CALL SCHMIDT (uwor,uwor2,lk,Z_min,yeari)
		ENDIF
! 
!     WRITE Outflow File 
! 
!     WRITE Mixing Energy Budget
! 
		IF(PMIXER) WRITE(13,670) jyear,jday,AECO,AEWS,AKETIL,AKE
! 
!     WRITE daily Nutrient Budget Control. Kg_day
! 

		IF(PNUTRIT) THEN
		  WRITE(unut,9091)jday,Balance_N/1000.0,Balance_N_Stack/1000.0,				&
     			InLake_Delta_N/1000.0,Nitrogen_T/1000.0,N_Out/1000.0,						&
     			N_Atm/1000.0,N_Ins/1000.0,N_Inf/1000.0,N_Ground/1000.0,					&
           (Set_Lost_N + Part_Lost_N)/1000.0,Sed_Al_N/1000.0,							&
     			Anox_Atm_N/1000.0, Balance_P/1000.0,Balance_P_Stack/1000.0,				&
           InLake_Delta_P/1000.0,Phosphorus_T/1000.0,P_Out/1000.0,P_Atm/1000.0,	&
           P_Ins/1000.0,P_Inf/1000.0,P_Ground/1000.0,										&
           (Set_Lost_P + Part_Lost_P)/1000.0,Sed_Al_P/1000.0  
! 
!     Re-inicialize the Nutrient Budget variables
! 
			Nitrogen_T_1   = Nitrogen_T 
			Phosphorus_T_1 = Phosphorus_T
	
			DO i_p =1,7
				Part_T_1(i_p) = Part_T(i_p)
			ENDDO

		ENDIF
!9091		FORMAT(i8,2X,23(F15.2,1X))

      ENDIF

		IF(PNUTRIT) THEN
			WRITE(unut,9091)jyear, jday,Balance_N/1000.0,Balance_N_Stack/1000.0,		&
     			InLake_Delta_N/1000.0,Nitrogen_T/1000.0,N_Out/1000.0,						&
				N_Atm/1000.0,N_Ins/1000.0,N_Inf/1000.0,N_Ground/1000.0,					&
				(Set_Lost_N + Part_Lost_N)/1000.0,Sed_Al_N/1000.0,							&
				Anox_Atm_N/1000.0, Balance_P/1000.0,Balance_P_Stack/1000.0,				&
				InLake_Delta_P/1000.0,Phosphorus_T/1000.0,P_Out/1000.0,P_Atm/1000.0,	&
				P_Ins/1000.0,P_Inf/1000.0,P_Ground/1000.0,									&
				(Set_Lost_P + Part_Lost_P)/1000.0,Sed_Al_P/1000.0  
! 
!     Re-inicialize the Nutrient Budget variables
! 
			Nitrogen_T_1   = Nitrogen_T 
			Phosphorus_T_1 = Phosphorus_T
	  
			DO i_p =1,7
	  			Part_T_1(i_p) = Part_T(i_p)
			ENDDO
		ENDIF
9091	FORMAT(2I6,1X,23(F15.2,1X))
! 
!     Set the clock to zero to start another day
! 
      iclock = 0
    	IF (ntot .eq. ndays .and. ITMFN .eq. 86400) THEN
! 
!     Create the LOG file as a final step detailing all parameters
!     used in the simulation that are variable
! 
!      CALL Log_File (partname, fpar,fin1)  

!MB  AT THIS STAGE REINITIALISE TO PREPARE FOR A NEW RUN
!MB  AS THE CODE DOES NOT ALLOW FOR CONTINUATION OF THE PRESENT RUN
!MB  CLOSE ALL OF THE FILES

		CLOSE (4)
		CLOSE (7)
		CLOSE (8)
		CLOSE (10)
		CLOSE (13)
		CLOSE (18)
		CLOSE (19)
		CLOSE (20)
		CLOSE (22)
		CLOSE (27)   !secchi depth file
		CLOSE (28)   !sediment release file
		RETURN	 
      ENDIF
! 
!*******************************************************************************************
!     CONTINUE daily loop
!********************************************DAILY TIME STEP LOOP ENDS**********************
! 
		GOTO 460 
      END SUBROUTINE SIMWQ			 
     
!************************************************************************
		SUBROUTINE SELFILE(FIN1,PARTNAME,INIFILEYN,WORKFILE)
!****************************************************************************
!MB     THIS GETS A .INI file and ALSO RETRIEVES THE PART OF IT THAT 
!MB     IS COMMON TO ALL THE OTHER INPUT FILENAMES
      USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
		CHARACTER*20 PARTNAME,FIN1
		CHARACTER*4 xs
		CHARACTER*20 WORKFILE
		LOGICAL*4 INIFILEYN
		LOGICAL FTEST
! 
!  quim New file reading structure
! 
		xs = '.IN1'
		PARTNAME = WORKFILE 

		WRITE(FIN1,150)PARTNAME,xs
150		FORMAT(A11,A4)
		INIFILEYN = .TRUE.
		RETURN
		END SUBROUTINE SELFILE

!************************************************************************
      SUBROUTINE FOPEN(PROMPT,LUN,FNAME,EXT,ST,FM) 
!***********************************************************************
!                                                                      
!     SUBROUTINE FOPEN                                                 
!                                                                      
!     DESCRIPTION:      This SUBROUTINE investigates the possibility of  
!                     openning a file correctly using the information  
!                     passed about the file's name, extension, status
!                     and FORMAT.                                      
!                                                                      
!     THE LOGIC IS:                                                    
!                                                                      
!         STEP 1.       IF fname = ' ' the user is prompt to enter a     
!                     filename.  The filename is converted to upper      
!                     CASE.                                            
!                                                                      
!         STEP 2.     The filename extension is checked against ext.   
!                     IF extension does not agree THEN the user is     
!                     given the option to exit the program or            
!                     CONTINUE.                                                                          
!                     IF the user chooses to CONTINUE THEN set           
!                     fname = ' ' and GOTO STEP 1.                     
!                                                                      
!         STEP 3.     The filename with correct extention is checked   
!                     IF it exists and appropriate action is taken     
!                     according to the file status.                    
!                                                                      
!     VARIABLES:                                                       
!                                                                      
!         LUN      (i*4)     file LOGICAL unit number                  
!         LENGTH   (i*4)     length of a CHARACTER string (ie. fname)  
!                                                                      
!         PROMPT   (CH**)    prompt for user to enter file name        
!         FNAME    (CH**     file name                                 
!         EXT      (CH*3)    file ext - 3 characters long              
!         ST       (CH**)    file status                               
!         FM       (CH**)    file FORMAT                               
!                                                                      
!         FEXIST   (L*1)     files exist?                              
!                                                                      
!     REQUIRED SUBPROGRAMS:                                            
!                                                                      
!         QUIT        This SUBROUTINE outputs a message as the user    
!                     exits the program.                               
!                                                                      
!     NON-STANDARD FUNCTION CALLS:                                     
!                                                                      
!         CHARNB      Lahey F77 intrinsic funtion which removes        
!                     trailing blanks from CHARACTER strings.          
!                                                                      
!***********************************************************************
!                                                                       
      USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

      INTEGER*4 LUN, IDX, i, LENGTH                                     
      CHARACTER*(*) PROMPT, FNAME, EXT, ST, FM                          
      CHARACTER*1 QUERY, OPTION                                         
      CHARACTER*11 FFORM                                                
      LOGICAL*1 FEXIST,SUCCESS
								       
      SUCCESS = .FALSE.                                                 
1000  IF (FNAME .eq. ' ') THEN            !Atenció això no és PODRÀ EXECUTAR                              
!MB  MODIFIED DUE TO CONVERSION (1 LINE) 
			WRITE (*,3010) 'CHARNB(PROMPT)'                               
3010     FORMAT (/, 1X, A, 1X)                                       
			READ (*, '(A)') FNAME                                       
      ENDIF                                                             
! 
! Convert lower CASE characters of FNAME, ST and FM, to upper CASE                  
!                                                                       
      LENGTH = LEN(FNAME)                                              
      DO 1005 i=1, LENGTH                                              
			IF (FNAME(i:i) .GE. 'a' .and. FNAME(i:i) .le. 'z') THEN     
				FNAME(i:i) = CHAR(ICHAR(FNAME(i:i)) - 32)             
			ENDIF                                                       
1005  CONTINUE
							  
      LENGTH = LEN(ST)                                              
      DO 1006, i=1, LENGTH                                              
			IF (ST(i:i) .GE. 'a' .and. ST(i:i) .le. 'z') THEN     
				 ST(i:i) = CHAR(ICHAR(ST(i:i)) - 32)             
			ENDIF                                                       
1006  CONTINUE     
						     
      LENGTH = LEN(FM)                                              
      DO 1007 i=1, LENGTH                                              
			IF (FM(i:i) .GE. 'a' .and. FM(i:i) .le. 'z') THEN     
				 FM(i:i) = CHAR(ICHAR(FM(i:i)) - 32)             
			ENDIF                                                       
1007  CONTINUE                                                          
!
! Check file extension                                                  
!     
      IF (EXT .ne. ' ') THEN                                                                  
			IDX = INDEX (FNAME, '.')                                          
      ENDIF                                                
           
! Query IF file exits and what FORMAT                                                          
!                                                                       
      INQUIRE (file=FNAME, EXIST=FEXIST, form=FFORM)                    
! 
! Convert FFORM to upper CASE IF required
! This is done because UNIX system RETURN variables as lower CASE
! where as DOS systems RETURN variables as upper CASE                 
!                                                                       
      LENGTH = LEN(FFORM)                                              
      DO 1009 i=1, LENGTH                                              
			IF (FFORM(i:i) .GE. 'a' .and. FFORM(i:i) .le. 'z') THEN     
				FFORM(i:i) = CHAR(ICHAR(FFORM(i:i)) - 32)             
			ENDIF                                                       
1009  CONTINUE                                                          
!                                                                       
! IF status 'NEW' and file already exits                                
! Give user option of over writing file                                 
! Or entering a new file name                                           
!                                                                       
1100  IF (FEXIST .and. ST .eq. 'NEW') THEN                              
			WRITE (*, 3020) FNAME                               
3020     FORMAT (/, 1X, 'File already exists : ', A,             & 
                 /, 1X, 'DO you want to over WRITE file (Y/N)? ')    
			READ (*,'(A)') QUERY                                        
								    
			IF (QUERY .eq. 'Y' .OR. QUERY .eq. 'Y') THEN                
				IF (FM .eq. FFORM) THEN                                  
					OPEN (UNIT=LUN, file=FNAME, status='UNKNOWN', form=FM)
					SUCCESS = .TRUE.                                      
				ELSE                                                     
					WRITE (*,3025) FNAME                          
3025           FORMAT (/, 1X, 'File has wrong FORMAT : ', A)        
					SUCCESS = .FALSE.                                     
				ENDIF                                                    
			ELSE                                                        
				SUCCESS = .FALSE.                                     
			ENDIF                                                       
								       
      ELSEIF (FEXIST .and. ST .eq. 'OLD') THEN                         
			IF (FM .eq. FFORM) THEN                                     
				OPEN (UNIT=LUN, file=FNAME, status=ST, form=FM)       
				SUCCESS = .TRUE.                                      
			ELSE                                                        
				WRITE (*,3025) FNAME                          
				SUCCESS = .FALSE.                                     
			ENDIF 			       
      ELSEIF (FEXIST .and. ST .eq. 'UNKNOWN') THEN                     
			IF (FM .eq. FFORM) THEN                                     
				OPEN (UNIT=LUN, file=FNAME, status='UNKNOWN', form=FM)
				SUCCESS = .TRUE.                                      
			ELSE                                                        
				WRITE (*,3025) FNAME                          
				SUCCESS = .FALSE.                                     
			ENDIF                                                       
								       
      ELSEIF (.NOT. FEXIST .and. ST .eq. 'OLD') THEN                   
			WRITE (*,3030) FNAME                                
3030     FORMAT (/, 1X, 'File does not exist : ', A)             
			SUCCESS = .FALSE.                                       
								       
      ELSEIF (.NOT. FEXIST .and. ST .eq. 'NEW') THEN                   
			OPEN (UNIT=LUN, file=FNAME, status=ST, form=FM)             
			SUCCESS = .TRUE.                                            
								       
      ELSEIF (.NOT. FEXIST .and. ST .eq. 'UNKNOWN') THEN               
			OPEN (UNIT=LUN, file=FNAME, status='UNKNOWN', form=FM)      
			SUCCESS = .TRUE.                                            
      ENDIF                                                             
								       
      IF (.NOT. SUCCESS) THEN                                           
			FNAME = ' '                                                 
1200     WRITE (*,3040)                                              
3040     FORMAT (//////, 12X, 'ENTER:',						&                      
						  ///, 25X, '1.)  CONTINUE',           &            
                      /, 25X, '2.)  QUIT PROGRAM',       &            
                    ///, 12X, 'OPTION?  ')                           
			READ (*,'(A)') OPTION                                       
								       
			IF (OPTION .lt. '1' .OR. OPTION .gt. '2') THEN
				WRITE (*,3050)                                       
3050        FORMAT (/, 5X, 'Invalid option, try again!')  
				GOTO 1200           
			ENDIF
						
			IF (OPTION .eq. '1') THEN                                   
				FNAME =' '                                            
				GOTO 1000                                             
								    
			ELSEIF (OPTION .eq. '2') THEN                              
				CALL QUIT                                             
				STOP
			ELSE                                                        
				GOTO 1200                                             
			ENDIF                                                       
      ENDIF                                                             
								       
      RETURN                                                            
      END SUBROUTINE FOPEN                                                               

!*******************************************************************************
		SUBROUTINE MOD5(fieldfile,seriefile,conturfile,blankfile,simfile,fizfile,		&
		temfile,perfile,TOUT,nvar,prof,pconta,n_profiles,conta,N_layers,neleva,FILEIN)
!*******************************************************************************
!                   
!    : Time series of linearly intepolated temperature at fixed elevations
!      e.x. (10,20,30,40)m above bottom.seriefile
!        PROFILE at selected days (Automatic or Set by the user)
!        of simulated selected variable (set by the user) 
!    	   All needed files to make contours with WINDSURF
!        based on simulated variable (T,Sal,Cla,...):
!    			BLANKING file
!             BOUNDARY file
!             CONTOUR  file
!           Written by Joaquim P. Losada and Geoff Schladow UC-Davis 2000 
!********************************************************************************
		USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER*4 nl,i,j,jtest3,neleva,ii
		INTEGER*4 print_lday,print_fday,print_jday
		INTEGER*4 nlay,N_layers
		INTEGER*4 flag                 !Printing Flag 
		INTEGER*4 conta(15) 
		INTEGER*4 cont,contador,nquim,ccont
		INTEGER*4 fday,lday	  
		INTEGER*4 nvar,mvar
		INTEGER*4 n_profiles, i_profile,i_assig
!	INTEGER*4 nnn(n_profiles+2),NNMAX     !Deactivated
		INTEGER*4 nnn(n_profiles),NNMAX
		INTEGER*4 yearfac,year,yearcounter,yearcounterj
		INTEGER*4 n_year,leap,leap_old       
		INTEGER*4 nlayer,jdayj
		INTEGER*4 yearfacj,yearj,yearfirstj,year_new
!		INTEGER*4 PARTSW,NMETSW 
      REAL*8 ddep(maxns),ttem(maxns),ssal(maxns)
      REAL*8 WWQUAl(maxns,28),ccf(maxns,7)
!IF		REAL*8 base,crl,
		REAL*8 deepest
!IF		REAL*8 ttemp,JJDAY
		REAL*8 deda(maxns),varm(maxns,28)
      REAL*8 vari(maxns),deptha(maxns) 
		REAL*8 tempera(neleva), prof(neleva)
		REAL*8 fdepth,des_dep(neleva)
		REAL*8 var(maxns)
!		REAL*8 depthy(500,n_profiles+2),varya(500,n_profiles+2)
		REAL*8 depthy(maxns,n_profiles),varya(maxns,n_profiles)     ! Deactivated Field
		REAL*8 TOUT	                   !Surface layer

!	common /DADES MODEL/ ttemp(maxns),JJDAY(maxns)

		CHARACTER*4  xs
		CHARACTER*20 conturfile
		CHARACTER*20 fizfile,temfile,perfile,blankfile,seriefile,fieldfile
      CHARACTER*20 simfile	 
		CHARACTER*20 xfile
			
		LOGICAL   FTEST
		LOGICAL*4 pconta 
!gbs
      CHARACTER*20 FILEIN	
		INTEGER*4 JJ
	 	
		INTEGER*4 JSTART,SYEAR,SYEARTEST,adddays,contday
!gbs*******************************************************************
!      initialise
      year = 0
		adddays=0

      OPEN(1,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD') 
      READ(1,*)PRODAT 
      READ(1,*)JSTART 
		CLOSE(1)
      YEAR=int(JSTART/1000)         
!gbs******************************************************************* 
!     DATA  yearfac/0/
!     DATA  yearfacj/0/

! 
!     TRIM MS-FORTRAN ROUTINE CARE MUST BE TAKEN WITH OTHER COMPILERS
!     *.SIM file names
! 
		xs = '.SIM' 
      xfile = TRIM(simfile)//xs	
6969  INQUIRE(file=xfile,EXIST=FTEST)  
      IF(.NOT.FTEST) THEN
			WRITE (*,*) 'INPUT .SIM file DOES NOT EXIST (MOD5)!'	  
			WRITE (*,*) 'NAME READ BY DLM'
			WRITE (*,*) 'PLEASE, ENTER NAME OF .SIM file: '
			READ(*,*) xfile
			GOTO 6969	  
		ENDIF
! 
!     Top and Bottom of the lake (m) above sea level 
! 
!gbs 17Jan06	OPEN(8,file=fizfile,form='FORMATTED',status='OLD')
!gbs 17Jan06	READ(8,*)base  
!gbs 17Jan06	READ(8,*)crl  
!gbs 17Jan06	deepest = crl - base
!gbs 17Jan06	CLOSE(8)		   
		IF(nlay.gt.maxns) THEN
			WRITE(*,1001)'ERROR IN MOD5: nlay',nlay,'MUST BE LESS THAN 800'
1001		FORMAT(I5,\,',')
			WRITE(*,*) 'Reduce the selected N_layers in *.WAT file'
			STOP
		ENDIF

		OPEN( 9,file=temfile,   status='UNKNOWN')	                !writing CONTOUR file	
		OPEN(10,file=xfile,     status='OLD',form='UNFORMATTED')    !.SIM file
		OPEN(11,file=perfile,   status='UNKNOWN')					!PROFILES file
		OPEN(12,file=conturfile,status='OLD',ACCESS='APPEND')		!CONTOUR file
		OPEN(13,file=seriefile, status='UNKNOWN',ACCESS='APPEND')	!TIME SERIE file
		OPEN(14,file=blankfile, status='OLD',ACCESS='APPEND')		!BLANKS file
		OPEN(15,file=fieldfile, status='OLD')						!MEASURED VARIABLES file
 		  
		i_profile    = 1 	 !Counter for Profile file
		i_assig      = 0
		contday      = 0     !Counter for first day
		depthmax     = 0     !Inicializes de max depth
      yearcounter  = 0
      yearcounterj = 0
		n_year       = 0
		leap         = 0
! 
!     READ simulated DATA on .SIM file 
! 
100   READ(10,END=5000)jday,nl,nchl,ZOOSW,PARTSW,NMETSW,				&
			(ddep(j),ttem(j),ssal(j),wwqual(j,1),wwqual(j,2),			&
			wwqual(j,3),wwqual(j,4),wwqual(j,5),wwqual(j,6),			&
			wwqual(j,7),wwqual(j,8),wwqual(j,9),wwqual(j,10),			&
			wwqual(j,11),wwqual(j,12),wwqual(j,13),wwqual(j,14),		&
			wwqual(j,15),wwqual(j,16),wwqual(j,17),wwqual(j,18),		&
			wwqual(j,19),wwqual(j,20),wwqual(j,21),wwqual(j,22),		&
			wwqual(j,23),wwqual(j,24),wwqual(j,25),wwqual(j,26),		&
			wwqual(j,27),wwqual(j,28),ccf(j,1),ccf(j,2),ccf(j,3),		&
			ccf(j,4),ccf(j,5),ccf(j,6),ccf(j,7),j=1,nl)
!gbs 17Jan06
		deepest=ddep(nl)
		nlay=int(deepest)

      IF(yearcounter.eq.0) THEN	
			yearcounter = 1
			n_year = 0
			year_new = year
			leap_old = leap	
		ENDIF
	
		IF(year.ne.year_new) THEN
			n_year = n_year + leap_old !*yearfac      		
			leap_old = leap
			year_new = year
		ENDIF
       print_jday = (jday-(jday/1000)*1000) + n_year

		cont  = 0   !Counter for PROFILE file
		ccont = 0   !Counter for CONTOUR file
! 
!     Assig days from which a profile will be extracted
! 
		IF(i_assig.eq.0) THEN
			i_assig = 1
			fday = jday-1		    	      !First day	
			lday = conta(n_profiles)         !Last day     
         print_fday = fday-(fday/1000)*1000 
		ENDIF

		IF(jday.eq.fday) THEN
			fdepth = ddep(nl)
		ENDIF
! 
!     Find the maximum depth
! 
		IF(depthmax.le.ddep(nl)) THEN
			depthmax = ddep(nl)
		ENDIF
! 
!     WRITE nodes for BOUNDARY file using Surfer	 print_jday+adddays,
!  
		WRITE(12,301) print_jday,',',ddep(nl)     !(jday-(jday/1000)*1000)+365*yearfac,',',ddep(nl)                                   !jday + 365*yearfac,',',ddep(nl)!-ddep(nl) 
301   FORMAT(I8,A1,F10.2)
! 
!     Writes nodes for BLANK file using Surfer
! 
		WRITE(14,331) print_jday,',',ddep(nl)                                   !jday + 365*yearfac,',',ddep(nl)!-ddep(nl) 
331   FORMAT(I8,A1,F10.2)
! 
!     Assignate the variable simulated to vector var. Number of variables 33
! 
		SELECT CASE(nvar)
			CASE(1)
				DO 400 i = 1,nl            !TEMPERATURE
					var(i) = ttem(i)
  400			CONTINUE
			CASE(2)					    !SALINITY
				DO 401 i = 1,nl           
					var(i) = ssal(i)
  401			CONTINUE
			CASE(3:9)                   !CHOLOROPHYLL
				DO 402 i = 1,nl           
					var(i) = wwqual(i,nvar-2)
  402			CONTINUE
			CASE(10:30)					!NUTRIENTS & Water quality
				DO i = 1,nl
					var(i) = wwqual(i,nvar-2)
				ENDDO
			CASE(31:37)				    !PARTICLES
				DO 403 i = 1,nl           
					var(i) = ccf(i,nvar-30)
  403			CONTINUE
			CASE DEFAULT
		END SELECT
! 
!     WRITE VARIABLE simulated on CONTOUR file
! 
		DO 666 j = 1,nl
			IF(ddep(j).lt.0.0)GOTO 666 !+++++++++++++++++++++++++++++++++++++++++++++
!			IF((nvar.ne.1).and.(ddep(j).lt.400))GOTO 666 
!          WRITE(9,101) print_jday+adddays,ddep(j),var(j)  ! print_jday will separate year  jday + 365*yearfac,ddep(j),var(j) 
101			FORMAT(I8,F10.2,F20.4)
666   CONTINUE
! 
!     Linar interpolation for Time Serie file 
!  
		DO 1500 i=1,nlay
			deptha(i)=i*(int(deepest)/(nlay))
			DO 1400 j=1,nl !nl=number of simulated layers
				IF(j.eq.1.and.deptha(i).lt.ddep(j))THEN 
					vari(i)=var(j)	      
					GOTO 1500
				ENDIF
				IF(j.eq.nl.and.deptha(i).gt.ddep(j))THEN	     
					IF(ccont.eq.0) THEN
						nquim = i !Nquim takes the value of the last layer (surface) 	
						deptha(nquim) = ddep(j) 
						vari(nquim) = var(j)		   
						ccont = 1
						GOTO 1500
					ENDIF
					vari(i)=0.0
					GOTO 1500
				ENDIF
  				IF(deptha(i).eq.ddep(j))THEN
					vari(i)=var(j)	      
					GOTO 1500
				ENDIF
				IF(deptha(i).ge.ddep(j).and.deptha(i).lt.ddep(j+1))THEN            
					vari(i)=var(j)+((var(j+1)-var(j))*((deptha(i)-ddep(j))/(ddep(j+1)-ddep(j))))	
					GOTO 1500
            ENDIF
1400     CONTINUE		 
1500   CONTINUE
!gbs
      
      DO i=1, neleva
			des_dep(i)=ddep(nl)-prof(i)
!	print*,i,des_dep(i),ddep(nl),prof(i)
!	pause
			IF(des_dep(i).eq.ddep(nl)) THEN
				tempera(i)=var(nl)
				flag = 1
			ELSE
				DO j=1,nl	       
					IF(des_dep(i).ge.ddep(j).and.des_dep(i).lt.ddep(j+1))THEN            
						tempera(i)=var(j)+((var(j+1)-var(j))*((des_dep(i)-ddep(j))/(ddep(j+1)-ddep(j))))
						flag = 1	    
					ENDIF
				ENDDO
			ENDIF
      ENDDO
	
!gbs 17Jan06 changed for writing from top to bottom	see bill modification 17Jan06 folder
!      DO 1200 j=1,nlay
!	 DO 3333 iii = 1,neleva            !flag control when to WRITE on
!	  IF(deptha(j).EQ.prof(iii)) THEN  !the file SERIET
!	   flag = 1
!	   tempera(iii) = vari(j)	
!wef 15sep05 **********************************************************
!	  ELSEIF((j .EQ. 1) .AND. (prof(iii) .LT. deptha(j))) THEN
!		tempera(iii) = vari(j)
!	  ELSEIF((prof(iii) .LT. deptha(j)) .AND. 
!     &				(prof(iii) .GT. deptha(j-1))) THEN
!	    tempera(iii) = var(j-1)+((var(j)-var(j-1))*
!     &               ((prof(iii)-deptha(j-1))/(deptha(j)-deptha(j-1))))
!	  ELSEIF((prof(iii) .GT. deptha(j)) .AND. 
!     &				(j .EQ. nlay)) THEN
!	    tempera(iii) = 0.d0
!	  ENDIF
!wef 15sep05 **********************************************************
!3333   CONTINUE
!1200  CONTINUE
! 		
!     WRITE linearly interpolated VARIABLE time serie at the assigned elevations
!     (For ex: 10,20,30,40m)
! 
		IF (flag.eq.1) THEN						  
         WRITE(13,404)print_jday+adddays,(tempera(jj),jj=1,neleva) !jday + 365*yearfac,(tempera(jj),jj=1,8)
404      FORMAT(I5,<neleva>F20.3)      
		ENDIF

      flag = 0
!gbs*************
		IF (MOD(FLOAT(YEAR),4.0).EQ.0.0) THEN
			IF(jday.eq.366) THEN
				year = year+1
				adddays=adddays+366
			ENDIF
		ELSE
			IF(jday.eq.365) THEN
				year = year+1
				adddays=adddays+365
			ENDIF
		ENDIF
!gbs********************

      GOTO 1447
		IF(jday.eq.lday.or.jday-1.eq.fday) THEN
3100	   READ(15,*,END=5144) nlayer
         READ(15,*)
			IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1) THEN       
				IF(nvar.gt.19) THEN
					WRITE(*,*)'ERROR MOD5: VARIABLE NUMBER must be less that 19'
					STOP
				ENDIF
				mvar = nvar
				DO 910 ii =1,nlayer
					READ(15,*) jdayj, deda(ii),(varm(ii,j),j=1,19)
 910			CONTINUE
			ELSE	 
				IF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1) THEN
 					IF(nvar.ge.20.and.nvar.le.26) THEN
                  WRITE(*,*)'ERROR MOD5: zoo and metals not simulated'
						STOP
					ENDIF
					IF(nvar.ge.27) THEN
						mvar = nvar-7
					ELSE
						mvar = nvar
					ENDIF		
					DO 920 ii =1,nlayer
						READ(15,*) jdayj, deda(ii),(varm(ii,j),j=1,26)	
 920           CONTINUE
				ELSE
               IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.eq.1) THEN	
						SELECT CASE(nvar)
							CASE(25:35) 
								WRITE(*,*)'ERROR MOD5: VARIABLE NUMBER must be less that 24'
								STOP
							CASE(20:21)
								WRITE(*,*)'ERROR MOD5: zoo not simulated, nvar 20,21'
								STOP
							CASE DEFAULT
						ENDSELECT
						mvar = nvar	 
		 
						DO 922 ii =1,nlayer
							READ(15,*) jdayj, deda(ii),(varm(ii,j),j=1,24)	   
 922						CONTINUE
					ELSE
						WRITE(*,*)'  WARNING: AT THE MOMENT MOD6 ONLY ACCEPTS'
						WRITE(*,*)'PARTICLES ON/OFF METALS ON/OFF and FITO ON/OFF'
						STOP
					ENDIF
				ENDIF
			ENDIF
! 
!     Get Julian day (1993109) ===> 109
!     Accounts for Leap year 
! 
         yearj = int(jdayj/1000)
			jtest3=MOD(yearj,4)

 			IF(jtest3.eq.0) THEN
            leap = 366
         ELSE
            leap = 365
         ENDIF
         IF(yearcounterj.eq.0) THEN	
				yearcounterj = 1
				n_year = 0
				year_new = year
				leap_old = leap	
			ENDIF
			IF(year.ne.year_new) THEN
				n_year = n_year + leap_old !*yearfac      		
				leap_old = leap
				year_new = year
			ENDIF
      
			IF(jdayj.eq.lday)GOTO 4100 
			IF(jdayj.eq.fday)GOTO 4100		
         GOTO 3100
4100     CONTINUE
! 
!     Linear interpolation to have the measured DATA at the same depths
!     that the simulated DATA (Patterson et al., 1984)
! 
!			CALL LINEAL(nl,nlayer,ddep,deda ,varm,vari)! 
!			WRITE measured variable of the last day on PROFILE file
! 
!!!      print_jdayj = (jdayj-(jdayj/1000)*1000) + n_year
!!!      WRITE(11,1210) int(print_jdayj),deda(nlayer) !+365*,deda(nlayer)
!!!		WRITE(11,*)nlayer
!!!		DO 1209 j=1,nlayer
!!!			WRITE(11,1220)deda(j),varm(j,mvar)
1209     CONTINUE 
		ENDIF

1447  CONTINUE
		GOTO 100

5000  CONTINUE
! 
!     NODES for completing the curve on BLANK file
!        
		print_lday = (lday - (lday/1000)*1000) + n_year
		WRITE(14,337)  print_lday,',',depthmax     !lday,',',depthmax             !-ddep(nl) 
		WRITE(14,337)1+print_fday,',',depthmax    !fday+1,',',depthmax           !-ddep(nl) 
		WRITE(14,337)1+print_fday,',',fdepth      !fday+1,',',fdepth             !-ddep(nl)   
337   FORMAT(I8,A1,F10.2)
	
		tout = var(nl) !Value of Variable at surface

!1331	FORMAT(<n_profiles+2>(I10,5X),/) ! Field file profiles desactivated
1331	FORMAT(<n_profiles>(I10,5X),/)
!1441 FORMAT(<n_profiles+2>(F6.2,2X,F16.3,2X,\),/) ! Field file profiles desactivated
1441  FORMAT(<n_profiles>(F6.2,2X,F16.3,2X,\),/)
	
		CLOSE(9)
      CLOSE(10)
		CLOSE(11)      
		CLOSE(12)
		CLOSE(13)
		CLOSE(14)
		CLOSE(15)
      RETURN
5144  CONTINUE
		WRITE(*,*)'WARNING. MOD5 END of File of fieldfile' 
		WRITE(*,*) 'Check IF First Day ',fday,' is in File ',fieldfile
		WRITE(*,*) 'Check IF Last  Day ', lday,' is in File ',fieldfile
		WRITE(*,*)'Please check # of simulated days'
		STOP
      END SUBROUTINE MOD5

!************************************************************************
		SUBROUTINE MOD6(nvar,simfile,fieldfile,fitfile,ERROR,ERROR_E,ERROR_H,ERROR_C)
!******************************************************************
!*  
!*               Written by JOAQUIM PEREZ (07/01/99)
!*
!*Calculates fiting function of the simulated days  
!*
!*INPUT:  nvar number of the variable to be performed 
!*        IDAY day(93xxx) of the measured DATA. Simulation DATA should contains
!*             this day
!*        simfile: binary file whith all days profile information
!*        fieldfile: file contains.PRO information to meke linear interpolation
!*                   when constructing OPTIMUS function
!*        var: vector whith the simulated DATA of the variable
!*   		vari vector whith the linearly interpolated DATA
!********************************************************************************
		USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER*4 nvar,mvar,nl,n_nk
 !     INTEGER*4 partsw,zoosw,nmetsw
		CHARACTER*20 simfile
		CHARACTER*20 fieldfile,fitfile
		CHARACTER*4 xs	
		CHARACTER*20 xfile
		REAL*8 varm(maxns),var(maxns),vari(maxns),VARFIELD(maxns)
		REAL*8 wwqual(maxns,28),ccf(maxns,7),ttem(maxns),ssal(maxns),ddep(maxns)
		REAL*8 deptha(maxns)
		REAL*8 cons,CONEXP,CONDIF,MEAN_FIELD,MEAN_SIMUL
		REAL*8 ERROR,ERROR_C,ERROR_E,ERROR_H,ERROR_VOL
		REAL*8 Scaling_Factor, error_total
		REAL*8 SUMA,SUMA_E,SUMA_H
		INTEGER*4 NLAY,I,II,J

		LOGICAL FTEST

		DATA Scaling_Factor/2.84291/
		DATA error_total/0.0/
! 
!     Build name file *.SIM and check IF it exists in the current directory
! 
		xs = '.SIM' 
      xfile = TRIM(simfile)//xs	
6969  INQUIRE(file=xfile,EXIST=FTEST)
      IF(.NOT.FTEST) THEN
			WRITE (*,*) 'INPUT .SIM file DOES NOT EXIST (MOD6)!'	  
			WRITE (*,*) 'NAME READ BY DLM', xfile 
			WRITE (*,*) 'PLEASE, ENTER NAME OF .SIM file: '
			READ(*,*) xfile
			GOTO 6969	  
		ENDIF
! 
!     OPEN Fitting and Simulation files
! 
      OPEN(2,file=fitfile,status='OLD',ACCESS='APPEND') 
		OPEN(10,file=xfile,form='UNFORMATTED',status='OLD')
		jday = 0	  
!		Searching for Initial Day (IDAY)
	
!		IF(jday.eq.IDAY) GOTO 200
!		GOTO 100 
!     READ on the FIELD file which contains the measured DATA
! 
 200	OPEN(11,file=fieldfile,status='OLD')
 300	READ(11,*,END=5000) nlay
      READ(11,*)
		IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1) THEN
			IF(nvar.gt.19) THEN
				WRITE(*,*)'ERROR MOD6: VARIABLE NUMBER must be less that 19'
				STOP
			ENDIF
			mvar = nvar
			DO 10 ii =1,nlay
				READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,19)
				varm(ii) = VARFIELD(mvar)
 10		CONTINUE
		ELSE	 
			IF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1) THEN
				SELECT CASE(nvar)
					CASE(20:21)
						WRITE(*,*)'ERROR MOD6: zoo not simulated, nvar 20,21'
						STOP
					CASE DEFAULT
				ENDSELECT
				IF(nvar.ge.27) THEN
					mvar = nvar-7
				ELSE
					mvar = nvar
				ENDIF
				DO 20 ii =1,nlay
					READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,26)
					varm(ii) = VARFIELD(mvar)
 20			CONTINUE
			ELSE
				IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.eq.1) THEN
					SELECT CASE(nvar)
						CASE(25:35) 
							WRITE(*,*)'ERROR MOD6: VARIABLE NUMBER must be less that 24'
							STOP
						CASE(20:21)
							WRITE(*,*)'ERROR MOD6: zoo not simulated, nvar 20,21'
							STOP
						CASE DEFAULT
					ENDSELECT
					mvar = nvar	 
					DO 22 ii =1,nlay
						READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,24)	             
                  varm(ii) = VARFIELD(mvar)
 22				CONTINUE
				ELSE
					WRITE(*,*)'WARNING: AT THE MOMMENT MOD6 ONLY ACCEPTS'
					WRITE(*,*)'PARTICLES ON/OFF METALS ON/OFF and FITO ON/OFF'
					STOP
				ENDIF
			ENDIF
		ENDIF      

      IF(jjday.le.jday) GOTO 207
! 
!     READ on the simulation file
! 
      n_nk = 0
100   READ(10,END=5000)jday,nl,nchl,ZOOSW,PARTSW,NMETSW,			&
			(ddep(j),ttem(j),ssal(j),wwqual(j,1),wwqual(j,2),		&
			wwqual(j,3),wwqual(j,4),wwqual(j,5),wwqual(j,6),		&
			wwqual(j,7),wwqual(j,8),wwqual(j,9),wwqual(j,10),		&
			wwqual(j,11),wwqual(j,12),wwqual(j,13),wwqual(j,14),	&
			wwqual(j,15),wwqual(j,16),wwqual(j,17),wwqual(j,18),	&
			wwqual(j,19),wwqual(j,20),wwqual(j,21),wwqual(j,22),	&
			wwqual(j,23),wwqual(j,24),wwqual(j,25),wwqual(j,26),	&
			wwqual(j,27),wwqual(j,28),ccf(j,1),ccf(j,2),ccf(j,3),	&
			ccf(j,4),ccf(j,5),ccf(j,6),ccf(j,7),j=1,nl)	
207   CONTINUE
      IF(jjday.eq.jday) GOTO 97	
		IF(jjday.gt.jday) THEN 
			n_nk = n_nk + 1
			GOTO 100
		ENDIF 	
		GOTO 300  
97    CONTINUE
!  
!     Check depths.Warning IF depth simulated < depth measured
!     IF so THEN substracts one layer to measured DATA and compares again
! 
!98	IF(ddep(nl).lt.deptha(nlay)) THEN
!			nlay = nlay-1
!			GOTO 98
!		ENDIF
! 
!     Checks depths.Warning IF depth measured < depth simulated
!     IF so THEN substracts one layer to simulated DATA and compares again
! 
98		IF(deptha(nlay).lt.ddep(nl)) THEN
			nl = nl-1
			GOTO 98
		ENDIF
! 
!     Fill vector var with the selected simulated variable
! 
		SELECT CASE(nvar)
			CASE(1)
				DO 400 i = 1,nl           !TEMPERATURE
					var(i) = ttem(i)
  400       CONTINUE
			CASE(2)					   !SALINITY
				DO 401 i = 1,nl           
					var(i) = ssal(i)
  401       CONTINUE
			CASE(3:9)                   !CHOLOROPHYLL
				DO 402 i = 1,nl           
					var(i) = wwqual(i,nvar-2)
  402       CONTINUE
         CASE(10:30)
				DO i = 1,nl
					var(i) = wwqual(i,nvar-2)
            ENDDO
			CASE(31:37)				   !PARTICLES
				DO 403 i = 1,nl           
					var(i) = ccf(i,nvar-30)
  403       CONTINUE
         CASE DEFAULT
		END SELECT
! 
!     Linear interpolation to have the simulated DATA at the same depths
!     that the field DATA (Jorgensen, 1983;Omlin,2000)
! 
!		CALL LINEAL(nlay,nl,deptha,ddep,var,vari)
! 
!     Calcule the error function IF simulated variable is interpolated
!     to measured DATA positions
! 
!		CALL OPTIMA(nlay,deptha,vari,varm,ERROR,ERROR_E,ERROR_H,ERROR_C,		&
!                 SUMA,SUMA_E,SUMA_H,JJDAY)
! 
!     Linear interpolation to have the measured DATA at the same depths
!     that the simulated DATA (Patterson et al., 1984)
! 
		CALL LINEAL(nl,nlay,ddep,deptha,varm,vari)
! 
!     Calcule the error function IF measured variable is interpolated
!     to simulated DATA positions
! 
		CALL OPTIMA(nl,ddep,var,vari,ERROR,ERROR_E,ERROR_H,ERROR_C,SUMA,SUMA_E,SUMA_H)! 
!     Add error to calculate the global error
! 
		error_total = error_total + error 
! 
!     WRITE fitting function
! 
		WRITE(2,33)JJDAY,ERROR,ERROR_E,ERROR_H,ERROR_C
 33	FORMAT(I8,3X,2F8.3,2X,F8.3,2X,F8.3)

      GOTO 300

		CLOSE(11)
      CLOSE(10)
		CLOSE(2)
		RETURN
5000  CONTINUE
! 
!     WRITE the global error
! 
		WRITE(2,43)'Global  ',error_total,ERROR_E,ERROR_H,ERROR_C
43		FORMAT(A8,2X,2F8.3,2X,F8.3,2X,F8.3)
!		WRITE(*,*)'ERROR.SIM LAST DAY FOUND ON MOD6 WITHOUT MATCHING IDAY'
!		WRITE(*,*)' jday=  ',jday,' IDAY=  ',IDAY
		CLOSE(10)
		CLOSE(11)
		CLOSE(2)
!		STOP	
		RETURN
		END SUBROUTINE MOD6
!****************************************************************
      SUBROUTINE LINEAL(nlay,nl,deptha,ddep,ttem,TEMPA)
!*****************************************************************
!*	Written by Joaquim Perez (08/01/99)
!*	Performs the linear interpolation of the input variable ttem
!*   	(at height ddep) onto heights deptha. 
!*    It colud be used to interpolate simulated
!*    variable at field depths, or to interpolate field DATA to 
!*    simulated depths
!*    
!*    ddep: vector with depths of the DATA ttem
!*    TEMPA: output vector with linearly interpolated simulated DATA
!*           onto heights deptha	
!*	ttem: input DATA vector 
!*****************************************************************
		USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER*4 I,J,nlay,nl
		REAL*8 deptha(maxns),ddep(maxns),ttem(maxns),tempa(maxns)
		DO 1500 i=1,nlay
			DO 1400 j=1,nl !nl=number of layers	                   
				IF(j.eq.1.and.deptha(i).lt.ddep(j))THEN 
					tempa(i)=ttem(j)	      
					GOTO 1500
				ENDIF
				IF(j.eq.nl.and.deptha(i).gt.ddep(j))THEN	     
					tempa(i) = -99
					GOTO 1500
				ENDIF

  				IF(deptha(i).eq.ddep(j))THEN
					tempa(i)=ttem(j)	      
					GOTO 1500
				ENDIF

				IF(deptha(i).ge.ddep(j).and.deptha(i).lt.ddep(j+1)) THEN
					tempa(i)=ttem(j)+((ttem(j+1)-ttem(j))*((deptha(i)-ddep(j))/(ddep(j+1)-ddep(j))))
					GOTO 1500
				ENDIF
1400		CONTINUE		 
1500  CONTINUE
		RETURN
		END SUBROUTINE LINEAL

!************************************************************************************************
      SUBROUTINE MOD7(simfile,numero,vector, FILEIN)
!************************************************************************************************* 
!     Ouput for MATLAB files
! 
		USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

      INTEGER*4 sio
		PARAMETER(sio = 3000)
		INTEGER*4 i,j,k
		INTEGER*4 nvar,IDAY,nl
		INTEGER*4 numero
		INTEGER*4 vector(37)
!IF		INTEGER*4 partsw,zoosw,nmetsw
      INTEGER*4 counter,prim,seco   
      INTEGER*4 secu,primu,tercu,alfa
		CHARACTER*20 simfile
		CHARACTER*20 xfile
		CHARACTER*20 fieldfile
		CHARACTER*4 xs	
      CHARACTER*20 diafile
		CHARACTER*4 tercoc,secoc,primc		
      REAL*8 wwqual(sio,28),ccf(sio,7),ttem(sio),ssal(sio),ddep(sio)
		REAL*8 deptha(sio)
		REAL*8 var(numero,sio)	
		LOGICAL FTEST
!gbs
      CHARACTER*20 FILEIN	
		INTEGER*4 year1, year2, year3, year4, year12, year13
		CHARACTER*4 yearc1, yearc2, yearc3, yearc4 	
		INTEGER*4 JSTART,YEAR
!gbs*******************************************************************
      OPEN(1,FILE=FILEIN,FORM='FORMATTED',STATUS='OLD') 
      READ(1,*)PRODAT 
      READ(1,*)JSTART 
      YEAR=JSTART/1000

		year1=(YEAR/1000)
		       
		year12=YEAR-(YEAR/1000)*1000
		year2=(year12/100)

		year13=YEAR-(YEAR/100)*100
		year3=(year13/10)
		 		
		year4=YEAR-(YEAR/10)*10  
			  
		yearc1   = char(year1+48)
		yearc2   = char(year2+48)
		yearc3   = char(year3+48)
		yearc4   = char(year4+48)
    
!gbs*******************************************************************	
		
		xs = '.SIM' 
      xfile = TRIM(simfile)//xs	
6969  INQUIRE(file=xfile,EXIST=FTEST) 
      IF(.NOT.FTEST) THEN
			WRITE (*,*) 'INPUT .SIM file DOES NOT EXIST (MOD7)!'	  
			WRITE (*,*) 'NAME READ BY DLM', xfile 
			WRITE (*,*) 'PLEASE, ENTER NAME OF .SIM file: '
			READ(*,*) xfile
			GOTO 6969	  
		ENDIF
! 
!     OPEN and READ on the *.SIM file
! 
		OPEN(10,file=xfile,form='UNFORMATTED',status='OLD')
	 		 
      i = 0
100		READ(10,END=5000)jday,nl,nchl,ZOOSW,PARTSW,NMETSW,			&
			(ddep(j),ttem(j),ssal(j),wwqual(j,1),wwqual(j,2),		&
			wwqual(j,3),wwqual(j,4),wwqual(j,5),wwqual(j,6),		&
			wwqual(j,7),wwqual(j,8),wwqual(j,9),wwqual(j,10),		&
			wwqual(j,11),wwqual(j,12),wwqual(j,13),wwqual(j,14),	&
			wwqual(j,15),wwqual(j,16),wwqual(j,17),wwqual(j,18),	&
			wwqual(j,19),wwqual(j,20),wwqual(j,21),wwqual(j,22),	&
			wwqual(j,23),wwqual(j,24),wwqual(j,25),wwqual(j,26),	&
			wwqual(j,27),wwqual(j,28),ccf(j,1),ccf(j,2),ccf(j,3),	&
			ccf(j,4),ccf(j,5),ccf(j,6),ccf(j,7),j=1,nl)
! 
!     Build the Matlab output file name using julian date method
! 
			jday    = jday-(jday/1000)*1000
			primu   = int(jday/100)
			alfa    = jday - int(jday/100)*100
			secu    = int(alfa/10)
			tercu   = alfa - int(alfa/10)*10
			primc   = char(primu+48)
			secoc   = char(secu+48)
			tercoc  = char(tercu+48)
			diafile = trim(adjustl('Y&day'))// trim(adjustl(yearc1))		&
     				 // trim(adjustl( yearc2))// trim(adjustl(yearc3))		&
     				 // trim(adjustl( yearc4))// trim(adjustl(primc))		&
					 // trim(adjustl( secoc ))// trim(adjustl(tercoc)) //'.dat'
               
	
			OPEN(11,file= diafile,status='unknown')
!gbs
      
			IF(mod(float(year),4.0).eq.0.0) THEN
				IF(jday.eq.366) THEN
					YEAR = YEAR +1
					year1=(YEAR/1000)
		       
					year12=YEAR-(YEAR/1000)*1000
					year2=(year12/100)

					year13=YEAR-(YEAR/100)*100
					year3=(year13/10)
		 			
					year4=YEAR-(YEAR/10)*10  
					  
					yearc1   = char(year1+48)
					yearc2   = char(year2+48)
					yearc3   = char(year3+48)
					yearc4   = char(year4+48)
				ENDIF
			ELSE
				IF(jday.eq.365) THEN
					YEAR = YEAR +1
					year1=(YEAR/1000)
						 
					year12=YEAR-(YEAR/1000)*1000
					year2=(year12/100)

					year13=YEAR-(YEAR/100)*100
					year3=(year13/10)
		 			
					year4=YEAR-(YEAR/10)*10  
					  
					yearc1   = char(year1+48)
					yearc2   = char(year2+48)
					yearc3   = char(year3+48)
					yearc4   = char(year4+48)
				ENDIF
			ENDIF
!gbs
! 
!     Assignates the variable simulated to vector var
! 
			DO k = 1,numero
				nvar=vector(k) 
				SELECT CASE(nvar)
					CASE(1)
						DO 400 i = 1,nl            !TEMPERATURE
							var(k,i) = ttem(i)
  400					CONTINUE
					CASE(2)					    !SALINITY
						DO 401 i = 1,nl           
							var(k,i) = ssal(i)
  401					CONTINUE
					CASE(3:9)                   !CHOLOROPHYLL
						DO 402 i = 1,nl           
							var(k,i) = wwqual(i,nvar-2)
  402					CONTINUE
					CASE(10:30)
						DO i = 1,nl
							var(k,i) = wwqual(i,nvar-2)
						ENDDO
					CASE(31:37)				     !PARTICLES
						DO 403 i = 1,nl           
							var(k,i) = ccf(i,nvar-30)
  403					CONTINUE
					CASE DEFAULT
				END SELECT
			ENDDO

			DO j = 1,nl      	
				WRITE(11,33)jday,ddep(j),(var(k,j),k=1,numero)
33				FORMAT(i8,f8.2,<numero>f18.4)
			ENDDO
			CLOSE(11)
		i=i+1
      GOTO 100
5000  CONTINUE 
      CLOSE(10)
		CLOSE(1)
		RETURN
		END SUBROUTINE MOD7
!**********************************************************************
		SUBROUTINE MOD8(nvar,fieldfile,prof,neleva)
!**********************************************************************
		USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER*4 I,II,J,nvar,mvar,nlay,neleva		
		CHARACTER*20 fieldfile
		REAL*8 varm(maxns),var(maxns),vari(maxns),VARFIELD(maxns)
		REAL*8 deptha(maxns),pruf(1),VAR_INT(4)
		REAL*8 prof(neleva),des_dep(neleva)
		CHARACTER*4 xs	
		CHARACTER*20 xfile
!IF      INTEGER*4 partsw,zoosw,nmetsw
		INTEGER*4 YEARCOUNTER,year,yearfac,dd(1500),kk,day, addday
		LOGICAL FTEST
! 
!     READ on the FIELD file which contains the measured DATA
! 
200	OPEN(11,file=fieldfile,status='OLD')
      REWIND(11)
!		OPEN(10,file='NO3_MEASURED_01_02.txt',status='UNKNOWN',
!		OPEN(10,file='Particles6_M_0102.txt',status='UNKNOWN',
		OPEN(10,file='Measured.txt',status='UNKNOWN',ACCESS='APPEND')

      YEARCOUNTER = 0 
		kk=2
		dd(1)=1
		addday=0
 300	READ(11,*,END=5000) nlay
      READ(11,*)
		IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1) THEN
			IF(nvar.gt.19) THEN
				WRITE(*,*)'ERROR MOD8: VARIABLE NUMBER must be less that 19'
				STOP
			ENDIF
			mvar = nvar
			DO 10 ii =1,nlay
				READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,19)
				varm(ii) = VARFIELD(mvar)
 10      CONTINUE
		ELSE	 
			IF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1) THEN
				SELECT CASE(nvar)
					CASE(20:21)
						WRITE(*,*)'ERROR MOD8: zoo not simulated, nvar 20,21' 
						STOP
					CASE DEFAULT
				ENDSELECT
				IF(nvar.ge.27) THEN
					mvar = nvar-7
				ELSE
					mvar = nvar
				ENDIF	
				DO 20 ii =1,nlay
					READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,26)	         
					varm(ii) = VARFIELD(mvar) 
 20         CONTINUE
            
			ELSE
				IF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.eq.1) THEN
					SELECT CASE(nvar)
						CASE(25:35) 
							WRITE(*,*)'ERROR MOD8: VARIABLE NUMBER must be less that 24'
							STOP
						CASE(20:21)
!							WRITE(*,*)'ERROR MOD8: zoo not simulated, nvar 20,21'
!							STOP
						CASE DEFAULT
					ENDSELECT
					mvar = nvar
					DO 22 ii =1,nlay
						READ(11,*) JJDAY, deptha(ii),(VARFIELD(j),j=1,24)	           
						varm(ii) = VARFIELD(mvar)
 22				CONTINUE
				ELSE
					WRITE(*,*)'WARNING: AT THE MOMENT MOD8 ONLY ACCEPTS'
					WRITE(*,*)'PARTICLES ON/OFF METALS ON/OFF and FITO ON/OFF'
					STOP
				ENDIF
			ENDIF
		ENDIF	
!gbs	
		DO i=1,neleva
			des_dep(i)=deptha(nlay)-prof(i)
!			print*,des_dep(i)
		ENDDO
!gbs
!     Linear interpolation to have the measured DATA at the set elevations
		CALL LINEAL(neleva,nlay,des_dep,deptha,varm,vari)
!gbs--------------------------------------------       
      year=int(jjday/1000)
		dd(kk)=jjday-1000*int(jjday/1000)  
		IF(dd(kk).lt.dd(kk-1)) THEN
			IF(mod(year,4).eq.0) THEN
				addday=addday+366
			ELSE
				addday=addday+365
			ENDIF	  
		ENDIF
		day=dd(kk)+addday	
		kk=kk+1
!gbs---------------------------------

		WRITE(10,111) DAY,(vari(i),i=1,neleva)   !JJDAY + 365*yearfac
111   FORMAT(I9,<neleva>f20.4)  !16 is based on current situation ELSE it is neleva
      GOTO 300

5000  CONTINUE
		CLOSE(11)
		CLOSE(10)
		RETURN
		END SUBROUTINE MOD8
!********************************************************

