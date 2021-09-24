 	   PROGRAM DLM1DWQ
!*******************************************************************************
	  USE MSFLIB          !SPECIFIC LIBRARY OF MS-FORTRAN FOR NAME
      USE PORTLIB         !SPECIFIC LIBRARY OF MS-FORTRAN FOR TIME 

      USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
!
!    NEW DEFINITIONS
!
	   INTEGER*4 ndia,lastday,jtest3
	   INTEGER*4 FITXER
	   INTEGER*4 COUNT1
!IF	   INTEGER*4 jjday      
	   INTEGER*4 prim, seco, dec,uni
	   INTEGER*4 nvar,iday
	   INTEGER*4 numero,vector(37)
	   INTEGER*4 neleva
	   INTEGER*4 nconta,n_profiles
	   INTEGER*4 n_layers
	   INTEGER*4 conta(15)
	   INTEGER*4 i,j,k,jk,jj,kk,lf,i_nexp,nexp

	   LOGICAL*4 NZOOSW
	   LOGICAL*4 NPARTSW
	   LOGICAL*4 METSW
	   LOGICAL*4 MATLAB
	   LOGICAL   FTEST
	   LOGICAL*4 PCONTA

!IF	   REAL(8) ::  TTEMP
	   REAL(8) ::  TOUT
	   REAL(8) ::  prof(500)
	   REAL(8) ::  ERROR,ERROR_VOL,ERROR_C,ERROR_E,ERROR_H
	   REAL(8) ::  MEAN_FIELD,MEAN_SIMUL
!	   common /DADES MODEL/ TTEMP(50),JJDAY(50)

       PARAMETER(FITXER = 15)
	   CHARACTER*40 file(10)
	   DATA file /'*.PAR','*.IN1','*.LOG','*.RUN','*.CLI', '*.QIN', '*.QOT','*.WQD','*.SIM','*.SMI'/
	   CHARACTER*1   TKOPT                                                
       CHARACTER*20  workfile
	   CHARACTER*12  fileout
	   CHARACTER*12  fileout2, fileout3, fileout4
	   CHARACTER*20  fileout5, fileout6, fileout7
	   CHARACTER*20  fileout8, fileout9, fileout10,fileout11
	   CHARACTER*20  fileout12,fileout13,fileout14,fileout15
       CHARACTER*20  fileout16,fileout17,fileout18,fileout19
	   CHARACTER*20  fileout20,fileout21,fileout22,fileout23
	   CHARACTER*20  fileout24, fileout25, fileout26 
	   CHARACTER*12  fileblank
	   CHARACTER*12  filecont
	   CHARACTER*12  fileserie
	   CHARACTER*30  fizfile,temfile,perfile,blank,contur
	   CHARACTER*30  fieldfile,HeatFile,HeatDay,seriefile,fitfile
	   CHARACTER*30  LakeFile,NutFile,MixFile,InFile,OutFile
	   CHARACTER*20  SchmidtFile, SchmidtFile2
	   CHARACTER*20  SenChlFile,SenDoFile

	   CHARACTER*4   primc, secoc, nvarc
	   CHARACTER*20  flnm(7)
	   CHARACTER*20  sename(16),SDname
	   CHARACTER*20  limname(1)
	   CHARACTER*5   LONGW
!IF	   CHARACTER*20  resname 
	   CHARACTER*20  filein	 

!    PROVISIONAL
	   REAL(8) ::  CONS,CONEXP,CONDIF
	   REAL(8) ::  time_begin,time_end
	   REAL(8) ::  caca_Part_I(7),caca_Part_F(7)

!------------------------------------------------
	   REAL*8 redPART,redTP,redTN, REDUCTION_PART, REDUCTION_TP,REDUCTION_TN
!     REAL*8,PARAMETER :: redPART = 1.0-0.40d0, 
!     &					redTP   = 1.0-0.40d0,
!     &                	redTN	= 1.0-0.40d0
!
!    NAME OF THE files	
!
	   fileout =   'Contour'
	   fileout2 =  'Profile'
	   fileout3 =  'Fit'
	   fileout4 =  'TimeS'
	   fileout5 =  'HeatDay'
	   fileout6 =  'Heat'
	   fileout7 =  'LakeNum'
	   fileout8 =  'WorkFile'
	   fileout9 =  'WorkFile2'
	   fileout10 =  'MixFile'
	   fileout11 =  'NutFile'
	   fileout12 = 'InFile'
	   fileout13 = 'OutFile'
	   fileout14 = 'Sen_Chl'
	   fileout15 = 'Sen_DO'

	   fileout16 = 'Sen_THP'
	   fileout17 = 'Sen_PP1'
	   fileout18 = 'Sen_POP'
	   fileout19 = 'Sen_RP'
	   fileout20 = 'Sen_NO3'
	   fileout21 = 'Sen_NH4'

	   fileout22 = 'Sen_PN1'
	   fileout23 = 'Sen_PON'
	   fileout24 = 'Sen_DON'
         
	   fileout25 = 'Lim_Nut_Lig'
   	
	   fileout26 = 'Secchi_Dep'

	   fileblank = 'Blank'
	   filecont =  'Bound'

!
!    ERASE OLD files		
   !	
	   COUNT1 = DELFILESQQ('*.bln')
	   COUNT1 = DELFILESQQ('*.dat')
	   COUNT1 = DELFILESQQ('*.txt')
      COUNT1 = DELFILESQQ('*.wqd')
      COUNT1 = DELFILESQQ('*.sim')
	   DO 301 j = 1,7	
	      COUNT1 = DELFILESQQ(file(j)) 
301	CONTINUE

!
!    Header file with file names and user settings
!
!	
!     OPEN(1,file='Tahoe_0102.hdr',status='OLD')
!     OPEN(1,file='Tahoe_2000.hdr',status='OLD')
!	   OPEN(1,file='Tahoe_9902.hdr',status='OLD')
!	   OPEN(1,file='Tahoe_22yrs.hdr',status='OLD')
  	   OPEN(1,file='TrTahoe.hdr',status='OLD')
!  	   OPEN(1,file='INPyramid.hdr',status='OLD')
!       OPEN(1,file='INPyramid0314.hdr',status='OLD')
!       OPEN(1,file='IG1012.hdr',status='OLD')

!gbs    Read IN *.HDR file WHICH CONTAINS ALL THE INFORMATION 

!	   READ (1,*)REDUCTION_PART, REDUCTION_TP, REDUCTION_TN
	   redPART = 1.0d0-REDUCTION_PART*0.01d0
	   redTP   = 1.0d0-REDUCTION_TP  *0.01d0
	   redTN	= 1.0d0-REDUCTION_TN  *0.01d0

      READ(1,*) CALIBRATION
	   READ(1,*) BULK	
	   READ(1,*) humidity

	   IF((bulk.eq.1).and.(humidity.ne.1)) THEN
	     WRITE(*,*) 'If bulk =1, THEN humidity =1'
	     STOP
	   ENDIF
!	   IF((bulk.eq.0).and.(humidity.ne.0)) THEN
!	      WRITE(*,*) 'If bulk =0, THEN humidity =0'
!	      STOP
!	   ENDIF

	   READ(1,*)ITMPR

	   IF(ITMPR.gt.86400) THEN
	      WRITE(*,*)'WARNING : ITMPR MUST BE EQUAL or LESS THAN 86400'
	      WRITE(*,*) ' ITMPR?'
	      READ(5,*)ITMPR
      ENDIF

!	   IF(ITMPR.ne.86400) THEN
!	      WRITE(*,*)'WARNING: IF ITMPR is not set to 86400 seconds, THEN:'
!	      WRITE(*,*)'1) Nutrient file will not be printed out '
!	      WRITE(*,*)'2) SW in the Heatdaily file may be printed 0.0'
!	      PAUSE
!     ENDIF	

!
!    TIME STEP FOR SUB-DAILY SIMULATION	
!
!	READ(1,*)ITIMES

!
!    NUMBER OF VARIABLES TO BE PRINT
!
      READ(1,*) numero   
	   IF(numero.gt.34)THEN
	      WRITE(*,*)'ERROR NVAR SHOULD BE LESS THAN 34'
	   ENDIF
!
!     ID NUMBER OF VARIABLES TO BE READ(T=1,sal=2,...,37)
!                                              	
	   READ(1,*) (vector(k), k = 1,numero)
!
!     READ SWITCH
!
      READ(1,5)RESNAM
      READ(1,'(A5)')LONGW 
      READ(1,*)NZOOSW  
      READ(1,*)NPARTSW 
      READ(1,*)METSW	 
      READ(1,*)PHEATR
	   READ(1,*)PINSRT
      READ(1,*)PMIXER
      READ(1,*)POUT
      READ(1,*)PLAKE
	   READ(1,*)PWORK
      READ(1,*)PNUTRIT
      READ(1,*)MATLAB
	
	   IF(.NOT. MATLAB.and.numero.gt.1) THEN
	     WRITE(*,*)'WARNING: IF NOT MATLAB OUTPUT,SET TO 1 # of VARIABLES TO BE PRINT'
   !	  PAUSE
	   ENDIF

	   IF(.NOT.PLAKE) THEN
	      WRITE(*,*) 'WARNING: LAKENUMBER SHOULD BE .TRUE. IF      &
         CALIBRATION PROCEDURE IS INTENDET USING	ROUTINE OPTIMA'
	   PAUSE
	   ENDIF
	
	   READ(1,*) nconta
	   IF(nconta.gt.10) THEN
	      WRITE(*,*)'ERROR. NUMBER OF DAYS GREATER THAN MAX. '
	      STOP	 
	   ELSE
	      READ(1,*) (conta(j),j=1,nconta)
	   ENDIF
	
	   DO i = 1,nconta-1
	      IF(conta(i+1).eq.conta(i)) THEN
	         WRITE(*,*)'Warning: repeated day. Please, check *.hdr file'
	      STOP
	      ENDIF

	      IF(conta(i+1).le.conta(i)) THEN
	      WRITE(*,*)'Warning: days not ordered. Please, check *.hdr file'
	      STOP
	      ENDIF
	   ENDDO

	   n_profiles = nconta
!
!     READ ELEVATIONS FOR OUTPUT (M)
!
	   READ(1,*)neleva
!	   IF(neleva.gt.10) THEN
!	      WRITE(*,*)'ERROR. NUMBER OF ELEVATIONS MUST BE LESS THAN 8'
!	      STOP
!	   ELSE
	      READ(1,*)(prof(j),j=1,neleva)
!	   ENDIF
!
!     READ NUMBER OF LAYERS OF OUTPUT
!
	   READ(1,*) n_layers
	   n_layers = int(n_layers)
!
!    HEADER files
!
      DO KK = 1,7
         READ(1,5)flnm(KK)
	      INQUIRE(file=flnm,EXIST=FTEST)
         IF(.NOT.FTEST) THEN
	         WRITE (*,*) 'INPUT file DOES NOT EXIST!'
	         WRITE (*,*)  flnm(KK)
	         STOP
	      ENDIF
5       FORMAT(20A)
	   ENDDO
!
!     NUMBER OF .PRO files	
!
	   READ(1,*) nexp 

	   DO 10 i = 1,nexp!
!-----------CPU time------------------!
         CALL CPU_TIME (time_begin)

	      READ(1,*) NDIA,filein,IDAY,ITIMES !,flnm(4),flnm(5),
!     &	 flnm(1),flnm(2),flnm(6)
!     &	 ,Z_red,DepthMax, Shape		                        ! ALAN PROJECT			   
	      INQUIRE(file=filein,EXIST=FTEST) 
         IF(.NOT.FTEST) THEN
	         WRITE (*,*) 'INPUT PRO file DOES NOT EXIST!'
	         WRITE (*,*) filein
	         STOP
	      ENDIF	 
!----------CONTROL OF MATLAB OUTPUT-------------
!        IF(MATLAB) GOTO 2099
         prim  = int(i/10)          !tenth digit
	      seco  = i - int(i/10)*10   !single digit
	      primc = char(prim+48)      !character of tenth digit
	      secoc = char(seco+48)      !character of single digit
!        temfile  = trim(adjustl(fileout )) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
!	      perfile  = trim(adjustl(fileout2)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      fitfile  = trim(adjustl(fileout3)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
!	      seriefile= trim(adjustl(fileout4)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      HeatDay  = trim(adjustl(fileout5)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      HeatFile = trim(adjustl(fileout6)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      LakeFile = trim(adjustl(fileout7)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      MixFile  = trim(adjustl(fileout10))// trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      NutFile  = trim(adjustl(fileout11))// trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      InFile   = trim(adjustl(fileout12))// trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      OutFile  = trim(adjustl(fileout13))// trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
   !	   blank    = trim(adjustl(fileblank))// trim(adjustl(primc))//trim(adjustl(secoc)) //'.bln'
   !	   contur   = trim(adjustl(filecont))//  trim(adjustl(primc))//trim(adjustl(secoc)) //'.bln'
	      SDname   = trim(adjustl(fileout24))// trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      SchmidtFile  = trim(adjustl(fileout8)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      SchmidtFile2 = trim(adjustl(fileout9)) // trim(adjustl(primc))//trim(adjustl(secoc)) // '.txt'
	      limname(1)= trim(adjustl(fileout23))//'.txt'
      
         IF(MATLAB) GOTO 2099
	 
!---------OPEN FITTING file--------------------------!
	      OPEN  (2,file=fitfile,status='unknown')
	      WRITE (2,6116)
6116     FORMAT(' DAY         ','ERROR   ','ERROR_E  ','ERROR_H  ',' ERROR_C ')
	      CLOSE (2)
!---------OPEN file HEATCONTEN.dat WRITING HEAT FLUX and HEAT CONTEN--------------
	      OPEN  (22,file=HeatDay,status='UNKNOWN')      
	      WRITE (22,5550)RESNAM
5550     FORMAT(' DAYLY SURFACE FLUXES W/M2 FOR ',A20)     
	      WRITE (22,5555)
5555     FORMAT(2X,'DAY',9X,'SWR',7X,'LWRin',6X,'LWRout',7X,'LWnet',7X,'EVAPO',6X,'CONDUC',  &
             5X,'SWR+LWR',4X,'EVA+COND',2X,'NET_FLUXES',1X,'T_Sur',1X,'T_Bot') 	
         CLOSE (22)

	      OPEN  (12,file=HeatFile,status='UNKNOWN')      
	      WRITE (12,5050)RESNAM
5050     FORMAT(' INSTANT SURFACE FLUXES W/M2 FOR ',A20)     
!	      WRITE (12,5851)
5851     FORMAT(4x,'ICLOCK',' LWnet ',2x,' EVAP ',2x,' Cond ',2x,' LWout',                 &
            2x,' LWin ',2x,' LWout ',2x,' SW ',2x,' Rad ',2x,'Total')
         CLOSE (12) 

2099     CONTINUE      
!---------Writes the header of the first line of Secchi Depth file------------------------
	      OPEN  (12,file=SDname,status='UNKNOWN')
	      WRITE (12,1265)
1265	   FORMAT(2X,'Year',2X,'Jday',X,'SecDepth',4X,'A_Tahoe',4X,         &
         'B_Tahoe',3X,'Kd_Tahoe',4X,'Bsedime',5X,'Chloro',8X,'POP',8X,    &
         'PON',2X,'NO3_Tahoe',2X,'NH4_Tahoe',2X,'THP_Tahoe',11X,          &
         'Parts1',11X,'Parts2',11X,'Parts3',11X,'Parts4',11X,             &
         'Parts5',11X,'Parts6',11X,'Parts7')
 	      CLOSE (12)
!
!        DATA file and WRITE ON BINARY files
!
         CALL MOD1(workfile,ndia,filein,flnm,LONGW,NZOOSW,NPARTSW,         &
     	      METSW,redPART,redTP,redTN)

	      OPEN(18,file=flnm(3),FORM='FORMATTED',status='OLD')
         READ(18,*) base
         READ(18,*) CRL, CRLNGTH					!Elevation of CREST OF THE RESERVOIR/LAKE and length of crest    
  	      REWIND(18)
	      CLOSE(18)
	  	
!gbs	   OPEN(18,file='Tahoe_new2.fiz',status='OLD')
!gbs	   READ(18,*)base 
!gbs	   READ(18,*) CRL, CRLNGTH
!gbs	   CLOSE(18)
!
!        PERFORM SIMULATION PHYSICS and WQ
!
!        CALL SIMWQ(workfile,HeatDay,HeatFile,LakeFile,              &
!     	             NutFile,MixFile,InFile,OutFile)

         CALL SIMWQ(workfile,HeatDay,HeatFile,LakeFile,SchmidtFile,SchmidtFile2, &                
          NutFile,MixFile,InFile,OutFile,sename,i_nexp,limname,SDname,redPART,redTP,redTN)	
!--------------MATLAB OUTPUT files-----------------------------------    
	      IF(MATLAB) THEN	
            CALL MOD7(workfile,numero,vector, FILEIN)
            GOTO 3099					
         ENDIF 
!    CONTOURS OF THE SIMULATED VARIABLES. CONTOUR, blank and PROFILE file
!    ARE ALL NEEDED TO BUID THE CONTOUR PLOT. ASUMES VARIABLE NVAR = FIRST
!    VARIABLE OF VECTOR
!
	      fizfile   = flnm(3)
	      fieldfile = flnm(7)
	
	      DO jj = 1,numero
	         nvar = vector(jj)	    
	         IF(jj.ge.10) THEN
	            dec=int(jj/10)
	            uni=mod(jj,10)
		         nvarc=char(dec+48)//char(uni+48)
	         ELSE
		         nvarc = char(jj+48)
	         ENDIF 
            temfile = trim(trim(adjustl(fileout)) // trim(adjustl(primc))        & 
                  //trim(adjustl(secoc)) //'_'//trim(adjustl(nvarc)) // '.txt')

	         perfile   = trim(adjustl(fileout2)) // trim(adjustl(primc))          &
                  //trim(adjustl(secoc)) //'_'//trim(adjustl(nvarc)) // '.txt'

	         seriefile = trim(adjustl(fileout4)) // trim(adjustl(primc))          &
                  //trim(adjustl(secoc)) // '_'//trim(adjustl(nvarc)) // '.txt'

	         blank     = trim(adjustl(fileblank))//trim(adjustl(primc))           &
                  //trim(adjustl(secoc))//'_'//trim(adjustl(nvarc)) // '.bln'

	         contur    = trim(adjustl(filecont))//trim(adjustl(primc))            &
                  //trim(adjustl(secoc))//'_'//trim(adjustl(nvarc)) // '.bln'

!	
!    WRITES ON THE FIRST LINE OF BOUNDARY files THE NUMBER OF PAIR POINTS and FLAG 1
!    VALID FOR SURFER and MATLAB NO AUTOMATIC
!
 	         OPEN  (12,file=contur,status='unknown')
	         WRITE (12,113)ndia,',',1
113         FORMAT(I6,A2,I2)
 	         CLOSE (12)
!
!    WRITES ON THE FIRST LINE OF blank files THE NUMBER OF PAIR POINTS and FLAG 1
!    VALID FOR SURFER and MATLAB NO AUTOMATIC
!
            OPEN  (12,file=blank,status='unknown')
	         WRITE (12,163)ndia +3,',',1
163         FORMAT(I6,A2,I2)
 	         CLOSE (12)
!
!quim WRITES ON THE FIRST LINE OF SERIETOT file THE DEPTHS AT WHICH THE VARIABLE 
!quim VALID FOR MICROCALC and MATLAB NO AUTOMATIC
!
	         OPEN  (12,file=seriefile,status='UNKNOWN')
!	         WRITE (12,263)(prof(lf),lf=1,neleva)      
263         FORMAT('jday ',<neleva>(F10.1,4X))
            CLOSE (12)
     
	         CALL MOD5(fieldfile,seriefile,contur,blank,workfile,fizfile,         &
     	         temfile,perfile,TOUT,NVAR,prof,PCONTA,n_profiles,conta,n_layers,neleva,FILEIN)
         ENDDO	! jj
!
!    CALCULATES FITTING FUNCTION FOR THE GENETIC ALGORITHM, 
	      CALL MOD6(NVAR,workfile,fieldfile,fitfile,ERROR,ERROR_E,ERROR_H,ERROR_C)!			         
!    WRITE SELECTED MEASURED VARIABLE INTERPOLATED AT THE depth OF THE SIMULATED VARIABLE
!
	      CALL MOD8(NVAR,fieldfile,prof,neleva)                 
3099     CONTINUE
!    REMOVE OLD files
	      DO 30 j = 1,6	
	         COUNT1 = DELFILESQQ(file(j))
30	      CONTINUE
!------CPU time------------------    
         OPEN(110,file='CPU_Time.tim',status='UNKNOWN',access='append')
         CALL CPU_TIME ( time_end )
         WRITE(110,*) 'CPU time (seconds): ',time_end - time_begin
         CLOSE(110)
10    CONTINUE

	   CLOSE(1)
	   CLOSE(4)
!	   CLOSE(2)
	   CLOSE(3)
      STOP
      END PROGRAM DLM1DWQ  				       
!********************************************************************************

      SUBROUTINE MOD1(BITINI,NDIA,FILEIN,FLNM,LONGW,NZOOSW,NPARTSW,    &
                METSW,reducn_part,reducn_TP,reducn_TN)

!
!     MOD1 DLM-WQ VERSION V1.7 TAHOE LAKE PROJECT UC-DAVIS 2000
!

!     PROGRAM PROCESSES ALL THE INITIAL DATA
!     UNITS: VOLUME=10**3 M**3, AREA=10**3 M**2, depth=METRES

! a(20000)=AREA OF EACH LAYER DETERMINED BY LINEAR INTERPOLATION
! AEM=
! AKH=
! ALPHA(I3)=HALF ANGLE OF STREAM
! an(100)=INTERPOLATION COEFFICIENT FOR VOLUME
! ar(100)=AREA OF EACH LAYER
! arlast=LAST AREA ENTRY
! bn(100)=INTERPOLATION COEFFICIENT FOR AREA
! base=BOTTOM ELEVATION OF RESERVOIR
! K_SOD=BIOLOGICAL SEDIMENT OXYGEN DEMAND (G M**-2 d**-1)
! bufinf=RUNTIME INFLOW DATA CHECKING
! bufsal=SALINITY RUNTIME INFLOW DATA CHECKING
! buftem=TEMPERATURE RUNTIME INFLOW DATA CHECKING
! bufwq=water quality RUNTIME INFLOW DATA CHECKING
! bufcf=PARTICLE RUNTIME INFLOW DATA CHECKING
! CDRAG(I3)=STREAMBED DRAG COEFFICIENT
! cf=ARRAY FOR 7 PARTICLE SIZES
! CFDOWN=DOWNFLOW PARTICLE ARRAY
! cfinf=INFLOW PARTICLE ARRAY
! CFINS=INFLOW PARTICLES
! CHANGE=LOGICAL EXPRESSION FOR CHANGES
! constbo=TEMPERATURE-DEPENDENT ORGANIC DECAY RATE (DAY**-1)
! constm=TEMPERATURE-DEPENDENT MORTALITY (DAY**-1)
! CONSTNH=TEMPERATURE-DEPENDENT NITRIFICATION (DAY**-1)
! constr=TEMPERATURE-DEPENDENT RESPIRATION (DAY**-1)
! COPDRP=CONVERSION RATE OF ORGANIC P TO DRP (DAY**-1)
! CONNH3=CONVERSION RATE OF ORGANIC N TO NH3 (DAY**-1)
! CRL=CREST ELEVATION OF RESERVOIR
! CK=CONVECTIVE OVERTURN
! CS=SHEAR EFFICIENCY
! CT=UNSTEADY EFFECTS
! d(20000)=depth INTERVALS OF 0.1M
! dadz(20000)=GRADIENTS OF AREA BETWEEN 0.1M LAYERS
! DDOWN(I3,maxpar)=
! densy(7)=DENSITY OF 3 PHYTOPLANKTON GROUPS, a TOTAL P, TOTAL N,PART FE and MN
! and BOD, and ZOOPLANKTON
! denscf=DENSITIES OF 7 PARTICLE ARRAYS
! DEPF(maxdata)=FISH depth FOR ZOOPLANKTON PART OF PROGRAM (M)
! depth=1-d depth ARRAY
! DIFF=DIFFUSIVITY OF SUBSTANCE
! DIINS(I3,maxpar)=INFLOW DENSITY
! DLWST(I3)=
! DMAX=MAXIMUM LAYER THICKNESS
! DMIN=MINIMUM LAYER THICKNESS
! DOLD(I3,maxpar)=
! DRW(I2)=OUTFLOW VOLUMES
! dtmax=SURFACE depth
! dv(20000)=GRADIENTS OF VOLUME BETWEEN 0.1M LAYERS
! ETA=WIND STIRRING
! etwat=BACKGROUND light ATTENUATION (M**-1)
! FERED=RATE OF REDUCTION OF FE III
! FEOXY=RATE OF OXIDATION OF FEII
! filein=INPUT file
! FISH(100)=NUMBER OF FISH/1000M3
! flnm(7)=INPUT files: X.OUT, X.IN, X.FIZ, X.MET, X.TAB, X.WAT, INITIAL,X.FLD
! FLOINF(I3)=
! FO=
! FSUM=
! FTEST=LOGICAL EXPRESSION FOR DETERMINING IF file EXISTS
! gromax=MAXIMUM GROWTH RATE OF PHYTOPLANKTON (DAY**-1)
! HFLOW(I3)=
! HSCP=HALF SATURATION CONSTANT FOR PHOSPHORUS UPTAKE (MG M**-3)
! hscn=HALF SATURATION CONSTANT FOR NITROGEN UPTAKE (MG M**-3)
! HSCSI=HALF SATURATION CONSTANT FOR SILICA UPTAKE (MG M**-3)
! HSCNHNO=PREFERENTIAL AMMONIUM UPTAKE FACTOR
! HSCZZ=HALF SATURATION CONSTANT FOR ZOOPLANKTON GRAZING (MG CHLA M**-3)
! HTSAVE=
! i,j,zz=VARIABLE INTEGERS
! ID=INDICATOR FOR INFLOW START TIME ERROR
! idbuf(I3)=
! ICNT(I3)=
! IDIF=NUMBER OF DIFFUSING SUBSTANCES
! IDIM1=NUMBER OF LAYERS
! IDTEST(I3)=TEST THAT ALL START DATES ARE OK
! ij=
! IO=9
! I2=9
! I3=100
! INPAR(I3,maxpar)=
! INUNIT=INPUT file INDICATOR
! iseqb(I3)=INDICATOR FOR INFLOW END TIME ERROR
! ISEQ1=
! ISEQ2=
! ISEQ3=
! ISTOR=NUMBER OF STORAGE TABLE LAYERS
! jday1=DAY TO START SIMULATION
! jyear1=YEAR TO START SIMULATION
! jcheck=END DATE OF SIMULATION
! JDLYD=DAY OF START DATE OF SIMULATION
! JDLYY=YEAR OF START DATE OF SIMULATION
! JD0LST=
! JD0STR=
! JLAST=JULIAN DATE FOR END OF SIMULATION
! JLASTD=DAY FOR END OF SIMULATION
! JLASTL=LAST DAY OF SIMULATION WHEN THERE IS a LEAP YEAR
! JLASTY=LAST YEAR OF SIMULATION WHEN THERE IS a LEAP YEAR
! JMAX=
! jstart=JULIAN START DATE FOR SIMULATION
! JSTRTD=DAY FOR START OF SIMULATION
! JSTRTL=
! JSTRTY=YEAR FOR START OF SIMULATION
! JTEST1=	Year (Normal)
! JTEST2=	Day
! JTEST3=	Year (Leap)
! karo=STORAGE AREA ADJUSTMENT IF > 0
! kar1=FIRST LAYER WITH a POSITIVE AREA
! KO=SWITH INDICATING ERROR IN water quality file
! kstoro=STORAGE VOLUME ADJUSTMENT IF < 0
! kstor1=FIRST LAYER WITH a POSITIVE STORAGE
! lanext=TEMPORARY VARIABLE FOR INTERPOLATING AREA
! latit=LATITUDE OF RESERVOIR
! LC=LENGTH OF RESERVOIR AT CREST
! light=light FOR SATURATION
! LK= DAYS COUNTER
! LONGW=TYPE OF LONGWAVE RADIATION MEASUREMENT
! lvnext=TEMPORARY VARIABLE FOR INTERPOLATING VOLUME
! LWIND=TYPE OF LONGWAVE RADIATION
! maxchl=NUMBER OF ARRAYS TAKEN BY CHLOROPHYLL a
! MSTEP=
! MNRED=RATE OF REDUCTION OF MN IV
! MNOXY=RATE OF OXIDATION OF MN II
! nchl=NUMBER OF ALGAL GROUPS (AS CHLOROPHYLL) CONSIDERED (1 or 3)
! ndata=NUMBER OF depth/AREA/VOLUME DATA POINTS
! NDAYS=NUMBER OF DAYS TO SIMULATE
! NEWDAZ=
! NOINS(I3)=
! NINMIN=MINIMUM INTERNAL NITROGEN CONCENTRATION (MG/M3)
! NINMAX=MAXIMUM INTERNAL NITROGEN CONCENTRATION (MG/M3)
! NLEAPL=
! ns=NUMBER OF LAYERS
! nstor=NUMBER OF DATA POINTS
! numdif=TOTAL NUMBER OF DIFFUSING SUBSTANCES
! NUMOUT=NUMBER OF OUTFLOWS
! numriv=NUMBER OF INFLOWS
! np=
! nx=
! nz=
! OLDSL=
! OLEN=BASIN LENGTH AT THE OUTLET
! OLEV(I2)=BASIN WIDTH AT THE OUTLET
! OS=FIXED LAYER SIZE FOR FISH
! OWID=WIDTH OF RESERVOIR
! OUNIT=OUTPUT file INDICATOR
! PARTSW=PARTICLE SWITCH
! PHEATR=PRINT HEATR file OPTION
! PHI(I3)=STREAMBED SLOPE
! PINMIN=MINIMUM INTERNAL PHOSPHORUS CONCENTRATION (MG/M3)
! PINMAX=MAXIMUM INTERNAL PHOSPHORUS CONCENTRATION (MG/M3)
! PHEATR=PRINT HEATR file OPTION
! PINSRT=PRINT INFLOW file OPTION
! PMIXER=PRINT MIXER file OPTION
! POUT=PRINT OUTFLOW file OPTION
! prodat=PROFILE STARTING DATE
! QDOWN=
! QINS(I3,maxpar)=INFLOW VOLUME
! QUERY=SWITCH FOR CORRECTION TO START DATE
! RAIN=DAILY RAINFALL
! RDSAT=SATURATION LEVEL FOR REACTIVE DISTANCE OF FISH (M)
! RESNAM=RESERVOIR NAME
! sal=1-d SALINITY ARRAY
! SALINF=INFLOW SALINITY
! SATFZ1,SATFZ1=GAIN RATE WHERE ZOOPLANKTON MINIMISE PREDATION
! (MICROGC/ZOOPLANKTON/DAY)
! SDOWN=DOWNFLOW SALINITY
! sedthp=SEDIMENT RELEASE RATE OF DRP (G/M2/d)
! sednh3=SEDIMENT RELEASE RATE OF NH3 (G/M2/d)
! sedno3=SEDIMENT RELEASE RATE OF NO3 (G/M2/d)
! SEDFE=SEDIMENT RELEASE RATE OF FE
! SEDMN=SEDIMENT RELEASE RATE OF MN
! sedtem=SEDIMENT TEMPERATURE MULTIPLIER FOR RELEASE RATE
! SINS(I3,maxpar)=INFLOW SALINITY
! SRAT=LONGWAVE RADIATION
! STARV1,STARV2=LOSS RATE WHERE ZOOPLANKTON MAXIMIZE INGESTION
! (MICROGC/ZOOPLANKTON/DAY)
! stlast=LAST VOLUME ENTRY
! stor(100)=CUMULATIVE VOLUME OF EACH LAYER
! SVPD=VAPOUR PRESSURE
! SW=SHORTWAVE RADIATION
! T4=TEMPERATURE
! TDOWN=DOWNFLOW TEMPERATURE
! teminf=INFLOW TEMPERATURE
! temp=1-d TEMPERATURE ARRAY
! thetabo=TEMPERATURE ADJUSTMENT FOR DECAY OF ORGANIC MATTER
! thetanh=TEMPERATURE ADJUSTMENT FOR NITRIFICATION
! thetat=TEMPERATURE ADJUSTMENT FOR GROWTH RATE
! THR=
! TI=
! TIMEI=
! TINS(I3,maxpar)=INFLOW TEMPERATURE
! TOTIN(I3)=
! UNMAX=MAXIMUM NITROGEN UPTAKE RATE (DAY**-1)
! UPMAX=MAXIMUM PHOSPHORUS UPTAKE RATE (DAY**-1)
! UAV=
! UF=
! UI=
! U6=WIND SPEED
! v(20000)=VOLUME OF EACH LAYER DETERMINED BY LINEAR INTERPOLATION
! wqdown=DOWNFLOW water quality ARRAY
! wqinf=INFLOW water quality ARRAY
! WQINS(I3,27,maxpar)=INFLOW water quality
! wqual=2-d water quality ARRAY (2ND DIMENSION CONSISTS OF (IN
! ORDER: CHLA AS GREEN ALGAE and DINOFLAGELLATES, CHLA
! AS BLUE-GREEN ALGAE, CHLA AS DIATOMS, OXYGEN, DISSOLVED REACTIVE
! PHOSPHORUS, INTERNAL PHOSPHORUS IN THREE ALGAL GROUPS,
! NITRATE-NITROGEN, AMMONIACAL NITROGEN, INTERNAL
! PHOSPHORUS IN BLUE-GREEN ALGAE,INTERNAL PHOSPHORUS IN DIATOMS,
! TOTAL PHOSPHORUS, INTERNAL NITROGEN IN ALL ALGAL GROUPS, TOTAL NITROGEN,
! SILICA, 2 GROUPS OF ZOOPLANKTON
! vcrl=
! VDA(ISTOR)=
! vmin=MINIMUM LAYER VOLUME
! vmax=MAXIMUM LAYER VOLUME
! WC=WIDTH OF RESERVOIR AT CREST
! wsl(100)=ELEVATION OF EACH LAYER
! X,y=TEMPORARY VARIABLES FOR CALCULATING STORAGE AT CREST
! ZOO=ZOOPLANKTON NUMBERS ARRAY
! ZOOSW=ZOOPLANKTON SWITCH
! ZRESP=base LEVEL ZOOPLANKTON RESPIRATION RATE (DAY-1)
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
      
      INTEGER*4 I2
      INTEGER*4 I3
		INTEGER*4 NDIA
!      INTEGER*4 maxpar, maxseg, maxinf
!
!    Changed for Lake Tahoe Project new definition of MAXINF
!
      PARAMETER(I2=25,I3=100) !maxpar=500, maxseg=10, maxinf=100)

      INTEGER*4 IO
 !     INTEGER*4 ICNT(I3), numseg(I3)
      INTEGER*4 ID
      INTEGER*4 idbuf(I3)
      INTEGER*4 IDIF
      INTEGER*4 IDIM1
      INTEGER*4 i,II
      INTEGER*4 ij,jhj
!IF      INTEGER*4 INPAR(I3,maxpar)
      INTEGER*4 INUNIT
      INTEGER*4 ISEQ1(3)
      INTEGER*4 ISEQ2(3)
      INTEGER*4 ISEQ3(3)
      INTEGER*4 ISEQB1(I3)
      INTEGER*4 ISEQB2(I3)
      INTEGER*4 ISEQB3(I3)
      INTEGER*4 ISTOR
      INTEGER*4 j,jj
!IF      INTEGER*4 jday1
      INTEGER*4 jcheck
      INTEGER*4 JDLDH
      INTEGER*4 JDLYD,JDD(I3)
      INTEGER*4 JDLYY,JYY(I3)
      INTEGER*4 JD0LST
      INTEGER*4 JD0STR
      INTEGER*4 JLASTD
      INTEGER*4 JLAST
	   INTEGER*4 JLASTL
      INTEGER*4 JLASTY
      INTEGER*4 JMAX
      INTEGER*4 jstart
      INTEGER*4 JSTRTD
      INTEGER*4 JSTRTL
      INTEGER*4 JSTRTY
      INTEGER*4 JTEST1
      INTEGER*4 JTEST2
      INTEGER*4 JTEST3
!IF      INTEGER*4 jyear1
      INTEGER*4 karo,kar1
      INTEGER*4 k, KO,KK,KOK
      INTEGER*4 kstoro,kstor1
      INTEGER*4 lanext
      INTEGER*4 LK
      INTEGER*4 lvnext
      INTEGER*4 LWIND
!IF	   INTEGER*4 NMETSW

!IF      INTEGER*4 maxchl
!IF      INTEGER*4 MSTEP
!IF      INTEGER*4 nchl
      INTEGER*4 ndata
!IF      INTEGER*4 NDAYS
      INTEGER*4 NEWDAZ
      INTEGER*4 NLEAPL
      INTEGER*4 np
!IF      INTEGER*4 NOINS(I3)
!IF      INTEGER*4 ns
      INTEGER*4 nstor
!IF      INTEGER*4 numdif
!IF      INTEGER*4 NUMOUT
      INTEGER*4 numriv
      INTEGER*4 nx
      INTEGER*4 nz
!IF      INTEGER*4 PARTSW
!IF      INTEGER*4 prodat
!IF      INTEGER*4 ZOOSW
      INTEGER*4 zz

      REAL(8), PARAMETER ::  CSRAT=0.99d0,TENM6=1.0d-6
!
!     MORE ROOM FOR TAHOE LAKE
!
	   PARAMETER(IDIM1=500,ISTOR=20000,IDIF=37)
!
!     MAXCHLA 2 GROUPS
!
      PARAMETER (IO=9)

!IF      REAL(8) ::  AEM
!IF      REAL(8) ::  AKH
!IF      REAL(8) ::  ALPHA(I3, maxseg)
!gbs  REAL(8) ::  ALPHA(I3
      REAL(8) ::  arlast,stlast,dtmax !IF BASE
      REAL(8) ::  bufinf(I3)
      REAL(8) ::  bufsal(I3)
      REAL(8) ::  buftem(I3)
      REAL(8) ::  bufwq(I3,28),bufwq_1(I3,28),bufwq_2(I3,28)
      REAL(8) ::  bufcf(I3,7),bufcf_1(I3,7)

!gbs  REAL(8) ::  CDRAG(I3)
!IF      REAL(8) ::  CDRAG(I3, maxseg)
!IF      REAL(8) ::  CFDOWN(I3,7,maxpar)
!IF      REAL(8) ::  cf(IDIM1,7)
!IF	   REAL(8) ::  cfinf(I3,7)
!IF	   REAL(8) ::  CFINS(I3,7,maxpar)
!IF      REAL(8) ::  CK
!IF      REAL(8) ::  CRL
!IF	   REAL(8) ::  CRLNGTH
!IF      REAL(8) ::  CS
!IF      REAL(8) ::  CT
!IF	   REAL(8) ::  coag
!IF      REAL(8) ::  DDOWN(I3,maxpar)
!IF      REAL(8) ::  depth(IDIM1)
!IF      REAL(8) ::  DIFF(IDIF)
!IF      REAL(8) ::  DIINS(I3,maxpar)
!IF      REAL(8) ::  DLWST(I3)
!IF      REAL(8) ::  DMAX
!IF      REAL(8) ::  DMIN
!IF      REAL(8) ::  DOLD(I3,maxpar)
!IF      REAL(8) ::  DRW(I2)
!IF      REAL(8) ::  ETA
!IF      REAL(8) ::  FLOINF(I3)
!IF      REAL(8) ::  FO
!IF      REAL(8) ::  FSUM
!IF      REAL(8) ::  HFLOW(I3)
!IF      REAL(8) ::  HTSAVE
      REAL(8) ::  latit
!IF      REAL(8) ::  LC
!IF      REAL(8) ::  OLDSL
!IF      REAL(8) ::  OLEN(I2)
!IF      REAL(8) ::  OLEV(I2)
!IF      REAL(8) ::  OWID(I2)
!gbs     REAL(8) ::  PHI(I3)
!!!wef added seglngth, bgnwdth, bgnele, and bgnarea to account for variables needed by plunge.for
!IF	   REAL*8 seglngth(I3, maxseg)
!IF	   REAL*8 bgnwdth(I3, maxseg)
!IF	   REAL*8 bgnele(I3, maxseg)
!!!
!IF      REAL(8) ::  QDOWN(I3,maxpar)
!IF      REAL(8) ::  QINS(I3,maxpar)
!IF      REAL(8) ::  RAIN
!IF	   REAL(8) ::  RH
!IF      REAL(8) ::  sal(IDIM1)
!IF      REAL(8) ::  SALINF(I3)
!IF      REAL(8) ::  SDOWN(I3,maxpar)
!IF      REAL(8) ::  SINS(I3,maxpar)
!IF      REAL(8) ::  SRAT
!IF      REAL(8) ::  SVPD
!IF      REAL(8) ::  SW
!IF      REAL(8) ::  T4
!IF	   REAL(8) ::  thetabs
!IF      REAL(8) ::  TDOWN(I3,maxpar)
!IF      REAL(8) ::  teminf(I3)
!IF      REAL(8) ::  temp(IDIM1)
!IF      REAL(8) ::  THR
!IF      REAL(8) ::  TI
!IF      REAL(8) ::  TIMEI
!IF      REAL(8) ::  TINS(I3,maxpar)
!IF      REAL(8) ::  TOTIN(I3)
!IF      REAL(8) ::  UAV
!IF      REAL(8) ::  UF
!IF      REAL(8) ::  UI
!IF      REAL(8) ::  U6
!IF      REAL(8) ::  vmax
!IF      REAL(8) ::  vmin
!IF      REAL(8) ::  wqual(IDIM1,28)
!IF      REAL(8) ::  wqdown(I3,28,maxpar)
!IF      REAL(8) ::  wqinf(I3,28)
!IF      REAL(8) ::  WQINS(I3,28,maxpar)
!
!    MORE ROOM FOR TAHOE LAKE 
!
      REAL(8) ::  maxdata
	   PARAMETER (maxdata = 500)
	   REAL(8) ::  wsl(20000),stor(20000),ar(20000),an(20000),bn(20000) !@
      REAL(8) ::  v(20000),dv(20000),d(20000),DVDA(20000)
!IF      REAL(8) ::a(20000),dadz(20000),
!IF      REAL(8) ::  vcrl
      REAL(8) ::  VDA(20000)
!IF     REAL(8) ::  WC
      REAL(8) ::  x,y

      LOGICAL*4 IDTEST(I3)
      LOGICAL*4 NPARTSW,METSW
      LOGICAL*4 NZOOSW

      CHARACTER*5  LONGW
      CHARACTER*20 flnm(7)
      CHARACTER*20 filein  
      CHARACTER*1  QUERY
!IF      CHARACTER*20 RESNAM
!IF      CHARACTER*20 RIVNAM(100)   
!
!    SPECIFIC water quality PARAMETERS
!
      REAL(8) ::  alpha_P 
!IF	   REAL(8) ::  K_SOD
!IF      REAL(8) ::  constbo
!IF	   REAL(8) ::  constnh
!IF      REAL(8) ::  constm(maxchl)
!IF      REAL(8) ::  constr(maxchl)
!    CHANGED FROM 7 TO 1
!    CHANGED FOR ALGAE SETTLING...

!IF	   REAL(8) ::  densy(1)
	   REAL(8) ::  Depthref
!IF	   REAL(8) ::  org_setl_vel
!IF	   REAL(8) ::  org_part_fac

!IF      REAL(8) ::  denscf(7)
!IF      REAL(8) ::  etwat
!IF	   REAL(8) ::  etca(maxchl)
!IF      REAL(8) ::  etpart(7)
!IF	   REAL(8) ::  hscnit
!IF	   REAL(8) ::  hscdo
!IF	   REAL(8) ::  KH_OSOD 
!IF      REAL(8) ::  light(maxchl)
!IF      REAL(8) ::  gromax(maxchl)
!IF      REAL(8) ::  hscn(maxchl)
	   REAL(8) ::  k1_Res         ! NEW SALTON SEA (Feb 2004)
	   REAL(8) ::  k2_Res         ! NEW SALTON SEA (Feb 2004)
       REAL(8) ::  porosity
!IF	   REAL(8) ::  thetabo
!IF      REAL(8) ::  thetanh
!IF      REAL(8) ::  thetat(maxchl)
!IF      REAL(8) ::  sedthp
!IF      REAL(8) ::  sednh3
!IF      REAL(8) ::  sedno3
!IF      REAL(8) ::  sedtem
	   REAL(8) ::  UCritic
!wef removed the hiding from next 4 variables
!IF      REAL*8 NINMIN(MAXCHL)
!IF      REAL*8 PINMIN(MAXCHL)
!IF      REAL*8 NINMAX(MAXCHL)
!IF      REAL*8 PINMAX(MAXCHL)
!wef removed the hiding from next 4 variables
!IF      REAL*8 SEDPO4
      REAL*8 SEDEFE
!IF      REAL*8 SEDMN  
!IF      REAL*8 ZRESP
!MB *********************************************
!MB MENU MODIFICATIONS BEGIN HERE

      INTEGER*4 ULIS, UFLD                                           
	   INTEGER*4 UIN1,URUN,UMET,URIV,UOUT,UPAR,UWQD
      CHARACTER*20 FLIS, FIN1, FFLD, FRUN,FMET,FRIV,FOUT, FPAR, FWQD                                          
      CHARACTER*20 BITINI       
      CHARACTER*4  EIN1
      CHARACTER*4  EPAR, ERUN, EMET,ERIV,EOUT,EWQD
      CHARACTER*80 MESS                                                 
      LOGICAL*1 SUCCESS                                                         
      INTEGER*4 FJDAY
	   INTEGER*4 caca 
!gbs***************************  
!********************************
	   INTEGER*4 :: ps, wq, SL, total_sublayers, sublayers
	   REAL(8) ::  in_depth(IDIM1)
	   REAL(8) ::  in_temp(IDIM1)
	   REAL(8) ::  in_sal(IDIM1)
	   REAL(8) ::  in_cf(IDIM1,7)
	   REAL(8) ::  in_wqual(IDIM1,100)
!************sat_oxygen**********************
	   REAL*8 DOsat, atm_press, part_press_wat_vap, DOE_theta,IZMF

!gbs*******************************  
!gbs 17Jan06	 
!IF	   REAL*8 zu,zq,ztair  
!
!    NEW MODEL VARIABLES 
!
	   LOGICAL FTEST
      CHARACTER*10 THESTAMP
!IF	   INTEGER*4 humidity
!IF	   REAL(8) ::  ypchla(maxchl)
!IF	   REAL(8) ::  ynchla(maxchl)
!IF	   REAL(8) ::  ysichla(maxchl)
!IF	   REAL(8) ::  czoo(maxchl)
!IF	   REAL(8) ::  halfn(maxchl)
!IF	   REAL(8) ::  halfp(maxchl)
!IF	   REAL(8) ::  cn1,cn2,cn3,cn4,cn5
!IF	   REAL(8) ::  cp1,cp2,cp3
	   REAL(8), PARAMETER :: one = 1.0d0
!IF	   REAL(8), PARAMETER :: alg_min = 0.19d0
	   REAL*8 reducn,Ured,NUred,CEred,total_reducn,total_year
	   REAL*8 redn_P, redn_N, reducn_part, reducn_TP, reducn_TN 
	   REAL*8 UP,NUP,CEP,UN,NUN,CEN
	   INTEGER*4 RIV_ID(I3),RIV
	   REAL*8 Uflow_prct(I3),NUflow_prct(I3)
!IF     SAVE FSUM,UI,UAV,TIMEI,THR,TI,AEM,FO,OLDSL,UF,HTSAVE,MSTEP,IDTEST
!IF     SAVE QDOWN,TDOWN,SDOWN,wqdown,CFDOWN,DDOWN,DOLD,QINS,TINS,SINS,         &
!IF      WQINS,CFINS,DIINS,HFLOW,TOTIN
      SAVE IDTEST                                                                                   
!
!    Logic File Units
! 
      PARAMETER (INUNIT=15)
      PARAMETER (UIN1=38, ULIS=39, UFLD=40,URUN=41)
      PARAMETER (UPAR=42,UWQD=43)
      PARAMETER (UMET=44, URIV=45,UOUT=46)
                                                   
 		NDAYS=NDIA				                                                                                                         
      EIN1 = '.IN1'           
      ERUN = '.RUN'
      EMET = '.CLI'
      ERIV = '.QIN'
      EOUT = '.QOT'
      EPAR = '.PAR'           
      EWQD = '.WQD'                                          
								       
      FLIS = ' '                                                        
      FFLD = ' '       
					   
3   	CONTINUE	  

!
!    Open files and start DATA preparation
!
7     CONTINUE
!
!    DATE and TIME CALLS MADE HERE TO GIVE ALL THE file NAMES an INDIVIDUAL
!     and UNIQUE NAME. EACH file NAME IS ALSO CONSTRUCTED HERE
!
	   CALL TIMESTAMP(THESTAMP)
      FFLD =   'W' // THESTAMP // '.BIN'
	   FPAR =   'W' // THESTAMP // '.PAR'
	   FIN1 =   'W' // THESTAMP // '.IN1'
	   FRUN =   'W' // THESTAMP // '.RUN'
       FMET =   'W' // THESTAMP // '.CLI'
       FRIV =   'W' // THESTAMP // '.QIN'
       FOUT =   'W' // THESTAMP // '.QOT'
	   FWQD =   'W' // THESTAMP // '.WQD'
	   BITINI = 'W' // THESTAMP
!
!    OPEN and READ .PRO file
!
      INQUIRE(file=filein,EXIST=FTEST) 
      IF(.NOT.FTEST) THEN
	      WRITE (*,*) 'INPUT .PRO file DOES NOT EXIST (MOD1)!'
	      WRITE (*,*)	
	      STOP
	   ENDIF

      OPEN(INUNIT,file=filein,FORM='FORMATTED',status='OLD')
 
      READ(INUNIT,*)prodat 
      READ(INUNIT,*)jstart 
      READ(INUNIT,*)vmin	
      READ(INUNIT,*)DMIN
      READ(INUNIT,*)DMAX
      READ(INUNIT,*)ns      
      READ(INUNIT,*)	      
!
!    FILTER FOR NOT YET AVAILABLE OPTIONS
!
11    IF(NZOOSW)THEN	! ZOOPLANKTON 
	      ZOOSW=1
     	   WRITE(*,*)'WARNING: TAHOE DLM-WQ DOES NOT SUPPORT ZOO'
	      STOP
      ELSE
	      ZOOSW=0
      ENDIF

      IF(NPARTSW)THEN	! PARTICLES
	      PARTSW=1
      ELSE
	      PARTSW=0
      ENDIF

      IF(METSW)THEN	! METALS
	      NMETSW=1
     	   WRITE(*,*)'WARNING: TAHOE DLM-WQ DOES NOT SUPPORT METALS'
	      STOP
      ELSE
	      NMETSW=0
      ENDIF
!
!    INICIALIZE WQ and PARTICLE VARIABLES
!    
	   DO i = 1,maxdata
	      DO j = 1, 100
            wqual(i,j)= 0.0d0
	      ENDDO  ! j
	      DO j = 1,7
	         cf(i,j)   = 0.0d0
	      ENDDO  ! j
	   ENDDO   ! i
!
!   	READ INITIAL CONDITIONS MOD5 i MOD4, MUST HAVE THE SAME READING STRUCURE
!	
      DO 14 i=1,ns
	 						! METALLS OFF, PARTICULES ON, ZOO OFF
	      IF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1)THEN
	         READ(INUNIT,*)depth(i),temp(i),sal(i),wqual(i,1),wqual(i,2),  &
            wqual(i,3),wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),        & 
            wqual(i,8),wqual(i,9),wqual(i,10),wqual(i,11),                 &
            wqual(i,12),wqual(i,13),wqual(i,14),wqual(i,15),wqual(i,16),   &
            wqual(i,17),wqual(i,18),wqual(i,19),wqual(i,20),wqual(i,21),   &
            wqual(i,22),wqual(i,23),wqual(i,24),wqual(i,25),wqual(i,26),   &
            wqual(i,27),wqual(i,28),cf(i,1),cf(i,2),cf(i,3),cf(i,4),       &
            cf(i,5),cf(i,6),cf(i,7)

				!METALLS OFF, PARTICULES OFF, ZOO OFF	
	      ELSEIF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1)THEN
	         READ(INUNIT,*)depth(i),temp(i),sal(i),wqual(i,1),wqual(i,2),  &
            wqual(i,3),wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),        &
            wqual(i,8),wqual(i,9),wqual(i,10),wqual(i,11),                 &
            wqual(i,12),wqual(i,13),wqual(i,14),wqual(i,15),wqual(i,16),   &
            wqual(i,17),wqual(i,18),wqual(i,19),wqual(i,20),wqual(i,21),wqual(i,22),wqual(i,23),wqual(i,24)
	      ENDIF
!
!    Set component oxygen concentrations
!
!2015/08	      wqual(i,4) = 0.0d0
!2015/08	      wqual(i,5) = 0.0d0
!2015/08	      wqual(i,6) = wqual(i,8)
!2015/08       wqual(i,6) = 0.0d0
!2015/08	      wqual(i,7) = 0.0d0
14    CONTINUE
!      STOP	   
!
!    NEW STRUCTURE ONLY ALLOWS 2 GROUPS OF ALGAE
!
      DO i=1,ns			   
!2013/04	      IF(wqual(i,3).ne.0.0d0)THEN 
!2013/04	         WRITE(*,*)'ERROR:TAHOE DLM-WQ ONLY SUPORTS 2 ALGAE GROUPS'
!2013/04	         STOP
!2013/04	      ENDIF
!2013/04	      IF(wqual(i,2).ne.0.0d0)THEN 
!2013/04	         nchl = 2
!2013/04	         GOTO 17                
!2013/04	      ENDIF
	      nchl = 7
!	      wqual(i,12) = 0.0d0
!	      wqual(i,18) = 0.0d0  
	   ENDDO
17    CONTINUE
!
!    Check for accpetable range of WQ parameters
!
	   DO i=1,ns			 
	      DO j=1,23		  
	         IF(wqual(i,j).lt.0.0d0)THEN
	            WRITE(*,*) 'Invalid WQ value in INITIAL PROFILE'
	            WRITE(*,*)i,j, wqual(i,j)     
	            STOP
	         ENDIF
	      ENDDO   ! i
	   ENDDO   ! j
      IF(PARTSW)THEN		 
	      DO i=1,ns			 
	         DO j=1,7		  
	            IF(cf(i,j).lt.0.0d0)THEN
	               WRITE(*,*) 'Invalid Particle value in INITIAL PROFILE'
	               WRITE(*,*)i,j, cf(i,j)     
	               STOP
	            ENDIF
	         ENDDO   ! i
	      ENDDO   ! j
	   ENDIF
      CLOSE(INUNIT)
!
!     OUTPUT PARAMETER VALUES SO THAT USER MAY MODIFY THEM
!
      WRITE(*,24) RESNAM
24    FORMAT(/,' AT PRESENT THE DLM-WQ MODEL IS SET UP FOR ',20A)

      IF(NMETSW.eq.1)THEN
	      IF(nchl.eq.3)THEN
	         WRITE(*,725) 
725         FORMAT(' TO SIMULATE THREE ALGAL GROUPS and METALS ',/)     
	      ELSE
	         WRITE(*,726)
726         FORMAT(' TO SIMULATE one ALGAL GROUP and METALS',/)
	      ENDIF
      ELSEIF(ZOOSW.eq.1.and.PARTSW.eq.1.and.NMETSW.ne.1)THEN
	      IF(nchl.eq.3)THEN
	         WRITE(*,25) 
25          FORMAT(' TO SIMULATE THREE ALGAL GROUPS, ZOOPLANKTON ','and PARTICLES ',/)     
	      ELSE
	         WRITE(*,26)
26          FORMAT(' TO SIMULATE one ALGAL GROUP, ZOOPLANKTON and',' PARTICLES',/)
	      ENDIF
      ELSEIF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1)THEN
	      IF(nchl.eq.3)THEN
	          WRITE(*,27) 
27           FORMAT(' TO SIMULATE THREE ALGAL GROUPS and' ,' PARTICLES',/)
	      ELSE
	         WRITE(*,28)
28          FORMAT(' TO SIMULATE one ALGAL GROUP and PARTICLES',/)
	      ENDIF     
      ELSEIF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1)THEN
	      IF(nchl.eq.3)THEN
	         WRITE(*,29) 
29          FORMAT(' TO SIMULATE THREE ALGAL GROUPS',/)
	      ELSE
	         WRITE(*,30)
30          FORMAT(' TO SIMULATE one ALGAL GROUP',/)
	      ENDIF
      ELSEIF(ZOOSW.eq.1.and.PARTSW.ne.1.and.NMETSW.ne.1)THEN
	      IF(nchl.eq.3)THEN
	         WRITE(*,31) 
31          FORMAT(' TO SIMULATE ZOOPLANKTON and 3 ALGAL GROUPS',/)
	      ELSE 
	         WRITE(*,32)
32          FORMAT (' TO SIMULATE ZOOPLANKTON and 1 ALGAL GROUP',/)
	      ENDIF
      ENDIF
!
!    RUN SPECIFICATIONS
!   	prodat jstart
!
      jyear1 = jstart/1000
      jday1  = jstart-jyear1*1000

      IF(prodat+1 .ne. jstart) THEN 
!	      WRITE(*,33) 
33       FORMAT(/,' *** PROFILE DATE or START DATE IS WRONG ***'/)
!         WRITE(*,34)prodat,jstart,NDAYS
34       FORMAT('INITIAL PROFILE DAY              ',4X,I8,/,      &
                'SIMULATION START DAY             ',4X,I8,/,      &
                'RUN DAYS                         ',6X,I4,/)
!	      STOP 
	   ENDIF
!
!    Check LONGWAVE ENTRY
!
	   IF (LONGW .eq. 'CLOUD'.or.LONGW .eq. 'CLOUD')THEN
	      LWIND = 1
	   ELSEIF(LONGW .eq. 'MLWIN'.or.LONGW .eq. 'MLWIN')THEN
	      LWIND = 2
	   ELSEIF(LONGW .eq. 'NETLW'.or.LONGW .eq. 'NETLW')THEN
	      LWIND = 3
	   ELSE
	      WRITE(*,*)' ERROR IN LONG WAVE TYPE STATED IN INPUT file '
	      STOP
	   ENDIF
      WRITE(*,*) 'GETTING STORAGE TABLE DATA FROM ',flnm(5)
      OPEN(30,file=flnm(5),FORM='FORMATTED',status='OLD')
      READ(30,*)
      READ(30,*)wsl(1),ar(1),stor(1)	
!						      
! CONVERT AREA UNITS TO MILLIONS OF SQUARE METRES
!
      ar(1)  = ar(1)/1000.0d0
      arlast = ar(1)
      stlast = stor(1)
      base   = wsl(1)
      wsl(1) = wsl(1)-base
      kstoro = 0
      karo   = 0

      IF(stor(1).le. 0.0d0) kstoro = kstoro + 1
      IF(ar(1)  .le. 0.0d0) karo   = karo + 1

      ndata	= 1
80    i	= ndata+1
      
      READ(30,*,END=84)wsl(i),ar(i),stor(i)
	   ar(i)	= ar(i)/1000.0d0  
      ndata	= ndata+1
!
! MAKE ALL ELEVATIONS RELATIVE TO RESERVOIR BOTTOM
!
      wsl(i)=wsl(i)-base
!
!  COUNT NUMBER OF ZERO ENTRIES FOR AREA and VOLUME
!
      IF(ar(i).le.0.0d0)   karo   = karo + 1
      IF(stor(i).le.0.0d0) kstoro = kstoro + 1

! CHECK CONSISTENCY OF INPUT DATA

      IF(wsl(i).lt.wsl(i-1).or.ar(i).lt.arlast.or.stor(i).lt.stlast)THEN
	      WRITE(*,82)i
82       FORMAT(' UNEXPECTED DATA READ ON RECORD ',I3,/,       & 
                ' DATA MUST INCREASE MONOTONICALLY')
	      STOP
      ENDIF
 !     d(i)=wsl(i)
 !     a(i)=ar(i)
 !     v(i)=stor(i)
      arlast	= ar(i)
      stlast	= stor(i)
      GOTO 80   
      
84    dtmax	= wsl(ndata)
!      nstor=ndata
!      GOTO  176
!****************************************************************      
      dtmax	= 10.0d0*dtmax
      nstor	= JFIX(dtmax+1.0d0/1000.0d0)
 	
! Calculate interpolation coefficients

      J=NDATA-1
      DO 100 I=2,J
	      IF(STOR(I).LE.0.0d0)GO TO 86
	      AN(I)=LOG10(STOR(I+1)/STOR(I))/LOG10(WSL(I+1)/WSL(I))
86       IF(AR(I).LE.0.0d0)GO TO 100
	      BN(I)=LOG10(AR(I+1)/AR(I))/LOG10(WSL(I+1)/WSL(I))
100   CONTINUE
      AN(1)=1.0d0
      BN(1)=1.0d0

! Note: the values of the exponents a,b for layer 1 are not used
! in any calculations as the area and volume are assumed to vary
! linearly from the bottom to the first layer

      AN(NDATA)=AN(NDATA-1)
      BN(NDATA)=BN(NDATA-1)

! DEVELOP TABLE OF VOLUMES, AREAS FOR depth INCREMENTS OF 0.1M
! kstor1 IS THE FIRST LAYER WITH a POSITIVE STORAGE
! kar1 IS THE FIRST LAYER WITH a POSITIVE AREA

      j		= 0
      kstor1	= kstoro
      kar1	= karo

      DO 175 i=1,nstor
	      d(i)=float(i)/10.0d0
 130     IF(j.eq.ndata)GOTO 140
	      j=j+1
	      IF(d(i).lt.wsl(j+1))GOTO 140
	      GOTO 130

! LINEARLY INTERPOLATE AREAS and VOLUMES FOR ALL DEPTHS BELOW THE
! FIRST CORRESPONDING NON-ZERO TABLE ENTRY

 140     IF(j.eq.1)THEN
	         lvnext=MAX0(2,kstor1)
	         lanext=MAX0(2,kar1)
	         v(i)=stor(lvnext)*d(i)/wsl(lvnext)
	         a(i)=ar(lanext)*d(i)/wsl(lanext)
	         GOTO 170
	      ENDIF

	      IF(j.lt.(kstor1))THEN  
	         v(i)=stor(kstor1)*d(i)/wsl(kstor1)
         ELSE
	         v(i)=stor(j)*(d(i)/wsl(j))**an(j)
         ENDIF

	      IF(j.lt.(kar1))THEN
	         a(i)=ar(kar1)*d(i)/wsl(kar1)
         ELSE
	         a(i)=ar(j)*(d(i)/wsl(j))**bn(j)
	      ENDIF

! RESET j FOR NEXT TABLE ENTRY i (j MUST VARY FROM 1 TO ndata FOR EACH
! PASS THROUGH THE LOOP IN ORDER TO CORRECTLY FIND THE LAYER)  
 170     j=j-1
 175  CONTINUE
 176  CONTINUE
! CALCULATE GRADIENTS BETWEEN EACH POINT FOR a and v

      DO 180 i=1,nstor-1
	      dv(i)   = (v(i+1)-v(i))
	      dadz(i) = (a(i+1)-a(i))
          area(i)=a(i)
          vol(i)=v(i)
!          WRITE(20,FMT='(I10, F15.8, F15.8, F15.8)') I,D(i),A(I),V(I)
180   CONTINUE      
	   dv(nstor)   = dv(nstor-1)
      dadz(nstor) = dadz(nstor-1)
          area(nstor)=a(nstor)
          vol(nstor)=v(nstor)
	   CLOSE(INUNIT)
!       PAUSE
	   DO i =1,nstor-1
	      IF(a(i).lt.0.0d0)THEN
            WRITE(*,*) 'Warning: error in the area table',i,a(i),i+1,a(i+1)	   
            PAUSE
         ENDIF
	   ENDDO
!
! ADD SAVE water quality and PARTICLES
!    TAHOE upgraded to include up to 100 streams
!		

!IF	   DATA FSUM,UI,UAV,TIMEI,THR,TI,AEM,FO,OLDSL,UF,HTSAVE/11*0.0d0/,MSTEP/0/,IDTEST/100*.FALSE./
      CLOSE(5)
      CLOSE(30)
      CLOSE(7)
! Initialize KO, THE INTEGER INDICATING ERROR IN water quality file
      KO=0
!
!    Initialize Inflow arrays
!    Initialize water quality and PARTICLE ARRAYS
!Bill 2012/06
      TOT_CNT=0
      PARCNT=0
      delta_vol=0.0d0

      DO 300 i=1,i3
	      TOTIN(i) = 0.0d0
	      HFLOW(i) = 0.0d0
	      DLWST(i) = 1000.0d0
	      ICNT(i)  = 0
	      NOINS(i) = 0   
	      DO 300 j=1,maxpar
	         QDOWN(i,j) = 0.0d0    
	         TDOWN(i,j) = 0.0d0    
	         SDOWN(i,j) = 0.0d0
	         DDOWN(i,j) = 0.0d0    
	         DOLD(i,j)  = 0.0d0    
	         QINS(i,j)  = 0.0d0    
	         TINS(i,j)  = 0.0d0    
	         SINS(i,j)  = 0.0d0
	         DIINS(i,j) = 0.0d0    
	         INPAR(i,j)=0   
	         DO 290 zz=1,28
	            wqdown(i,zz,j) = 0.0d0
	            WQINS(i,zz,j)  = 0.0d0
290         CONTINUE
	         DO 295 zz=1,7
	            CFDOWN(i,zz,j) = 0.0d0
	            CFINS(i,zz,j)  = 0.0d0
295         CONTINUE
!Bill 2012/06
            INTIME(J)=0.0D0
            INSINS(J)=0.0D0
            INDIINS(J)=0.0D0
            INQINS (J)=0.0D0
            INTINS (J)=0.0D0
            INPINS (J)=0.0D0
            INWQ3 (J)=0.0D0 
300   CONTINUE

! READ IN INFLOW, OUTFLOW and MET DATA FOR one EXTRA DAY, SINCE
! AVERAGING THESE DATA

      NDAYS = NDAYS + 1

!  READ STORAGE TABLE DATA
      DO 305 j=1,nstor
	      vda (j)  = v(j)	
	      dvda(j) = dv(j)     
305   CONTINUE
!
!  READ PHYSICAL DATA
!
      WRITE(*,*)'GETTING PHYSICAL DATA FROM ',flnm(3) ! ARXIU .FIZ
      OPEN(8,file=flnm(3),FORM='FORMATTED',status='OLD')
      READ(8,*) base
      READ(8,*) CRL, CRLNGTH					!Elevation of CREST OF THE RESERVOIR/LAKE and length of crest
      READ(8,*) lc
      READ(8,*) wc
      READ(8,*) latit
      READ(8,*) NUMOUT	                !OUTLET #
      READ(8,*) (OLEV(i),i=1,NUMOUT)	    !OUTLET ELEVATION (M)
      READ(8,*) (OLEN(i),i=1,NUMOUT)
      READ(8,*) (OWID(i),i=1,NUMOUT)
	   READ(8,*)CK							 !CK convective overturn
	   READ(8,*)ETA						 !ETA wind stirring
	   READ(8,*)CT							 !CT unsteady effects
	   READ(8,*)CS							 !CS shear efficiency
	   READ(8,*)AKH						 !AKH billowing
      READ(8,*) numriv		
!
!MB READING THE RIVER NAMES and segment details
!
	   DO 707 I=1, NUMRIV
         READ(8,*)
		   READ(8,*)(RIVNAM(I),numseg(I))	
		   READ(8,*)(ALPHA(I,K),K=1,numseg(I))		!1/2 angle	
		   READ(8,*)(seglngth(I,K),K=1,numseg(I))	!segment length	
		   READ(8,*)(bgnwdth(I,K),K=1,numseg(I))	!beginning width	
		   READ(8,*)(bgnele(I,K),K=1,numseg(I))	!beginning bottom elevation  
		   READ(8,*)(CDRAG(I,K),K=1,numseg(I))		!drag coefficient  
707	CONTINUE 
!gbs17Jan06 added sensor heights
		READ(8,*)
		READ(8,*)
		READ(8,*) zu, zq, ztair	 
!wef	READ(8,*)CK							 !CK convective overturn
!wef	READ(8,*)ETA						 !ETA wind stirring
!wef	READ(8,*)CT							 !CT unsteady effects
!wef	READ(8,*)CS							 !CS shear efficiency
!wef	READ(8,*)AKH						 !AKH billowing
	   CLOSE(8)
! 
! HERE ARE SOME OF THE PARAMETERS FROM THE PHYSICAL file
!
      numdif=37
	  DIFF(1) = 1.0D-7   	! DIFF(1) = 1.4D-08, 5.4D-08
      DIFF(2) = 1.25D-09     !DIFF(2) = 1.25D-09
      DIFF(3) = 0.143192D-11 
      DIFF(4) = 0.143192D-11 
      DIFF(5) = 0.143192D-11 
      DIFF(6) = 1.25D-09
      DIFF(7) = 1.25D-09
      DIFF(8) = 1.25D-09
      DIFF(9) = 1.25D-09
      DIFF(10)= 1.25D-09
      DIFF(11)= 0.143192D-11
      DIFF(12)= 1.25D-09
      DIFF(13)= 0.143192D-11
      DIFF(14)= 0.143192D-11
      DIFF(15)= 0.143192D-11
      DIFF(16)= 0.143192D-11
      DIFF(17)= 1.25D-09
      DIFF(18)= 1.25D-09    
      DIFF(19)= 0.143192D-11
      DIFF(20)= 0.143192D-11
      DIFF(21)= 0.143192D-11
      DIFF(22)= 0.143192D-11
      DIFF(23)= 1.25D-09
      DIFF(24)= 0.1D-11
      DIFF(25)= 0.1D-11
      DIFF(26)= 0.143192D-11
      DIFF(27)= 0.143192D-11
      DIFF(28)= 0.143192D-11
      DIFF(29)= 0.143192D-11
      DIFF(30)= 0.143192D-11
      DIFF(31)= 0.143192D-11
      DIFF(32)= 0.143192D-11
      DIFF(33)= 0.143192D-11
      DIFF(34)= 0.143192D-11
      DIFF(35)= 0.143192D-11
      DIFF(36)= 0.143192D-11
      DIFF(37)= 0.143192D-11
!
!    READ IN THE water quality DATA
!             	                     	 		
	   WRITE(*,*) 'GETTING water quality DATA FROM ', flnm(6)
	   WRITE(*,*) 'NUMBER OF ALGAE SPECIES,nchl =  ', nchl
	   OPEN(21,file=flnm(6),FORM='FORMATTED',status ='OLD')
	   READ(21,*)
	   READ(21,*)
	   READ(21,*)(gromax(i),i = 1,nchl)
	   READ(21,*)
	   READ(21,*)(constr(i),i = 1,nchl)
	   READ(21,*)
	   READ(21,*)(constm(i),i = 1,nchl)    
	   READ(21,*)
	   READ(21,*)(thetat(i),i = 1,nchl)
	   READ(21,*)
	   READ(21,*)(light(i),i = 1,nchl)
	   READ(21,*)
	   READ(21,*) etwat
	   READ(21,*)
	   READ(21,*)(etca(i),i = 1,nchl)     !NEW
      READ(21,*)
	   READ(21,*)(etpart(i),i = 1,nchl)   !NEW
	   READ(21,*)
	   READ(21,*)
	   READ(21,*)(ycchla(i),i = 1,nchl)   !NEW
       READ(21,*)
	   READ(21,*)(ypchla(i),i = 1,nchl)   !NEW
	   READ(21,*)
	   READ(21,*)(ynchla(i),i = 1,nchl)   !NEW 
	   READ(21,*)
	   READ(21,*)(ysichla(i),i = 1,nchl)  !NEW
	   READ(21,*)
	   READ(21,*)(phyto_setl_vel(i),i = 1,nchl) !NEW only for Chlorophyll
	   READ(21,*)
	   READ(21,*) org_setl_vel			  !NEW for PON and PON
	   READ(21,*)
	   READ(21,*)(phyto_part_fac(i),i = 1,nchl) !NEW only for Chlorophyll
	   READ(21,*)
	   READ(21,*) org_part_fac             !NEW for PON and PON
	   READ(21,*)
	   READ(21,*)
	   READ(21,*) K_SOD  !biolc
	   READ(21,*)
	   READ(21,*) constbo
	   READ(21,*)
	   READ(21,*) constnh
	   READ(21,*)
	   READ(21,*) cn1   !NEW
	   READ(21,*)
	   READ(21,*) cn2   !NEW
	   READ(21,*)
	   READ(21,*) cn3   !NEW 
	   READ(21,*)
	   READ(21,*) cn4   !NEW
	   READ(21,*)
	   READ(21,*) cn5   !NEW
	   READ(21,*)
	   READ(21,*) cp1   !NEW
	   READ(21,*)
	   READ(21,*) cp2   !NEW
	   READ(21,*)
	   READ(21,*) cp3   !NEW
	   READ(21,*)
       READ(21,*)(halfc(i), i =1,nchl) !NEW
	   READ(21,*)
	   READ(21,*)(halfn(i), i =1,nchl) !NEW
	   READ(21,*)
	   READ(21,*)(halfp(i), i = 1,nchl) !NEW
	   READ(21,*)
       READ(21,*)(halfsi(i), i =1,nchl) !NEW
	   READ(21,*)
	   READ(21,*)(hscn(i), i = 1,nchl)  !NEW
	   READ(21,*)
 	   READ(21,*) hscnit      !NEW
       READ(21,*)       
	   READ(21,*) hscdo       !NEW
	   READ(21,*)
	   READ(21,*) KH_OSOD      !hsccdo      !NEW
	   READ(21,*)
	   READ(21,*) densy(1)    !NEW
	   READ(21,*)
	   READ(21,*)
	   READ(21,*) thetanh
	   READ(21,*)
	   READ(21,*) thetabo
	   READ(21,*)
	   READ(21,*) thetabs
	   READ(21,*)
	   READ(21,*)
	   READ(21,*) sedthp
	   READ(21,*)
	   READ(21,*) sednh3
	   READ(21,*)
	   READ(21,*) sedno3
	   READ(21,*)
	   READ(21,*) sedtem
	   READ(21,*)
	   READ(21,*)
	   READ(21,*)(czoo(i),i = 1,nchl)    !NEW
	   READ(21,*)
	   READ(21,*)
	   READ(21,*)(denscf(i),i = 1,7)
	   READ(21,*)
	   READ(21,*) coag
	   CLOSE(21)
!2015/9      IF(nchl.eq.1)THEN         ! ONLY 1 GROUP
!2015/9	      DO 322 i=2,maxchl		  
!2015/9	         gromax(i) = 0.0d0	  
!2015/9	         constr(i) = 0.0d0   
!2015/9	         constm(i) = 0.0d0
!2015/9	         thetat(i) = 0.0d0
!2015/9	         light(i)  = 0.0d0
!2015/9		      hscn(i)   = 0.0d0	! NEW
!2015/9		      ycchla(i)= 0.000      !NEW
!2015/9              ypchla(i) = 0.0d0   ! NEW
!2015/9		      ynchla(i) = 0.0d0   ! NEW
!2015/9		      ysichla(i)= 0.0d0   ! NEW 
!2015/9322      CONTINUE
!2015/9      ENDIF	
	   GOTO 5178
! !
! !    Parameter Filter for SALTON SEA Project
! !
      ! DO 325 i=1,nchl
         ! IF(gromax(i).lt.0.2.or.gromax(i).gt.8.0.or.constr(i).lt.0.02            &
            ! .or.constr(i).gt.0.8.or.constm(i).lt.0.003.or.                       &
            ! constm(i).gt.0.17.or.thetat(i).lt.1.0.or.thetat(i).gt.1.14) KO = 1
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Phytoplankton Parameters'		 
            ! WRITE(*,*) 'i,gromax(i),constr(i),constm(i),thetat(i)'
            ! WRITE(*,*) i,gromax(i),constr(i),constm(i),thetat(i)
	         ! KO  = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF
 
         ! IF(light(i).lt.48.0.or.light(i).gt.194.0.or. etwat.lt.0.0.or.           &
            ! etwat.gt.3.0.or.etca(i).lt.0.0.or.etca(i).gt.2.0                    &
            ! .or.etpart(i).lt.0.00.or.etpart(i).gt.0.00)KO=1 
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Light parameters'		 
	         ! WRITE(*,*) 'light(i),etwat,etca(i),etpart(i)'             
            ! WRITE(*,*) light(i),etwat,etca(i),etpart(i)
	         ! KO  = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF
         ! IF(ypchla(i).lt.0.5.or.ypchla(i).gt.1.0.or.ynchla(i).lt.0.0           &
            ! .or.ynchla(i).gt.0.0.or.ysichla(i).lt.0.0.or.ysichla(i).gt.0.0) KO = 1

         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Phyto Nutrient Ratios'		 
            ! WRITE(*,*) 'i,ypchla(i),ynchla(i),ysichla(i)'
	         ! WRITE(*,*)  i,ypchla(i),ynchla(i),ysichla(i)
	         ! KO  = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF

         ! IF(halfn(i).lt.0.0.or.halfn(i).gt.0.0.or.halfp(i).lt.1.0.or.               &
            ! halfp(i).gt.5.0.or.hscn(i).lt.0.0.or.hscn(i).gt.0.0)  KO=1
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Phyto Nutrients'		 
            ! WRITE(*,*) 'halfn(i),hscn(i),halfp(i)'
			   ! WRITE(*,*)  halfn(i),hscn(i),halfp(i)
	         ! KO = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF

	      ! IF(czoo(i).lt.0.0d0.or.czoo(i).gt.0.0)KO=1
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*)'Error in *.WAT. Zoo feeding rater'		 
            ! WRITE(*,*) czoo(i)
	         ! KO  = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF
! 325   CONTINUE

      ! DO i = 1,nchl
	      ! IF(phyto_setl_vel(i).lt.0.0.or.phyto_setl_vel(i).gt.30.0.or.            &
            ! phyto_part_fac(i).lt.250000.or.phyto_part_fac(i).gt.300000) KO = 1    
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Settling Phyto'		 
            ! WRITE(*,*) 'i,phyto_setl_vel(i),phyto_part_fac(i)'
			   ! WRITE(*,*)  i,phyto_setl_vel(i),phyto_part_fac(i)
	         ! KO = 0
	         ! KOK = 1
		      ! PAUSE
         ! ENDIF
      ! ENDDO
	   ! IF(org_setl_vel.lt.0.05.or.org_setl_vel.gt.3.3.or.               &
         ! org_part_fac.lt.7000.or.org_part_fac.gt.9000) KO = 1 

	   ! IF(KO.eq.1)THEN
		   ! WRITE(*,*) 'Error in *.WAT. Settling Parts'		 
		   ! WRITE(*,*) 'org_setl_vel,org_part_fac'
		   ! WRITE(*,*)  org_setl_vel,org_part_fac
	
	      ! KO  = 0
	      ! KOK = 1
	      ! PAUSE
	   ! ENDIF      

      ! IF(cn1.lt.0.0.or.cn1.gt.0.0.or.cn2.lt.0.0.or.cn2.gt.0.0.or.cn3.lt.0.0.or.  &
        ! cn3.gt.0.0.or.cn4.lt.0.0.or.cn4.gt.0.0.or.cn5.lt.0.0.or.cn5.gt.0.0) KO = 1

      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. Nutrients Cycle'		 
         ! WRITE(*,*) 'cn1,cn2,cn3,cn4,cn5'
         ! WRITE(*,*)  cn1,cn2,cn3,cn4,cn5
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF

      ! IF(cp1.lt.0.6.or.cp1.gt.1.8.or.cp2.lt.0.00000545.or.cp2.gt.0.00001635      &
        ! .or.cp3.lt.0.0.or.cp3.gt.0.0)KO=1

      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. Nutrients Cycle'		 
         ! WRITE(*,*) 'cp1,cp2,cp3'
         ! WRITE(*,*)  cp1,cp2,cp3
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF

      ! IF(thetanh.lt.0.0.or.thetanh.gt.0.0.or.thetabo.lt.1.0.or.               &
         ! thetabo.gt.1.14.or.thetabs.lt.1.0.or.thetabs.gt.1.14)KO=1

      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. Nutrient const. Temperature'		 
         ! WRITE(*,*) 'thetanh,thetabo,thetabs'
         ! WRITE(*,*)  thetanh,thetabo,thetabs
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF
	 
      ! IF(hscnit.lt.0.0.or.hscnit.gt.4.0.or.hscdo.lt.0.0.or.                   &
         ! hscdo.gt.2.0.or.KH_OSOD.lt.0.0.or.KH_OSOD.gt.4.0)KO =1	
      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. DO parameters'		 
         ! WRITE(*,*) 'hscnit,hscdo,KH_OSOD'
	      ! WRITE(*,*)  hscnit,hscdo,KH_OSOD
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF
	
	   ! IF(sedthp.lt.60.0.or.sedthp.gt.2810.0.or.sednh3.lt.0.0.or.                 &
         ! sednh3.gt.0.0.or.sedtem.lt.1.0.or.sedtem.gt.1.14)KO=1

      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. Sediment parameters'		 
         ! WRITE(*,*) 'sedthp,sednh3,sedtem'
	      ! WRITE(*,*)  sedthp,sednh3,sedtem
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF
	 
	   ! IF(densy(1).lt.1001.or.densy(1).gt.1100)KO = 1
      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. BOD density'		 
         ! WRITE(*,*) 'densy(1)'
	      ! WRITE(*,*)  densy(1)
	      ! KO  = 0
	      ! KOK = 1
		   ! PAUSE
      ! ENDIF
      
      ! IF(K_SOD.lt.0.0d0.or.K_SOD.gt.15.0.or.constbo.lt.0.0d0.or.constbo.gt.0.05.or.    &
         ! constnh.lt.0.0d0.or.constnh.gt.0.5)    KO=1

      ! IF(KO.eq.1)THEN
         ! WRITE(*,*) 'Error in *.WAT. BOD constants'		 
         ! WRITE(*,*) 'K_SOD,constbo,constnh'
	      ! WRITE(*,*)  K_SOD,constbo,constnh
	      ! KO = 0
	      ! KOK = 1
	      ! PAUSE
      ! ENDIF

	   ! IF(NPARTSW)THEN
	      ! DO i = 1,7
	         ! IF(denscf(i).lt.0.0.or.denscf(i).gt.0.0) KO = 1
         ! ENDDO
	    	! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Particle density'		 
            ! WRITE(*,*) 'denscf'
	         ! WRITE(*,*)  denscf
	         ! KO  = 0
	         ! KOK = 1
	         ! PAUSE
         ! ENDIF	    
		   ! IF(coag.lt.0.0d0.or.coag.gt.0.0) KO = 1
         ! IF(KO.eq.1)THEN
            ! WRITE(*,*) 'Error in *.WAT. Particle coagulation'		 
            ! WRITE(*,*) 'coag'
	         ! WRITE(*,*)  coag
	         ! KO  = 0
	         ! KOK = 1
	         ! PAUSE
         ! ENDIF
      ! ENDIF
      ! IF(KOK.eq.1)THEN
	      ! WRITE(*,*) 'WARNING: INVALID OPTION SPECIFIED IN ',flnm(6)
	      ! STOP
      ! ENDIF

! !
! !   
!
      DO j = 1, nchl
	      DO i = 1, ns		
            IF(wqual(i,j).lt.alg_min(j)) THEN
	            WRITE(*,*)'Warning: Chla < Cha min in Layer',i,'Chla',wqual(i,j)
	            WRITE(*,*)'It will be replaced by minimum value :', alg_min(j)
	            PAUSE
	            wqual(i,j) = alg_min(j)
	         ENDIF
	      ENDDO   ! i
	   ENDDO   ! j 

!
!    Check the Initial Profiles PhytoN and PhytoP
!
      DO j = 1,nchl
	      DO i = 1,ns
!	         IF(abs(wqual(i,10+j)-(wqual(i,j)*ypchla(j))).gt.1.0D-4)THEN
!               WRITE(*,*)'Warning: IP < IP min in Layer',i
!	            STOP
!	   wqual(i,10+j) = wqual(i,j)*ypchla(j)
!	         ENDIF
!	  IF(wqual(i,16+j).ne.wqual(i,j)*ynchla(j))THEN
!         wqual(i,16+j) = wqual(i,j)*ynchla(j)
!	  ENDIF	
	      ENDDO   !i
	   ENDDO   !j 

! CONVERT light FROM MICROE/M2/S TO WATTS

!      DO 333 i=1,maxchl			!ATENCIó ES PARLA DE FOTOADAPTACIó
!	 light(i)=0.2*light(i)		!a SCHLADOW, PERO NO HO HE VIST (05/05/98)
!333   CONTINUE

5178  CONTINUE	

!
!eg  Conversion from SW to PAR units 
!eg 	0.45 W/M2 TOTAL SW = 1 W/M2 PAR
!
	   DO i = 1, nchl
!	      light(i) = light(i) * .45		
	   ENDDO

!
!    ADJUST HEIGHTS TO MEASURE FROM base (M)
!
      CRL=CRL-base	
      DO 334 i=1,NUMOUT
	      OLEV(i) = OLEV(i)-base	  
334   CONTINUE

! CALCULATE vmax, vmin, vcrl .NOTE: vmax DEFAULTS TO 2*vmin

      vmin=v(nstor)*vmin
      vmax=2.0*vmin    !DAVIS ESTABA a 10*vmin

!  CALCULATE STORAGE AT CREST LEVEL, vcrl

      X	= CRL*10.0d0
      y	= MOD(X,one)																							! Problems with Definition
      ij=x-y
      IF (ij .gt. nstor) THEN
	      y	= y+float(ij-nstor)
	      ij	= nstor
      ENDIF
      vcrl= v(ij)+y*dv(ij)
!  READ INITIAL CONDITIONS and ENHANCE RESOLUTION OF INITIAL PROFILE BY
!  INTRODUCING LAYERS BETWEEN THOSE FROM THE INPUT DATA and FILLING THE
!  VALUES WITH THE MEAN VALUE OF THE SURROUNDING LAYERS (i.E. LINEARLY
!  INTERPOLATE BETWEEN INPUT DATA VALUES). AS LONG AS SUFFICIENT SPACE
!  EXISTS, CONTINUE REDOUBLING THE NUMBER OF LAYERS.
! include water quality and PARTICLES

!  Maximum Layer Number 100. Increase Initial Number x 2
!  ************************************************************
!*************************************************************
!      IF(ns*2 .gt. 200) GOTO 360                                                                               !DAVIS SPECIFIC TAHOE

!!  Creates (2NS-1) more layers 

!335    DO 340 i=1,ns
!	  nx=2*ns-2*i+1	! ... 13, 11, 9, 7, 5, 3, 1.
!	  np=ns+1-i
!	  depth(nx)=depth(np)
!	  sal(nx)=sal(np)
!	  temp(nx)=temp(np)
!	    DO 338 zz=1,28
!	     wqual(nx,zz)=wqual(np,zz)
!338       CONTINUE
!	  DO 339 zz=1,7
!	     cf(nx,zz)=cf(np,zz)
!339       CONTINUE
!340    CONTINUE

!!  include water quality

!      DO 350 i=2,ns
!	 nx=2*i-2  ! PARELLS
!	 np=2*i-3  ! SENARS
!	 nz=2*i-1  ! SENARS
!	 temp(nx)	= (temp(np)+temp(nz))	/ 2.0d0
!	 depth(nx)	= (depth(np)+depth(nz))	/ 2.0d0
!	 sal(nx)	= (sal(np)+sal(nz))		/ 2.0d0
!
!	 DO 344 zz=1,28
!	    wqual(nx,zz)=(wqual(np,zz)+wqual(nz,zz))/2.0d0
!344      CONTINUE
!	 DO 345 zz=1,7
!	    cf(nx,zz)=(cf(np,zz)+cf(nz,zz))/2.0d0
!345      CONTINUE
!350   CONTINUE  

!      ns = 2*ns-1					! ACTUALITZA EL NOBRE DE CAPES TOTALS                     
!      IF (ns*2 .le. 200) GOTO 335 ! REPET DESDOBLAMENT FINS AL MAX PERMES
!360   CONTINUE
!************************************************************************
!******************************GBS***************************************
! Rename the initial profile starting with in_

	   DO k=1,ns
		   in_depth(k)=depth(k)
		   in_temp(k)=temp(k)
		   in_sal(k)=sal(k)
!           for water quality parameters
		   DO wq =1,28
			   in_wqual(k,wq)=wqual(k,wq)
		   ENDDO
!           for particles
		   DO ps =1,7
			   in_cf(k,ps)=cf(k,ps)
		   ENDDO
	   ENDDO !END of k look
!*******************************************************
! Make a finer distribution of initial profile by linear interpolation

	   total_sublayers=1
	   depth(1)=in_depth(1)
	   temp(1)=in_temp(1)
	   sal(1)=in_sal(1)
!           for water quality parameters
	   DO wq =1,28
	      wqual(1,wq)=in_wqual(1,wq)
	   ENDDO
!           for particles

	   DO ps =1,7
	      cf(1,ps)=in_cf(1,ps)
	   ENDDO

	   DO i=2,ns		
		   sublayers=int(in_depth(i)-in_depth(i-1))	    
! Linear interpolation from top to bottom
		DO SL=total_sublayers+1, total_sublayers+sublayers		   
		   depth(SL)=depth(SL-1)+(in_depth(i)-in_depth(i-1))/sublayers
			temp(SL)=temp(SL-1)+(in_temp(i)-in_temp(i-1))/sublayers
			sal(SL)=sal(SL-1)+(in_sal(i)-in_sal(i-1))/sublayers
!           for water quality parameters
			DO wq =1,28
	         wqual(SL,wq)=wqual(SL-1,wq)+(in_wqual(i,wq)-in_wqual(i-1, wq))/sublayers
			ENDDO !END of wq loop
!       for particles
			DO ps =1,7
	         cf(SL,ps)=cf(SL-1,ps)+(in_cf(i,ps)-in_cf(i-1, ps))/sublayers
			ENDDO !END of ps loop	       
	   ENDDO !END of SL loop
		total_sublayers=total_sublayers+sublayers
!		depth(total_sublayers)=in_depth(i)
!		temp(total_sublayers)=in_temp(i)
!		sal(total_sublayers)=in_sal(i)
!           for water quality parameters
!		DO wq =1,28
!			wqual(total_sublayers,wq)=in_wqual(i,wq)
!			ENDDO
!           for particles
!		DO ps =1,7
!			cf(total_sublayers,ps)=in_cf(i,ps)
!		ENDDO
	ENDDO !END of initial layer i loop
	ns=total_sublayers

!	DO i=1,ns
!	 WRITE (99,360) depth(i),temp(i),sal(i),wqual(i,1),wqual(i,2),
!    & wqual(i,3),wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),
!    & wqual(i,8),wqual(i,9),wqual(i,10),wqual(i,11),
!   & wqual(i,12),wqual(i,13),wqual(i,14),wqual(i,15),
!    & wqual(i,16),wqual(i,17),wqual(i,18),wqual(i,19),
!    & wqual(i,20),wqual(i,21),wqual(i,22),wqual(i,23),
!    & wqual(i,24),wqual(i,25),wqual(i,26),wqual(i,27), wqual(i,28),
!     &cf(i,1),cf(i,2),cf(i,3),cf(i,4),
!     & cf(i,5),cf(i,6),cf(i,7)
!	ENDDO
!360   FORMAT(5f7.2, 7f15.2)
!      CLOSE (99)
!	PAUSE
!****************************************************************
!************************GBS*************************************	   
!
!    WRITE INITIAL CONDITIONS file TO BINARY files
!    include water quality and PARTICLES
! 
      WRITE(*,*)'WRITING INITIAL CONDITIONS TO ',FIN1
      CALL FOPEN (' ',UIN1,FIN1,EIN1,'NEW','UNFORMATTED')           
	   WRITE (UIN1)NUMOUT,numriv,numdif,nchl,nstor,ns,ZOOSW,PARTSW,NMETSW,     &
         (temp(i),sal(i),depth(i),cf(i,1),cf(i,2),cf(i,3),cf(i,4),            &
         cf(i,5),cf(i,6),cf(i,7),wqual(i,1),wqual(i,2),                       &
         wqual(i,3),wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),              &
         wqual(i,8),wqual(i,9),wqual(i,10),wqual(i,11),                       &
         wqual(i,12),wqual(i,13),wqual(i,14),wqual(i,15),                     &
         wqual(i,16),wqual(i,17),wqual(i,18),wqual(i,19),                     &
         wqual(i,20),wqual(i,21),wqual(i,22),wqual(i,23),                     &
         wqual(i,24),wqual(i,25),wqual(i,26),wqual(i,27),wqual(i,28),i=1,ns)
      CLOSE(UIN1)
!
!     FIND FIRST RECORDS IN DAILY INFLOW, OUTFLOW and MET files
!
      CALL FOPEN (' ',URUN, FRUN,ERUN,'NEW','UNFORMATTED')
      CALL FOPEN (' ',UMET, FMET, EMET,'NEW','UNFORMATTED')
      CALL FOPEN (' ',URIV, FRIV, ERIV,'NEW','UNFORMATTED')
      CALL FOPEN (' ',UOUT, FOUT,EOUT, 'NEW','UNFORMATTED')
      OPEN(60,file=flnm(1),FORM='FORMATTED',status='OLD')   !OUTFLOW DATA FILE
      OPEN(7,file=flnm(2),FORM='FORMATTED',status='OLD')    !RIVFLOW DATA FILE
      OPEN(15,file=flnm(4),FORM='FORMATTED',status='OLD')   !METEOROLOGICAL DATA FILE
  !    jstart=1000*jyear1+jday1-1
      jstart=1000*jyear1+jday1     
      WRITE(*,*)'PROCESS DAILY DATA ... SEARCHING FOR DAY ',jstart     
      READ(7,*)  
	  READ(7,*)
      READ(7,FMT='(I10)') RIVFLO_TIMESTEP      

400   READ(7,*,END=900)ISEQB1(1),ISEQB2(1),ISEQB3(1),idbuf(1),bufinf(1)
      ISEQ1(1)=ISEQB1(1)   
      JLAST=ISEQB1(1)	
      IF(ISEQB1(1) .lt. jstart) GOTO 400
      READ(60,*)
      READ(60,FMT='(I10)') OUTFLO_TIMESTEP
410   READ(60,*,END=920)ISEQ2(1)
      JLAST=ISEQ2(1)
      IF(ISEQ2(1).lt. jstart) GOTO 410
      READ(15,*)
      READ(15,FMT='(I10)') MET_TIMESTEP
420   READ(15,*,END=940)ISEQ3(1)
      JLAST=ISEQ3(1)
      IF(ISEQ3(1) .lt. jstart) GOTO 420
!  HERE, ONCE FIRST DESIRED RECORD FROM EACH file IS FOUND
!  2015/9-------------------------------------------------------------------------
        RIV_SUBTIMESTEP=86400/RIVFLO_TIMESTEP
        OUT_SUBTIMESTEP=86400/OUTFLO_TIMESTEP
        MET_SUBTIMESTEP=86400/MET_TIMESTEP
!2015/9------------------------------------------------------------------------------
      BACKSPACE 7
      BACKSPACE 60
      BACKSPACE 15
      IF(ISEQ1(1) .eq. ISEQ2(1) .and. ISEQ2(1) .eq. ISEQ3(1).and.ISEQ1(1) .eq. jstart) GOTO 500
      JMAX = MAX0(ISEQ1(1),ISEQ2(1))
      JMAX = MAX0(JMAX,ISEQ3(1))
      IF(ISEQ1(1) .gt. jstart)WRITE(IO,430)ISEQ1(1)
      IF(ISEQ2(1) .gt. jstart)WRITE(IO,440)ISEQ2(1)
      IF(ISEQ3(1) .gt. jstart)WRITE(IO,450)ISEQ3(1)
430   FORMAT(' INFLOW DATA BEGINS ON ',I10)
440   FORMAT(' OUTFLOW DATA BEGINS ON ',I10)
450   FORMAT(' METEOROLOGICAL DATA BEGINS ON ',I10)
      WRITE(*,460)JMAX
460   FORMAT(/,' START SIMULATION ON ',I10,' (y/n)? ')
      READ(6,470)QUERY
470   FORMAT(A1)

      IF(QUERY .eq. 'N' .or. QUERY .eq. 'n') THEN
	      STOP
      ENDIF
      jstart = JMAX
      GOTO 400
!*******************************************************************
!  DAILY LOOP COMMENCES WHEN ALL INPUT DATA STARTS ON THE SAME DAY
!*******************************************************************
500   CONTINUE
!------------------------------------------------------------------
!2015/09	   OPEN(61, file='Urban_NonUrban_Flows1.inf', FORM='FORMATTED',Status='Old')		
!2015/09	   READ(61,*)
!2015/09	   DO riv = 1,64
!2015/09	      READ(61,*) RIV_ID(riv),Uflow_prct(riv),NUflow_prct(riv)
!2015/09	      Uflow_prct(riv) = Uflow_prct(riv)/100.0d0 
!2015/09	      NUflow_prct(riv)=	NUflow_prct(riv)/100.0d0 
!	print*, riv,RIV_ID(riv),Uflow_prct(riv),NUflow_prct(riv)
!	PAUSE
!2015/09       ENDDO
!2015/09    Uflow_prct(riv) = 0.01d0 
!2015/09    NUflow_prct(riv)= 0.01d0 
!------------------------------------------------------------------- 
      jcheck=jstart-1	
      DO 800 LK=1,NDAYS
	      JLAST=jcheck
	      jcheck=jcheck+1	
	      WRITE(*,*)'Processing day ',jcheck
!  CHECK FOR END OF YEAR
!  JTEST1 IS THE YEAR, JTEST2 IS THE DAY
!  JTEST3 INDICATES a LEAP YEAR WHEN IT EQUALS 0

	      JTEST1=jcheck/1000
	      JTEST2=jcheck-1000*JTEST1	
	      JTEST3=MOD(JTEST1,4)
	      IF(JTEST2 .gt. 365 .and. JTEST3 .ne. 0) THEN
	         JTEST1=JTEST1+1
	         JTEST2=1
	         jcheck=1000*JTEST1+JTEST2
	      ENDIF
	      IF(JTEST2 .gt. 366) THEN
	         JTEST1=JTEST1+1
	         JTEST2=1
	         jcheck=1000*JTEST1+JTEST2
	      ENDIF
!	WRITE(*, 9011)LK,jtest2,jtest1,reducn_part,reducn_TP,reducn_TN
! 9011 FORMAT(2i6, i10, 3f10.2) 
!******************************************************************************

!  READ DAILY INFLOW, OUTFLOW and MET DATA
! include water quality and PARTICLES. DON'T READ THE 1ST LINE
!**************SUB_DAILY DATA READING**********************************************
    DO 531 II=1,RIV_SUBTIMESTEP
	   DO 530 i=1,numriv
	      IF(NMETSW.eq.1)THEN
	         READ(7,*,END=900)ISEQB1(i),ISEQB2(i),ISEQB3(i),idbuf(i),bufinf(i),      &
             buftem(i),bufsal(i),bufwq(i,1),bufwq(i,2),        &
             bufwq(i,3),bufwq(i,8),bufwq(i,9),bufwq(i,10),     &
             bufwq(i,11),bufwq(i,12),bufwq(i,13),bufwq(i,14),  &
             bufwq(i,15),bufwq(i,16),bufwq(i,17),bufwq(i,18),  &
             bufwq(i,19),bufwq(i,20),bufwq(i,21),bufwq(i,24),  &
             bufwq(i,25),bufwq(i,26),bufwq(i,27),bufwq(i,28)   
			   bufwq(i,22) = 0.0d0
			   bufwq(i,23) = 0.0d0
			   bufcf(i,1)  = 0.0d0
			   bufcf(i,2)	= 0.0d0
			   bufcf(i,3)	= 0.0d0
			   bufcf(i,4)	= 0.0d0
			   bufcf(i,5)	= 0.0d0       
			   bufcf(i,6)	= 0.0d0
			   bufcf(i,7)	= 0.0d0

	      ELSEIF(ZOOSW.eq.1.and.PARTSW.eq.1.and.NMETSW.ne.1)THEN
	         READ(7,*,END=900)ISEQB1(i),ISEQB2(i),ISEQB3(i),idbuf(i),bufinf(i),      &
             buftem(i),bufsal(i),bufwq(i,1),bufwq(i,2),        &
             bufwq(i,3),bufwq(i,8),bufwq(i,9),bufwq(i,10),     &
             bufwq(i,11),bufwq(i,12),bufwq(i,13),bufwq(i,14),  &
             bufwq(i,15),bufwq(i,16),bufwq(i,17),bufwq(i,18),  &
             bufwq(i,19),bufwq(i,20),bufwq(i,21),bufwq(i,22),  &
             bufwq(i,23),bufcf(i,1),bufcf(i,2),bufcf(i,3),     &
             bufcf(i,4),bufcf(i,5),bufcf(i,6),bufcf(i,7)
			   bufwq(i,24)	= 0.0d0
			   bufwq(i,25)	= 0.0d0
			   bufwq(i,26)	= 0.0d0
			   bufwq(i,27)	= 0.0d0
			   bufwq(i,28)	= 0.0d0

	      ELSEIF(ZOOSW.eq.1.and.PARTSW.ne.1.and.NMETSW.ne.1)THEN
	         READ(7,*,END=900)ISEQB1(i),ISEQB2(i),ISEQB3(i),idbuf(i),bufinf(i),      &
             buftem(i),bufsal(i),bufwq(i,1),bufwq(i,2),        &
             bufwq(i,3),bufwq(i,8),bufwq(i,9),bufwq(i,10),     &
             bufwq(i,11),bufwq(i,12),bufwq(i,13),bufwq(i,14),  &
             bufwq(i,15),bufwq(i,16),bufwq(i,17),bufwq(i,18),  &
             bufwq(i,19),bufwq(i,20),bufwq(i,21),bufwq(i,22),  &
             bufwq(i,23)
			   bufwq(i,24) = 0.0d0
			   bufwq(i,25) = 0.0d0
			   bufwq(i,26) = 0.0d0
			   bufwq(i,27)	= 0.0d0
			   bufwq(i,28)	= 0.0d0
			   bufcf(i,1)	= 0.0d0
			   bufcf(i,2)	= 0.0d0
			   bufcf(i,3)	= 0.0d0
			   bufcf(i,4)	= 0.0d0
			   bufcf(i,5)	= 0.0d0
			   bufcf(i,6)	= 0.0d0
			   bufcf(i,7)	= 0.0d0
	
	      ELSEIF(ZOOSW.ne.1.and.PARTSW.eq.1.and.NMETSW.ne.1)THEN
	         READ(7,*,END=900)ISEQB1(i),ISEQB2(i),ISEQB3(i),idbuf(i),bufinf(i),        &
     		   buftem(i),bufsal(i),bufwq(i,1),bufwq(i,2),         &
     		   bufwq(i,3),bufwq(i,4),bufwq(i,5),bufwq(i,6),       &
     		   bufwq(i,7),bufwq(i,8),bufwq(i,9),bufwq(i,10),      &      ! Reading SRP bufwq(i,10) instead of THP bufwq(i,10)
     		   bufwq(i,11),bufwq(i,12),bufwq(i,13),bufwq(i,14),   &
     		   bufwq(i,15),bufwq(i,16),bufwq(i,17),bufwq(i,18),   &
     		   bufwq(i,19),bufwq(i,20),bufwq(i,21),bufwq(i,22),   &
     		   bufwq(i,23),bufwq(i,24),bufwq(i,25),bufwq(i,26),bufwq(i,27),bufwq(i,28),bufcf(i,1),bufcf(i,2),     &
     		   bufcf(i,3),bufcf(i,4),bufcf(i,5),bufcf(i,6),bufcf(i,7)

!gbs*************Change (Reduction or increase) ************************
!---------------------Stream DO-----------------------------------------
!	         DOsat = 14.5532-0.38217*buftem(i)+0.0054258*buftem(i)*                &
!                     buftem(i)-(bufsal(i)/1.80655)*(1.665d-4-5.866d-6*buftem(i)   &
!                     +9.796d-8*buftem(i)*buftem(i))

!------effect of altitude----------
!            atm_press=dexp(5.25*dlog(1.0-(CRL+BASE)/(44.3*1000.0)))
!	         part_press_wat_vap=dexp(11.8571-3840.70/(buftem(i)+273.15)-          &
!                               216961.0/(buftem(i)+273.15)**2.0)
!            DOE_theta=0.000975-1.426d-5*buftem(i)+6.436d-8*buftem(i)**2.0

!	         bufwq(i,8)=DOsat*atm_press*(1.0-part_press_wat_vap/                  &
!             atm_press)*(1.0-DOE_theta*atm_press)/((1.0-part_press_wat_vap)*(1.0-DOE_theta))
		
 	      ELSEIF(ZOOSW.ne.1.and.PARTSW.ne.1.and.NMETSW.ne.1)THEN
	         READ(7,*,END=900)ISEQB1(i),ISEQB2(i),ISEQB3(i),idbuf(i),bufinf(i),         &
               buftem(i),bufsal(i),bufwq(i,1),bufwq(i,2),         &
               bufwq(i,3),bufwq(i,8),bufwq(i,9),bufwq(i,10),      &
               bufwq(i,11),bufwq(i,12),bufwq(i,13),bufwq(i,14),   &
               bufwq(i,15),bufwq(i,16),bufwq(i,17),bufwq(i,18),   &
               bufwq(i,19),bufwq(i,20),bufwq(i,21)
			   bufwq(i,22) = 0.0d0
			   bufwq(i,23) = 0.0d0
			   bufwq(i,24) = 0.0d0
			   bufwq(i,25) = 0.0d0
			   bufwq(i,26) = 0.0d0
			   bufwq(i,27) = 0.0d0
			   bufwq(i,28) = 0.0d0
			   bufcf(i,1)  = 0.0d0
			   bufcf(i,2)	= 0.0d0
			   bufcf(i,3)	= 0.0d0
			   bufcf(i,4)	= 0.0d0
			   bufcf(i,5)	= 0.0d0         
			   bufcf(i,6)	= 0.0d0
			   bufcf(i,7)	= 0.0d0
	      ENDIF

!2015/9	      DO 527 j=5,6			!AQUESTES COLUMNES NO TENEN
!2015/9		         bufwq(i,j)	= 0.0d0       !ASSIGNADA CAP VARIABLE   
!2015/9	527      CONTINUE                      
!2015/9		      bufwq(i,4)		= bufwq(i,8)	!OXIGEN DISOLT(i,8) a 
 
!  CHECK FOR VALID OPTIONS IN INTERNAL NUTRIENT CONCENTRATIONS
! IN INFLOWS. IF NOT THEN CONVERT THE INTERNAL CONCENTRATION
! OF PHOSPHORUS and NITROGEN FROM UNITS OF MASS P or N DIVIDED BY
! MASS CHLA TO a MASS OF P or N/VOLUME water
530      CONTINUE    ! i the number of river     
	      DO 545 i=1,numriv
!           CHECK DATE ON DATA CARD
	         IF(ISEQB1(i).ne. jcheck) THEN
	            WRITE(*,535)i,jcheck
535            FORMAT(' INFLOW DATA ERROR (WRONG DATE FOR RECORD) ',I5,' DAY ',I10)
	            PAUSE
!	            STOP
	         ENDIF
!  CHECK RIVER ID
		      ID=idbuf(i)
	         IF(ID .lt. 1 .or. ID .gt. numriv) THEN
	            WRITE(*,536)jcheck,ID
536            FORMAT(' INFLOW DATA ERROR (RIVER ID OUT OF BOUNDS)',I5,' ID = ',I10)
	            PAUSE
	            STOP
		      ENDIF
		      IF(IDTEST(ID)) THEN
	            WRITE(*,537)ID,jcheck
537            FORMAT(' INFLOW DATA ERROR (MULTIPLE ENTRIES FOR',' RIVER ',I5,' ON DAY ',I10,')')
	            PAUSE
	            STOP
		      ENDIF
!  HERE IF INFLOW DATA IS OK
! include water quality and PARTICLES
		      FLOINF(ID)=bufinf(i)
		      teminf(ID)=buftem(i)
	         SALINF(ID)=bufsal(i)
			 
	         DO 538 zz=1,28
	            wqinf(ID,zz)=bufwq(i,zz)
538         CONTINUE            
	         DO 539 zz=1,7
	            cfinf(ID,zz)=bufcf(i,zz)
539         CONTINUE
		      IDTEST(ID)=.TRUE.
545      CONTINUE
!
!  RESET IDTEST FOR NEXT DAY
!
	      DO 550,i=1,numriv
	         IDTEST(i)=.FALSE.
550       CONTINUE
                    
	      JDLYY=jcheck/1000         !year
	      JDLYD=jcheck-JDLYY*1000   !day
	      DO i=1,numriv
	         JYY(I)=JDLYY
	         JDD(I)=JDLYD
	      ENDDO
!-------------- writing in binary file------------------------------------------------  
	      WRITE(URIV)(JYY(i),JDD(i),ISEQB2(i),ISEQB3(i),FLOINF(i),teminf(i),SALINF(i),     &
                wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),wqinf(i,5),        &
                wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),wqinf(i,10),       &
                wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),wqinf(i,15),   &
                wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),wqinf(i,20),   &
                wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),wqinf(i,25),   &
                wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),cfinf(i,2),     &
                cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),cfinf(i,7), i=1,numriv)     
 
531   CONTINUE
!************************************************************************************************************************************
!  GET OUTFLOW DATA
!****************************************************************************************************************************************
        DO 562 II=1,OUT_SUBTIMESTEP
            READ(60,*,END=920)ISEQ2(1),ISEQ2(2),ISEQ2(3),(DRW(i),i=1,NUMOUT)
	      IF(ISEQ2(1) .ne. jcheck)THEN
!	         WRITE(*,560)jcheck,ISEQ2(1)
!560         FORMAT(' OUTFLOW DATA ERROR ... EXPECTED DAY ',I10,' BUT READ DAY ',I10)
!	         PAUSE
!	         STOP
          ENDIF
!-------------- writing in binary file------------------------------------------------           
           WRITE(UOUT)JDLYY,JDLYD,ISEQ2(2),ISEQ2(3),(DRW(i),i=1,NUMOUT)
562     CONTINUE
  
!***********************************************************************************************************************************
!    GET and check meteorological DATA
!***************************************************************************************************************************************
        DO 563 II=1,MET_SUBTIMESTEP
        IF(humidity.eq.1) THEN
	         READ(15,*,END=940)ISEQ3(1),ISEQ3(2),ISEQ3(3),SW,SRAT,T4,RH,U6,RAIN,msecchi_dep			
	         IF(RH.gt.1.00) RH=1.0 		 
			 !RH=(RH)/(exp(2.3026*(7.5*T4/(T4+237.3)+0.7858))) !Vapour pressure to RH! ! SCT commented outand replaced with direct code from heat evaporation
	 		 SVPD = (rh)*(10.0**(9.286-(2322.38/(T4+273.15))))
	         IF(U6.le.0.30) U6=3.0d-1	
	         IF(LWIND.eq.1 .and. SRAT.gt.1) THEN
	            WRITE(*,*)'ERROR: Incoherent LW definition at *.MET input DATA'
	            WRITE(*,*)'Check LW switch in *.HDR file'
	            STOP
	         ENDIF
	         IF((RH.gt.1.0).or.(RH.lt.0.0)) THEN
	            WRITE(*,*)'Check relative humidity (DECIMAL) in *.MET' 
	            STOP
            ENDIF
	      ELSE
		      IF(LWIND.eq.1 .and. SRAT.gt.1) THEN
	            WRITE(*,*)'ERROR: Incoherent LW definition at *.MET input DATA'
	            WRITE(*,*)'Check LW switch in *.HDR file'
	            STOP
            ENDIF
            READ(15,*,END=940)ISEQ3(1),ISEQ3(2),ISEQ3(3),SW,SRAT,T4,SVPD,U6,RAIN,msecchi_dep                  
	      ENDIF
!	      SRAT=SRAT*1.10
!	      T4=T4+4.5
!	      U6=U6*0.93
	      IF((Rain.lt.0.0).or.(U6.lt.0.0).or.(SW.lt.0.0)) THEN
	         WRITE(*,*)'Check for negative value in file *.MET' 
	         STOP
         ENDIF
!        PRINT*, ISEQ3(1),ISEQ3(2),ISEQ3(3)
!        PAUSE
	      IF(ISEQ3(1) .ne. jcheck)THEN
!	         WRITE(*,570)jcheck,ISEQ3(1)
!570         FORMAT(' METEOROLOGICAL DATA ERROR ... EXPECTED DAY ',I10,' BUT READ DAY ',I10)
!	         PAUSE
!	         STOP
          ENDIF
!-------------- writing in binary file------------------------------------------------          
          IF(humidity.eq.1) THEN
              WRITE(UMET)JDLYY,JDLYD,ISEQ3(2),ISEQ3(3),SW,SRAT,T4,RH,U6,RAIN,msecchi_dep
          ELSE
              WRITE(UMET)JDLYY,JDLYD,ISEQ3(2),ISEQ3(3),SW,SRAT,T4,SVPD,U6,RAIN,msecchi_dep
          ENDIF
563  CONTINUE    
!****************************************************************************************************
          
!  WRITE DAILY RUN TIME DATA TO THE BINARY file 'D1RUN.dat'
!  include water quality and PARTICLES
!2015/09
!2015/09	      IF(humidity.eq.1) THEN
!2015/09              WRITE(URUN)JDLYY,JDLYD,(FLOINF(i),teminf(i),SALINF(i),         &
!2015/09              wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),wqinf(i,5),        &
!2015/09              wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),wqinf(i,10),       &
!2015/09              wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),wqinf(i,15),   &
!2015/09              wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),wqinf(i,20),   &
!2015/09              wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),wqinf(i,25),   &
!2015/09              wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),cfinf(i,2),     &
!2015/09              cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),cfinf(i,7),        &
!2015/09              i=1,numriv),(DRW(i),i=1,NUMOUT),SW,SRAT,T4,RH,U6,RAIN
!2015/09          ELSE
!2015/09              WRITE(URUN)JDLYY,JDLYD,(FLOINF(i),teminf(i),SALINF(i),         &
!2015/09              wqinf(i,1),wqinf(i,2),wqinf(i,3),wqinf(i,4),wqinf(i,5),        &
!2015/09              wqinf(i,6),wqinf(i,7),wqinf(i,8),wqinf(i,9),wqinf(i,10),       &
!2015/09              wqinf(i,11),wqinf(i,12),wqinf(i,13),wqinf(i,14),wqinf(i,15),   &
!2015/09              wqinf(i,16),wqinf(i,17),wqinf(i,18),wqinf(i,19),wqinf(i,20),   &
!2015/09              wqinf(i,21),wqinf(i,22),wqinf(i,23),wqinf(i,24),wqinf(i,25),   &
!2015/09              wqinf(i,26),wqinf(i,27),wqinf(i,28),cfinf(i,1),cfinf(i,2),     &
!2015/09              cfinf(i,3),cfinf(i,4),cfinf(i,5),cfinf(i,6),cfinf(i,7),        &
!2015/09              i=1,numriv),(DRW(i),i=1,NUMOUT),SW,SRAT,T4,SVPD,U6,RAIN
!2015/09          ENDIF 
!2015/09
800   CONTINUE ! total number of days	      						        	
!*****************************************************************************
!*********************END OF DAILY LOOP***************************************
!*****************************************************************************
!
!    CHANGED DIM maxpar
!
!2015/9      DO 802 i=1,numriv
!2015/9	      DO 802 j=1,maxpar
!2015/9	         DO 801 KK=5,7
!2015/9	            wqdown(i,KK,j)=0.0d0
!2015/9	            WQINS(i,KK,j)=0.0d0
!2015/9 801         CONTINUE
!2015/9	      wqdown(i,8,j)=wqdown(i,4,j)
!2015/9	      WQINS(i,8,j)=WQINS(i,4,j)
!2015/9 802   CONTINUE

      CLOSE(10)
      CLOSE(60)
      CLOSE(7)
      CLOSE(15)
      CLOSE(22)

! REVERT BACK TO THE ACTUAL NUMBER OF DAYS TO SIMULATE

      NDAYS=NDAYS-1

!  WRITE PARAMETER file TO 'D1TOUT.dat'
! include water quality and PARTICLES

!810   OPEN(12,file='D1TOUT.dat',FORM='UNFORMATTED',status='UNKNOWN')
    
      CALL FOPEN (' ',UPAR,FPAR,EPAR,'NEW','UNFORMATTED')  
      WRITE(*,*)'WRITING PARAMETER DATA TO ',FPAR
      WRITE(UPAR)NDAYS,LWIND,(DIFF(I),I=1,NUMDIF),                            &
         CRL,CRLNGTH,LC,WC,LATIT,(OLEV(I),OLEN(I),OWID(I),                    &
         I=1,NUMOUT),(numseg(I),I=1,numriv),                                  &
         ((ALPHA(I,k),seglngth(I,k),bgnwdth(I,k),bgnele(I,k),CDRAG(I,k),      &
         I=1,numriv),k=1,10),CK,ETA,                                          &
         CT,CS,AKH,(V(I),DV(I),A(I),DVDA(I),VDA(I),DADZ(I),I=1,NSTOR),VMAX,   &
         VMIN,VCRL,DMAX,DMIN,(((WQDOWN(I,ZZ,J),WQINS(I,ZZ,J),i=1,numriv),     &
         ZZ=1,28),j=1,5),(((CFDOWN(I,ZZ,J),CFINS(I,ZZ,J),i=1,numriv),         &
         ZZ=1,7),j=1,5),((QDOWN(I,J),QINS(I,J),SDOWN(I,J),SINS(I,J),          &
         TDOWN(I,J),TINS(I,J),DDOWN(I,J),DOLD(I,J),DIINS(I,J),                &
         INPAR(I,J),J=1,MAXPAR),HFLOW(I),TOTIN(I),DLWST(I),ICNT(I),           &
         NOINS(I),I=1,NUMRIV),MSTEP,FSUM,UI,UAV,TIMEI,THR,TI,AEM,FO,          &
         OLDSL,UF,HTSAVE,JSTART,RESNAM, zu, zq,ztair
!gbs 17Jan06 added sensor heights zu, zq, ztair
      CLOSE(UPAR)

      CALL FOPEN (' ',UWQD,FWQD,EWQD,'NEW','UNFORMATTED')           

      WRITE(*,*) 'WRITING PARAMETER DATA TO ',FWQD
      WRITE(UWQD)(gromax(i),constr(i),constm(i),thetat(i),light(i),            &
         ycchla(i),ynchla(i),ypchla(i),ysichla(i),halfc(i),halfn(i),halfp(i),  &
         halfsi(i),hscn(i),czoo(i),etca(i),i=1,nchl),densy(1),etwat,           &
         K_SOD,KH_OSOD,constbo,constnh,cn1,cn2,cn3,cn4,cn5,cp1,cp2,cp3,        &
         thetabo,thetanh,thetabs,sedthp,sednh3,sedno3,sedtem,                  &
         coag,(denscf(i),etpart(i),i=1,7),(phyto_setl_vel(i),phyto_part_fac(i),i=1,nchl), &
         org_setl_vel,org_part_fac    !!    DAVIS: ADDED  
      CLOSE(UWQD)

	   CLOSE(URUN)
       CLOSE(UMET)
       CLOSE(URIV)
       CLOSE(UOUT)
	   CLOSE(UIN1)
	   CLOSE(UPAR)
	   CLOSE(UWQD)
	   CLOSE(ULIS)
	   CLOSE(UFLD)

	   WRITE(*,*)'THE RELEVANT DATA files ARE AS FOLLOWS : '
	   WRITE(*,*)FIN1
	   WRITE(*,*)FRUN
       WRITE(*,*)FMET
       WRITE(*,*)FRIV
       WRITE(*,*)FOUT
	   WRITE(*,*)FPAR
	   WRITE(*,*)FWQD
      RETURN
!
!     HERE IF INPUT DATA ENDS BEFORE EXPECTED
!
900   WRITE(*,910)JLAST
910   FORMAT(' INFLOW DATA file ENDS ON ',I10,/)
      GOTO 960

920   WRITE(*,930)JLAST
930   FORMAT(' OUTFLOW DATA file ENDS ON ',I10,/)
      GOTO 960

940   WRITE(*,950)JLAST
950   FORMAT(' METEOROLOGICAL DATA file ENDS ON ',I10,/)

!  CALCULATE HOW MANY DAYS OF DATA WE HAVE

960   JLASTY	= JLAST/1000
      JLASTD	= JLAST-1000*JLASTY
      JSTRTY	= jstart/1000
      JSTRTD	= jstart-1000*JSTRTY
      NLEAPL	= JLASTY/4
      JLASTL	= MOD(JLASTY,4)
      IF(JLASTL .gt. 0)NLEAPL=NLEAPL+1
      JD0LST	= 366*NLEAPL+365*(JLASTY-1-NLEAPL)+JLASTD
      NLEAPL	= JSTRTY/4
      JSTRTL	= MOD(JSTRTY,4)
      IF(JSTRTL .gt. 0)NLEAPL=NLEAPL+1
      JD0STR	= 366*NLEAPL+365*(JSTRTY-1-NLEAPL)+JSTRTD
      NEWDAZ	= JD0LST-JD0STR+1
      NEWDAZ	= MAX0(NEWDAZ,0)  
      WRITE(*,970)NEWDAZ
970   FORMAT(' THERE IS ONLY ENOUGH DATA TO SIMULATE ',I10,' DAYS')
      WRITE(*,'(a)') ' DO YOU WISH TO PROCEED (Y/N)? : '
      READ(*,470)QUERY

      IF(QUERY .eq. 'N' .or. QUERY .eq. 'n')THEN
	      STOP
      ENDIF
      NDAYS=NEWDAZ

      GOTO 9999

8000  WRITE (MESS,3800)                                    
3800  FORMAT ('READ ERROR IN file')
      CALL ERRORS (MESS)    
!
! RE-Initialize INPUT file NAMES
! 
      DO 8010 i=1,6                                                             
	      flnm(i) = ' '                                                        
8010  CONTINUE                                                                  
      CLOSE(ULIS)
      FLIS = ' '                                                                
      GOTO 3         
!
9999  CONTINUE
!	WRITE(*,*)'INITIAL DATA file IS ',FINI

	   CLOSE(URUN)
       CLOSE(UMET)
       CLOSE(URIV)
       CLOSE(UOUT)
	   CLOSE(UIN1)
	   CLOSE(UPAR)
	   CLOSE(UWQD)
	   CLOSE(ULIS)
	   CLOSE(UFLD)
      RETURN 
      END SUBROUTINE MOD1	   

!***********************************************************************
      SUBROUTINE ERRORS (MESSGE) 
!                                                                  
!***********************************************************************
!                                                                                                                                        
!     SUBROUTINE ERRORS                                                                                          
!                                                                                                                                        
!     DESCRIPTION:      WRITES ERROR MESSAGES TO THE SCREEN.                     
!                                                                                                                                        
!      THE LOGIC IS:                                                   
!                                                                      
!         STEP 1.     OUTPUT MESSAGE TO THE SCREEN.                    
!                     PROMPT THE USER TO PRESS ENTER TO CONTINUE.      
!                     return.                                          
!                                                                      
!     VARIABLES:                                                       
!                                                                      
!         CH       (CH*1)    READS CARRIAGE return TO CONTINUE PROGRAM 
!                                                                                                                                        
!***********************************************************************
!
      CHARACTER*80 MESSGE

      WRITE (*,3000) MESSGE                                             
3000  FORMAT(/, 1X, 78A)

      WRITE (*,3020)                                                    
3020  FORMAT (/, 1X, 'PRESS <ENTER> TO CONTINUE')                       
      READ (*,2010)                                                  
2010  FORMAT (A1)       
!
! EXIT PROGRAM GRACEFULLY
!
      CALL QUIT
      END SUBROUTINE ERRORS                                                       
 !************************************************************************
      SUBROUTINE QUIT                                                   
!                                                                       
!***********************************************************************
!                                                                  
!     SUBROUTINE QUIT                                              
!                                                                  
!     DESCRIPTION:    THIS SUBROUTINE OUTPUTS a MESSAGE AS THE USER
!                     EXITS THE PROGRAM.                           
!                                                                  
!***********************************************************************
!                                                                       
      WRITE (*,2)                                                    
2     FORMAT(//15X,'________________________________________________'   &
             //15X,'                   BYE!'                            &
             //15X,'________________________________________________' )

      STOP
      RETURN    
      END SUBROUTINE QUIT
!---------------------------------------------------------------------------
