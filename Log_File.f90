!***********************************************************************
      SUBROUTINE Log_File (partname, fpar,fin1)
!***********************************************************************
! 
!     Subroutine to WRITE the LOG file 
!     Set up for 2 algae groups and adapted to LAKE TAHOE project
!     written by Joaquim Losada and Geoff Schladow. UC-DAVIS 2000.
! 
!***********************************************************************      
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
		
      INTEGER*4    ulog,JSTART
		INTEGER*4 LWIND
		INTEGER*4    i,j,K,kk
		CHARACTER*20 partname
		CHARACTER*20 logname
		CHARACTER*20 fpar,fin1
		CHARACTER*80 qp(41)  
      CHARACTER*80 qpx(41)
		REAL*8       xinf(maxinf),xout(maxout)
		REAL*8 latit
		PARAMETER(ulog = 70)    
! 
!     Definition of the WRITE staments
! 
		qp(1)='Maximum algal growth rate/day  '
		qp(2)='Maximum algal respiratory rate/day '
		qp(3)='Maximum algal mortality rate/day '
		qp(4)='temperature multiplier for growth/respiration/death '
		qp(5)='Light saturation for phytoplankton/microE/m2/s '
		qp(6)='Light attenuation, background/m-1 '
		qp(7)='Specific extinction coefficient Chla/m2/Chla(mg/L) '
		qp(8)='Phosphorus to chlorophyll mass ratio  '
		qp(9)='Nitrogen to chlorophyll mass ratio  '
		qp(10)='Silica to chlorophyll mass ratio  '
		qp(11)='Setting velocity for phytoplankton (m/day) '
		qp(12)='Setting velocity for detritus POP & PON (m/day) '
		qp(13)='Phytoplankton transfer function (#/Cha(microg)'
		qp(14)='POP and PON transfer function (#/(microg) '
      qp(15)='Oxygen demand of sub-euphotic sediments/g/m2/day '
      qp(16)='Decomposition rate of BOD/day '
		qp(17)='Half saturation constnat eff DO on denitrif.(g/m3) '
      qp(18)='PON ==> NH4 rate/day  '
      qp(19)='PON ==> doN rate/day  '
      qp(20)='DON ==> NH4 rate/day  '
      qp(21)='NH4 ==> NO3 rate/day  '
      qp(22)='NO3 ==> N2 rate/day '
      qp(23)='POP ==> THP rate/day  '
      qp(24)='POP ==> RP rate/day  ' 
      qp(25)='RP ==> THP rate/day ' 
      qp(26)='Half sat constant for N nutrient limitation  '
      qp(27)='Half sat constant for P nutrient limitation  '
		qp(28)='Half sat constant for Ammonia preferential uptake factor '
      qp(29)='Half sat cons for limit of reac by DO for nitrif/microg/l'
		qp(30)='Half sat constant for limit of reac by DO for BOD /microg/l  '
		qp(31)='Half sat constant for limit of reac by DO for sed proc. '
      qp(32)='Density of BOD for settling kg/m3 '
      qp(33)='Nitrification temperature multiplier '
      qp(34)='Organic decomposition temperature multiplier '
		qp(35)='Biological and chemical sediment oxygen demand '
		qp(36)='Release rate of phosphorus THP from sediments/microg/m2/day '
		qp(37)='Release rate of nitrogen NH4 from sediments/microg/m2/day '
		qp(38)='Temperature multiplier for sediment nutrient release '
		qp(39)='Zoo feeding rate on Cholophyll '
		qp(40)='Density of 7 particle size groups/Kg/m3 '
		qp(41)='Coagulation efficiency factor for particles '
! 
!     Range of Values
! 
		qpx(1)='(1.0-2.5 range, 1.75 mean):'
		qpx(2)='(0.05-0.15 range, 0.10 mean):'
		qpx(3)='(0.003-0.17 range,  0.0865 mean):'
		qpx(4)='(1.0-1.1 range, 1.08 mean):'
		qpx(5)='(100.0-500.0 range, 300.0 mean)'
		qpx(6)='(0.05-0.65 range, 0.45 mean) '
		qpx(7)='0.08(Geoff)  Scalar'
		qpx(8)='(0.5-1.0 range, 0.75 mean)'
		qpx(9)='(7.0-15.0 range, 11 mean) '
		qpx(10)='(40.0-50.0 range, 45.0 mean) '
		qpx(11)='(0.5-3.0 range for Tahoe, 3.26 mean)('
		qpx(12)='(0.05-0.5 range for Tahoe, 0.275 mean)'
		qpx(13)='transfer function (#/(microg) '
		qpx(14)='transfer function (#/(microg) '
		qpx(15)='(0.02-15.0 range, 5.0 mean)'
		qpx(16)='(0.005-0.05 range, 0.01 mean)'
      qpx(17)='(0.005-0.2 range, 0.01 mean)'
		qpx(18)='( range,  mean) (NEW)'
		qpx(19)='( range,  mean) (NEW)'
		qpx(20)='( range,  mean) (NEW)'
		qpx(21)='( range,  mean) (NEW)'
		qpx(22)='(0.01-0.100 range, 0.05 mean) '
		qpx(23)='(0.001-0.1 range, 0.05 mean) '
		qpx(24)='(range, mean) '
		qpx(25)='(0.01-0.100 range, 0.05 mean) '
		qpx(26)='(10-400)(NEW)vector'
		qpx(27)='(4-80)(NEW)vector'
		qpx(28)='(20.0-120.0 range, 25.0 mean)'
		qpx(29)='/microg/l  (NEW 1 value. Keep fixed)'
		qpx(30)='microg/l  (NEW 1 value. Keep fixed'
		qpx(31)='microg/l  (NEW 1 value. Keep fixed'
		qpx(32)='(1001-1100 range, 1.040 mean)('
		qpx(33)='(1.02-1.14 range, 1.08 mean)'
		qpx(34)='(1.02-1.14 range, 1.08 mean)'
		qpx(35)='(NEW. 1 Values. Keep fixed)'
		qpx(36)='(0.0-5000.0 range, 2500.0 mean)'
		qpx(37)='(0.0-10000.0 range, 5000.0 mean)'
		qpx(38)='(1.02-1.14 range, 1.08 mean):'
		qpx(39)='(microgChl/l)'
		qpx(40)='(1005-1200 range, 1125 mean)'
		qpx(41)='(-0.05-0.05 range, 0.005 mean) '
! 
!     Set the LOG file name
! 	
		logname = 'LogFile.dat'
! 
!     Open the LOG file 
! 
      OPEN(ulog,FILE=LOGNAME,FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(*,*)'WRITING LOG FILE : ',logname
! 
!     Get and WRITE the relevant parts of the PARMETER FILE to LOG FILE
! 
      OPEN(4,FILE=FPAR,FORM='UNFORMATTED',STATUS='OLD')

		READ(4)NDAYS,LWIND,(DIFF(I),I=1,NUMDIF),CRL,CRLNGTH,LC,WC,LATIT,		&
				(OLEV(I),OLEN(I),OWID(I),I=1,NUMOUT),									&
				((ALPHA(I,k),seglngth(I,k),bgnwdth(I,k),								&
				bgnele(I,k),CDRAG(I,k),I=1,NUMINF),k=1,maxseg),						&
				CK,ETA,CT, CS,AKH,(VV(I),DVV(I),A(I),DVVDA(I),VVDASH(I),			&
				DADZ(I),I=1,NUMSTO),VMAX,VMIN,VCRL,DMAX,DMIN,						&
				(((WqDOWN(I,KK,J),WqINS(I,KK,J),I=1,NUMINF),KK=1,28),J=1,5),	&
				(((CFINS(I,KK,J),CFDOWN(I,KK,J),I=1,NUMINF),KK=1,7),J=1,5),		&
				((qDOWN(I,J),qINS(I,J),SDOWN(I,J),SINS(I,J),							&
				TDOWN(I,J),TINS(I,J),DDOWN(I,J),DOLD(I,J),DIINS(I,J),				&
				INPAR(I,J),J=1,MAXPAR),HFLOW(I),TOTIN(I),DLWST(I),ICNT(I),		&
				NOINS(I),I=1,NUMINF),MSTEP,FSUM,UI,UAV,TIMEI,THR,TI,AEM,FO,		&
				OLDSL,UF,HTSAVE,JSTART,RESNAM, zu, zq, ztair

!gbs17Jan06 added meteorological sensors height zu, zq and ztair

      CLOSE(4)

		WRITE(ulog,*)
		WRITE(ulog,*)'INITIAL PARAMETER DATA'
		WRITE(ulog,*)
		WRITE(ulog,5543)RESNAM,JSTART,dmax,DMIN
		WRITE(ulog,5538)NDAYS,LWIND,CRL,LC,WC,LATIT
		WRITE(ulog,5539)

		DO i = 1,numout
			WRITE(ulog,5540)OLEV(i),OLEN(i),OWID(i)
		END DO

		WRITE(ulog,5541)
		DO i = 1,numinf
!			WRITE(ulog,5542)alpha(i),PHI(i),cdrag(i)
		ENDDO
! 
!     Get and WRITE INITIAL DATA FILE to LOG FILE
! 
      OPEN (7,FILE=FIN1,FORM='UNformatTED',STATUS='OLD')
		READ (7)numout,numinf,numdif,nchl,NUMSTO,ns,ZOOSW,PARTSW,NMETSW,		&
			(temp(i),sal(i),depth(i),cf(i,1),cf(i,2),cf(i,3),cf(i,4),			&
			cf(i,5),cf(i,6),cf(i,7),wqual(i,1),wqual(i,2),wqual(i,3),			&
			wqual(i,4),wqual(i,5),wqual(i,6),wqual(i,7),								&
			wqual(i,8),wqual(i,9),wqual(i,10),wqual(i,11),							&
			wqual(i,12),wqual(i,13),wqual(i,14),wqual(i,15),						&
			wqual(i,16),wqual(i,17),wqual(i,18),wqual(i,19),						&
			wqual(i,20),wqual(i,21),wqual(i,22),wqual(i,23),						&
			wqual(i,24),wqual(i,25),wqual(i,26),wqual(i,27),						&
			wqual(i,28),i=1,ns)
      CLOSE(7)

		WRITE(ulog,*)
		WRITE(ulog,*)'INITIAL PROFILE DATA'
      WRITE(ulog,5531)numout,numinf,numdif,nchl,NUMSTO,ns
      WRITE(ulog,5532)
	
		DO i = 1,ns

			WRITE(ulog,5533)depth(i),temp(i),sal(i),wqual(i,1),					&
				wqual(i,2),wqual(i,3),wqual(i,8),wqual(i,9),wqual(i,10),			&					
				wqual(i,11),wqual(i,12),wqual(i,13),wqual(i,14),					&
				wqual(i,15),wqual(i,16),wqual(i,17),wqual(i,18),					&
				wqual(i,19),wqual(i,20),wqual(i,21),									&
				cf(i,1),cf(i,2),cf(i,3),cf(i,4),cf(i,5),cf(i,6),cf(i,7)
		ENDDO

		WRITE(ulog,*)

5531  FORMAT(/,'NUMBER OF OUTLETS',5X,I5,/,										&
               'NUMBER OF INFLOWING STREAMS',5X,I5,/,							&
               'NUMBER OF DIFFUSING SUBSTANCES',5X,I5,/,						&
               'NUMBER OF ALGAL GROUPS',5X,I5,/,								&
               'NUMBER OF DATA POINTS',5X,I5,/,									&
               'NUMBER OF PROFILE LAYERS',5X,I5)		

5532  FORMAT(/,5x,'DEPTH',7x,'TEMP',8X,'sal',12x,'CHLA1',8x,				&
        'CHLA2',8X,'CHLA3',8X 'DO',11x,'BOD',10x,'THP',10x,'PIN1',9x,	&
        'PIN2',9x,'POP',9X,'RP',9x,'NO3',11x,'NH4',11X,'NIN1',9x,			&
        'NIN2',9x,'PON',9x,'DON',12X,'Si',14X,									&
        'PART1',8X,'PART2',10X,'PART3',10X,'PART4',10X,'PART5',10X,		&
        'PART6',10X,'PART7')

		WRITE(ulog,*)
5533  FORMAT(20(F10.3,3X),7(F15.1))

5538  FORMAT(/,'DAYS SIMULATED',3X,I5,/,											&	
       'TYPE OF LONGWAVE RADIATION',3X,I5,/,										&
       'CREST ELEVATION',3X,F7.2,/,													&
       'LENGTH AT CREST',3X,F12.2,/,												&
       'WIDTH AT CREST',3X,F12.2,/,													&
       'LATITUDE',3X,F7.2)
		
5539  FORMAT(/,'OUTLET ELEVATION',3X,'LENGTH AT OUTLET',3X,'WIDTH AT OUTLET')
5540  FORMAT(3(F10.2,10X))
5541  FORMAT(/,'HALF ANGLE',3X,'STREAM BED SLOPE',3X,'DRAG COEFFICIENT')
5542  FORMAT(3(F10.3,8X))

5543  FORMAT(/,'RESERVOIR NAME',3X,A20,/,											&
       'START DATE',3X,I7,/,															&
       'MAX LAYER THICKNESS',3X,F7.2,/,											&
       'MIN LAYER THICKNESS',3X,F7.2)
! 
!    	WRITE water quality parameters
! 
	   WRITE(ulog,*)  
	   WRITE(ulog,*)'WATER QUALITY DATA'
	   WRITE(ulog,*)qp(1)//qpx(1)
	   WRITE(ulog,*)(gromax(i),i = 1,nchl)
	   WRITE(ulog,*)qp(2)//qpx(2)
	   WRITE(ulog,*)(constr(i),i = 1,nchl)
	   WRITE(ulog,*)qp(3)//qpx(3)
	   WRITE(ulog,*)(constm(i),i = 1,nchl)    
	   WRITE(ulog,*)qp(4)//qpx(4)
	   WRITE(ulog,*)(thetat(i),i = 1,nchl)
	   WRITE(ulog,*)qp(5)//qpx(5)
	   WRITE(ulog,*)(light(i),i = 1,nchl)
	   WRITE(ulog,*)qp(6)//qpx(6)
	   WRITE(ulog,*)etwat
	   WRITE(ulog,*)qp(7)//qpx(7)
	   WRITE(ulog,*)(etca(i),i = 1,nchl)     !NEW
	   WRITE(ulog,*)qp(8)//qpx(8)
	   WRITE(ulog,*)(ypchla(i),i = 1,nchl)   !NEW
	   WRITE(ulog,*)qp(9)//qpx(9)
	   WRITE(ulog,*)(ynchla(i),i = 1,nchl)   !NEW 
	   WRITE(ulog,*)qp(10)//qpx(10)
	   WRITE(ulog,*)(ysichla(i),i = 1,nchl)  !NEW
	   WRITE(ulog,*)qp(11)//qpx(11)
	   WRITE(ulog,*)(phyto_setl_vel(i),i = 1,nchl) !NEW
	   WRITE(ulog,*)qp(12)//qpx(12)
	   WRITE(ulog,*)org_setl_vel     			 !NEW 
	   WRITE(ulog,*)qp(13)//qpx(13)
	   WRITE(ulog,*)(phyto_part_fac(i),i = 1,nchl) !NEW 
	   WRITE(ulog,*)qp(14)//qpx(14)
	   WRITE(ulog,*)org_part_fac             !NEW
	   WRITE(ulog,*)qp(15)//qpx(15)
	   WRITE(ulog,*)K_SOD
	   WRITE(ulog,*)qp(16)//qpx(16)
	   WRITE(ulog,*)constbo
	   WRITE(ulog,*)qp(17)//qpx(17)
	   WRITE(ulog,*)constnh
	   WRITE(ulog,*)qp(18)//qpx(18)
	   WRITE(ulog,*)cn1                       !NEW
	   WRITE(ulog,*)qp(19)//qpx(19)
	   WRITE(ulog,*)cn2                       !NEW
	   WRITE(ulog,*)qp(20)//qpx(20)
	   WRITE(ulog,*)cn3                       !NEW 
	   WRITE(ulog,*)qp(21)//qpx(21)
	   WRITE(ulog,*)cn4                       !NEW
	   WRITE(ulog,*)qp(22)//qpx(22)
	   WRITE(ulog,*)cn5                       !NEW
	   WRITE(ulog,*)qp(23)//qpx(23)
	   WRITE(ulog,*)cp1                       !NEW
	   WRITE(ulog,*)qp(24)//qpx(24)
	   WRITE(ulog,*)cp2                       !NEW
	   WRITE(ulog,*)qp(25)//qpx(25)
	   WRITE(ulog,*)cp3                       !NEW
	   WRITE(ulog,*)qp(26)//qpx(26)
	   WRITE(ulog,*)(halfn(i), i =1,nchl)  !NEW
	   WRITE(ulog,*)qp(27)//qpx(27)
	   WRITE(ulog,*)(halfp(i), i = 1,nchl) !NEW
	   WRITE(ulog,*)qp(28)//qpx(28)
	   WRITE(ulog,*)(hscn(i), i = 1,nchl)  !NEW
	   WRITE(ulog,*)qp(29)//qpx(29)
 	   WRITE(ulog,*)hscnit                 !NEW
      WRITE(ulog,*)qp(30)//qpx(30)       
	   WRITE(ulog,*)hscdo                  !NEW
	   WRITE(ulog,*)qp(31)//qpx(31)
	   WRITE(ulog,*)KH_OSOD                 !NEW
	   WRITE(ulog,*)qp(32)//qpx(32)
	   WRITE(ulog,*)densy(1)               !NEW
	   WRITE(ulog,*)qp(33)//qpx(33)
	   WRITE(ulog,*)thetanh
	   WRITE(ulog,*)qp(34)//qpx(34)
	   WRITE(ulog,*)thetabo
	   WRITE(ulog,*)qp(35)//qpx(35)
	   WRITE(ulog,*) thetabs
	   WRITE(ulog,*)qp(36)//qpx(36)
	   WRITE(ulog,*)sedthp
	   WRITE(ulog,*)qp(37)//qpx(37)
	   WRITE(ulog,*)sednh3
	   WRITE(ulog,*)qp(38)//qpx(38)
	   WRITE(ulog,*)sedtem
	   WRITE(ulog,*)qp(39)//qpx(39)
	   WRITE(ulog,*)(czoo(i),i = 1,nchl)    !NEW
	   WRITE(ulog,*)qp(40)//qpx(40)
	   WRITE(ulog,*)(denscf(i),i = 1,7)
	   WRITE(ulog,*)qp(41)//qpx(41)
	   WRITE(ulog,*)coag
! 
!     WRITE multiplicative factors for wind, inflow and outflow
! 
		WRITE(ulog,*)
		WRITE(ulog,5600)'WIND FACTOR : ',xwind

		DO i = 1,numinf
			WRITE(ulog,5601)'INFLOW FACTOR ',i,': ',xinf(i)
		ENDDO

		DO i = 1,numout
			WRITE(ulog,5602)'OUTFLOW FACTOR ',i,': ',xout(i)
		ENDDO
5600  FORMAT(A20,F7.2)
5601  FORMAT(A20,I5,A5,7F7.2)
5602  FORMAT(A20,I5,A5,7F7.2)

		CLOSE(ulog)
		RETURN
		END SUBROUTINE Log_File