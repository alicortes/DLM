!************************************************************************
		SUBROUTINE ZOOPLANKTON
!************************************************************************
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
		
		REAL(8), PARAMETER :: pi= 3.14159265359d0
		REAL(8), PARAMETER :: alfa = 3.0d0, beta= 4.0d0
		INTEGER*4 i, j, k, nzoo
		PARAMETER (nzoo = 1)
		REAL(8) zoop(maxns,	nzoo)
		REAL(8) zoogrowth_rate(maxns,nzoo)
		REAL(8) grazing_algae (maxns, nzoo)
		REAL(8) pref_algae (nzoo, nchl) 
		REAL(8) grazing_detrius (maxns, nzoo)
		REAL(8) pref_detrius (nzoo)
		REAL(8) grazing_max (nzoo)
		REAL(8)	half_sat_const (nzoo)
		REAL(8)	pref_zoo2algae (nzoo, nchl), pref_zoo2det (nzoo)
		REAL(8) temp_opt, temp_ref, kmax, temp_const
		REAL(8)	total_feed (maxns, nzoo)
		REAL(8)	phyto(maxns,nchl)
		REAL(8)	POC (maxns)
		REAL(8) pred1,pred2,sum_pref
		REAL(8) zoop_growth(maxns, nzoo), zoop_death(maxns,nzoo)
		REAL(8) mortality (nzoo), predation (maxns, nzoo)
		PARAMETER (pred1 = 0.15d0, pred2 = 40.0d0)
	
!gbs-------Zooplankton feeding------------------------------------
		DO i = 1,ns
			POC(i)=0.0d0
			DO k = 1,nchl
!gbs-------------Respiration and mortality------------------------
				POC(i) = POC(i)+death(i,k)
				PHYTO(i,k)= wqual(i,k)
			ENDDO	 !j

			DO j=1,nzoo
				zoogrowth_rate(i,j) = 1.0d0

				IF (zoop(i,j).le.0.00001d0) THEN
					zoop(i,j)=0.0001d0
				ELSE
					zoop(i,j)=zoop(i,j)
				ENDIF
			ENDDO
		ENDDO !i
!gbs---------------------------------------------------------------
!		IF(iclock.ge.sunrise.and.iclock+nosecs.le.sunset) THEN ! Day conditions 
!			GOTO 1
!		ELSE
!			GOTO 99
!		ENDIF
    1 CONTINUE
!gbs-------------------------------------------------------------------
		DO j = 1, nzoo
			grazing_max (j) = 0.45d0      !per day
			half_sat_const (j) = 120.0d0
			pref_zoo2det (j)=0.01
			mortality(j) = 0.05
			DO k = 1, nchl
				pref_zoo2algae (j, k)=0.25
			ENDDO
		ENDDO

!		zoo = alfa*sin(pi*jday/180 + beta)      
!     zoop = 0.01D0 
!gbs	Growth of zooplankton--------------------------------------------- 
		DO i = 1,ns
!-------------Temperature multiplier---------------------------  
			temp_opt=5.6d0   !	   

         IF(temp(i).le.temp_opt) THEN	   !!  temp(i)
				kmax =0.1*DEXP(-0.004*((temp(i)-temp_opt)**2.0))
			ELSEIF(temp(i).gt.temp_opt) THEN
				kmax =0.1*DEXP(-0.004*((temp_opt-temp(i))**2.0))			     !28 33
			ENDIF
			DO j = 1,nzoo  
!--------------preference of zooplnkton to eat to diff phytoplankton---
            sum_pref=0.0d0
				total_feed(i,j)=0.0d0
            DO k=1,nchl
					sum_pref=pref_zoo2algae(j,k)*phyto(i,k)
				ENDDO
				DO k=1,nchl

					pref_algae (j,k)=pref_zoo2algae(j,k)*phyto(i,k)/		&
									(sum_pref+pref_zoo2det(j)*POC(i))
					total_feed(i,j)=sum_pref+pref_zoo2det(j)*POC(i)
				ENDDO
				pref_detrius (j)=pref_zoo2det(j)*POC(i)/						&
                        (sum_pref+pref_zoo2det(j)*POC(i))

!-------------- grazing rate determination ------------------------------

				grazing_algae (i,j)=0.0d0
				DO k=1,nchl 	      
	   			grazing_algae (i,j)= grazing_max(j)*pref_algae (j,k)*phyto(i,k)/		&
                   (half_sat_const(j)+total_feed(i,j))
               grazing_detrius(i,j) = grazing_max(j)*pref_detrius (j)*death(i,k)/	&
                  (half_sat_const(j)+total_feed(i,j))    
				ENDDO

!---------------- ZOOP Growth -------------------------------------------	      
				zoop_growth(i,j) = zoop(i,j)*zoogrowth_rate(i,j)*kmax*						&
                  (grazing_algae (i,j)+pref_detrius (j))
            zoop(i,j)=zoop(i,j)+zoop_growth(i,j)*(DFLOAT(nosecs)/86400.0) 	
!---------------- ZOOP Mortality ----------------------------------------------
            temp_ref=10.0
				temp_const=0.05
				zoop_death(i,j) = mortality(j)*DEXP(temp_const*(temp(i)-temp_ref))*zoop(i,j)
				zoop(i,j)=zoop(i,j)-zoop_death(i,j)*(DFLOAT(nosecs)/86400.0)

!---------------- ZOOP Predation ----------------------------------------------
!		  j=1 copepods and j = 2 cladocerans
      		IF (j.eq.1) THEN
					predation(i,j)=pred1*(zoop(i,j)**2.0)/(pred2+zoop(i,j))
				ELSEIF(j.eq.2) THEN
					predation(i,j)=pred1*(zoop(i,j)**3.0)/((pred2)**2.0+zoop(i,j)**2.0)
				ENDIF
				zoop(i,j)=zoop(i,j)-predation(i,j)*(DFLOAT(nosecs)/86400.0)
!-------------------------------------------------------------------
!	 IF(depth(i).ge.450.0.and.depth(i).le.460.0d0) THEN
!	 	 print*,zoop_death(i,j),zoop_growth(i,j),predation(i,j)
!	ENDIF
!--------------------------------------------------------------------
			ENDDO	  !j or zooplankton loop
		ENDDO  !i or number of layer loop
!gbs ---------------end of zooplankton growth, death and predation -------------
!gbs --------------- Balance of phytoplankton adn nutrients --------------------
		DO i = 1, ns
			DO j = 1, nchl
	!			wqual (i,3) = wqual (i,3)- grazing_algae (i,1)
	!  			death  (i,j) = death  (i,j) - grazing_detrius(i,1)
    !            print*,'gbszoo', i,j,death(i,j)
			ENDDO
		ENDDO
  99	CONTINUE	
		RETURN
		END SUBROUTINE zooplankton
