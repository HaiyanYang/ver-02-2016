C********************************************************************************************
C********************* SUBROUTINE UEL for phantom-node method *******************************
C********************************************************************************************
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &       PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &       KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     &       LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
	  INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*)
      DIMENSION SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL)
      DIMENSION DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*)
      DIMENSION JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
      DIMENSION PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
C     Variables to be updated:
C     RHS         Residual Vector
C     AMATRX      K matrix
C
C     Variables passed in:
C     COORDS      Nodal coordinate array: COORD(J,N) is jth coord of nth node
C     U           Total accumulated DOF array.  Contains accumulated displacements,
C                 ordered as (u_i^1, u_i^2)
C     DU          Incremental displacements, ordered as
C                 DU(2*(N-1)+I,1) = u_i^n
C     NNODE       No. nodes on element.
C     NDOFEL      No. degrees of freedom for element
C     NPROPS      No. real valued element properties
C     MLVARX      Dimensioning variable.
C     NRHS        No. RHS vectors.
C     MCRD        Largest of max value of COORDINATES parameter or active DOF <3.
C
C
C======================	For 2D solid phantom elements:============================
C
C	NNODE=8		8-point element, assume crack starts in the middle
C
C				4(8)----3(7)	(8)----(7)	 4 ---- 3
C				|	 	|	-->  | ---- | +  | ---- |
C				1(5)----2(6)	 1 ---- 2	(5)----(6)
C								
C	NDOFEL=8 	when no damage is detected. 4*2 
C	NDOFEL=16	when damage is detected. 8*2	
C	MLVARX=1
C	NRHS=1
C   MCRD=2		(X,Y)
C	NSVARS = 49
C	
	DOUBLE PRECISION KMATRX0(8,8),KMATRX1(8,8),KMATRX2(8,8),KMATRX3(8,8)
	DOUBLE PRECISION KU0(8),KU1(8),KU2(8),KU3(8),UR(8),UG(8),XE(8)
	DOUBLE PRECISION UE1(8),UM1(8),UE2(8),UM2(8),UE3(16),UM3(8),VM3(8)
	DOUBLE PRECISION STRAIN0(6),STRAIN1(6),STRAIN2(6)
	DOUBLE PRECISION STRESS0(6),STRESS1(6),STRESS2(6),STRESS3(2,2)
	DOUBLE PRECISION TMATRX3(8,16),TMATRX3T(16,8)
	DOUBLE PRECISION KT3(8,16),K3G(16,16),F3G(16)
	DOUBLE PRECISION TMATRX2(8,8),TMATRX2T(8,8),KT2(8,8),K2G(8,8),F2G(8)
	DOUBLE PRECISION TMATRX1(8,8),TMATRX1T(8,8),KT1(8,8),K1G(8,8),F1G(8)
	DOUBLE PRECISION XYZE0(3,8),XYZE1(3,8),XYZE2(3,8),XYZE3(2,4)
	DOUBLE PRECISION XYZEA(2),XYZEB(2),RA(4),RB(4)
	DOUBLE PRECISION NMATRXA(2,8),NMATRXB(2,8)
	DOUBLE PRECISION DM(2),CFSTAT(2),U0_EFF(2),UF_EFF(2)
	DOUBLE PRECISION E1,E2,E3,G12,G23,G13,V12,V23,V13,THETA,XT,XC,YT,YC
	DOUBLE PRECISION S,G1T,G1C,GNC,GSC,BKETA
	DOUBLE PRECISION DMEQ(2),KNNU(2),KSSU(2),KNNUEQ(2),KSSUEQ(2)
	INTEGER CNCR(8),CNCG(8),CNC1(8),CNC2(8),CNC3(16),IINC,ISTEP0,IINC0
	CHARACTER(len=15):: ELEMENT
	DATA IUEL/0/, FSTAT/0.D0/, CFSTAT/0.D0, 0.D0/, DM/0.D0, 0.D0/
	DATA U0_EFF/0.D0, 0.D0/, UF_EFF/0.D0, 0.D0/
	DATA DMEQ/0.D0, 0.D0/, KNNU/0.D0, 0.D0/, KSSU/0.D0, 0.D0/
	DATA KNNUEQ/0.D0, 0.D0/, KSSUEQ/0.D0, 0.D0/
c	DATA XYZCA/0.D0, 0.D0/, XYZCB/0.D0, 0.D0/
	SAVE IUEL,FSTAT,CFSTAT,DM,U0_EFF,UF_EFF,STRESS0,IINC,ISTEP0,IINC0
	SAVE NMATRXA,NMATRXB,XYZEA,XYZEB,DMEQ,KNNU,KSSU,KNNUEQ,KSSUEQ
C SVARS are not saved by ABAQUS!! need to use DATA & SAVE as defined above
c to store updated information
	
	IF (IUEL .EQ. 0) THEN
		DO I=1,NSVARS
			SVARS(I)=0.D0
		END DO
		IUEL = 1
c		WRITE(7,*) 'FIRST CALL'
c		WRITE(7,*) 'NSVARS:', NSVARS
	END IF
C
C		DO K=1,5
C		WRITE(7,*) 'SVARS(K):', K, SVARS(K)
C		END DO
C	
C----- Read material properties-----------------------
C-- Read from UEL PROPERTY line in input file --------
c-- floating var. first, then integer, 8 var. per line
      	E1   = PROPS(1)
      	E2   = PROPS(2)
      	E3   = PROPS(3)
      	G12  = PROPS(4)
      	G23  = PROPS(5)
      	G13  = PROPS(6)
      	V12  = PROPS(7)
      	V23  = PROPS(8)
      	V13  = PROPS(9)
      	THETA = PROPS(10) !- material orientation(fibre dir.)
      	XT = PROPS(11)
      	XC = PROPS(12)
      	YT = PROPS(13)
      	YC = PROPS(14)
      	S = PROPS(15)
      	G1T = PROPS(16)
      	G1C = PROPS(17)
      	GNC = PROPS(18)
      	GSC = PROPS(19)
      	BKETA = PROPS(20) !- ETA for BK formula of mixed-mode Gc
c       	 		      	
       WRITE(7,*) 'FSTAT:', FSTAT
c      			
!----- Initialize AMATRX and RHS -------------
       DO  I = 1,NDOFEL
         RHS(I,1)=0.D0
         DO  J = 1,NDOFEL
           AMATRX(J,I) = 0.D0
         END DO
       END DO     
!-Real U and Ghost U	
       	DO I = 1,4
       		UR(I*2-1)=U((I-1)*4+1)
       		UR(I*2)=U((I-1)*4+2)
       		UG(I*2-1)=U((I-1)*4+3)
       		UG(I*2)=U((I-1)*4+4)
       	END DO       	
      	!- CONNECTIVITY MATRIX
      	DO I = 1,4
       		CNCR(I*2-1)=(I-1)*4+1
       		CNCR(I*2)=(I-1)*4+2
       		CNCG(I*2-1)=(I-1)*4+3
       		CNCG(I*2)=(I-1)*4+4
       	END DO 
c       	DO K=1,8
c		WRITE(7,*) 'UR:', K, UR(K)
c		WRITE(7,*) 'UG:', K, UG(K)
c		END DO    	      
c       
C===== WHEN THE ELEMENT IS INTACT: normal 4-node quadrilateral elm ==========              
      IF (FSTAT .EQ. 0.D0) THEN
      	JTYPE0=204
      	NNODE0=4
      	NDOFEL0=8
      	!-ghost nodes tied to real nodes
       		DO I = 1,4
       	 		U((I-1)*4+3)=U((I-1)*4+1)
       	 		U((I-1)*4+4)=U((I-1)*4+2)
       		END DO
			DO I=1,8
				UG(I)=UR(I)
			END DO
		!-----Initialize and value XYZE matrix -------
      	DO I=1,3
        	DO J=1,8
          	XYZE0(I,J)=0.D0
        	end do
      	end do
      	DO I=1,2
        	DO J=1,4
          	XYZE0(I,J)=COORDS(I,J)
        	end do
      	end do
C      	!----- Import values from U to U0 ------------
C      	DO I=1,8 !-Only the first 4 real nodes
C      		UE0(I)=UR(I)
C      	END DO
		!----- Initialize STRESS and STRAIN ----------
		DO I=1,6
			STRESS0(I)=0.0
			STRAIN0(I)=0.0
		END DO
		!----- Get K matrix, KU, stress and strain -----------------
		CALL NORMALELM(E1,E2,E3,G12,G13,G23,V12,V13,V23,THETA,
     1	JTYPE0,NNODE0,KMATRX0,KU0,NDOFEL0,XYZE0,UR,STRAIN0,STRESS0)
     	!----- Check for failure initiation ------------------------
     	DO I=1,3
     	WRITE(7,*) 'STRESS0:', I, STRESS0(I)
     	END DO
     	NST=3
     	CALL HASHIN(NST,STRESS0,YT,S,FSTAT2)
C       	SVARS(1)=FSTAT
C       		     	
     	IF (FSTAT2 .EQ. 0.D0) THEN !- no failure, export AMATRX & RHS
C       		SVARS(2)=0.D0 !-DAMAGE = 0
     	!- Assemble to global AMATRX and RHS
     		DO I = 1,8
     		DO J=1,8
     	AMATRX(CNCR(I),CNCR(J))=AMATRX(CNCR(I),CNCR(J))+KMATRX0(I,J)
     		END DO
     		RHS(CNCR(I),1)=RHS(CNCR(I),1)-KU0(I)
     		END DO
       	!-------- Store stress&strain in SDVs -------
       		DO I=1,6
       		SVARS(5+I)=STRESS0(I)
       		SVARS(11+I)=STRAIN0(I)
       		END DO
       	ELSE !- failure initiated, FSTAT2=1
       		FSTAT=FSTAT2
			SVARS(1)=FSTAT
			IINC=KINC
			ISTEP0=KSTEP
			IINC0=KINC
       	!-------- Determine crack surface -----------
       	!--Determine XYZCA=sum(Ra_I*XYZ_I),XYZCB=sum(Rb_I*XYZ_I)
       	!--then ABCCA=sum(Ra_I*ABC_I),ABCCB=sum(Rb_I*ABC_I)
       	!-- ABCCA(2): crack tip A coords in ref. geo.
       	!-- ABCCB(2): crack tip B coords in ref. geo.
       	!-- Here assume crack lies in the middle
C				4(8)----3(7)	(8)----(7)	 4 ---- 3
C				|	 	|	-->  | ---- | +  | ---- |
C				1(5)----2(6)	 1 ---- 2	(5)----(6)
			RA(1)=0.5
			RA(2)=0.D0
			RA(3)=0.D0
			RA(4)=0.5
			RB(1)=0.D0
			RB(2)=0.5
			RB(3)=0.5
			RB(4)=0.D0
			DO I=1,2
				DO J=1,4
				NMATRXA(I,J)=0.D0
				NMATRXB(I,J)=0.D0
				END DO
			END DO
			DO I=1,2
				DO J=1,4
				NMATRXA(I,I+(J-1)*2)=RA(J)
				NMATRXB(I,I+(J-1)*2)=RB(J)
				END DO
			END DO
			!- Nodal coords vector
      		XE(1)=COORDS(1,1)
      		XE(2)=COORDS(2,1)
      		XE(3)=COORDS(1,2)
      		XE(4)=COORDS(2,2)
      		XE(5)=COORDS(1,3)
      		XE(6)=COORDS(2,3)
      		XE(7)=COORDS(1,4)
      		XE(8)=COORDS(2,4)  
      		!- Ghost nodes' coords
      		CALL MATRIX_MUL(NMATRXA,XE,XYZEA,2,8,1)
      		CALL MATRIX_MUL(NMATRXB,XE,XYZEB,2,8,1)
       	END IF
C       	
      END IF !-FSTAT=0
        WRITE(7,*) 'FSTAT:', FSTAT
C
C===== WHEN THE ELEMENT IS DAMAGED: a mesh of 2 normal quadri elm ==
c===== with 1 interface elm ========================================      	            
      IF (FSTAT .GE. 1.D0) THEN
C		
C------- Elm A information-
      	JTYPE1=204
      	NNODE1=4
      	NDOFEL1=8
C				4(8)----3(7)	(8)----(7)	 4 ---- 3
C				|	 	|	-->  | ---- | +  | ---- |
C				1(5)----2(6)	 1 ---- 2	(5)----(6)
     	!- Nodal disp. vector of elm 1
      	UE1(1)=UR(1) !-U(1)
      	UE1(2)=UR(2) !-U(2)
      	UE1(3)=UR(3) !-U(5)
      	UE1(4)=UR(4) !-U(6)
      	UE1(5)=UG(5) !-U(11)
      	UE1(6)=UG(6) !-U(12)
      	UE1(7)=UG(7) !-U(15)
      	UE1(8)=UG(8) !-U(16)
      	!- CONNECTIVITY MATRIX      	
		CNC1(1)=CNCR(1)
		CNC1(2)=CNCR(2)
		CNC1(3)=CNCR(3)
		CNC1(4)=CNCR(4)
		CNC1(5)=CNCG(5)
		CNC1(6)=CNCG(6)
		CNC1(7)=CNCG(7)
		CNC1(8)=CNCG(8)
      	!- Initialize coords matrix of elm 1
      	DO I=1,3
        	DO J=1,8
          	XYZE1(I,J)=0.D0
        	end do
      	end do
      	!- Real nodes' coords
      	DO I=1,2
          	XYZE1(I,1)=COORDS(I,1)
          	XYZE1(I,2)=COORDS(I,2)
C          	WRITE(7,*) 'XYZE1(I,1):', I, XYZE1(I,1)
C          	WRITE(7,*) 'XYZE1(I,2):', I, XYZE1(I,2)
          	XYZE1(I,3)=XYZEB(I)
          	XYZE1(I,4)=XYZEA(I)
C          	WRITE(7,*) 'XYZE1(I,3):', I, XYZE1(I,3)
C          	WRITE(7,*) 'XYZE1(I,4):', I, XYZE1(I,4)
      	end do
      	!- Initialize UM1 and TMATRX1
      	DO I=1,8
C      		UM1(I)=0.D0
      		DO J=1,8
      			TMATRX1(I,J)=0.D0
      		END DO
      	END DO
      	!- Assign value to TMATRX1
      	DO I=1,4
      		TMATRX1(I,I)=1.D0
      	END DO
      	DO I=1,8
      		TMATRX1(5,I)=NMATRXB(1,I)
      		TMATRX1(6,I)=NMATRXB(2,I)
      		TMATRX1(7,I)=NMATRXA(1,I)
      		TMATRX1(8,I)=NMATRXA(2,I)
      	END DO
      	!- Calculate UM1:
      	CALL MATRIX_MUL(TMATRX1,UE1,UM1,8,8,1)
c      	DO I=1,8
C      	WRITE(7,*) 'UE1(I):', I, UE1(I)
c      	WRITE(7,*) 'UM1(I):', I, UM1(I)
c      	END DO
		!- Calculate elm1's K matrix and F vector    	
		CALL NORMALELM(E1,E2,E3,G12,G13,G23,V12,V13,V23,THETA,
     1	JTYPE1,NNODE1,KMATRX1,KU1,NDOFEL1,XYZE1,UM1,STRAIN1,STRESS1)
C     	
		!- Transfer KMATRIX1 and KU1 into matrix w.r.t UE1 -----
		CALL MATRIX_TRANSPOSE(TMATRX1,TMATRX1T,8,8)
		CALL MATRIX_MUL(KMATRX1,TMATRX1,KT1,8,8,8)
		CALL MATRIX_MUL(TMATRX1T,KT1,K1G,8,8,8)
		CALL MATRIX_MUL(TMATRX1T,KU1,F1G,8,8,1)
C		DO I=1,8
C     	WRITE(7,*) 'F1G:', I, F1G(I)
C     	END DO
     	!- Assemble to global AMATRX and RHS
     	DO I = 1,8
     	DO J=1,8
     	AMATRX(CNC1(I),CNC1(J))=AMATRX(CNC1(I),CNC1(J))+K1G(I,J)
     	END DO
     	RHS(CNC1(I),1)=RHS(CNC1(I),1)-F1G(I)
     	END DO
     	!- Update stresses & strains
     	DO I=1,6
       		SVARS(17+I)=STRESS1(I)
       		SVARS(23+I)=STRAIN1(I)
       	END DO      	
C      	
C------- Elm B information-      	
      	JTYPE2=204
      	NNODE2=4
      	NDOFEL2=8
C				4(8)----3(7)	(8)----(7)	 4 ---- 3
C				|	 	|	-->  | ---- | +  | ---- |
C				1(5)----2(6)	 1 ---- 2	(5)----(6)  	
      	!- Nodal disp. vector of elm 2
      	UE2(1)=UG(1) !-U3
      	UE2(2)=UG(2) !-U4
      	UE2(3)=UG(3) !-U7
      	UE2(4)=UG(4) !-U8
      	UE2(5)=UR(5) !-U9
      	UE2(6)=UR(6) !-U10
      	UE2(7)=UR(7) !-U13
      	UE2(8)=UR(8) !-U14
      	!- CONNECTIVITY MATRIX      	
		CNC2(1)=CNCG(1)
		CNC2(2)=CNCG(2)
		CNC2(3)=CNCG(3)
		CNC2(4)=CNCG(4)
		CNC2(5)=CNCR(5)
		CNC2(6)=CNCR(6)
		CNC2(7)=CNCR(7)
		CNC2(8)=CNCR(8)
      	!- Initialize coords matrix of elm 2      	
      	DO I=1,3
        	DO J=1,8
          	XYZE2(I,J)=0.D0
        	end do
      	end do
      	!- Real nodes' coords
      	DO I=1,2
          	XYZE2(I,1)=XYZEA(I)
          	XYZE2(I,2)=XYZEB(I)
C          	WRITE(7,*) 'XYZE2(I,1):', I, XYZE2(I,1)
C          	WRITE(7,*) 'XYZE2(I,2):', I, XYZE2(I,2)      	
          	XYZE2(I,3)=COORDS(I,3)!-real node
          	XYZE2(I,4)=COORDS(I,4)!-real node
C          	WRITE(7,*) 'XYZE2(I,3):', I, XYZE2(I,3)
C          	WRITE(7,*) 'XYZE2(I,4):', I, XYZE2(I,4)
      	end do      	
      	!- Initialize UM2 and TMATRX2
      	DO I=1,8
C      		UM2(I)=0.D0
      		DO J=1,8
      			TMATRX2(I,J)=0.D0
      		END DO
      	END DO
      	!- Assign value to TMATRX2
      	DO I=5,8
      		TMATRX2(I,I)=1.D0
      	END DO
      	DO I=1,8
      		TMATRX2(1,I)=NMATRXA(1,I)
      		TMATRX2(2,I)=NMATRXA(2,I)
      		TMATRX2(3,I)=NMATRXB(1,I)
      		TMATRX2(4,I)=NMATRXB(2,I)
      	END DO
      	!- Calculate U2:
      	CALL MATRIX_MUL(TMATRX2,UE2,UM2,8,8,1)
c      	DO I=1,8
C      	WRITE(7,*) 'UE2(I):', I, UE2(I)
c      	WRITE(7,*) 'UM2(I):', I, UM2(I)
c      	END DO 
		!- Calculate elm2's K matrix and F vector       	
     	CALL NORMALELM(E1,E2,E3,G12,G13,G23,V12,V13,V23,THETA,
     1	JTYPE2,NNODE2,KMATRX2,KU2,NDOFEL2,XYZE2,UM2,STRAIN2,STRESS2)
		!- Transfer KMATRIX2 and KU2 into matrix w.r.t U -----
		CALL MATRIX_TRANSPOSE(TMATRX2,TMATRX2T,8,8)
		CALL MATRIX_MUL(KMATRX2,TMATRX2,KT2,8,8,8)
		CALL MATRIX_MUL(TMATRX2T,KT2,K2G,8,8,8)
		CALL MATRIX_MUL(TMATRX2T,KU2,F2G,8,8,1)
C		DO I=1,8
C     	WRITE(7,*) 'F2G:', I, F2G(I)
C     	END DO
     	!- Assemble to global AMATRX and RHS
     	!- Assemble to global AMATRX and RHS
     	DO I = 1,8
     	DO J=1,8
     	AMATRX(CNC2(I),CNC2(J))=AMATRX(CNC2(I),CNC2(J))+K2G(I,J)
     	END DO
     	RHS(CNC2(I),1)=RHS(CNC2(I),1)-F2G(I)
     	END DO
     	!- Update stresses & strains     	
     	DO I=1,6
       		SVARS(29+I)=STRESS2(I)
       		SVARS(35+I)=STRAIN2(I)
       	END DO
C ------ Elm C (cohesive) information-
      	NNODE3=4
      	NDOFEL3=8
      	JTYPE3=404
C				4(8)----3(7)	(8)----(7)	 4 ---- 3
C				|	 	|	-->  | ---- | +  | ---- |
C				1(5)----2(6)	 1 ---- 2	(5)----(6)
      	DO I=1,2
        	DO J=1,4
          	XYZE3(I,J)=0.D0
        	end do
      	end do
      	!- Nodal coords of elm 3
      	DO I=1,2
          	XYZE3(I,1)=XYZE1(I,4)
C          	WRITE(7,*) 'XYZE3(I,1):', I, XYZE3(I,1)
      	    XYZE3(I,2)=XYZE1(I,3)
C      	    WRITE(7,*) 'XYZE3(I,2):', I, XYZE3(I,2)
          	XYZE3(I,3)=XYZE2(I,2)
C          	WRITE(7,*) 'XYZE3(I,3):', I, XYZE3(I,3)
          	XYZE3(I,4)=XYZE2(I,1)
C          	WRITE(7,*) 'XYZE3(I,4):', I, XYZE3(I,4)
      	end do
      	!- Nodal disp. vector of elm 2      	
      	DO I=1,8
      		UE3(I)=UE1(I)
      		UE3(8+I)=UE2(I)
      	END DO
      	!- CONNECTIVITY MATRIX      	
      	DO I=1,8
      		CNC3(I)=CNC1(I)
      		CNC3(8+I)=CNC2(I)
      	END DO
      	!- Initialize UM3 and TMATRX
      	DO I=1,8
      		UM3(I)=0.D0
      		VM3(I)=0.D0
      		DO J=1,16
      			TMATRX3(I,J)=0.D0
      		END DO
      	END DO
      	!- TMATRX3: UM3=TMATRX3*UE3
      	DO I=1,8
      		TMATRX3(1,I)=NMATRXA(1,I)
      		TMATRX3(2,I)=NMATRXA(2,I)
      		TMATRX3(3,I)=NMATRXB(1,I)
      		TMATRX3(4,I)=NMATRXB(2,I)
      		TMATRX3(5,8+I)=NMATRXB(1,I)
      		TMATRX3(6,8+I)=NMATRXB(2,I)
      		TMATRX3(7,8+I)=NMATRXA(1,I)
      		TMATRX3(8,8+I)=NMATRXA(2,I)
      	END DO
      	!- UM3=TMATRX*U
      	CALL MATRIX_MUL(TMATRX3,UE3,UM3,8,16,1)
      	!- VM3=TMATRX*V
C      	CALL MATRIX_MUL(TMATRX3,V,VM3,8,16,1)
      	!- SDVs from previous step
C      	CFSTAT(1)=SVARS(2)
C     	CFSTAT(2)=SVARS(3)
C     	DM(1)=SVARS(4)
C     	DM(2)=SVARS(5)
C     	STRESS3(1,1)=SVARS()
C     	STRESS3(2,1)=SVARS()
C     	STRESS3(1,2)=SVARS()
C     	STRESS3(2,2)=SVARS()
C     	U0_EFF(1)=SVARS(46)
C     	U0_EFF(2)=SVARS(47)
C     	UF_EFF(1)=SVARS(48)
C     	UF_EFF(2)=SVARS(49)
     	DO I=1,2
     	DO J=1,2
     		STRESS3(I,J)=0.D0
     	END DO
     	END DO
C     	WRITE(7,*) 'DM1:', DM(1)
C     	WRITE(7,*) 'DM2:', DM(2)
c		WRITE(7,*) 'IINC:', IINC
c		WRITE(7,*) 'KINC:', KINC
      	!- Calculate KMATRIX3,KU3,STRESS,DAMAGE
     	CALL INTERFACE(STRESS0,STRESS1,STRESS2,E2,G12,YT,S,GNC,
     &		GSC,BKETA,UM3,VM3,STRESS3,CFSTAT,DM,U0_EFF,UF_EFF,
     &		KU3,KMATRX3,XYZE3,IINC,KINC,IINC0,DTIME,DMEQ,KNNU,
     &		KSSU,KNNUEQ,KSSUEQ,KSTEP,ISTEP0)
     	IINC=KINC
c     	WRITE(7,*) 'CFSTAT1:', CFSTAT(1)
c     	WRITE(7,*) 'CFSTAT2:', CFSTAT(2)
c     	WRITE(7,*) 'DM1:', DM(1)
c     	WRITE(7,*) 'DM2:', DM(2)
c     	WRITE(7,*) 'U0_EFF(1)', U0_EFF(1)
c     	WRITE(7,*) 'U0_EFF(2)', U0_EFF(2)
c     	WRITE(7,*) 'UF_EFF(1)', UF_EFF(1)
c     	WRITE(7,*) 'UF_EFF(2)', UF_EFF(2)
     	!- Update to SDVs
     	SVARS(2)=CFSTAT(1)
     	SVARS(3)=CFSTAT(2)
     	SVARS(4)=DM(1)
     	SVARS(5)=DM(2)
     	SVARS(42)=STRESS3(1,1)
     	SVARS(43)=STRESS3(2,1)
     	SVARS(44)=STRESS3(1,2)
     	SVARS(45)=STRESS3(2,2)
C     	SVARS(46)=U0_EFF(1)
C     	SVARS(47)=U0_EFF(2)
C     	SVARS(48)=UF_EFF(1)
C     	SVARS(49)=UF_EFF(2)
		!- Transfer KMATRIX3 and KU3 into matrix w.r.t UE3 -----
		CALL MATRIX_TRANSPOSE(TMATRX3,TMATRX3T,8,16)
		CALL MATRIX_MUL(KMATRX3,TMATRX3,KT3,8,8,16)
		CALL MATRIX_MUL(TMATRX3T,KT3,K3G,16,8,16)
		CALL MATRIX_MUL(TMATRX3T,KU3,F3G,16,8,1)			
C		DO I=1,8
C    	WRITE(7,*) 'F3G:', I, F3G(I)
C     	END DO
		!- Assemble to global AMATRX and RHS
     	DO I = 1,16
     	DO J = 1,16
     	AMATRX(CNC3(I),CNC3(J))=AMATRX(CNC3(I),CNC3(J))+K3G(I,J)
     	END DO
     	RHS(CNC3(I),1)=RHS(CNC3(I),1)-F3G(I)
     	END DO
C     	
       	DO I=2,3
       	WRITE(7,*)'STRESS1:', I, STRESS1(I)
       	WRITE(7,*)'STRESS2:', I, STRESS2(I)
       	WRITE(7,*)'STRESS3:', I-1, STRESS3(I-1,1)
     	END DO   	
c
      END IF !-FSTAT=1
C
C		DO K=1,5
C		WRITE(7,*) 'SVARS(K):',K, SVARS(K)
C		END DO

       RETURN
       END
C********************************************************************************************
C********************* END SUBROUTINE UEL ***************************************************
C********************************************************************************************
C
C
C
c-----------------SUBROUTINE ELM_TYPE ---------------------------
C===== Subroutine to determine element type and =================
c===== dimension, no.strain, no.integration pnt =================

	  SUBROUTINE ELM_TYPE(JTYPE,ELEMENT,NDIM,J2DTYPE,NST,NINTP)
c	  
      INCLUDE 'ABA_PARAM.INC'
c     
	  CHARACTER(len=15):: ELEMENT	
c 	JTYPE: 		elm number in input file
c	ELEMENT: 	elm name
c	NDIM: 		geometry dimension
c	J2DTYPE: 	2D elm type(plain stress/strain)
c	NST:		no. of strains/stresses (3 for 2D, 6 for 3D)
c	NINTP:		no. of integration points	  
	  IF (JTYPE .LT. 200) THEN
C		1D BAR ELEMENT
		ndim = 1
		NST = 1 !-no. of stresses/strains
		IF (JTYPE .EQ. 102) THEN
			ELEMENT = 'bar'
			NINTP = 1
		END IF
      ELSE IF (JTYPE .GT. 200 .AND. JTYPE .LT. 300) THEN
C       2D element
C       JTYPE= 203 & 213 -- 3-node triangle
C            = 206 & 216 -- 6-node triangle
C            = 204 & 214 -- 4-node quadrilateral
C            = 208 & 218 -- 8-node quadrilateral
		ndim = 2
		NST = 3
        IF (JTYPE .LT. 210) THEN
C         Plane strain element
          J2DTYPE=0
        ELSE
C         Plane stress element
          J2DTYPE=1
        END IF
C      
        IF (JTYPE .EQ. 203 .OR. JTYPE .EQ. 206 .OR.
     1      JTYPE .EQ. 213 .OR. JTYPE .EQ. 216) THEN
     		ELEMENT ='triangle'
        	IF (JTYPE .EQ. 203 .OR. JTYPE .EQ. 213) NINTP = 3
c          	IF (JTYPE .EQ. 206 .OR. JTYPE .EQ. 216) NINTP = 6
        ELSE
        	ELEMENT ='quadrilateral'
        	IF (JTYPE .EQ. 204 .OR. JTYPE .EQ. 214) NINTP = 4
c      		IF (JTYPE .EQ. 208 .OR. JTYPE .EQ. 218) NINTP = 9
       	END IF
C       	
      ELSE IF (JTYPE .GT. 300 .AND. JTYPE .LT. 400) THEN
C       3D element
C       JTYPE= 304 -- 4-node tetrahedron
c			 = 306 -- 6-node prism
C            = 308 -- 8-node hexahedron
C            = 320 -- 20-node hexahedron
		ndim = 3
		NST = 6
		IF (JTYPE .EQ. 304) THEN
		  ELEMENT ='tetrahedron'
		  NINTP = 1 !-constant stress elm
		ELSE IF (JTYPE .EQ. 306) THEN
          ELEMENT ='prism'
          NINTP = 6
        ELSE
          ELEMENT ='hexahedron'
          IF (JTYPE.EQ.308) NINTP = 8
c          IF (JTYPE.EQ.320) NINTP = 27
        END IF
C        
      ELSE IF (JTYPE .GT. 400) THEN
C		INTERFACE element
c		JTYPE=404 -- 4-node 2D interface
c			 =406 -- 6-node 3D interface (triangular)
c			 =408 -- 8-node 3D interface (rectangular)
        IF (JTYPE .EQ. 404) THEN 
        	ELEMENT = 'int2d4'
        	NINTP = 2
        	NST = 2
        ELSE IF (JTYPE .EQ. 406) THEN
          	ELEMENT ='int3d6'
		  	NINTP = 3
		  	NST = 3
		ELSE IF (JTYPE .EQ. 408) THEN
		  	ELEMENT ='int3d8'
		  	NINTP = 4
		  	NST = 3
		END IF
      END IF !END DEERMINING ELM TYPE
C     
      RETURN
      END
c-------------------------------------------------------------------------
C
C----------------------- Subroutine DEEMAT--------------------------------
C D-matrix formulation
C Returns the elastic dee-matrix for isotropic/orthotropic material
C-------------------------------------------------------------------------
      subroutine deemat(E1,E2,E3,G12,G23,G13,V12,V23,V13,
     & ndim,J2DTYPE,dee)
C
      INCLUDE 'ABA_PARAM.INC'
C
      parameter (zero=0.d0,one=1.d0,two=2.d0,one_half=1.d0/2.d0)      
      DIMENSION DEE(6,6)
C
      DO I=1,6
        DO J=1,6
          DEE(I,J) = zero
        END DO
      END DO
C
        v21= e2*v12/e1
        v31= e3*v13/e1
        v32= e3*v23/e2
C
        del= one-v12*v21-v13*v31-v23*v32-two*v21*v32*v13
        if (ndim .eq. 2) then
          if (j2dtype .eq. 0) then
            dee(1,1)= e1/(one-v12*v21)
            dee(1,2)= e1*v21/(one-v12*v21)
            dee(2,2)= e2/(one-v12*v21)
            dee(2,1)= dee(1,2)
            dee(3,3)= g12
          else
            dee(1,1)= e1*(one-v23*v32)/del
            dee(1,2)= e1*(v21+v23*v31)/del
            dee(2,2)= e2*(one-v13*v31)/del
            dee(2,1)= dee(1,2)
            dee(3,3)= g12
          end if
        else if (ndim .eq. 3) then
          del = del/e1/e2/e3
          dee(1,1)= (one-v23*v32)/e2/e3/del
          dee(1,2)= (v21+v23*v31)/e2/e3/del
          dee(1,3)= (v31+v21*v32)/e2/e3/del
          dee(2,1)= dee(1,2)
          dee(2,2)= (one-v13*v31)/e1/e3/del
          dee(2,3)= (v32+v12*v31)/e1/e3/del
          dee(3,1)= dee(1,3)
          dee(3,2)= dee(2,3)
          dee(3,3)= (one-v12*v21)/e1/e2/del
          dee(4,4)= g12
          dee(5,5)= g23
          dee(6,6)= g13
        end if
C
C
      return
      end
C
C
C------------------- Subroutine deemat_global ----------------------------
C To transform material matrix from local to global coordinate	
C-------------------------------------------------------------------------
      subroutine DEEMAT_GLOBAL(NST,DEE,DEEG,THETA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      parameter (zero=0.d0,two=2.d0,four=4.d0,pi=3.14159265359D0)      
      DIMENSION DEE(6,6), DEEG(6,6)
      DOUBLE PRECISION M, N
C
      DO I=1,6
        DO J=1,6
          deeg(I,J)=ZERO
        END DO
      END DO
C
      m=cos(pi*theta/180.d0)
      n=sin(pi*theta/180.d0)

      IF (NST .EQ. 3) THEN
        deeg(1,1) = m*m*m*m*dee(1,1) + two*m*m*n*n*(dee(1,2)
     &            + two*dee(3,3)) + n*n*n*n*dee(2,2)
        deeg(1,2) = n*n*m*m*(dee(1,1) + dee(2,2) - four*dee(3,3))
     &            + (n*n*n*n+m*m*m*m)*dee(1,2)
        deeg(2,1) = deeg(1,2)
        deeg(2,2) = n*n*n*n*dee(1,1) + two*m*m*n*n*(dee(1,2)
     &            + two*dee(3,3)) + m*m*m*m*dee(2,2)
        deeg(1,3) = n*m*(m*m*(dee(1,1) - dee(1,2) - two*dee(3,3))
     &            + n*n*(dee(1,2) - dee(2,2) + two*dee(3,3)))
        deeg(3,1) = deeg(1,3)
        deeg(2,3) = n*m*(n*n*(dee(1,1) - dee(1,2) - two*dee(3,3))
     &            + m*m*(dee(1,2) - dee(2,2) + two*dee(3,3)))
        deeg(3,2) = deeg(2,3)
        deeg(3,3) = m*m*n*n*(dee(1,1)+dee(2,2)-2*dee(1,2))
     &			  +(m*m-n*n)**2*dee(3,3)
      ELSE IF (NST .EQ. 6) THEN
        deeg(1,1) = m*m*m*m*dee(1,1) + two*m*m*n*n*(dee(1,2)
     &            + two*dee(4,4)) + n*n*n*n*dee(2,2)
        deeg(1,2) = n*n*m*m*(dee(1,1) + dee(2,2) - four*dee(4,4))
     &            + (n*n*n*n+m*m*m*m)*dee(1,2)
        deeg(2,1) = deeg(1,2)
        deeg(2,2) = n*n*n*n*dee(1,1) + two*m*m*n*n*(dee(1,2)
     &            + two*dee(4,4)) + m*m*m*m*dee(2,2)
        deeg(1,3) = m*m*dee(1,3) + n*n*dee(2,3)
        deeg(3,1) = deeg(1,3)
        deeg(1,4) = n*m*(m*m*(dee(1,1) - dee(1,2) - two*dee(4,4))
     &            + n*n*(dee(1,2) - dee(2,2) + two*dee(4,4)))
        deeg(4,1) = deeg(1,4)
        deeg(2,3) = n*n*dee(1,3) + m*m*dee(2,3)
        deeg(3,2) = deeg(2,3)
        deeg(2,4) = n*m*(n*n*(dee(1,1) - dee(1,2) - two*dee(4,4))
     &            + m*m*(dee(1,2) - dee(2,2) + two*dee(4,4)))
        deeg(4,2) = deeg(2,4)
        deeg(3,3) = dee(3,3)
        deeg(3,4) = m*n*(dee(1,3) - dee(2,3))
        deeg(4,3) = deeg(3,4)
        deeg(5,5) = m*m*dee(5,5) + n*n*dee(6,6)
        deeg(5,6) = m*n*(dee(6,6) - dee(5,5))
        deeg(6,5) = deeg(5,6)
        deeg(6,6) = n*n*dee(5,5) + m*m*dee(6,6)
        deeg(4,4) = n*n*m*m*(dee(1,1) - two*dee(1,2) + dee(2,2))
     &            + (n*n - m*m)*(n*n - m*m)*dee(4,4)
      END IF
C
      RETURN
      End
c
c
C================== Subroutine INITIALGAUSS ==============================
C Initialize Gauss integrating points
C Returns the local coordinates of the integrating points and weights
C=========================================================================
      subroutine initialgauss(element,xyzg,wt,nig)
C
      INCLUDE 'ABA_PARAM.INC'
C
      dimension xyzg(3,8),wt(8),w(3),v(9),s(8,3)
      character(len=15):: element
      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +  five=5.d0,six=6.d0,eight=8.d0,anine=9.d0)
C
      root3= one/sqrt(three)
      r15= sqrt(three)/sqrt(five)
      one_half=one/two
      one_third= one/three
      two_third= two/three
      one_sixth= one/six
      w(1)=five/anine
      w(2)=eight/anine
      w(3)=w(1)
C
!---- Initialize weight array, position matrix ------------------------
      do i=1,8
         wt(i)=zero
	   do j=1,3
	      s(i,j)=0.0; xyzg(j,i)=0.0
	   end do
	end do
!----------------------------------------------------------------------
      n=0
      do i=1,3
        do j=1,3
          n=n+1
          v(n)=w(i)*w(j)
        end do
      end do
C
!----- assign value according to element type & integration point number 
C
C	  
      if (element .eq. 'bar') then
        if (nig .eq. 1) then
          s(1,1) = zero
          wt(1)  = two
        else if (nig .eq. 2) then
          s(1,1)= -root3
          s(2,1)=  root3
          do i=1,nig
            wt(i) = one
          end do
        else if (nig .eq. 3) then
          s(1,1)= -r15
          s(2,1)=  zero
          s(3,1)=  r15
          do i=1,nig
            wt(i) = w(i)
          end do
        end if    
      else if (element .eq. 'triangle') then
        if (nig .eq. 1) then
          s(1,1)= one_third
          s(1,2)= one_third
          wt(1) = one_half
        else if (nig .eq. 3) then
          s(1,1)= one_half
          s(1,2)= one_half
          s(2,1)= one_half
          s(2,2)= zero
          s(3,1)= zero
          s(3,2)= one_half
          wt(1) = one_sixth
          wt(2) = wt(1)
          wt(3) = wt(1)
        end if
      else if (element .eq. 'quadrilateral') then
        if (nig .eq. 1) then
          s(1,1)= zero
          s(1,2)= zero
          wt(1) = four
        else if (nig .eq. 4) then
          s(1,1)= -root3
          s(1,2)= -root3
          s(2,1)=  root3
          s(2,2)= -root3
          s(3,1)= -root3
          s(3,2)=  root3
          s(4,1)=  root3
          s(4,2)=  root3
          do i=1,nig
            wt(i) = one
          end do
        end if     
      else if (element .eq. 'prism') then
        if (nig .eq. 1) then
          s(1,1)= one_third
          s(1,2)= one_third
          s(1,3)= zero
          wt(1) = two
        else if (nig .eq. 6) then
          s(1,1)= one_sixth
          s(1,2)= one_sixth
          s(1,3)=-root3
          s(2,1)= two_third 
          s(2,2)= one_sixth
          s(2,3)=-root3
          s(3,1)= one_sixth 
          s(3,2)= two_third
          s(3,3)=-root3
          s(4,1)= one_sixth
          s(4,2)= one_sixth
          s(4,3)= root3
          s(5,1)= two_third
          s(5,2)= one_sixth
          s(5,3)= root3
          s(6,1)= one_sixth
          s(6,2)= two_third
          s(6,3)= root3
          do i=1,nig
            wt(i) = one_third
          end do
        end if
      else if (element .eq. 'hexahedron') then
        if (nig .eq. 1) then
          s(1,1)= zero
          s(1,2)= zero
          s(1,3)= zero
          wt(1) = eight
        else if (nig .eq. 8) then
          s(1,1)= -root3
          s(1,2)= -root3
          s(1,3)= -root3
          s(2,1)=  root3
          s(2,2)= -root3
          s(2,3)= -root3
          s(3,1)= -root3
          s(3,2)=  root3
          s(3,3)= -root3
          s(4,1)=  root3
          s(4,2)=  root3
          s(4,3)= -root3
          s(5,1)= -root3
          s(5,2)= -root3
          s(5,3)=  root3
          s(6,1)=  root3
          s(6,2)= -root3
          s(6,3)=  root3
          s(7,1)= -root3
          s(7,2)=  root3
          s(7,3)=  root3
          s(8,1)=  root3
          s(8,2)=  root3
          s(8,3)=  root3
          do i=1,nig
            wt(i) = one
          end do
        end if
      else if(element .eq. 'int2d4') then
C
      else if (element .eq. 'int3d6') then
C
      else if (element .eq. 'int3d8') then
C
      end if
!------------------------------------------------------------------------------
      do i=1,3
        do j=1,8
          xyzg(i,j)=s(j,i)
        end do
      end do
C
      return
      end
C==============================================================================
c==============================================================================
C==============================================================================
C
C
C============================================================================
C=========================== SUBROUTINE SHAPE ===============================
C------ Shape functions and derivatives for normal element ------------------
C============================================================================
       SUBROUTINE SHAPE(ELEMENT,NNODE,PQ,I,F,DF)
C
      	INCLUDE 'ABA_PARAM.INC'
C      
	   	dimension PQ(3,8),F(8),DF(8,3)
		CHARACTER(len=15):: ELEMENT
C
!-------initialize F and DF--------------
      	do k=1,8
	   		F(k)=0.0
	   		do j=1,3
	   		DF(k,j)=0.0
	   		end do
	  	end do
!----------------------------------------
C
      	xi=PQ(1,I)
      	eta=PQ(2,I)
      	zeta=PQ(3,I)
C      	
      	if(element.eq.'bar') then
C      	
      	else if(element.eq.'triangle') then
C      	
		else if(ELEMENT.eq.'quadrilateral') then
		    etam= 0.25*(1.0-eta)
        	etap= 0.25*(1.0+eta)
        	xim = 0.25*(1.0-xi)
        	xip = 0.25*(1.0+xi)
c			IF(NNODE.EQ.4) THEN
        	F(1)=0.25*(1.0-xi)*(1.0-eta)
        	F(2)=0.25*(1.0+xi)*(1.0-eta)
        	F(3)=0.25*(1.0+xi)*(1.0+eta)
        	F(4)=0.25*(1.0-xi)*(1.0+eta)
          	DF(1,1) = -etam
          	DF(2,1) =  etam
          	DF(3,1) =  etap
          	DF(4,1) = -etap
          	DF(1,2) = -xim
          	DF(2,2) = -xip
          	DF(3,2) =  xip
          	DF(4,2) =  xim
c			END IF
		else if(element.eq.'tetrahedron') then
C		
		else if(element.eq.'prism') then
C		
		else if(element.eq.'hexahedron') then
C		
		else if(element.eq.'int2d4') then
C		
		else if(element.eq.'int3d6') then
C		
		else if(element.eq.'int3d8') then
C		
		end if
C
       RETURN
       END
C================================================================
c================================================================
C================================================================
c
C============================================================================
C=========================== SUBROUTINE NMATRX ==============================
C============================================================================
       SUBROUTINE NMATRIX(ELEMENT,NNODE,NDOFEL,PQ,NDIM,NMATRX)
c
      	INCLUDE 'ABA_PARAM.INC'
c      
	   	DOUBLE PRECISION PQ(NDIM),NMATRX(NDIM,NDOFEL),F(NNODE)
		CHARACTER(len=15):: ELEMENT
C
!-------initialize TMATRX --------------
      	do k=1,NDOFEL
	   		do j=1,NDIM
	   		NMATRX(j,k)=0.0
	   		end do
	  	end do
!----------------------------------------
		IF(NDIM .EQ. 2) THEN
      	xi=PQ(1)
      	eta=PQ(2)
      	ELSE IF(NDIM .EQ. 3) THEN
      	xi=PQ(1)
      	eta=PQ(2)      	
      	zeta=PQ(3)
      	ELSE
      	xi=PQ(1)      	
      	END IF
C      	
C      	if(element .eq. 'bar') then
C      	
C      	else if(element .eq. 'triangle') then
C      	
C		else if(ELEMENT .eq. 'quadrilateral') then
        	F(1)=0.25*(1.0-xi)*(1.0-eta)
        	F(2)=0.25*(1.0+xi)*(1.0-eta)
        	F(3)=0.25*(1.0+xi)*(1.0+eta)
        	F(4)=0.25*(1.0-xi)*(1.0+eta)
        	DO I=1,4
C        		WRITE(7,*)'F(I):', I, F(I)
        		NMATRX(1,2*I-1)=F(I)
        		NMATRX(2,2*I)=F(I)
        	END DO
C		else if(element .eq. 'tetrahedron') then
C		
C		else if(element .eq. 'prism') then
C		
C		else if(element .eq. 'hexahedron') then
C		
C		else if(element .eq. 'int2d4') then
C		
C		else if(element .eq. 'int3d6') then
C		
C		else if(element .eq. 'int3d8') then
C		
C		end if
C
C      	do k=1,NDOFEL
C	   		do j=1,NDIM
C	   		WRITE(7,*) 'NMATRX(J,K):', J, K, NMATRX(J,K)
C	   		end do
C	  	end do
       RETURN
       END
C================================================================
c================================================================
C================================================================
C
C--------Subroutine MATRIX_MUL----------
C matrix multiplication, a*b=c
C---------------------------------------
      subroutine MATRIX_MUL(a,b,c,l,m,n)
C
      INCLUDE 'ABA_PARAM.INC'
C
      dimension a(l,m), b(m,n), c(l,n)
C
      do i=1,l
        do j=1,n
          c(i,j) = 0.d0 !-initialize
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do
C
      return
      end
C

C-------------------------- SUBROUTINE DETERMINANT -----------------------
C Returns the determinant of a 1x1 2x2 3x3 jacobian matrix		
C-------------------------------------------------------------------------
      subroutine determinant(Ajacobi,DetJ,n)
C
      INCLUDE 'ABA_PARAM.INC'
C
      dimension Ajacobi(3,3)
C
      if (n .eq. 1) then
        DetJ= 1.d0
      else if (n .eq. 2) then
        DetJ= Ajacobi(1,1)*Ajacobi(2,2) - Ajacobi(1,2)*Ajacobi(2,1)
      else if (n .eq. 3) then
        DetJ= Ajacobi(1,1)*(Ajacobi(2,2)*Ajacobi(3,3)-
     +            Ajacobi(3,2)*Ajacobi(2,3))
        DetJ= DetJ- Ajacobi(1,2)*(Ajacobi(2,1)*Ajacobi(3,3)-
     +            Ajacobi(3,1)*Ajacobi(2,3))
        DetJ= DetJ+ Ajacobi(1,3)*(Ajacobi(2,1)*Ajacobi(3,2)-
     +            Ajacobi(3,1)*Ajacobi(2,2))
      end if
C
      return
      end
C
C     
C================== Subroutine BEEMAT ====================================
C B-matrix formulation
C bee-matrix for 2d elasticity (nst=3) or for 3-d (nst=6)
C=========================================================================
C
      subroutine beemat(bee,deriv,nod,nst,NDOFEL)
C
      INCLUDE 'ABA_PARAM.INC'
C
      dimension deriv(8,3), bee(6,NDOFEL)
C
      do i=1,6
        do j=1,NDOFEL
          bee(i,j) = 0.d0
        end do
      end do
C
      if (nst .eq. 3) then
        do m=1,nod
          k= 2*m
          l=k-1
          x=deriv(m,1)
          y=deriv(m,2)
          bee(1,l)= x
          bee(3,k)= x
          bee(2,k)= y
          bee(3,l)= y
        end do
      else if (nst .eq. 6) then
        do m=1,nod
          n= 3*m
          k=n-1
          l=k-1
          x=deriv(m,1)
          y=deriv(m,2)
          z=deriv(m,3)
          bee(1,l)=x
          bee(4,k)=x
          bee(6,n)=x
          bee(2,k)=y
          bee(4,l)=y
          bee(5,n)=y
          bee(3,n)=z
          bee(5,k)=z
          bee(6,l)=z
        end do
      end if
C
      return
      end
C================================================================
c================================================================
C================================================================
C
C-------------------------- Subroutine INVERT ----------------------------
C Inverts a square matrix onto itself					 
C-------------------------------------------------------------------------
      subroutine invert(Ajacobi,DetJ,n)
C
      INCLUDE 'ABA_PARAM.INC'
C
      dimension Ajacobi(3,3),Ajac(3,3)
C
      do i=1,n
        do j=1,n
          Ajac(i,j)=Ajacobi(i,j)
        end do
      end do
C
      if(n .eq. 2) then
        Ajacobi(1,1)=Ajac(2,2)
        Ajacobi(2,1)=-Ajac(2,1)
        Ajacobi(1,2)=-Ajac(1,2)
        Ajacobi(2,2)=Ajac(1,1)
      else if(n .eq. 3) then
        Ajacobi(1,1)=Ajac(2,2)*Ajac(3,3)-Ajac(3,2)*Ajac(2,3)
        Ajacobi(2,1)=Ajac(3,1)*Ajac(2,3)-Ajac(2,1)*Ajac(3,3)
        Ajacobi(3,1)=Ajac(2,1)*Ajac(3,2)-Ajac(3,1)*Ajac(2,2)
        Ajacobi(1,2)=Ajac(3,2)*Ajac(1,3)-Ajac(1,2)*Ajac(3,3)
        Ajacobi(2,2)=Ajac(1,1)*Ajac(3,3)-Ajac(3,1)*Ajac(1,3)
        Ajacobi(3,2)=Ajac(3,1)*Ajac(1,2)-Ajac(1,1)*Ajac(3,2)
        Ajacobi(1,3)=Ajac(1,2)*Ajac(2,3)-Ajac(2,2)*Ajac(1,3)
        Ajacobi(2,3)=Ajac(2,1)*Ajac(1,3)-Ajac(1,1)*Ajac(2,3)
        Ajacobi(3,3)=Ajac(1,1)*Ajac(2,2)-Ajac(2,1)*Ajac(1,2)
      end if
C
      do i=1,n
        do j=1,n
          Ajacobi(i,j)=Ajacobi(i,j)/DetJ
        end do
      end do
C
      return
      end
      

C--------Subroutine MATRIX_TRANSPOSE------
C matrix transpose
C-----------------------------------------
      subroutine MATRIX_TRANSPOSE(a,at,m,n)
C
      INCLUDE 'ABA_PARAM.INC'
      dimension a(m,n), at(n,m)
C
      do i=1,m
        do j=1,n
          at(j,i) = a(i,j)
        end do
      end do
C
      return
      end
C      
C--------- Subroutine transfer_strain-------------------------------------
C Transfer strains form global to material coord sys	
C-------------------------------------------------------------------------
      SUBROUTINE TRANSFER_STRAIN(NST,STRAIN,THETA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STRAIN(6),STMD(6),T(6,6)     
      PARAMETER (pi=3.14159265359D0,ZERO=0.D0,ONE=1.D0,TWO=2.D0)
C
      C=cos(pi*theta/180.D0)
      S=sin(pi*theta/180.D0)
      DO I=1,6
        DO J=1,6
          T(I,J) = ZERO
        END DO
      END DO
      IF (NST .EQ. 3) THEN
        T(1,1)=C*C
        T(1,2)=S*S
        T(1,3)=C*S
        T(2,1)=S*S
        T(2,2)=C*C
        T(2,3)=-C*S
        T(3,1)=-TWO*C*S
        T(3,2)=TWO*C*S
        T(3,3)=C*C-S*S
      ELSE
        T(1,1)=C*C
        T(1,2)=S*S
        T(1,4)=C*S
        T(2,1)=S*S
        T(2,2)=C*C
        T(2,4)=-C*S
        T(3,3)=ONE
        T(4,1)=-TWO*C*S
        T(4,2)=TWO*C*S
        T(4,4)=C*C-S*S
        T(5,5)=C
        T(5,6)=-S
        T(6,5)=S
        T(6,6)=C
      END IF
      CALL MATRIX_MUL(T,STRAIN,STMD,6,6,1)
C
      DO I=1,NST
        STRAIN(I) = STMD(I)
      END DO
C
      RETURN
      End SUBROUTINE TRANSFER_STRAIN
C
C
C----------------------------------------------------------------
c--- Subroutine KMATRIX_N to calculate K-matrix ---
C--- and nodal force vector of normal elm 
c----------------------------------------------------------------
	SUBROUTINE NORMALELM(E1,E2,E3,G12,G13,G23,V12,V13,V23,THETA,
     1	JTYPE,NNODE,KMATRX,KU,NDOFEL,XYZE,U,STRAIN,STRESS)
c     
	  INCLUDE 'ABA_PARAM.INC'
C
	DOUBLE PRECISION KMATRX(NDOFEL,NDOFEL),KU(NDOFEL),U(NDOFEL)
	DOUBLE PRECISION W(8),PQ(3,8)
	DOUBLE PRECISION FN(8),DN(8,3),XYZE(3,8)
	DOUBLE PRECISION JAC(3,3),GN(8,3)
	DOUBLE PRECISION BEE(6,NDOFEL),BEET(NDOFEL,6)
	DOUBLE PRECISION DEE(6,6),DEEG(6,6)
	DOUBLE PRECISION BTD(NDOFEL,6),BTDB(NDOFEL,NDOFEL)
	DOUBLE PRECISION STRESS(6),STRAIN(6)
	DOUBLE PRECISION E1,E2,E3,G12,G13,G23,V12,V13,V23,THETA,DETJ
C	
	CHARACTER(len=15):: ELEMENT
C	
!------ Determine Element Type -------------------------------------
c 	JTYPE: 		elm number in input file
c	ELEMENT: 	elm name
c	NDIM: 		geometry dimension
c	J2DTYPE: 	2D elm type(plain stress/strain)
c	NST:		no. of strains/stresses (3 for 2D, 6 for 3D)
c	NINTP:		no. of integration points
	CALL ELM_TYPE(JTYPE,ELEMENT,NDIM,J2DTYPE,NST,NINTP)
C	
!----- Initialize AMATRX and RHS -------------
       DO  I = 1,NDOFEL
         KU(I)=0.D0
         DO  J = 1,NDOFEL
           KMATRX(J,I) = 0.D0
         END DO
       END DO
!----- Initialize STRESS and STRAIN ----------
		DO I=1,6
			STRESS(I)=0.0
			STRAIN(I)=0.0
		END DO
C
!----- Calculate material local stiffness matrix DEE ----------------------------------------
		CALL DEEMAT(E1,E2,E3,G12,G23,G13,V12,V23,V13,ndim,J2DTYPE,DEE)
C
!-------Calculate material stiffness matrix in global coordinates DEEG-----------------------
		CALL DEEMAT_GLOBAL(NST,DEE,DEEG,THETA)
C		
!------ Set up integration points and weights -----------------------------------------------
       	CALL INITIALGAUSS(ELEMENT,PQ,W,NINTP)
C
C======== Calculate AMATRX = K (exclude interface elm) ======================================
       	DO IG=1,NINTP !-CALCULATE AT EACH INT POINT
C
			CALL SHAPE(ELEMENT,NNODE,PQ,IG,FN,DN) !-to get N and DN
			CALL MATRIX_MUL(XYZE,DN,JAC,3,8,3) !-JAC=XYZE*DN
     		CALL DETERMINANT(JAC,DETJ,NDIM) !-get DET(JAC)
     		CALL INVERT(JAC,DETJ,NDIM) ! invert JAC onto itself
     		CALL MATRIX_MUL(DN,JAC,GN,8,3,3) !-GN=DN*inv(JAC)
     		CALL BEEMAT(BEE,GN,NNODE,NST,NDOFEL) !-calculate B matrix 6*NDOFEL
     		CALL MATRIX_TRANSPOSE(BEE,BEET,6,NDOFEL) !-BEET = transpose(BEE)
     		CALL MATRIX_MUL(BEET,DEEG,BTD,NDOFEL,6,6) !-BTD = BEET * DEEG
     		CALL MATRIX_MUL(BTD,BEE,BTDB,NDOFEL,6,NDOFEL) !-BTDB = BEET * DEEG * BEE
C
       		DO I=1,NDOFEL
          		DO J=1,NDOFEL
              		KMATRX(I,J) = KMATRX(I,J)+BTDB(I,J)*DETJ*W(IG)
          		end do
        	end do
C       	
!-------- Calculate strain in global coords: epsilon = B * U --------------------------------
        	DO I=1,6
           	DO J=1,NDOFEL
              	strain(I) = strain(I)+BEE(I,J)*U(J)/NINTP !-average strain in the elm
           	END DO
        	END DO
C     	
       	END DO !-LOOPED OVER ALL INT POINTS. IG=NINTP
C
C======== Calculate KU ===============================================================
		CALL MATRIX_MUL(KMATRX,U,KU,NDOFEL,NDOFEL,1) !- KU = AMATRX*U = -RHS
C
C======== Transform STRAIN into local coords ======
		CALL TRANSFER_STRAIN(NST,STRAIN,THETA)
C		
C======== Calculate STRESS in local coords ========
		CALL MATRIX_MUL(DEE,STRAIN,STRESS,6,6,1)
C
       RETURN
       END SUBROUTINE NORMALELM
C
C
C----------------------------------------------------------                  
C---------Subroutine HASHIN--------------------------------
C HASHIN FAILURE CRITERIA FOR MATRIX TENSILE FAILURE
C----------------------------------------------------------
      subroutine hashin(nst,STRESS,YT,S,FSTAT)
c
      include 'aba_param.inc'
C      
      DOUBLE PRECISION stress(6), sigma(6), FF
C
		IF(NST .EQ. 3) THEN
        SIGMA(1) = STRESS(1)
        SIGMA(2) = STRESS(2)
        SIGMA(3) = 0.D0
        SIGMA(4) = STRESS(3)
        SIGMA(5) = 0.D0
        SIGMA(6) = 0.D0
      	ELSE
        	DO I=1,6
          	SIGMA(I) = STRESS(I)
        	end do
		END IF
C
	FF=(max(0.d0,SIGMA(2))/YT)**2+(SIGMA(4)/S)**2+(SIGMA(6)/S)**2
C     
	 IF(FF .GE. 1.D0) THEN
	 	FSTAT=1.D0
	 ELSE
	 	FSTAT=0.D0	
	 END IF
C	 
      Return  
      END subroutine hashin
C
C
C      
C
C=========================== SUBROUTINE INTERFACE ===================
C
      SUBROUTINE INTERFACE(STRESS0,STRESS1,STRESS2,EN,ES,YT,S,GNC,
     &	GSC,BKETA,U,V,TRACT,CFSTAT,DM,U0_EFF,UF_EFF,KU,KMATRX,COORDS
     & 	,IINC,KINC,IINC0,DTIME,DMEQ,KNNU,KSSU,KNNUEQ,KSSUEQ,KSTEP,
     &	ISTEP0)
C
	INCLUDE 'ABA_PARAM.INC'
c
      DOUBLE PRECISION KU(8),KMATRX(8,8),COORDS(2,4),U(8),V(8)
	  DOUBLE PRECISION DM(2),CFSTAT(2),U0_EFF(2),UF_EFF(2)
      DOUBLE PRECISION W(2),PQ(2),F(4),DF(4)       
	  DOUBLE PRECISION RNM(2),TANGENT(2),QMATRX(2,2),QMATRXT(2,2)
	  DOUBLE PRECISION UREL(2),VREL(2),FMATRX(2,8),FMATRXT(8,2)
	  DOUBLE PRECISION GAP(2),GAP0(2),VGAP(2),TRACT(2,2),DMATRX(2,2)
	  DOUBLE PRECISION QF(2,8),DQF(2,8),FTQT(8,2),FTQTDQF(8,8)
	  DOUBLE PRECISION FTQTTAU(8),STRESS0(6),S0N,S0T,OFSET(2)
	  DOUBLE PRECISION EN,ES,YT,S,GNC,GSC,BKETA,DET,KNN0,KNNC
	  DOUBLE PRECISION KSS0,KNN,KSS,T_EFF,U_EFF,GN,GS,B,GMC
	  DOUBLE PRECISION DM2,DFAIL,TRACT2(2),DMEQ(2)
	  DOUBLE PRECISION KNNU(2),KSSU(2),KNNUEQ(2),KSSUEQ(2)
	  DOUBLE PRECISION STRESS1(6),STRESS2(6),DTIME,MDTIME
	  INTEGER OSCI
c
	NINTP=2
	NNODE=4
	DFAIL=1.D0
C -- Initial penalty stiffness
c	KNN0 = 10**3*EN
c	KNN0 = EN
	KNN0 = 0.01*EN
C	KSS0 = 10**3*ES
	KSS0 = KNN0
	KNNC = 10**2*EN
C	
	MDTIME=1E-20
	OSCI=0
	IF(IINC .EQ. KINC) THEN
		OSCI=1
	ELSE !-Previous increment has converged
		DO I=1,2 !-update DMEQ, KNNUEQ, KSSUEQ
		DMEQ(I)=DM(I) !-equilibrium DM
		KNNUEQ(I)=KNNU(I)
		KSSUEQ(I)=KSSU(I)
		END DO
	END IF
C
C	The above arrays are used inside the subroutine
C	NINTP		No. integration points
C	PQ(NINTP)	Local coordinates of integration points
C	W(NINTP)	Weights of integration points
C	F(NNODE)	Shape function
C	DF(NNODE)	Derivative of shape function w.r.t local coordinate
C
c		DO I=1,8
c		WRITE(7,*) 'U INTERFACE:', I, U(I)
c		END DO
c		
       DO  I = 1,8
         KU(I)=0.D0
         DO  J = 1,8
           KMATRX(J,I) = 0.D0
         END DO
       END DO
C
C      Set up integration points and weights
       CALL KINTPT(PQ,W,NNODE,NINTP)
C
C		WRITE(7,*) 'W(1):', W(1)
C		WRITE(7,*) 'W(2):', W(2)
       DO IG=1,NINTP !-Looping over all integration points
C
C       Shape functions and derivatives
       	CALL KSHAPE(NNODE,PQ(IG),F)
	   	CALL DSHAPE(NNODE,PQ(IG),DF)
C	   	WRITE(7,*) 'F(1):', F(1)
C	   	WRITE(7,*) 'F(2):', F(2)
C	   	WRITE(7,*) 'DF(1):', DF(1)
C
C        Compute Normal and tangent vectors to boundary at int pt.
C        Plane of interface defined as average of upper and lower
C        surfaces (hopefully coincident, but you never know...)
	   	TANGENT(1)=.5D0*(COORDS(1,2)-COORDS(1,1)
     &					+COORDS(1,3)-COORDS(1,4))
	   	TANGENT(2)=.5D0*(COORDS(2,2)-COORDS(2,1)
     &					+COORDS(2,3)-COORDS(2,4))
C     	WRITE(7,*) 'TAN(1):', TANGENT(1)
C     	WRITE(7,*) 'TAN(2):', TANGENT(2)
c	
	   	CALL KUNITV(TANGENT,DET)
C	   	WRITE(7,*) 'TAN(1):', TANGENT(1)
C     	WRITE(7,*) 'TAN(2):', TANGENT(2)
C     	WRITE(7,*) 'DET:', DET
	   	DET=.5D0*DET !-JAC = 0.5 * TANGENT
c
	   	RNM(1)=-TANGENT(2)
	   	RNM(2)=TANGENT(1)
C	    WRITE(7,*) 'RNM(1):', RNM(1)
C     	WRITE(7,*) 'RNM(2):', RNM(2)
C		Q Matrix: global coords to local coords
C		Q = TRANSPOSE[RNM,TANGENT]
		DO J=1,2
			QMATRX(1,J)=RNM(J)
			QMATRX(2,J)=TANGENT(J)
		END DO
C
C       Relative displacement&velocity, in global coords
       	UREL(1) = 0.D0
       	UREL(2) = 0.D0
	   	VREL(1)=0.D0
	   	VREL(2)=0.D0
C		FMATRX: UREL (at each int pnt) = F*U
		DO J = 1,8
			DO K = 1,2
			FMATRX(K,J)=0.D0
			END DO
		END DO	   	
	   	DO N = 1,2
	   	FMATRX(N,N)= -F(1)
	   	FMATRX(N,N+2)= -F(2)
	   	FMATRX(N,N+4)= F(2)
	   	FMATRX(N,N+6)= F(1)
	   	END DO
C	   	
	   	CALL MATRIX_MUL(FMATRX,U,UREL,2,8,1)
	   	CALL MATRIX_MUL(FMATRX,V,VREL,2,8,1)
c	   	WRITE(7,*) 'UREL(1):', UREL(1)
c       	WRITE(7,*) 'UREL(2):', UREL(2)
c
C
C       GAP(1) is normal separation, GAP(2) is separation in tan dirn
c		GAP = Q * UREL; VGAP = Q * VREL
        GAP0(1) =  RNM(1)*UREL(1)+RNM(2)*UREL(2)
       	GAP0(2) =  TANGENT(1)*UREL(1)+TANGENT(2)*UREL(2)
	   	VGAP(1)= RNM(1)*VREL(1)+RNM(2)*VREL(2)
	   	VGAP(2)= TANGENT(1)*VREL(1)+TANGENT(2)*VREL(2)
		WRITE(7,*) 'GAP0(1):', GAP0(1)
       	WRITE(7,*) 'GAP0(2):', GAP0(2)
C -- STRESS0 in local coords --
		S0N=STRESS0(1)*RNM(1)**2+STRESS0(2)*RNM(2)**2+
     &	2*STRESS0(3)*RNM(1)*RNM(2) !-normal
     	S0T=(STRESS0(1)-STRESS0(2))*RNM(1)*RNM(2)
     &	+STRESS0(3)*(RNM(2)**2-RNM(1)**2) !-shear
C -- Initial offsets for cohesive laws 
		IF(S0N .GT. 0.D0) THEN
     		OFSET(1)=S0N/KNN0
     	ELSE
     		OFSET(1)=0.D0
     	END IF
		IF(DMEQ(IG) .EQ. 0.D0) THEN !-first time	
     		OFSET(2)=S0T/KSS0
		ELSE !-Update ofset(2)
			OFSET(2)=SIGN(S0T/KSS0,GAP0(2))
		END IF
c		WRITE(7,*) 'SON:', S0N
c     	WRITE(7,*) 'SOT:', S0T
c     	WRITE(7,*) 'OFSET1:', OFSET(1)
c     	WRITE(7,*) 'OFSET2:', OFSET(2)
c -- Cohesive Law at integration pnt I ----------------------	
		IF(CFSTAT(IG) .EQ. 0.D0) THEN !-first time
			GAP(1)=OFSET(1)
			GAP(2)=OFSET(2)
			KNN = KNN0
			TRACT(1,IG)=KNN*GAP(1)
			KSS = KSS0
			TRACT(2,IG)=KSS*GAP(2)
C			WRITE(7,*) 'KNN', KNN
C			WRITE(7,*) 'KSS', KSS
c			WRITE(7,*) 'TRACT(1)', TRACT(1,IG)
c			WRITE(7,*) 'TRACT(2)', TRACT(2,IG)	
C 			-- Effective displacement			
			U_EFF=SQRT(GAP(1)**2+GAP(2)**2)
c			WRITE(7,*) 'U_EFF:', U_EFF		
		! -- Effective stress
			T_EFF=SQRT(TRACT(1,IG)**2+TRACT(2,IG)**2)
c			WRITE(7,*) 'T_EFF:', T_EFF
		! -- Work of normal and shear stress
			GN=0.5*TRACT(1,IG)*GAP(1)
			GS=0.5*TRACT(2,IG)*GAP(2)
c			WRITE(7,*) 'GN:', GN
c			WRITE(7,*) 'GS:', GS		
		! -- Mixed-mode ratio (B-K)
			B=GS/(GN+GS)
c			WRITE(7,*) 'B:', B
		! -- Mixed-mode fracture energy
			GMC=GNC+(GSC-GNC)*B**BKETA
c			WRITE(7,*) 'GMC:', GMC
			U0_EFF(IG)=U_EFF
			UF_EFF(IG)=2*GMC/T_EFF+U0_EFF(IG)
			CFSTAT(IG)=1.D0
		ELSE
			GAP(1)=GAP0(1)+OFSET(1)
			GAP(2)=GAP0(2)+OFSET(2)
c			WRITE(7,*) 'GAP(1):', GAP(1)
c      		WRITE(7,*) 'GAP(2):', GAP(2)
C 			-- Effective displacement			
			U_EFF=SQRT(max(0.d0,GAP(1))**2+GAP(2)**2)
c			WRITE(7,*) 'U_EFF:', U_EFF	
		END IF
c
C		WRITE(7,*) 'DM:', DM(IG)
C		WRITE(7,*) 'IINC:', IINC
C		WRITE(7,*) 'IINC0:', IINC0
C		WRITE(7,*) 'DTIME:', DTIME
		IF((DM(IG) .LT. DFAIL).AND.((IINC .NE. IINC0)
     &	.OR. (KSTEP .NE. ISTEP0)).AND.(DTIME .GT. MDTIME)) THEN
		! -- Degradation: linear
			DM2=0.D0
			IF(UF_EFF(IG) .LE. U0_EFF(IG)) THEN
				DM2=DFAIL
			ELSE
				IF(U_EFF .NE. 0.D0) THEN
				WRITE(7,*) 'U_EFF:', U_EFF
				WRITE(7,*) 'U0_EFF:', U0_EFF(IG)
				WRITE(7,*) 'UF_EFF:', UF_EFF(IG)
			DM2=UF_EFF(IG)*(U_EFF-U0_EFF(IG))/(U_EFF*(UF_EFF(IG)-U0_EFF(IG)))
				WRITE(7,*) 'DMTEMP:', DM2
				DM2=MIN(DFAIL,DM2)
				DM2=MAX(0.D0,DM2)
     			END IF
			END IF
			WRITE(7,*) 'DM2:', DM2
			WRITE(7,*) 'DMEQ:', DMEQ(IG)
C
			IF(DM2 .GE. DMEQ(IG)) THEN
				DM(IG)=DM2
C				Update stresses				
				IF(GAP0(1) .GE. 0.D0) THEN
				KNN = KNN0*(1.D0-DM(IG))
				TRACT(1,IG)=KNN*GAP(1)
				ELSE
				KNN = KNNC
				TRACT(1,IG)=KNN*GAP0(1)
				END IF
				KSS = KSS0*(1.D0-DM(IG))
				TRACT(2,IG)=KSS*GAP(2)				
c				Separate DMU for two modes for unloading
				IF(GAP0(1) .GT. 0.D0) THEN				
				KNNU(IG)=MIN(KNNC,TRACT(1,IG)/GAP0(1))
				ELSE
				KNNU(IG)=KNNC
				END IF
				IF(GAP0(2) .NE. 0.D0) THEN
				KSSU(IG)=MIN(KNNC,TRACT(2,IG)/GAP0(2))
				ELSE
				KSSU(IG)=KNNC
				END IF
			ELSE !-Unloading
				DM(IG)=DMEQ(IG)
				KNNU(IG)=KNNUEQ(IG)
				KSSU(IG)=KSSUEQ(IG)
				IF(GAP0(1) .GT. 0.D0) THEN
				KNN = KNNU(IG)
				ELSE
				KNN = KNNC
				END IF
				KSS = KSSU(IG)
				TRACT(1,IG)=KNN*GAP0(1)
				TRACT(2,IG)=KSS*GAP0(2)
				WRITE(7,*) 'KNNUEQ', KNNUEQ(IG)
				WRITE(7,*) 'KSSUEQ', KSSUEQ(IG)
				WRITE(7,*) 'TRACT(1)', TRACT(1,IG)
				WRITE(7,*) 'TRACT(2)', TRACT(2,IG)
			END IF
		ELSE IF((IINC .EQ. IINC0) .AND. (KSTEP .EQ. ISTEP0)) THEN
			KNN = KNN0
			TRACT(1,IG)=KNN*OFSET(1)
			KSS = KSS0
			TRACT(2,IG)=KSS*OFSET(2)		
		ELSE ! DTIME<MDTIME or DM=DFAIL
C			DM(IG)=DMEQ(IG)
			KNNU(IG)=KNNUEQ(IG)
			KSSU(IG)=KSSUEQ(IG)
			IF(GAP0(1) .GT. 0.D0) THEN
				KNN = KNNU(IG)
			ELSE
				KNN = KNNC
			END IF
			KSS = KSSU(IG)
			TRACT(1,IG)=KNN*GAP0(1)
			TRACT(2,IG)=KSS*GAP0(2)	
		END IF
c-------------------------------------------------------------------------
c
c -- Calculate K matrix and RHS contribution of integration pnt I --------
C -- K = (FMATRX'*Q'*D*Q*FMATRX)*W(I), F=K*U=(FMATRX'*Q'*TRACT)*W(I) ----- 
		DO J=1,2
		DO K=1,2
			DMATRX(J,K)=0.D0
		END DO
		END DO
		DMATRX(1,1)=KNN
		DMATRX(2,2)=KSS
C
		CALL MATRIX_MUL(QMATRX,FMATRX,QF,2,2,8)
		CALL MATRIX_TRANSPOSE(QF,FTQT,2,8)
		CALL MATRIX_MUL(DMATRX,QF,DQF,2,2,8)
		CALL MATRIX_MUL(FTQT,DQF,FTQTDQF,8,2,8)
C		
       	DO J=1,8
          	DO K=1,8
            	KMATRX(J,K) = KMATRX(J,K)+FTQTDQF(J,K)*DET*W(IG)
          	end do
        end do		
C
		TRACT2(1)=TRACT(1,IG)
		TRACT2(2)=TRACT(2,IG)
		CALL MATRIX_MUL(FTQT,TRACT2,FTQTTAU,8,2,1)
c		
		DO J=1,8
			KU(J)=KU(J)+FTQTTAU(J)*W(IG)*DET
		END DO
C
C
       END DO !-Looped over all integration pnts
C
		
c		DO J=1,8
c		WRITE(7,*) 'KU3:', KU(J)
c		END DO
       RETURN
       END
C
C
C
C========================== SUBROUTINE INTPTS =================================
C
       SUBROUTINE KINTPT(PQ,W,KNODE,NINTP)
C
	  INCLUDE 'ABA_PARAM.INC'
C       
       DOUBLE PRECISION PQ(2),W(2),CN
C
C      Subroutine to initialize integration points PQ and weights W
C      NINTP = No. integration points (DEFAULT=2)
C      KNODE = No. of nodes on element (DEFAULT=4)
C
C      Integration points for 2D interface elements
C
C				4----3
C				|	 |
C				|	 |
C				1----2
c
	IF (NINTP .EQ. 2) THEN
		CN=0.5773502691896260D0
		PQ(1)=-CN
		PQ(2)=CN
		W(1)=1.D0
		W(2)=1.D0
	END IF
C
C
       RETURN
       END
C
C=========================== SUBROUTINE SHAPE ===============================
C
       SUBROUTINE KSHAPE(KNODE,PQ,F)
C
	  INCLUDE 'ABA_PARAM.INC'
C
       DOUBLE PRECISION F(4)
C
C  
	IF(KNODE .EQ. 4) THEN
		F(1)=.5D0*(1-PQ)
		F(2)=.5D0*(1+PQ)
		F(3)=F(2)
		F(4)=F(1)
	END IF
C
       RETURN
       END
C
C============================== SUBROUTINE DSHAPE ====================
C
       SUBROUTINE DSHAPE(KNODE,PQ,DF)
C
	  INCLUDE 'ABA_PARAM.INC'
C
C
       DOUBLE PRECISION DF(4)
C
C
	IF(KNODE .EQ. 4) THEN
		DF(1)=-.5D0
		DF(2)=.5D0
		DF(3)=DF(2)
		DF(4)=DF(1)
	END IF
C
       RETURN
       END
C
C========================== SUBROUTINE KUNITV ===================
C
       SUBROUTINE KUNITV(A,AMAG)
C
	  INCLUDE 'ABA_PARAM.INC'
C       
	DOUBLE PRECISION A(2), AMAG
C
C     Normalize vector A and return its magnitude as AMAG
C
      AMAG = DSQRT(A(1)*A(1)+A(2)*A(2))
	A(1) = A(1)/AMAG
	A(2) = A(2)/AMAG
C
       RETURN
	 END