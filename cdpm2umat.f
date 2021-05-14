c*****************************************************************
c
c Copyright (c) 2020 Peter Grassl
c
c Permission is hereby granted, free of charge, to any person obtaining a copy
c of this software and associated documentation files (the "Software"), to deal
c in the Software without restriction, including without limitation the rights
c to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
c copies of the Software, and to permit persons to whom the Software is
c furnished to do so, subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in all
c copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
c IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
c FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
c AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
c LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c SOFTWARE.
c**********************************************************************
c
c
      
      Subroutine umat50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,cma,qmat,elsizv,idelev)
c
c

c*****************************************************************
c Implementation of the stress evaluation for CDPM2 in UMAT50V in LS-DYNA
c Support page: http://petergrassl.com/Research/DamagePlasticity/CDPMLSDYNA/index.html
c
c Key references: 
c 1) P. Grassl, D. Xenos, U. Nyström, R. Rempling, K. Gylltoft.: "CDPM2: A damage-plasticity approach to modelling the failure of concrete". International Journal of Solids and Structures. Volume 50, Issue 24, pp. 3805-3816, 2013"
c 2) P. Grassl, U. Nyström, R. Rempling and K. Gylltoft, "A damage-plasticity model for the dynamic failure of concrete", 8th International Conference on Structural Dynamics, Leuven, Belgium, 2011.
c 3) P. Grassl and M. Jirasek: "Damage-plastic model for concrete failure". International Journal of Solids and Structures. Vol. 43, pp. 7166-7196, 2006."
c
c
c Authors of this subroutine: Dimitrios Xenos and Peter Grassl
c
c Compatible with LSDYNA Release 9.3.1
c
c*****************************************************************

      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer idelev(*)
      logical failels(*)
      character*5 etype
c
      real eps(6);
c       sig(nlq,6)   --- nominal stress in voigt notation xx,yy,zz,xy,yz,xz
c       epsps(nlq)   --- cummulative plastic strain (kappa)
c       hisv(nlq,35) - history variables
c       hisv(nlq,1) --- kappa
c       hisv(nlq,2) --- equivalent strain
c       hisv(nlq,3) --- plastic strain xx
c       hisv(nlq,4) --- plastic strain yy
c       hisv(nlq,5) --- plastic strain zz
c       hisv(nlq,6) --- plastic strain xy
c       hisv(nlq,7) --- plastic strain zy
c       hisv(nlq,8) --- plastic strain xz
c       hisv(nlq,9) --- kappa tension kdt
c       hisv(nlq,10) -- kappa tension 1 kdt1
c       hisv(nlq,11) -- kappa tension 2 kdt2
c       hisv(nlq,12) -- kappa compression kdc
c       hisv(nlq,13) -- kappa compression 1 kdc1
c       hisv(nlq,14) -- kappa compression 2 kdc2
c       hisv(nlq,15) -- damage variable tension omegaT (in LS-DYNA omegaT is written as wt)
c       hisv(nlq,16) -- damage variable tension omegaC (in LS-DYNA omegaC is written as wc)
c       hisv(nlq,17) -- strain rate factor used for the dynamic formulation based on quasi-static analysis  
c       hisv(nlq,18) -- alphac is the compression factor given in paper in IJSS by Grassl et al. in equation 46
c       hisv(nlq,19) -- equivalent strain tension eqstrT
c       hisv(nlq,20) -- equivalent strain compression eqstrC
c       hisv(nlq,21) -- total strain along xx
c       hisv(nlq,22) -- total strain along yy
c       hisv(nlq,23) -- total strain along zz
c       hisv(nlq,24) -- total strain along xy
c       hisv(nlq,25) -- total strain along yz
c       hisv(nlq,26) -- total strain along xz
c       hisv(nlq,27) -- equivalent strain (without rate factor influence)

c We assume that elen contains the length associated with one 
c integration point. This needs to be improved for triangular elements.
      real elen,trsh
      common/aux34loc/trsh(nlq),elen(nlq)
c       cm(1)   YM (Youngs modulus)
c	cm(2)	PR (Poissons ratio)
c	cm(3)	ECC (Eccentricity)
c	cm(4)	QH0 (Initial hardening)
c	cm(5)	FT (Uniaxial tension strength)
c	cm(6)	FC (Uniaxial compression strength)
c	cm(7)	HP (Hardening parameter)
c	cm(8)	AH (Hardening ductility measure)
c	cm(9)	BH (Hardening ductility measure)
c	cm(10)	CH (Hardening ductility measure)
c	cm(11)	DH (Hardening ductility measure)
c	cm(12)	AS (Damage ductility measure)
c	cm(13)	DF (Dilation constant)
c       cm(14)  ERT (Energy strain rate type)
c		     = 0.0: Rate does not affect fracture energy        
c                    = 1.0: Rate effect on fracture energy 
c                    = 2.0: Square of rate effect on fracture energy
c       cm(15)	TYPE (tensile damage type)
c		     = 0.0: Linear softening     
c                    = 1.0: Bi-Linear softening
c                    = 2.0: Exponential softening
c	cm(16)	BS (Damage: ductility parameter) 
c	cm(17)	WF (Damage: disp threshold 0)
c       cm(18)	WF1 (Damage: disp threshold 1)
c       cm(19)	FT1 (Damage: stress threshold 1)
c       cm(20)	SRT (Strength strain rate type)
c		     = 0.0: No rate effects        
c                    = 1.0: Model code 2010 first branch only 
c                    = 2.0: Model code 2010 first and second branch
c       cm(21)	failflag 
c		     = 0.0: if not ALL  gausspoints of the element are damaged        
c                    = 1.0: if ALL gausspoints of the element are damaged
c     cm(22)	efc (Damage Compression: Strain/displacement threshold )
c     cm(23)	damageflag (damage flag)
c     cm(24)	pflag (print flag)
c     bqs(nlq)  - pressure at each gauss point
      integer maxnip,nnm1,lft,llt
c     maxnip    ---- variable denoting max number of integration points
      integer mx,i,j
c
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c
c     ym   --------------- Young's modulus
c     pr   --------------- Poisson's ratio
c     ecc  --------------- eccentricity parameter used in the plasticity law. Its calibration is described in Jirasek & Bazant (2002)
c     qh0  --------------- value of the derivative of the 1st hardening function at kappa=0 as described in eq. (30) of the IJSS paper by Grassl et al.
c     ft   --------------- tensile strength
c     fc   --------------- compressive strength
c     hp   --------------- value of the derivative of the 2nd hardening function as described in eq. (31) of the IJSS paper by Grassl et al.
c     ah,bh,ch,dh -------- hardening parameters used in eq. (33) of the IJSS paper by Grassl et al.
c     as,bs,df ----------- softening  parameters used in eqs. (56) and (50) of the IJSS paper by Grassl et al.
c     eratetype ---------  parameter denoting whether and how strain rate dependence is taken into account for fracture energy
c     type --------------- softening type used in the formulation of the tensile damage variable
c     wf   --------------- max crack opening displacement used in the tensile damage law
c     wf1  --------------- max crack opening displacement used in the bilinear tensile damage law
c     ft1  --------------- tensile stress threshold used in the billinear tensile damage law
c     sratetype ---------  parameter denoting whether and how strain rate dependence is taken into account for strength
c     failflg ------------ flag denoting when an element should be deleted
c     efc  --------------- compressive strain/displacement threshold used as a parameter in the compressive damage law
c     m0   --------------- parameter used in the plasticity law calculated in eq.(20) of the IJSS paper by Grassl et al.
c     damageflag ------------ flag to denote whether isotropic or anisotropic damage law is used
c     printflag ---------- flag to print input only once
c
c Variables used for the evaluation of the plasticity algorithm
      real sig(6),sigEff(6),oldStrain(6),convStrain(6),
     $     deltaTotStrain(6),
     $     tempTotStrain(6),elStrain(6),princStress(3),
     $     princDir(3,3),totStrain(6),plastStrain(6),
     $     sum,tempTheta,
     $     strain(6), sigVTrial,rhoTrial,thetaTrial,apexStress,
     $     tempkappaP,yieldval
      integer subincCounter,rtype,subincFlag,converged,l,k
c       totstrain      -----------  total strain vector equal to the sum of  plastic and elastic strains
c       sigVTrial      ----------- trial volumetric stress 
c       rhoTrial       ----------- trial deviatoric stress
c       thetaTrial     ----------- trial Lode angle
c       apexStress     ----------- apexstress. Used only in the vertex case of the plasticity algorithm
c       yieldval   - value of the yield function
c       subincCounter   - counter of subincrementations performed in the plasticity algorithm
c       subincFlag   - flag denoting whether or not subincrementation is taking place
c                      =0 no subincrementation is taking place
c                      =1 subincrementation is taking place
c       rtype      - return type of the yield surface
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged       - integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     sigEff(6)        ---------- effective stress (no damage included)
c     oldStrain(6)     ---------- old elastic strain 
c     convStrain(6)    ---------- last elastic strain vector for which the plasticity algorithm has converged
c     deltaTotStrain(6)   ---------- difference between last and new strain vector
c     tempTotStrain(6)  ---------- temporary elastic strain
c     elStrain(6)      ----------  elastic strain
c     princStress(3)   ----------  array containing principal stresses
c     princDir(3,3)    ----------  matrix containing principal eigenvectors stored columnwise
c     plastStrain(6)   ----------  plastic strain vector
c     sum              ----------  variable used for summation
c     tempTheta        ---------- temporary Lode angle
c     strain(6)        ---------- strain array
c     l,k              ---------- counters used in various iterations
c
c Variables used in the damage algorithm
      real alpha,omegaT,omegaC,effStressT(6),effStressC(6),
     $     rateFactor,cdpm2u_computeRateFactor,strainrate(6),epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     length,deltaElStrain(6),stressOld(6), damagedGP    
c     alpha             --------- variable used to identify amount of contribution of compressive stresses and is defined in eq. (46) of the IJSS paper by Grassl et al.
c     omegaT            --------- tensile damage variable
c     omegaC            --------- compressive damage variable
c     effStressT        --------- effective tensile part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     effStressC        --------- effective compressive part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     rateFactor       ---------- variable used to incorporate impact effects on the constitutive law
c     cdpm2u_computeRateFactor --------- function used to calculate the rate factor
c     strainrate(6)    ---------- rate of the strain tensor used for the calculation of the rateFacto
c     epsilonT         ---------- tensile equivalent strain
c     epsilonC         ---------- compressive equivalent strain
c     kappaDT          ---------- history parameter kappaDT        
c     kappaDT1         ---------- history parameter kappaDT1       
c     kappaDT2         ---------- history parameter kappaDT2       
c     kappaDC          ---------- history parameter kappaDC        
c     kappaDC1         ---------- history parameter kappaDC1       
c     kappaDC2         ---------- history parameter kappaDC2       
c     length           ---------- characteristic length used for the formulation of the tensile damage function based on the crack-band approach
c     deltaElStrain(6) ---------- elastic strains of the previous step
c     damagedGP        ---------- number of damaged gausspoints within analysed element 
c
c ----------------------------------- Variable initialisation -----------------------------------------------
c
c      mx=48*(mxt(lft)-1)
c
c     Material constants
c
      ym=cm(1)
      pr=cm(2)
      ecc=cm(3)
      qh0=cm(4)
      ft=cm(5)
      fc=cm(6)
      hp=cm(7)
      ah=cm(8)
      bh=cm(9)
      ch=cm(10)
      dh=cm(11)
      as=cm(12)
      df=cm(13)
      eratetype=cm(14)
      type=cm(15)
      bs=cm(16)
      wf=cm(17)
      wf1=cm(18)      
      ft1=cm(19)
      sratetype=cm(20)
      failflg=cm(21)
      efc=cm(22)
      damageflag=cm(23)
      printflag=cm(24)
c This is the global tolerance. All other tolerances are made relative to this global tolerance
      gTol=1.e-6
      damagedGP=0.

c This function needs allows us to run the code without having to enter all 
c default parameters
      call cdpm2u_giveDefaultValuesIfNotGiven()
      
      cm(25) = printflag;
      
      m0 = 3. * ( fc** 2. - ft**2. ) / ( fc * ft ) * ecc / ( ecc + 1. );      
      do i=lft,llt
         tempkappaP=hsvs(i,1)
c Write strain increment vector
         eps(1) = d1(i)
         eps(2) = d2(i)
         eps(3) = d3(i)
         eps(4) = d4(i)
         eps(5) = d5(i)
         eps(6) = d6(i)

         do l=1,6
            totStrain(l)   = eps(l) + hsvs(i,l+20)
            strainrate(l)  = eps(l)/dt1siz(i)
            plastStrain(l) =  hsvs(i,l+2)
            oldStrain(l)  = hsvs(i,l+20)
            convStrain(l) = oldStrain(l)
            tempTotStrain(l) = totStrain(l)
            deltaTotStrain(l) = eps(l)
         enddo

         subincCounter=0
         subincFlag=0
         
         converged=1
                 
         do while ( converged .eq. 1 .or. subincFlag .eq. 1 ) 
            
            do l=1,6
               elStrain(l) = tempTotstrain(l) - plastStrain(l)              
            enddo   

            call cdpm2u_computeStressesfromStrains(sigEff,elStrain,
     $           ym,pr)
            
            call cdpm2u_computeTrialCoordinates(sigEff,sigVTrial,
     $           rhoTrial,tempTheta)
            thetaTrial=tempTheta
            call cdpm2u_computeYieldValue(yieldval,sigVTrial,rhoTrial,
     $           thetaTrial,tempKappaP)

            apexStress=0.

            if (yieldval .gt. 0.) then
               call cdpm2u_checkForVertexCase(apexStress,sigVTrial,
     $              tempKappaP,rtype)
               if (rtype.eq.1 .or. rtype .eq. 2) then
                  call cdpm2u_performVertexReturn(sigEff,
     $                 apexStress,tempKappaP,rtype,converged)
               end if
               
               if (rtype.eq.0) then
                  call cdpm2u_performRegularReturn(sigEff, 
     $                 tempKappaP,converged,ym,pr,gTol)
               end if
            else
               converged=0
               do l=1,6        
                  plastStrain(l)=hsvs(i,l+2)
               enddo
               goto 925
            end if                     
            
            if ( converged .eq. 1 ) then
               subincCounter=subincCounter+1
               if ( subincCounter .gt. 10 ) then
                  write(*,*) '*** Perform Plasticity return with' 
                  write(*,*) 'subincrementation methodology'
                  write(*,*) 'No convergence reached !***'
                  stop
               else if (subincCounter .gt. 9 .and. 
     $                 tempKappaP .lt. 1.0 ) then
                  tempKappaP=1.
               end if
               subIncFlag = 1
               do l=1,6
                  deltaTotStrain(l)=deltaTotStrain(l)*0.5
                  tempTotStrain(l)=convStrain(l)+deltaTotStrain(l)
               enddo
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 0) then
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $          ym,pr)
               do l=1,6        
                  plastStrain(l)=totStrain(l)-elStrain(l)
               enddo
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 1) then
c               write(*,*) '*** Subincrementation required',subincCounter               
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $              ym,pr)
               do l=1,6
                  plastStrain(l)=tempTotStrain(l)-elStrain(l)
                  convStrain(l)=tempTotStrain(l)
                  deltaTotStrain(l)=totStrain(l)-convStrain(l)
                  tempTotStrain(l)=totStrain(l)
               enddo
               subincCounter = 0
               subincFlag=0
               converged=1            
            end if
         end do
         
                     
 925     continue
         if (damageflag .eq. 3.0) then
            omegaT=0.0
            omegaC=0.0            
            epsilonT=0.0
            kappaDT=0.0
            kappaDT1=0.0
            kappaDT2=0.0
            kappaDC=0.0
            kappaDC1=0.0
            kappaDC2=0.0
            omegaT=0.0
            omegaC=0.0
            rateFactor=0.0
            alpha=0.0
            epsilonT=0.0
            epsilonC=0.0
            do l=1,6
c               sig(i,l)=sigEff(l)
               sig(l)=sigEff(l)
            enddo            
            goto 152
         end if

c Initialize parameters used in the damage algorithm              
         rateFactor=hsvs(i,17)

         epsilonT=hsvs(i,19)
         epsilonC=hsvs(i,20)
         kappaDT=hsvs(i,9)
         kappaDT1=hsvs(i,10)
         kappaDT2=hsvs(i,11)
         kappaDC=hsvs(i,12)
         kappaDC1=hsvs(i,13)
         kappaDC2=hsvs(i,14)
         omegaT=hsvs(i,15)
         omegaC=hsvs(i,16)
         alpha=hsvs(i,18)
         call cdpm2u_computeAlpha(effStressT,effStressC,sigEff,alpha)
         length=elen(i)
         sum=0.
         do l=1,6
c     Compute norm of increment of plastic strains       
            sum=sum+(plastStrain(l)-hsvs(i,l+2))**2
            deltaElStrain(l)=(hsvs(i,l+20)-hsvs(i,l+2))
         enddo
         call cdpm2u_computeStressesfromStrains(stressOld,deltaElStrain,
     $        ym,pr)
         sum=sqrt(sum)
         rateFactor=hsvs(i,17)
         call cdpm2u_computeDamage(omegaC,omegaT,strainrate,
     $        rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,
     $        kappaDT2,kappaDC,kappaDC1,kappaDC2,sigEff,sum,
     $        tempKappaP,length,stressOld,hsvs(i,18),hsvs(i,27))
         
         do l=1,6
            if (damageflag .eq. 0.0) then
               sig(l)=(1.-omegaT)*effStressT(l)+
     $              (1.-omegaC)*effStressC(l)               
            else if (damageflag .eq. 1.0) then 
                sig(l)=(1.-omegaT)*sigEff(l)
           else if (damageflag .eq. 2.0) then
               sig(l)=1.-(1.-omegaT*(1-alpha))*
     $              (1-omegaC*alpha)*sigEff(l)
            end if
         enddo

 152     continue
         
c     Write the history variable at the end of the routine
         epsps(i)= tempkappaP 
         hsvs(i,1)= tempkappaP 
         hsvs(i,2)= epsilonT

         hsvs(i,3)= plastStrain(1)
         hsvs(i,4)= plastStrain(2)
         hsvs(i,5)= plastStrain(3)
         hsvs(i,6)= plastStrain(4)
         hsvs(i,7)= plastStrain(5)
         hsvs(i,8)= plastStrain(6)
         
         hsvs(i,9)=  kappaDT
         hsvs(i,10)= kappaDT1
         hsvs(i,11)= kappaDT2
         hsvs(i,12)= kappaDC
         hsvs(i,13)= kappaDC1
         hsvs(i,14)= kappaDC2
         hsvs(i,15)= omegaT
         hsvs(i,16)= omegaC
         hsvs(i,17)=rateFactor
         hsvs(i,18)=alpha
         hsvs(i,19)=epsilonT
         hsvs(i,20)=epsilonC

         hsvs(i,21)=totStrain(1)
         hsvs(i,22)=totStrain(2)
         hsvs(i,23)=totStrain(3)
         hsvs(i,24)=totStrain(4)
         hsvs(i,25)=totStrain(5)
         hsvs(i,26)=totStrain(6)

         if ( hsvs(i,15).gt. 0.9995 .and. hsvs(i,16).gt. 0.9995 ) then
            damagedGP=damagedGP+1.0
         end if

         if (failflg .gt. 0.0) then
            damagedGP=damagedGP/FLOAT(maxnip)
            if (damagedGP .ge. failflg) then
               failur=.true.
                  failels(i)=1
                  do l=1,6
                     sig(l)=0.0
                  enddo
            end if
         end if
         
c     Write sig components
         sig1(i) = sig(1)
         sig2(i) = sig(2)
         sig3(i) = sig(3)
         sig4(i) = sig(4)
         sig5(i) = sig(5)
         sig6(i) = sig(6)
      enddo
      return
      end
      
      subroutine cdpm2u_giveDefaultValuesIfNotGiven()
c     Subroutine to check if all model parameters have values. If a parameter does not have any values a default value is provided.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c
      real epsilon
      if( pr .lt. 0) then
         pr = 0.2
      end if
      if (wf1 .le. 0.) then
         wf1= 0.15*wf
      end if
      if (ft1 .le. 0.) then
         ft1= 0.3*ft
      end if
      if (qh0 .le. 0.) then
         qh0=0.3
      end if
      if (hp .lt. 0.) then
         hp=0.5
      end if
      if (ah .le. 0.) then
         ah=8.e-2
      end if
      if (bh .le. 0.) then
         bh=3.e-3
      end if
      if (ch .le. 0.) then
         ch=2.
      end if
      if (dh .le. 0.) then
         dh=1.e-6
      end if
      if (df .le. 0.) then
         df=0.85
      end if
      if (eratetype .le. 0.) then
         eratetype=0.
      end if
      if (efc .le. 0. ) then
         efc=100.e-6
      end if
      if (as .le. 0.) then
         as=15.
      end if
      if (bs .le. 0.) then
         bs=1.
      end if
      if (sratetype .le. 0.) then
         sratetype=0.
      end if
      if (damageflag .le. 0.) then
         damageflag=0.
      end if
      if (gTol .le. 0.) then
         gTol=1.e-6
      end if
      if (ecc .le. 0.) then
         epsilon=ft*((1.16*fc)**2.-fc**2.)/(1.16*fc*(fc**2.-ft**2.))
         ecc=(1.+epsilon)/(2.-epsilon)
      end if

c Write output
c     Print all input variables
c
      if (printflag .lt. 0.5) then
      write(*,*) '*****************************************************'
      write(*,*) '*****************************************************'
      write(*,*) '  LS-DYNA is using the usermaterial CDPM2 (UMAT50V)  '
      write(*,*) '  -------------------------------------------------  '
      write(*,*) '                                                     '
      write(*,*) '  Material parameters:                               '
      write(*,2) '   YM (Youngs modulus)............... = ',ym
      write(*,2) '   PR (Poissons ratio)............... = ',pr
      write(*,2) '   ECC (Eccentricity)................ = ',ecc
      write(*,2) '   QH0 (Initial hardening)........... = ',qh0
      write(*,2) '   FT (Uniaxial tension strength).... = ',ft
      write(*,2) '   FC (Uniaxial compression strength) = ',fc
      write(*,2) '   HP (Hardening modulus)............ = ',hp
      write(*,2) '   AH (Hardening ductility measure).. = ',ah
      write(*,2) '   BH (Hardening ductility measure).. = ',bh
      write(*,2) '   CH (Hardening ductility measure).. = ',ch
      write(*,2) '   DH (Hardening ductility measure).. = ',dh
      write(*,2) '   AS (Damage ductility measure)..... = ',as
      write(*,2) '   DF (Dilation constant)............ = ',df
      write(*,2) '   ERATETYPE (Strength rate type).... = ',eratetype     
      write(*,*) '     EQ.0.0: constant fracture energy            '
      write(*,*) '     EQ.1.0: DIF times fracture energy'
      write(*,*) '     EQ.2.0: square of DIF times fracture energy'
      write(*,2) '   TYPE (Damage type)................ = ',type
      write(*,*) '     EQ.0.0: Linear softening           '
      write(*,*) '     EQ.1.0: Bi-Linear softening        '
      write(*,*) '     EQ.2.0: Exponential softening      '
      write(*,*) '     EQ.4.0: No damage                  '
      write(*,2) '   BS (Damage: ductility parameter)   = ',bs
      write(*,2) '   WF (Damage: disp threshold 0)..... = ',wf
      write(*,2) '   WF1 (Damage: disp threshold 1).... = ',wf1
      write(*,2) '   FT1 (Damage: stress threshold 1).. = ',ft1
      write(*,2) '   EFC (strain threshold in comp).... = ',efc
      write(*,2) '   SRATETYPE (Strength rate type).... = ',sratetype     
      write(*,*) '     EQ.0.0: No rate effects            '
      write(*,*) '     EQ.1.0: MC 2010 first branch'
      write(*,*) '     EQ.2.0: MC 2010 first and second branch'
      write(*,*) '   DAMAGEFLAG (damage flag)......... =  ',damageflag
      write(*,*) '     EQ.0.0: standard model with two damage params'
      write(*,*) '     EQ.1.0: iso model with one damage param'
      write(*,*) '     EQ.2.0: standard model with one damage param'
      write(*,*) '     EQ.3.0: no damage'
      write(*,*) '*****************************************************'
      write(*,*) '  History variables:                                 '
      write(*,*) '   Var  #1 = kappa'
      write(*,*) '   Var  #2 = equivalient strain'
      write(*,*) '   Var  #3 = plastic strain direction 1'
      write(*,*) '   Var  #4 = plastic strain dierction 2'
      write(*,*) '   Var  #5 = plastic strain direction 3'
      write(*,*) '   Var  #6 = hardening tension kdt'
      write(*,*) '   Var  #7 = hardening tension kdt1'
      write(*,*) '   Var  #8 = hardening tension kdt2'
      write(*,*) '   Var  #9 = hardening compression kdc'
      write(*,*) '   Var #10 = hardening compression kdc1'
      write(*,*) '   Var #11 = hardening compression kdc2'
      write(*,*) '   Var #12 = damage function tension wt'
      write(*,*) '   Var #13 = damage function compression wc'
      write(*,*) '   Var #14 = (internal book keeping)'
      write(*,*) '   Var #15 = compression factor alphac'
      write(*,*) '   Var #16 = ratefactor alphar'
      write(*,*) '   Var #17 = elastic strain direction 1'
      write(*,*) '   Var #18 = elastic strain direction 2'
      write(*,*) '   Var #19 = elastic strain direction 3'
      write(*,*) '   Var #20 = equivalent strain tension'
      write(*,*) '   Var #21 = equivalent strain compression'
      write(*,*) '   Var #22-#27 = undamaged stresses'
      write(*,*) '*****************************************************'
      write(*,*) '  Check that the bulk modulus and shear modulus are  '
      write(*,*) '  correctly calculated. They MUST be set on the      '
      write(*,*) '  input card.                                        '
      write(*,1) '  With E = ',ym,' and pr = ',pr
      write(*,2) '  BLK = E/(3*(1-2pr)) = ',ym/(3.0*(1.-2.*pr))
      write(*,2) '  SHR = E/(2*(1+pr))  = ',ym/(2.*(1.+pr))
      write(*,*) '*****************************************************'
      write(*,*) '*****************************************************'
 1        format(1x,A,1pE9.3,A,1pE9.3)
 2            format(1x,A,1pE12.5)
      printflag=1.0
      end if
c
      return
      end    


      subroutine cdpm2u_computeStrainsfromStresses(stress,strain,ym,pr)
c     Subroutine to calculate elastic strains from the stress tensor. Performs operation epsilon = D : sigma
      real stress(6),strain(6),ym,pr
c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      strain(1)=(stress(1) - pr * stress(2) - pr * stress(3))/ym
      strain(2)=(-pr*stress(1) + stress(2) - pr * stress(3))/ym
      strain(3)=(-pr*stress(1) - pr * stress(2) + stress(3))/ym
      strain(4)=(2. * (1+pr) * stress(4) )/ym
      strain(5)=(2. * (1+pr) * stress(5) )/ym
      strain(6)=(2. * (1+pr) * stress(6) )/ym
      return
      end
      
      subroutine cdpm2u_computeStressesfromStrains(stress,strain,ym,pr)
c     Subroutine to calculate strains from the elastic strain tensor. Performs operation sigma = C : epsilon

c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      real factor,  stress(6),strain(6),ym,pr
      factor = ym/((1.+pr)*(1.-2.*pr))
      stress(1)=factor*((1.-pr)*strain(1) + pr * strain(2) + 
     $     pr * strain(3))
      stress(2)=factor*(pr*strain(1) + (1.-pr) * strain(2) + 
     $     pr * strain(3))
      stress(3)=factor*(pr*strain(1) + pr * strain(2) + 
     $     (1.-pr)*strain(3))
      stress(4)=factor*(((1.-2.*pr)/2.) * strain(4) )
      stress(5)=factor*(((1.-2.*pr)/2.) * strain(5) )
      stress(6)=factor*(((1.-2.*pr)/2.) * strain(6) )
      return
      end
c ---------------------------------------------------------VERTEX RETURN FUNCTIONS ---------------------------------------------------

      subroutine cdpm2u_checkForVertexCase(apexStress,sigV,tempkappa,
     $     rtype)
c     Subroutine that check whether the current stress state requires plasticity return to the vertex, at which
c     derivative of plastic potential and yield surface are discontinuous.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c
      real apexStress,sigV,tempkappa,qh2,cdpm2u_qh2fun 
      integer rtype
c     sigV              -----------  volumetric stress <Input>
c     tempKappa         -----------  cummulative plastic strain  <Input>
c     apexStress        -----------  sigmaV of the yield surface for the current 
c                                    tempkappa and rho=theta=0 <Output>
c       rtype           -----------  return type of the yield surface  <Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     qh2               -----------  variables containing the results of the hardening functions 
c     cdpm2u_qh2fun     -----------  function to calculate the hardening function  given in eq. (31) of IJSS paper by P. Grassl et al.

      if ( sigV .gt. 0. ) then
         rtype = 1
         if (tempKappa .lt. 1.) then
            apexStress = 0.
         else
            qh2=cdpm2u_qh2fun(tempKappa,hp)
            apexStress=qh2*fc/m0
         end if        
      else if ( sigV .lt. 0. .and. tempKappa .lt. 1.) then
         rtype = 2
         apexStress = 0.
      else 
         rtype = 0
         apexStress=0.
      end if
      
      return 
      end
      
      
      subroutine cdpm2u_performVertexReturn(stress,apexStress,
     $     kappa,rtype,converged)
c     Subroutine that performs plasticity return close whenever the stress state needs
c     to be returned to the apex. If the stress state is not an actual vertex case
c     rtype=0 is returned.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      real sigV,rho,apexStress,kappa,stress(6),theta,
     $     yieldValue,yieldValueMid,sig2,dSig,sigMid,sigAnswer,
     $     ratioPotential,kappa0,tempKappaP, cdpm2u_computeTempKappa,
     $     cdpm2u_computeRatioPotential,ratioTrial,yldTol
      integer i,j,k,rtype,maxiter,converged
c     apexStress     ------------------ variable containing the apex to return <Input>
c     kappa          ------------------ cummulative plastic strain  <Input>
c     rtype          ------------------  return type of the yield surface  <Input/Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged      ------------------ integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     kappa0,tempKappa---------------- variables for cummulative plastic strains
c     yieldValue,yieldValueMid ------- variables with yield values 
c     sigV           ----------------- volumetric stress 
c     rho            ----------------- deviatoric stress 
c     sig2,sigMid,sigAnswer,dSig ----- variables containing volumetric stress units
c     ratioPotential ----------------- variable containing the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     ratioTrial     ----------------- ratio of rho/sigmaV
c     yldTol         ----------------- tolerance used in bisection method solver
c     cdpm2u_computeTempKappa--------- function to calculate tempKappa according to eq.(32) of the IJSS paper by Grassl et al.
c     cdpm2u_computeRatioPotential --- function to calculate the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     maxiter         ---------------- parameter denoting max number of iterations performed using the bisection method
      yldTol=gTol
      yieldValue = 0.
      yieldValueMid = 0.
      sig2 = 0.
      kappa0=kappa
      tempKappaP=kappa
      maxiter=250

      call cdpm2u_computeTrialCoordinates(stress,sigV,rho,theta)

      sig2 = apexStress
      
      tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigV,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValue,sigV, 0., 0., tempKappaP)
      
      tempKappaP =
     $     cdpm2u_computeTempKappa(kappa0, sigV, rho, sig2,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValueMid,sig2, 0.,0.,
     $     tempKappaP)
      
      if ( yieldValue * yieldValueMid .ge. 0. )  then
         converged=1
         rtype = 0  
         goto 501
      end if
      
      if ( yieldValue .lt. 0.0 ) then
         dSig = sig2 - sigV
         sigAnswer = sig2
      else 
         dSig = sigV - sig2
         sigAnswer = sig2
      end if
      
      do  j = 1, maxiter
         dSig = 0.5 * dSig
         sigMid = sigAnswer + dSig
         tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigMid,
     $        ym,pr)
         
         call cdpm2u_computeYieldValue(yieldValueMid,sigMid, 0., 0.,
     $        tempKappaP)
         
        if ( yieldValueMid .le. 0. ) then
            sigAnswer = sigMid
         end if
         if (abs(yieldValueMid) .lt. yldTol .and. 
     $        yieldValueMid .le. 0.) then

            ratioPotential =
     $           cdpm2u_computeRatioPotential(sigAnswer, tempKappaP)
            
            ratioTrial = rho / ( sigV - sigAnswer );
            
            if ( ( ( ( ratioPotential .ge. ratioTrial ) .and. 
     $           rtype .eq. 1 ) ) .or.
     $           ( ( ratioPotential .le. ratioTrial ) .and. 
     $           rtype .eq. 2  ) ) then
               goto 500
            else    
               converged=1
               rtype = 0           
               goto 501
            endif
         endif
      enddo
 500  do k = 1, 3
         stress(k) = sigAnswer
         stress(k+3) = 0.
      enddo
      kappa=tempKappaP                          
      converged=0
 501  continue
      return
      end


      real function cdpm2u_computeTempKappa(kappaInitial,sigV1,rho,
     $     sigV2,ym,pr)
c     Function to calculate the tempKappa whenever requested from the performVertexReturn function.
c     TempKappa is calculated according to eq. (32) of the IJSS paper by P. Grassl et al.

      real sigV1,rho,sigV2,kappaInitial,ym,pr,
     $     equivalentDeltaPlasticStrain,kM,gM,ducMeas, 
     $     cdpm2u_computeDucMeas
c     sigV1 -------------- volumetric stress in the previous stress state <Input>
c     sigV2 -------------- volumetric stress in the current stress state <Input>
c     rho   -------------- deviatoric stress  <Input>
c     kappaInitial ------- previous kappaP (cummulative plastic strain) <Input>
c     ym    -------------- Young's modulus  <Input>
c     pr    -------------- Poisson's ratio <Input>
c     kM,gM -------------- bulk and shear moduli
c     equivalentDeltaPlasticStrain  -----  Increase of the plastic strains
c     ducMeas------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas------ function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM = ym / ( 2. * ( 1. + pr ) )

      equivalentDeltaPlasticStrain = sqrt( 1. / 9. *  (( sigV1 - sigV2 ) 
     $     /  kM )** 2.  + (rho / ( 2. * gM ))** 2. )
                                         
      ducMeas = cdpm2u_computeDucMeas(sigV2, 0., 3.141592653589793/3.)

      cdpm2u_computeTempKappa=kappaInitial+equivalentDeltaPlasticStrain/
     $     ducMeas
      return
      end

      real function cdpm2u_computeRatioPotential(sig ,kappa)
c     Function to calculate the ratio of the derivatives of the plastic potential, given in eq.(22) of the IJSS paper by P. Grassl et al. with respect to the deviatoric and volumetric stress respectively.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      real AGParam,BGParam,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,R,mQ,kappa,
     $     sig,rho,
     $     dgdsig,dgdrho,Al,Bl
      integer j
c     sig          ---------------- volumetric stress <Input>
c     rho          ---------------- deviatoric stress <Input>
c     kappa        ---------------- cummulative plastic strain kappaP <Input>
c     dgdsig,dgdrho --------------- derivatives of the plastic potential
c     AGParam,BGParam ------------- components of the plastic potential function
c     qh1,qh2          ------------ variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun - functions to calculate the hardening functions  given in eqs. (30) in IJSS paper by Grassl et al.
c     Al,Bl        ---------------- components of the function given in eq.(23) of the IJSS paper by P. Grassl et al.
c     R,mQ         ---------------- variables 
      rho=0.
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =    qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2.*df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(1.5) *rho/fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * (1. - qh1 ) * Bl + 6. ) +
     $     m0 * (qh1**2.) / ( sqrt(6.) * fc )
      
      cdpm2u_computeRatioPotential= dgdrho/dgdsig*3.*(1.-2.*pr ) / 
     $     ( 1. + pr )
      return 
      end

c ------------------------------------ Plasticity Algorithm functions (regular return functions) ---------------------------------
      subroutine cdpm2u_performRegularReturn(stress,kappa,converged,ym,
     $     pr,yieldTol)
c     Subroutine to perform regular plasticity return
      
      real stress(6),tempKappa,ym,pr, resid(4), normalisedResid(4),PI,
     $     jacobian(4,4), inverseJac(4,4),increment(4),
     $     deltaLambda,normOfResiduals, unknowns(4),
     $     trialSig,trialRho,trialTheta,kappa,kappaP,tempKappaP,
     $     yieldTol,sig,rho,ddkappadDeltaLambdadInv(2),
     $     dgdInv(2),ddgddInv(2,2), dkappadDeltaLambda,
     $     dfdkappa, ddgdInvDkappa(2), ddkappadDeltaLambdadKappa,
     $     ddkappaddDeltaLambdadInv(2),kM,gM,dfdInv(2),sum,
     $     stressPrincipal(3),princDir(3,3)
      integer i,j,iterations,totIter,converged,error
      
c     stress(6)      ------------------  effective stress components are in order xx,yy,zz,xy,yz,xz  <Input/Output>
c     kappa          ------------------  cummulative plastic strain (kappa) <Input>
c     ym             ------------------  Young's modulus  <Input>
c     pr             ------------------  Poisson's ratio <Input>
c     yieldTol       ------------------ tolerance of the N-R solver <Input>
c     converged      ------------------ integer showing whether solution has converged <Output>
c                                        = 0 solution has converged
c                                        = 1 solution has not converged
c     resid(4)       ------------------ array containing the residuals of the N-R solver
c     normalisedResid(4) -------------- array containing the normalised residuals of the N-R solver
c     PI             ------------------ parameter defined as the pi=3.14....
c     jacobian(4,4)  ------------------ matrix containing the jacobian of the problem in order to calculate the solution
c     inversJac(4,4) ------------------ matrix containing the inverse matrix of the jacobian of the problem in order to calculate the solution
c     increment(4)   ------------------ array containing the increments of the unknowns at each N-R iteration
c     deltaLambda    ------------------ plastic multiplier
c     normOfResiduals ----------------- norm the array resid(4)
c     unknowns(4)    ------------------ array with the unknowns in order sigV,rho,kappa,deltaLambda
c     trialSig       ------------------ trial volumetric stress (initial guess)
c     trialRho       ------------------ trial deviatoric stress (initial guess)
c     trialTheta     ------------------ trial Lode angle (initial guess)
c     kappaP          ------------------ cummulative plastic strain kappa (initial guess)
c     tempKappa      ------------------ temporary cummulative plastic strain kappa (each iteration)
c     sig            ------------------ temporary volumetric stress (each iteration)
c     rho            ------------------ temporary deviatoric stress (each iteration)
c     ddkappadDeltaLambdadInv(2) ------ derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric stress
c     dfdInv(2)      ------------------ derivative of the yield function with respect to the volumetric and deviatoric stress
c     dgdInv(2)      ------------------ derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     ddgddInv(2,2)  ------------------ second derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     dkappadDeltaLambda -------------- derivative of kappa with respect to the plastic multiplier
c     dfdkappa       ------------------ derivative of the yield function with respect to kappa
c     ddgdInvDkappa(2) ---------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and to kappa
c     ddkappadDeltaLambdadKappa ------- derivative of kappa with respect to the plastic multiplier and kappa
c     ddkappaddDeltaLambdadInv(2) ----- derivative of kappa with respect to the plastic multiplier and deviatoric and volumetric stress
c     kM             ------------------ bulk modulus
c     gM             ------------------ shear modulus
c     sum            ------------------ variable used for summation
c     stressPrincipal(3)  ------------- array containing the principal stresses 
c     princDir(3,3)  ------------------ matrix containing eigenvectors of the effective stress tensor stored columnwise 
c     i,j,iterations ------------------ integers used as counters
c     totIter        ------------------ maximum number of iterations of the N-R algorithm
c     error          ------------------ integer indicating whether the inversion of the jacobian matrix was successful
      iterations=0
      totIter=100

      PI=3.1415926535897932384626433832795029
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM =  ym / ( 2. * ( 1. + pr ) )
      normOfResiduals=1.
      do i=1,4
         resid(i)=0.
         normalisedResid(i)=0.
         unknowns(i)=0.
         increment(i)=0.
      enddo
      deltaLambda=0.
      call cdpm2u_computePrincValues(stress,stressPrincipal,0,princDir)
      call cdpm2u_computeTrialCoordinates(stress,trialSig,trialRho,
     $     trialTheta)     
      kappaP=kappa
      tempKappaP=kappa
      sig=trialSig
      rho=trialRho
      unknowns(1)=trialSig
      unknowns(2)=trialRho
      unknowns(3)=tempKappaP
      unknowns(4)=0.
      
      call cdpm2u_computeYieldValue(resid(4),sig, rho,trialTheta,
     $     tempKappaP)
      normOfResiduals=1.
      do while (normOfResiduals .gt. yieldTol)
c      write(*,*) 'normOfResiduals = ',normOfResiduals
         iterations=iterations+1
         if (iterations .eq. totIter) then
            converged=1
            goto 600
         end if 
         normalisedResid(1)=resid(1)/kM
         normalisedResid(2)=resid(2)/2./gM
         normalisedResid(3)=resid(3)
         normalisedResid(4)=resid(4)
         normOfResiduals=sqrt(normalisedResid(1)**2.+
     $        normalisedResid(2)**2.+normalisedResid(3)**2. +
     $        normalisedResid(4)**2.)
         
         if (isnan(normOfResiduals)) then
            converged=1
            goto 600
         end if

         if (normOfResiduals .gt. yieldTol) then
c     ----------------- compute jacobian ---------------------------------
           call cdpm2u_computedfdInv(dfdInv,sig,rho,trialTheta,
     $           tempKappaP) 
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgddInv(ddgddInv,sig,rho,trialTheta,
     $          tempKappaP)
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP)
           call cdpm2u_computedfdKappa(dfdkappa,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $          trialTheta,tempKappaP)
           call cdpm2u_computeddKappadDeltaLambdadKappa(
     $          ddkappadDeltaLambdadKappa,sig,rho,tempKappaP,trialTheta)

           call cdpm2u_computeddKappadDeltaLambdadInv(
     $          ddKappaddDeltaLambdadInv,sig,rho,tempKappaP,trialTheta)

           jacobian(1,1) = 1. + kM * deltaLambda *   ddgddInv(1, 1)
           jacobian(1, 2) = kM * deltaLambda * ddgddInv(1, 2)
           jacobian(1, 3) = kM * deltaLambda * ddgdInvdKappa(1)
           jacobian(1, 4) = kM * dgdInv(1)
           
           jacobian(2, 1) = 2. *gM *deltaLambda *ddgddInv(2, 1)
           jacobian(2, 2) = 1. + 2. *gM *deltaLambda *  ddgddInv(2, 2)
           jacobian(2, 3) = 2. *gM *deltaLambda * ddgdInvdKappa(2)
           jacobian(2, 4) = 2. *gM *dgdInv(2)
           
           jacobian(3, 1) = deltaLambda * ddKappaddDeltaLambdadInv(1)
           jacobian(3, 2) = deltaLambda * ddKappaddDeltaLambdadInv(2)
           jacobian(3, 3) = deltaLambda * ddkappadDeltaLambdadKappa - 1.
           jacobian(3, 4) = dkappadDeltaLambda
           
           jacobian(4, 1) = dfdInv(1)
           jacobian(4, 2) = dfdInv(2)
           jacobian(4, 3) = dfdKappa
           jacobian(4, 4) = 0.

           call cdpm2u_computeInverseJac(inverseJac,jacobian,error)
           if (error.eq.-1) then
              converged=1
              goto 600
           end if
           do i=1,4
              sum = 0.
              do j = 1,4
                 sum =sum+ inverseJac(i, j) * resid(j)
              enddo
              increment(i) =0.-sum
              unknowns(i)=unknowns(i)+increment(i)
           enddo
           if (unknowns(4) .le. 0.) then
              unknowns(4)=0.
           end if
           if (unknowns(2) .le. 0.) then
              unknowns(2)=0.
           end if
           if (unknowns(3)-kappaP .le. 0.) then
              unknowns(3)=kappaP
           end if
           sig = unknowns(1)
           rho = unknowns(2)
           tempKappaP=unknowns(3)
           deltaLambda=unknowns(4)
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP)
           resid(1) = sig - trialSig + kM *deltaLambda * dgdInv(1)
           resid(2) = rho - trialRho +  2.* gM *deltaLambda * dgdInv(2)
           resid(3) = -tempKappaP +kappaP+deltaLambda*dkappadDeltaLambda
           call cdpm2u_computeYieldValue(resid(4),sig,rho,trialTheta,
     $          tempKappaP)
        end if
      end do
      converged=0
      
      
      stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(trialTheta)
      stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta - 2. * PI/ 3.)
      stressPrincipal(3) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta + 2. * PI / 3.)

      call cdpm2u_transformStressVectorTo(stress,princDir,
     $     stressPrincipal)
      
      kappa=tempKappaP
 600  continue
      return
      end
      
c ------------------------------------ Plasticity Algorithm functions (general functions) ----------------------------------------
      subroutine cdpm2u_computedfdInv(dfdInv,sig,rho,theta,kappa) 
c     Subroutine to calculate the derivative of the yield function with respect to the volumetric and deviatoric stresses respectively 
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real rFunction,qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     dfdsig,dfdrho,sig,dfdInv(2),Bl
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdInv(2)        -------------  derivatives of yield function with respect to volumeric and deviatoric stress <Output>
c     dfdsig           -------------  derivative of the yield function with respect to volumetric stress
c     dfdrho           -------------  derivative of the yield function with respect to deviatoric stress
c     Al,Bl            -------------  variables corresponding to components of the yield function
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.       
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      
      
      Al =( 1. - qh1 ) * ( sig / fc + rho / ( sqrt(6.) * fc ) )** 2. +
     $     sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      
      dfdsig= 4. * ( 1. - qh1 ) / fc * Al * Bl + qh2* qh1**2. * m0 / fc
      dfdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6. )+ 
     $     rFunction * m0 * qh2 * qh1**2./ ( sqrt(6.) * fc )
      
      dfdInv(1) = dfdsig
      dfdInv(2) = dfdrho
      return
      end
      
      subroutine cdpm2u_computedgdInv(dgdInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the plastic potential function with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,cdpm2u_qh2fun,kappa,
     $     Bl,AGParam,BGParam,R,mQ,dgdsig,dgdrho,sig,dgdInv(2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2)        -------------  derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     dgdsig           -------------  derivative of the plastic potential with respect to volumetric stress
c     dgdrho           -------------  derivative of the plastic potential with respect to deviatoric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.      
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) /  ( log(AGParam)+
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6.) + 
     $     m0 * qh1**2. / ( sqrt(6.) * fc )
      
      dgdInv(1) = dgdsig
      dgdInv(2) = dgdrho
      return 
      end 
      
      subroutine cdpm2u_computeddgddInv(ddgddInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the derivatives of the plastic potential with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa, Bl,AGParam,
     $     BGParam,R,mQ,sig, dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho,
     $     ddgddSig,ddgddRho,ddgdSigdRho,ddgdRhodSig,ddgddInv(2,2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2,2)      -------------  derivatives of derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     ddgddsig         -------------  second derivative of the plastic potential with respect to volumetric stress
c     ddgddrho         -------------  second derivative of the plastic potential with respect to deviatoric stress
c     ddgdsigdrho      -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     ddgdrhodsig      -------------  derivative of the plastic potential with respect to deviatoric and volumetric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho, ---  variables corresponding to derivatives of components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.            
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) / ( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      dMQDSig = AGParam / ( BGParam * fc ) * exp(R)
      dAlDSig = 2. * ( 1. - qh1 ) * Bl / fc
      dBlDSig = 1. / fc
      dAlDRho = 2. * ( 1. - qh1 ) * Bl / ( fc * sqrt(6.) ) + 
     $     sqrt(3. / 2.) / fc;
      dBlDRho = 1. / ( fc * sqrt(6.) )
      
      ddgddSig = 4. * ( 1. - qh1 ) / fc * ( dAlDSig * Bl + Al * 
     $     dBlDSig ) +qh1**2. * dMQDSig / fc
      ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * Bl + 6. ) +Al * dBlDRho * 4. *
     $     ( 1. - qh1 ) / ( sqrt(6.) * fc )
      ddgdSigdRho = 4. * (1. - qh1 )/fc *( dAlDRho * Bl + Al * dBlDRho )
      ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - 
     $     qh1 ) * Bl + 6. ) + Al / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * dBlDSig )
      
      ddgddInv(1, 1) = ddgddSig
      ddgddInv(1, 2) = ddgdSigdRho
      ddgddInv(2, 1) = ddgdRhodSig
      ddgddInv(2, 2) = ddgddRho
      
      return
      end
      
      subroutine cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,
     $     sig,rho,theta,kappa)
c     Subroutine to calculate the derivative of the hardening variable with respect to plastic multiplier
      real rho,sig,theta,kappa,dkappadDeltaLambda,equivalentDGDStress,
     $     ductilityMeasure,dgdInv(2),cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     dkappadDeltaLambda -----------  derivative of the hardening variable with respect to plastic multiplier <Output>
c     dgdInv(2)        -------------  derivative of plastic potential with respect to volumeric and deviatoric stress
c     ductilityMeasure -------------  ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- norm of the derivative of the plastic potential with respect to volumetric and deviatoric stress
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      equivalentDGDStress = sqrt( 1. / 3.*dgDInv(1)** 2.+dgdInv(2)** 2.)
      ductilityMeasure = cdpm2u_computeDucMeas(sig, rho,theta)
      dkappadDeltaLambda = equivalentDGDStress / ductilityMeasure
      return
      end
      
      subroutine cdpm2u_computedfdKappa(dfdkappa,sig,rho,theta,kappa)
c     Subroutine to calculate the derivative of the yield function with respect to the cummulative plastic strain(kappa)
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      real dfdkappa,sig,rho,theta,kappa,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,dfdqh1,
     $     dfdqh2,dqh1dkappa,dqh2dkappa,cdpm2u_dqh1dkappaFun,
     $     cdpm2u_dqh2dkappaFun,Al,
     $     Bl,rFunction
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdkappa         -------------  derivative of the yield function with respect to cummulative plastic strain (kappa) <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     Al,Bl            -------------  variables corresponding to components of the yield function
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      Al = ( 1. - qh1 ) * ( ( sig / fc + rho / ( sqrt(6.) * 
     $     fc ) )) **2.  + sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      dfdqh1 = -2. *Al *(Bl** 2.) + 2. * qh1 * qh2 *   m0 * ( sig / fc +
     $     rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *qh1*(qh2**2.)

      dfdqh2 = (qh1**2.) * m0 * ( sig / fc + rho * rFunction /
     $     (sqrt(6.) * fc)) -  2. *qh2 *(qh1** 2.)
      dfdkappa =  dqh1dkappa * dfdqh1 + dqh2dkappa * dfdqh2
      
      if ( dfdkappa .gt. 0. ) then
         dfdkappa = 0.
      end if
      
      return
      end      
      
      subroutine cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $     theta,kappa)
c     Subroutine to calculate the derivative of the plastic potential function with respect to the volumetric and deviatoric stresses and the cummulative plastic strain (kappa)      
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      
      real qh1,cdpm2u_qh1fun,qh2,cdpm2u_qh2fun,dqh1dkappa,
     $     cdpm2u_dqh1dkappaFun,dqh2dkappa,
     $     cdpm2u_dqh2dkappaFun, AGParam,BGParam,R,mQ,dAGParamdKappa,
     $     BGParamTop,BGParamBottom,dBGParamTopDKappa,
     $     dBGParamBottomDKappa,dBGParamDKappa,RTop,RBottom,dRTopDKappa,
     $     dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard,
     $     dDGDSigDKappa,ddgdInvdKappa(2),kappa,sig,rho,theta,
     $     dDGDRhoDKappa
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     ddgdInvdKappa(2) -------------  derivative of the plastic potential with respect to deviatoric and volumetric strains and the cummulative plastic strain (kappa) <Output>
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     dDGDSigDKappa    -------------  derivative of the plastic potential with respect to volumetric stress and cummulative plastic strain (kappa)
c     dDGDRhoDKappa    -------------  derivative of the plastic potential with respect to deviatoric stress and cummulative plastic strain (kappa)
c     Al,Bl,AGParam,BGParam,R,mQ,dAGParamdKappa,BGParamTop,BGParamBottom,dBGParamTopDKappa,  dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard           -------------  variables corresponding to components and their derivatives of the plastic potential
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)
      
      AGParam = ft * qh2 * 3 / fc + m0 / 2
      BGParam =  qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2+ m0 / 2) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      dAGParamDKappa = dqh2dkappa * 3. * ft / fc
      BGParamTop = qh2 / 3. * ( 1. + ft / fc );
      BGParamBottom = ( log(AGParam) + log(df + 1.) - 
     $     log(2 * df - 1.) - log(3. * qh2 + m0 / 2) )
      dBGParamTopDKappa = dqh2dkappa / 3.
      dBGParamBottomDKappa = -3. * dqh2dkappa / ( 3 * qh2 + m0 / 2. )
      dBGParamDKappa =( dBGParamTopDKappa * BGParamBottom - BGParamTop * 
     $     dBGParamBottomDKappa ) / (BGParamBottom**2.)
      RTop = ( sig - ft / 3. * qh2 )
      RBottom = fc * BGParam
      dRTopDKappa = -ft / 3. * dqh2dkappa
      dRBottomDKappa = fc * dBGParamDKappa
      dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / 
     $     (RBottom** 2.)
      dMQDKappa = dAGParamDKappa * exp(R) + AGParam *dRDKappa *exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * (Bl** 2.) + sqrt(3. / 2.) * rho / fc
      dAlDYieldHard = -(Bl** 2.)
      
      dDGDSigDKappa =  ( -4. * Al * Bl / fc + 4. * ( 1 - qh1 ) / fc * 
     $     dAlDYieldHard * Bl ) * dqh1dkappa +
     $     dqh1dkappa * 2 * qh1 * mQ / fc + qh1 * dMQDKappa / fc
      dDGDRhoDKappa =
     $     ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 )* 
     $     Bl + 6. ) - 4. * Al / ( sqrt(6.) * fc ) * Bl + m0/( sqrt(6.)* 
     $     fc ) ) * 2 * qh1 * dqh1dkappa
      
      ddgdInvdKappa(1) = dDGDSigDKappa
      ddgdInvdKappa(2) = dDGDRhoDKappa
      
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadKappa(
     $     ddkappadDeltaLambdadKappa,sig,rho,kappa,theta)
c     Subroutine to compute the derivative of the cummulative plastic strain(kappa) with respect to the plastic multiplier and the hardening variable kappa
      real equivalentDGDStress,dEquivalentDGDStressDKappa,ducMeas,
     $     ddkappadDeltaLambdadKappa,dgdInv(2),ddgdInvdKappa(2),
     $     sig,rho,kappa,theta,cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddkappadDeltaLambdadKappa ----  derivative of the kappa with respect to plastic multiplier and kappa <Output> 
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     ddgdInvdKappa(2) ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and cummulative plastic strain (kappa)
c     dgdInv(2)        ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDKappa ---  scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress and kappa which is the norm of the derivative of the increment of the plastic strains with respect to kappa 
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa, sig, rho,theta,
     $     kappa)      
      equivalentDGDStress =sqrt( 1./ 3.*( dgdInv(1)**2.)+ dGDInv(2)**2.)
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivalentDGDStressDKappa = (2./3.*dgdInv(1) * ddgdInvdKappa(1) +
     $     2. * dgdInv(2) *ddgdInvdKappa(2) ) / 2./equivalentDGDStress
      ddkappadDeltaLambdadKappa= dEquivalentDGDStressDKappa/ducMeas
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadInv(
     $     ddKappaddDeltaLambdadInv, sig,rho,kappa,theta)
c     Subroutine to compute the derivative of the cummulative plastic strain (kappa) with respect to the plastic multiplier and the volumetric and deviatoric stress
      real equivDGDStress, dgdInv(2),ddgddInv(2, 2),
     $     dEquivDGDStressDInv(2),ducMeas,
     $     ddKappaddDeltaLambdadInv(2),sig,rho,kappa,theta,
     $     dDuctilityMeasureDInv(2),cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddKappaddDeltaLambdadInv -----  derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric strain <Output> 
c     ddgdd(2,2)       -------------  second derivative of the plastic potential with respect to volumetric and deviatoric stress 
c     dgdInv(2)        -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     dDuctilityMeasureDInv(2) -----  derivative of the ductility measure with respect to volumetric and deviatoric stress
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDInv(2) -- scalar function of the second derivative of the plastic potential with respect to volumetric and deviatoric stress which is the derivative of the norm of the increment of the plastic strains with respect to deviatoric and volumetric plastic strains 
      
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgddInv(ddgddInv, sig, rho,theta, kappa)
      
      equivDGDStress = sqrt( 1. / 3. * (dGDInv(1)** 2.) +
     $     dGDInv(2)**2. )
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivDGDStressDInv(1) =  ( 2. / 3.*dgdInv(1)*ddgddInv(1,1) + 
     $     2. * dgdInv(2)*ddgddInv(2, 1) ) / ( 2. * equivDGDStress)
      dEquivDGDStressDInv(2) = ( 2. / 3.*dgdInv(1)*ddgddInv(1,2) +
     $     2. * dgdInv(2)*ddgddInv(2, 2) ) /( 2. * equivDGDStress )
      call cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, rho,
     $     theta, kappa)
      
      ddKappaddDeltaLambdadInv(1) = ( dEquivDGDStressDInv(1) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(1) ) / (ducMeas** 2.)
      ddKappaddDeltaLambdadInv(2) = ( dEquivDGDStressDInv(2) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(2) ) / (ducMeas** 2.)
      
      return
      end
      
      
      subroutine cdpm2u_computeTrialCoordinates(stress, sigV,rho, theta)
c     Subroutine which returns volumetric and deviatoric stress and the Lode angle
c     based on the given stress tensor
      real stress(6),rho,sigV,theta,tempdevSig(6),j2,j3
      integer i,j
      
c     stress(6)     -------------  Stress array. stress={sigXX,sigYY,sigZZ,sigXY,sigYZ,sigXZ} <Input>
c     sigV          -------------  Volumetric stress. <Output>
c     rho           -------------  Norm of deviatoric stress. <Output>
c     theta         -------------  Lode angle.  <Output>
c     tempdevSig(6) -------------  Array containing  deviatoric stress tensor
c     j2            -------------  Second invariant of the deviatoric stress tensor
c     j3            -------------  Third invariant of the deviatoric stress tensor
c     i,j           -------------  counters used in iterations
      
      sigV= (stress(1)+stress(2)+stress(3))/3.
      do j=1,3
         tempdevSig(j)=stress(j)-sigV
         tempdevSig(j+3)=stress(j+3)
      enddo
      j2=0.5*(tempdevSig(1)*tempdevSig(1)+tempdevSig(2)*tempdevSig(2)+ 
     $     tempdevSig(3)*tempdevSig(3))+ tempdevSig(4)*tempdevSig(4)+ 
     $     tempdevSig(5)*tempdevSig(5)+tempdevSig(6)*tempdevSig(6)
      
      if(j2 .eq. 0.) then
         theta=0.
         rho=0.
      else   
         rho=sqrt(2.*j2)
         j3= (1./3.) * ( tempdevSig(1)*tempdevSig(1)*tempdevSig(1) + 
     $        3.*tempdevSig(1) *tempdevSig(4)*tempdevSig(4)+
     $        3.*tempdevSig(1) * tempdevSig(6) * tempdevSig(6) + 
     $        6. * tempdevSig(5) * tempdevSig(4)  *  tempdevSig(6) +
     $        3. * tempdevSig(2) * tempdevSig(4)**2.+
     $        3 * tempdevSig(3) *tempdevSig(6)*tempdevSig(6)+
     $        tempdevSig(2)*tempdevSig(2)*tempdevSig(2)+ 
     $        3. * tempdevSig(2) * tempdevSig(5)*tempdevSig(5)+
     $        3. * tempdevSig(3)*tempdevSig(5)*tempdevSig(5)+ 
     $        tempdevSig(3)*tempdevSig(3)*tempdevSig(3))
         theta=(3.*sqrt(3.)/2.)*j3/(j2**(3./2.))
      end if
      if (theta .gt. 1.) then
         theta=1.
      else  if (theta .lt. -1.) then
         theta=-1.
      end if
      theta=1./3.*acos(theta)
      return
      end
      
      subroutine cdpm2u_computeYieldValue(answer,sigV, rho,theta,kappa)
c     Function to evaluate the yield function f based on the given stress High -Westergaard coordinates and kappa.
c     The equation is given in eq. (18) of IJSS paper by P. Grassl et al. Returns a real number 
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      
      real rFunction,qh1,qh2, Al,sigV,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     answer

c     sigV             -------------  volumetric stress <Input>
c     rho              -------------  deviatoric stress <Input>
c     theta            -------------  Lode angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     answer           -------------  variable containing the result of the yield function <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.
c     Al               -------------  variable used to siplify and facilitate the calculation of the yield surface
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     1     ( 2. * ecc - 1. )**2.  ) /
     2     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     3     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     4     cos(theta)**2.+ 5. * ecc**2. - 4. * ecc) )                                                              

      Al = ( 1. - qh1 ) * (sigV / fc + rho / (sqrt(6.) * fc) ) ** 2. +
     &     sqrt(1.5) * rho / fc
      
      answer=Al**2. + qh1**2. * qh2* m0 * ( sigV / fc + rho *
     &     rFunction / ( sqrt(6.) * fc ) ) - qh1**2. * qh2**2.      
      return
      end
      
      
      real function cdpm2u_qh1fun(kappa,qh0,hp)
c     Function to calculate the first hardening function given in eq. (30) of the IJSS paper 
c     by Grassl et al.

      real hp,qh0, kappa,answer

c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     qh0    ---------- Initial hardening parameter <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer=qh0
      else if ((kappa .gt. 0.0) .and. (kappa .lt. 1.0)) then
         answer=
     $        (1.0 - qh0 - hp) *(kappa**3.0)-
     $        ( 3.0 * ( 1.0 - qh0 ) - 3.0 * hp ) * (kappa**2.0) +
     $        ( 3.0 * ( 1.0 - qh0 ) - 2.0 * hp ) * kappa + qh0
      else 
         answer=1.
      endif      
      cdpm2u_qh1fun=answer      
      return
      end
      
      real function cdpm2u_qh2fun(kappa,hp)
c     Function to calculate the first hardening function given in eq. (31) of the IJSS paper 
c     by Grassl et al.
      
      real kappa,answer,hp
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      
      if ( kappa .le. 0. ) then
         answer=1.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=1.
      else 
         answer= 1.+(kappa-1.)*hp
      endif      
      cdpm2u_qh2fun=answer     
      return 
      end
      
      real function cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
c     Function to calculate the derivative of the first hardening function, given in eq. (30) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real kappa,answer,hp,qh0
c     kappa  ---------- Cummulative plastic strain <Input> 
c     qh0    ---------- Initial hardening parameter <Input>
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer= 3. * ( 1 - qh0 ) - 2. * hp
      else if ( kappa .ge. 0. .and. kappa .lt. 1. ) then
         answer=   3. * ( 1. - qh0 - hp ) * (kappa**2.)
     $        - 2. * ( 3. * ( 1. - qh0 ) - 3. * hp ) * kappa
     $        + ( 3. * ( 1. - qh0 ) - 2. * hp )
      else 
         answer=  0.
      endif
      cdpm2u_dqh1dkappaFun=answer
      return
      end
      
      real function cdpm2u_dqh2dkappaFun(kappa,hp)
c     Function to calculate the derivative of the second hardening function, given in eq. (31) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real kappa,hp,answer
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      if ( kappa .le. 0. ) then
         answer=0.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=0.
      else 
         answer=hp
      endif
      cdpm2u_dqh2dkappaFun=answer
      return 
      end  


      real function cdpm2u_computeDucMeas( sigV,rho,theta)
c     Function to calculate ductility measure according to eq. (33) of IJSS paper by P. Grassl et al.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag


      real thetaConst,x,  sigV,rho,theta ,eh,fh,answer
c     sigV ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta---------------- Lode angle <Input>
c     eh,fh---------------- hardening parameters
c     thetaConst,x -------- variables
c     answer--------------- variable containing the answer of the function

      thetaConst = (2. * cos(theta))**2.
      x = -( sigV + fc / 3 ) / fc
      if ( x .lt. 0. ) then
         eh = bh - dh
         fh = ( bh - dh ) * ch / ( ah - bh )
         answer = ( eh * exp(x / fh) + dh ) / thetaConst
      else 
         answer = ( ah + ( bh - ah ) * exp( -x / ( ch ) ) ) / thetaConst
      endif
      cdpm2u_computeDucMeas=answer
      return
      end
      
      
      
      subroutine cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, 
     $     rho,theta, kappa)
c     Subroutine to compute the derivative of the ductility measure with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real dDuctilityMeasureDInv(2), sig,rho, kappa,theta,x,dXDSig,
     $     EHard,FHard,dDuctilityMeasureDX,theta1
c     sig  ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta --------------- Lode angle <Input>
c     dDuctilityMeasureDInv(2) --  derivative of the ductility measure with respect to volumetric and deviatoric stress <Output>
c     eh,fh --------------- hardening parameters
c     theta1,x ------------ variables
c     EHard,FHard --------- hardening parameters
c     dXDSig -------------- derivative of the function R, given in eq. (34) in IJSS paper by P.Grassl et al., with respect to volumetric stress
c     dDuctilityMeasureDX - derivative of the ductility measure with respect to the function R

      theta1 = (2. * cos(theta))** 2.
      x = ( -( sig + fc / 3. ) ) / fc
      
      if ( x .lt. 0. ) then
         dXDSig = -1. / fc
         EHard = bh - dh
         FHard = ( bh - dh ) * ch / ( ah - bh )
         
         dDuctilityMeasureDX = EHard  / FHard *exp(x / FHard) / theta1
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      else 
         dXDSig = -1. / fc
         dDuctilityMeasureDX = -( bh - ah ) / ( ch ) / theta1 *
     $        exp( -x / ( ch ) )
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      endif
      return
      end
c     --------------------------------------------- Functions used in the damage algorithm ---------
      subroutine cdpm2u_computeAlpha(stressT,stressC,stress,alpha)
      
c     Subroutine calculating the alpha according to eq.46 of the IJSS paper by Grassl et al. and splitting stress tensor in tensile and compressive part in the principal effective stress coordinate system. Finally principal effective tensile and compressive stress tensors are rotated back to the original system

      real stress(6),stressT(6),stressC(6)
      real princDir(3, 3),princStress(3),princStressT(3),alpha,
     $     princStressC(3), squareNormOfPrincipalStress,alphaTension
      integer i
c     stress(6)       -------------  effective stress tensor (no damage included) <Input>
c     stressT(6)      -------------  effective tensile stress tensor (no damage included) <Output>
c     alpha           ------------- alphaC used in the material model <Output>
c     stressC(6)      -------------  effective compressive stress tensor (no damage included) 
c     alphaTension    ------------- 1-alphaC <Output>
c     princDir(3,3)   -------------  matrix containing the principal directions of original stress tensor stored columnwise
c     princStress(3)  ------------- array containing principal stresses {sigma1,sigma2,sigma3} 
c     princStressT(3) ------------- array containing principal tensile stresses {sigmaT1,sigmaT2,sigmaT3} 
c     princStressC(3) ------------- array containing principal compressive stresses {sigmaC1,sigmaC2,sigmaC3} 
c     squareNormOfPrincipalStress - norm of princStress 

      call cdpm2u_computePrincValues(stress,princStress,0,princDir)
      
c     Split the principal values in a tension and a compression part
      do i = 1,3
         if ( princStress(i) .ge. 0. ) then
            princStressT(i) = princStress(i)
            princStressC(i) = 0.
         else 
            princStressC(i) = princStress(i)
            princStressT(i) = 0.
         endif
      enddo
      
c     Transform the tension and compression principal stresses back to the original coordinate system

      call cdpm2u_transformStressVectorTo(stressT, princDir,
     $     princStressT)
      call cdpm2u_transformStressVectorTo(stressC, princDir,
     $     princStressC)
      
c     Determine the two factors from the stress
      squareNormOfPrincipalStress = 0.

         squareNormOfPrincipalStress = princStress(1)**2.+
     $     princStress(2)**2.+princStress(3)**2.

      
      alphaTension = 0.
      
      if ( squareNormOfPrincipalStress .gt. 0. ) then
         do i = 1,3
            alphaTension = alphaTension+ princStressT(i) *
     $           ( princStressT(i) + princStressC(i) ) /
     $           squareNormOfPrincipalStress
         enddo
      endif

      alpha= 1. - alphaTension

      return
      end

      real function cdpm2u_computeRateFactor(alpha,strainrate,fc,
     $     tol,sratetype)
c     Function to incorporate impact effects in the constitutive law and calculate rateFactor. 
c     All functions used are based on Model Code 2010. 
      
      real strainrate(6),princDir(3,3),princStrainRate(3),max,min,tol,
     $     alphaS,gammaS, deltaS,betaS,strainRateTension0,
     $     strainRateCompression0,rate,ratioT,ratioC,rateFactorTension,
     $     rateFactorCompression,alpha,rateFactor,
     $     strainRateRatioCompression,sratetype
      integer k
c     alpha                     --------------   alpha is the variable used in CDPM2U to evaluate contribution of compressive stresses to principal stress tensor  <Input>
c     totstrain(6)             --------------   array containing strain increments {xx,yy,zz,xy,yz,xz} <Input>
c     deltaTime                 --------------    time step  length <Input>
c     oldStrainRate                 --------------    strain rate  <Input/Output>
c     rateFactor                --------------   rateFactor used to incorporate impact effects on the constitutive law <Input>
c     fc                        --------------   concrete compressive strength <Input>
c     tol                       --------------   tolerance used to identify whether it is a tensile or compressive strain state <Input>
c     princDir(3,3)             --------------   matrix containing the principal directions of the strain rate tensor. Eigenvectors stored columnwise
c     princStrainRate(6)        --------------   array containing the 3 eigenvalues of the strain rate tensor. {rate1,rate2,rate3}
c     max                       --------------   max(tensile) principal strain rate
c     min                       --------------   min (compressive) principal strain rate
c     ratioT                    --------------   max strainrate/sig0Tension (see MC90)
c     ratioC                    --------------   max strainrate/sig0Compression (see MC90)
c     rateFactorTension         --------------   contribution of tensile strain rates to the rate Factor
c     rateFactorCompression     --------------   contribution of compressive strain rates to the rate Factor
c     strainRate                 --------------    strain rate  
      rateFactorTension=1.
      rateFactorCompression=1.
      
      call cdpm2u_computePrincValues(strainRate,princStrainRate,1,
     $     princDir)
      
      max= -1.e-20
      min = 1.e20

      do k=1,3
         if (max .lt. princStrainRate(k)) then
            max=princStrainRate(k)           
         end if
         if (min .gt. princStrainRate(k)) then
            min=princStrainRate(k)
         end if
      enddo
      
      if ( 1. - alpha .gt.tol ) then 
         rate = max
      else
         rate =  min
      endif

      ratioT= rate/1.e-6

      if ( rate .lt. 1.e-6 ) then
         rateFactorTension = 1.
      else if ( 1.e-6 .lt. rate .and. rate .lt. 10.) then
         rateFactorTension = ratioT**0.018
      else 
         rateFactorTension = 0.0062 * (ratioT **(1./3.))
      endif

      
      ratioC= rate/(-30.e-6)

      if ( rate .gt. -30.e-6 ) then
         rateFactorCompression = 1.
      else if (-30.e-6 .gt. rate .and. rate .gt. -30) then
         rateFactorCompression = ratioC**(0.014)        
      else 
         rateFactorCompression =  0.012*(ratioC**(1./3.))
      endif

      
      rateFactor = ( 1. - alpha ) * rateFactorTension + 
     $     alpha * rateFactorCompression
      cdpm2u_computeRateFactor=rateFactor

      return
      end

      subroutine cdpm2u_computeDamage(omegaC,omegaT,strainRate,
     $     rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,kappaDT2,
     $     kappaDC,kappaDC1,kappaDC2,stress,deltaPlasticStrainNormT,
     $     tempKappaP,len,stressOld,oldAlpha,epsilon)
c     Subroutine to perform the damage return. Both compressive and tensile damage variables are calculated.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag


      real omegaC,omegaT,tempRateFactor,rateFactor,epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     stress(6),deltaPlasticStrainNormT,deltaPlasticStrainNormC,
     $     tempKappaP,len,omegaOldC,omegaOldT,alpha,sigElastic,
     $     rhoElastic, thetaElastic,pHelp,help,tempEquivStrain,qhelp,
     $     tempEquivStrainT,tempEquivStrainC,fsT,fsC,ducMeas,Rs,
     $     strainRate(6),yieldTolDamage,cdpm2u_computeRateFactor,
     $     deltaPlasticStrainNormTNew,
     $     cdpm2u_computeDeltaPlasticStrainNormT,
     $     deltaPlasticStrainNormCNew,
     $     cdpm2u_computeDeltaPlasticStrainNormC,
     $     alphaZero,e0,rFunction,
     $     cdpm2u_computeDamageT,cdpm2u_computeDamageC,
     $     cdpm2u_computeEquivalentStrain,stressOld(6),
     $     minEquivStrain,oldAlpha,epsilon
      integer unloadingFlag,step1flag
c     omegaC          --------------- compressive damage variable                        <Input/Output>
c     omegaT          --------------- tensile damage variable                            <Input/Output>
c     totstrain(6)    --------------- array containing the rate of the strain tensor     <Input>
c     strainRate      ---------------  strain rate     <Input/Output>
c     deltaTime       ---------------  time step length     <Input>
c     rateFactor      --------------- variable to incorporate impact effects on the constitutive law  <Input/Output>
c     alpha           --------------- variable showing the contribution of damage on the strains      <Input/Output>
c     epsilonT        --------------- tensile equivalent strain                <Input/Output>
c     epsilonC        --------------- compressive equivalent strain            <Input/Output>
c     epsilon         ---------------  old equivalent strain (no rate influnce)            <Input/Output>
c     kappaDT         --------------- history parameter kappaDT                <Input/Output>
c     kappaDT1        --------------- history parameter kappaDT1               <Input/Output>
c     kappaDT2        --------------- history parameter kappaDT2               <Input/Output>
c     kappaDC         --------------- history parameter kappaDC                <Input/Output>
c     kappaDC1        --------------- history parameter kappaDC1               <Input/Output>
c     kappaDC2        --------------- history parameter kappaDC2               <Input/Output>
c     stress(6)       --------------- stress tensor                            <Input>
c     deltaPlasticStrainNormT ------- norm of the increment of plastic strains <Input>
c     tempKappaP      --------------- cummulative plastic strain (kappaP)      <Input>
c     deltaElStrain(6)  ---------------  increment of the elastic strain vector  <Input>
c     oldAlpha        ---------------  old alpha  <Input>
c     len             --------------- characteristic length used to combine damage law with the crack-band approach <Input>
c     deltaPlasticStrainNormTNew ---- norm of the increment of the cummulative plastic strain used in the damage algorithm for tension
c     cdpm2u_computeDeltaPlasticStrainNormT- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the tensile damage variable
c     deltaPlasticStrainNormCNew ---- norm of the increment of the cummulative plastic strains used in the damage algorithm for compression
c     cdpm2u_computeDeltaPlasticStrainNormC- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the compressive damage variable
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     alphaZero       --------------- parameter used for the calculation of the damage ductility measure based on eqs.(56-57) of the IJSS paper by Grassl et al.
c     e0              --------------- parameter equal to ft/E
c     cdpm2u_computeDamageT,cdpm2u_computeDamageC - functions used to calculate tensile and compressive damage variables respectively
c     yieldTolDamage  --------------- tolerance used in the damage algorithm for cased where Hp=0
c     Rs              --------------- variable used in the calculation of the ductility measure for damage given in eq.(57) of the IJSS paper by Grassl et al.
c     xs              --------------- ductility measure for damage given in eq.(56) of the IJSS paper by Grassl et al.
c     fsT,fsC         --------------- loading function for the tensile and compressive damage variables
c     cdpm2u_computeRateFactor ------------- function used for the calculation of the rate factor to include impact effects on the constitutive law
c     omegaOldT,omegaOldC ----------- old (previous step) tensile and compressive damage variables
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression)
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     minEquivStrain    -------------- when unloading is occuring corresponds to the minEquivStrain before reloading occurs
c     step1flag    -------------- flag indicating whether we are at the first step
c                                        =1 it is analysis 1st step
c                                        =0  it not analysis 1st step

      e0=ft/ym
      yieldTolDamage=gTol*10.0
      omegaOldC=omegaC
      omegaOldT=omegaT
      deltaPlasticStrainNormC=deltaPlasticStrainNormT
      if (rateFactor .eq. 0.) then
         step1flag=1
         rateFactor=1.
      else
         step1flag=0
      end if
      call cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,tempEquivStrain,epsilon,ym,pr,
     $     gtol)
      
c-------------------Compute tensile and compressive equivalent strains-------------------------------------
      if (strrateflg .gt. 0.0 .and. omegaC .eq. 0. .and. 
     $     omegaT.eq. 0. ) then
         tempRateFactor=cdpm2u_computeRateFactor(alpha,strainrate,
     $        fc,gTol,strrateflg)
      else 
         tempRateFactor=rateFactor
      end if
      tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/
     $     tempRateFactor
      
      if (unloadingFlag .eq. 0) then
         tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $        tempRateFactor
      else
          tempEquivStrainC=epsilonC+ oldAlpha*(minEquivStrain-epsilon)/
     $        tempRateFactor + alpha*(tempEquivStrain-minEquivStrain)/
     $         tempRateFactor
      end if
   
c     Note rate factor is calculated only once at the onset of damage
      if ( ( tempEquivStrainT .gt. e0 .or. tempEquivStrainC .gt. e0
     $     ) .and. ( ( omegaT .eq. 0. ) .and. (omegaC .eq. 0. ) ).and.
     $     strrateflg .gt. 0.0 .and. step1flag .ne. 1) then
         tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/rateFactor
         if (unloadingFlag .eq. 0) then
            tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $           rateFactor
         else
            tempEquivStrainC=epsilonC+oldAlpha*(minEquivStrain-
     $           epsilon)/rateFactor+(tempEquivStrain-minEquivStrain)*
     $           alpha/rateFactor 
         end if
      else
         rateFactor = tempRateFactor 
      endif
      
      fsT = (tempEquivStrainT - kappaDT)/e0
      fsC = (tempEquivStrainC - kappaDC)/e0

      epsilon=tempEquivStrain
      epsilonT=tempEquivStrainT
      epsilonC=tempEquivStrainC
c -------------------- Compute Ductility Measure Damage --------------------------------------------------
      call cdpm2u_computeTrialCoordinates(stress, sigElastic,
     $rhoElastic,thetaElastic)
      Rs = 0.
      alphaZero= 1./sqrt(6.) 
      if ( sigElastic .lt. 0. ) then
         if ( rhoElastic .gt. 1.e-16 ) then
            Rs = -sigElastic /(alphaZero*rhoElastic)
         else 
            Rs = -sigElastic * 1.e16 / alphaZero
         endif
      else 
         Rs = 0.
      endif
      ducMeas = 1. + ( as - 1. ) * Rs ** bs
      
c----- Check which damage surfaces (tensile/compressive) are active --------------------------------------
      if (fsT .lt. -yieldTolDamage .and. 
     $     fsC .lt. -yieldTolDamage) then
c     no increase of the damage variables required
      else if (  fsT .ge. -yieldTolDamage .and. 
     $        fsC .lt. -yieldTolDamage) then
c     only tensile damage surface active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT, deltaPlasticStrainNormT,kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / ducMeas / rateFactor
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / 
     $        ducMeas
         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT,rateFactor)
      else if (  fsT .lt.  -yieldTolDamage .and. 
     $        fsC .ge.  -yieldTolDamage) then
c     only compressive damage surface active
         deltaPlasticStrainNormCNew =
     $        cdpm2u_computeDeltaPlasticStrainNormC(alpha, 
     $        tempEquivStrainC,
     $        deltaPlasticStrainNormC, kappaDC,rhoElastic,
     $        tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas*rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC

         omegaC = 
     $        cdpm2u_computeDamageC(kappaDC, 
     $        kappaDC1, kappaDC2, omegaOldC,rateFactor)
      else if (  fsT .ge.-yieldTolDamage  .and. 
     $        fsC .ge.-yieldTolDamage ) then
c Both compressive and tensile damage surfaces are active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT,deltaPlasticStrainNormT, kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / (ducMeas * rateFactor)
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / ducMeas

         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT,rateFactor)
c     only compressive damage surface active
         deltaPlasticStrainNormCNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormC(
     $        alpha,tempEquivStrainC,deltaPlasticStrainNormC, 
     $        kappaDC,rhoElastic,tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas * rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC
         omegaC = cdpm2u_computeDamageC(kappaDC, kappaDC1, 
     $        kappaDC2, omegaOldC,rateFactor)
      endif
    
      return
      end

      subroutine cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,equivStrainNew,equivStrainOld,
     $     ym,pr,gtol)
c     Function to check if there is unloading and reloading occuring within one step (usually happens during cyclic loading). If such a process is happening the algorithm returns a flag and an approximate value of the minimum equivalent strain. Moreover it returns the current equivalent strain.

      real stress(6),minEquivStrain,stress1(6),ym,pr,
     $     deltaStress(6),stressPlus(6),stressMinus(6),
     $     equivStrainOld,equivStrainNew,equivStrain1,stressOld(6),
     $     equivStrainPlus,equivStrainMinus,sigV,rho,theta,
     $     cdpm2u_computeEquivalentStrain,gtol
      integer i,j,unloadingFlag
c     stress(6)   ------------- effective stress tensor(no damage included) <Input>
c     stressOld(6) ---------- effective stress vector in  previous step <Input>
c     equivStrainOld ---------- equivalent strain of the previous loading step <input>
c     ym             ---------- Young's modulus <Input>
c     pr             ---------- Poisson's ratio <Input>
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression) <Output>
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     equivStrainNew ---------- equivalent strain of the current loading step  <Output>
c     minEquivStrain   -------- minimum equivalent strain <Output>
c     stressPlus    ----------   stress vector representing previous step's elastic strain vector plus 0.99*deltaStrain
c     stressMinus   ----------   stress vector representing current step's elastic strain vector minus 0.01*deltaStrain
c     sigV          ---------- volumentric stress
c     rho           ---------- deviatoric stress
c     theta         ---------- Lode angle
c     equivStrain1  ---------- equivalent strain of various strain vectors
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al

      call cdpm2u_computeTrialCoordinates(stress, sigV, rho,theta)
      equivStrainNew= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
    
      do i=1,6
        deltaStress(i)=stress(i)-stressOld(i)
      enddo
      do i=1,6
         stressPlus(i)=stressOld(i)+0.01*deltaStress(i)
         stressMinus(i)=stressOld(i)+0.99*deltaStress(i)
      enddo

      call cdpm2u_computeTrialCoordinates(stressPlus, sigV, rho,theta)
      equivStrainPlus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      call cdpm2u_computeTrialCoordinates(stressMinus, sigV, rho,theta)
      equivStrainMinus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      unloadingFlag=0
      minEquivStrain=equivStrainOld
      if ( (equivStrainPlus .lt. equivStrainOld .and.
     $     equivStrainMinus .lt. equivStrainNew) .and. 
     $     (abs(equivStrainPlus - equivStrainOld) .gt. gtol/10. .and. 
     $     abs(equivStrainMinus - equivStrainNew) .gt. gtol/10.)) then
         unloadingFlag=1
         write(*,*) '*** Unloading and reloading occurs.'
         write(*,*) 'within a single step. Subincrementation performed' 
         do i=1,100
            do j=1,6
               stress1(j)=stressOld(j)+deltaStress(j)*FLOAT(i)/100.0 
            enddo
            call cdpm2u_computeTrialCoordinates(stress1,sigV,rho,theta)
            equivStrain1= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
            if( equivStrain1 .le. minEquivStrain) then
               minEquivStrain=equivStrain1
            else
               goto 625
            end if
         enddo
      end if
      
 625  continue
      return
      end
      
      real function  cdpm2u_computeEquivalentStrain(sigElastic,
     $     rhoElastic, thetaElastic)
c     Function to calculate equivalent strain used in the damage part according to eq. 37 of the IJSS paper by Grassl et al.

      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real  answer, rFunction,thetaElastic,rhoElastic,
     $     sigElastic, pHelp,qHelp,help,e0
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     pHelp,help,qhelp -------------- variables used in the calculation of the equivalent strain
c     rFunction       --------------- parameter given by eq.(19) of the IJSS paper by Grassl et al.      
c     e0              --------------- parameter equal to ft/E
      e0=ft/ym
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(thetaElastic)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(thetaElastic) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     $     cos(thetaElastic)**2.+ 5. * ecc**2. - 4. * ecc) )
      pHelp = -m0 * ( rhoElastic * rFunction / ( sqrt(6.) * fc ) + 
     $     sigElastic / fc )
      qHelp = -3. / 2. * (rhoElastic** 2.) / (fc**2.)
      help = -0.5 * pHelp + sqrt((pHelp** 2.) / 4. - qHelp)
c     negative help values are not of interest and create problems since we compute the square root 
      answer = 0.
      if ( help .gt. 0. ) then
         answer= help * e0
      endif
      cdpm2u_computeEquivalentStrain=answer
      return 
      end
      
      real function  cdpm2u_computeDamageT(kappa,kappaOne,kappaTwo,le, 
     $      omegaOld,rateFactor)
c     Function to calculate damage due to tension.

      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real kappa, kappaOne, kappaTwo, omegaOld,residual,le,
     $     residualDerivative,omega,tol,help,e0,yieldTolDamage,
     $     rateFactor,wf1Mod,wfMod,residualStrength
      integer iter,newtonIter
c     kappa           ----------------------- damage history parameter kappaDT <input>
c     kappaOne        ----------------------- damage history parameter kappaDT1 <input>
c     kappaTwo        ----------------------- damage history parameter kappaDT2 <input
c     omegaOld        ----------------------- old damage variable (previous step) omegaT <input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     residualDerivative--------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaT
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaT for exponential softening law
c     help            -----------------------  help variable used in the billinear softening law
c     iter            ----------------------- counter of performed N-R iterations in the exponential softening law
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaT in the exponential softening law
c     yieldTolDamage  ----------------------- tolerance used in the damage algorithm if Hp=0

      wfMod = wf
      wf1Mod = wf1
      ftTemp = ft*(1-yieldTolDamage)
      if(strrateflg .eq. 2) then
         wfMod = wf/rateFactor
         wf1Mod = wf1/rateFactor
      else if(strrateflg .eq. 3) then
         wfMod = wf/(rateFactor*rateFactor)
         wf1Mod = wf1/(rateFactor*rateFactor)
      end if

      newtonIter=100
      tol=gTol/100.0
      e0=ft/ym
      yieldTolDamage=gTol*10.0
      residualStrength=1.e-3*ft
      if ( kappa  .gt. e0*(1-yieldTolDamage) ) then
        if ( type .eq. 0. ) then
c     Linear damage law
           omega = ( ym * kappa * wfMod - ftTemp * wfMod +
     $      ftTemp * kappaOne * le ) /
     $          ( ym * kappa * wfMod - ftTemp * le * kappaTwo )
           help = le * kappaOne + le * omega * kappaTwo
           if ( help .ge. 0. .and. help .lt. wfMod .and. 
     $          (1-omega)*ym*kappa .gt. residualStrength) then
              goto 185
           endif
c           write(*,*) 'Residual stress computed\n'
           omega = 1.-1.e-3*ft/(ym*kappa)            
        else if ( type .eq. 1. ) then
c     Bilinear damage law
           omega = ( ym * kappa * wf1Mod - ftTemp * wf1Mod - ( ft1 -
     $      ftTemp ) * 
     $          kappaOne * le ) /( ym * kappa * wf1Mod + 
     $          ( ft1 - ftTemp ) * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo
            if ( help .ge. 0. .and. help .lt. wf1Mod .and. 
     $           (1-omega)*ym*kappa .gt. residualStrength) then
               goto 185
            endif
            
            omega = ( ym * kappa * ( wfMod - wf1Mod ) -
     $       ft1 * ( wfMod - wf1Mod ) +
     $           ft1 * kappaOne * le  - ft1 * wf1Mod ) / ( ym * kappa * 
     $           ( wfMod - wf1Mod )  - ft1 * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo

            if ( help .gt. wf1Mod .and. help .lt. wfMod  .and. 
     $           (1-omega)*ym*kappa .gt. residualStrength) then
               goto 185
            endif
c            write(*,*) 'Residual stress computed\n'
            omega = 1.-1.e-3*ftTemp/(ym*kappa)
            
         else if ( type .eq. 2. ) then
c     Exponential: Iterative solution with N-R procedure
            omega = 1.
            residual=0.
            residualDerivative = 0.
            iter = 0

 135        continue
            iter=iter+1
            residual = ( 1 - omega ) * ym * kappa - ftTemp *
     $           exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod)
            residualDerivative = -ym * kappa + ftTemp * le * 
     $           kappaTwo / wfMod * exp(-le * ( omega * kappaTwo + 
     $           kappaOne ) / wfMod)
            omega =omega- residual / residualDerivative;
            if ( iter .gt. newtonIter ) then
               write(*,*) '*** Algorithm for tensile damage-
     $No convergence reached after 100 iterations ***'
               stop
            end if
            if (abs(residual/ftTemp) .ge. 1.e-8) then
               goto 135
            end if
         end if
      else 
         omega = 0.;
      endif
      
       
       if ( omega .gt. 1. ) then
          omega = 1.
       endif
       
       if ( omega .lt. 0. .or. omega .lt. omegaOld) then
          omega=omegaOld
       endif
 185   cdpm2u_computeDamageT= omega
       return 
       end  

      real function cdpm2u_computeDamageC(kappa, kappaOne, 
     $     kappaTwo, omegaOld, rateFactor)
c     Function to calculate damage due to compression
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      
      real kappa, kappaOne, kappaTwo, omegaOld,residual,dResidualDOmega,
     $     exponent,omega,tol,kappaDC,e0,yieldTolDamage,rateFactor,
     $     efcMod
      integer nite,newtonIter
c     kappa           ----------------------- damage history parameter kappaDC <Input>
c     kappaOne        ----------------------- damage history parameter kappaDC1 <Input>
c     kappaTwo        ----------------------- damage history parameter kappaDC2 <Input
c     omegaOld        ----------------------- old damage variable (previous step) omegaC <Input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     dResidualDOmega ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     exponent        ----------------------- exponent used in the formulation of the damage law (in the current version assumed =1.)
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaC
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaC
c     nite            ----------------------- counter of performed N-R iterations
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaC

      ftTemp = ft*(1-yieldTolDamage)
      efcMod = efc
      if(strrateflg .eq. 2) then
         efcMod = efc/rateFactor
      else if(strrateflg .eq. 3) then
         efcMod = efc/(rateFactor*rateFactor)
      end if

      
      if (damageflag .eq. 1.0) then
         omega=0.0
         kappaOne=0.0
         kappaTwo=0.0
         kappa=0.0
         goto 188
      end if
      
      newtonIter=100
      tol=gTol/100.0
      yieldTolDamage=gTol*10.

      omega=1.
      nite = 0
      residual = 0.
      dResidualDOmega = 0.
      exponent = 1.
      e0=ft/ym
      if ( kappa .gt. e0*(1-yieldTolDamage)) then
 187     continue 
         nite=nite+1;
         residual =  ( 1. - omega ) * ym * kappa - ftTemp * exp( - ( 
     $        kappaOne + omega * kappaTwo ) / efcMod )
         dResidualDOmega =-ym * kappa + ftTemp *
     $   kappaTwo / efcMod * exp( -( 
     $        kappaOne + omega * kappaTwo ) / efcMod )
         omega = omega- residual / dResidualDOmega
         if ( nite .gt. newtonIter ) then
            write(*,*) '*** Algorithm for compressive damage-
     $No convergence reached.'
c Set omega=omegaOld =', omegaOld,' ***'
c            omega = omegaOld
            stop
         endif
         if( abs(residual/ftTemp) .ge. tol )   goto 187
         
       else 
          omega = 0.
       endif

       if ( omega .gt. 1 ) then
          omega = 1.
       endif
       
       if ( omega .lt. 0. .or. omega .lt. omegaOld ) then
          omega = omegaOld
       endif
       
 188   cdpm2u_computeDamageC=omega
       return
       end

      
      real function  cdpm2u_computeDeltaPlasticStrainNormT(tempKappaD,
     $     plastStrNorm, kappaD,e0,yieldTolDamage)      
c     Function returning the norm of the increment of the plastic strain tensor.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch

      real e0,yieldTolDamage,answer,tempKappaD,plastStrNorm, kappaD,
     $     factor
c     tempKappaD      -------------------    temporary (current) KappaDt                   <Input>
c     kappaD          -------------------    old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------    plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------    variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------    tolerance used when Hp=0                      <Input>
c     factor          -------------------    factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     answer          -------------------    calculated norm             
      factor = 0.
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0  * ( 1. - yieldTolDamage )) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer = plastStrNorm*factor
      else 
         answer=plastStrNorm 
      end if
      cdpm2u_computeDeltaPlasticStrainNormT=answer
      return 
      end

      real function  cdpm2u_computeDeltaPlasticStrainNormC( 
     $     alpha,tempKappaD,
     $     plastStrNorm, kappaD,rho,tempKappa,yieldTolDamage)
c     Function returning the norm of the increment of the plastic strain tensor multiplied by alphaC and betaC according to eq. (48) of the IJSS paper by Grassl et al.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch      
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag

      real   tempKappaD,plastStrNorm, kappaD,factor,qh2,alpha,rho,e0,
     $     cdpm2u_qh2fun,extraFactor,tempKappa,kappa,
     $     yieldTolDamage,answer

c     tempKappaD      -------------------  temporary (current) KappaDt                   <Input>
c     kappaD          -------------------  old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------  plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------  variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------  tolerance used when Hp=0                      <Input>
c     rho             -------------------  deviatoric stress                             <Input>
c     alpha           -------------------  variable calculated based on eq. 46 of the IJSS paper by P. Grassl et al.    <Input>
c     tempKappa       -------------------  temporary (current step) cummulative plastic strain (kappaP)   <Input>
c     factor          -------------------  factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     qh2             -------------------  variables containing the results of the hardening function 
c     cdpm2u_qh2fun          -------------------  functions to calculate the hardening functions  given in eq. (31) in IJSS paper by Grassl et al.
c     extraFactor     -------------------  variable betaC given in eq. (50) of the IJSS paper by Grassl et al.
c     answer          -------------------  calculated answer         
      e0=ft/ym
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer=plastStrNorm*factor
      else 
         answer=plastStrNorm
      end if
      
      qh2=cdpm2u_qh2fun(tempKappa,hp)
      
      if (rho<1.e-16) then
         extraFactor =ft * qh2 * sqrt(2. / 3.) / 1.e-16 / sqrt( 1. + 
     $        2.*  (df** 2.) )
      else 
         extraFactor =ft * qh2 * sqrt(2. / 3.) / rho / sqrt( 1. + 
     $        2.* (df** 2.) )
      endif
      answer=answer*extraFactor*alpha
      cdpm2u_computeDeltaPlasticStrainNormC=answer
      return
      end
      

c ---------------------------------------------General functions--------------------------------

      subroutine cdpm2u_computePrincValues(ax,w,flag,z)
c     Subroutine to solve the eigenvalues and eigenvectors of real symmetric matrix by jacobi method.
c     Eigenvalues are stored in descending order (from max to min) in array w and eigenvectors
c     are stored columnwise in matrix z.

      real ax,w,z,matrix
      integer flag
      dimension ax(6),w(3), z(3,3),matrix(3,3)
      real ssum, aa, co, si, tt, tol, sum, aij, aji;
      integer ite, i, j, k, ih
      
      real swap
      integer ii,jj,kk
c     ax(6)    ---------------------            Input tensor. It is a 6x1 array and components are given 
c                                               in order xx,yy,zz,xy,yz,xz  <Input>
c     flag     ---------------------            Flag denoting whether it is a stress or a strain tensor in 
c                                                order to convert voigt shear strains(gammas) to tensorial
c                                                shear strains(epsilons).
c                                                = 0: Stress tensor
c                                                = 1: Strain tensor
c     w(3)     --------------------             3x1 Array containing eigenvectors of ax stored in 
c                                               descending order(max to min)  <Output>
c     z(3,3)   --------------------             Matrix containing all eigenvectors stored columnwise <Output>
c     matrix(3,3) -----------------             Stress/Strain tensor that is being processed
c     tol      --------------------             Tolerance of the algorithm (default 10 figures,1e-10)
c     ite,i,j,k,ih ----------------             Counters used in various loops
c     sum,ssum --------------------             Summation variables used in various addition procedures
c     aa,co,si,tt,aij,aji ---------             Variables used in the algorithm 

      tol=1.e-10
c     Reconstruct tensor based on given array
      matrix(1,1)=ax(1)
      matrix(2,2)=ax(2)
      matrix(3,3)=ax(3)
      if (flag .eq. 0) then
         matrix(2,1)=ax(4)
         matrix(1,2)=ax(4)
         matrix(3,1)=ax(6)
         matrix(1,3)=ax(6)
         matrix(3,2)=ax(5)
         matrix(2,3)=ax(5)
      else
         matrix(2,1)=ax(4)/2.
         matrix(1,2)=ax(4)/2.
         matrix(3,1)=ax(6)/2.
         matrix(1,3)=ax(6)/2.
         matrix(3,2)=ax(5)/2.
         matrix(2,3)=ax(5)/2.
      endif
      
c     Initialise w,z and check if zero stress state
      do i=1,3
         w(i) = matrix(i, i)
      enddo
      sum=0.
      do i=1,3
         do j=1,3
            sum =sum+ abs( matrix(i, j) )
            z(i, j) = 0.0
         enddo
         z(i, i) = 1.0
      enddo
      if ( sum .le. 0.0 ) then
         goto 900
      endif
            
c     Reduce to matrix diagonal
      ite=0

 272  continue
      ssum = 0.0
      do j=2,3
         ih = j - 1
         do i=1,ih
            if ( abs( matrix(i, j) ) / sum  .gt. tol ) then
               ssum =ssum+ abs( matrix(i, j) )
c     CALCULATE ROTATION ANGLE
               aa = atan2( matrix(i, j) * 2.0, w(i) - w(j) ) /  2.0
               si = sin(aa)
               co = cos(aa)
               
c     MODIFY "I" AND "J" COLUMNS OF "matrix" and "z"
               do k=1,i-1
                  tt = matrix(k, i)
                  matrix(k, i) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (i,i)
               tt = w(i)
               w(i) = co * tt + si *matrix(i, j)
               aij = -si * tt + co *matrix(i, j)
               tt = z(i, i)
               z(i, i) = co * tt + si *z(i, j)
               z(i, j) = -si * tt + co *z(i, j)
               
               do k=i+1,j-1
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (j,j)
               tt = matrix(i, j)
               aji = co * tt + si *w(j)
               w(j) = -si * tt + co *w(j)
               
               tt = z(j, i)
               z(j, i) = co * tt + si *z(j, j)
               z(j, j) = -si * tt + co *z(j, j)
               
               do k=j+1,3
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(j, k)
                  matrix(j, k) = -si * tt + co *matrix(j, k)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     MODIFY DIAGONAL TERMS
               w(i) = co * w(i) + si * aji
               w(j) = -si * aij + co *w(j)
               matrix(i, j) = 0.0
            else 
c     matrix(I,J) MADE ZERO BY ROTATION
            endif
         enddo
      enddo
      
      ite=ite+1
c     CHECK FOR CONVERGENCE
      if ( ite .gt. 50 ) then
         write(*,*)  '*** Compute principal values.
     $Too many iterations! ***'
         stop
      endif
      
      if ( abs(ssum) / sum .gt. tol  ) goto 272
      
      do ii=1,2
         do jj=1,2
            if ( w(jj + 1) > w(jj) ) then
c     swap eigenvalues and eigenvectors
               swap = w(jj + 1);
               w(jj + 1) = w(jj);
               w(jj) = swap;
               do kk=1,3
                  swap = z(kk, jj + 1);
                  z(kk, jj + 1) = z(kk, jj);
                  z(kk, jj) = swap;
               enddo
            endif
         enddo
      enddo
 900  continue
      return
      end
      
      subroutine cdpm2u_computeInverseJac(jacinv,jac,error)
c     Subroutine to calculate the inverse of a 4x4 matrix
      real jacinv(4,4),jac(4,4),tmp(4,4)
      real  piv, linkomb,dtol
      integer i,j,k,error
c     jac(4,4)    ---------  matrix whose inverse will be calculated <Input>
c     jacinv(4,4) ---------  inverse matrix  <Output>
c     error       --------- integer showing whether solution has converged <Output>
c                            = 0 solution has converged
c                            =-1 solution has not converged
c     tmp(4,4)    --------- temporary matrix
c     piv,lincomb --------- variables used in gaussian elimination
c     i,j,k       --------- counters
c     dtol        --------- tolerance of the algorithm

      dtol=1.e-20           
      do i=1,4
         do j=1,4
            tmp(i,j)=jac(i,j)
            jacinv(i,j)=0.
         enddo
         jacinv(i,i)=1.
      enddo
      
      do i=1,3
         piv = tmp(i, i)
         if (abs(piv) .lt. dtol) then
            error=-1
            goto 452
         endif

         do j=i+1,4
            linkomb = tmp(j, i) / tmp(i, i)
            do k=i,4
               tmp(j, k) = tmp(j, k) - tmp(i, k) * linkomb
            enddo
            do k=1,4
               jacinv(j, k) =jacinv(j, k)-jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=4,2,-1
         piv = tmp(i, i)
         do j=i-1,1,-1
            linkomb = tmp(j, i) / piv
            do k=i,1,-1
               tmp(j, k) =  tmp(j, k) -tmp(i, k) * linkomb
            enddo
            do k=4,1,-1
               jacinv(j, k) =jacinv(j, k) - jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=1,4
         do j=1,4
            jacinv(i, j) = jacinv(i, j) / tmp(i, i)
         enddo
      enddo 
      error=0
 452  continue
      return
      end
      
      
      
      subroutine cdpm2u_transformStressVectorTo(stress,princDir,
     $     princStress)
c     Rotates principal stress tensor back to the original Coordinate system. It calculates transformation 
c     matrix  and then multiplies with principal stress tensor.

      real stress(6),princDir(3,3),princStress(3),transpose(3,3),
     $     transformationMtrx(6,6),sum,princStressTensor(6)
      integer i,j
c     stress (6)             ----------------- stress vector <Output>
c     princDir(3,3)          ----------------- matrix containing eigenvectors stored columnwise <Input>
c     princStress(3)         ----------------- stress vector containing principal stresses =
c                                              {sigma1,sigma2,sigma3}  <Input>
c     princStressTensor(6)   ----------------- stress vector at principal axis coordinate system 
c                                              ={sigma1,sigma2,sigma3,0,0,0}
c     transpose(3,3)         ----------------- matrix containing eigenvectors stored rowise
c     sum                    ----------------- variable used in summation during matrix multiplication
c     transformationMtrx(6,6)----------------- transformation matrix used to transform principal stress
c                                              vector to stress vector in original CS
c     i and j                -----------------  integers used as counters


      do i=1,3
         princStressTensor(i)= princStress(i)
         princStressTensor(i+3)= 0.
         do j=1,3
            transpose(i,j)=princDir(j,i)
         enddo
      enddo
      
      
      transformationMtrx(1,1)=transpose(1,1)*transpose(1,1)
      transformationMtrx(1,2)=transpose(2,1)*transpose(2,1)
      transformationMtrx(1,3)=transpose(3,1)*transpose(3,1)
      transformationMtrx(1,4)=2.*transpose(1,1)*transpose(2,1)
      transformationMtrx(1,5)=2.*transpose(2,1)*transpose(3,1)
      transformationMtrx(1,6)=2.*transpose(1,1)*transpose(3,1)

      transformationMtrx(2,1)=transpose(1,2)*transpose(1,2)
      transformationMtrx(2,2)=transpose(2,2)*transpose(2,2)
      transformationMtrx(2,3)=transpose(3,2)*transpose(3,2)
      transformationMtrx(2,4)=2.*transpose(1,2)*transpose(2,2)
      transformationMtrx(2,5)=2.*transpose(2,2)*transpose(3,2)
      transformationMtrx(2,6)=2.*transpose(1,2)*transpose(3,2)

      transformationMtrx(3,1)=transpose(1,3)*transpose(1,3)
      transformationMtrx(3,2)=transpose(2,3)*transpose(2,3)
      transformationMtrx(3,3)=transpose(3,3)*transpose(3,3)
      transformationMtrx(3,4)=2.*transpose(1,3)*transpose(2,3)
      transformationMtrx(3,5)=2.*transpose(2,3)*transpose(3,3)
      transformationMtrx(3,6)=2.*transpose(1,3)*transpose(3,3)

      transformationMtrx(4,1)=transpose(1,1)*transpose(1,2)
      transformationMtrx(4,2)=transpose(2,1)*transpose(2,2)
      transformationMtrx(4,3)=transpose(3,1)*transpose(3,2)
      transformationMtrx(4,4)=transpose(1,1)*transpose(2,2)+
     $     transpose(2,1)*transpose(1,2)
      transformationMtrx(4,5)=transpose(2,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(2,2)
      transformationMtrx(4,6)=transpose(1,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(1,2)

      transformationMtrx(5,1)=transpose(1,2)*transpose(1,3)
      transformationMtrx(5,2)=transpose(2,2)*transpose(2,3)
      transformationMtrx(5,3)=transpose(3,2)*transpose(3,3)
      transformationMtrx(5,4)=transpose(1,2)*transpose(2,3)+
     $     transpose(2,2)*transpose(1,3)
      transformationMtrx(5,5)=transpose(2,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(2,3)
      transformationMtrx(5,6)=transpose(1,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(1,3)
    

      transformationMtrx(6,1)=transpose(1,1)*transpose(1,3)
      transformationMtrx(6,2)=transpose(2,1)*transpose(2,3)
      transformationMtrx(6,3)=transpose(3,1)*transpose(3,3)
      transformationMtrx(6,4)=transpose(1,1)*transpose(2,3)+
     $     transpose(2,1)*transpose(1,3)
      transformationMtrx(6,5)=transpose(2,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(2,3)
      transformationMtrx(6,6)=transpose(1,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(1,3)

 
      do i=1,6
         sum=0.
         do j=1,6
            sum=sum+  transformationMtrx(i,j)*princStressTensor(j)
         enddo
         stress(i)=sum
      enddo
c
      return
      end


      subroutine utan50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat)
c
c*****************************************************************
c Implementation of the stiffness for CDPM2 in UMAT50V in LS-DYNA
c Support page: http://petergrassl.com/Research/DamagePlasticity/CDPMLSDYNA/index.html
c
C Only elastic stiffness so far.
c
c      Authors of this subroutine: Peter Grassl
c
c Last updated: 20 April 2016
c
c*****************************************************************

      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
      
      real es(6,6)
      real ym,pr
c
c Define elastic stiffness to start with. Needs to be improved later.
c
c
      do k=1,6
         do j=1,6
            es(j,k) =0.
         enddo
      enddo

      ym=cm(1)
      pr=cm(2)      

      const = ym/((1.+pr)*(1.-2.*pr))
      es(1,1)= const*(1.-pr)
      es(2,2)= es(1,1)
      es(3,3) = es(1,1)
      es(4,4) = const*(1.-2.*pr)/2.
      es(5,5) = es(4,4)
      es(6,6) = es(4,4)
      es(1,2) = const*pr
      es(1,3) = es(1,2)
      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(2,3) = es(1,2)
      es(3,2) = es(2,3)

      do i=lft,llt
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo
      
      return
      end 
