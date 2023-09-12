
! ==============================================================================================
! Program Small_Stars
! Description: This program forms planetary systems around small stars (M and Brown Dwarfs).
! Reference: Miguel, Y. et al. 2020, MNRAS 491, 1998–2009
! ==============================================================================================

implicit none

! ------------------------- Included Files -------------------------
include "parameters.par"  ! Declaration of the variables and constants

! ------------------------- File Operations -------------------------

! Output file below:
 open(330, file="LHS_3154-GaussMstar-Chamaleon.sal")

!Generacion de los satellites iniciales they semimajor axis is chosen random

!y del disco inicial
do ii=1, nsistems

!inicializo algunos valores
  n=1
  time=1.d0           !yr
  mtg0=0.0
  mts0=0.0

!We use the mass of the star as a gaussian as the observations
  mstar=0.1118 !LHS 3154 in solar masses

! In the standard scenario, the disk mass is found as a relation with the stellar mass. 
!Following Williams & Cieza 2011 this relation is: 
  md= 0.01*mstar!solar masses

! the disipation timescale is log-uniform calculated
    y1=dlog10(tau_min)
    y2=dlog10(tau_max)
    yu=dble(ran(semilla))*(y2-y1)+y1
    taugas=10.**(yu)                    !yr

! Some parameters derived from the Mdisk, they are for the planet trapping (check papers by Cridland et al.)
  if (md >= 1.d-4 .and. md < 1d-3)then
    md_max_pt = 5.d-2 !in Earth Masses
    md_min_pt = 5.d-3 !in Earth Masses
  else
    md_max_pt = 0.5d-1 !in Earth Masses
    md_min_pt = 5.d-2 !in Earth Masses
  endif

!Stellar lluminosity
  elestar=0.00119 !for LHS 3154 in solar luminosity

!Location of the iceline
!to get it, I replaced the T for 170K and used the relation between the luminosity and the stellar mass
  rice= 0.075*(mstar/0.1)         !AU new estimation, using Ida+2016 Irr profile.   
  mstar=emesol*mstar              !stellar mass in gr

!Calculus of the inner radius following Ormel+2017.
!This is the initial one, then it will change with time, because Mg^dot changes exponentially.
  Mg_dot_0=1.d-10
  Mg_dot=1.d-10 
  rmin= (((180.**4.*(0.5*rsun)**12.)/(4.*g*mstar*(Mg_dot*emesol/yearsec)**2.))**(1./7.)) !cm
  amin=rmin                 !inner location for the seeds generation ! [cm]

!Calculus of the separation radius between viscous and radiative heating regimes (Ida, Guillot, Morby 2016)
  a_vis_irr= 1.8d0 * elestar**(-20.d0/33.) * (mstar/emesol)**(31.d0/33.) * &
               (alfa/1.d-3)**(-14.d0/33.) * (Mg_dot/1.d-8)**(28.d0/33.) ! AU

!==================================================================
!                         Grid definition
!==================================================================

!For the solids and gas surface density I will use a simple non-exponential profile that goes with a^-1
!They are also the same profiles used in Testi+2016
  sigmag0=md*(2.0d0-gama)*emesol/(2.*pi*(rc*aucm)**2.)    !g/cm2
  sigmad0=sigmag0/Md_ratio                                !gr/cm2 using a gas to dust ratio of 100
  do i=1,ngrilla
    if(i == 1)then
      a_gr(i)=rmin/aucm      !au
      av_gr(i)=rmin          !cm
    else
      if(a_gr(i) < sep) a_gr(i)= a_gr(i-1)+grillajup  !en AU
      if(a_gr(i) >= sep) a_gr(i)= a_gr(i-1)+grillajup_out  !en AU
      av_gr(i)= a_gr(i)*aucm        !en cm
    endif
    if(a_gr(i) <= rice) fice=1.d0
    if(a_gr(i) > rice) fice=2.d0
    sigma_s(i)=sigmad0*fice*(a_gr(i)/rc)**(-gama)* &
      exp(-(a_gr(i)/rc)**(2.-gama))         !Initial solid surface density [g/cm2]
    sigma_g(i)=sigmag0*(a_gr(i)/rc)**(-gama)* &
      exp(-(a_gr(i)/rc)**(2.-gama))         !Initial gas surface density [g/cm2]
    sigma_g_ini(i)=sigma_g(i)  
    if(a_gr(i) < a_vis_irr)then             
      T(i)=200. * (mstar/emesol)**(3.d0/10.) * (alfa/1.d-3)**(-1.d0/5.) * &
           (Mg_dot/1.d-8)**(2.d0/5) * a_gr(i)**(-9.d0/10.) !K (Ida, Guillot, Morby 2016)
    else
      T(i)=150.*elestar**(2./7)*(mstar/emesol)**(-1./7)*a_gr(i)**(-3./7)  !K (Ida, Guillot, Morby 2016)
    endif
    if(a_gr(i)<rice) m_iso_disk_tot= 2.*pi*av_gr(i)*sigma_s(i)*grillajup*aucm + m_iso_disk_tot
    if(a_gr(i) >= rmax/aucm)then

!Calculus of total number of points in the grid
      ind=i
      exit
    endif
    mtg0=mtg0+sigma_g(i)*2*pi*av_gr(i)*grilla !total initial gaseous disc
    mts0=mts0+sigma_s(i)*2*pi*av_gr(i)*grilla !total initial solid disc
    m_iso(i)=2.*pi*av_gr(i)*grilla*sigma_s(i)  !isolation mass in gr

  enddo

!==================================================================
!                         Initial Embryos
!==================================================================
  a(1)=rmin/aucm 
  do i = 1, MaxSed

! Initial position of each seed [cm]
! It is calculated with the prescription of Ida & Lin VI (2010)
    if(i > 1)a(i)=a(i-1)+Delta_a
    acm(i)=a(i)*aucm                         !in cm
    a_ini(i)=a(i)                            !initial embryos location [AU]
    rplanet(i)=size_ini*1.d5                 !radio del embryon de satelite en cm 
    emepla(i)=ro_sat*ctpi*rplanet(i)**3.     !all the embryos have the same initial mass in g
    eme_core(i)=emepla(i)
    erreh(i)=a(i)*(emepla(i)/(3.d0*mstar))**0.33333333 !radio de hill del embryo en au        

!calculus of the isolation mass, for this I have to search the sigma_s closest to the planet location
    if(a(i) < sep )incr=grillajup
    if(a(i) >= sep )incr=grillajup_out
    do kk= 1,ind
        if(dabs(a_gr(kk)-a(i)) < incr/2.d0)then
            m_iso(i)=2.*pi*acm(i)*sigma_s(kk)*10.*erreh(i)*aucm ! the isolation mass in gr
            exit
        endif
    enddo
    Delta_a=10.d0*(m_iso(i)/(3.*mstar))**(1./3.)*a(i)  !au 

  endif

! calculates the initial water percetentage
    r = ((a(i)-rice)*10.)
    w1(i) = ((2./sqrt(pi)*sign(1.d0,r)*sqrt(1.-exp(-r**2.))*(sqrt(pi)/2.+31./200.* &
      exp(-r**2.)-341./8000.*exp(-2.*r**2.)))*25.+25.)/100. !water percentage
    indgas(i)=0    
    in(i)=i
    tot_sed=i

!1) Calculate the initial eccentricity of each planet
!with a log-uniform distribution between 0.001 and 0.01
    semilla=semilla*3
    z1=dlog10(0.001d0)
    z2=dlog10(0.01d0)
    ecc(i)=10.**(dble(ran(semilla))*(z2-z1)+z1)

    if(a(i) >= rmax/aucm)then
      tot_sed=i-1
      in(i)=i-1
      exit
    endif
  enddo

!==================================================================
!                      Accretion Calculation
!==================================================================
 200     continue               ! lazo por paso de tiempo
  if(time > taugas)goto 1111 !if the time is larger than the gas dissipation timescale it doesn't accrete or migrate anymore
!in 1111 check resonance trapping and possible collisions

! Dissipation of gas changes some things-------------------------
  mtg=0.0
  mts=0.0         
  det=exp(-time/taugas)

!     calculus of the gas that remains on the disc. It is disipating due to
!     viscosity is falling into the star (taugas is this time-scale)
!     calculus of the solids that remain in the disc. They dessapera because the embryos
!     accretion and because it falls into Jupiter due to gas drag effect

  do i=1,ind
    sigma_g(i)=sigma_g_ini(i)*det
    if(sigma_g(i) <= 1.d-10)sigma_g(i)=1.d-10 
  enddo

! The inner disk also moves because the gas accretion into the star decreases exponentially
  Mg_dot=Mg_dot_0*det
  if(Mg_dot <= 1.d-30)Mg_dot=1.d-30
  rmin= (((180.**4.*(0.5*rsun)**12.)/(4.*g*mstar*(Mg_dot*emesol/yearsec)**2.))**(1./7.)) !cm

! The radius that separates the viscous and irradiated regime also changes because Mg_dot changes
  a_vis_irr= 1.8d0 * elestar**(-20.d0/33.) * (mstar/emesol)**(31.d0/33.) * &
               (alfa/1.d-3)**(-14.d0/33.) * (Mg_dot/1.d-8)**(28.d0/33.) ! AU
  if(a_vis_irr <= 1.d-10)a_vis_irr=1.d-10

!---
  do k=1,tot_sed
    i=in(k)
    if(a(i) <= rmin/aucm)cycle !if the satellite is inside the inner radius it does not accrete anymore 

!     Solid accretion

!    mean solids in the feeding zone
    zai= a(i)-5.d0*(2.*emepla(i)/(3.*mstar))**(1./3.)*a(i) !AU, mutual Hill radius (Ida-Lin 2010)
    zaf= a(i)+5.d0*(2.*emepla(i)/(3.*mstar))**(1./3.)*a(i) !AU
    cp=0
    nonzero=0
    sigmadp= 0.d0
    sigmagp= 0.d0

! mean density in the feeding zone
    do j= 1,ind
      if(a_gr(j) >= zai.and.a_gr(j) < zaf)then
        cp= cp+1
        if(sigma_s(j) > 1.d-8)nonzero=nonzero+1
        sigmadp=sigma_s(j)+sigmadp
        sigmagp=sigma_g(j)+sigmagp
      endif                        
    enddo
    sigmadp=sigmadp/cp !g cm-2
    sigmagp=sigmagp/cp !g cm-2
    if(cp == 0)then

!     In case the feeding zone is too small, finds a grid point
      if(a(i) < sep )incr=grillajup
      if(a(i) >= sep )incr=grillajup_out
      do kk= 1,ind
        if(dabs(a_gr(kk)-a(i)) < incr/2.d0)then
          sigmadp=sigma_s(kk)
          sigmagp=sigma_g(kk)
          exit
        endif
      enddo
    endif

    eme_sca(i)=2.*3.16*mstar*rplanet(i)/acm(i) !stop accreting mass in g
    if(emepla(i) >= eme_sca(i))then
      goto 3333        !stops accretion (Ida&Lin 2004,2008)
    endif
    if(sigmadp <= 1.e-10)then
      goto 3333        !It stops accretion if there are no more solids
    endif

!=======================================================================

! standard accretion of Ida & Lin 2004
    tacc1(i)=1.2e5*(10./sigmadp)*(a(i)**0.5)*(emepla(i)/emet)**(1./3.)*(mstar/emesol)**(-1./6.)* &
        ((sigmagp/2.4e3)**(-1./5.)*(a(i)**(1./20.))*(1.e20/1.e18)**(1./15.))**2       
    dmdt=emepla(i)/tacc1(i)
    emedot_crit=dmdt*1.d6/emet !coefficient that we will use for the critical mass
    dm=dmdt*dt       !gr

! Estimates the total mass in the FZ to know how much mass it can accrete
    max_mass_fz=sigmadp*pi*((zaf*aucm)**2.-(zai*aucm)**2.)

! If the mass that tries to accrete is larger than the one available, then it accretes the available mass
    if(dm > max_mass_fz) dm = max_mass_fz

! We remove from the disk what the planet accreted
    if(a(i) < sep )incr=grillajup
    if(a(i) >= sep )incr=grillajup_out
    if(cp == 0)then 
      do kk= 1,ind-1
        if(dabs(a_gr(kk)-a(i)) < incr/2.d0)then 
          sigma_s(kk)=sigma_s(kk)-  &
            (dm/(pi*(av_gr(kk+1)**2.-av_gr(kk)**2.)))
          if(sigma_s(kk) <= 1.d-10)then                        
            sigma_s(kk)=1.d-10
          endif
          exit
        endif
      enddo
    else
      do jj=1,ind-1
        if(a_gr(jj) > zaf)exit
        if(a_gr(jj) >= zai.and.a_gr(jj) <= zaf)then
          sigma_s(jj)=sigma_s(jj)- &
            (dm/(pi*(av_gr(jj+1)**2.-av_gr(jj)**2.)*nonzero))
          if(sigma_s(jj) <= 1.d-10)then
            sigma_s(jj)=1.d-10                                     
          endif
        endif
      enddo
    endif
    r = ((a(i)-rice)*10.)
    w2 = ((2./sqrt(pi)*sign(1.d0,r)*sqrt(1.-exp(-r**2.))*(sqrt(pi)/2.+31./200.*exp(-r**2.)-341./8000.*exp(-2.*r**2.)))*25.+25.)/100. !water percentage
    w_int = (w1(i)*eme_core(i) + w2*dm)/(eme_core(i)+dm)
    w1(i) = w_int
    eme_core(i)=eme_core(i)+dm
    emepla(i)=emepla(i)+dm
    rplanet(i)=(emepla(i)/(ctpi*ro_sat))**0.3333333 !cm

!==================================================================
!                          Gas Accretion
!==================================================================
    em=emepla(i)/emet           !em is the mass in Earth masses
    crit_gas(i)=10.d0*emedot_crit**0.25 !critical mass (IL04, Ikoma2000)
    if(em > crit_gas(i))then
      if(sigmagp > 1.e-20)then !It doesn't accrete gas if there isn't gas
        tkh=1.d3*(emepla(i)/(100.d0*emet))**(-3.)         !Ida&Lin, last paper 2013
        tasa=em/tkh !Earth mass per year

! Now I calculate the limit to gas accretion (it can't be larger than the gas accretion onto the star, Ida+2013)
        nu= alfa*(0.05*a(i)**0.25*acm(i))**2 * &
            ((g*emestar/acm(i)**3.)**0.5)
        m_dot_disk= 3.*pi*sigmagp*nu
        m_dot_disk= m_dot_disk/det
        if (tasa*emet > m_dot_disk)tasa = m_dot_disk/emet
        em=em+tasa*dt !Earth mass

! Thermal condition to stop gas accretion (locally) (ida & lin 2004)
        eme_th=0.95d3*a(i)**(3./4.)*elestar**(3./8.)*(emesol/mstar)**0.5 !Earth masses

! Global condition to stop the gas accretion (ida & lin 2004 - papers I to III)
        eme_glob=pi*acm(i)**2.*sigmagp/emet !Earth masses
        if(em < eme_th .or. em < eme_glob)then

! calculus of the new mass and Hill radius of the planet
          emegas(i)=emegas(i)+tasa*dt !Earth masses
          emepla(i)=em*emet    !gr
          erreh(i)=a(i)*(emepla(i)/(3.d0*mstar))**0.33333333
        endif
      endif
    endif
 3333       continue

!New planetary radius after the gas accretion
    rplanet(i)=(emepla(i)/(ctpi*ro_sat))**0.3333333 !cm
  enddo

!=======================================================================
!=======================================================================
!Dynamical effects 
!=======================================================================

!Planetesimal's migration due to gas drag effect
  call migra_planetes(ind,a_gr,sigma_s,T,sigma_g,av_gr,time,dt,mstar) 

!=======================================================================

!Planetary type I migration
   call migrarI(a,in,det,emepla,acm,elestar,&
               mstar,sigmag0,time,tot_sed,rmin,a_vis_irr,Mg_dot,dt,gama,rc) 

!=======================================================================

!Planetary type II migration
  rm=10.d0*exp(2.*time/taugas) !AU (Miguel+2011a)
  call migrar(a,rm,in,acm,rplanet,ind,a_gr,elestar, &
              emepla,det,av_gr,sigma_g,time,taugas,sigmag0, &
              mstar,tot_sed,rmin,a_vis_irr,Mg_dot,dt,gama,rc)

!=======================================================================
1111 continue
!Resonance trapping
  call resonance_trapping(a,in,acm,emepla,mstar,tot_sed)

  !=======================================================================
  if (time <= taugas)then ! the trapping only works if there is still gas in the disk
   call migration_trapping(a,in,acm,a_gr,emepla,av_gr,mstar,&
     tot_sed,rmin,a_vis_irr,md_max_pt,md_min_pt,rice,a_ini,erreh,&
     Mg_dot,elestar)
  endif
  do k=1,tot_sed
      i=in(k)
      if(a(i) < rmin/aucm)then
        a(i) = rmin/aucm
        acm(i)=a(i)*aucm
        erreh(i)=a(i)*(emepla(i)/(3.d0*mstar))**0.33333 !radio de hill del planeta 
      endif
  enddo

!=======================================================================
!updates time
  if(tot_sed>1)call indexx(tot_sed,a,in) ! order the embryos according to a
  time=time+dt 
  if(time < tfin)goto 200
  do k=1,tot_sed
    i=in(k)   
    write(330,'(15(e15.4,1x))')a(i),emepla(i)/emet,emegas(i),m_iso(i)/emet, &
        alfa,taugas,md,mstar/emesol,rc,a_ini(i),w1(i)*100.,ecc(i),rmin/aucm,a_vis_irr,rice
  enddo
  call flush(20)
enddo
stop
end

!==================================================================
!==================================================================
!==================================================================
!                          SUBROUTINES
!==================================================================
!==================================================================
!==================================================================

!this routine traps the planets in the ice line and in the place
!where the viscosity and irradiation regimes change
!(Ver Alex Cridland+2018)
subroutine migration_trapping(a,in,acm,a_gr,emepla,av_gr,mstar,&
  tot_sed,rmin,a_vis_irr,md_max_pt,md_min_pt,rice,a_ini,erreh,&
  Mg_dot,elestar)
implicit none
include "parameters.par"  !declaration of the variables and constants
do jj=1,tot_sed
    i=in(jj)

! Check first if its in the right mass range for the appropiate Mdisk
    if(emepla(i)/emet < md_max_pt .and. md_min_pt < emepla(i)/emet)then

! Calculate the eme_gap, the planet is only trapped in type I mig.
      if(a(i) < a_vis_irr)then             
        Tcal=200. * (mstar/emesol)**(3.d0/10.) * (alfa/1.d-3)**(-1.d0/5.) * &
          (Mg_dot/1.d-8)**(2.d0/5) * a(i)**(-9.d0/10.) !K (Ida, Guillot, Morby 2016)
      else
        Tcal=150.*elestar**(2./7)*(mstar/emesol)**(-1./7)*a(i)**(-3./7)
      endif  
      omega1=(g*mstar/acm(i)**3.)**0.5               
      eme_gap1=80.*alfa*kb*Tcal/(me*mu*omega1**2.*acm(i)**2.)*mstar
      if(emepla(i) > eme_gap1) cycle        !no es atrapado si migra por migII

!If the satellite has an initial semimajor axis larger that rice and now is below it its trapped there
      if (rice < a_vis_irr)then
        if(a_ini(i) >  a_vis_irr)then
          if(a(i)<a_vis_irr) a(i) = a_vis_irr
        elseif(a_ini(i) <  a_vis_irr .and. a_ini(i) > rice)then
           if(a(i) < rice) a(i) = rice
        elseif(a_ini(i) < rice)then
          if(a(i) < rmin/aucm) a(i) = rmin/aucm
        endif
      else
        if(a_ini(i) > rice)then
          if(a(i) < rice) a(i) = rice
        elseif(a_ini(i) <  rice .and. a_ini(i) > a_vis_irr)then
           if(a(i) < a_vis_irr) a(i) = a_vis_irr
        elseif(a_ini(i) < a_vis_irr)then
          if(a(i) < rmin/aucm) a(i) = rmin/aucm
        endif
      endif
    else
      if(a(i) < rmin/aucm) a(i) = rmin/aucm
    endif
    acm(i)=a(i)*aucm
    erreh(i)=a(i)*(emepla(i)/(3.d0*mstar))**0.33333 
  enddo
  return
  end

!==================================================================
!==================================================================

!(Ver Ida y Lin, The Astrophysical Journal, 604:388, 2004)
subroutine migrar(a,rm,in,acm,rplanet,ind,a_gr,elestar, &
  emepla,det,av_gr,sigma_g,time,taugas,sigmag0,mstar,&
  tot_sed,rmin,a_vis_irr,Mg_dot,dt,gama,rc)
implicit none
include "parameters.par"  !declaration of the variables and constants

!type II Migration in the disk-dominated regime(Hasegawa & Ida, 2013)
rmcm=rm*aucm                       !rm in cm
sigmag=sigmag0/rm         
sigmag=sigmag0*(rm/rc)**(-gama)* &
  exp(-(rm/rc)**(2.-gama))         !Initial gas surface density [g/cm2]
sigmag=sigmag*det
om=dsqrt(g*mstar/(rmcm)**3.) 
if(rm < a_vis_irr)then             
  Tcalrm=200. * (mstar/emesol)**(3.d0/10.) * (alfa/1.d-3)**(-1.d0/5.) * &
  (Mg_dot/1.d-8)**(2.d0/5) * rm**(-9.d0/10.) !K (Ida, Guillot, Morby 2016)
else
  Tcalrm=150.*elestar**(2./7)*(mstar/emesol)**(-1./7)*rm**(-3./7) !K (Ida, Guillot, Morby 2016)

endif
hm=(2.*kb*Tcalrm/(mu*me*om**2.))**0.5!cm
do j=1,tot_sed
  i=in(j)
  if(acm(i) > rmin)then   !it stops migration if it reaches the iner boundary 
    if(a(i) < a_vis_irr)then             
      Tcal=200. * (mstar/emesol)**(3.d0/10.) * (alfa/1.d-3)**(-1.d0/5.) * &
        (Mg_dot/1.d-8)**(2.d0/5) * a(i)**(-9.d0/10.) !K (Ida, Guillot, Morby 2016)
    else
      Tcal=150.*elestar**(2./7)*(mstar/emesol)**(-1./7)*a(i)**(-3./7)
    endif  
    omega1=(g*mstar/acm(i)**3.)**0.5            
    eme_gap1=80.*alfa*kb*Tcal/(me*mu*omega1**2.*acm(i)**2.)*mstar 
    if(emepla(i) > eme_gap1)then !starts type II migration if it opened up a gap 

!condition for one regime or the other 
      eme_disk=0.0
      do kk= 1,ind
        if(a_gr(kk) <= a(i) .and. a(i) < sep)then
          eme_disk=eme_disk+ &
            2.*pi*av_gr(kk)*grilla*sigma_g(kk)
        elseif(a_gr(kk) <= a(i) .and. a(i) >= sep)then
          eme_disk=eme_disk+ &
            2.*pi*av_gr(kk)*grilla_out*sigma_g(kk)
        endif
      enddo
      if(emepla(i).lt.eme_disk)then !disk-dominated type II mig regime
        tau2=a(i)/100.*taugas !in years
      else             !satellite-dominated type II mig regime
        sg=(rm-a(i))/dabs(rm-a(i))
        tau2=1./(3.*sg*alfa*sigmag*rmcm**2.*om/(emepla(i)*omega1)* &
          (hm/acm(i))**2.*om) ! in sec
        tau2=tau2/yearsec !in years
      endif
      dadt=a(i)/tau2   !AU/years
      dadt=dadt*cmigII 
      da=dadt*dt       !AU
      a(i)=a(i)-dadt*dt
      acm(i)=a(i)*aucm
    endif
  endif
enddo
return
end

!==================================================================
!==================================================================
subroutine migrarI(a,in,det,emepla,acm,elestar,&
    mstar,sigmag0,time,tot_sed,rmin,a_vis_irr,Mg_dot,dt,gama,rc)
implicit none
include "parameters.par"  !declaration of the variables and constants
do j=1,tot_sed
  i=in(j)

  if(acm(i) > rmin)then   ! It doesn't migrate if it reached the inner boundary 
    if(a(i) < a_vis_irr)then             
        Tcal=200. * (mstar/emesol)**(3.d0/10.) * (alfa/1.d-3)**(-1.d0/5.) * &
          (Mg_dot/1.d-8)**(2.d0/5) * a(i)**(-9.d0/10.) !K (Ida, Guillot, Morby 2016)
    else
        Tcal=150.*elestar**(2./7)*(mstar/emesol)**(-1./7)*a(i)**(-3./7)
    endif  
    omega1=(g*mstar/acm(i)**3.)**0.5               
    eme_gap1=80.*alfa*kb*Tcal/(me*mu*omega1**2.*acm(i)**2.)*mstar
    if(emepla(i) > eme_gap1) cycle        !it doesn't migrate with type I if it migrates with type II 

!gas density in the planet's location
    sigmag=sigmag0*(a(i)/rc)**(-gama)* &
        exp(-(a(i)/rc)**(2.-gama))         !Initial gas surface density [g/cm2]
    sigmag=sigmag*det

    cs=5./3.*Tcal*(kb/(mu*me)) !this is the square of cs
    tau_migI=1./(2.7+1.1)*cs/(acm(i)*omega1)**2.* &
        (mstar/emepla(i))* &
        (mstar/(acm(i)**2.*sigmag))/(omega1*31536000.) 
    dadt=acm(i)/tau_migI
    dadt=-dadt*constmigI
    da=dadt*dt               
    acm(i)=acm(i)+dadt*dt !cm
    a(i)=acm(i)/aucm !AU 
  endif
enddo
end

!==================================================================
!==================================================================
subroutine migra_planetes(ind,a_gr,sigma_s,T,sigma_g,av_gr, &
  time,dt,mstar)
implicit none
include "parameters.par"  !declaration of the variables and constants

!calculus from Thommes et al (2003) and Mosqueira & Estrada(2003), with changes according to the stellar mass

do j=1,ind-1
  tau1=1.e6*(sat_disc_size)*(100./T(j))**(3./2.)*(1.d4/sigma_g(j))*(mstar/(0.1*emesol)) !years 
  dsigmadt_n(j)=(sigma_s(j)/tau1)!gr/cm2/yr
enddo
do j=1,ind-1
  sigma_s(j)=sigma_s(j)-dsigmadt_n(j)*dt+dsigmadt_n(j+1)*dt
  if(sigma_s(j).lt.1d-10)sigma_s(j)=0.0
enddo
return
end

!==================================================================
!==================================================================
subroutine resonance_trapping(a,in,acm,emepla,mstar,tot_sed)
implicit none
include "parameters.par"  !declaration of the variables and constants

  do j=1,tot_sed
    i=in(j)
    difk_ch=150.        !random value to start  
    do kk=1,tot_sed
      k=in(kk)
      if(k.ne.i)difk(k)=a(i)-a(k)
      if(difk(k) > 0.0 .and. difk(k) < difk_ch)then
          difk_ch=difk(k)
          numero=k 
      endif
    enddo
    if(difk_ch < 150.)then
      dif=abs(a(i)-a(numero))
      rhill_sum=((2.*emepla(i))/(3.*mstar))**(1./3.)*a(i)
      if(dif < 5.*rhill_sum)then
        if(a(i) < a(numero))then
          algo=a(i)+dif/2.
          a(i)=algo-rhill_sum*5./2.
          a(numero)=algo+rhill_sum*5./2.
          acm(i)=a(i)*aucm
          acm(numero)=a(numero)*aucm
        else
          algo=a(numero)+dif/2.
          a(i)=algo+rhill_sum*5./2.
          a(numero)=algo-rhill_sum*5./2.    
          acm(i)=a(i)*aucm
          acm(numero)=a(numero)*aucm
        endif
      endif
    endif
  enddo
return
end


!==================================================================
!==================================================================
!====================================
      DOUBLE PRECISION FUNCTION RND(DSEED)
      DOUBLE PRECISION DSEED,D2A32
      D2A32=2.D0**32
      DSEED=DMOD(3125.*DSEED,D2A32)
      RND=DSEED/D2A32
      RETURN
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software

!=====================================
!==================================================================
!==================================================================

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
      INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
      IF(L.GT.1)THEN
      L=L-1
      INDXT=INDX(L)
      Q=ARRIN(INDXT)
      ELSE
      INDXT=INDX(IR)
      Q=ARRIN(INDXT)
      INDX(IR)=INDX(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
      INDX(1)=INDXT
      RETURN
      ENDIF
      ENDIF
      I=L
      J=L+L
20    IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
      IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
      ENDIF
      IF(Q.LT.ARRIN(INDX(J)))THEN
      INDX(I)=INDX(J)
      I=J
      J=J+J
      ELSE
      J=IR+1
      ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
      GO TO 10
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software

!==================================================================
!==================================================================
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
          IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
          NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software

!==================================================================
!==================================================================
FUNCTION ran(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ran
INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
  am=nearest(1.0,-1.0)/IM
  iy=ior(ieor(888889999,abs(idum)),1)
  ix=ieor(777755555,abs(idum))
  idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0) iy=iy+IM
ran=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran

!C  (C) Copr. 1986-92 Numerical Recipes Software