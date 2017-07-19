!             Program block bhf-matrix.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95  
!
!                    Begin bhf-matrix code
!
SUBROUTINE  setup_hfmatrix
  USE constants
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  IMPLICIT NONE
  INTEGER :: number_orbits
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: coeffs  

  !     reserve space in memory for various arrays
  ALLOCATE ( ra (n_rel), wra (n_rel));  ALLOCATE ( rgkk (n_rel), wgkk (n_rel))
  ALLOCATE ( hol (n_rel, 0:lmax, 0:nmax) )
  !     set up mesh points in lab frame
  CALL rel_mesh                   
  !     setup ho wave functions
  CALL ho_wfunction                
  !    Expansion coefficients for sp wave functions. 
  ALLOCATE(coeffs(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE ( bhf_hol (n_rel,all_orbit%total_orbits))
  ALLOCATE ( wave_function (n_rel,all_orbit%total_orbits))
  coeffs = 0.0_dp; wave_function = 0.0_dp
  !  perform the HF calculation
  CALL brueckner_hartree_fock(coeffs,all_orbit%total_orbits)
  !  update the g-matrix
  CALL setupg_bhf(coeffs,all_orbit%total_orbits)
  SELECT CASE (type_of_interaction)
  CASE ('core-diagrams')
     CALL core_diagrams
  CASE ('open-diagrams')
     IF (n_body_int == 'onebody') CALL onebody_contribution
  END SELECT
  DEALLOCATE(coeffs) ;    DEALLOCATE ( bhf_hol,wave_function )
  DEALLOCATE ( rgkk, wgkk, ra, wra) ; DEALLOCATE ( hol)

END SUBROUTINE setup_hfmatrix
!
!           Set up the BHF G-mtx in the lab-frame 
!
SUBROUTINE setupg_bhf(coeffs,ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE hoosc_gmatrix
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs
  REAL(DP), ALLOCATABLE :: bhf_coeff(:,:)
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  INTEGER ::  p_parity, ang_mom, isospin_z, n_confs, ie, bra, ket
  REAL(DP), ALLOCATABLE :: gna(:,:,:),  temp(:,:,:)

  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations 
           CALL number_gmatrix_confs &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           n_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(n_confs+n_confs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,gmatrix_configs)
           ALLOCATE(temp(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(gna(n_confs, n_confs,n_startenergy_g))
           ALLOCATE(bhf_coeff(n_confs, n_confs))
           gna=0.0_dp; bhf_coeff = 0.0_dp; temp = 0.0_dp
           CALL fetch_matrix(isospin_z,p_parity,ang_mom,gna,gmatrix_configs)
           CALL bhf_coefficients(ang_mom,isospin_z,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
           DO ie=1,n_startenergy_g
              temp(:,:,ie) = MATMUL(gna(:,:,ie),TRANSPOSE(bhf_coeff(:,:)))
              gna(:,:,ie) = MATMUL(bhf_coeff(:,:),temp(:,:,ie))
           ENDDO
           !     update table of matrix elements 
           CALL update(isospin_z,p_parity,ang_mom,gna,gmatrix_configs)
           !     free space in heap
           DEALLOCATE(gmatrix_configs%config_ab)
           DEALLOCATE(bhf_coeff)
           DEALLOCATE(gna, temp)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE setupg_bhf
!
!                    Update out G-mtx
!
SUBROUTINE update(it,ip,ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, ijd, it, ip, ij, ia, ib, ic,id,ie
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(IN)  :: gna

  ijd=ij+ij
  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL replace_g(ia,ib, ic, id, it, ip, ij ,gna(j,i,:))
        WRITE(9,'(7I4,5X,10(5X,E12.6))') it, ip, ijd, ia, ib, ic, id, (gna(j,i,ie),ie=1,n_startenergy_g)
     ENDDO
  ENDDO

END SUBROUTINE update
!
!                    Get just the G-mtx to be used in various HF iterations
!
SUBROUTINE fetch_matrix(it,ip,ij,gna,gmatrix_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE relcm_gmatrix
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER :: i,j, it, ip, ij, ia, ib, ic,id
  REAL(DP) :: norm, dij
  REAL(DP), DIMENSION(gmatrix_configs%number_confs, &
       gmatrix_configs%number_confs,n_startenergy_g), INTENT(INOUT)  :: gna
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  DO i=1,gmatrix_configs%number_confs
     ia=gmatrix_configs%config_ab(i*2-1)
     ib=gmatrix_configs%config_ab(i*2)
     DO j=i,gmatrix_configs%number_confs
        ic=gmatrix_configs%config_ab(j+j-1)
        id=gmatrix_configs%config_ab(j+j)
        CALL pphhmtx(ia,ib,ic,id,ij,ans)
        norm=1.0_dp/dij(ia,ib)/dij(ic,id)
        gna(j,i,:)=ans(:)*norm
        gna(i,j,:) = ans(:)*norm
     ENDDO
  ENDDO

END SUBROUTINE fetch_matrix
!
!                 Set up h.o. wf for rel cm system and lab frame
!                 It computes also a complex
!
SUBROUTINE ho_wfunction
  USE constants
  USE wave_functions
  IMPLICIT NONE
  INTEGER :: n, l, i, j
  REAL(DP) :: ph, sum_rel
  REAL(DP)  :: cx(0:200), factor, z_lab,xp

  DO n=0,nmax
     ph=(-1.D0)**n
     DO l=0,lmax
        factor = 0.5D0*((n+l+2)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        sum_rel=0.0_dp
        DO i=1,n_rel
           !  real ho wave function
           z_lab= ra(i)*ra(i)*oscl*oscl; cx = 0.0_dp
           CALL laguerre_general( n, l+0.5D0, z_lab, cx )
           xp = EXP(-z_lab*0.5D0)*((ra(i)*oscl)**l)*cx(n)
           hol(i,l,n) = xp*ph*factor*(oscl**(1.5D0)) ! lab wf
           sum_rel=sum_rel+ wra(i)*(hol(i,l,n)*ra(i))**2
        ENDDO
        WRITE(6,'(21H Norm cm ho wf n,l : ,2I3,2X,3F12.7)') n, l, sum_rel
     ENDDO
  ENDDO

END SUBROUTINE ho_wfunction
!
!     Brueckner-Hartree-Hock self-consistency
!
SUBROUTINE brueckner_hartree_fock(coeffs, ncoeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: sp_configs
  INTEGER, INTENT(IN) :: ncoeffs
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(INOUT)  :: coeffs  
  INTEGER :: a, c, ket , bra, nrot, hole_states, hf_iter, max_hfiter
  INTEGER :: la, lc, na, nc, k, i
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_vect
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: hf_mtx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hf_eigen
  REAL(DP), DIMENSION(n_rel) ::  sum_wf
  REAL(DP) :: e_kin, hf, sigma
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: int_factor, kin_energy
  ALLOCATE (int_factor(n_rel), kin_energy(n_rel) )
  int_factor(:)=wra(:)*ra(:)*ra(:)

  !     initialize bhf harmonic oscillator sp wave function
  bhf_hol=0.0_dp; max_hfiter = 100
  ALLOCATE (hf_mtx(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_vect(all_orbit%total_orbits,all_orbit%total_orbits))
  ALLOCATE (hf_eigen(all_orbit%total_orbits))
  hf_vect = 0.0_dp
  !  for first iteration coeffs has only diagonal values equal 1
  DO bra = 1, all_orbit%total_orbits
     coeffs(bra,bra) = 1.0_dp
  ENDDO
  hf_iter = 0; sigma = 1.0_dp
  DO WHILE ((hf_iter <= max_hfiter).AND.(ABS(sigma) > 1E-8))
     !     set up bhf matrix to be diagonalized
     hf_mtx = 0.0_dp; hf_vect = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        DO ket= bra, all_orbit%total_orbits
           c = ket
           IF(all_orbit%ll(a) /= all_orbit%ll(c)) CYCLE
           IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
           IF(all_orbit%itzp(a) /= all_orbit%itzp(c) ) CYCLE
           hf = 0.0_dp               
           SELECT CASE (physical_system) 
           CASE('nuclear_physics')
              kin_energy(:)=0.5*ra(:)*ra(:)*hbarc*hbarc/p_mass(all_orbit%itzp(c))
           CASE('atomic_physics')
              kin_energy(:)=0.5*ra(:)*ra(:)
           END SELECT
           CALL  diagram_HF(a, c, coeffs, hf_iter, hf, ncoeffs)
           !  compute the kinetic energy or unperturbed one-body H
           e_kin = 0.0_dp
           SELECT CASE (type_of_renormv)
           CASE ('no-core')
              IF ( a == c) THEN
                 DO k=1,n_rel
                    e_kin=e_kin+2*hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                         hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                         int_factor(k)*kin_energy(k)
                 ENDDO
              ENDIF
           CASE ('v-nrg')
              IF ( a == c) THEN
                 DO k=1,n_rel
                    e_kin=e_kin+2*hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                         hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                         int_factor(k)*kin_energy(k)
                 ENDDO
              ENDIF
           CASE ('vlowk')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           CASE ('v-krg')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           CASE ('g-matrix')
              DO k=1,n_rel
                 e_kin=e_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))* &
                      hol(k,all_orbit%ll(c),all_orbit%nn(c))* &
                      int_factor(k)*a_factor*kin_energy(k)
              ENDDO
           END SELECT
           hf_mtx(bra,ket) = e_kin+hf
           hf_mtx(ket,bra) = hf_mtx(bra,ket)
        ENDDO
     ENDDO
     !     obtain the BHF coefficients 
     CALL matrix_diag(hf_mtx,all_orbit%total_orbits, all_orbit%total_orbits,&
          hf_eigen, hf_vect,nrot)
     !     set up the new bhf harmonic  oscillator wave function
     !     Note memory stride
     sigma = 0.0_dp
     DO bra = 1, all_orbit%total_orbits
        a = bra
        la = all_orbit%ll(a)
        na = all_orbit%nn(a)
        sum_wf = 0.0_dp
        DO ket= 1, all_orbit%total_orbits
           c = ket
           lc = all_orbit%ll(c)
           nc = all_orbit%nn(c)
           sum_wf(:)=sum_wf(:)+hf_vect(ket,bra)*hol(:,lc,nc)
           coeffs(a,c) = hf_vect(bra,ket)
        ENDDO
        bhf_hol(:,a)=sum_wf(:)
        ! set up of new single-particle energy 
        sigma = sigma +ABS(all_orbit%e(a) - hf_eigen(bra))
        all_orbit%e(a) = hf_eigen(bra)                    
     ENDDO
     sigma = sigma/all_orbit%total_orbits
     WRITE(6,*) 'iteration nr and sigma', hf_iter, sigma
     DO i=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(i) /= 'hole') CYCLE         
        WRITE(6,'(7HNumber:,6(I3,2X),2X,E20.10)') i, all_orbit%nn(i), all_orbit%ll(i), &
             all_orbit%jj(i), &
             all_orbit%nshell(i), all_orbit%itzp(i), all_orbit%e(i)
     ENDDO
     hf_iter = hf_iter+1
  ENDDO
  DO i=1, all_orbit%total_orbits
     IF ( keep_originalenergies == 'no') THEN 
        IF (all_orbit%model_space(i) == 'inside')  THEN
           all_orbit%evalence(i) = all_orbit%e(i)
        ELSE
           all_orbit%evalence(i) = 0.0_dp
        ENDIF
     ELSEIF ( keep_originalenergies == 'yes') THEN
        all_orbit%e(i) = all_orbit%e_original(i)
     ENDIF
     WRITE(8,'(7HNumber:,5(I4,2X),2X,E12.6,2X,E12.6,2X,A10,2X,A10)') i, all_orbit%nn(i), all_orbit%ll(i), &
          all_orbit%jj(i), all_orbit%itzp(i), all_orbit%e(i), all_orbit%evalence(i), &
          all_orbit%orbit_status(i), all_orbit%model_space(i)
  ENDDO
  DEALLOCATE ( hf_mtx)
  DEALLOCATE ( hf_eigen)
  DEALLOCATE ( hf_vect)
  !     new harmonic oscillator wave function
  wave_function=bhf_hol
  DEALLOCATE (int_factor, kin_energy )

END SUBROUTINE brueckner_hartree_fock
!
!    The Hartree-Fock diagram
!
SUBROUTINE diagram_HF(a,c,coeffs, iteration,onebody_diagram_HF, ncoeffs)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c, iteration, ncoeffs
  INTEGER :: j_min, j_max, jph, h, h1, h2, i
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP) :: val, ang_mom_factor, w
  REAL(DP), INTENT(INOUT) :: onebody_diagram_HF
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  onebody_diagram_HF=0.0_dp
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO h1=1, all_orbit%total_orbits
        IF  (all_orbit%jj(h) /= all_orbit%jj(h1)) CYCLE
        IF  (all_orbit%ll(h) /= all_orbit%ll(h1)) CYCLE
        IF  (all_orbit%itzp(h) /= all_orbit%itzp(h1)) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF  (all_orbit%jj(h) /= all_orbit%jj(h2)) CYCLE
           IF  (all_orbit%ll(h) /= all_orbit%ll(h2)) CYCLE
           IF  (all_orbit%itzp(h) /= all_orbit%itzp(h2)) CYCLE
           SELECT CASE (type_of_renormv)
           CASE ('no-core')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h2,h)*coeffs(h1,h)
              ENDDO
           CASE ('v-nrg')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h2,h)*coeffs(h1,h)
              ENDDO
           CASE ('vlowk')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           CASE ('v-krg')
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 onebody_diagram_HF=onebody_diagram_HF+ans(1)*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           CASE ('g-matrix')
              w = all_orbit%e(h2)+(all_orbit%e(c)+all_orbit%e(a))/2
              DO jph=j_min,j_max
                 ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
                 CALL pphhmtx(a,h1,c,h2,jph,ans); IF ( ans(1) == 0.0_dp ) CYCLE
                 CALL interpolate(w,e_start_g,ans,val)
                 onebody_diagram_HF=onebody_diagram_HF+val*ang_mom_factor*coeffs(h1,h)*coeffs(h2,h)
              ENDDO
           END SELECT
        ENDDO
     ENDDO
  ENDDO

END  SUBROUTINE diagram_HF

SUBROUTINE bhf_coefficients(ang_mom,itz,n_confs,ncoeffs,gmatrix_configs,bhf_coeff,coeffs)
  USE single_particle_orbits
  USE configurations
  USE constants
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: gmatrix_configs
  INTEGER, INTENT(IN) :: ncoeffs, n_confs, ang_mom, itz
  REAL(DP), DIMENSION(ncoeffs,ncoeffs), INTENT(IN)  :: coeffs  
  REAL(DP), DIMENSION(n_confs,n_confs), INTENT(INOUT)  :: bhf_coeff
  REAL(DP) :: dij, fnorm
  REAL(DP) :: e_dir, e_exc, e
  INTEGER :: ket , bra, a,b, la, ja, lb, jb, p,q, lq, lp, jp, jq

  fnorm = 0.d0
  DO bra = 1, n_confs
     a = gmatrix_configs%config_ab(bra+bra-1)
     b = gmatrix_configs%config_ab(bra+bra)
     la=all_orbit%ll(a); jb=all_orbit%jj(b)
     ja=all_orbit%jj(a); lb=all_orbit%ll(b)
     DO ket = 1, n_confs
        p = gmatrix_configs%config_ab(ket+ket-1)
        q = gmatrix_configs%config_ab(ket+ket)
        lp = all_orbit%ll(p); jp=all_orbit%jj(p)
        lq = all_orbit%ll(q); jq=all_orbit%jj(q)
        IF ( itz /= 0) THEN
           fnorm = 1.d0/ dij(a,b) /  dij(p,q)
        ELSEIF( itz == 0 ) THEN
           fnorm = 1.d0
        ENDIF
        e_dir=0.0_dp; e_exc=0.D0; e = 0.0_dp
        ! direct term
        IF ( (la == lp).AND.( lb == lq ).AND.( ja == jp ).AND.( jb == jq )) THEN
           e_dir = coeffs(p,a)*coeffs(q,b)
        ENDIF
        ! exchange term
        IF ( (la == lq).AND.( lb == lp ).AND.( ja == jq ).AND.( jb == jp ) .AND. itz /= 0 ) THEN
           e_exc = coeffs(p,b)*coeffs(q,a)
        ENDIF
        e= (e_dir-e_exc*((-1.0_dp)**((2*ang_mom-jp-jq)/2)))*fnorm
        bhf_coeff(bra,ket) = e
     ENDDO
  ENDDO

END SUBROUTINE bhf_coefficients

!
!     This function returns the matrix for V
!
SUBROUTINE pphhmtx(ja,jb,jc,jd,jt,gmtpn)
  USE stored_bare_interaction
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  LOGICAL TRIAG
  INTEGER, INTENT(IN) :: ja,jb,jc,jd,jt
  INTEGER :: iph
  REAL(DP), INTENT(INOUT), DIMENSION(n_startenergy_g) :: gmtpn
  REAL(DP) :: dij, delta, xxx

  gmtpn=0.
  IF(2*all_orbit%nn(ja)+all_orbit%ll(ja) +2*all_orbit%nn(jb)+all_orbit%ll(jb) > nlmax ) RETURN
  IF(2*all_orbit%nn(jc)+all_orbit%ll(jc) +2*all_orbit%nn(jd)+all_orbit%ll(jd) > nlmax ) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  IF((-1)**(all_orbit%ll(ja)+all_orbit%ll(jb)) /=  &
       (-1)**(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((ja == jb).AND.(MOD(jt,2)/=0)) RETURN       
  IF((jc == jd).AND.(MOD(jt,2)/=0)) RETURN       
  IF(triag(all_orbit%jj(ja),all_orbit%jj(jb),2*jt)) RETURN
  IF(triag(all_orbit%jj(jc),all_orbit%jj(jd),2*jt)) RETURN
  CALL mtx_elements(ja,jb,jc,jd,jt,gmtpn)

END SUBROUTINE pphhmtx
!
!     Calculates the crosscoupled matrix element type 1
!
SUBROUTINE cross_coupled_mtxel1(ja,jb,jc,jd,jtot,cross)
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE 
  REAL(DP),  DIMENSION (n_startenergy_g) :: ans
  INTEGER :: ja,jb,jc,jd,jtot,jbra_min, jket_min, jbra_max ,jket_max, &
       jt, j_min, j_max, iph
  REAL(DP) :: fnorm, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_g), INTENT(OUT) :: cross

  cross=0.
  fnorm=SQRT(2.*jtot+1.)*iph((all_orbit%jj(ja)+all_orbit%jj(jd))/2+jtot)
  IF(iph(all_orbit%ll(ja)+all_orbit%ll(jb)) & 
       /=iph(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  jbra_min=ABS(all_orbit%jj(ja)-all_orbit%jj(jb))
  jbra_max=all_orbit%jj(ja)+all_orbit%jj(jb)
  jket_min=ABS(all_orbit%jj(jc)-all_orbit%jj(jd))
  jket_max=all_orbit%jj(jc)+all_orbit%jj(jd)
  j_max=MIN(jbra_max,jket_max)/2
  j_min=MAX(jbra_min,jket_min)/2
  IF(j_min > j_max) RETURN
  DO jt=j_min,j_max
     ang_mom_factor=sjs(all_orbit%jj(jc),all_orbit%jj(ja),2*jtot, &
          all_orbit%jj(jb),all_orbit%jj(jd),2*jt)* &
          (2.*jt+1.)*fnorm*iph(jt)
     CALL pphhmtx(ja,jb,jc,jd,jt,ans)
     cross=cross+ans*ang_mom_factor
  ENDDO

END SUBROUTINE cross_coupled_mtxel1
!
!     Calculates the crosscoupled matrix element type 2
!
SUBROUTINE cross_coupled_mtxel2(ja,jb,jc,jd,jtot,cross)
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE 
  REAL(DP), DIMENSION (n_startenergy_g) :: ans
  INTEGER :: ja,jb,jc,jd,jtot, jbra_min, jket_min, jbra_max ,jket_max, &
       jt, j_min, j_max, iph
  REAL(DP) :: fnorm, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_g), INTENT(OUT) :: cross

  cross=0.
  fnorm=SQRT(2.*jtot+1.)*iph((all_orbit%jj(ja)+all_orbit%jj(jd))/2+ &
       jtot+all_orbit%jj(jb))
  IF(iph(all_orbit%ll(ja)+all_orbit%ll(jb)) & 
       /=iph(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  jbra_min=ABS(all_orbit%jj(ja)-all_orbit%jj(jb))
  jbra_max=all_orbit%jj(ja)+all_orbit%jj(jb)
  jket_min=ABS(all_orbit%jj(jc)-all_orbit%jj(jd))
  jket_max=all_orbit%jj(jc)+all_orbit%jj(jd)
  j_max=MIN(jbra_max,jket_max)/2
  j_min=MAX(jbra_min,jket_min)/2
  IF(j_min > j_max) RETURN
  DO jt=j_min,j_max
     ang_mom_factor=sjs(all_orbit%jj(jc),all_orbit%jj(jb),2*jtot, &
          all_orbit%jj(ja),all_orbit%jj(jd),2*jt)* &
          (2.*jt+1.)*fnorm
     CALL pphhmtx(ja,jb,jc,jd,jt,ans)
     cross=cross+ans*ang_mom_factor
  ENDDO

END SUBROUTINE cross_coupled_mtxel2
!
!     This function returns the matrix for V
!     Plain brute force
!
SUBROUTINE decompose(ja,jb,jc,jd,jt,Vdeco)
  USE stored_bare_interaction
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  LOGICAL TRIAG
  INTEGER, INTENT(IN) :: ja,jb,jc,jd,jt
  INTEGER :: iph, ia, ib, ic, id, JPmin, JPmax
  INTEGER ::  jja, jjb, jjc, jjd, la, lb, lc, ld, k, LP, LPmin, LPmax, L, LLmin, LLmax, S, SP, JP
  REAL(DP), INTENT(INOUT), DIMENSION(n_startenergy_g,0:2) :: Vdeco
  REAL(DP), DIMENSION(n_startenergy_g) :: gmtx
  REAL(DP) :: dij, delta, w6, w6p, w9, w9p, w9pp, w9ppp

  Vdeco = 0.
  IF(2*all_orbit%nn(ja)+all_orbit%ll(ja) +2*all_orbit%nn(jb)+all_orbit%ll(jb) > nlmax ) RETURN
  IF(2*all_orbit%nn(jc)+all_orbit%ll(jc) +2*all_orbit%nn(jd)+all_orbit%ll(jd) > nlmax ) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  IF((-1)**(all_orbit%ll(ja)+all_orbit%ll(jb)) /=  &
       (-1)**(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((ja == jb).AND.(MOD(jt,2)/=0)) RETURN       
  IF((jc == jd).AND.(MOD(jt,2)/=0)) RETURN       
  IF(triag(all_orbit%jj(ja),all_orbit%jj(jb),2*jt)) RETURN
  IF(triag(all_orbit%jj(jc),all_orbit%jj(jd),2*jt)) RETURN

  ! The single-particle j-values multiplied by two
  jja = all_orbit%jj(ja)
  jjb = all_orbit%jj(jb)
  jjc = all_orbit%jj(jc)
  jjd = all_orbit%jj(jd)
  ! The single-particle l-values
  la = all_orbit%ll(ja)
  lb = all_orbit%ll(jb)
  lc = all_orbit%ll(jc)
  ld = all_orbit%ll(jd)
  ! LS-jj variables for bra and ket side
  LPmin =  ABS(la-lb); LPmax = la+lb; LLmin = ABS(lc-ld); LLmax = lc+ld
  Vdeco = 0.
  ! Start decomposition
  DO k = 0, 2
     DO LP = LPmin, LPmax
        DO SP = 0, 1
           IF (triag(LP,SP, jt)) CYCLE
           !  First LS-jj transformation coefficient
           w9p=snj(2*la,1,jja,2*lb,1,jjb,2*LP,2*SP,2*JT)
           w9p=w9p*SQRT((jja+1.)*(jjb+1.)*(2*LP+1.)*(2*SP+1.))
           DO L = LLmin, LLmax
              IF ( triag (L, LP, k)) CYCLE
              DO S = 0, 1
                 IF ( SP /= S ) CYCLE
                 IF ( triag (L, S, jt)) CYCLE
                 IF ( triag (S, SP, k)) CYCLE
                 !   Second LS-jj transformation coefficient
                 w9=snj(2*lc,1,jjc,2*ld,1,jjd,2*L,2*S,2*JT)
                 w9=w9*SQRT((jjc+1.)*(jjd+1.)*(2*L+1.)*(2*S+1.))*w9p
                 ! First 6j symbol
                 w6=sjs(2*LP,2*SP, 2*jt,2*S,2*L,2*k)*w9
                 JPmin = ABS(LP-SP); JPmax = LP+SP
                 DO JP = JPmin, JPmax
                    ! Second 6j symbol
                    w6p=sjs(2*LP,2*SP, 2*JP,2*S,2*L,2*k)*(2*JP+1)*((-1.0)**JP)*w6
                    ! Loops over the new single-particle states
                    DO ia = 1, all_orbit%total_orbits
                       IF(all_orbit%ll(ia) /= la) CYCLE
                       IF(all_orbit%itzp(ia) /= all_orbit%itzp(ja) ) CYCLE
                       IF(all_orbit%nn(ia) /= all_orbit%nn(ja) ) CYCLE
                       IF ( triag(2*la, 1,all_orbit%jj(ia))) CYCLE      
                       DO ib = 1, all_orbit%total_orbits
                          IF(all_orbit%ll(ib) /= lb) CYCLE
                          IF(all_orbit%itzp(ib) /= all_orbit%itzp(jb) ) CYCLE
                          IF ( triag(2*lb, 1,all_orbit%jj(ib)) )CYCLE    
                          IF(all_orbit%nn(ib) /= all_orbit%nn(jb) ) CYCLE  
                          w9pp=snj(2*la,1,all_orbit%jj(ia),2*lb,1,all_orbit%jj(ib),2*LP,2*SP,2*JP)
                          w9pp=w9pp*w6p*SQRT((all_orbit%jj(ia)+1.)*(all_orbit%jj(ib)+1.)*(2*LP+1.)*(2*SP+1.))
                          DO ic = 1, all_orbit%total_orbits
                             IF(all_orbit%ll(ic) /= lc) CYCLE
                             IF(all_orbit%itzp(ic) /= all_orbit%itzp(jc) ) CYCLE
                             IF ( triag(2*lc, 1,all_orbit%jj(ic))) CYCLE      
                             IF(all_orbit%nn(ic) /= all_orbit%nn(jc) ) CYCLE
                             DO id = 1, all_orbit%total_orbits
                                IF(all_orbit%ll(id) /= ld) CYCLE
                                IF(all_orbit%itzp(id) /= all_orbit%itzp(jd) ) CYCLE
                                IF ( triag(2*ld, 1,all_orbit%jj(id))) CYCLE      
                                IF(all_orbit%nn(id) /= all_orbit%nn(jd) ) CYCLE
                                w9ppp=snj(2*lc,1,all_orbit%jj(ic),2*ld,1,all_orbit%jj(id),2*L,2*S,2*JP)
                                w9ppp=w9ppp*SQRT((all_orbit%jj(ic)+1.)*(all_orbit%jj(id)+1.)*(2*L+1.)*(2*S+1.))
                                CALL pphhmtx(ia,ib,ic,id,jp,gmtx)
                                Vdeco(:,k) = Vdeco(:,k) + gmtx(:)*w9ppp*w9pp
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     Vdeco(:,k) =  Vdeco(:,k)*((-1.0)**jt)*(2*k+1.0)/dij(ja,jb)/dij(jc,jd)
  ENDDO

END SUBROUTINE decompose




