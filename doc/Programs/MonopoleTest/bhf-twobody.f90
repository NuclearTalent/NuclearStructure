!             Program block bhf-twobody.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO, Norway 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90/F95 
!             LAST UPGRADE : Oct 2007
!                            
!             Program to calculate effective interactions 
!             for nuclei within the spherical shell
!             model. The program runs in the proton-neutron
!             formalism. 
!             Syntax throughout is that of free source format,
!             additional compiler options may be necessary
!             The code computes effective interactions for
!             the nuclear shell model and allows one to 
!             perform Hartree-Fock calculations as well.
!
!
!      Sets up onebody contributions
!
SUBROUTINE onebody_contribution
  USE single_particle_orbits
  USE onebody_diagrams
  USE constants
  USE wave_functions
  IMPLICIT NONE
  REAL(DP), DIMENSION(n_startenergy_veff) :: onebody_diagram_1, &
       onebody_diagram_2,onebody_diagram_3, onebody_diagram_5, &
       onebody_diagram_6,onebody_diagram_7,onebody_diagram_8, &
       onebody_diagram_9,onebody_diagram_4, onebody_diagram_10, &
       onebody_diagram_11,onebody_diagram_12,onebody_diagram_13, &
       onebody_diagram_14,onebody_diagram_15,onebody_diagram_16, &
       onebody_diagram_17,onebody_diagram_18,onebody_diagram_19, &
       onebody_diagram_20,onebody_diagram_21, onebody_diagram_monopole
  REAL(DP), DIMENSION(1,1,n_startenergy_veff) :: sbox
  REAL(DP) :: sum_kin, particle_mass
  REAL(DP), DIMENSION(all_orbit%total_orbits,all_orbit%total_orbits) :: kinenergy
  INTEGER :: a, c, iph, k

  ALLOCATE (  one_body_terms &
       (all_orbit%total_orbits, &
       all_orbit%total_orbits,n_startenergy_veff) )
  ALLOCATE (  one_body_folded &
       (all_orbit%total_orbits, &
       all_orbit%total_orbits) )
  kinenergy = 0.0_dp
  WRITE(6,*) 'a, c,  onebody part, kinetic energy and  total sp energy'
  DO a = 1, all_orbit%total_orbits
     IF (all_orbit%model_space(a) == 'outside') CYCLE         
     DO c= a, all_orbit%total_orbits
        !        CALL one_body_diagram_phase(a,c,phase)
        sbox=0.0_dp
        IF (all_orbit%model_space(c) == 'outside') CYCLE         
        IF(iph(all_orbit%ll(a)) /= iph(all_orbit%ll(c))) CYCLE
        IF ( all_orbit%jj(a) /= all_orbit%jj(c)) CYCLE
        IF(all_orbit%itzp(a) /= all_orbit%itzp(c) ) CYCLE
        IF((itzmin == 1).AND.(all_orbit%itzp(a) == -1)) CYCLE
        IF((itzmax == -1).AND.(all_orbit%itzp(a) == 1) ) CYCLE
        SELECT CASE (physical_system) 
        CASE('nuclear_physics')
           particle_mass = p_mass(all_orbit%itzp(a))
           sum_kin=0.0_dp
           IF ( hf_iterations == 0 ) THEN
              DO k=1,n_rel
                 sum_kin=sum_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))*hol(k,all_orbit%ll(c),all_orbit%nn(c))  &
                      *wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ELSE
              DO k=1,n_rel
                 sum_kin=sum_kin+wave_function(k,a)*wave_function(k,c)*wra(k)*(ra(k)**4)*0.5*hbarc*hbarc/particle_mass
              ENDDO
           ENDIF
        CASE('atomic_physics')
           sum_kin=0.0_dp
           IF ( hf_iterations == 0 ) THEN
              DO k=1,n_rel
                 sum_kin=sum_kin+hol(k,all_orbit%ll(a),all_orbit%nn(a))*hol(k,all_orbit%ll(c),all_orbit%nn(c))  &
                      *wra(k)*(ra(k)**4)*0.5
              ENDDO
           ELSE
              DO k=1,n_rel
                 sum_kin=sum_kin+wave_function(k,a)*wave_function(k,c)*wra(k)*(ra(k)**4)*0.5
              ENDDO
           ENDIF
        END SELECT
        SELECT CASE (type_of_renormv)
        CASE ('no-core')
           IF ( a == c) THEN
              sum_kin=2*sum_kin
           ELSE
              sum_kin = 0.0_dp
           ENDIF
        CASE ('v-nrg')
           IF ( a == c) THEN
              sum_kin=2*sum_kin
           ELSE
              sum_kin = 0.0_dp
           ENDIF
        CASE ('vlowk')
           sum_kin = a_factor*sum_kin
        CASE ('v-krg')
           sum_kin = a_factor*sum_kin
        CASE ('g-matrix')
           sum_kin = sum_kin*a_factor
        END SELECT
        kinenergy(c,a) = sum_kin
        onebody_diagram_1 = 0.0_dp
        SELECT CASE ( order_of_interaction )
        CASE ('first')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_monopole(a,c,onebody_diagram_monopole)
           one_body_terms(a,c,:)=onebody_diagram_1(:)
        CASE ('second')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_2(a,c,onebody_diagram_2)
           CALL diagram_3(a,c,onebody_diagram_3)
           one_body_terms(a,c,:)= &
                (onebody_diagram_1(:)+ &
                onebody_diagram_2(:)+ &
                onebody_diagram_3(:))
        CASE ('third')
           CALL diagram_1(a,c,onebody_diagram_1)
           CALL diagram_2(a,c,onebody_diagram_2)
           CALL diagram_3(a,c,onebody_diagram_3)
           CALL diagram_4(a,c,onebody_diagram_4)
           CALL diagram_5(a,c,onebody_diagram_5)
           CALL diagram_6(a,c,onebody_diagram_6)
           CALL diagram_7(a,c,onebody_diagram_7)
           CALL diagram_8(a,c,onebody_diagram_8)
           CALL diagram_9(a,c,onebody_diagram_9)
           CALL diagram_10(a,c,onebody_diagram_10)
           CALL diagram_11(a,c,onebody_diagram_11)
           CALL diagram_12(a,c,onebody_diagram_12)
           CALL diagram_13(a,c,onebody_diagram_13)
           CALL diagram_14(a,c,onebody_diagram_14)
           CALL diagram_15(a,c,onebody_diagram_15)
           CALL diagram_16(a,c,onebody_diagram_16)
           CALL diagram_17(a,c,onebody_diagram_17)
           CALL diagram_18(a,c,onebody_diagram_18)
           CALL diagram_19(a,c,onebody_diagram_19)
           CALL diagram_20(a,c,onebody_diagram_20)        
           CALL diagram_21(a,c,onebody_diagram_21)
           one_body_terms(a,c,:)= &
                (onebody_diagram_1(:)+ &
                onebody_diagram_2(:)+ &
                onebody_diagram_3(:)+onebody_diagram_4(:)+onebody_diagram_5(:)+ &
                onebody_diagram_6(:)+onebody_diagram_7(:)+onebody_diagram_8(:)+ &
                onebody_diagram_9(:)+onebody_diagram_10(:)+onebody_diagram_11(:)+ &
                onebody_diagram_12(:)+onebody_diagram_13(:)+onebody_diagram_14(:)+ &
                onebody_diagram_15(:)+onebody_diagram_16(:)+onebody_diagram_17(:)+ &
                onebody_diagram_18(:)+onebody_diagram_19(:)+onebody_diagram_20(:)+ &
                onebody_diagram_21(:))
        END SELECT
        one_body_terms(c,a,:)=one_body_terms(a,c,:)
        sbox(1,1,:)=one_body_terms(a,c,:)

        WRITE(6,'(2I4,2X,4F12.6)') a , c, one_body_terms(a,c,n_startenergy_veff/2+1), &
             onebody_diagram_monopole(n_startenergy_veff/2+1), &
             kinenergy(c,a),one_body_terms(a,c,n_startenergy_veff/2+1)+kinenergy(c,a)
     ENDDO
  ENDDO

END SUBROUTINE onebody_contribution
