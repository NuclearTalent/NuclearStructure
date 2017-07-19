!             Program block heff-diagrams.f
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F90  
!             LAST UPGRADE : October 2007
!
!             This program block sets up all types of diagrams for
!             one-body, two-body, three-body, effective operators and closed core diagrams.
!             Three-body diagrams are only to second order in the interaction
!
!     Begin one-body diagrams
!
!
!
!
!     HF diagram
!
SUBROUTINE diagram_monopole(a,c,onebody_diagram_1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, i, k
  REAL(DP), DIMENSION(n_startenergy_veff) :: w 
  REAL(DP) :: val, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g,0:2) :: ans
  REAL(DP), DIMENSION(n_startenergy_veff) :: sum
  REAL(DP), DIMENSION(0:2) :: kdiagram

  onebody_diagram_1=0.; kdiagram = 0.0
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
     ang_mom_factor = (all_orbit%jj(h)+1.)
     CALL monopole(a,h,ans)
     DO k = 0, 2
        sum = 0.0
        DO i=1, n_startenergy_veff
           call interpolate(w(i),e_start_g,ans(:,k),val)
           sum(i) = sum(i) + val*ang_mom_factor
           onebody_diagram_1(i)=onebody_diagram_1(i)+val*ang_mom_factor
        ENDDO
        kdiagram(k) = kdiagram(k)+sum(n_startenergy_veff/2+1)
     ENDDO
  ENDDO
  WRITE(6,'(2X,2I3,2X,3F12.8)') a, a, kdiagram
END  SUBROUTINE diagram_monopole

!
!     HF diagram
!
SUBROUTINE monopole(a,i,answer)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, i
  INTEGER :: j_min, j_max, jph, k
  REAL(DP) :: ang_mom_factor, sumJ
  REAL(DP), DIMENSION(n_startenergy_g,0:2), INTENT(OUT) :: answer
  REAL(DP), DIMENSION(n_startenergy_g,0:2) :: ans
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1

  answer = 0.;   sumJ = 0.
  j_min=ABS((all_orbit%jj(a)-all_orbit%jj(i))/2)
  j_max=(all_orbit%jj(a)+all_orbit%jj(i))/2
  DO jph=j_min,j_max
     ang_mom_factor=(2.*jph+1.)
     sumJ= sumJ + ang_mom_factor
     CALL decompose(a,i,a,i,jph,ans)
!     DO k = 0, 2
!        answer(:,k) = answer(:,k)+ans(:,k)*ang_mom_factor
        answer = answer+ ans*ang_mom_factor
!     ENDDO
  ENDDO
  answer = answer/sumJ
END  SUBROUTINE monopole
!
!     HF diagram
!
SUBROUTINE diagram_1(a,c,onebody_diagram_1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: j_min, j_max, jph, h, i
  REAL(DP), DIMENSION(n_startenergy_veff) :: w 
  REAL(DP) :: val, ang_mom_factor
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  onebody_diagram_1=0.
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
     DO jph=j_min,j_max
        ang_mom_factor=(2.*jph+1.)/(all_orbit%jj(a)+1.)
        CALL pphhmtx(a,h,c,h,jph,ans) ; IF ( ans(1) == 0.0_dp) CYCLE
        DO i=1, n_startenergy_veff
           call interpolate(w(i),e_start_g,ans,val)
           onebody_diagram_1(i)=onebody_diagram_1(i)+val*ang_mom_factor
        ENDDO
     ENDDO
  ENDDO
END  SUBROUTINE diagram_1



!
!     2p1h diagram
!
SUBROUTINE diagram_2(a,c,onebody_diagram_2)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h, p1, p2, i, nshell1, nshell2, iosc
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL DENCHECK

  onebody_diagram_2=0.
  DO h=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE         
     nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) == 'hole') CYCLE         
           DO p2=1, all_orbit%total_orbits
              IF (all_orbit%orbit_status(p2) == 'hole') CYCLE         
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              iosc=nshell2-nshell1
              IF(dencheck(iosc))CYCLE
              de=all_orbit%evalence(c)+all_orbit%e(h)- &
                   all_orbit%e(p1)-all_orbit%e(p2)+wcn
              CALL pphhmtx(a,h,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              CALL pphhmtx(p1,p2,c,h,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
              DO  i=1, n_startenergy_veff
                 CALL interpolate(w(i),e_start_g,ans1,val1)
                 CALL interpolate(w(i),e_start_g,ans2,val2)
                 onebody_diagram_2(i)=onebody_diagram_2(i)+ &
                      factr*0.5*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_2
!
!     2h1p diagram
!
SUBROUTINE diagram_3(a,c,onebody_diagram_3)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, p, i, h1, h2, nshell1, nshell2, iosc
  REAL(DP) :: val1, val2, factr
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2 , de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_3
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL DENCHECK

  onebody_diagram_3=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) == 'hole') CYCLE         
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/(all_orbit%jj(a)+1.)
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE         
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE         
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
              iosc=nshell2-nshell1
              IF(dencheck(iosc)) CYCLE
              de=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p)-all_orbit%evalence(a)+wcn
              CALL pphhmtx(h1,h2,c,p,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              CALL pphhmtx(a,p,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO i=1, n_startenergy_veff
                 CALL  interpolate(w1(i),e_start_g,ans1,val1)
                 CALL  interpolate(w2(i),e_start_g,ans2,val2)
                 onebody_diagram_3(i)=onebody_diagram_3(i)- & 
                      factr*0.5d0*val1*val2/de(i)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_3



SUBROUTINE diagram_4(a,c,onebody_diagram_4)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h, p1, p2, p3, p4, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w, de 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_4
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_4=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              IF (dencheck(iosc1)) CYCLE
              CALL pphhmtx(a,h,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              den1=all_orbit%evalence(c)+all_orbit%e(h)-all_orbit%e(p1)-all_orbit%e(p2)
              DO p3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                 DO p4=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p4) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(p3)+all_orbit%nshell(p4)
                    iosc2=nshell1-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%evalence(c)+all_orbit%e(h)-  &
                         all_orbit%e(p3)-all_orbit%e(p4)
                    CALL pphhmtx(p1,p2,p3,p4,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p3,p4,c,h,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w(i),e_start_g,ans1,val1)
                       CALL interpolate(w(i),e_start_g,ans2,val2)
                       CALL interpolate(w(i),e_start_g,ans3,val3)
                       onebody_diagram_4(i)=onebody_diagram_4(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_4



SUBROUTINE diagram_5(a,c,onebody_diagram_5)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, p3, h2, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_5
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_5=0.
  DO p3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p3))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              nshell2=all_orbit%nshell(a)+all_orbit%nshell(p3)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p3)
              CALL pphhmtx(h1,h2,c,p3,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 DO p2=1, all_orbit%total_orbits                    
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(p1)+& 
                         all_orbit%nshell(p2)
                    iosc2=nshell1-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    CALL pphhmtx(a,p3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=all_orbit%e(h1)+all_orbit%e(h2)+ &
                         all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_5(i)=onebody_diagram_5(i)- &
                            0.25*val1*factr*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_5



SUBROUTINE diagram_6(a,c,onebody_diagram_6)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, p3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de      
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_6
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_6=0.
  DO p2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p2))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
     nshell2=all_orbit%nshell(a)+all_orbit%nshell(p2)
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              DO p1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 CALL pphhmtx(a,p2,h1,h2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                 DO p3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(c)+nshell1
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p3)+ &
                         all_orbit%nshell(a)
                    iosc2=nshell4-nshell3
                    den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p3)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(h1,h2,p1,p3,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p3,c,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                    w2=all_orbit%e(h2)+all_orbit%e(h1)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w3=w2
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_6(i)=onebody_diagram_6(i)- &
                            0.25d0*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_6


SUBROUTINE diagram_7(a,c,onebody_diagram_7)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p, h2, h3, h4, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_7
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_7=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE   
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO h1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE   
           DO h2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p)- &
                   all_orbit%evalence(a)
              CALL pphhmtx(a,p,h1,h2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO h3=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE 
                 DO h4=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h4)
                    iosc2=nshell2-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%e(p)- &
                         all_orbit%evalence(a)
                    CALL pphhmtx(h1,h2,h3,h4,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(h3,h4,c,p,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w2=all_orbit%e(h3)+all_orbit%e(h4)-all_orbit%evalence(a)- &
                         all_orbit%e(p)+w1
                    w3=all_orbit%e(h3)+all_orbit%e(h4)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_7(i)=onebody_diagram_7(i)- &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_7



SUBROUTINE diagram_8(a,c,onebody_diagram_8)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, iosc1, iosc2, i
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_8
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_8=0.
  DO h3=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h3))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h3))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h3)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h3)+all_orbit%evalence(c)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    iosc2=nshell3-nshell2
                    den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(a,h3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(h1,h2,c,h3,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%evalence(c)+all_orbit%e(h3)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=w2+all_orbit%e(h3)+all_orbit%evalence(c)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_8(i)=onebody_diagram_8(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_8



SUBROUTINE diagram_9(a,c,onebody_diagram_9)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j_min, j_max, jtot, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP) :: val1, val2, val3, factr,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_9
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_9=0.
  DO h2=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h2))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
     DO jtot=j_min,j_max
        factr=(2.*jtot+1.)/((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              iosc1=nshell2-nshell1
              den1=all_orbit%e(h2)+all_orbit%evalence(c)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              IF(dencheck(iosc1)) CYCLE
              DO h1=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                 DO h3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h1)+all_orbit%nshell(c)+ &
                         all_orbit%nshell(h3)
                    nshell4=all_orbit%nshell(a)+nshell2
                    iosc2=nshell4-nshell3
                    den2=all_orbit%e(h1)+all_orbit%e(h3)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    IF(dencheck(iosc2)) CYCLE
                    CALL pphhmtx(h1,h3,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(a,h2,h1,h3,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL pphhmtx(p1,p2,c,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    w1=all_orbit%e(h1)+all_orbit%e(h3)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h2)+all_orbit%evalence(a)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+w1
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_9(i)=onebody_diagram_9(i)+ &
                            0.25*factr*val1*val2*val3/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_9



SUBROUTINE diagram_10(a,c,onebody_diagram_10)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_10
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_10=0.
  DO h1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        IF(all_orbit%jj(h1).ne.all_orbit%jj(h2)) CYCLE
        IF(iph(all_orbit%ll(c)+all_orbit%ll(h2)) /= &
             iph(all_orbit%ll(a)+all_orbit%ll(h1))) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(h1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
        DO h3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
           par1=iph(all_orbit%ll(h1)+all_orbit%ll(h3))
           par2=iph(all_orbit%ll(h2)+all_orbit%ll(h3))
           IF(par1 /= par2) CYCLE 
           j2min=ABS((all_orbit%jj(h3)-all_orbit%jj(h1))/2)
           j2max=(all_orbit%jj(h3)+all_orbit%jj(h1))/2
           DO jtot1=j1min,j1max
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(h1))+1.)
                 DO p1=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
                    DO p2=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                       nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h3)
                       nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                       iosc1=nshell2-nshell1
                       nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)+ &
                            all_orbit%nshell(c)
                       nshell4=nshell2+all_orbit%nshell(a)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc1).or.dencheck(iosc2)) CYCLE
                       den1=all_orbit%e(h1)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)
                       den2=all_orbit%e(h2)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)
                       CALL pphhmtx(a,h1,c,h2,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h2,h3,p1,p2,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p1,p2,h1,h3,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)+all_orbit%e(h3)- &
                            all_orbit%e(p1)-all_orbit%e(p2)+wcn
                       w2=all_orbit%e(h2)+all_orbit%e(h3)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h3)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_10(i)=onebody_diagram_10(i)- &
                               0.5*factr1*factr2*  &
                               val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_10

SUBROUTINE diagram_11(a,c,onebody_diagram_11)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_11
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_11=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE    
     DO p2=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
        IF(all_orbit%jj(p1) /= all_orbit%jj(p2)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        DO p3=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(p3))
           par2=iph(all_orbit%ll(p2)+all_orbit%ll(p3))
           IF(par1 /= par2) CYCLE 
           j2min=ABS((all_orbit%jj(p3)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p3)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO h1=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                    DO h2=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                       nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p3)
                       iosc1=nshell2-nshell1
                       nshell3=all_orbit%nshell(c)+nshell1
                       nshell4=all_orbit%nshell(p3)+all_orbit%nshell(p2)+ &
                            all_orbit%nshell(a)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc1).or.dencheck(iosc2)) CYCLE
                       den1=all_orbit%e(h1)+all_orbit%e(h2)-  &
                            all_orbit%e(p1)-all_orbit%e(p3)
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)
                       CALL pphhmtx(a,p2,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p1,p3,h1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%e(p3)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_11(i)=onebody_diagram_11(i) &
                               +0.5*factr1*factr2*val1* &
                               val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_11


SUBROUTINE diagram_12(a,c,onebody_diagram_12)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_12
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_12=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)
        nshell2=all_orbit%nshell(p1)
        iosc1=nshell2-nshell1
        IF(dencheck(iosc1)) CYCLE
        den1=all_orbit%e(h1)-all_orbit%e(p1)

        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE    
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(h2))
           par2=iph(all_orbit%ll(h1)+all_orbit%ll(h2))
           IF(par1 /= par2) CYCLE
           j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(h2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,h1,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO p2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                    DO p3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                       nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell4=all_orbit%nshell(p3)+all_orbit%nshell(p2)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)
                       CALL pphhmtx(p1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p2,p3,h1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_12(i)=onebody_diagram_12(i)+0.5 &
                               *factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_12


SUBROUTINE diagram_13(a,c,onebody_diagram_13)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, p3, h2, par1, par2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, iph
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_13
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_13=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
        nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
        iosc1=nshell2-nshell1
        den1=all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%evalence(a)-all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        DO h2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE    
           par1=iph(all_orbit%ll(p1)+all_orbit%ll(h2))
           par2=iph(all_orbit%ll(h1)+all_orbit%ll(h2))
           IF(par1 /= par2) CYCLE
           j2min=ABS((all_orbit%jj(h2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(h2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,p1,c,h1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO p2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
                    DO p3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE 
                       nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                       nshell4=all_orbit%nshell(p2)+all_orbit%nshell(p3)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                            all_orbit%e(p2)-all_orbit%e(p3)
                       CALL pphhmtx(h1,h2,p2,p3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(p2,p3,p1,h2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%evalence(c)-all_orbit%evalence(a)+wcn
                       w3=w2
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_13(i)=onebody_diagram_13(i)+ &
                               0.5*factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_13




SUBROUTINE diagram_14(a,c,onebody_diagram_14)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2
  REAL(DP) :: val1, val2, val3, factr2, factr1,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_14
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_14=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE    
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)
        nshell2=all_orbit%nshell(p1)
        iosc1=nshell2-nshell1
        den1=all_orbit%e(h1)-all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        DO p2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
           j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
           DO jtot1=j1min,j1max
              CALL pphhmtx(a,h1,c,p1,jtot1,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
              DO jtot2=j2min,j2max
                 factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE 
                    DO h3=1, all_orbit%total_orbits
                       IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE    
                       nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h3)
                       nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                       iosc2=nshell4-nshell3
                       IF(dencheck(iosc2)) CYCLE
                       den2=all_orbit%e(h3)+all_orbit%e(h2)- &
                            all_orbit%e(p1)-all_orbit%e(p2)
                       CALL pphhmtx(p1,p2,h2,h3,jtot2,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                       CALL pphhmtx(h2,h3,h1,p2,jtot2,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                       w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
                       w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                       w3=all_orbit%e(h1)+all_orbit%e(h2)+ &
                            all_orbit%e(h3)-all_orbit%e(p1)+wcn
                       de=(den1+wcn)*(den2+wcn)
                       DO i=1, n_startenergy_veff
                          CALL interpolate(w1(i),e_start_g,ans1,val1)
                          CALL interpolate(w2(i),e_start_g,ans2,val2)
                          CALL interpolate(w3(i),e_start_g,ans3,val3)
                          onebody_diagram_14(i)=onebody_diagram_14(i)- &
                               0.5*factr1*factr2*val1*val2*val3/de(i)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_14
!
!
!
SUBROUTINE diagram_15(a,c,onebody_diagram_15)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  j1min, j1max, h1, p1, p2, h3, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j2min, j2max, jtot1, jtot2, ie1, ie2
  REAL(DP) :: val1, val2, val3, factr2, factr1
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1,w2,w3,den1, den2, &
       de,xmtx1,xmtx2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_15
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_15=0.
  DO p1=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE 
     DO h1=1, all_orbit%total_orbits
        IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE 
        IF(all_orbit%jj(p1) /= all_orbit%jj(h1)) CYCLE
        j1min=ABS((all_orbit%jj(a)-all_orbit%jj(p1))/2)
        j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
        nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
        nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
        iosc1=nshell2-nshell1
        den1=wcn+all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%evalence(a)- &
             all_orbit%e(p1)
        IF(dencheck(iosc1)) CYCLE
        xmtx1=0.
        DO jtot1=j1min,j1max
           factr1=(2.*jtot1+1.)/((all_orbit%jj(a))+1.)
           CALL pphhmtx(a,p1,c,h1,jtot1,ans1); IF(ans1(1).eq.0.) CYCLE 
           w1=all_orbit%evalence(c)+all_orbit%e(h1)+wcn
           DO ie1=1, n_startenergy_veff
              CALL interpolate(w1(ie1),e_start_g,ans1,val1)
              xmtx1(ie1)=xmtx1(ie1)+factr1*val1
           ENDDO
        ENDDO
        xmtx2=0.
        DO p2=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE 
           j2min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
           j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
           DO jtot2=j2min,j2max
              factr2=(2.*jtot2+1.)/((all_orbit%jj(p1))+1.)
              DO h2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE       
                 DO h3=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE 
                    nshell3=all_orbit%nshell(h3)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=wcn+all_orbit%e(h3)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    CALL pphhmtx(h1,p2,h2,h3,jtot2,ans2); IF(ans2(1).eq.0.) CYCLE 
                    CALL pphhmtx(h2,h3,p1,p2,jtot2,ans3); IF(ans3(1).eq.0.) CYCLE 
                    w2=all_orbit%evalence(c)+all_orbit%e(h1)+all_orbit%e(h2)+ &
                         all_orbit%e(h3)-all_orbit%evalence(a)-all_orbit%e(p1)+wcn
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%e(h3)- &
                         all_orbit%evalence(a)+wcn
                    de=den1*den2
                    DO ie2=1, n_startenergy_veff
                       CALL interpolate(w2(ie2),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie2),e_start_g,ans3,val3)
                       xmtx2(ie2)=xmtx2(ie2)+factr2*val2*val3/de(ie2)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        onebody_diagram_15=onebody_diagram_15-0.5*xmtx1*xmtx2
     ENDDO
  ENDDO

END SUBROUTINE diagram_15
!
!
!
SUBROUTINE diagram_16(a,c,onebody_diagram_16)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_16
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_16=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(a,h1,h,p1,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                         all_orbit%e(p1)-all_orbit%e(p2)
                    den2=all_orbit%e(h2)+all_orbit%e(h)- &
                         all_orbit%evalence(a)-all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h,h2,c,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=w2+all_orbit%e(h)-all_orbit%e(p2)
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_16(i)=onebody_diagram_16(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_16


SUBROUTINE diagram_17(a,c,onebody_diagram_17)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_17
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_17=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(h,p1,c,h1,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)+ &
                         all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)+ &
                         all_orbit%nshell(a)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(h)
                    nshell4=all_orbit%nshell(p2)+all_orbit%nshell(a)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)+all_orbit%evalence(c)-all_orbit%evalence(a)
                    den2=all_orbit%e(h2)+all_orbit%e(h)-all_orbit%evalence(a)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(a,p2,h,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=w1+all_orbit%e(h)-all_orbit%e(p2)
                    w3=all_orbit%e(h)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_17(i)=onebody_diagram_17(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_17


SUBROUTINE diagram_18(a,c,onebody_diagram_18)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_18
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_18=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(p,h1,c,p1,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(a,h2,p,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+  &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(p)+all_orbit%jj(h1)+ &
                         all_orbit%jj(p2)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%e(p2)+wcn
                    w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_18(i)=onebody_diagram_18(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_18



SUBROUTINE diagram_19(a,c,onebody_diagram_19)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_19
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_19=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              CALL cross_coupled_mtxel1(a,p1,p,h1,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                    nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                    iosc1=nshell2-nshell1
                    nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2).or.dencheck(iosc1))CYCLE
                    den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                         all_orbit%e(p2)
                    den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(p,p2,c,h2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h1)+all_orbit%jj(p2)+ &
                         all_orbit%jj(p)-jtot)
                    w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)- &
                         all_orbit%evalence(a)+wcn
                    w2=w1+all_orbit%evalence(a)-all_orbit%e(p2)
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_19(i)=onebody_diagram_19(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_19



SUBROUTINE diagram_20(a,c,onebody_diagram_20)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, p, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg,  den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_20
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_20=0.
  DO p=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) == 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(p))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(p))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(c)+all_orbit%nshell(h1)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%evalence(c)+all_orbit%e(h1)-all_orbit%e(p)- &
                   all_orbit%e(p1)
              CALL cross_coupled_mtxel1(p,p1,c,h1,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h2)+all_orbit%nshell(c)
                    nshell4=all_orbit%nshell(p)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(a,h2,p,p2,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(h1,p2,p1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(p)+all_orbit%jj(p2))/2+  &
                         all_orbit%jj(h1)+all_orbit%jj(p2)+ &
                         all_orbit%jj(p)-jtot)
                    w1=all_orbit%e(h2)+all_orbit%evalence(c)+wcn
                    w2=w1+all_orbit%e(h1)-all_orbit%e(p)
                    w3=all_orbit%evalence(c)+all_orbit%e(h1)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_20(i)=onebody_diagram_20(i)+ &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_20



SUBROUTINE diagram_21(a,c,onebody_diagram_21)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER ::  h1, p1, p2, h, h2, &
       nshell1, nshell2, nshell3, nshell4, iosc1, iosc2, i, &
       j_min, j_max, jtot, iph
  REAL(DP) :: val1, val2, val3, factr, sg, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff) :: w1, w2, w3, de
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: onebody_diagram_21
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL DENCHECK

  onebody_diagram_21=0.
  DO h=1, all_orbit%total_orbits
     IF(all_orbit%orbit_status(h) /= 'hole' ) CYCLE
     j_min=ABS((all_orbit%jj(a)-all_orbit%jj(h))/2)
     j_max=(all_orbit%jj(a)+all_orbit%jj(h))/2
     DO jtot=j_min,j_max
        factr=sqrt(2.*jtot+1.)*((all_orbit%jj(a))+1.)
        DO p1=1, all_orbit%total_orbits
           IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
           DO h1=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(h1)
              nshell2=all_orbit%nshell(p1)+all_orbit%nshell(a)
              iosc1=nshell2-nshell1
              IF(dencheck(iosc1)) CYCLE
              den1=all_orbit%e(h)+all_orbit%e(h1)-all_orbit%evalence(a)- &
                   all_orbit%e(p1)
              CALL cross_coupled_mtxel1(a,p1,h,h1,jtot,ans1) ; IF ( ans1(1) == 0.0_dp) CYCLE
              DO p2=1, all_orbit%total_orbits
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 DO h2=1, all_orbit%total_orbits
                    IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                    nshell3=all_orbit%nshell(h)+all_orbit%nshell(h2)
                    nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
                    iosc2=nshell4-nshell3
                    IF(dencheck(iosc2)) CYCLE
                    den2=all_orbit%e(h)+all_orbit%e(h2)-all_orbit%evalence(a)- &
                         all_orbit%e(p2)
                    CALL cross_coupled_mtxel1(h1,p2,p1,h2,jtot,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
                    CALL cross_coupled_mtxel1(h,h2,c,p2,jtot,ans3) ; IF ( ans3(1) == 0.0_dp) CYCLE
                    sg=iph((all_orbit%jj(p1)+all_orbit%jj(h1)+ &
                         all_orbit%jj(h2)+all_orbit%jj(a)+ &
                         all_orbit%jj(h)+all_orbit%jj(p2))/2+ &
                         all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                         all_orbit%jj(h)-jtot)
                    w1=all_orbit%e(h)+all_orbit%e(h1)+wcn
                    w2=w1+all_orbit%e(h2)-all_orbit%evalence(a)
                    w3=all_orbit%evalence(c)+all_orbit%e(h2)+all_orbit%e(h)- &
                         all_orbit%evalence(a)+wcn
                    de=(den1+wcn)*(den2+wcn)*factr
                    DO i=1, n_startenergy_veff
                       CALL interpolate(w1(i),e_start_g,ans1,val1)
                       CALL interpolate(w2(i),e_start_g,ans2,val2)
                       CALL interpolate(w3(i),e_start_g,ans3,val3)
                       onebody_diagram_21(i)=onebody_diagram_21(i)- &
                            val1*val2*val3*sg/de(i)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE diagram_21
!
!
!     Begin two-body diagrams
!
!
!
!     G-matrix contribution
!
SUBROUTINE diag1(a,b,c,d,jtot,dg1)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, b, c, d, jtot
  INTEGER :: i
  REAL(DP) :: val
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT)  :: dg1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w

  dg1=0.
  CALL pphhmtx(a,b,c,d,jtot,ans)
  IF ( ans (1) == 0.0_dp) RETURN
  IF ( type_of_interaction == 'coupled-cluster') THEN
     w=all_orbit%e(c)+all_orbit%e(d)+wcn
  ELSEIF ( type_of_interaction == 'open-diagrams') THEN
     w=all_orbit%evalence(c)+all_orbit%evalence(d)+wcn
  ENDIF
  DO i=1, n_startenergy_veff
     CALL interpolate(w(i),e_start_g,ans,val)
     dg1(i)=val
  ENDDO

END SUBROUTINE diag1





!
!     Begin closed core diagrams
!
!
!     first-order contribution to the core-energy
!
SUBROUTINE core_E_first_order(first_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, ie
  REAL(DP) :: val, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w
  REAL(DP), INTENT(OUT) :: first_order
  REAL(DP), DIMENSION(n_startenergy_g) :: ans

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           CALL pphhmtx(h1,h2,h1,h2,jtot,ans)
           IF ((ans(1) == 0.)) CYCLE
           w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           DO ie=1,n_startenergy_veff
              CALL interpolate(w(ie),e_start_g,ans,val)
              diagram(ie)=diagram(ie)+fact*val
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  first_order=diagram(n_startenergy_veff/2+1)*0.5

END SUBROUTINE core_E_first_order
!
!  The second order diagram, 2p-2h
!
SUBROUTINE core_E_second_order(second_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, ie, nshell1, nshell2, idiff
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 CALL pphhmtx(p1,p2,h1,h2,jtot,ans1)
                 CALL pphhmtx(h1,h2,p1,p2,jtot,ans2)
                 IF ( ans1(1) == 0.) CYCLE; IF ( ans2(1) == 0.) CYCLE
                 de=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w(ie),e_start_g,ans1,val1)
                    CALL interpolate(w(ie),e_start_g,ans2,val2)
                    diagram(ie)=diagram(ie)+val1*val2*fact/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)*0.25

END SUBROUTINE core_E_second_order


!
!  The second order diagram in  a recoupled way, 2p-2h
!
SUBROUTINE core_E_second_orderrecoupl(second_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min, iph
  INTEGER :: h1, h2, p1, p2, ie, nshell1, nshell2, idiff
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, sg
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(p1))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
        DO jtot=j_min,j_max
           DO h2=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 nshell1=all_orbit%nshell(h1)+all_orbit%nshell(h2)
                 nshell2=all_orbit%nshell(p1)+all_orbit%nshell(p2)
                 idiff=nshell2-nshell1
!                 IF(dencheck(idiff)) CYCLE
                 CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans1)
                 CALL cross_coupled_mtxel1(h1,h2,p1,p2,jtot,ans2)
                 sg=iph( (all_orbit%jj(p1)+all_orbit%jj(p2) +  &
                      3*all_orbit%jj(h1)+3*all_orbit%jj(h2))/2)
                 IF ( ans1(1) == 0.) CYCLE; IF ( ans2(1) == 0.) CYCLE
                 de=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 DO ie=1,n_startenergy_veff
                    CALL interpolate(w(ie),e_start_g,ans1,val1)
                    CALL interpolate(w(ie),e_start_g,ans2,val2)
                    diagram(ie)=diagram(ie)+sg*val1*val2/de(ie)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)*0.25

END SUBROUTINE core_E_second_orderrecoupl
!
!  The second order diagram with two HF insertions
!
SUBROUTINE core_E_second_order_HF(second_order)
  USE constants
  USE single_particle_orbits              
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, p1, ie, jtot
  REAL(DP), INTENT(OUT) :: second_order
  REAL(DP) :: val1, val2, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, de
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits
        IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
        IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        de=all_orbit%e(h2)+wcn-all_orbit%e(p1)
        DO h1=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.))
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              CALL cross_coupled_mtxel1(h2,h3,p1,h3,jtot,ans2)
              IF ( ans2(1) == 0.) CYCLE
              DO ie=1,n_startenergy_veff
                 CALL interpolate(w1(ie),e_start_g,ans1,val1)
                 CALL interpolate(w2(ie),e_start_g,ans2,val2)
                 diagram(ie)=diagram(ie)+val1*val2*fact/de(ie)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  second_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE core_E_second_order_HF
!
!  The third order diagram, 4p-2h
!
SUBROUTINE core_E_third_order1(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, p3, p4, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn-& 
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(h1,h2,p1,p2,jtot,ans1)
                 IF ((ans1(1) == 0.)) CYCLE
                 DO p3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                    DO p4=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(p4) == 'hole' ) CYCLE
                       de2=all_orbit%e(h1)+all_orbit%e(h2)+wcn-&  
                            all_orbit%e(p3)-all_orbit%e(p4)
                       CALL pphhmtx(p1,p2,p3,p4,jtot,ans2)
                       IF(ans2(1) == 0.) CYCLE
                       CALL pphhmtx(p3,p4,h1,h2,jtot,ans3)
                       IF(ans3(1) == 0.) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w(ie),e_start_g,ans1,val1)
                          CALL interpolate(w(ie),e_start_g,ans2,val2)
                          CALL interpolate(w(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+fact*val1*val2*val3/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)*0.125

END SUBROUTINE core_E_third_order1
!
!  The third order diagram, 2p-4h
!
SUBROUTINE core_E_third_order2(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min
  INTEGER :: h1, h2, p1, p2, h3, h4, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(h2))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
        w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
        DO jtot=j_min,j_max
           fact=2*jtot+1.
           DO p1=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&  
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL pphhmtx(p1,p2,h1,h2,jtot,ans1)
                 IF ((ans1(1) == 0.)) CYCLE
                 DO h3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    DO h4=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                       de2=all_orbit%e(h3)+all_orbit%e(h4)+wcn -&
                            all_orbit%e(p1)-all_orbit%e(p2)
                       w2 = w1 + all_orbit%e(h3)+all_orbit%e(h4) -&
                            all_orbit%e(p1)-all_orbit%e(p2)
                       w3 = de2
                       CALL pphhmtx(h1,h2,h3,h4,jtot,ans2)
                       IF(ans2(1) == 0.) CYCLE
                       CALL pphhmtx(h3,h4,p1,p2,jtot,ans3)
                       IF(ans3(1) == 0.) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+fact*val1*val2*val3/de1(ie)/de2(ie) 
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)*0.125

END SUBROUTINE core_E_third_order2
!
!  The third order diagram, 3p-3h
!
SUBROUTINE core_E_third_order3(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: jtot, j_max, j_min, nshell4, iph
  INTEGER :: h1, h2, h3, p1, p2, p3, ie
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact, sg, deno
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.
  DO h1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        j_min=ABS((all_orbit%jj(h1)-all_orbit%jj(p1))/2)
        j_max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
        DO jtot=j_min,j_max
           fact=SQRT(2*jtot+1.)
           DO h2=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 de1=all_orbit%e(h1)+all_orbit%e(h2)+wcn -&
                      all_orbit%e(p1)-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(p1,p2,h1,h2,jtot,ans1) 
                 IF ( ans1(1) == 0.0_dp) CYCLE
                 DO h3=1, all_orbit%total_orbits 
                    IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
                    DO p3=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(p3) == 'hole' ) CYCLE
                       de2=all_orbit%e(h1)+all_orbit%e(h3)+wcn -&
                            all_orbit%e(p1)-all_orbit%e(p3)
                       w2 = all_orbit%e(h3)+all_orbit%e(h1) +wcn
                       w3 = w2 + all_orbit%e(h2)-all_orbit%e(p1)
                       sg=iph((all_orbit%jj(p2)+all_orbit%jj(h2)+ &
                            all_orbit%jj(h3)+all_orbit%jj(h1)+  &
                            all_orbit%jj(p1)+all_orbit%jj(p3))/2+ &
                            all_orbit%jj(h2)+all_orbit%jj(h1)+ &
                            all_orbit%jj(h3)-jtot)
                       CALL cross_coupled_mtxel1(h1,h3,p1,p3,jtot,ans2) 
                       IF ( ans2(1) == 0.0_dp) CYCLE
                       CALL cross_coupled_mtxel1(p3,h2,h3,p2,jtot,ans3) 
                       IF ( ans3(1) == 0.0_dp) CYCLE
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+sg*fact*val1*val2*val3/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE core_E_third_order3
!
!  The third order TDA diagram with two HF insertions
!
SUBROUTINE tda_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)-all_orbit%nshell(p2)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+wcn-all_orbit%e(p2)
                 CALL cross_coupled_mtxel1(h2,p2,p1,h3,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(h3,h4,p2,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+all_orbit%e(h3)+wcn
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE tda_E_third_order_HF
!
!  The third order RPA-1 diagram with two HF insertions
!
SUBROUTINE rpa1_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)+all_orbit%nshell(h2)-&
                      all_orbit%nshell(p2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+all_orbit%e(h2)+wcn-all_orbit%e(p2)-all_orbit%e(p1)
                 CALL cross_coupled_mtxel1(h2,h3,p1,p2,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(p2,h4,h3,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+w2-all_orbit%e(p1)
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE rpa1_E_third_order_HF
!
!  The third order RPA-2 diagram with two HF insertions 
!
SUBROUTINE rpa2_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO h1=1, all_orbit%total_orbits   ! outer hole line
     IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
     DO p1=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
        DO h2=1, all_orbit%total_orbits 
           IF (all_orbit%jj(h2) /= all_orbit%jj(p1)) CYCLE
           IF (all_orbit%ll(h2) /= all_orbit%ll(p1)) CYCLE
           IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
           idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p1)
           IF(dencheck(idiff1)) CYCLE
           de1=all_orbit%e(h2)+wcn-all_orbit%e(p1)
           w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h1,h2,h1,p1,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              DO p2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
                 IF (all_orbit%jj(h3) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(h3) /= all_orbit%ll(p2)) CYCLE
                 idiff2=all_orbit%nshell(h3)+all_orbit%nshell(h2)-&
                      all_orbit%nshell(p2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h3)+all_orbit%e(h2)+wcn-all_orbit%e(p2)-all_orbit%e(p1)
                 CALL cross_coupled_mtxel1(p1,p2,h2,h3,jtot,ans2)
                 IF ( ans2(1) == 0.) CYCLE
                 DO h4=1, all_orbit%total_orbits   ! outer hole line
                    IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
                    fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h4)+1.))
                    CALL cross_coupled_mtxel1(h3,h4,p2,h4,jtot,ans3)
                    IF ( ans3(1) == 0.) CYCLE
                    w3=all_orbit%e(h4)+w2-all_orbit%e(p1)
                    DO ie=1,n_startenergy_veff
                       CALL interpolate(w1(ie),e_start_g,ans1,val1)
                       CALL interpolate(w2(ie),e_start_g,ans2,val2)
                       CALL interpolate(w3(ie),e_start_g,ans3,val3)
                       diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE rpa2_E_third_order_HF
!
!  The third order diagram with 3  HF insertions, 4 holes, 2 particles
!
SUBROUTINE hf31_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, p1, p2, ie, idiff1, idiff2, jtot, jtot2, j_min, j_max
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p2=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p2) == 'hole' ) CYCLE
     DO h2=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
        IF (all_orbit%jj(h2) /= all_orbit%jj(p2)) CYCLE
        IF (all_orbit%ll(h2) /= all_orbit%ll(p2)) CYCLE
        idiff1=all_orbit%nshell(h2)-all_orbit%nshell(p2)
        IF(dencheck(idiff1)) CYCLE
        de1=all_orbit%e(h2)+wcn-all_orbit%e(p2)
        DO h4=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
           w1=all_orbit%e(h4)+all_orbit%e(h2)+wcn
           CALL cross_coupled_mtxel1(h2,h4,p2,h4,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              DO p1=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
                 IF (all_orbit%jj(p1) /= all_orbit%jj(p2)) CYCLE
                 IF (all_orbit%ll(p1) /= all_orbit%ll(p2)) CYCLE
                 w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                 idiff2=all_orbit%nshell(h2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h2)+wcn-all_orbit%e(p1)
                 j_min=ABS((all_orbit%jj(p2)-all_orbit%jj(p1))/2)
                 j_max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
                 DO jtot2 = j_min, j_max
                    CALL cross_coupled_mtxel1(h3,p2,h3,p1,jtot2,ans2)
                    IF ( ans2(1) == 0.) CYCLE
                    DO h1=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                       fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.) &
                            *(all_orbit%jj(h4)+1.)/(all_orbit%jj(p2)+1.))
                       CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot2,ans3)
                       IF ( ans3(1) == 0.) CYCLE
                       w3=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)+val1*val2*val3*fact/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE hf31_E_third_order_HF
!
!  The third order diagram with 3  HF insertions, 5 holes, 1 particle
!
SUBROUTINE hf32_E_third_order_HF(third_order)
  USE constants
  USE single_particle_orbits
  IMPLICIT NONE
  INTEGER :: h1, h2, h3, h4, h5, p1, ie, idiff1, idiff2, jtot, jtot2, j_min, j_max
  REAL(DP), INTENT(OUT) :: third_order
  REAL(DP) :: val1, val2, val3, fact
  REAL(DP), DIMENSION(n_startenergy_veff ) :: diagram, w1, w2, w3, de1, de2
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2, ans3
  LOGICAL dencheck

  diagram=0.; jtot = 0
  DO p1=1, all_orbit%total_orbits 
     IF(all_orbit%orbit_status(p1) == 'hole' ) CYCLE
     DO h4=1, all_orbit%total_orbits 
        IF(all_orbit%orbit_status(h4) /= 'hole' ) CYCLE
        IF (all_orbit%jj(h4) /= all_orbit%jj(p1)) CYCLE
        IF (all_orbit%ll(h4) /= all_orbit%ll(p1)) CYCLE
        idiff1=all_orbit%nshell(h4)-all_orbit%nshell(p1)
        IF(dencheck(idiff1)) CYCLE
        de1=all_orbit%e(h4)+wcn-all_orbit%e(p1)
        DO h5=1, all_orbit%total_orbits 
           IF(all_orbit%orbit_status(h5) /= 'hole' ) CYCLE
           w1=all_orbit%e(h4)+all_orbit%e(h5)+wcn
           CALL cross_coupled_mtxel1(h4,h5,p1,h5,jtot,ans1)
           IF ( ans1(1) == 0.) CYCLE
           DO h3=1, all_orbit%total_orbits 
              IF(all_orbit%orbit_status(h3) /= 'hole' ) CYCLE
              DO h2=1, all_orbit%total_orbits 
                 IF(all_orbit%orbit_status(h2) /= 'hole' ) CYCLE
                 IF (all_orbit%jj(p1) /= all_orbit%jj(h2)) CYCLE
                 IF (all_orbit%ll(p1) /= all_orbit%ll(h2)) CYCLE
                 w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
                 idiff2=all_orbit%nshell(h2)-all_orbit%nshell(p1)
                 IF(dencheck(idiff2)) CYCLE
                 de2=all_orbit%e(h2)+wcn-all_orbit%e(p1)
                 j_min=ABS((all_orbit%jj(h2)-all_orbit%jj(h4))/2)
                 j_max=(all_orbit%jj(h2)+all_orbit%jj(h4))/2
                 DO jtot2 = j_min, j_max
                    CALL cross_coupled_mtxel1(h3,h2,h3,h4,jtot2,ans2)
                    IF ( ans2(1) == 0.) CYCLE
                    DO h1=1, all_orbit%total_orbits 
                       IF(all_orbit%orbit_status(h1) /= 'hole' ) CYCLE
                       fact=SQRT((all_orbit%jj(h1)+1.)*(all_orbit%jj(h3)+1.) &
                            *(all_orbit%jj(h5)+1.)/(all_orbit%jj(p1)+1.))
                       CALL cross_coupled_mtxel1(h1,p1,h1,h2,jtot2,ans3)
                       IF ( ans3(1) == 0.) CYCLE
                       w3=all_orbit%e(h2)+all_orbit%e(h1)+wcn
                       DO ie=1,n_startenergy_veff
                          CALL interpolate(w1(ie),e_start_g,ans1,val1)
                          CALL interpolate(w2(ie),e_start_g,ans2,val2)
                          CALL interpolate(w3(ie),e_start_g,ans3,val3)
                          diagram(ie)=diagram(ie)-val1*val2*val3*fact/de1(ie)/de2(ie)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  third_order=diagram(n_startenergy_veff/2+1)

END SUBROUTINE hf32_E_third_order_HF

!
!
!     Begin three-body diagrams
!
!
!     Sets up three-body contrib with one-body and two-body diagrams only

SUBROUTINE two_2_three(a,b,c,d,e,f,j_ab,j_de,jtot,two_d)
  USE constants
  USE single_particle_orbits
  USE stored_qbox
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, j_ab, j_de, j_ac, jtot, iph, &
       j_ac_min, j_ac_max, j_bc_min, j_bc_max, i
  REAL(DP) :: term, val, dij   
  REAL(DP), DIMENSION(n_startenergy_g) :: q1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w
  REAL(DP),  DIMENSION(n_startenergy_veff ), INTENT(OUT) :: two_d
  REAL(DP), DIMENSION(n_startenergy_veff) :: dg1,dg2,dg3,dg4, &
       dg5,dg6,dg7,dg8,dg9
  LOGICAL triag

  two_d=0.0_dp; dg1=0.0_dp; dg2=0.0_dp; dg3=0.0_dp; dg4=0.0_dp; dg5=0.0_dp; dg6=0.0_dp; dg7=0.0_dp;
  dg8=0.0_dp; dg9=0.0_dp
  j_bc_min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j_bc_max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j_ac_min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j_ac_max=(all_orbit%jj(a)+all_orbit%jj(c))/2

  ! The three-particle starting energy
  w=all_orbit%evalence(d)+all_orbit%evalence(e)+all_orbit%evalence(f)+wcn
  IF(c == f) THEN
     IF(j_ab == j_de) THEN
        CALL pphhmtx(a,b,d,e,j_ab,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg1(i)=dg1(i)+val
        ENDDO
     ENDIF
  ENDIF

  IF(c == d) THEN
     term=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
          all_orbit%jj(f),jtot,2*j_ab)*SQRT((2.*j_ab+1.)*(2.*j_de+1.))&
          *iph((all_orbit%jj(e)+all_orbit%jj(f)+all_orbit%jj(d)+jtot)/2) &
          *iph((ABS(all_orbit%jj(d)-jtot)/2)+j_ab)
     CALL pphhmtx(a,b,e,f,j_ab,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg2(i)=dg2(i)+val*term
     ENDDO
  ENDIF

  IF(c == e) THEN
     term=-sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
          all_orbit%jj(f),2*j_ab)*SQRT((2.*j_ab+1.)*(2.*j_de+1.)) &
          *iph((all_orbit%jj(e)+all_orbit%jj(f))/2+j_ab+j_de)
     CALL pphhmtx(a,b,d,f,j_ab,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg3(i)=dg3(i)+val*term
     ENDDO
  ENDIF

  IF(b == f) THEN
     term=-sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)*SQRT((2.*j_ab+1.)*(2.*j_de+1.))&
          *iph((all_orbit%jj(b)+all_orbit%jj(c))/2+j_de+j_ab)
     CALL pphhmtx(a,c,d,e,j_de,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg4(i)=dg4(i)+val*term
     ENDDO
  ENDIF

  IF(b == d) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(d),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(e),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(d),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(d)+all_orbit%jj(e))/2-j_de+1)
        CALL pphhmtx(a,c,e,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg5(i)=dg5(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(b == e) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac, &
             jtot,all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1)
        CALL pphhmtx(a,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg6(i)=dg6(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(a == f) THEN
     term=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,all_orbit%jj(c), &
          jtot,2*j_de)*SQRT((2.*j_ab+1.)*(2.*j_de+1.)) &
          *iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(a)+jtot)/2)&
          *iph(ABS(all_orbit%jj(a)-jtot)/2+j_de)
     CALL pphhmtx(b,c,d,e,j_de,q1)
     DO i=1, n_startenergy_veff
        CALL interpolate(w(i),e_start_g,q1,val)
        dg7(i)=dg7(i)+val*term
     ENDDO
  ENDIF

  IF(a == d) THEN
     DO j_ac=j_bc_min,j_bc_max
        IF(triag(2*j_ac,all_orbit%jj(a),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(d),jtot)) CYCLE
        term=sjs(all_orbit%jj(b),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(a),2*j_ab) &
             *sjs(all_orbit%jj(e),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(d),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(d)+all_orbit%jj(e)+all_orbit%jj(b)+ &
             all_orbit%jj(a))/2-j_ab-j_de)
        CALL pphhmtx(b,c,e,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg8(i)=dg8(i)+val*term
        ENDDO
     ENDDO
  ENDIF

  IF(a == e) THEN
     DO j_ac=j_bc_min,j_bc_max
        IF(triag(2*j_ac,all_orbit%jj(a),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(b),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(a),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac,jtot, &
             all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1) &
             *iph((all_orbit%jj(b)+all_orbit%jj(a))/2-j_ab+1)
        CALL pphhmtx(b,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg9(i)=dg9(i)+val*term
        ENDDO
     ENDDO
  ENDIF
  two_d=dg1+dg2+dg3+dg4+dg5+dg6+dg7+dg8+dg9

END SUBROUTINE two_2_three


!
!
!     Begin three-body diagrams
!
!
!     Sets up three-body contrib with one-body and two-body diagrams only

SUBROUTINE two_2_threex(a,b,c,d,e,f,j_ab,j_de,jtot,two_d)
  USE constants
  USE single_particle_orbits
  USE stored_qbox
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, j_ab, j_de, j_ac, jtot, iph, j_ac_min, j_ac_max, i
  REAL(DP) :: term, val, dij
  REAL(DP), DIMENSION(n_startenergy_g) :: q1
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w
  REAL(DP),  DIMENSION(n_startenergy_veff ), INTENT(OUT) :: two_d
  REAL(DP), DIMENSION(n_startenergy_veff) :: dg1
  LOGICAL triag

  two_d = 0.0_dp; dg1 = 0.0_dp
  j_ac_min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j_ac_max=(all_orbit%jj(a)+all_orbit%jj(c))/2

  ! The three-particle starting energy
  w=all_orbit%evalence(d)+all_orbit%evalence(e)+all_orbit%evalence(f)+wcn
  
  IF(b == e) THEN
     DO j_ac=j_ac_min,j_ac_max
        IF(triag(2*j_ac,all_orbit%jj(b),jtot)) CYCLE
        IF(triag(2*j_ac,all_orbit%jj(e),jtot)) CYCLE
        term=sjs(all_orbit%jj(a),all_orbit%jj(c),2*j_ac,jtot, &
             all_orbit%jj(b),2*j_ab) &
             *sjs(all_orbit%jj(d),all_orbit%jj(f),2*j_ac, &
             jtot,all_orbit%jj(e),2*j_de) &
             *SQRT((2.*j_ab+1.)*(2.*j_de+1.))*(2.*j_ac+1.) &
             *iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+1)
        CALL pphhmtx(a,c,d,f,j_ac,q1)
        DO i=1, n_startenergy_veff
           CALL interpolate(w(i),e_start_g,q1,val)
           dg1(i)=dg1(i)+val*term
        ENDDO
     ENDDO
  ENDIF
  two_d = dg1

END SUBROUTINE two_2_threex



!
!              Three-body diagram 1 with particle intermediate state
!

SUBROUTINE three_body_1a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(b)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2)
     CALL pphhmtx(a,p,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)
        sixj2=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_de, &
             all_orbit%jj(f),jtot,2*j_bc)
        fact2=fact1*sixj1*sixj2*(2.*j_bc+1.)
        CALL pphhmtx(b,c,p,f,j_bc,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_1a
!
!              Three-body diagram 1 with hole intermediate state
!
SUBROUTINE three_body_1b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))*iph((all_orbit%jj(p)+ &
          all_orbit%jj(b)+all_orbit%jj(f)+all_orbit%jj(c))/2+all_orbit%jj(p))
     CALL pphhmtx(a,p,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(p)+all_orbit%e(f)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(p)+all_orbit%e(f)+all_orbit%e(d)+all_orbit%e(e)- &
          all_orbit%e(b)-all_orbit%e(c)+wcn
     w2=all_orbit%e(p)+all_orbit%e(f)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c), jtot,2*j_bc)
        sixj2=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_de, &
             all_orbit%jj(f),jtot,2*j_bc)
        fact2=fact1*sixj1*sixj2*(2.*j_bc+1.)
        CALL pphhmtx(b,c,p,f,j_bc,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_1b



SUBROUTINE three_body_2a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab-j_de))
     CALL pphhmtx(p,b,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     IF(ANS1(2) == 0.)CYCLE
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)-all_orbit%e(b)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)
        sixj2=sjs(all_orbit%jj(p),all_orbit%jj(b),2*j_de,jtot, &
             all_orbit%jj(f),2*j_ac)
        fact2=fact1*sixj1*sixj2*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,f,j_ac,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_2a




SUBROUTINE three_body_2b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, &
       sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab- &
          j_de+all_orbit%jj(p)))
     CALL pphhmtx(p,b,d,e,j_de,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(f)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(p)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)
        sixj2=sjs(all_orbit%jj(p),all_orbit%jj(b),2*j_de,jtot, &
             all_orbit%jj(f),2*j_ac)
        fact2=fact1*sixj1*sixj2*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,f,j_ac,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        IF(ans2(2) == 0.) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE three_body_2b


SUBROUTINE three_body_3a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, iph,  &
       nshell1, nshell2, idiff, jtot
  REAL(DP) :: val1, val2, fact                         
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2 
  REAL(DP), DIMENSION(n_startenergy_veff) :: den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de)* &
          sjs(all_orbit%jj(p),all_orbit%jj(f),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)
     CALL pphhmtx(a,b,p,f,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     CALL pphhmtx(p,c,d,e,j_de,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
     den=wcn+all_orbit%e(d)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(c)
     w2=all_orbit%e(d)+all_orbit%e(e)+wcn
     w1=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+ &
          all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        CALL interpolate(w2(ie),e_start_g,ans2,val2)
        three_diagram(ie)=three_diagram(ie)+fact*val1*val2/den(ie)
     ENDDO
  ENDDO

END SUBROUTINE three_body_3a




SUBROUTINE three_body_3b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, iph,  &
       nshell1, nshell2, idiff, jtot
  REAL(DP) :: val1, val2, fact                          
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(a)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(f)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_ab+j_de+ &
          all_orbit%jj(p))* &
          sjs(all_orbit%jj(p),all_orbit%jj(f),2*j_ab,jtot, &
          all_orbit%jj(c),2*j_de)
     CALL pphhmtx(a,b,p,f,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     CALL pphhmtx(p,c,d,e,j_de,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
     den=wcn+all_orbit%e(p)+all_orbit%e(f)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(p)+all_orbit%e(f)+wcn
     w2=all_orbit%e(f)-all_orbit%e(a)+all_orbit%e(d)+all_orbit%e(e)- &
          all_orbit%e(b)+all_orbit%e(p)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        CALL interpolate(w2(ie),e_start_g,ans2,val2)
        three_diagram(ie)=three_diagram(ie)+fact*val1*val2/den(ie)
     ENDDO
  ENDDO

END SUBROUTINE three_body_3b


SUBROUTINE three_body_4a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(f)+all_orbit%nshell(e)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2)
     CALL pphhmtx(a,b,d,p,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(f)+all_orbit%e(e)-all_orbit%e(p)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ef=j1min,j1max
        sixj1=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_ef)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             all_orbit%jj(f),jtot,2*j_ef)
        fact2=fact1*sixj1*sixj2*(2.*j_ef+1.)
        CALL pphhmtx(p,c,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_4a


SUBROUTINE three_body_4b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(p)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+all_orbit%jj(p))
     CALL pphhmtx(a,b,d,p,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+all_orbit%e(d)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(b)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_ef=j1min,j1max
        sixj1=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_ef)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             all_orbit%jj(f),jtot,2*j_ef)
        fact2=fact1*sixj1*sixj2*(2.*j_ef+1.)
        CALL pphhmtx(p,c,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_4b


SUBROUTINE three_body_5a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(c)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+j_ab-j_de))
     CALL pphhmtx(a,b,p,e,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(c)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_df=j1min,j1max
        sixj1=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_df)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
             all_orbit%jj(f),2*j_df)
        fact2=fact1*sixj1*sixj2*(2.*j_df+1.)
        CALL pphhmtx(p,c,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_5a




SUBROUTINE three_body_5b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j1max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph(ABS((all_orbit%jj(c)-all_orbit%jj(f))/2+ &
          j_ab-j_de+all_orbit%jj(p)))
     CALL pphhmtx(a,b,p,e,j_ab,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(a)-all_orbit%e(b)
     w1=all_orbit%e(e)+all_orbit%e(p)+wcn
     w2=all_orbit%e(e)+all_orbit%e(p)+all_orbit%e(f)+ &
          all_orbit%e(d)-all_orbit%e(a)-all_orbit%e(b)+wcn
     DO ie=1,n_startenergy_veff
        CALL interpolate(w1(ie),e_start_g,ans1,val1)
        vert1(ie)=val1
     ENDDO
     DO j_df=j1min,j1max
        sixj1=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ab, &
             jtot,all_orbit%jj(c),2*j_df)
        sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
             jtot,all_orbit%jj(f),2*j_df)
        fact2=fact1*sixj1*sixj2*(2.*j_df+1.)
        CALL pphhmtx(p,c,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
        DO ie=1,n_startenergy_veff
           CALL interpolate(w2(ie),e_start_g,ans2,val2)
           three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_5b



SUBROUTINE three_body_6a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot, j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+j_de)
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(a)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,e,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                jtot,all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_df, &
                all_orbit%jj(e),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)*iph(j_bc+j_df)
           CALL pphhmtx(a,p,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_6a



SUBROUTINE three_body_6b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_df,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(b)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+j_de+all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(p)+all_orbit%e(e)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+all_orbit%e(e)+ &
          all_orbit%e(p)-all_orbit%e(b)-all_orbit%e(c)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,e,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                jtot,all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_df, &
                all_orbit%jj(e),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)*iph(j_bc+j_df)
           CALL pphhmtx(a,p,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_6b


SUBROUTINE three_body_7a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(d)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_de+j_ab)
     den=wcn+all_orbit%e(d)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(b)+wcn
     w2=all_orbit%e(f)+all_orbit%e(d)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,e,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
                all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ac,jtot, &
                all_orbit%jj(b),2*j_df)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)
           CALL pphhmtx(p,b,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_7a


SUBROUTINE three_body_7b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_df, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3 
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(d)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(d)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(c)+all_orbit%nshell(a)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(e)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(c)+all_orbit%jj(f))/2+j_de+j_ab+ &
          all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(e)-all_orbit%e(c)-all_orbit%e(a)
     w1=all_orbit%e(p)+all_orbit%e(e)+wcn
     w2=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(c)-all_orbit%e(a)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,p,e,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_df=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,jtot, &
                all_orbit%jj(f),2*j_df)
           sixj3=sjs(all_orbit%jj(p),all_orbit%jj(e),2*j_ac,jtot, &
                all_orbit%jj(b),2*j_df)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_df+1.)
           CALL pphhmtx(p,b,d,f,j_df,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_7b


SUBROUTINE three_body_8a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_ef,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(e)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(a)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+jtot)
     den=wcn+all_orbit%e(e)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(a)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(a)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,d,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_ef, &
                all_orbit%jj(d),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_bc+j_ef)
           CALL pphhmtx(a,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_8a



SUBROUTINE three_body_8b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_bc, iph, j_ef,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(b)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(b)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(b)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(e)+all_orbit%jj(f)+ &
          all_orbit%jj(c))/2+jtot+all_orbit%jj(p))
     den=wcn+all_orbit%e(p)+all_orbit%e(d)-all_orbit%e(b)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)+ &
          all_orbit%e(p)-all_orbit%e(b)-all_orbit%e(c)+wcn
     DO j_bc=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab, &
             all_orbit%jj(c),jtot,2*j_bc)*(2.*j_bc+1.)
        CALL pphhmtx(b,c,d,p,j_bc,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(a),all_orbit%jj(p),2*j_ef, &
                all_orbit%jj(d),jtot,2*j_bc)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_bc+j_ef)
           CALL pphhmtx(a,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_8b



SUBROUTINE three_body_9a(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF( all_orbit%orbit_status(p) /= 'particle') CYCLE
     nshell1=all_orbit%nshell(e)+all_orbit%nshell(f)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(b)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+j_ab)
     den=wcn+all_orbit%e(e)+all_orbit%e(f)-all_orbit%e(p)-all_orbit%e(b)
     w1=all_orbit%e(d)+all_orbit%e(e)+all_orbit%e(f)- &
          all_orbit%e(b)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+wcn
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,d,p,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de, &
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ac, &
                all_orbit%jj(b),jtot,2*j_ef)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_ac+j_ef)
           CALL pphhmtx(b,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE three_body_9a


SUBROUTINE three_body_9b(a,b,c,d,e,f,j_ab,j_de,jtot,three_diagram)
  USE constants
  USE single_particle_orbits
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER :: a, b, c, d, e, f, ie, p, j_ab, j_de, j_ac, j_ef, iph,  &
       nshell1, nshell2, idiff, j1min, j1max, jtot,  j2min, j2max
  REAL(DP) :: val1, val2, fact1, fact2, sixj1, sixj2, sixj3
  REAL(DP), DIMENSION(n_startenergy_veff), INTENT(OUT) :: three_diagram 
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  REAL(DP), DIMENSION(n_startenergy_veff) :: vert1, den, w1, w2
  LOGICAL dencheck

  three_diagram=0.
  j1min=ABS((all_orbit%jj(a)-all_orbit%jj(c))/2)
  j1max=(all_orbit%jj(a)+all_orbit%jj(c))/2
  j2min=ABS((all_orbit%jj(e)-all_orbit%jj(f))/2)
  j2max=(all_orbit%jj(e)+all_orbit%jj(f))/2
  DO p=1,all_orbit%total_orbits
     IF(all_orbit%orbit_status(p) /= 'hole') CYCLE
     nshell1=all_orbit%nshell(a)+all_orbit%nshell(c)
     nshell2=all_orbit%nshell(p)+all_orbit%nshell(d)
     idiff=nshell1-nshell2
     IF(dencheck(idiff))CYCLE
     fact1=-SQRT((2.*j_ab+1.)*(2.*j_de+1.))* &
          iph((all_orbit%jj(b)+all_orbit%jj(c)+all_orbit%jj(f)+ &
          all_orbit%jj(e))/2+j_ab+all_orbit%jj(p))
     den=wcn+all_orbit%e(d)+all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)
     w1=all_orbit%e(d)+all_orbit%e(p)+wcn
     w2=all_orbit%e(f)+all_orbit%e(e)+all_orbit%e(d)+ &
          all_orbit%e(p)-all_orbit%e(a)-all_orbit%e(c)+wcn         
     DO j_ac=j1min,j1max
        sixj1=sjs(all_orbit%jj(a),all_orbit%jj(b),2*j_ab,jtot, &
             all_orbit%jj(c),2*j_ac)*(2.*j_ac+1.)
        CALL pphhmtx(a,c,d,p,j_ac,ans1) ; IF (ans1(1) == 0.0_dp) CYCLE 
        DO ie=1,n_startenergy_veff
           CALL interpolate(w1(ie),e_start_g,ans1,val1)
           vert1(ie)=val1
        ENDDO
        DO j_ef=j2min,j2max
           sixj2=sjs(all_orbit%jj(d),all_orbit%jj(e),2*j_de,&
                all_orbit%jj(f),jtot,2*j_ef)
           sixj3=sjs(all_orbit%jj(d),all_orbit%jj(p),2*j_ac, &
                all_orbit%jj(b),jtot,2*j_ef)
           fact2=fact1*sixj1*sixj2*sixj3*(2.*j_ef+1.)*iph(j_ac+j_ef)
           CALL pphhmtx(b,p,e,f,j_ef,ans2) ; IF (ans2(1) == 0.0_dp) CYCLE
           DO ie=1,n_startenergy_veff
              CALL interpolate(w2(ie),e_start_g,ans2,val2)
              three_diagram(ie)=three_diagram(ie)+fact2*vert1(ie)*val2/ &
                   den(ie)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE three_body_9b
!
!     Begin one-body effective operator diagrams
!
!
!     effective operator diagram with ph intermdediate state, first order in G
!
SUBROUTINE effective_operator_1(a,c,effoper_diagram_1)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, l, p,  nshell1, nshell2, idiff, iphase, iph
  REAL(DP) :: val, factr
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_1
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  LOGICAL dencheck
  effoper_diagram_1=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
        nshell1=all_orbit%nshell(a)+all_orbit%nshell(p)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(h)
        idiff=nshell1-nshell2
        IF(dencheck(idiff)) CYCLE
        de=all_orbit%evalence(c)+all_orbit%e(h)-all_orbit%evalence(a)-all_orbit%e(p)+wcn
        IF(bare_operator(h,p) == 0.) CYCLE
        iphase=iph((3*all_orbit%jj(h)+all_orbit%jj(p))/2-lambda)
        factr=iphase*bare_operator(h,p)/SQRT(2.*lambda+1.)
        CALL cross_coupled_mtxel1(a,p,c,h,lambda,ans)
        IF(ans(2) == 0.) CYCLE
        w=all_orbit%evalence(c)+all_orbit%e(h)+wcn
        DO l=1,n_startenergy_veff
           CALL interpolate(w(l),e_start_g,ans,val)
           effoper_diagram_1(l)=effoper_diagram_1(l)+val*factr/de(l)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_1
!
!     effective operator diagram with ph intermdediate state, first order in G
!     hermitian conjugate of previous diagram
!
SUBROUTINE effective_operator_2(a,c,effoper_diagram_2)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, l, p,  nshell1, nshell2, idiff, iphase, iph
  REAL(DP) :: val, factr
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  LOGICAL dencheck

  effoper_diagram_2=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
        nshell1=all_orbit%nshell(a)+all_orbit%nshell(p)
        nshell2=all_orbit%nshell(c)+all_orbit%nshell(h)
        idiff=nshell1-nshell2
        IF(dencheck(idiff)) CYCLE
        IF(bare_operator(p,h) == 0.) CYCLE
        iphase=iph((3*all_orbit%jj(h)+all_orbit%jj(p))/2-lambda)
        factr=iphase*bare_operator(p,h)/SQRT(2.*lambda+1.)
        CALL cross_coupled_mtxel1(a,h,c,p,lambda,ans)
        IF(ans(2) == 0.) CYCLE
        de=all_orbit%evalence(a)+all_orbit%e(h)-all_orbit%evalence(c)-all_orbit%e(p)+wcn
        w=all_orbit%evalence(a)+all_orbit%e(h)+wcn
        DO l=1,n_startenergy_veff
           CALL interpolate(w(l),e_start_g,ans,val)
           effoper_diagram_2(l)=effoper_diagram_2(l)+val*factr/de(l)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_2
!
!
!
SUBROUTINE effective_operator_3(a,c,effoper_diagram_3)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: w1, w2, de
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_3
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_3=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(c)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%evalence(a)-all_orbit%e(p1)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,p2,p1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1+all_orbit%e(h2)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_3(l)=effoper_diagram_3(l)+& 
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_3
!
!
!
SUBROUTINE effective_operator_4(a,c,effoper_diagram_4)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_4
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_4=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(a)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p1)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,h2,h1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              w2=w1+all_orbit%e(h2)-all_orbit%evalence(c)
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_4(l)=effoper_diagram_4(l)+&
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_4
!
!
!
subroutine effective_operator_5(a,c,effoper_diagram_5)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff )  :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_5
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_5=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE    
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE 
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)-& 
                   all_orbit%evalence(a)-all_orbit%e(p1)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-& 
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,h2,p1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3* &
                   all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                   all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_5(l)=effoper_diagram_5(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_5
!
!
!
subroutine effective_operator_6(a,c,effoper_diagram_6)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_6
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_6=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p1)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,p2,h1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3* &
                   all_orbit%jj(h2)+all_orbit%jj(p1)+ &
                   all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_6(l)=effoper_diagram_6(l)+& 
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_6
!
!
!
subroutine effective_operator_7(a,c,effoper_diagram_7)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_7
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_7=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE         
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE            
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(p2,h2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h2)+all_orbit%evalence(a)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,p1,c,h1,lambda,ans1)
              CALL cross_coupled_mtxel1(h1,h2,p1,p2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(p2,h2)/(2.*lambda+1.)
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+ &
                   all_orbit%e(h2)-all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_7(l)=effoper_diagram_7(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_7
!
!
!
SUBROUTINE effective_operator_8(a,c,effoper_diagram_8)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph
  REAL(DP) :: val1, val2, factr, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_8
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck, triag

  effoper_diagram_8=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE 
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE   
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              IF(bare_operator(h2,p2) == 0.) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p2),all_orbit%jj(h2))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(p1),all_orbit%jj(h1))) CYCLE
              IF(triag(2*lambda,all_orbit%jj(c),all_orbit%jj(a))) CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell4=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h2)+all_orbit%evalence(c)- &
                   all_orbit%evalence(a)-all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              CALL cross_coupled_mtxel1(a,h1,c,p1,lambda,ans1)
              CALL cross_coupled_mtxel1(p1,p2,h1,h2,lambda,ans2)
              IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
              factr=iph((3*all_orbit%jj(h1)+3*all_orbit%jj(h2)+&
                   all_orbit%jj(p1)+all_orbit%jj(p2))/2) &
                   *bare_operator(h2,p2)/(2.*lambda+1.)
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w1=all_orbit%evalence(c)-all_orbit%evalence(a)+w2
              de= (den1 + wcn)*(den2 + wcn)
              DO l=1,n_startenergy_veff
                 CALL interpolate(w1(l),e_start_g,ans1,val1)
                 CALL interpolate(w2(l),e_start_g,ans2,val2)
                 effoper_diagram_8(l)=effoper_diagram_8(l)+ &
                      factr*val1*val2/de(l)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_8
!
!
!
SUBROUTINE effective_operator_9(a,c,effoper_diagram_9)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_9
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_9=0.
  DO p3=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE 
     DO p2=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE        
        IF(bare_operator(p3,p2) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p3)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(h)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%e(h)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p3)
              w1=all_orbit%e(h)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(h)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(h)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p3), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,h,p1,p3,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,p2,c,h,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((3*all_orbit%jj(p3)+3*all_orbit%jj(h)&
                      +all_orbit%jj(p2)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p1))/2)*sixjj*iph(lambda-jt)&
                      *bare_operator(p3,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_9(l)=effoper_diagram_9(l)+&
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_9
!
!
!
SUBROUTINE effective_operator_10(a,c,effoper_diagram_10)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_10
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_10=0.
  DO h1=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE         
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h1,h2) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO p2=1, all_orbit%total_orbits                  
              IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%e(h2)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%evalence(a)+wcn
              w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(h1),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-0.5D0*iph((all_orbit%jj(h2)+&
                      all_orbit%jj(a))/2+jt)*sixjj*(2.*jt+1.) &
                      *bare_operator(h1,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_10(l)=effoper_diagram_10(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_10
!
!
!
SUBROUTINE effective_operator_11(a,c,effoper_diagram_11)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p, h1, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_11
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_11=0.
  DO h3=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h2,h3) == 0.) CYCLE
        DO p=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(h3)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(p)
              idiff1=nshell1-nshell3
              idiff2=nshell4-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p)-all_orbit%evalence(a)
              den2=all_orbit%e(h1)+all_orbit%e(h3)- &
                   all_orbit%e(p)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h3))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h3))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h3)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local      
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(h3),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,p,h1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,h3,c,p,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((3*all_orbit%jj(h1)+3* &
                      all_orbit%jj(h2)+all_orbit%jj(p)&
                      +2*all_orbit%jj(a)+&
                      all_orbit%jj(h3))/2)*sixjj*iph(lambda-jt) &
                      *bare_operator(h2,h3)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_11(l)=effoper_diagram_11(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_11
!
!
!
SUBROUTINE effective_operator_12(a,c,effoper_diagram_12)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p1, p2,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_12
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_12=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     DO p2=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        IF(bare_operator(p2,p1) == 0.0) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell3-nshell2
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p2)
              den2=all_orbit%e(h2)+all_orbit%e(h1)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p1), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,p1,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=iph((all_orbit%jj(p1)+all_orbit%jj(a))/2+jt)
                 factr=-0.5D0*factr*(2.*jt+1.)*sixjj*bare_operator(p2,p1)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_12(l)=effoper_diagram_12(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_12
!
!
!
SUBROUTINE effective_operator_13(a,c,effoper_diagram_13)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_13
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_13=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE	 
     DO p3=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE	   
        IF(bare_operator(h,p3) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO p2=1, all_orbit%total_orbits
              IF(all_orbit%orbit_status(p2) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(p3)
              idiff1=nshell1-nshell2
              idiff2=nshell1-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h)- &
                   all_orbit%e(p3)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(c)+all_orbit%e(h)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p3), &
                      all_orbit%jj(h),2*jt)
                 IF(sixjj == 0.) CYCLE
                 factr=0.5D0*iph(jt+(all_orbit%jj(a)+ &
                      all_orbit%jj(p3))/2)*sixjj &
                      *(2.*jt+1.)*bare_operator(h,p3)
                 CALL pphhmtx(a,p3,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_13(l)=effoper_diagram_13(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_13
!
!
!
SUBROUTINE effective_operator_14(a,c,effoper_diagram_14)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h, p1, l, p2, p3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_14
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_14=0.
  DO h=1,all_orbit%total_orbits
     IF (all_orbit%orbit_status(h) /= 'hole') CYCLE
     DO p3=1, all_orbit%total_orbits
        IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
        IF(bare_operator(p3,h) == 0.) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE 
           DO p2=1, all_orbit%total_orbits
              IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE 
              nshell1=all_orbit%nshell(h)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p3)
              idiff1=nshell1-nshell2
              idiff2=nshell1-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h)+all_orbit%evalence(a)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h)- &
                   all_orbit%e(p3)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p3))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(h))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(p3))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(a)+all_orbit%e(h)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h), &
                      all_orbit%jj(p3),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL pphhmtx(a,h,p1,p2,jt,ans1)
                 IF(ans1(2) == 0.) CYCLE
                 CALL pphhmtx(p1,p2,c,p3,jt,ans2)
                 IF(ans2(2) == 0.) CYCLE
                 factr=iph(jt+(all_orbit%jj(a)+3*all_orbit%jj(h)&
                      +2*all_orbit%jj(p3))/2)*sixjj &
                      *0.5*(2.*jt+1.)*bare_operator(p3,h)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_14(l)=effoper_diagram_14(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_14
!
!
!
SUBROUTINE effective_operator_15(a,c,effoper_diagram_15)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p, h1, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_15
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_15=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
     DO h3=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
        IF(bare_operator(h3,p) == 0.0D0) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(h3)
              idiff1=nshell2-nshell1
              idiff2=nshell3-nshell1
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(a)-all_orbit%e(p)
              den2=all_orbit%evalence(c)-all_orbit%e(p)+ &
                   all_orbit%e(h3)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h3))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(p))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(h3))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h1)+wcn
              w2=w1+all_orbit%evalence(c)+all_orbit%e(h3)- &
                   all_orbit%e(p)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p), &
                      all_orbit%jj(h3),2*jt)
                 if(sixjj == 0.0) cycle
                 CALL pphhmtx(a,p,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,h3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=0.5*iph(jt+(all_orbit%jj(a)+ &
                      all_orbit%jj(p))/2)*sixjj*(2.*jt+1.)*bare_operator(h3,p)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_15(l)=effoper_diagram_15(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_15
!
!
!
SUBROUTINE effective_operator_16(a,c,effoper_diagram_16)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, p, l, h2, h3,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_16
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_16=0.
  DO p=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p) /= 'particle') CYCLE
     DO h3=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
        IF(bare_operator(p,h3) == 0.) CYCLE
        DO h1=1,all_orbit%total_orbits
           IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h3)
              idiff1=nshell2-nshell1
              idiff2=nshell3-nshell1
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%evalence(c)-all_orbit%e(p)
              den2=all_orbit%evalence(a)-all_orbit%e(p)+ &
                   all_orbit%e(h3)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h3))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h3))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%evalence(a)+all_orbit%e(h2)+ &
                   all_orbit%e(h3)+all_orbit%e(h1)- &
                   all_orbit%e(p)-all_orbit%evalence(c)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h3), &
                      all_orbit%jj(p),2*jt)
                 IF(sixjj == 0.) CYCLE
                 factr=iph(jt+(2*all_orbit%jj(p)+ &
                      3*all_orbit%jj(h3)+all_orbit%jj(a))/2)*sixjj &
                      *0.5*(2.*jt+1.)*bare_operator(p,h3)
                 CALL pphhmtx(a,h3,h1,h2,jt,ans1)
                 CALL pphhmtx(h1,h2,c,p,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_16(l)=effoper_diagram_16(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_16
!
!
!
SUBROUTINE effective_operator_17(a,c,effoper_diagram_17)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, h2, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_17
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_17=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(h2,p2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(c)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(a)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1+all_orbit%evalence(c)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.0) CYCLE
                 CALL cross_coupled_mtxel2(a,p1,h1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,p2,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p2)+3*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(h2,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_17(l)=effoper_diagram_17(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_17
!
!
!
SUBROUTINE effective_operator_18(a,c,effoper_diagram_18)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_18
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_18=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE      
        IF(bare_operator(p2,h2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(p2)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(h2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)- &
                   all_orbit%e(p1)-all_orbit%evalence(c)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w1=w2+all_orbit%evalence(a)-all_orbit%evalence(c)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(a),all_orbit%jj(c), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,p1,h1,p2,jt,ans1)
                 CALL cross_coupled_mtxel2(h1,h2,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(-lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +2*all_orbit%jj(c)+3*all_orbit%jj(p2)+&
                      3*all_orbit%jj(h2))/2)*sixjj*bare_operator(p2,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_18(l)=effoper_diagram_18(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_18
!
!
!
SUBROUTINE effective_operator_19(a,c,effoper_diagram_19)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_19
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_19=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE            
        IF(bare_operator(h2,p2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(c)+all_orbit%nshell(h2)
              nshell4=all_orbit%nshell(a)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)- &
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(c)-all_orbit%e(p2)+wcn
              w2=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(h2), &
                      all_orbit%jj(p2),2*jt)
                 IF(sixjj == 0.) CYCLE
                 CALL cross_coupled_mtxel2(a,h1,p1,h2,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,p2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +all_orbit%jj(p2)+3*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(h2,p2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_19(l)=effoper_diagram_19(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_19
!
!
!
SUBROUTINE effective_operator_20(a,c,effoper_diagram_20)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, nshell4, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr, den1, den2, sixjj
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_20
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_20=0.
  DO p2=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        IF(bare_operator(p2,h2) == 0.0D0) CYCLE
        DO p1=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p1)
              nshell3=all_orbit%nshell(a)+all_orbit%nshell(h2)
              nshell4=all_orbit%nshell(c)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell4-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)-& 
                   all_orbit%e(p1)-all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)- &
                   all_orbit%e(p2)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              j2max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(JMIN_LOCAL > JMAX_LOCAL) CYCLE
              w2=all_orbit%e(h1)+all_orbit%e(h2)+ &
                   all_orbit%evalence(a)-all_orbit%e(p2)+wcn
              w1=all_orbit%e(h1)+all_orbit%evalence(a)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 sixjj=sjs(all_orbit%jj(c),all_orbit%jj(a), &
                      2*lambda,all_orbit%jj(p2), &
                      all_orbit%jj(h2),2*jt)
                 IF(sixjj == 0.0D0) CYCLE
                 CALL cross_coupled_mtxel2(a,h1,p1,p2,jt,ans1)
                 CALL cross_coupled_mtxel2(p1,h2,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=-iph(-lambda-jt+(all_orbit%jj(p1)+ &
                      3*all_orbit%jj(h1)+2*all_orbit%jj(a) &
                      +3*all_orbit%jj(p2)+5*all_orbit%jj(h2))/2)*sixjj &
                      *bare_operator(p2,h2)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_20(l)=effoper_diagram_20(l)&
                         +factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_20
!
!
!
SUBROUTINE effective_operator_21(a,c,effoper_diagram_21)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_21
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_21=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.) CYCLE
     fact=-0.5/(all_orbit%jj(c)+1.)*bare_operator(a,p1)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE	    
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE  
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,c,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_21(l)=effoper_diagram_21(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_21
!
!
!
SUBROUTINE effective_operator_22(a,c,effoper_diagram_22)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_22
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_22=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.) CYCLE
     fact=-bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.) * &
          iph((all_orbit%jj(p1)-all_orbit%jj(a)))
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell2
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(a,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_22(l)=effoper_diagram_22(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_22
!
!
!
SUBROUTINE effective_operator_23(a,c,effoper_diagram_23)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_23
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_23=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((-all_orbit%jj(a)+all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(a)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_23(l)=effoper_diagram_23(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_23
!
!
!
SUBROUTINE effective_operator_24(a,c,effoper_diagram_24)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_24
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_24=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_24(l)=effoper_diagram_24(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_24
!
!
!
SUBROUTINE effective_operator_23folded(a,c,effoper_diagram_23f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_23f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_23f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.25/(all_orbit%jj(a)+1.)* &
          iph((-all_orbit%jj(a)+all_orbit%jj(h1)))
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(a)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_23f(l)=effoper_diagram_23f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_23folded
!
!
!
SUBROUTINE effective_operator_24folded(a,c,effoper_diagram_24f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_24f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_24f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.25/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h2)+all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1-nshell3
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)+all_orbit%e(h2)-all_orbit%e(p2)- &
                   all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=all_orbit%evalence(c)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_24f(l)=effoper_diagram_24f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_24folded

!
!
!
SUBROUTINE effective_operator_25(a,c,effoper_diagram_25)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, p3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_25
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_25=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.) CYCLE
     fact=bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(p1)-all_orbit%jj(a)))
     DO p2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        DO p3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p3)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(a)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(a)-all_orbit%e(p3)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h1))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p3))/2
              j3min=ABS(all_orbit%jj(p1)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(p1)+all_orbit%jj(h1))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p3))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p2,p3,p1,h1,jt,ans1)
                 CALL pphhmtx(a,h1,p2,p3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_25(l)=effoper_diagram_25(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_25
!
!
!
SUBROUTINE effective_operator_26(a,c,effoper_diagram_26)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, p2, l, p3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_26
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_26=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.) CYCLE
     fact=bare_operator(a,p1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO p2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
        DO p3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p3) /= 'particle') CYCLE
           DO h1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(p2)+all_orbit%nshell(p3)
              nshell3=all_orbit%nshell(h1)+all_orbit%nshell(c)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%evalence(c)-all_orbit%e(p3)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h1))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p3))/2
              j3min=ABS(all_orbit%jj(p1)-all_orbit%jj(h1))/2
              j1max=(all_orbit%jj(p1)+all_orbit%jj(h1))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p3))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%evalence(c)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)  		  
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,h1,p2,p3,jt,ans1)
                 CALL pphhmtx(p2,p3,c,h1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_26(l)=effoper_diagram_26(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_26
!
!
!
SUBROUTINE effective_operator_27(a,c,effoper_diagram_27)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_27
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_27=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(h1)-all_orbit%jj(a))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_27(l)=effoper_diagram_27(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_27
!
!
!
SUBROUTINE effective_operator_28(a,c,effoper_diagram_28)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_28
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_28=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_28(l)=effoper_diagram_28(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_28
!
!
!
SUBROUTINE effective_operator_27folded(a,c,effoper_diagram_27f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_27f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_27f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=-bare_operator(h1,c)*0.25/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(h1)-all_orbit%jj(a))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans1)
                 CALL pphhmtx(a,h2,p1,p2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_27f(l)=effoper_diagram_27f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_27folded
!
!
!
SUBROUTINE effective_operator_28folded(a,c,effoper_diagram_28f)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_28f
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_28f=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%model_space(h1) == 'outside') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=-bare_operator(a,h1)*0.25/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO p1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(h2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(h2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,c,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_28f(l)=effoper_diagram_28f(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_28folded
!
!
!
SUBROUTINE effective_operator_29(a,c,effoper_diagram_29)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_29
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_29=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(a,p1) == 0.0D0) CYCLE
     fact=-0.5/(all_orbit%jj(c)+1.)*bare_operator(a,p1)* &
          iph((all_orbit%jj(c)-all_orbit%jj(p1))/2)
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE	    
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE  
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(c)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+all_orbit%evalence(c)&
                   -all_orbit%e(p1)+wcn
              w2=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,c,p2,jt,ans1)
                 CALL pphhmtx(p1,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_29(l)=effoper_diagram_29(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_29
!
!
!
SUBROUTINE effective_operator_30(a,c,effoper_diagram_30)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: h1, h2, l, p2, p1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, iph, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_30
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_30=0.
  DO p1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(p1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(p1))) CYCLE
     IF(bare_operator(p1,c) == 0.0D0) CYCLE
     fact=-bare_operator(p1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(a)-all_orbit%jj(p1)))
     DO h1=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
        DO p2=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(p2) /= 'particle') CYCLE
           DO h2=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
              nshell1=all_orbit%nshell(p1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h1)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(p2)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=all_orbit%e(h1)+all_orbit%e(h2)-all_orbit%e(p1)- &
                   all_orbit%e(p2)
              den2=all_orbit%evalence(a)-all_orbit%e(p1)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p2))/2
              j2min=ABS(all_orbit%jj(p2)-all_orbit%jj(p1))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(h2))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(h2))/2
              j2max=(all_orbit%jj(p2)+all_orbit%jj(p1))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p2))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h1)+all_orbit%e(h2)+wcn
              w2=w1+all_orbit%evalence(a)-all_orbit%e(p1)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,h2,p1,p2,jt,ans1)
                 CALL pphhmtx(a,p2,h1,h2,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_30(l)=effoper_diagram_30(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_30
!
!
!
SUBROUTINE effective_operator_31(a,c,effoper_diagram_31)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, h2, l, h3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_31
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_31=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(a) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(a) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(h1,c) == 0.) CYCLE
     fact=bare_operator(h1,c)*0.5/(all_orbit%jj(a)+1.)* &
          iph((all_orbit%jj(a)-all_orbit%jj(h1))/2)
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        DO h3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
           DO p1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(a)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h3)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(a)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=-(all_orbit%e(p1)+all_orbit%evalence(a)-all_orbit%e(h3)- &
                   all_orbit%e(h2))
              den2=all_orbit%e(h1)-all_orbit%evalence(a)
              j1min=ABS(all_orbit%jj(a)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h3))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h3))/2
              j3max=(all_orbit%jj(a)+all_orbit%jj(p1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              w1=w2+all_orbit%e(h1)-all_orbit%evalence(a)
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h2,h3,h1,p1,jt,ans1)
                 CALL pphhmtx(a,p1,h2,h3,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_31(l)=effoper_diagram_31(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_31
!
!
!
SUBROUTINE effective_operator_32(a,c,effoper_diagram_32)
  USE bare_operator_value
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a, c
  INTEGER :: p1, h2, l, h3, h1,  nshell1, nshell2, &
       idiff1, idiff2, nshell3, j1min, j2min, &
       j3min, j1max, j2max, j3max, jt, jmin_local, jmax_local, iph
  REAL(DP) :: val1, val2, factr,  fact, den1, den2
  REAL(DP), DIMENSION(n_startenergy_veff ) :: de, w1, w2
  REAL(DP), DIMENSION(n_startenergy_veff ), INTENT(OUT) :: effoper_diagram_32
  REAL(DP), DIMENSION(n_startenergy_g) :: ans1, ans2
  LOGICAL dencheck

  effoper_diagram_32=0.
  DO h1=1, all_orbit%total_orbits
     IF (all_orbit%orbit_status(h1) /= 'hole') CYCLE
     IF((all_orbit%jj(c) /= all_orbit%jj(h1)).AND. &
          (all_orbit%ll(c) /= all_orbit%ll(h1))) CYCLE
     IF(bare_operator(a,h1) == 0.) CYCLE
     fact=bare_operator(a,h1)*0.5/(all_orbit%jj(c)+1.)* &
          iph((all_orbit%jj(c)-all_orbit%jj(h1))/2)
     DO h2=1,all_orbit%total_orbits
        IF (all_orbit%orbit_status(h2) /= 'hole') CYCLE
        DO h3=1, all_orbit%total_orbits
           IF (all_orbit%orbit_status(h3) /= 'hole') CYCLE
           DO p1=1,all_orbit%total_orbits
              IF (all_orbit%orbit_status(p1) /= 'particle') CYCLE
              nshell1=all_orbit%nshell(h1)-all_orbit%nshell(c)
              nshell2=all_orbit%nshell(h2)+all_orbit%nshell(h3)
              nshell3=all_orbit%nshell(p1)+all_orbit%nshell(c)
              idiff1=nshell1
              idiff2=nshell2-nshell3
              IF(dencheck(idiff1)) CYCLE
              IF(dencheck(idiff2)) CYCLE
              den1=-(all_orbit%e(p1)+all_orbit%evalence(c)-all_orbit%e(h3)- &
                   all_orbit%e(h2))
              den2=all_orbit%e(h1)-all_orbit%evalence(c)
              j1min=ABS(all_orbit%jj(c)-all_orbit%jj(p1))/2
              j2min=ABS(all_orbit%jj(h2)-all_orbit%jj(h3))/2
              j3min=ABS(all_orbit%jj(h1)-all_orbit%jj(p1))/2
              j1max=(all_orbit%jj(h1)+all_orbit%jj(p1))/2
              j2max=(all_orbit%jj(h2)+all_orbit%jj(h3))/2
              j3max=(all_orbit%jj(c)+all_orbit%jj(p1))/2
              jmin_local=MAX(j1min,j2min,j3min)
              jmax_local=MIN(j1max,j2max,j3max)
              IF(jmin_local > jmax_local) CYCLE
              w1=all_orbit%e(h2)+all_orbit%e(h3)+wcn+ &
                   all_orbit%e(h1)-all_orbit%evalence(c)
              w2=all_orbit%e(h2)+all_orbit%e(h3)+wcn
              de= (den1 + wcn)*(den2 + wcn)
              DO jt = jmin_local, jmax_local
                 CALL pphhmtx(h1,p1,h2,h3,jt,ans1)
                 CALL pphhmtx(h2,h3,c,p1,jt,ans2)
                 IF( (ans1(1) == 0.0_dp).OR.(ans2(1) == 0.0_dp) ) CYCLE
                 factr=fact*(2.*jt+1.)
                 DO l=1,n_startenergy_veff
                    CALL interpolate(w1(l),e_start_g,ans1,val1)
                    CALL interpolate(w2(l),e_start_g,ans2,val2)
                    effoper_diagram_32(l)=effoper_diagram_32(l)+ &
                         factr*val1*val2/de(l)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE effective_operator_32

