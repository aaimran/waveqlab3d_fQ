module initial_stress_condition

  use common, only : wp

  implicit none

contains
  subroutine initial_stress_tensor(s0,x,problem,Mat)

    ! s0 = (/sxx0,syy0,szz0,sxy0,sxz0,syz0 /)
    !           1    2    3    4    5    6
    
    implicit none

    real(kind = wp), intent(in) :: x(:)
    real(kind = wp), optional, intent(in) :: Mat(3)
    real(kind = wp) :: Rx,Ry,Rt
    real(kind = wp), intent(out) :: s0(:)
    character(*),intent(in) :: problem
    real(kind = wp) :: Omega, b11, b33, b13, sigma22, Pf, u0, psi
    real(kind = wp),parameter::pi=3.14159265359_wp
    

    select case(problem)

    case('TPV26')
       if (x(2) .le. 15d0) then
          Omega = 1d0
       elseif (x(2) .ge. 15d0 .and. x(2) .le. 20d0) then
          Omega = (20d0 - x(2))/5d0
       elseif (x(2) .ge. 20d0) then

          Omega = 0d0

       end if

       b11 = 0.926793d0
       b33 = 1.073206d0
       b13 = -0.169029d0

       Pf = 9.8d0*x(2)

       sigma22 = -2.67d0*9.8d0*x(2)

       s0(1) = Omega*(b33*(sigma22 + Pf) - Pf) + (1d0-Omega)*sigma22 + Pf
       s0(2) = sigma22 + Pf
       s0(3) = Omega*(b11*(sigma22 + Pf) - Pf) + (1d0-Omega)*sigma22 + Pf
       s0(4) = 0d0
       s0(5) = Omega*(b13*(sigma22 + Pf))
       s0(6) = 0d0

    case('TPV27')
       if (x(2) .le. 15d0) then
          Omega = 1d0
       elseif (x(2) .ge. 15d0 .and. x(2) .le. 20d0) then
          Omega = (20d0 - x(2))/5d0
       elseif (x(2) .ge. 20d0) then

          Omega = 0d0

       end if

       b11 = 0.926793d0
       b33 = 1.073206d0
       b13 = -0.169029d0

       Pf = 9.8d0*x(2)

       sigma22 = -2.67d0*9.8d0*x(2)

       s0(1) = Omega*(b33*(sigma22 + Pf) - Pf) + (1-Omega)*sigma22 + Pf
       s0(2) = sigma22 + Pf
       s0(3) = Omega*(b11*(sigma22 + Pf) - Pf) + (1-Omega)*sigma22 + Pf
       s0(4) = 0d0
       s0(5) = Omega*(b13*(sigma22 + Pf))
       s0(6) = 0d0

    case('TPV28')
       s0(1) = -60d0
       s0(2) = 0d0
       s0(3) = 60d0
       s0(4) = 0d0
       s0(5) = 29.38d0
       s0(6) = 0d0


    case('TPV29','TPV30')

       if (x(2) .le. 17d0) then
          Omega = 1d0
       elseif (x(2) .ge. 17d0 .and. x(2) .le. 22d0) then
          Omega = (22d0 - x(2))/5d0
       elseif (x(2) .ge. 22d0) then
          Omega = 0d0
       end if

       b11 = 1.025837d0
       b33 = 0.974162d0
       b13 =-0.158649d0

       Pf = 9.8d0*x(2)

       sigma22 = -2.67d0*9.8d0*x(2)

       s0(1) = Omega*(b33*(sigma22 + Pf) - Pf) + (1d0-Omega)*sigma22 + Pf
       s0(2) = sigma22 + Pf
       s0(3) = Omega*(b11*(sigma22 + Pf) - Pf) + (1d0-Omega)*sigma22 + Pf
       s0(4) = 0d0
       s0(5) = Omega*(b13*(sigma22 + Pf))
       s0(6) = 0d0



    case('TPV101', 'TPV102')
       s0(1) = -120d0
       s0(2) = 0d0
       s0(3) = -120d0
       s0(4) = 0d0
       s0(5) = 75d0
       s0(6) = 0d0

    case('TPV31','TPV32')

       u0 = 32.03812032_wp

       s0(1) = -60_wp*Mat(2)/u0
       s0(2) = 0_wp
       s0(3) = -60_wp*Mat(2)/u0
       s0(4) = 0_wp
       s0(5) = 30_wp*Mat(2)/u0
       s0(6) = 0_wp  

    case('TPV33')

       if(x(3) < -9.8_wp) then
           Rx = (-x(3)-9.8_wp)/10.0_wp
       elseif(x(3) > 1.1_wp) then
           Rx = (x(3)-1.1_wp)/10.0_wp
       else
           Rx = 0.0_wp
       end if

       if(x(2) < 2.3_wp) then
           Ry = (-x(2)+2.3_wp)/10.0_wp
       elseif(x(2) > 8.0_wp) then
           Ry = (x(2)-8.0_wp)/10.0_wp
       else
           Ry = 0.0_wp
       end if 

       Rt = min(1.0_wp,sqrt(Rx**2+Ry**2))

       s0(1) = -60.0_wp
       s0(2) = 0.0_wp
       s0(3) = 0.0_wp
       s0(4) = 0.0_wp
       s0(5) = 30.0_wp*(1_wp-Rt)
       s0(6) = 0.0_wp 

    case('TPV34')

       ! Bring in material properties for point, use shear modulus
       u0 = 32.03812032_wp

       s0(1) = -60_wp*Mat(2)/u0
       s0(2) = 0_wp
       s0(3) = -60_wp*Mat(2)/u0
       s0(4) = 0_wp
       s0(5) = 30_wp*Mat(2)/u0
       s0(6) = 0_wp  

    case('SCITS2016')

       ! Bring in material properties for point, use shear modulus
       u0 = 32.03812032_wp

       s0(1) = -60_wp*Mat(2)/u0
       s0(2) = 0_wp
       s0(3) = -60_wp*Mat(2)/u0
       s0(4) = 0_wp
       s0(5) = 30_wp*Mat(2)/u0
       s0(6) = 0_wp                        
       
    case('rate-weakening')
       
       psi = 50.0_wp/180.0_wp*pi
       s0(1) = -126.0_wp
       s0(5) = 42.0_wp
       s0(3) = (1.0_wp - 2.0_wp*s0(5)/(s0(1)*tan(2.0_wp*psi)))*s0(1)
       s0(2) = 0.5_wp*(s0(1) + s0(3))
       s0(4) = 0.0_wp
       s0(6) = 0.0_wp


    case('mg_a1')
       
       s0(1) = -60.0_wp
       s0(2) = 0.0_wp
       s0(3) = -60.0_wp
       s0(4) = 0.0_wp
       s0(5) = 29.38_wp
       s0(6) = 0_wp
       
    case('mg_b1a')
       
       s0(1) = -60.0_wp
       s0(2) = 0.0_wp
       s0(3) = -60.0_wp
       s0(4) = 0.0_wp
       s0(5) = 31.0_wp
       s0(6) = 0_wp
       
    case('TPV10')
       
       s0(1) = -7.378*x(2)
       s0(2) = 0.0_wp 
       s0(3) = - 7.378*x(2)
       s0(4) = 0.0_wp
       s0(5) = -0.55_wp*s0(1)
       s0(6) = 0d0

    case('wasatch100')
                                                                                                                                                                                                                                              

       Pf = 9.8d0*(x(2)+2.149d0)
       sigma22=-2.7*9.8d0*(x(2)+2.149d0)

       !print *, x(2)
       
       !Pf = 9.8d0*(x(2)+0d0*2.781d0)
       !sigma22=-2.7*9.8d0*(x(2)+0d0*2.781d0)
!       Omega=3.3*x(2)                                                                                                                                                                                  

       s0(1) = sigma22+Pf-0.4d0*sigma22
       s0(2) = sigma22+Pf
!       s0(3)=sigma22+Pf+0.25d0*(-0.4d0*sigma22)                                                                                                                                                        
       s0(3)=(s0(1)+s0(2))*0.5d0
!       s0(1) = sigma22+Pf+Omega                                                                                                                                                                        
!       s0(2) = sigma22+Pf                                                                                                                                                                              
!       s0(3)=sigma22+Pf+0.25d0*Omega                                                                                                                                                                   
       s0(4) = 0d0
       s0(5) = 0d0
       s0(6) = 0d0


       
    
    case default

       stop 'invalid problem in intial stress tensor'

    end select

  end subroutine initial_stress_tensor

end module initial_stress_condition
