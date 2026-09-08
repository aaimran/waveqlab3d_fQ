!> This module contains routines that compute the metrics
!> of the grid using finite differences.

module metrics

  use common, only : wp
  implicit none

contains

  ! These routines are taken from k.0_wpd. They were subroutines
  ! originally. They were made into pure functions so that the
  ! three do loops in the grid.f90 file coulde be replaced by
  ! an equivalent concurrent do loop.

  pure function D6x(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mx)

    implicit none

    integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mx
    real(kind = wp),intent(in) :: hx        ! Spatial step
    real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,1:3), intent(in) :: u
    integer,intent(in) :: x, y, z

    integer :: n
    real(kind = wp), dimension(3) :: u_x, D6x

    n = size(u_x)

    if (x .eq. 1) then
       u_x(1:n) =  0.036577936277544_wp*u(6, y, z, 1:n) &
            & + 0.104488069284042_wp*u(5, y, z, 1:n) &
            & - 0.450398306578272_wp*u(4, y, z, 1:n) &
            & - 0.141512858744873_wp*u(3, y, z, 1:n) &
            & + 2.033378678700676_wp*u(2, y, z, 1:n) &
            & - 1.582533518939116_wp*u(1, y, z, 1:n)
    end if

    if (x .eq. 2) then
       u_x(1:n) = -0.014903449191300_wp*u(6, y, z, 1:n) &
            & - 0.069112065532624_wp*u(5, y, z, 1:n) &
            & + 0.258816087376832_wp*u(4, y, z, 1:n) &
            & + 0.287258622978251_wp*u(3, y, z, 1:n) &
            & - 0.462059195631158_wp*u(1, y, z, 1:n)

    end if
    if (x .eq. 3) then
       u_x(1:n) = -0.018129342917256_wp*u(6, y, z, 1:n) &
            & - 0.022902190275815_wp*u(5, y, z, 1:n) &
            & + 0.606235523609147_wp*u(4, y, z, 1:n) &
            & - 0.636451095137907_wp*u(2, y, z, 1:n) &
            & + 0.071247104721830_wp*u(1, y, z, 1:n)
    end if

    if (x .eq. 4) then
       u_x(1:n) =  0.013435342414630*u(7, y, z, 1:n) &
            & - 0.051642265516119_wp*u(6, y, z, 1:n) &
            & + 0.520262285050482_wp*u(5, y, z, 1:n) &
            & - 0.306681191361148_wp*u(3, y, z, 1:n) &
            & - 0.290087484386815_wp*u(2, y, z, 1:n) &
            & + 0.114713313798970_wp*u(1, y, z, 1:n)
    end if



    if (x .eq. 5) then
       u_x(1:n) =  0.018281071473911_wp*u(8, y, z, 1:n) &
            & - 0.164529643265203_wp*u(7, y, z, 1:n) &
            & + 0.769199413962647_wp*u(6, y, z, 1:n) &
            & - 0.707905442575989_wp*u(4, y, z, 1:n) &
            & + 0.015764336127392_wp*u(3, y, z, 1:n) &
            & + 0.105400944933782_wp*u(2, y, z, 1:n) &
            & - 0.036210680656541_wp*u(1, y, z, 1:n)
    end if

    if (x .eq. 6) then
       u_x(1:n) = 0.016437980868017_wp*u(9, y, z, 1:n) &
            & - 0.147941827812150_wp*u(8, y, z, 1:n) &
            & + 0.739709139060752_wp*u(7, y, z, 1:n) &
            & - 0.691649024426814_wp*u(5, y, z, 1:n) &
            & + 0.063183694641876_wp*u(4, y, z, 1:n) &
            & + 0.011220896474665_wp*u(3, y, z, 1:n) &
            & + 0.020437334208704_wp*u(2, y, z, 1:n) &
            & - 0.011398193015050_wp*u(1, y, z, 1:n)
    end if


    if (x .eq. mx-5) then
       u_x(1:n) = - 0.016437980868017_wp*u(mx-8, y, z, 1:n) &
            & + 0.147941827812150_wp*u(mx-7, y, z, 1:n) &
            & - 0.739709139060752_wp*u(mx-6, y, z, 1:n) &
            & + 0.691649024426814_wp*u(mx-4, y, z, 1:n) &
            & - 0.063183694641876_wp*u(mx-3, y, z, 1:n) &
            & - 0.011220896474665_wp*u(mx-2, y, z, 1:n) &
            & - 0.020437334208704_wp*u(mx-1, y, z, 1:n) &
            & + 0.011398193015050_wp*u(mx, y, z, 1:n)
    end if


    if (x .eq. mx-4) then
       u_x(1:n) = - 0.018281071473911_wp*u(mx-7, y, z, 1:n) &
            & + 0.164529643265203_wp*u(mx-6, y, z, 1:n) &
            & - 0.769199413962647_wp*u(mx-5, y, z, 1:n) &
            & + 0.707905442575989_wp*u(mx-3, y, z, 1:n) &
            & - 0.015764336127392_wp*u(mx-2, y, z, 1:n) &
            & - 0.105400944933782_wp*u(mx-1, y, z, 1:n) &
            & + 0.036210680656541_wp*u(mx, y, z, 1:n)
    end if


    if (x .eq. mx-3) then
       u_x(1:n) = - 0.013435342414630*u(mx-6, y, z, 1:n) &
            & + 0.051642265516119_wp*u(mx-5, y, z, 1:n) &
            & - 0.520262285050482_wp*u(mx-4, y, z, 1:n) &
            & + 0.306681191361148_wp*u(mx-2, y, z, 1:n) &
            & + 0.290087484386815_wp*u(mx-1, y, z, 1:n) &
            & - 0.114713313798970_wp*u(mx, y, z, 1:n)
    end if


    if (x .eq. mx-2) then
       u_x(1:n) = 0.018129342917256_wp*u(mx-5, y, z, 1:n) &
            & + 0.022902190275815_wp*u(mx-4, y, z, 1:n) &
            & - 0.606235523609147_wp*u(mx-3, y, z, 1:n) &
            & + 0.636451095137907_wp*u(mx-1, y, z, 1:n) &
            & - 0.071247104721830_wp*u(mx, y, z, 1:n)
    end if



    if (x .eq. mx-1) then
       u_x(1:n) = 0.014903449191300_wp*u(mx-5, y, z, 1:n) &
            & + 0.069112065532624_wp*u(mx-4, y, z, 1:n) &
            & - 0.258816087376832_wp*u(mx-3, y, z, 1:n) &
            & - 0.287258622978251_wp*u(mx-2, y, z, 1:n) &
            & + 0.462059195631158_wp*u(mx, y, z, 1:n)

    end if



    if (x .eq. mx) then
       u_x(1:n) =   - 0.036577936277544_wp*u(mx-5, y, z, 1:n) &
            & - 0.104488069284042_wp*u(mx-4, y, z, 1:n) &
            & + 0.450398306578272_wp*u(mx-3, y, z, 1:n) &
            & + 0.141512858744873_wp*u(mx-2, y, z, 1:n) &
            & - 2.033378678700676_wp*u(mx-1, y, z, 1:n) &
            & + 1.582533518939116_wp*u(mx, y, z, 1:n)
    end if



    if ((x .ge. 7) .and. (x .le. mx-6)) then

       u_x(1:n) = -1.0_wp/60.0_wp*u(x-3, y, z, 1:n) + 0.15_wp*u(x-2, y, z, 1:n) &
            & - 0.75_wp*u(x-1, y, z, 1:n) + 0.75_wp*u(x+1, y, z, 1:n) - 0.15_wp*u(x+2, y, z, 1:n) &
            & + 1.0_wp/60.0_wp*u(x+3, y, z, 1:n)

       !u_x(1:n) = 0.016666666666667_wp*u(x-4, y, z, 1:n) - 0.133333333333333_wp*u(x-3, y, z, 1:n) + 0.5_wp*u(x-2, y, z, 1:n) &
       !     & - 1.333333333333333_wp*u(x-1, y, z, 1:n) + 0.583333333333333_wp*u(x, y, z, 1:n) + 0.4_wp*u(x+1, y, z, 1:n)&
       !     & - 0.033333333333333_wp*u(x+2, y, z, 1:n)
                   
    end if


    u_x = (1.0_wp/hx)*u_x

    D6x = u_x

  end function D6x

  pure function D6y(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, my)

    implicit none

    integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, my
    real(kind = wp),intent(in) :: hx        ! Spatial step
    real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,1:3), intent(in) :: u
    integer,intent(in) :: x, y, z

    integer :: n
    real(kind = wp), dimension(3) :: u_x, D6y

    n = 3

    if (y .eq. 1) then

       u_x(1:n) =  0.036577936277544_wp*u(x, 6, z, 1:n) &
            & + 0.104488069284042_wp*u(x, 5, z, 1:n) &
            & - 0.450398306578272_wp*u(x, 4, z, 1:n) &
            & - 0.141512858744873_wp*u(x, 3, z, 1:n) &
            & + 2.033378678700676_wp*u(x, 2, z, 1:n) &
            & - 1.582533518939116_wp*u(x, 1, z, 1:n)
    end if

    if (y .eq. 2) then
       u_x(1:n) = -0.014903449191300_wp*u(x, 6, z, 1:n) &
            & - 0.069112065532624_wp*u(x, 5, z, 1:n) &
            & + 0.258816087376832_wp*u(x, 4, z, 1:n) &
            & + 0.287258622978251_wp*u(x, 3, z, 1:n) &
            & - 0.462059195631158_wp*u(x, 1, z, 1:n)

    end if
    if (y .eq. 3) then
       u_x(1:n) = -0.018129342917256_wp*u(x, 6, z, 1:n) &
            & - 0.022902190275815_wp*u(x, 5, z, 1:n) &
            & + 0.606235523609147_wp*u(x, 4, z, 1:n) &
            & - 0.636451095137907_wp*u(x, 2, z, 1:n) &
            & + 0.071247104721830_wp*u(x, 1, z, 1:n)
    end if

    if (y .eq. 4) then
       u_x(1:n) =  0.013435342414630*u(x, 7, z, 1:n) &
            & - 0.051642265516119_wp*u(x, 6, z, 1:n) &
            & + 0.520262285050482_wp*u(x, 5, z, 1:n) &
            & - 0.306681191361148_wp*u(x, 3, z, 1:n) &
            & - 0.290087484386815_wp*u(x, 2, z, 1:n) &
            & + 0.114713313798970_wp*u(x, 1, z, 1:n)
    end if



    if (y .eq. 5) then
       u_x(1:n) =  0.018281071473911_wp*u(x, 8, z, 1:n) &
            & - 0.164529643265203_wp*u(x, 7, z, 1:n) &
            & + 0.769199413962647_wp*u(x, 6, z, 1:n) &
            & - 0.707905442575989_wp*u(x, 4, z, 1:n) &
            & + 0.015764336127392_wp*u(x, 3, z, 1:n) &
            & + 0.105400944933782_wp*u(x, 2, z, 1:n) &
            & - 0.036210680656541_wp*u(x, 1, z, 1:n)
    end if

    if (y .eq. 6) then
       u_x(1:n) = 0.016437980868017_wp*u(x, 9, z, 1:n) &
            & - 0.147941827812150_wp*u(x, 8, z, 1:n) &
            & + 0.739709139060752_wp*u(x, 7, z, 1:n) &
            & - 0.691649024426814_wp*u(x, 5, z, 1:n) &
            & + 0.063183694641876_wp*u(x, 4, z, 1:n) &
            & + 0.011220896474665_wp*u(x, 3, z, 1:n) &
            & + 0.020437334208704_wp*u(x, 2, z, 1:n) &
            & - 0.011398193015050_wp*u(x, 1, z, 1:n)
    end if


    if (y .eq. my-5) then
       u_x(1:n) = - 0.016437980868017_wp*u(x, my-8, z, 1:n) &
            & + 0.147941827812150_wp*u(x, my-7, z, 1:n) &
            & - 0.739709139060752_wp*u(x, my-6, z, 1:n) &
            & + 0.691649024426814_wp*u(x, my-4, z, 1:n) &
            & - 0.063183694641876_wp*u(x, my-3, z, 1:n) &
            & - 0.011220896474665_wp*u(x, my-2, z, 1:n) &
            & - 0.020437334208704_wp*u(x, my-1, z, 1:n) &
            & + 0.011398193015050_wp*u(x, my, z, 1:n)
    end if


    if (y .eq. my-4) then
       u_x(1:n) = - 0.018281071473911_wp*u(x, my-7, z, 1:n) &
            & + 0.164529643265203_wp*u(x, my-6, z, 1:n) &
            & - 0.769199413962647_wp*u(x, my-5, z, 1:n) &
            & + 0.707905442575989_wp*u(x, my-3, z, 1:n) &
            & - 0.015764336127392_wp*u(x, my-2, z, 1:n) &
            & - 0.105400944933782_wp*u(x, my-1, z, 1:n) &
            & + 0.036210680656541_wp*u(x, my, z, 1:n)
    end if


    if (y .eq. my-3) then
       u_x(1:n) = - 0.013435342414630*u(x, my-6, z, 1:n) &
            & + 0.051642265516119_wp*u(x, my-5, z, 1:n) &
            & - 0.520262285050482_wp*u(x, my-4, z, 1:n) &
            & + 0.306681191361148_wp*u(x, my-2, z, 1:n) &
            & + 0.290087484386815_wp*u(x, my-1, z, 1:n) &
            & - 0.114713313798970_wp*u(x, my, z, 1:n)
    end if


    if (y .eq. my-2) then
       u_x(1:n) = 0.018129342917256_wp*u(x, my-5, z, 1:n) &
            & + 0.022902190275815_wp*u(x, my-4, z, 1:n) &
            & - 0.606235523609147_wp*u(x, my-3, z, 1:n) &
            & + 0.636451095137907_wp*u(x, my-1, z, 1:n) &
            & - 0.071247104721830_wp*u(x, my, z, 1:n)
    end if



    if (y .eq. my-1) then
       u_x(1:n) = 0.014903449191300_wp*u(x, my-5, z, 1:n) &
            & + 0.069112065532624_wp*u(x, my-4, z, 1:n) &
            & - 0.258816087376832_wp*u(x, my-3, z, 1:n) &
            & - 0.287258622978251_wp*u(x, my-2, z, 1:n) &
            & + 0.462059195631158_wp*u(x, my, z, 1:n)

    end if



    if (y .eq. my) then
       u_x(1:n) =   - 0.036577936277544_wp*u(x, my-5, z, 1:n) &
            & - 0.104488069284042_wp*u(x, my-4, z, 1:n) &
            & + 0.450398306578272_wp*u(x, my-3, z, 1:n) &
            & + 0.141512858744873_wp*u(x, my-2, z, 1:n) &
            & - 2.033378678700676_wp*u(x, my-1, z, 1:n) &
            & + 1.582533518939116_wp*u(x, my, z, 1:n)
    end if



    if ((y .ge. 7) .and. (y .le. my-6)) then

       u_x(1:n) = -1.0_wp/60.0_wp*u(x, y-3, z, 1:n) + 0.15_wp*u(x, y-2, z, 1:n) &
            & - 0.75_wp*u(x, y-1, z, 1:n) + 0.75_wp*u(x, y+1, z, 1:n) - 0.15_wp*u(x, y+2, z, 1:n) &
            & + 1.0_wp/60.0_wp*u(x, y+3, z, 1:n)

       !u_x(1:n) = 0.016666666666667_wp*u(x, y-4, z, 1:n) - 0.133333333333333_wp*u(x, y-3, z, 1:n) + 0.5_wp*u(x, y-2, z, 1:n) &
       !     & - 1.333333333333333_wp*u(x, y-1, z, 1:n) + 0.583333333333333_wp*u(x, y, z, 1:n) + 0.4_wp*u(x, y+1, z, 1:n)&
       !     & - 0.033333333333333_wp*u(x, y+2, z, 1:n)
    end if



    u_x = (1.0_wp/hx)*u_x

    D6y = u_x

  end function D6y


  pure function D6z(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mz)

    implicit none

    integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mz
    real(kind = wp),intent(in) :: hx        ! Spatial step
    real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,1:3), intent(in) :: u
    integer,intent(in) :: x, y, z

    integer :: n
    real(kind = wp), dimension(3) :: u_x, D6z

    n = 3

    if (z .eq. 1) then
       u_x(1:n) =  0.036577936277544_wp*u(x, y, 6, 1:n) &
            & + 0.104488069284042_wp*u(x, y, 5, 1:n) &
            & - 0.450398306578272_wp*u(x, y, 4, 1:n) &
            & - 0.141512858744873_wp*u(x, y, 3, 1:n) &
            & + 2.033378678700676_wp*u(x, y, 2, 1:n) &
            & - 1.582533518939116_wp*u(x, y, 1, 1:n)
    end if

    if (z .eq. 2) then
       u_x(1:n) = -0.014903449191300_wp*u(x, y, 6, 1:n) &
            & - 0.069112065532624_wp*u(x, y, 5, 1:n) &
            & + 0.258816087376832_wp*u(x, y, 4, 1:n) &
            & + 0.287258622978251_wp*u(x, y, 3, 1:n) &
            & - 0.462059195631158_wp*u(x, y, 1, 1:n)

    end if
    if (z .eq. 3) then
       u_x(1:n) = -0.018129342917256_wp*u(x, y, 6, 1:n) &
            & - 0.022902190275815_wp*u(x, y, 5, 1:n) &
            & + 0.606235523609147_wp*u(x, y, 4, 1:n) &
            & - 0.636451095137907_wp*u(x, y, 2, 1:n) &
            & + 0.071247104721830_wp*u(x, y, 1, 1:n)
    end if

    if (z .eq. 4) then
       u_x(1:n) =  0.013435342414630*u(x, y, 7, 1:n) &
            & - 0.051642265516119_wp*u(x, y, 6, 1:n) &
            & + 0.520262285050482_wp*u(x, y, 5, 1:n) &
            & - 0.306681191361148_wp*u(x, y, 3, 1:n) &
            & - 0.290087484386815_wp*u(x, y, 2, 1:n) &
            & + 0.114713313798970_wp*u(x, y, 1, 1:n)
    end if



    if (z .eq. 5) then
       u_x(1:n) =  0.018281071473911_wp*u(x, y, 8, 1:n) &
            & - 0.164529643265203_wp*u(x, y, 7, 1:n) &
            & + 0.769199413962647_wp*u(x, y, 6, 1:n) &
            & - 0.707905442575989_wp*u(x, y, 4, 1:n) &
            & + 0.015764336127392_wp*u(x, y, 3, 1:n) &
            & + 0.105400944933782_wp*u(x, y, 2, 1:n) &
            & - 0.036210680656541_wp*u(x, y, 1, 1:n)
    end if

    if (z .eq. 6) then
       u_x(1:n) = 0.016437980868017_wp*u(x, y, 9, 1:n) &
            & - 0.147941827812150_wp*u(x, y, 8, 1:n) &
            & + 0.739709139060752_wp*u(x, y, 7, 1:n) &
            & - 0.691649024426814_wp*u(x, y, 5, 1:n) &
            & + 0.063183694641876_wp*u(x, y, 4, 1:n) &
            & + 0.011220896474665_wp*u(x, y, 3, 1:n) &
            & + 0.020437334208704_wp*u(x, y, 2, 1:n) &
            & - 0.011398193015050_wp*u(x, y, 1, 1:n)
    end if


    if (z .eq. mz-5) then
       u_x(1:n) = - 0.016437980868017_wp*u(x, y, mz-8, 1:n) &
            & + 0.147941827812150_wp*u(x, y, mz-7, 1:n) &
            & - 0.739709139060752_wp*u(x, y, mz-6, 1:n) &
            & + 0.691649024426814_wp*u(x, y, mz-4, 1:n) &
            & - 0.063183694641876_wp*u(x, y, mz-3, 1:n) &
            & - 0.011220896474665_wp*u(x, y, mz-2, 1:n) &
            & - 0.020437334208704_wp*u(x, y, mz-1, 1:n) &
            & + 0.011398193015050_wp*u(x, y, mz, 1:n)
    end if


    if (z .eq. mz-4) then
       u_x(1:n) = -0.018281071473911_wp*u(x, y, mz-7, 1:n) &
            & + 0.164529643265203_wp*u(x, y, mz-6, 1:n) &
            & - 0.769199413962647_wp*u(x, y, mz-5, 1:n) &
            & + 0.707905442575989_wp*u(x, y, mz-3, 1:n) &
            & - 0.015764336127392_wp*u(x, y, mz-2, 1:n) &
            & - 0.105400944933782_wp*u(x, y, mz-1, 1:n) &
            & + 0.036210680656541_wp*u(x, y, mz, 1:n)



    end if


    if (z .eq. mz-3) then
       u_x(1:n) = - 0.013435342414630*u(x, y, mz-6, 1:n) &
            & + 0.051642265516119_wp*u(x, y, mz-5, 1:n) &
            & - 0.520262285050482_wp*u(x, y, mz-4, 1:n) &
            & + 0.306681191361148_wp*u(x, y, mz-2, 1:n) &
            & + 0.290087484386815_wp*u(x, y, mz-1, 1:n) &
            & - 0.114713313798970_wp*u(x, y, mz, 1:n)
    end if


    if (z .eq. mz-2) then
       u_x(1:n) = 0.018129342917256_wp*u(x, y, mz-5, 1:n) &
            & + 0.022902190275815_wp*u(x, y, mz-4, 1:n) &
            & - 0.606235523609147_wp*u(x, y, mz-3, 1:n) &
            & + 0.636451095137907_wp*u(x, y, mz-1, 1:n) &
            & - 0.071247104721830_wp*u(x, y, mz, 1:n)
    end if



    if (z .eq. mz-1) then
       u_x(1:n) = 0.014903449191300_wp*u(x, y, mz-5, 1:n) &
            & + 0.069112065532624_wp*u(x, y, mz-4, 1:n) &
            & - 0.258816087376832_wp*u(x, y, mz-3, 1:n) &
            & - 0.287258622978251_wp*u(x, y, mz-2, 1:n) &
            & + 0.462059195631158_wp*u(x, y, mz, 1:n)

    end if



    if (z .eq. mz) then
       u_x(1:n) =   - 0.036577936277544_wp*u(x, y, mz-5, 1:n) &
            & - 0.104488069284042_wp*u(x, y, mz-4, 1:n) &
            & + 0.450398306578272_wp*u(x, y, mz-3, 1:n) &
            & + 0.141512858744873_wp*u(x, y, mz-2, 1:n) &
            & - 2.033378678700676_wp*u(x, y, mz-1, 1:n) &
            & + 1.582533518939116_wp*u(x, y, mz, 1:n)
    end if



    if ((z .ge. 7) .and. (z .le. mz-6)) then
       u_x(1:n) = -1.0_wp/60.0_wp*u(x, y, z-3, 1:n) + 0.15_wp*u(x, y, z-2, 1:n) &
            & - 0.75_wp*u(x, y, z-1, 1:n) + 0.75_wp*u(x, y, z+1, 1:n) - 0.15_wp*u(x, y, z+2, 1:n) &
            & + 1.0_wp/60.0_wp*u(x, y, z+3, 1:n)

       !u_x(1:n) = 0.016666666666667_wp*u(x, y, z-4, 1:n) - 0.133333333333333_wp*u(x, y, z-3, 1:n) + 0.5_wp*u(x, y, z-2, 1:n) &
       !     & - 1.333333333333333_wp*u(x, y, z-1, 1:n) + 0.583333333333333_wp*u(x, y, z, 1:n) + 0.4_wp*u(x, y, z+1, 1:n)&
       !     & - 0.033333333333333_wp*u(x, y, z+2, 1:n)
    end if





    u_x = (1.0_wp/hx)*u_x

    D6z = u_x

  end function D6z



pure function Dfx6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mx, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mx
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dfx6

  n = size(u_x)

  if (x .eq. 1) then
    u_x(m:p) = -1.593253459452621_wp*u(x, y, z, m:p) + 2.093702214477299*u(x+1, y, z, m:p) &
               - 0.275607596716317_wp*u(x+2, y, z, m:p) - 0.302855902188632_wp*u(x+3, y, z, m:p) &
               + 0.023993033880123_wp*u(x+4, y, z, m:p) + 0.054021710000147_wp*u(x+5, y, z, m:p)
  end if

  if (x .eq. 2) then
    u_x(m:p) = -0.460835552271137_wp*u(x-1, y, z, m:p) - 0.029551468764759_wp*u(x, y, z, m:p)&
                   + 0.393228064437077_wp*u(x+1, y, z, m:p) + 0.105980141988699_wp*u(x+2, y, z, m:p)&
                   + 0.030739159126096_wp*u(x+3, y, z, m:p) - 0.039560344515975_wp*u(x+4, y, z, m:p)


                   !- 0.460835552271137*u_old[0] - 0.029551468764759*u_old[1]\
                   !+ 0.393228064437077*u_old[2] + 0.105980141988699*u_old[3]\
                   !+ 0.030739159126096*u_old[4] - 0.039560344515975*u_old[5]
                   

  end if
  if (x .eq. 3) then
     u_x(m:p) = 0.096333139747688_wp*u(x-2, y, z, m:p) - 0.629259919521358_wp*u(x-1, y, z, m:p)&
                  - 0.306540940575249_wp*u(x, y, z, m:p) + 1.258766829227482_wp*u(x+1, y, z, m:p)&
                  - 0.589452745232692_wp*u(x+2, y, z, m:p) + 0.197069524204597_wp*u(x+3, y, z, m:p)&
                  - 0.026915887850467_wp*u(x+4, y, z, m:p)
  end if

  if (x .eq. 4) then

     u_x(m:p) = 0.092596590361100_wp*u(x-3, y, z, m:p) - 0.240558007429252_wp*u(x-2, y, z, m:p)&
                  - 0.243684445410829_wp*u(x-1, y, z, m:p) - 0.318292975073623_wp*u(x, y, z, m:p)&
                  + 0.946962507015985_wp*u(x+1, y, z, m:p) - 0.330443502642251_wp*u(x+2, y, z, m:p)&
                  + 0.106765523632994_wp*u(x+3, y, z, m:p) - 0.013345690454124_wp*u(x+4, y, z, m:p)
                  
  end if



  if (x .eq. 5) then
    u_x(m:p) =  -0.015609194255155_wp*u(x-4, y, z, m:p) + 0.033485250803205_wp*u(x-3, y, z, m:p)&
                   + 0.041897570295664_wp*u(x-2, y, z, m:p) - 0.450132219590077_wp*u(x-1, y, z, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.406282074938901_wp*u(x+1, y, z, m:p)&
                   - 0.550950133911491_wp*u(x+2, y, z, m:p) + 0.146920035709731_wp*u(x+3, y, z, m:p)&
                   - 0.018365004463716_wp*u(x+4, y, z, m:p)
  end if

  if (x .eq. 6) then
     
    u_x(m:p) = - 0.016163101657400_wp*u(x-5, y, z, m:p) + 0.040531470682718_wp*u(x-4, y, z, m:p)&
                   - 0.005338567444210_wp*u(x-3, y, z, m:p) - 0.025427590054054_wp*u(x-2, y, z, m:p)&
                   - 0.357460227911152_wp*u(x-1, y, z, m:p) - 0.572337434946968_wp*u(x, y, z, m:p)&
                   + 1.313958528183954_wp*u(x+1, y, z, m:p) - 0.492734448068983_wp*u(x+2, y, z, m:p)&
                   + 0.131395852818395_wp*u(x+3, y, z, m:p) - 0.016424481602299_wp*u(x+4, y, z, m:p)
  end if


  if (x .eq. mx-5) then
    u_x(m:p) = 0.032848963204599_wp*u(x-2, y, z, m:p) - 0.394187558455186_wp*u(x-1, y, z, m:p)&
                 - 0.572337434946968_wp*u(x, y, z, m:p) + 1.257688453771463_wp*u(x+1, y, z, m:p)&
                 - 0.406675342006739_wp*u(x+2, y, z, m:p) + 0.120254802392339_wp*u(x+3, y, z, m:p)&
                 - 0.054367610816546_wp*u(x+4, y, z, m:p) + 0.016775726857039_wp*u(x+5, y, z, m:p)
  end if


  if (x .eq. mx-4) then
     u_x(m:p) = 0.036730008927433_wp*u(x-2, y, z, m:p) - 0.399693508760138_wp*u(x-1, y, z, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.303115093827634_wp*u(x+1, y, z, m:p)&
                   - 0.402190050630647_wp*u(x+2, y, z, m:p) + 0.047235828122731_wp*u(x+3, y, z, m:p)&
                   + 0.008331008040049_wp*u(x+4, y, z, m:p)
  end if


  if (x .eq. mx-3) then
    u_x(m:p) = -0.020661154127888_wp*u(x-2, y, z, m:p) - 0.327107204245208_wp*u(x-1, y, z, m:p)&
                - 0.318292975073623_wp*u(x, y, z, m:p) + 0.624133691970994_wp*u(x+1, y, z, m:p)&
                + 0.118346128899150_wp*u(x+2, y, z, m:p) - 0.076418487423426_wp*u(x+3, y, z, m:p)
  end if


  if (x .eq. mx-2) then
    u_x(m:p) = - 0.008748664712219_wp*u(x-3, y, z, m:p) + 0.061405392278187_wp*u(x-2, y, z, m:p)&
                  - 0.491468255323896_wp*u(x-1, y, z, m:p) - 0.306540940575249_wp*u(x, y, z, m:p)&
                  + 0.885608401570530_wp*u(x+1, y, z, m:p) - 0.140255933237354_wp*u(x+2, y, z, m:p)
  end if



  if (x .eq. mx-1) then
     u_x(m:p) = 0.029492540133095_wp*u(x-4, y, z, m:p) + 0.021790841692085_wp*u(x-3, y, z, m:p)&
                  - 0.215422101432619_wp*u(x-2, y, z, m:p) - 0.279404147185598_wp*u(x-1, y, z, m:p)&
                  - 0.029551468764759_wp*u(x, y, z, m:p) + 0.473094335557797_wp*u(x+1, y, z, m:p)

  end if



  if (x .eq. mx) then
     u_x(m:p) = - 0.052048915548040_wp*u(x-5, y, z, m:p) - 0.044953974933765_wp*u(x-4, y, z, m:p)&
                + 0.366971721882124_wp*u(x-3, y, z, m:p) + 0.189297839436616_wp*u(x-2, y, z, m:p)&
                - 2.039450367044345_wp*u(x-1, y, z, m:p) + 1.580183696207410_wp*u(x, y, z, m:p)
  end if



  if ((x .ge. 7) .and. (x .le. mx-6)) then

     u_x(m:p) = 0.033333333333333_wp*u(x-2, y, z, m:p) - 0.4_wp*u(x-1, y, z, m:p) &
                     - 0.583333333333333_wp*u(x, y, z, m:p)   + 1.333333333333333_wp*u(x+1, y, z, m:p)&
                     - 0.5_wp*u(x+2, y, z, m:p) + 0.133333333333333_wp*u(x+3, y, z, m:p)&
                     - 0.016666666666667_wp*u(x+4, y, z, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dfx6 = u_x

end function Dfx6



pure function Dbx6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mx, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mx
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dbx6

  n = size(u_x)

  if (x .eq. 1) then
    u_x(m:p) = - 1.580183696207410_wp*u(x, y, z, m:p) + 2.039450367044345_wp*u(x+1, y, z, m:p) &
               - 0.189297839436616_wp*u(x+2, y, z, m:p) - 0.366971721882124_wp*u(x+3, y, z, m:p) &
               + 0.044953974933765_wp*u(x+4, y, z, m:p) + 0.052048915548040_wp*u(x+5, y, z, m:p)
  end if

  if (x .eq. 2) then
    u_x(m:p) = -0.473094335557797_wp*u(x-1, y, z, m:p) + 0.029551468764759_wp*u(x, y, z, m:p)&
               + 0.279404147185598_wp*u(x+1, y, z, m:p) + 0.215422101432619_wp*u(x+2, y, z, m:p)&
               - 0.021790841692085_wp*u(x+3, y, z, m:p) - 0.029492540133095_wp*u(x+4, y, z, m:p)


  end if
  if (x .eq. 3) then
     u_x(m:p) = 0.140255933237354_wp*u(x-2, y, z, m:p) - 0.885608401570530_wp*u(x-1, y, z, m:p) &
                + 0.306540940575249_wp*u(x, y, z, m:p) + 0.491468255323896_wp*u(x+1, y, z, m:p) &
                - 0.061405392278187_wp*u(x+2, y, z, m:p) + 0.008748664712219_wp*u(x+3, y, z, m:p)
  end if

  if (x .eq. 4) then

     u_x(m:p) = 0.076418487423426_wp*u(x-3, y, z, m:p) - 0.118346128899150_wp*u(x-2, y, z, m:p) &
                -0.624133691970994_wp*u(x-1, y, z, m:p) + 0.318292975073623_wp*u(x, y, z, m:p) &
                +0.327107204245208_wp*u(x+1, y, z, m:p) + 0.020661154127888_wp*u(x+2, y, z, m:p)
                  
  end if



  if (x .eq. 5) then
    u_x(m:p) =  - 0.008331008040049_wp*u(x-4, y, z, m:p) - 0.047235828122731_wp*u(x-3, y, z, m:p)&
                + 0.402190050630647_wp*u(x-2, y, z, m:p) - 1.303115093827634_wp*u(x-1, y, z, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.399693508760138_wp*u(x+1, y, z, m:p)&
                - 0.036730008927433_wp*u(x+2, y, z, m:p)
  end if

  if (x .eq. 6) then
     
    u_x(m:p) = - 0.016775726857039_wp*u(x-5, y, z, m:p) + 0.054367610816546_wp*u(x-4, y, z, m:p) &
               - 0.120254802392339_wp*u(x-3, y, z, m:p) + 0.406675342006739_wp*u(x-2, y, z, m:p)&
               - 1.257688453771463_wp*u(x-1, y, z, m:p) + 0.572337434946968_wp*u(x, y, z, m:p) &
               + 0.394187558455186_wp*u(x+1, y, z, m:p) - 0.032848963204599_wp*u(x+2, y, z, m:p)
  end if


  if (x .eq. mx-5) then
    u_x(m:p) = 0.016424481602299_wp*u(x-4, y, z, m:p) - 0.131395852818395_wp*u(x-3, y, z, m:p)&
               + 0.492734448068983_wp*u(x-2, y, z, m:p) - 1.313958528183954_wp*u(x-1, y, z, m:p)&
               + 0.572337434946968_wp*u(x, y, z, m:p) + 0.357460227911152_wp*u(x+1, y, z, m:p)&
               + 0.025427590054054_wp*u(x+2, y, z, m:p) + 0.005338567444210_wp*u(x+3, y, z, m:p)&
               - 0.040531470682718_wp*u(x+4, y, z, m:p) + 0.016163101657400_wp*u(x+5, y, z, m:p)
  end if


  if (x .eq. mx-4) then
     u_x(m:p) = 0.018365004463716_wp*u(x-4, y, z, m:p) - 0.146920035709731_wp*u(x-3, y, z, m:p)&
                + 0.550950133911491_wp*u(x-2, y, z, m:p) - 1.406282074938901_wp*u(x-1, y, z, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.450132219590077_wp*u(x+1, y, z, m:p)&
                - 0.041897570295664_wp*u(x+2, y, z, m:p) - 0.033485250803205_wp*u(x+3, y, z, m:p)&
                + 0.015609194255155_wp*u(x+4, y, z, m:p)

  end if


  if (x .eq. mx-3) then
    u_x(m:p) = 0.013345690454124_wp*u(x-4, y, z, m:p) - 0.106765523632994_wp*u(x-3, y, z, m:p)&
                + 0.330443502642251_wp*u(x-2, y, z, m:p) - 0.946962507015985_wp*u(x-1, y, z, m:p)&
                + 0.318292975073623_wp*u(x, y, z, m:p) + 0.243684445410829_wp*u(x+1, y, z, m:p)&
                + 0.240558007429252_wp*u(x+2, y, z, m:p) - 0.092596590361100_wp*u(x+3, y, z, m:p)
  end if


  if (x .eq. mx-2) then
    u_x(m:p) = 0.026915887850467_wp*u(x-4, y, z, m:p) - 0.197069524204597_wp*u(x-3, y, z, m:p)&
              + 0.589452745232692_wp*u(x-2, y, z, m:p) - 1.258766829227482_wp*u(x-1, y, z, m:p)&
              + 0.306540940575249_wp*u(x, y, z, m:p) + 0.629259919521358_wp*u(x+1, y, z, m:p)&
              - 0.096333139747688_wp*u(x+2, y, z, m:p)
  end if



  if (x .eq. mx-1) then
     u_x(m:p) = 0.039560344515975_wp*u(x-4, y, z, m:p) - 0.030739159126096_wp*u(x-3, y, z, m:p) &
               - 0.105980141988699_wp*u(x-2, y, z, m:p) - 0.393228064437077_wp*u(x-1, y, z, m:p) &
               + 0.029551468764759_wp*u(x, y, z, m:p) + 0.460835552271137_wp*u(x+1, y, z, m:p)

  end if



  if (x .eq. mx) then
     u_x(m:p) = - 0.054021710000147_wp*u(x-5, y, z, m:p) - 0.023993033880123_wp*u(x-4, y, z, m:p)&
                + 0.302855902188632_wp*u(x-3, y, z, m:p) + 0.275607596716317_wp*u(x-2, y, z, m:p)&
                - 2.093702214477299_wp*u(x-1, y, z, m:p) + 1.593253459452621_wp*u(x, y, z, m:p)
  end if



  if ((x .ge. 7) .and. (x .le. mx-6)) then

     

     u_x(m:p) = 0.016666666666667_wp*u(x-4, y, z, m:p) - 0.133333333333333_wp*u(x-3, y, z, m:p) + 0.5_wp*u(x-2, y, z, m:p) &
          & - 1.333333333333333_wp*u(x-1, y, z, m:p) + 0.583333333333333_wp*u(x, y, z, m:p) + 0.4_wp*u(x+1, y, z, m:p)&
          & - 0.033333333333333_wp*u(x+2, y, z, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dbx6 = u_x

end function Dbx6



pure function Dfy6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, my, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, my
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dfy6

  n = size(u_x)

  if (y .eq. 1) then
    u_x(m:p) = -1.593253459452621_wp*u(x, y, z, m:p) + 2.093702214477299*u(x, y+1, z, m:p) &
               - 0.275607596716317_wp*u(x, y+2, z, m:p) - 0.302855902188632_wp*u(x, y+3, z, m:p) &
               + 0.023993033880123_wp*u(x, y+4, z, m:p) + 0.054021710000147_wp*u(x, y+5, z, m:p)
  end if

  if (y .eq. 2) then
    u_x(m:p) = -0.460835552271137_wp*u(x, y-1, z, m:p) - 0.029551468764759_wp*u(x, y, z, m:p)&
                   + 0.393228064437077_wp*u(x, y+1, z, m:p) + 0.105980141988699_wp*u(x, y+2, z, m:p)&
                   + 0.030739159126096_wp*u(x, y+3, z, m:p) - 0.039560344515975_wp*u(x, y+4, z, m:p)
                   

  end if
  if (y .eq. 3) then
     u_x(m:p) = 0.096333139747688_wp*u(x, y-2, z, m:p) - 0.629259919521358_wp*u(x, y-1, z, m:p)&
                  - 0.306540940575249_wp*u(x, y, z, m:p) + 1.258766829227482_wp*u(x, y+1, z, m:p)&
                  - 0.589452745232692_wp*u(x, y+2, z, m:p) + 0.197069524204597_wp*u(x, y+3, z, m:p)&
                  - 0.026915887850467_wp*u(x, y+4, z, m:p)
  end if

  if (y .eq. 4) then

     u_x(m:p) = 0.092596590361100_wp*u(x, y-3, z, m:p) - 0.240558007429252_wp*u(x, y-2, z, m:p)&
                  - 0.243684445410829_wp*u(x, y-1, z, m:p) - 0.318292975073623_wp*u(x, y, z, m:p)&
                  + 0.946962507015985_wp*u(x, y+1, z, m:p) - 0.330443502642251_wp*u(x, y+2, z, m:p)&
                  + 0.106765523632994_wp*u(x, y+3, z, m:p) - 0.013345690454124_wp*u(x, y+4, z, m:p)
                  
  end if



  if (y .eq. 5) then
    u_x(m:p) =  -0.015609194255155_wp*u(x, y-4, z, m:p) + 0.033485250803205_wp*u(x, y-3, z, m:p)&
                   + 0.041897570295664_wp*u(x, y-2, z, m:p) - 0.450132219590077_wp*u(x, y-1, z, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.406282074938901_wp*u(x, y+1, z, m:p)&
                   - 0.550950133911491_wp*u(x, y+2, z, m:p) + 0.146920035709731_wp*u(x, y+3, z, m:p)&
                   - 0.018365004463716_wp*u(x, y+4, z, m:p)
  end if

  if (y .eq. 6) then
     
    u_x(m:p) = - 0.016163101657400_wp*u(x, y-5, z, m:p) + 0.040531470682718_wp*u(x, y-4, z, m:p)&
                   - 0.005338567444210_wp*u(x, y-3, z, m:p) - 0.025427590054054_wp*u(x, y-2, z, m:p)&
                   - 0.357460227911152_wp*u(x, y-1, z, m:p) - 0.572337434946968_wp*u(x, y, z, m:p)&
                   + 1.313958528183954_wp*u(x, y+1, z, m:p) - 0.492734448068983_wp*u(x, y+2, z, m:p)&
                   + 0.131395852818395_wp*u(x, y+3, z, m:p) - 0.016424481602299_wp*u(x, y+4, z, m:p)
  end if


  if (y .eq. my-5) then
    u_x(m:p) = 0.032848963204599_wp*u(x, y-2, z, m:p) - 0.394187558455186_wp*u(x, y-1, z, m:p)&
                 - 0.572337434946968_wp*u(x, y, z, m:p) + 1.257688453771463_wp*u(x, y+1, z, m:p)&
                 - 0.406675342006739_wp*u(x, y+2, z, m:p) + 0.120254802392339_wp*u(x, y+3, z, m:p)&
                 - 0.054367610816546_wp*u(x, y+4, z, m:p) + 0.016775726857039_wp*u(x, y+5, z, m:p)
  end if


  if (y .eq. my-4) then
     u_x(m:p) = 0.036730008927433_wp*u(x, y-2, z, m:p) - 0.399693508760138_wp*u(x, y-1, z, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.303115093827634_wp*u(x, y+1, z, m:p)&
                   - 0.402190050630647_wp*u(x, y+2, z, m:p) + 0.047235828122731_wp*u(x, y+3, z, m:p)&
                   + 0.008331008040049_wp*u(x, y+4, z, m:p)
  end if


  if (y .eq. my-3) then
    u_x(m:p) = -0.020661154127888_wp*u(x, y-2, z, m:p) - 0.327107204245208_wp*u(x, y-1, z, m:p)&
                - 0.318292975073623_wp*u(x, y, z, m:p) + 0.624133691970994_wp*u(x, y+1, z, m:p)&
                + 0.118346128899150_wp*u(x, y+2, z, m:p) - 0.076418487423426_wp*u(x, y+3, z, m:p)
  end if


  if (y .eq. my-2) then
    u_x(m:p) = - 0.008748664712219_wp*u(x, y-3, z, m:p) + 0.061405392278187_wp*u(x, y-2, z, m:p)&
                  - 0.491468255323896_wp*u(x, y-1, z, m:p) - 0.306540940575249_wp*u(x, y, z, m:p)&
                  + 0.885608401570530_wp*u(x, y+1, z, m:p) - 0.140255933237354_wp*u(x, y+2, z, m:p)
  end if



  if (y .eq. my-1) then
     u_x(m:p) = 0.029492540133095_wp*u(x, y-4, z, m:p) + 0.021790841692085_wp*u(x, y-3, z, m:p)&
                  - 0.215422101432619_wp*u(x, y-2, z, m:p) - 0.279404147185598_wp*u(x, y-1, z, m:p)&
                  - 0.029551468764759_wp*u(x, y, z, m:p) + 0.473094335557797_wp*u(x, y+1, z, m:p)

  end if



  if (y .eq. my) then
     u_x(m:p) = - 0.052048915548040_wp*u(x, y-5, z, m:p) - 0.044953974933765_wp*u(x, y-4, z, m:p)&
                + 0.366971721882124_wp*u(x, y-3, z, m:p) + 0.189297839436616_wp*u(x, y-2, z, m:p)&
                - 2.039450367044345_wp*u(x, y-1, z, m:p) + 1.580183696207410_wp*u(x, y, z, m:p)
  end if



  if ((y .ge. 7) .and. (y .le. my-6)) then

     u_x(m:p) = 0.033333333333333_wp*u(x, y-2, z, m:p) - 0.4_wp*u(x, y-1, z, m:p) &
                     - 0.583333333333333_wp*u(x, y, z, m:p)   + 1.333333333333333_wp*u(x, y+1, z, m:p)&
                     - 0.5_wp*u(x, y+2, z, m:p) + 0.133333333333333_wp*u(x, y+3, z, m:p)&
                     - 0.016666666666667_wp*u(x, y+4, z, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dfy6 = u_x

end function Dfy6



pure function Dby6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, my, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, my
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dby6

  n = size(u_x)

  if (y .eq. 1) then
    u_x(m:p) = - 1.580183696207410_wp*u(x, y, z, m:p) + 2.039450367044345_wp*u(x, y+1, z, m:p) &
               - 0.189297839436616_wp*u(x, y+2, z, m:p) - 0.366971721882124_wp*u(x, y+3, z, m:p) &
               + 0.044953974933765_wp*u(x, y+4, z, m:p) + 0.052048915548040_wp*u(x, y+5, z, m:p)
  end if

  if (y .eq. 2) then
    u_x(m:p) = -0.473094335557797_wp*u(x, y-1, z, m:p) + 0.029551468764759_wp*u(x, y, z, m:p)&
               + 0.279404147185598_wp*u(x, y+1, z, m:p) + 0.215422101432619_wp*u(x, y+2, z, m:p)&
               - 0.021790841692085_wp*u(x, y+3, z, m:p) - 0.029492540133095_wp*u(x, y+4, z, m:p)
                   

  end if
  if (y .eq. 3) then
     u_x(m:p) = 0.140255933237354_wp*u(x, y-2, z, m:p) - 0.885608401570530_wp*u(x, y-1, z, m:p) &
                + 0.306540940575249_wp*u(x, y, z, m:p) + 0.491468255323896_wp*u(x, y+1, z, m:p) &
                - 0.061405392278187_wp*u(x, y+2, z, m:p) + 0.008748664712219_wp*u(x, y+3, z, m:p)
  end if

  if (y .eq. 4) then

     u_x(m:p) = 0.076418487423426_wp*u(x, y-3, z, m:p) - 0.118346128899150_wp*u(x, y-2, z, m:p) &
                -0.624133691970994_wp*u(x, y-1, z, m:p) + 0.318292975073623_wp*u(x, y, z, m:p) &
                +0.327107204245208_wp*u(x, y+1, z, m:p) + 0.020661154127888_wp*u(x, y+2, z, m:p)
                  
  end if



  if (y .eq. 5) then
    u_x(m:p) =  - 0.008331008040049_wp*u(x, y-4, z, m:p) - 0.047235828122731_wp*u(x, y-3, z, m:p)&
                + 0.402190050630647_wp*u(x, y-2, z, m:p) - 1.303115093827634_wp*u(x, y-1, z, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.399693508760138_wp*u(x, y+1, z, m:p)&
                - 0.036730008927433_wp*u(x, y+2, z, m:p)
  end if

  if (y .eq. 6) then
     
    u_x(m:p) = - 0.016775726857039_wp*u(x, y-5, z, m:p) + 0.054367610816546_wp*u(x, y-4, z, m:p) &
               - 0.120254802392339_wp*u(x, y-3, z, m:p) + 0.406675342006739_wp*u(x, y-2, z, m:p)&
               - 1.257688453771463_wp*u(x, y-1, z, m:p) + 0.572337434946968_wp*u(x, y, z, m:p) &
               + 0.394187558455186_wp*u(x, y+1, z, m:p) - 0.032848963204599_wp*u(x, y+2, z, m:p)
  end if


  if (y .eq. my-5) then
    u_x(m:p) = 0.016424481602299_wp*u(x, y-4, z, m:p) - 0.131395852818395_wp*u(x, y-3, z, m:p)&
               + 0.492734448068983_wp*u(x, y-2, z, m:p) - 1.313958528183954_wp*u(x, y-1, z, m:p)&
               + 0.572337434946968_wp*u(x, y, z, m:p) + 0.357460227911152_wp*u(x, y+1, z, m:p)&
               + 0.025427590054054_wp*u(x, y+2, z, m:p) + 0.005338567444210_wp*u(x, y+3, z, m:p)&
               - 0.040531470682718_wp*u(x, y+4, z, m:p) + 0.016163101657400_wp*u(x, y+5, z, m:p)
  end if


  if (y .eq. my-4) then
     u_x(m:p) = 0.018365004463716_wp*u(x, y-4, z, m:p) - 0.146920035709731_wp*u(x, y-3, z, m:p)&
                + 0.550950133911491_wp*u(x, y-2, z, m:p) - 1.406282074938901_wp*u(x, y-1, z, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.450132219590077_wp*u(x, y+1, z, m:p)&
                - 0.041897570295664_wp*u(x, y+2, z, m:p) - 0.033485250803205_wp*u(x, y+3, z, m:p)&
                + 0.015609194255155_wp*u(x, y+4, z, m:p)

  end if


  if (y .eq. my-3) then
    u_x(m:p) = 0.013345690454124_wp*u(x, y-4, z, m:p) - 0.106765523632994_wp*u(x, y-3, z, m:p)&
                + 0.330443502642251_wp*u(x, y-2, z, m:p) - 0.946962507015985_wp*u(x, y-1, z, m:p)&
                + 0.318292975073623_wp*u(x, y, z, m:p) + 0.243684445410829_wp*u(x, y+1, z, m:p)&
                + 0.240558007429252_wp*u(x, y+2, z, m:p) - 0.092596590361100_wp*u(x, y+3, z, m:p)
  end if


  if (y .eq. my-2) then
    u_x(m:p) = 0.026915887850467_wp*u(x, y-4, z, m:p) - 0.197069524204597_wp*u(x, y-3, z, m:p)&
              + 0.589452745232692_wp*u(x, y-2, z, m:p) - 1.258766829227482_wp*u(x, y-1, z, m:p)&
              + 0.306540940575249_wp*u(x, y, z, m:p) + 0.629259919521358_wp*u(x, y+1, z, m:p)&
              - 0.096333139747688_wp*u(x, y+2, z, m:p)
  end if



  if (y .eq. my-1) then
     u_x(m:p) = 0.039560344515975_wp*u(x, y-4, z, m:p) - 0.030739159126096_wp*u(x, y-3, z, m:p) &
               - 0.105980141988699_wp*u(x, y-2, z, m:p) - 0.393228064437077_wp*u(x, y-1, z, m:p) &
               + 0.029551468764759_wp*u(x, y, z, m:p) + 0.460835552271137_wp*u(x, y+1, z, m:p)

  end if



  if (y .eq. my) then
     u_x(m:p) = - 0.054021710000147_wp*u(x, y-5, z, m:p) - 0.023993033880123_wp*u(x, y-4, z, m:p)&
                + 0.302855902188632_wp*u(x, y-3, z, m:p) + 0.275607596716317_wp*u(x, y-2, z, m:p)&
                - 2.093702214477299_wp*u(x, y-1, z, m:p) + 1.593253459452621_wp*u(x, y, z, m:p)
  end if



  if ((y .ge. 7) .and. (y .le. my-6)) then

     

     u_x(m:p) = 0.016666666666667_wp*u(x, y-4, z, m:p) - 0.133333333333333_wp*u(x, y-3, z, m:p) + 0.5_wp*u(x, y-2, z, m:p) &
          & - 1.333333333333333_wp*u(x, y-1, z, m:p) + 0.583333333333333_wp*u(x, y, z, m:p) + 0.4_wp*u(x, y+1, z, m:p)&
          & - 0.033333333333333_wp*u(x, y+2, z, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dby6 = u_x

end function Dby6




pure function Dfz6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mz, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mz
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dfz6

  n = size(u_x)

  if (z .eq. 1) then
    u_x(m:p) = -1.593253459452621_wp*u(x, y, z, m:p) + 2.093702214477299*u(x, y, z+1, m:p) &
               - 0.275607596716317_wp*u(x, y, z+2, m:p) - 0.302855902188632_wp*u(x, y, z+3, m:p) &
               + 0.023993033880123_wp*u(x, y, z+4, m:p) + 0.054021710000147_wp*u(x, y, z+5, m:p)
  end if

  if (z .eq. 2) then
    u_x(m:p) = -0.460835552271137_wp*u(x, y, z-1, m:p) - 0.029551468764759_wp*u(x, y, z, m:p)&
                   + 0.393228064437077_wp*u(x, y, z+1, m:p) + 0.105980141988699_wp*u(x, y, z+2, m:p)&
                   + 0.030739159126096_wp*u(x, y, z+3, m:p) - 0.039560344515975_wp*u(x, y, z+4, m:p)
                   

  end if

  if (z .eq. 3) then
     u_x(m:p) = 0.096333139747688_wp*u(x, y, z-2, m:p) - 0.629259919521358_wp*u(x, y, z-1, m:p)&
                  - 0.306540940575249_wp*u(x, y, z, m:p) + 1.258766829227482_wp*u(x, y, z+1, m:p)&
                  - 0.589452745232692_wp*u(x, y, z+2, m:p) + 0.197069524204597_wp*u(x, y, z+3, m:p)&
                  - 0.026915887850467_wp*u(x, y, z+4, m:p)
  end if

  if (z .eq. 4) then

     u_x(m:p) = 0.092596590361100_wp*u(x, y, z-3, m:p) - 0.240558007429252_wp*u(x, y, z-2, m:p)&
                  - 0.243684445410829_wp*u(x, y, z-1, m:p) - 0.318292975073623_wp*u(x, y, z, m:p)&
                  + 0.946962507015985_wp*u(x, y, z+1, m:p) - 0.330443502642251_wp*u(x, y, z+2, m:p)&
                  + 0.106765523632994_wp*u(x, y, z+3, m:p) - 0.013345690454124_wp*u(x, y, z+4, m:p)
                  
  end if



  if (z .eq. 5) then
    u_x(m:p) =  -0.015609194255155_wp*u(x, y, z-4, m:p) + 0.033485250803205_wp*u(x, y, z-3, m:p)&
                   + 0.041897570295664_wp*u(x, y, z-2, m:p) - 0.450132219590077_wp*u(x, y, z-1, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.406282074938901_wp*u(x, y, z+1, m:p)&
                   - 0.550950133911491_wp*u(x, y, z+2, m:p) + 0.146920035709731_wp*u(x, y, z+3, m:p)&
                   - 0.018365004463716_wp*u(x, y, z+4, m:p)
  end if

  if (z .eq. 6) then
     
    u_x(m:p) = - 0.016163101657400_wp*u(x, y, z-5, m:p) + 0.040531470682718_wp*u(x, y, z-4, m:p)&
                   - 0.005338567444210_wp*u(x, y, z-3, m:p) - 0.025427590054054_wp*u(x, y, z-2, m:p)&
                   - 0.357460227911152_wp*u(x, y, z-1, m:p) - 0.572337434946968_wp*u(x, y, z, m:p)&
                   + 1.313958528183954_wp*u(x, y, z+1, m:p) - 0.492734448068983_wp*u(x, y, z+2, m:p)&
                   + 0.131395852818395_wp*u(x, y, z+3, m:p) - 0.016424481602299_wp*u(x, y, z+4, m:p)
  end if


  if (z .eq. mz-5) then
    u_x(m:p) = 0.032848963204599_wp*u(x, y, z-2, m:p) - 0.394187558455186_wp*u(x, y, z-1, m:p)&
                 - 0.572337434946968_wp*u(x, y, z, m:p) + 1.257688453771463_wp*u(x, y, z+1, m:p)&
                 - 0.406675342006739_wp*u(x, y, z+2, m:p) + 0.120254802392339_wp*u(x, y, z+3, m:p)&
                 - 0.054367610816546_wp*u(x, y, z+4, m:p) + 0.016775726857039_wp*u(x, y, z+5, m:p)
  end if


  if (z .eq. mz-4) then
     u_x(m:p) = 0.036730008927433_wp*u(x, y, z-2, m:p) - 0.399693508760138_wp*u(x, y, z-1, m:p)&
                   - 0.593528379527062_wp*u(x, y, z, m:p) + 1.303115093827634_wp*u(x, y, z+1, m:p)&
                   - 0.402190050630647_wp*u(x, y, z+2, m:p) + 0.047235828122731_wp*u(x, y, z+3, m:p)&
                   + 0.008331008040049_wp*u(x, y, z+4, m:p)
  end if


  if (z .eq. mz-3) then
    u_x(m:p) = -0.020661154127888_wp*u(x, y, z-2, m:p) - 0.327107204245208_wp*u(x, y, z+1, m:p)&
                - 0.318292975073623_wp*u(x, y, z, m:p) + 0.624133691970994_wp*u(x, y, z+1, m:p)&
                + 0.118346128899150_wp*u(x, y, z+2, m:p) - 0.076418487423426_wp*u(x, y, z+3, m:p)
  end if


  if (z .eq. mz-2) then
    u_x(m:p) = - 0.008748664712219_wp*u(x, y, z-3, m:p) + 0.061405392278187_wp*u(x, y, z-2, m:p)&
                  - 0.491468255323896_wp*u(x, y, z-1, m:p) - 0.306540940575249_wp*u(x, y, z, m:p)&
                  + 0.885608401570530_wp*u(x, y, z+1, m:p) - 0.140255933237354_wp*u(x, y, z+2, m:p)
  end if



  if (z .eq. mz-1) then
     u_x(m:p) = 0.029492540133095_wp*u(x, y, z-4, m:p) + 0.021790841692085_wp*u(x, y, z-3, m:p)&
                  - 0.215422101432619_wp*u(x, y, z-2, m:p) - 0.279404147185598_wp*u(x, y, z-1, m:p)&
                  - 0.029551468764759_wp*u(x, y, z, m:p) + 0.473094335557797_wp*u(x, y, z+1, m:p)

  end if



  if (z .eq. mz) then
     u_x(m:p) = - 0.052048915548040_wp*u(x, y, z-5, m:p) - 0.044953974933765_wp*u(x, y, z-4, m:p)&
                + 0.366971721882124_wp*u(x, y, z-3, m:p) + 0.189297839436616_wp*u(x, y, z-2, m:p)&
                - 2.039450367044345_wp*u(x, y, z-1, m:p) + 1.580183696207410_wp*u(x, y, z, m:p)
  end if



  if ((z .ge. 7) .and. (z .le. mz-6)) then

     u_x(m:p) = 0.033333333333333_wp*u(x, y, z-2, m:p) - 0.4_wp*u(x, y, z-1, m:p) &
                     - 0.583333333333333_wp*u(x, y, z, m:p)   + 1.333333333333333_wp*u(x, y, z+1, m:p)&
                     - 0.5_wp*u(x, y, z+2, m:p) + 0.133333333333333_wp*u(x, y, z+3, m:p)&
                     - 0.016666666666667_wp*u(x, y, z+4, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dfz6 = u_x

end function Dfz6



pure function Dbz6(hx, u, x, y, z, mbx, pbx, mby, pby, mbz, pbz, mz, m, p)

  implicit none

  integer, intent(in) ::  mbx, pbx, mby, pby, mbz, pbz, mz
  real(kind = wp),intent(in) :: hx        ! Spatial step
  integer,intent(in) :: x, y, z, m, p
  real(kind = wp), dimension(mbx:pbx,mby:pby,mbz:pbz,m:p), intent(in) :: u
  

  integer :: n
  real(kind = wp), dimension(m:p) :: u_x, Dbz6

  n = size(u_x)

  if (z .eq. 1) then
    u_x(m:p) = - 1.580183696207410_wp*u(x, y, z, m:p) + 2.039450367044345_wp*u(x, y, z+1, m:p) &
               - 0.189297839436616_wp*u(x, y, z+2, m:p) - 0.366971721882124_wp*u(x, y, z+3, m:p) &
               + 0.044953974933765_wp*u(x, y, z+4, m:p) + 0.052048915548040_wp*u(x, y, z+5, m:p)
  end if

  if (z .eq. 2) then
    u_x(m:p) = -0.473094335557797_wp*u(x, y, z-1, m:p) + 0.029551468764759_wp*u(x, y, z, m:p)&
               + 0.279404147185598_wp*u(x, y, z+1, m:p) + 0.215422101432619_wp*u(x, y, z+2, m:p)&
               - 0.021790841692085_wp*u(x, y, z+3, m:p) - 0.029492540133095_wp*u(x, y, z+4, m:p)
                   

  end if
  if (z .eq. 3) then
     u_x(m:p) = 0.140255933237354_wp*u(x, y, z-2, m:p) - 0.885608401570530_wp*u(x, y, z-1, m:p) &
                + 0.306540940575249_wp*u(x, y, z, m:p) + 0.491468255323896_wp*u(x, y, z+1, m:p) &
                - 0.061405392278187_wp*u(x, y, z+2, m:p) + 0.008748664712219_wp*u(x, y, z+3, m:p)
  end if

  if (z .eq. 4) then

     u_x(m:p) = 0.076418487423426_wp*u(x, y, z-3, m:p) - 0.118346128899150_wp*u(x, y, z-2, m:p) &
                -0.624133691970994_wp*u(x, y, z-1, m:p) + 0.318292975073623_wp*u(x, y, z, m:p) &
                +0.327107204245208_wp*u(x, y, z+1, m:p) + 0.020661154127888_wp*u(x, y, z+2, m:p)
                  
  end if



  if (z .eq. 5) then
    u_x(m:p) =  - 0.008331008040049_wp*u(x, y, z-4, m:p) - 0.047235828122731_wp*u(x, y, z-3, m:p)&
                + 0.402190050630647_wp*u(x, y, z-2, m:p) - 1.303115093827634_wp*u(x, y, z-1, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.399693508760138_wp*u(x, y, z+1, m:p)&
                - 0.036730008927433_wp*u(x, y, z+2, m:p)
  end if

  if (z .eq. 6) then
     
    u_x(m:p) = - 0.016775726857039_wp*u(x, y, z-5, m:p) + 0.054367610816546_wp*u(x, y, z-4, m:p) &
               - 0.120254802392339_wp*u(x, y, z-3, m:p) + 0.406675342006739_wp*u(x, y, z-2, m:p)&
               - 1.257688453771463_wp*u(x, y, z-1, m:p) + 0.572337434946968_wp*u(x, y, z, m:p) &
               + 0.394187558455186_wp*u(x, y, z+1, m:p) - 0.032848963204599_wp*u(x, y, z+2, m:p)
  end if


  if (z .eq. mz-5) then
    u_x(m:p) = 0.016424481602299_wp*u(x, y, z-4, m:p) - 0.131395852818395_wp*u(x, y, z-3, m:p)&
               + 0.492734448068983_wp*u(x, y, z-2, m:p) - 1.313958528183954_wp*u(x, y, z-1, m:p)&
               + 0.572337434946968_wp*u(x, y, z, m:p) + 0.357460227911152_wp*u(x, y, z+1, m:p)&
               + 0.025427590054054_wp*u(x, y, z+2, m:p) + 0.005338567444210_wp*u(x, y, z+3, m:p)&
               - 0.040531470682718_wp*u(x, y, z+4, m:p) + 0.016163101657400_wp*u(x, y, z+5, m:p)
  end if


  if (z .eq. mz-4) then
     u_x(m:p) = 0.018365004463716_wp*u(x, y, z-4, m:p) - 0.146920035709731_wp*u(x, y, z-3, m:p)&
                + 0.550950133911491_wp*u(x, y, z-2, m:p) - 1.406282074938901_wp*u(x, y, z-1, m:p)&
                + 0.593528379527062_wp*u(x, y, z, m:p) + 0.450132219590077_wp*u(x, y, z+1, m:p)&
                - 0.041897570295664_wp*u(x, y, z+2, m:p) - 0.033485250803205_wp*u(x, y, z+3, m:p)&
                + 0.015609194255155_wp*u(x, y, z+4, m:p)

  end if


  if (z .eq. mz-3) then
    u_x(m:p) = 0.013345690454124_wp*u(x, y, z-4, m:p) - 0.106765523632994_wp*u(x, y, z-3, m:p)&
                + 0.330443502642251_wp*u(x, y, z-2, m:p) - 0.946962507015985_wp*u(x, y, z-1, m:p)&
                + 0.318292975073623_wp*u(x, y, z, m:p) + 0.243684445410829_wp*u(x, y, z+1, m:p)&
                + 0.240558007429252_wp*u(x, y, z+2, m:p) - 0.092596590361100_wp*u(x, y, z+3, m:p)
  end if


  if (z .eq. mz-2) then
    u_x(m:p) = 0.026915887850467_wp*u(x, y, z-4, m:p) - 0.197069524204597_wp*u(x, y, z-3, m:p)&
              + 0.589452745232692_wp*u(x, y, z-2, m:p) - 1.258766829227482_wp*u(x, y, z-1, m:p)&
              + 0.306540940575249_wp*u(x, y, z, m:p) + 0.629259919521358_wp*u(x, y, z+1, m:p)&
              - 0.096333139747688_wp*u(x, y, z+2, m:p)
  end if



  if (z .eq. mz-1) then
     u_x(m:p) = 0.039560344515975_wp*u(x, y, z-4, m:p) - 0.030739159126096_wp*u(x, y, z-3, m:p) &
               - 0.105980141988699_wp*u(x, y, z-2, m:p) - 0.393228064437077_wp*u(x, y, z-1, m:p) &
               + 0.029551468764759_wp*u(x, y, z, m:p) + 0.460835552271137_wp*u(x, y, z+1, m:p)

  end if



  if (z .eq. mz) then
     u_x(m:p) = - 0.054021710000147_wp*u(x, y, z-5, m:p) - 0.023993033880123_wp*u(x, y, z-4, m:p)&
                + 0.302855902188632_wp*u(x, y, z-3, m:p) + 0.275607596716317_wp*u(x, y, z-2, m:p)&
                - 2.093702214477299_wp*u(x, y, z-1, m:p) + 1.593253459452621_wp*u(x, y, z, m:p)
  end if



  if ((z .ge. 7) .and. (y .le. mz-6)) then

     

     u_x(m:p) = 0.016666666666667_wp*u(x, y, z-4, m:p) - 0.133333333333333_wp*u(x, y, z-3, m:p) + 0.5_wp*u(x, y, z-2, m:p) &
          & - 1.333333333333333_wp*u(x, y, z-1, m:p) + 0.583333333333333_wp*u(x, y, z, m:p) + 0.4_wp*u(x, y, z+1, m:p)&
          & - 0.033333333333333_wp*u(x, y, z+2, m:p)
                 
  end if


  u_x = (1.0_wp/hx)*u_x

  Dbz6 = u_x

end function Dbz6


end module metrics
