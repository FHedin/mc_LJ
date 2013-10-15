!
! Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
! All rights reserved.
!
! The 3-clause BSD license is applied to this software.
! see LICENSE.txt
! 
!

subroutine aziz_ne_ne(r,pot)
    ! HFD-B potential from Aziz, Chem. Phys. 130 (1989) p 187
    ! reads r in angstroems gives potential energy in cm-1.
    
    implicit none
    
    !local
    real(8) :: x
    real(8) :: v_star
    real(8) :: x2,x4,x6,x8,x10,f

    !local constants
    real(8),parameter :: a_s=8.9571795d+5
    real(8),parameter :: alpha_s=13.86434671d0
    real(8),parameter :: c_6=1.21317545d0
    real(8),parameter :: c_8=0.53222749d0
    real(8),parameter :: c_10=0.24570703d0
    real(8),parameter :: beta_s=-0.12993822d0
!     real(8),parameter :: beta=-0.0136d0
    real(8),parameter :: d=1.36d0
    real(8),parameter :: epsi=42.25d0*0.695039d0
    real(8),parameter :: r_m=3.091d0
    
    !received
    real(8),intent(in)  :: r
    real(8),intent(out) :: pot
    
    !code starts here

    x=r/r_m

    f=dexp(-(d/x-1.d0)**2)
    if (x.ge.d) f=1.d0

    x2=x*x
    x4=x2*x2
    x6=x2*x4
    x8=x4*x4
    x10=x4*x6

    v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-f*(c_6/x6+c_8/x8+c_10/x10)

    pot=epsi*v_star

    return
end subroutine aziz_ne_ne


subroutine aziz_ar_ne(r,pot)
    !HFD-B potential from Barrow, Aziz, JCP 89, 6189 (1988)
    !reads r in angstroems gives potential energy in cm-1.

    implicit none
    
    !local
    real(8) :: x
    real(8) :: v_star
    real(8) :: x2,x4,x6,x8,x10,x12,f

    !local constants
    real(8),parameter :: a_s=1.651205d+5
    real(8),parameter :: alpha_s=9.69290567d0
    real(8),parameter :: c_6=1.09781826d0
    real(8),parameter :: c_8=0.34284623d0
    real(8),parameter :: c_10=0.30103922d0
    real(8),parameter :: c_12=0.74483225
    real(8),parameter :: beta_s=-2.27380851d0
!     real(8),parameter :: beta=-0.1868d0
    real(8),parameter :: d=1.44d0
    real(8),parameter :: epsi=67.59d0*0.695039d0
    real(8),parameter :: r_m=3.4889d0
    
    !received
    real(8),intent(in)  :: r
    real(8),intent(out) :: pot
    
    !code starts here

    x=r/r_m

    f=dexp(-(d/x-1.d0)**2)
    if (x.ge.d) f=1.d0

    x2=x*x
    x4=x2*x2
    x6=x2*x4
    x8=x4*x4
    x10=x4*x6
    x12=x6*x6

    v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-f*(c_6/x6+c_8/x8+c_10/x10+c_12/x12)

    pot=epsi*v_star

    return
end subroutine aziz_ar_ne


subroutine aziz_ar_ar(r,pot)

    !HFD-B potential from Aziz, JCP 92, 1030 (1990)
    !reads r in angstroems gives potential energy in cm-1.

    implicit none
    
    !local
    real(8) :: x
    real(8) :: v_star
    real(8) :: x2,x4,x6,x8,x10,f

    !local constants
    real(8),parameter :: a_s=1.14211845d+5
    real(8),parameter :: alpha_s=9.00053441d0
    real(8),parameter :: c_6=1.09971113d0
    real(8),parameter :: c_8=0.54511632d0
    real(8),parameter :: c_10=0.39278653d0
    real(8),parameter :: beta_s=-2.60270226d0
!     real(8),parameter :: beta=-0.1840d0
    real(8),parameter :: d=1.04d0
    real(8),parameter :: epsi=143.25d0*0.695039d0
    real(8),parameter :: r_m=3.761d0
    
    !received
    real(8),intent(in)  :: r
    real(8),intent(out) :: pot
    
    !code starts here

    x=r/r_m

    f=dexp(-(d/x-1.d0)**2)
    if (x.ge.d) f=1.d0

    x2=x*x
    x4=x2*x2
    x6=x2*x4
    x8=x4*x4
    x10=x4*x6

    v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-f*(c_6/x6+c_8/x8+c_10/x10)

    pot=epsi*v_star

    return
end subroutine aziz_ar_ar







