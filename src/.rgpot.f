        subroutine aziz_ne_ne(r,pot)

c HFD-B potential from Aziz, Chem. Phys. 130 (1989) p 187
c reads r in angstroems gives potential energy in cm-1.


        implicit none
        real*8 pot,x,r
        real*8 a_s,alpha_s,c_6,c_8,c_10
        real*8 beta_s,beta
        real*8 epsilon,d,r_m,v_star
        real*8 x2,x4,x6,x8,x10,f

        a_s=8.9571795d+5
        alpha_s=13.86434671d0
        c_6=1.21317545d0
        c_8=0.53222749d0
        c_10=0.24570703d0
        beta_s=-0.12993822d0
        beta=-0.0136d0
        d=1.36d0

        epsilon=42.25d0*0.695039d0
        r_m=3.091d0

        x=r/r_m

        f=dexp(-(d/x-1.d0)**2)
        if (x.ge.d) f=1.d0

        x2=x*x
        x4=x2*x2
        x6=x2*x4
        x8=x4*x4
        x10=x4*x6

        v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-
     &         f*(c_6/x6+c_8/x8+c_10/x10)

        pot=epsilon*v_star

        return
        end




        subroutine aziz_ar_ne(r,pot)

c HFD-B potential from Barrow, Aziz, JCP 89, 6189 (1988)
c reads r in angstroems gives potential energy in cm-1.


        implicit none
        real*8 pot,x,r
        real*8 a_s,alpha_s,c_6,c_8,c_10,c_12
        real*8 beta_s,beta
        real*8 epsilon,d,r_m,v_star
        real*8 x2,x4,x6,x8,x10,x12,f

        a_s=1.651205d+5
        alpha_s=9.69290567d0
        c_6=1.09781826d0
        c_8=0.34284623d0
        c_10=0.30103922d0
        c_12=0.74483225
        beta_s=-2.27380851d0
        beta=-0.1868d0
        d=1.44d0

        epsilon=67.59d0*0.695039d0
        r_m=3.4889d0

        x=r/r_m

        f=dexp(-(d/x-1.d0)**2)
        if (x.ge.d) f=1.d0

        x2=x*x
        x4=x2*x2
        x6=x2*x4
        x8=x4*x4
        x10=x4*x6
        x12=x6*x6

        v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-
     &         f*(c_6/x6+c_8/x8+c_10/x10+c_12/x12)

        pot=epsilon*v_star

        return
        end



        function aziz_ar_ar(r,pot)

c HFD-B potential from Aziz, JCP 92, 1030 (1990)
c reads r in angstroems gives potential energy in cm-1.


        implicit none
        real*8 pot,x,r
        real*8 a_s,alpha_s,c_6,c_8,c_10
        real*8 beta_s,beta
        real*8 epsilon,d,r_m,v_star
        real*8 x2,x4,x6,x8,x10,f

        a_s=1.14211845d+5
        alpha_s=9.00053441d0
        c_6=1.09971113d0
        c_8=0.54511632d0
        c_10=0.39278653d0
        beta_s=-2.60270226d0
        beta=-0.1840d0
        d=1.04d0

        epsilon=143.25d0*0.695039d0
        r_m=3.761d0

        x=r/r_m

        f=dexp(-(d/x-1.d0)**2)
        if (x.ge.d) f=1.d0

        x2=x*x
        x4=x2*x2
        x6=x2*x4
        x8=x4*x4
        x10=x4*x6

        v_star=a_s*dexp(-alpha_s*x+beta_s*x**2)-
     &         f*(c_6/x6+c_8/x8+c_10/x10)

        pot=epsilon*v_star

        return
        end







