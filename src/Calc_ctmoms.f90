    module ctmoms

!  Returns sample estimators of the first four central moments from a set {xdt} of kd data
!  using k-statistics (https://mathworld.wolfram.com/k-Statistic.html);
!  output are single precision real; undefined values are reported as NaN.

    contains

      subroutine moments(kd,xdt,xm1,usig2,usig,uskew,ukurt)
        use iparamt
        implicit none
        
        integer(4) :: kd, i
        real :: xdt(kd)
        real :: ykd, ykd1, ykd2, ykd3
        real :: xm1, xm2, xm3, xm4 
        real :: xsig2, usig2, xsig, usig
        real :: xskew, uskew, xkurt, ukurt
        
        if(kd.le.0)then
           xm1 = rnan 
           usig2 = rnan 
           usig = rnan 
           uskew = rnan 
           ukurt = rnan 
         
        else
           ykd = real(kd)      
           ykd1 = real(kd-1)
           ykd2 = real(kd-2)
           ykd3 = real(kd-3)
           
           xm1 = 0.0                
           do i = 1,kd
              xm1 = xm1 + xdt(i)
           enddo
           xm1 = xm1/ykd                                           !<x>           

!     standard deviation        
           if(kd.ge.2)then
              xm2 = 0.0
              do i = 1,kd
                 xm2 = xm2 + (xdt(i)-xm1)**2
              enddo
              xm2 = xm2/ykd                                        !sample mom
              xsig2 = xm2                                          
              xsig = sqrt(xsig2)                                   
              usig2 = ykd*xm2/ykd1                                 !variance
              if(usig2.lt.1.e-20) usig2 = 0.0
              usig = sqrt(usig2)                                   !SD
           else
              usig2 = rnan 
              usig = rnan 
           endif
           
!     coeffient of skewness 
           if(usig2.gt.0.0.and.kd.ge.3)then
              xm3 = 0.0
              do i = 1,kd
                 xm3 = xm3 + (xdt(i)-xm1)**3
              enddo
              xm3 = xm3/ykd                                        !sample mom 
              xskew = xm3/sqrt(xm2**3)                             
              uskew = xskew*sqrt(ykd*ykd1)/ykd2                    !g1 
           else
              uskew = rnan 
           endif
         
!     coefficient of excess kurtosis 
           if(usig2.gt.0.0.and.kd.ge.4)then
              xm4 = 0.0
              do i = 1,kd
                 xm4 = xm4 + (xdt(i)-xm1)**4
              enddo
              xm4 = xm4/ykd                                        !sample mom 
              xkurt = xm4/(xm2**2) - 3.0                           
              ukurt = ((ykd+1.0)*xkurt + 6.0)*ykd1/(ykd2*ykd3)     !g2 
           else
              ukurt = rnan 
           endif
           
        endif 
      
        return
      end subroutine moments

    end module ctmoms
    
