    module cquants
  
    contains

!   -----------------------------------------------------------------------------      
!     Returns quantiles of data distributions
!     Uses quicksort_real_F77.F (https://github.com/jasonanema/Quicksort_Fortran77) with compare.f90
!     IN: ikp = 0 / 1 / 2 (interpolation type: round to nearest integer / preserve median / NIST); 
!         noq (# required quantiles); ndat (# data); vec (data vector)
!     OUT: qqs(noq) (equally spaced quantiles); qlo(3) (q_{P=0.02;0.01;0.001}); qhi(3) (q_{P=0.98,0.99,0.999})

      subroutine quantile(ikp,noq,ndat,vec,qqs,qlo,qhi)
        implicit none      

        integer(4) :: ikp, i, noq, ndat, nor
        real(4) :: vec(ndat),step, xs, xnor, xd
        real(4), parameter :: eplo(3) = (/ 0.020, 0.010, 0.001 /), ephi(3) = (/ 0.980, 0.990, 0.999 /)
        real(4) :: qqs(noq), qlo(3), qhi(3)

        external compar
        integer(2) compar
        
        if(noq.gt.ndat)then
           print*,' Error from subrt < percentile >: '
           print*,' # percentile > # data. Stop'
           print*,' '
           stop
        endif
      
        call quicksort_real_F77(vec,ndat,compar)
        
!     Equally spaced quantiles    
        step = 1.0/real(noq+1)
        do i = 1,noq
           xs = step*i
           if(ikp.eq.0.or.ikp.eq.1) xnor = xs*ndat + 0.5 
           if(ikp.eq.2) xnor = xs*ndat + xs 
           nor = int(xnor)
           xd = xnor-nor
           if(nor.eq.0)then
              qqs(i) = vec(1)
           elseif(nor.eq.ndat)then
              qqs(i) = vec(ndat)
           else
              if(ikp.eq.0) qqs(i) = vec(nor) 
              if(ikp.eq.1.or.ikp.eq.2) qqs(i) = vec(nor) + xd*(vec(nor+1)-vec(nor))     
           endif
        enddo

!     Low-end tail quantiles qlo(1:3)
        do i=1,3
           xnor = eplo(i)*real(ndat) + eplo(i)
           nor = int(xnor)
           xd = xnor-nor
           if(nor.eq.0)then
              qlo(i) = vec(1)
           elseif(nor.eq.ndat)then
              qlo(i) = vec(ndat)
           else
              qlo(i) = vec(nor) + xd*(vec(nor+1)-vec(nor))
           endif
        enddo
        
!     High-end tail quantiles qhi(1:3)
        do i=1,3
           xnor = ephi(i)*real(ndat) + ephi(i)
           nor = int(xnor)
           xd = xnor-nor
           if(nor.eq.0)then
              qhi(i) = vec(1)
           elseif(nor.eq.ndat)then
              qhi(i) = vec(ndat)
           else
              qhi(i) = vec(nor) + xd*(vec(nor+1)-vec(nor))
           endif
        enddo
      
      end subroutine quantile

    end module cquants
    
