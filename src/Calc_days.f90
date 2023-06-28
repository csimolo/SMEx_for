    module calendar

    contains

!   -----------------------------------------------------------------------------    
!   Returns end day of every month for leap/non-leap years

      integer(4) function edom(sj,sk)
        use iparamt 
        implicit none
        integer(4) :: sj, sk

        if(ileap.eq.0)then
           edom = 30

        elseif(ileap.eq.1.or.ileap.eq.2)then
         
           if(sj.eq.4.or.sj.eq.6.or.sj.eq.9.or.sj.eq.11)then
              edom = 30
            
           elseif(sj.eq.2)then
              edom = 28
              if(ileap.eq.2)then
                 if(mod(sk,4).eq.0.and.mod(sk,100).gt.0) edom = 29  
                 if(mod(sk,400).eq.0) edom = 29  
              endif

           else
              edom = 31
           endif

        endif

      end function edom

!   -----------------------------------------------------------------------------      
!   Returns total # days in a year

      integer(4) function nid(sk)
        implicit none
        integer(4) :: sk, sj

        nid = 0
        do sj = 1,12
           nid = nid + edom(sj,sk)
        enddo

      end function nid
    
!   -----------------------------------------------------------------------------     
!   Returns 1st day of every season; ss = 1,..4 from MAM to DJF
      
      integer(4) function iws(ss,sk)
        implicit none
        integer(4) :: ss, sk, sj, sm(1:4)

        sm = (/2,5,8,11/)
        iws = 1
        if(ss.ge.1.and.ss.le.4)then
           do sj = 1,sm(ss)
              iws = iws + edom(sj,sk)
           enddo
        endif

      end function iws

    end module calendar
