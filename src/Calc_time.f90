    module timing

    contains
!   Returns cpu and wall clock time in seconds 
      subroutine mtime (tcpu, twall,sk)
        implicit none

        real(4) :: tcpu
        real(8) :: twall
        character(8) :: date
        character(10) :: time
        character(5) :: zone
        integer(4) :: values(8), j, sk
        integer(4) :: month(12)
      
        call cpu_time(tcpu)
        call date_and_time(date, time, zone, values)

        month  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        sk = values(1)
        if(mod(sk,4).eq.0.and.mod(sk,100).gt.0) month(2)=29
        if(mod(sk,400).eq.0) month(2)=29
        
        if(values(2).gt.1)then
           do j = 1, values(2)-1
              values(3) = values(3) + month(j)
           enddo
        endif
        values(3) = values(3)-1
!      print*,' no.days since 1st Jan: ', values(3)

        twall = dble(values(3))*86400.d0 + dble(values(5))*3600.d0 +   &
             dble(values(6))*60.0d0 +dble(values(7)) +  dble(values(8))*1.d-3        
        
        return
      end subroutine mtime
      
    end module timing
    
