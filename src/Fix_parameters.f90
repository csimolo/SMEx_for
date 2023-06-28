    module iparamt
      integer(4) :: ileap, dctime
      integer(4), parameter :: nidx(0:2) = (/ 360, 365, 366 /) 
      integer(4), parameter:: ledom(12) = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)  
      real, parameter :: kelv = 273.15000, missv_1 = -1.e+32, missv_2 = 1.e+20
      real, parameter :: rnan = transfer(Z'FFC00000', 1.)
      real(8), parameter :: dnan = transfer(X'FFF8000000000000', 1.d0)
      save
    end module iparamt
