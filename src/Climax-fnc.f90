!   Compute day-of-the-year normals for grid-point daily temperature data
!   from input netcdf files
!   ---------------------------------------------------------------------
!   March 2022 
!   Claudia Simolo - ISAC-CNR
!   Institute of Atmospheric Sciences and Climate 
!   National Research Council of Italy
!   c.simolo@isac.cnr.it  
!
!   Use:
!   netcdf library (https://github.com/Unidata/netcdf-fortran);
!   dfftpack library (https://github.com/fortran-lang/fftpack);
!   auxlib.a (only: Calc_time.f90, Fix_parameters.f90, Calc_days.f90)
!
!   Input: 
!   - list of time-ordered nc files to be processed, with the same variables' and grid structure; 
!   files include daily data of one or more years and cover the full reference period for normals 
!   - list of boundary years of each file (YYYY_i YYYY_f) ;
!   - nc names of dimensions (time, longitude, latitude);
!   - nc names of required variables (longitude, latitude, temperature);
!   - data calendar (via the 'ileap' index);
!   - reference period;
!   - name and path of output files
!
!   Output: 
!   - filtered normals for every day of the year and grid point, printed as single precision real
!   on ascii files named 'fout_lon_lat'; 'fout' (max a20) is defined as input, 'lon'='xxxyy' with 
!   xxx.yy in (0,360) and 'lat'='+/-xxyy' with +/-xx.yy in (-90,+90); missing values are set NaN;
!   - list of output files with grid-point coordinates ('ls_fout');
!   - list of grid-points with undefined normals ('fout.log')
!
!   Settings:  
!   degC = conversion to degrees Celsius 
!   Nlatx = # latitudes processed at once;
!   dw = half width (days) of the moving window about each day of the year;
!   nhar = # Fourier components for filtered normals

    program Climax_nc

      use iparamt
      use calendar
      use timing

      implicit none
      include '/usr/local/include/netcdf.inc'
      
!     Input and storing arrays 
      character(200) :: ls_ncfiles, ls_ncdates
      character(20) :: dmtime, dmlon, dmlat, vrlon, vrlat, vrmain
      character(200), allocatable :: fname(:)
      integer(4) :: iset, nset, ios, i,j,k, gap, yG1(10), yG2(10)
      integer(4) :: ncid, retval, ndims, nvars, ngatts, unlimdimid
      integer(4), parameter :: nmax_dims = 100, nmax_vars = 100, nmax_atts = 100
      integer(4) :: ncount(nmax_dims), start(nmax_dims)                           
      character(50) :: Dname(nmax_dims),Vname(nmax_vars)                                
      integer(4) :: dim, Dlen(nmax_dims), Did(nmax_dims)
      integer(4) :: var, Vtype, VND, Vid(nmax_vars), VDids(nf_max_var_dims),Vatt
      integer(4), allocatable :: y0(:), yF(:), numbf(:), Ntm(:)
      integer(4) :: vidlon,vidlat,vidmain, itm, ilat, ilon, Ntmf, Nlat, Nlon
      real(8), allocatable :: lons_in(:), tlons_in(:), lats_in(:)
      real, allocatable :: tas_in(:,:,:,:)

!     Output & settings 
      integer, parameter :: Nlatx =64
      integer, parameter :: dw = 15, nhar = 4
      logical, parameter :: degC = .true.
      integer(4) :: y1sp, y2sp, ochar 
      real(4) :: const
      integer(4) :: y1, yN, nft, jjmin, ngp
      integer(4) :: nbk, jbk, ix1, ix2, ixx, nxx
      integer(4) :: id, j1,j2, imax, jj, iw, ii, ki, ims, lng
      real, allocatable :: xa(:,:), rawnor(:)
      real(8), allocatable :: idft(:), rft(:), Ak(:), Bk(:), idfnor(:)
      character(200) :: opath
      character(20) :: fout                                                  
      character(5) :: chlon, chlat
      integer(4) :: flon, flat, fchar
      character(50) :: fileout
      real(4) :: t0c, tfc
      real(8) :: t0w, tfw
      integer(4) :: yy0, yyf

      call mtime(t0c, t0w,yy0)

      print*,' '
      print*,'Start program: Climax_nc'
      print*,' '
      print*,'--> SET INPUT '   
      print*,' '

!     Enter list of input nc files (with path)
      read(*,'(a)') ls_ncfiles
!     Enter list of boundary years of input files (YYYY_i YYYY_f) 
      read(*,'(a)') ls_ncdates

!     Enter nc names of dimensions (time, longitude, latitude)
      read*, dmtime
      read*, dmlon
      read*, dmlat
!     Enter nc names of required variables (longitude, latitude, main variable)
      read*, vrlon
      read*, vrlat
      read*, vrmain
!     Enter data calendar: ileap = 0,1,2 for 360_days, 365_days (non-leap), 365/366_days (leap)
      read*, ileap
      
!     Enter reference period for normals
      read*, y1sp, y2sp

!     Enter path (max a200) of output files
      read(*,'(a)') opath
      do i = 1,200
         if(opath(i:i).eq.' ') exit
      enddo
      ochar = i-1
!     Enter base name (max a20) of output files 
      read(*,'(a)') fout
      do i = 1,20
         if(fout(i:i).eq.' ') exit
      enddo
      fchar = i-1

      write(*,'(1x,a)')'Input stored in file inp_nn'
      print*,' '

!     Switch from degrees K to °C
      const = 0.0
      if(degC) const = kelv

!     Count input files      
      nset = 0
      open(13,file = ls_ncfiles, status = 'old')
      do
         read(13,'(a)',iostat=ios)
         if(ios.lt.0) exit
         nset = nset+1
      enddo
      close(13)
      write(*,'(1x,a,i3)') '# input file(s): ',nset
      print*,''
      
!     Allocate arrays 
      allocate(fname(nset))
      allocate(y0(nset))
      allocate(yF(nset))
      allocate(numbf(nset))
      allocate(Ntm(nset))
      
      open(13,file = ls_ncfiles, status = 'old')
      open(14,file = ls_ncdates, status = 'old')
      do iset = 1,nset
         read(13,'(a)',iostat=ios) fname(iset)
         if(ios.lt.0) exit
         read(14,*) y0(iset), yF(iset)

         if(iset.gt.1)then
            if(y0(iset).lt.y0(iset-1))then
               print*,'Error: input file(s) not time ordered'
               print*,'Stop. Check input list!'
               stop
            endif
         endif

      enddo
      close(13)
      close(14)

      if(y1sp.lt.y0(1).or.y2sp.gt.yF(nset))then
         print*,'Error: reference period for normals not included in input file(s)'
         print*,'Stop. Check input list!'
         stop
      endif

!     Check time gaps and overlaps in input files
      gap = 0
      if(nset.gt.1)then
         do iset = 2,nset
            if(y0(iset).gt.yF(iset-1)+1)then
               gap = gap+1
               yG1(gap) = yF(iset-1)+1
               yG2(gap) = y0(iset)-1
            elseif(y0(iset).lt.yF(iset-1)+1)then
               write(*,'(1x,2(a,i3))')'Warning: overlap between set ',iset-1,' and ',iset
               print*,''
            endif
            if(yF(iset).gt.y2sp) exit
         enddo
      endif
      
      if(gap.gt.0)then
         do i =1,gap
            if(y2sp.lt.yG1(i).or.y1sp.gt.yG2(i))then
               cycle
            else
               write(*,'(1x,a,2I5)')'Error: gap in the period of normals: ',yG1(i),yG2(i)
               print*,'Stop. Check input list!'
               stop
            endif
         enddo
      endif
            
      print*,'--> GET DATA '   
      print*,' '
      
!     Cycle on input files to read dimensions and variables       
      do iset = 1,nset
                  
!     Open files             
         write(*,'(3x,2a)')'OPEN: ', fname(iset) 
         retval = nf_open(fname(iset),nf_nowrite,numbf(iset))
         if(retval.ne.nf_noerr) call handle_err(retval)
         ncid = numbf(iset)
         write(*,'(3x,a,i10)')'ncid =',ncid
         print*,' '
         retval = nf_inq(ncid, ndims, nvars, ngatts, unlimdimid)
         if(retval.ne.nf_noerr) call handle_err(retval)

!     Check dimensions
         do dim = 1,ndims
            retval = nf_inq_dim(ncid, dim, Dname(dim), Dlen(dim))
            if(retval.ne.nf_noerr) call handle_err(retval)
            retval = nf_inq_dimid(ncid, Dname(dim), Did(dim))
            if(retval.ne.nf_noerr) call handle_err(retval)
            
            write(*,'(3x,a,i2,2(2x,a),1x,i5)')'Dim #', Did(dim), Dname(dim),' lenght:', Dlen(dim)
            if(Dname(dim).eq.dmtime) Ntm(iset) = Dlen(dim)
            if(Dname(dim).eq.dmlon) Nlon = Dlen(dim)
            if(Dname(dim).eq.dmlat) Nlat = Dlen(dim)
         enddo
         print*,' '

!     Check variables type and shape         
         do var = 1,nvars
            retval = nf_inq_var(ncid, var, Vname(var), Vtype, VND, VDids, Vatt)
            if(retval.ne.nf_noerr) call handle_err(retval)
            retval = nf_inq_varid(ncid, Vname(var), Vid(var))
            if(retval.ne.nf_noerr) call handle_err(retval)

            if(Vname(var).eq.vrlon.or.Vname(var).eq.vrlat.or.Vname(var).eq.vrmain)then
               write(*,'(3x,a,i2,2x,a)')'Var #',Vid(var), Vname(var)
               if(Vtype.eq.4) write(*,'(3x,a,i3,a)')'Vtype = ', Vtype, ' - int '
               if(Vtype.eq.5) write(*,'(3x,a,i3,a)')'Vtype = ', Vtype, ' - float '
               if(Vtype.eq.6) write(*,'(3x,a,i3,a)')'Vtype = ', Vtype, ' - double '
               
               if(VND.eq.0) write(*,'(3x,a,i3,a)')'VND = ', VND, ' - scalar '
               if(VND.eq.1) write(*,'(3x,a,i3,a)')'VND = ', VND, ' - vector '
               if(VND.eq.2) write(*,'(3x,a,i3,a)')'VND = ', VND, ' - matrix '
               if(VND.ge.3) write(*,'(3x,a,i3,a)')'VND = ', VND, ' - array '

               write(*,'(3x,a,i2,2x,a,10i4)') '# dims = ',VND,' VDids: ', (VDids(dim), dim = 1,VND)
               print*,''
            endif

            if(Vname(var).eq.vrlon) vidlon = Vid(var) 
            if(Vname(var).eq.vrlat) vidlat = Vid(var) 
            if(Vname(var).eq.vrmain) vidmain = Vid(var) 
         enddo
         
         if(yF(iset).gt.y2sp) exit
         if(iset.eq.nset) exit
      enddo !iset

      nset = iset
      Ntmf = maxval(Ntm(1:nset))

      y1 = minval(y0(1:nset))
      yN = maxval(yF(1:nset))

!     Partion space into latitude blocks 
      if(mod(Nlat,Nlatx).eq.0)then
         nbk = Nlat/Nlatx
      else
         nbk = int(real(Nlat)/real(Nlatx)) +1
      endif
      
!     Allocate arrays
      allocate(lons_in(Nlon))
      allocate(tlons_in(Nlon))
      allocate(lats_in(Nlat))
      
      nft = nidx(ileap)
      jjmin = y2sp-y1sp+1

      allocate(xa(nft,y1:yN)) 
      allocate(rawnor(nft))
      allocate(idfnor(nft))
      allocate(rft(nft))
      allocate(Ak(nft))
      allocate(Bk(nft))
      allocate(idft(nft))
      
!     Store coordinates into program arrays 
      ncid = numbf(1)

      write(*,'(3x,a,i2,2x,a)')'Get var #',vidlon,vrlon   !lon   
      retval = nf_get_var_double(ncid,vidlon,lons_in)  
      if(retval.ne.nf_noerr) call handle_err(retval)      
      do ilon = 1,Nlon
         tlons_in(ilon) = lons_in(ilon)
         if(lons_in(ilon).lt.0.0) tlons_in(ilon) = lons_in(ilon)+360.0
      enddo

      write(*,'(3x,a,i2,2x,a)')'Get var #',vidlat,vrlat   !lat
      retval = nf_get_var_double(ncid,vidlat,lats_in) 
      if(retval.ne.nf_noerr) call handle_err(retval)      
      
      print*,' '
      print*,'--> COMPUTING day-of-the-year NORMALs'

      open(20, file = opath(1:ochar)//'ls_'//fout(1:fchar), status = 'unknown')
      open(22, file = opath(1:ochar)//fout(1:fchar)//'.log', status ='unknown')  
      
!     Cycle on latitude blocks
      ngp = 0
      do jbk = 1,nbk
         ix1 = 1+(jbk-1)*Nlatx
         ix2 = min(ix1+Nlatx-1,Nlat)
         nxx = ix2-ix1+1   
         allocate(tas_in(Nlon,nxx,Ntmf,nset))
         
         print*,' '
         write(*,'(3x,a,i2,2x,a)')'Get var #',vidmain,vrmain
         write(*,'(3x,a,i2,a,i3,a)')'block #',jbk,':',nxx,' lats'
         
         do iset = 1,nset
            ncid = numbf(iset)
            write(*,'(3x,a,i3,2x,i10)')'file #',iset,ncid
            start(1) = 1
            start(2) = ix1
            start(3) = 1
            ncount(1) = Nlon
            ncount(2) = nxx
            ncount(3) = Ntm(iset)
            retval = nf_get_vara_real(ncid,vidmain,start,ncount,tas_in(:,1:nxx,:,iset))
            if(retval.ne.nf_noerr) call handle_err(retval)
         enddo !iset
         print*,' '
         
!     Cycle on grid points in each block
         do ilat = ix1,ix2 
            write(*,'(3x,a,1x,f8.4)')'>> lat circle ', lats_in(ilat)
            flat = nint(lats_in(ilat)*100)
            if(flat.lt.0) write(chlat,'(a1,i4.4)')'-',abs(flat)
            if(flat.ge.0) write(chlat,'(a1,i4.4)')'+',flat

            ixx = ilat - (jbk-1)*Nlatx
            
            do ilon = 1,Nlon 
               flon = nint(tlons_in(ilon)*100)
               write(chlon,'(i5.5)')flon
            
               ngp = ngp + 1
               fileout = fout(1:fchar)//'_'//chlon//'_'//chlat
               write(20,*) fileout,lons_in(ilon), lats_in(ilat)
            
!     Store grid-point data into xa(id,k). If ileap=2, xa(1:366,k) and xa(60,k)= NaN for k nonleap
               do iset = 1,nset               
                  itm = 0                   
               
                  do k = y0(iset),yF(iset)
                     id = 0                 

                     j1 = 1
                     j2 = 12
                     do j = j1,j2
                     
                        imax = edom(j,k)
                        if(ileap.eq.2) imax = ledom(j)
                        do i = 1,imax                                                   
                           id = id +1
                           if(i.gt.edom(j,k))then
                              xa(id,k) = rnan
                              cycle
                           endif

                           itm = itm + 1                     
                           xa(id,k) = tas_in(ilon,ixx,itm,iset) - const    
                           if( tas_in(ilon,ixx,itm,iset).eq.missv_1.or.tas_in(ilon,ixx,itm,iset)  &  
                                .eq.missv_2 ) xa(id,k) = rnan
                        enddo !i
                        
                     enddo   !j
                  enddo         !k
               
                  if(itm.ne.Ntm(iset))then
                     print*,' Error: last value itm =',itm,' .ne.',Ntm(iset)
                     print*,' Stop'
                     stop
                  endif

               enddo !iset

!     Compute raw grid-point day-of-the-year normals, as the average in y1sp-y2sp using 2dw+1 day windows   
               do id = 1,nft
                  jj = 0                    
                  rawnor(id) = 0.0
               
                  do k = y1sp,y2sp
                     
                     do iw = id-dw,id+dw 
                        if(iw.le.0)then     
                           ii = iw + nft
                           if(k.eq.y1) cycle
                           ki = k-1
                        elseif(iw.gt.nft)then 
                           ii = iw-nft
                           if(k.eq.yN) cycle
                           ki = k+1
                        else
                           ii = iw
                           ki = k
                        endif

                        if(ii.eq.60.and.nid(ki).ne.nft)then
                           if(id.gt.60) ii = id-dw-1
                           if(id.lt.60) ii = id+dw+1
                        endif
                     
                        if( isnan(xa(ii,ki))) cycle
                        rawnor(id) = rawnor(id)+ xa(ii,ki)
                        jj = jj+1
                     enddo !iw

                  enddo  !k
               
                  if(jj.ge.jjmin)then
                     rawnor(id) = rawnor(id)/real(jj)
                  else
                     rawnor(id) = rnan
                  endif
                  
               enddo  !id
            
!     Filter out noise: compute idft of raw normals by FFT 
               open(24, file = opath(1:ochar)//fileout, status = 'unknown')

               ims = 0
               do ii = 1,nft
                  rft(ii) = dble(rawnor(ii))
                  if( isnan(rawnor(ii)) ) ims = ims+1
               enddo
               if(ims.ge.1) write(22,*) lons_in(ilon), lats_in(ilat), rnan
               
               call dFour(nft,rft,nhar,lng,Ak,Bk,idft)
               do ii = 1, nft
                  idfnor(ii) = idft(ii)
                  write(24,*) ii, real(idfnor(ii)) 
                  if(ii.gt.lng) cycle 
               enddo
               close(24)
               
            enddo            !ilon
         enddo               !ilat

         deallocate(tas_in)
      enddo  !jbk

      close(20)                                                      
      close(22)
      
      print*,' '
      write(*,'(1x,a,1x,i10)')'--> END PROCESSING: total # gridpoints =',ngp
      write(*,'(5x,a,1x,a)')'Output printed on: ',opath(1:ochar)//fout(1:fchar)//'_*'
      print*,' '
      
!   Close input files
      do iset = 1,nset
         write(*,'(3x,a,i10)')'CLOSE file #',numbf(iset)
         retval = nf_close(numbf(iset))
         if(retval.ne.nf_noerr) call handle_err(retval)         
      enddo
      print*,' '
      
!  Deallocate arrays
      deallocate(fname)
      deallocate(y0)
      deallocate(yF)
      deallocate(numbf)
      deallocate(Ntm)

      deallocate(lons_in)
      deallocate(tlons_in)
      deallocate(lats_in)
      
      deallocate(xa)
      deallocate(rawnor)
      deallocate(idfnor)
      deallocate(rft)
      deallocate(Ak)
      deallocate(Bk)
      deallocate(idft)

      print*,'End program: Climax_nc'
      print*,' '

      call mtime (tfc, tfw, yyf)
      if(yyf.gt.yy0)then
         print*,'Execution time not available'
      else
         print*,'Execution time: '
         write(*,'(1x,a,f12.2)')'CPU  time (s) =', tfc-t0c
         write(*,'(1x,a,f12.2)')'Wall time (s) =', tfw-t0w
         print*,' '
      endif
      
      stop
    end program Climax_nc
    
!   Functions & Subroutines  ----------------------------------------------------
    
!   -----------------------------------------------------------------------------  
!   Returns errors in reading nc files    

    subroutine handle_err(errcode)
      implicit none 
      include '/usr/local/include/netcdf.inc'
      
      integer(4) :: errcode
        
      print *, 'Error: ',nf_strerror(errcode)
      stop 2
    end subroutine handle_err

!   -----------------------------------------------------------------------------      
!   Returns Fourier components and inverse discrete transform

!   IN: nft, rft(1:nft), nhar; nft, nhar unchanged on exit; rft changed to fourier coeffs
!   OUT: lng, Ak(1:lng), Bk(1:lng), idft(1:nft) 
    
    subroutine dFour(nft,rft,nhar,lng,Ak,Bk,idft)
      implicit none
      integer(4):: nft, nws, lng, nhar, ii, it
      real(8):: rft(nft), Ak(nft), Bk(nft), wsave(2*nft+15)
      real(8):: rmean, xarg, xsum, idft(nft)

      Ak = 0.d0
      Bk = 0.d0
      
      call dffti(nft,wsave)
      call dfftf(nft,rft,wsave)

      rmean = rft(1)/dble(nft)
      if(mod(nft,2).eq.0)then   !nft even
         lng = int(real(nft)/2.)
         Ak(lng) = rft(nft)
      else                      !nft odd
         lng = int(real(nft+1)/2.)
      endif
      
!     non-normalized Fourier coeffs k = 1, n/2 -1 n even; k = 1,(n-1)/2, n odd
      do ii = 2,lng
         Ak(ii-1) = rft(2*ii-2)
         Bk(ii-1) = -rft(2*ii-1)              
      enddo
      
!     compute idft up to nhar harmonics (nhar .le. nlg)
      if(nhar.gt.lng)then
         print*,' >>  * Error * from subroutine dFour: change nhar!'
         print*,' >>  Stop '
         print*,' '
         stop
      endif

      xarg = 2.d0*acos(-1.d0)/dble(nft)
      do it = 0,nft-1
         xsum = 0.d0
         do ii = 1, nhar
            xsum = xsum + (Ak(ii)*dcos(ii*xarg*it) + &
                 Bk(ii)*dsin(ii*xarg*it))
         enddo
         idft(it+1) = rmean + 2.d0*xsum/dble(nft)
      enddo
      
      return
    end subroutine dFour

!   -----------------------------------------------------------------------------        
