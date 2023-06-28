!   Compute annual/seasonal extremes (min and max values) and frequencies of exceeding reference quantiles 
!   (the 1st, the 99th and 99.9th) for grid-point daily data read from input netcdf files
!   ------------------------------------------------------------------------------------------------------
!   September 2022 
!   Claudia Simolo - ISAC-CNR
!   Institute of Atmospheric Sciences and Climate 
!   National Research Council of Italy
!   c.simolo@isac.cnr.it  
!
!   Use: 
!   netcdf library (https://github.com/Unidata/netcdf-fortran);
!   quicksort_real_F77.F (https://github.com/jasonanema/Quicksort_Fortran77);
!   auxlib.a (only: Calc_time.f90, Fix_parameters.f90, Calc_days.f90, compare.f, quicksort_real_F77.F, Calc_quantiles.f90)
!
!   Input:
!   - list of time-ordered nc files to be processed, with the same variables' and grid structure; files include  
!   daily data of one or more years and cover the period for calculation of the extremes and reference quantiles;
!   - list of boundary years of each file (YYYY_i YYYY_f) ;
!   - nc names of dimensions (time, longitude, latitude);
!   - nc names of required variables (longitude, latitude, temperature);
!   - data calendar (via the 'ileap' index);
!   - period of interest for the extremes;
!   - period for reference quantiles;
!   - data selection (annual/seasonal, absolute values/anomalies);
!   - name and path of output files 
!
!   Output:
!   - extremes values (min and max) as a function of time for every grid point, printed as single precision real 
!   on ascii files named 'f1out-smp_lon_lat'; 'f1out' (max a20) is defined as input, 'smp' (a3, input) is a label 
!   specifying data selection (year/season), 'lon' = 'xxxyy' (xxx.yy in 0,360) and 'lat' = '+/-xxyy' (+/- xx.yy 
!   in -90,+90); 
!   - exceedance frequencies (relative to the reference quantiles) as a function of time for every grid point, 
!   printed as above on ascii files named 'f2out-smp_lon_lat'; missing values are set NaN;
!   - grid-point reference quantiles (1st;99th;99.9th pct) with coordinates ;
!   - list of output files with grid-point coordinates ('ls_f1out-smp'; 'ls_f2out-smp') ; 
!   - list of grid points with undefined exceedance frequencies ('f2out-smp.log') 
!
!   Computational settings:  
!   degC = conversion to degree Celsius ; 
!   Nlatx = # latitudes processed at once

    program EXtremes_nc
  
      use iparamt
      use calendar
      use cquants  
      use timing
      
      implicit none
      include '/usr/local/include/netcdf.inc'      

!     Input and storing arrays       
      character(200) :: ls_ncfiles, ls_ncdates
      character(20) :: dmtime, dmlon, dmlat, vrlon, vrlat, vrmain
      character(200), allocatable :: fname(:)
      integer(4) :: nset, iset,ios, i,j,k, gap, yG1(10), yG2(10)
      integer(4) :: ncid, retval, ndims, nvars, ngatts, unlimdimid
      integer(4), parameter :: nmax_dims = 100, nmax_vars = 100, nmax_atts = 100
      integer(4) :: ncount(nmax_dims), start(nmax_dims)                           
      character(50) :: Dname(nmax_dims),Vname(nmax_vars)                                
      integer(4) :: dim, Dlen(nmax_dims), Did(nmax_dims)
      integer(4) :: var, Vtype, VND, Vid(nmax_vars), VDids(nf_max_var_dims),Vatt 
      integer(4), allocatable :: y0(:),yF(:),numbf(:), Ntm(:)
      integer(4) :: vidlon,vidlat,vidmain, itm, ilat, ilon, Ntmf, Nlat, Nlon
      real(8), allocatable :: lons_in(:), tlons_in(:), lats_in(:)
      real, allocatable :: tas_in(:,:,:,:)

!     Options & output settings
      integer, parameter :: Nlatx =64
      logical, parameter :: degC =.true.    
      real(4) :: const
      logical :: cdev
      integer(4) :: y1sp, y2sp, isn, H1,H2, nrchar, enchar, o1char, o2char
      integer(4) :: y1, yN, nnr, xlen, jjmin, ngp, inr
      integer(4) :: nbk, jbk, ix1, ix2, ixx, nxx
      integer(4) :: id, j1,j2, iq, iw, ii1, ii2, jj
      real, allocatable :: norm(:), xa(:,:), xh(:)
      integer(4), parameter :: nq = 1  
      real(4):: qxs(nq), qxlo(3),qxhi(3), mthlo(2),mthup(2) 
      real(4), parameter :: xlo = -1.e+20, xup = 1.e+20
      character(200) :: nrpath, nrname, o1path, o2path
      character(20) :: f1out, f2out                                                  
      character(3) :: smp
      character(5) :: chlon, chlat
      integer(4) :: f1char, f2char,flon, flat
      character(50) :: fileout1, fileout2
      integer(4) :: dcnlo(2), dcnup(2), dtot
      real(4) :: vmin, vmax, eplo(2), epup(2)
      real(4) :: t0c, tfc
      real(8) :: t0w, tfw
      integer(4) :: yy0, yyf

      call mtime(t0c, t0w,yy0)

      print*,' '
      print*,'Start program: EXtremes_nc'
      print*,''
      print*,'--> SET INPUT '   
      print*,' '

!     Enter list of input files (with path)
      read(*,'(a)') ls_ncfiles
!     Enter list of boundary years of of input files (YYYY_i YYYY_f) 
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
      
!     Enter boundary years for calculation of extremes
      read*, y1sp, y2sp
!     Select data samples via 'isn' index: 0/1/2/3/4 for year/mam/jja/son/djf
      read*, isn
      if(isn.eq.0) smp ='fyr'
      if(isn.eq.1) smp ='mam'
      if(isn.eq.2) smp ='jja'
      if(isn.eq.3) smp ='son'
      if(isn.eq.4) smp ='djf'
      if(isn.lt.0.or.isn.gt.4)then
         print*,'Error: isn out of range > Stop'
         print*,'Choose isn in [0,4]! '
         stop
      endif
!     Enter boundary years for calculation of reference quantiles 
      read*, H1, H2

!     Enter cdev = T for anomalies, else F
      read*,cdev
!     Enter path for normals -not used if cdev = F'
      read(*,'(a)')nrpath
      do i = 1,200
         if(nrpath(i:i).eq.' ') exit
      enddo
      nrchar = i-1
!     Enter base name of normals 
      read(*,'(a)') nrname
      do i = 1,50
         if(nrname(i:i).eq.' ') exit
      enddo
      enchar = i-1

!     Enter path (max a200) of output files #1
      read(*,'(a)') o1path
      do i = 1,200
         if(o1path(i:i).eq.' ') exit
      enddo
      o1char = i-1
!     Define base name (max a20) of output files #1 
      read(*,'(a)') f1out
      do i = 1,20
         if(f1out(i:i).eq.' ') exit
      enddo
      f1char = i-1

!     Enter path (max a200) of output files #2
      read(*,'(a)') o2path
      do i = 1,200
         if(o2path(i:i).eq.' ') exit
      enddo
      o2char = i-1
!     Define base name (max a20) of output files #2 
      read(*,'(a)') f2out
      do i = 1,20
         if(f2out(i:i).eq.' ') exit
      enddo
      f2char = i-1

      write(*,'(1x,a)')'Input stored in file inp_ee'
      print*,' '

!     Switch from degree K to °C   
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
      print*,' '

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
         print*,'Error: period for extremes not included in input file(s)'
         print*,'Stop. Check input list!'
         stop
      endif

      if(H1.lt.y0(1).or.H2.gt.yF(nset))then
         print*,'Error: period for quantiles not included in input file(s)'
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
         enddo
      endif
      
      if(gap.gt.0)then

         do i =1,gap
            if(y2sp.lt.yG1(i).or.y1sp.gt.yG2(i))then
               cycle
            else
               write(*,'(1x,a,2I5)')'Error: gap in the period of extremes: ',yG1(i),yG2(i)
               print*,'Stop. Check input list!'
               stop
            endif
         enddo
         do i =1,gap
            if(H2.lt.yG1(i).or.H1.gt.yG2(i))then
               cycle
            else
               write(*,'(1x,a,2I5)')'Error: gap in the period of quantiles: ',yG1(i),yG2(i)
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
         write(*,'(3x,a,i10)')'ncid = ',ncid
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
         
      enddo !iset
      
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
      
      nnr = nidx(ileap)
      allocate(norm(nnr))
      allocate(xa(nnr,y1:yN)) 
      
      xlen = (H2-H1+1)*nnr
      allocate(xh(xlen))
      
      if(isn.eq.0) jjmin = int(0.5*real(xlen))
      if(isn.gt.0) jjmin = int(0.5*real(xlen)/4.)
      
!     Store coordinates into prog arrays       
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
      print*,'--> COMPUTING EXTREMEs '

      open(20, file = o1path(1:o1char)//'ls_'//f1out(1:f1char)//'-'//smp, status ='unknown')  
      open(22, file = o2path(1:o2char)//'ls_'//f2out(1:f2char)//'-'//smp, status ='unknown')  
      open(24, file = o2path(1:o2char)//'THs-'//f2out(1:f2char)//'-'//smp, status ='unknown')
      open(26, file = o2path(1:o2char)//f2out(1:f2char)//'-'//smp//'.log', status ='unknown')

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
      
!     Cycle on grid points inside each block
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
               fileout1 = f1out(1:f1char)//'-'//smp//'_'//chlon//'_'//chlat
               fileout2 = f2out(1:f2char)//'-'//smp//'_'//chlon//'_'//chlat
               write(20,*) fileout1,lons_in(ilon),lats_in(ilat)
               write(22,*) fileout2,lons_in(ilon),lats_in(ilat)

!     Read grid-point normals from files
               if(cdev)then
                  open(28, file = nrpath(1:nrchar)//nrname(1:enchar)//'_'//chlon//'_'//chlat, status='old')
                  do inr = 1,nnr
                     read(28,*,iostat=ios) id, norm(inr)
                     if(ios.lt.0) exit 
                  enddo
                  close(28)
               endif
               
!     Store grid-point data into xa(id,k). If ileap=2 and k non-leap xa(366,k) is nan
               do iset = 1,nset
                  itm = 0            
               
                  do k = y0(iset),yF(iset)
                     id = 0          
                  
                     j1 = 1
                     j2 = 12
                     do j = j1,j2                     
                        do i = 1, edom(j,k)                            
                           itm = itm + 1
                           id = id +1
                        
                           xa(id,k) = tas_in(ilon,ixx,itm,iset) - const
                           if( tas_in(ilon,ixx,itm,iset).eq.missv_1.or.tas_in(ilon,ixx,itm,iset)  &  
                                .eq.missv_2 ) xa(id,k) = rnan
                        
!     Subtract day-of-the-year normals
                           if(.not.cdev) cycle
                           if(nid(k).ne.nidx(ileap).and.id.ge.60)then
                              xa(id,k) = xa(id,k) - norm(id+1)
                           else
                              xa(id,k) = xa(id,k) - norm(id)
                           endif

                        enddo !i
                     enddo   !j
                  enddo         !k
               
                  if(itm.ne.Ntm(iset))then
                     print*,'Error: last value itm =',itm,' .ne.',Ntm(iset)
                     print*,'Stop'
                     stop
                  endif

               enddo !iset

!     Set xa(id,k) to nan for k within gaps   
               if(gap.gt.0)then
                  do i=1,gap
                     do k = yG1(i),yG2(i)
                        xa(1:nnr,k) = rnan
                     enddo
                  enddo
               endif
            
!     Compute reference quantiles of annual/seasonal grid-point distributions        
               jj = 0
               do k = H1,H2
                  if(k.lt.y1) cycle
                  if(k.gt.yN) exit
               
!     Full year sample                  
                  ii1 = 1
                  ii2 = nid(k)

!     MAM / JJA / SON sample                  
                  if(isn.ge.1.and.isn.le.3)then
                     ii1 = iws(isn,k)
                     ii2 = iws(isn+1,k)-1 
!     DJF sample (Jan-Feb)    
                  elseif(isn.eq.4)then
                     ii2 = iws(1,k)-1
                  endif

                  do iw = ii1,ii2
                     if(isnan(xa(iw,k))) cycle
                     jj = jj+1          
                     xh(jj) = xa(iw,k)
                  enddo

                  if(isn.lt.4) cycle

!     DJF sample (Dec)
                  if(k.eq.y1) cycle
                  do iw = iws(4,k-1),nid(k-1) 
                     if(isnan(xa(iw,k-1))) cycle
                     jj = jj+1
                     xh(jj) = xa(iw,k-1)
                  enddo         !iw
                  
               enddo            !k
            
               if(jj.ge.jjmin)then
                  call quantile(2,nq,jj,xh,qxs,qxlo,qxhi)
                  mthlo(1) = qxlo(2)
                  mthlo(2) = qxlo(3)
                  mthup(1) = qxhi(2)
                  mthup(2) = qxhi(3)
               
               else
                  mthlo(1:2) = rnan
                  mthup(1:2) = rnan
                  write(26,*) lons_in(ilon), lats_in(ilat), rnan
               endif

               write(24,*) lons_in(ilon), lats_in(ilat), mthlo(1), (mthup(iq),iq=1,2)

!     Compute grid-point annual/seasonal extremes 
               open(30, file = o1path(1:o1char)//fileout1, position = 'append')
               open(32, file = o2path(1:o2char)//fileout2, position = 'append')

               do k = y1sp,y2sp
                  if(k.lt.y1) cycle
                  if(k.gt.yN) exit

                  vmin = xup
                  vmax = xlo
                  dtot = 0
                  dcnup = 0
                  dcnlo = 0
               
!     Full year sample                  
                  ii1 = 1
                  ii2 = nid(k)

!     MAM / JJA / SON sample                  
                  if(isn.ge.1.and.isn.le.3)then
                     ii1 = iws(isn,k)
                     ii2 = iws(isn+1,k)-1 
!     DJF sample (Jan-Feb)    
                  elseif(isn.eq.4)then
                     ii2 = iws(1,k)-1
                  endif

                  do iw = ii1,ii2
                     if(isnan(xa(iw,k))) cycle
                  
                     vmin = min(vmin,xa(iw,k))
                     vmax = max(vmax,xa(iw,k))

                     dtot = dtot+1
                     do iq=1,2
                        if(xa(iw,k).le.mthlo(iq)) dcnlo(iq) = dcnlo(iq)+1
                        if(xa(iw,k).ge.mthup(iq)) dcnup(iq) = dcnup(iq)+1
                     enddo

                  enddo

                  if(isn.eq.4)then
!     DJF sample (Dec)
                     if(k.eq.y1) cycle
                     do iw = iws(4,k-1),nid(k-1) 
                        if(isnan(xa(iw,k-1))) cycle

                        vmin = min(vmin,xa(iw,k-1))
                        vmax = max(vmax,xa(iw,k-1))

                        dtot = dtot+1
                        do iq=1,2
                           if(xa(iw,k-1).le.mthlo(iq)) dcnlo(iq) = dcnlo(iq)+1
                           if(xa(iw,k-1).ge.mthup(iq)) dcnup(iq) = dcnup(iq)+1
                        enddo
                     enddo         !iw
                  endif
            
                  write(30,*) k, vmin, vmax
                  write(32,*) k, real(dcnlo(1))/real(dtot),(real(dcnup(iq))/real(dtot), iq=1,2)
               enddo            !k
               close(30)
               close(32)
            
            enddo            !ilon
         enddo               !ilat

         deallocate(tas_in)
      enddo  !jbk
         
      close(20)                                                      
      close(22)
      close(24)
      
      print*,' '
      write(*,'(1x,a,1x,i10)')'--> END PROCESSING: total # gridpoints =',ngp
      write(*,'(5x,a,1x,a)')'Output printed on: ',o1path(1:o1char)//f1out(1:f1char)//'-'//smp//'_*'
      write(*,'(5x,a,1x,a)')'Output printed on: ',o2path(1:o2char)//f2out(1:f2char)//'-'//smp//'_*'
      print*,' '
      
!     Close input files
      do iset = 1,nset
         write(*,'(3x,a,i10)')'CLOSE file #',numbf(iset)
         retval = nf_close(numbf(iset))
         if(retval.ne.nf_noerr) call handle_err(retval)         
      enddo
      print*,' '
      
!     Deallocate arrays
      deallocate(fname)
      deallocate(y0)
      deallocate(yF)
      deallocate(numbf)
      deallocate(Ntm)

      deallocate(lons_in)
      deallocate(tlons_in)
      deallocate(lats_in)
            
      deallocate(xa)
      deallocate(xh)
      deallocate(norm)
      
      print*,'End program: EXtremes_nc'
      print*,' '

      call mtime (tfc, tfw, yyf)
      if(yyf.gt.yy0)then
         print*,'Execution time not available'
      else
         print*,'Total time for execution '
         write(*,'(1x,a,f12.2)')'CPU  time (s) =', tfc-t0c
         write(*,'(1x,a,f12.2)')'Wall time (s) =', tfw-t0w
         print*,' '
      endif

    end program EXtremes_nc
    
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

