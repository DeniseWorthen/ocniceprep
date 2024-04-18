program ocniceprep

  use ESMF
  use netcdf
  use init_mod   ,     only : nxt, nyt, nlevs, nxr, nyr, outvars, readnml, readcsv
  use init_mod   ,     only : wgtsdir, griddir, ftype, fsrc, fdst, input_file, angvar
  use init_mod   ,     only : do_ocnprep, debug, logunit
  use arrays_mod ,     only : b2d, b3d, rgb2d, rgb3d, dstlon, dstlat, setup_packing
  use arrays_mod ,     only : nbilin2d, nbilin3d, bilin2d, bilin3d
  use utils_mod  ,     only : getfield, packarrays, remap, dumpnc, nf90_err
  use utils_esmf_mod , only : createRH, remapRH, ChkErr, rotremap

  implicit none

  type(ESMF_VM)      :: vm
  character(len=160) :: gridfile
  character(len=160) :: wgtsfile
  character(len=160) :: fout

  ! dimensions, units and variables from source file used in creation of
  ! output netcdf
  real(kind=8), allocatable, dimension(:) :: Layer          !< the vertical grid center
  real(kind=8), allocatable, dimension(:) :: angsrc         !< the rotation angle at the Ct points for the src grid
  real(kind=8), allocatable, dimension(:) :: angdst         !< the rotation angle at the Ct points for the dst grid

  ! work arrays for output netcdf
  real(kind=8), allocatable, dimension(:,:)   :: out2d !< 2D destination grid output array
  real(kind=8), allocatable, dimension(:,:,:) :: out3d !< 3D destination grid output array

  ! Time information
  real(kind=8)       :: timestamp
  character(len= 40) :: timeunit
  character(len= 20) :: vname, vunit
  character(len=120) :: vlong
  integer(kind=4)    :: istep1, myear, mmonth, mday, msec

  character(len=120) :: meshfsrc, meshfdst

  integer :: nvalid
  integer :: n,nn,rc,ncid,varid
  integer :: idimid,jdimid,kdimid,edimid,timid
  integer :: idx1,idx2
  ! debug
  integer :: i,j

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

  ! --------------------------------------------------------
  ! initialize ESMF
  ! --------------------------------------------------------

  call ESMF_Initialize(rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! --------------------------------------------------------
  ! read the nml file and a file containing the list of
  ! variables to be remapped
  ! --------------------------------------------------------

  call readnml
  call readcsv(nvalid)

  ! nn = 0
  ! do j = 1,nyr
  !    do i = 1,nxr
  !       nn = nn+1
  !       if(j.eq.nyr)print *,i,nn
  !    end do
  ! end do
  ! i=114945 : (105,320)
  ! --------------------------------------------------------
  ! create a regrid RH from source to destination
  ! --------------------------------------------------------

  meshfsrc = trim(griddir)//fsrc(3:5)//'/'//'mesh.'//trim(fsrc)//'.nc'
  meshfdst = trim(griddir)//fdst(3:5)//'/'//'mesh.'//trim(fdst)//'.nc'
  print '(a)',trim(meshfsrc),trim(meshfdst)
  call createRH(trim(meshfsrc),trim(meshfdst),rc)

  ! --------------------------------------------------------
  ! read the master grid file and obtain the rotation angle
  ! on the source and destination grids
  ! in ocnpost, anglet is retrieved from CICE's own history
  ! file, here it is retrieved from the tripole grid file, which
  ! has the opposite sense for anglet (same as MOM6).
  ! the rotation formulas assume the same sense as MOM6, so
  ! in ocnpost, this requires for cice that sinrot = -sin(anglet)
  ! here, we need -sin(-anglet), which is sin(anglet), so no
  ! sign change is required
  ! --------------------------------------------------------

  allocate(angsrc(nxt*nyt)); angsrc = 0.0
  allocate(angdst(nxr*nyr)); angdst = 0.0

  gridfile = trim(griddir)//fsrc(3:5)//'/'//'tripole.'//trim(fsrc)//'.nc'
  call nf90_err(nf90_open(trim(gridfile), nf90_nowrite, ncid), 'open: '//trim(gridfile))
  call getfield(trim(gridfile), 'anglet', dims=(/nxt,nyt/), field=angsrc)
  call nf90_err(nf90_close(ncid), 'close: '//trim(gridfile))

  gridfile = trim(griddir)//fdst(3:5)//'/'//'tripole.'//trim(fdst)//'.nc'
  call nf90_err(nf90_open(trim(gridfile), nf90_nowrite, ncid), 'open: '//trim(gridfile))
  call getfield(trim(gridfile), 'anglet', dims=(/nxr,nyr/), field=angdst)
  ! spatial domain of output resolution
  !if (do_ocnprep) then
  !   call getfield(trim(gridfile), 'latCt', dims=(/nxr,nyr/), field=lath)
  !   call getfield(trim(gridfile), 'lonCt', dims=(/nxr,nyr/), field=lonh)
  !   call getfield(trim(gridfile), 'lonCu', dims=(/nxr,nyr/), field=lonq)
  !   call getfield(trim(gridfile), 'latCv', dims=(/nxr,nyr/), field=latq)
  !end if
  call nf90_err(nf90_close(ncid), 'close: '//trim(gridfile))
  ! --------------------------------------------------------
  ! get the 3rd (vertical or ncat) dimension
  ! --------------------------------------------------------

  call nf90_err(nf90_open(trim(input_file), nf90_nowrite, ncid), 'open: '//trim(input_file))
  if (do_ocnprep) then
     call nf90_err(nf90_inq_dimid(ncid, 'Layer', varid), 'get dimension Id: Layer'//trim(input_file))
     call nf90_err(nf90_inquire_dimension(ncid, varid, len=nlevs), 'get dimension Id: Layer'//trim(input_file))
  else
     call nf90_err(nf90_inq_dimid(ncid, 'ncat', varid), 'get dimension Id: ncat'//trim(input_file))
     call nf90_err(nf90_inquire_dimension(ncid, varid, len=nlevs), 'get dimension Id: ncat'//trim(input_file))
  endif
  do n = 1,nvalid
     if (debug) then
        write(logunit,'(a12,i4,a10,3(a6),a2,i4)')trim(outvars(n)%var_name)//', ',outvars(n)%var_dimen, &
             ', '//trim(outvars(n)%var_remapmethod),', '//trim(outvars(n)%var_grid),                   &
             ', '//trim(outvars(n)%var_pair),', '//trim(outvars(n)%var_pair_grid),', ',n
     end if
     if (do_ocnprep) then
        call nf90_err(nf90_inq_varid(ncid, trim(outvars(n)%var_name), varid), 'get variable Id: '//trim(outvars(n)%var_name))
        call nf90_err(nf90_get_att(ncid, varid,  'long_name', outvars(n)%long_name), 'get variable attribute: long_name '//trim(outvars(n)%var_name))
        call nf90_err(nf90_get_att(ncid, varid,      'units', outvars(n)%units), 'get variable attribute: units '//trim(outvars(n)%var_name)        )
     end if
  end do
  call nf90_err(nf90_close(ncid), 'close: '//trim(input_file))

  if (do_ocnprep) then
     ! timestamp
     call nf90_err(nf90_inq_varid(ncid, 'Time', varid), 'get variable Id: Time '//trim(input_file))
     call nf90_err(nf90_get_var(ncid, varid, timestamp), 'get variable: timestamp '//trim(input_file))
     call nf90_err(nf90_get_att(ncid, varid,    'units', timeunit), 'get variable attribute : units '//trim(input_file))
     ! layer
     allocate(Layer(nlevs)) ; Layer = 0.0
     call nf90_err(nf90_inq_varid(ncid, 'Layer', varid), 'get variable Id: Layer '//trim(input_file))
     call nf90_err(nf90_get_var(ncid, varid, layer), 'get variable: Layer '//trim(input_file))
  else
     call nf90_err(nf90_inq_attribute(ncid, nf90_global, 'istep1', istep1), 'get global attribute istep1 '//trim(input_file))
     call nf90_err(nf90_inq_attribute(ncid, nf90_global,  'myear',  myear), 'get global attribute myear '//trim(input_file))
     call nf90_err(nf90_inq_attribute(ncid, nf90_global, 'mmonth', mmonth), 'get global attribute mmonth '//trim(input_file))
     call nf90_err(nf90_inq_attribute(ncid, nf90_global,   'mday',   mday), 'get global attribute mday '//trim(input_file))
     call nf90_err(nf90_inq_attribute(ncid, nf90_global,   'msec',   msec), 'get global attribute msec '//trim(input_file))
  end if


 ! --------------------------------------------------------
 ! create packed arrays for mapping and remap packed arrays
 ! to the destination grid
 ! --------------------------------------------------------

  call setup_packing(nvalid,outvars)

  ! 2D bilin
  if (allocated(bilin2d)) then
     call packarrays(trim(input_file), trim(wgtsdir)//fsrc(3:5)//'/', cos(angsrc), sin(angsrc),         &
          b2d, dims=(/nxt,nyt/), nflds=nbilin2d, fields=bilin2d)
     rgb2d = 0.0
     call remapRH(src_field=bilin2d, dst_field=rgb2d)

     if (debug) then
        write(logunit,'(a)')'remap 2D fields bilinear with RH '
        write(logunit,'(a)')'packed min/max values, mapped min/max values'
        do n = 1,nbilin2d
           write(logunit,'(i4,a10,3(a2,a6),4g14.4)')n,trim(b2d(n)%var_name),'  ',                       &
                trim(b2d(n)%var_grid),'  ',trim(b2d(n)%var_pair),'  ', trim(b2d(n)%var_pair_grid),      &
                minval(bilin2d(:,n)), maxval(bilin2d(:,n)),minval(rgb2d(:,n)), maxval(rgb2d(:,n))
        end do
        call dumpnc(trim(ftype)//'.'//trim(fsrc)//'.bilin2d.nc', 'bilin2d', dims=(/nxt,nyt/),           &
             nflds=nbilin2d, field=bilin2d)
        call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin2d.nc', 'rgbilin2d', dims=(/nxr,nyr/),       &
             nflds=nbilin2d, field=rgb2d)
     end if
  end if

  ! 3D bilin
  if (allocated(bilin3d))then
     call packarrays(trim(input_file), trim(wgtsdir)//fsrc(3:5)//'/', cos(angsrc), sin(angsrc),         &
          b3d, dims=(/nxt,nyt,nlevs/), nflds=nbilin3d, fields=bilin3d)
     rgb3d = 0.0
     !call remapRH(src_field=bilin3d, dst_field=rgb3d)

     if (debug) then
        write(logunit,'(a)')'remap 3D fields bilinear with RH'
        write(logunit,'(a)')'packed min/max values,mapped min/max values'
        do n = 1,nbilin3d
           write(logunit,'(i4,a10,3(a2,a6),4g14.4)')n,trim(b3d(n)%var_name),'  ',                       &
                trim(b3d(n)%var_grid),'  ',trim(b3d(n)%var_pair),'  ', trim(b3d(n)%var_pair_grid),      &
                minval(bilin3d(:,:,n)), maxval(bilin3d(:,:,n)),minval(rgb3d(:,:,n)), maxval(rgb3d(:,:,n))
        end do
        call dumpnc(trim(ftype)//'.'//trim(fsrc)//'.bilin3d.nc', 'bilin3d', dims=(/nxt,nyt,nlevs/),     &
             nk=nlevs, nflds=nbilin3d, field=bilin3d)
        call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin3d.nc', 'rgbilin3d', dims=(/nxr,nyr,nlevs/), &
             nk=nlevs, nflds=nbilin3d, field=rgb3d)
     end if
  end if

  !--------------------------------------------------------
  ! rotate on Ct from EW->IJ and remap back to native staggers
  !--------------------------------------------------------

  if (allocated(bilin2d)) then
     call rotremap(trim(wgtsdir)//fdst(3:5)//'/', b2d, cos(angdst), sin(angdst), dims=(/nxr,nyr/), &
          nflds=nbilin2d, fields=rgb2d)

     call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin2d.Bu.ij.nc', 'rgbilin2d', dims=(/nxr,nyr/),       &
          nflds=nbilin2d, field=rgb2d)
  end if

    if (allocated(bilin3d)) then
       call rotremap(trim(wgtsdir)//fdst(3:5)//'/', b3d, cos(angdst), sin(angdst), dims=(/nxr,nyr,nlevs/), &
            nflds=nbilin3d, fields=rgb3d)

       call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin3d.Bu.ij.nc', 'rgbilin3d', dims=(/nxr,nyr,nlevs/),       &
            nk=nlevs, nflds=nbilin3d, field=rgb3d)
    end if
#ifdef test
  ! --------------------------------------------------------
  ! write the mapped fields
  ! --------------------------------------------------------

  allocate(out2d(nxr,nyr)); out2d = 0.0
  allocate(out3d(nxr,nyr,nlevs)); out3d = 0.0

  fout = trim(ftype)//'.'//trim(fdst)//'.nc'
  if (debug) write(logunit, '(a)')'output file: '//trim(fout)

  gridfile = trim(griddir)//fdst(3:5)//'/'//'tripole.'//trim(fdst)//'.nc'
  if (ocn_prep) then
     call setup_ocnrestart(trim(input_file),trim(fout),trim(gridfile))
  else
     call setup_icerestart(trim(input_file),trim(fout),trim(gridfile))
  end if

  call nf90_err(nf90_create(trim(fout), nf90_clobber, ncid), 'create: '//trim(fout))
  if (ocn_prep) then
     call nf90_err(nf90_def_dim(ncid, 'lonh', nxr,  idimid), 'define dimension: lonh')
     call nf90_err(nf90_def_dim(ncid, 'lath', nyr,  jdimid), 'define dimension: lath')
     call nf90_err(nf90_def_dim(ncid, 'lonq', nxr, qidimid), 'define dimension: lonq')
     call nf90_err(nf90_def_dim(ncid, 'latq', nxr, qjdimid), 'define dimension: latq')
     call nf90_err(nf90_def_dim(ncid, 'Layer',  nlevs, kdimid), 'define dimension: Layer')
     call nf90_err(nf90_def_dim(ncid, 'time', nf90_unlimited, timid), 'define dimension: time')
     ! define the time variable
     call nf90_err(nf90_def_var(ncid, 'time', nf90_double, (/timid/), varid), 'define variable: time')
     call nf90_err(nf90_put_att(ncid, varid,    'units', trim(timeunit)), 'put variable attribute: units')
     call nf90_err(nf90_put_att(ncid,  varid, 'calendar', trim(timecal)), 'put variable attribute: calendar')
     ! spatial grid
     call nf90_err(nf90_def_var(ncid, 'lonh', nf90_double,  (/idimid/), varid), 'define variable: lonh')
     call nf90_err(nf90_put_att(ncid, varid, 'units', 'degrees_east'),  'put variable attribute: units')
     call nf90_err(nf90_def_var(ncid, 'lath', nf90_double,  (/jdimid/), varid), 'define variable: lath' )
     call nf90_err(nf90_put_att(ncid, varid, 'units', 'degrees_north'), 'put variable attribute: units')
     call nf90_err(nf90_def_var(ncid, 'lonq', nf90_double, (/qidimid/), varid), 'define variable: lonq')
     call nf90_err(nf90_put_att(ncid, varid, 'units', 'degrees_east'),  'put variable attribute: units')
     call nf90_err(nf90_def_var(ncid, 'latq', nf90_double, (/qjdimid/), varid), 'define variable: latq' )
     call nf90_err(nf90_put_att(ncid, varid, 'units', 'degrees_north'), 'put variable attribute: units')
     ! vertical grid
     call nf90_err(nf90_def_var(ncid, 'Layer', nf90_double,  (/kdimid/), varid), 'define variable: Layer')
     call nf90_err(nf90_put_att(ncid, varid,    'units', 'm'), 'put variable attribute: units')
  else
     call nf90_err(nf90_def_dim(ncid, 'ni', nxr, idimid), 'define dimension: ni')
     call nf90_err(nf90_def_dim(ncid, 'nj', nyr, jdimid), 'define dimension: nj')
     call nf90_err(nf90_def_dim(ncid, 'ncat',  nlevs, kdimid), 'define dimension: ncat')
  end ifu

  if (allocated(b2d)) then
     do n = 1,nbilin2d
        vname = trim(b2d(n)%var_name)
        vunit = trim(b2d(n)%units)
        vlong = trim(b2d(n)%long_name)
        call nf90_err(nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,timid/), varid), 'define variable: '// vname)
        if (do_ocnprep) then
           call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units')
           call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name')
        end if
     enddo
  end if
  if (allocated(c2d)) then
     do n = 1,nconsd2d
        vname = trim(c2d(n)%var_name)
        vunit = trim(c2d(n)%units)
        vlong = trim(c2d(n)%long_name)
        call nf90_err(nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,timid/), varid), 'define variable: '// vname)
        if (do_ocnprep) then
           call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units' )
           call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name' )
        end if
     enddo
  end if
  if (allocated(b3d)) then
     do n = 1,nbilin3d
        vname = trim(b3d(n)%var_name)
        vunit = trim(b3d(n)%units)
        vlong = trim(b3d(n)%long_name)
        call nf90_err(nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,kdimid,timid/), varid), 'define variable: '// vname)
        if (do_ocnprep) then
           call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units' )
           call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name' )
        end if
     enddo
  end if
  if (.not. ocn_prep) then
     call nf90_err(nf90_put_att(ncid, nf90_global, 'istep1', istep1), 'put global attribute istep1')
     call nf90_err(nf90_put_att(ncid, nf90_global,  'myear',  myear), 'put global attribute myear')
     call nf90_err(nf90_put_att(ncid, nf90_global, 'mmonth', mmonth), 'put global attribute mmonth')
     call nf90_err(nf90_put_att(ncid, nf90_global,   'mday',   mday), 'put global attribute mday')
     call nf90_err(nf90_put_att(ncid, nf90_global,   'msec',   msec), 'put global attribute msec')
  end if
  call nf90_err(nf90_enddef(ncid), 'enddef: '// trim(fout))

  if (ocn_prep) then
     ! dimensions
     call nf90_err(nf90_inq_varid(ncid, 'lonh', varid), 'get variable Id: lonh')
     call nf90_err(nf90_put_var(ncid,   varid, lonCt),     'put variable: lonh')
     call nf90_err(nf90_inq_varid(ncid, 'lath', varid), 'get variable Id: lath')
     call nf90_err(nf90_put_var(ncid,   varid, latCt),     'put variable: lath')
     call nf90_err(nf90_inq_varid(ncid, 'lonq', varid), 'get variable Id: lonq')
     call nf90_err(nf90_put_var(ncid,   varid, lonCu),     'put variable: lonq')
     call nf90_err(nf90_inq_varid(ncid, 'latq', varid), 'get variable Id: latq')
     call nf90_err(nf90_put_var(ncid,   varid, latCv),     'put variable: latq')
     ! time
     call nf90_err(nf90_inq_varid(ncid, 'time', varid), 'get variable Id: time')
     call nf90_err(nf90_put_var(ncid, varid, timestamp), 'put variable: time')
     ! vertical
     call nf90_err(nf90_inq_varid(ncid, 'z_l', varid), 'get variable Id: z_l')
     call nf90_err(nf90_put_var(ncid, varid, z_l)    , 'put variable: z_l')
  end if

  if (allocated(rgb2d)) then
     do n = 1,nbilin2d
        out2d(:,:) = reshape(rgb2d(:,n), (/nxr,nyr/))
        vname = trim(b2d(n)%var_name)
        call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable Id: '//vname)
        call nf90_err(nf90_put_var(ncid,   varid, out2d), 'put variable: '//vname)
     end do
  end if
  if (allocated(rgc2d)) then
     do n = 1,nconsd2d
        out2d(:,:) = reshape(rgc2d(:,n), (/nxr,nyr/))
        vname = trim(c2d(n)%var_name)
        call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable Id: '//vname)
        call nf90_err(nf90_put_var(ncid,   varid, out2d), 'put variable: '//vname)
     end do
  end if
  if (allocated(rgb3d)) then
     do n = 1,nbilin3d
        out3d(:,:,:) = reshape(rgb3d(:,:,n), (/nxr,nyr,nlevs/))
        vname = trim(b3d(n)%var_name)
        call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable Id: '//vname)
        call nf90_err(nf90_put_var(ncid,   varid, out3d), 'put variable: '//vname)
     end do
  end if
  call nf90_err(nf90_close(ncid), 'close: '// trim(fout))
  write(logunit,'(a)')trim(fout)//' done'
#endif
  stop

end program ocniceprep
