module restarts_mod

  use netcdf
  use init_mod   , only : nlevs, nxr, nyr
  use init_mod   , only : debug, logunit
  use arrays_mod , only : b2d, b3d, rgb2d, rgb3d
  use arrays_mod , only : nbilin2d, nbilin3d, bilin2d, bilin3d
  use utils_mod  , only : nf90_err

  implicit none

contains
  !----------------------------------------------------------
  ! set up ice restart file
  !----------------------------------------------------------
  subroutine setup_icerestart(fin, fout)

    character(len=*), intent(in) :: fin, fout

    integer :: istep1, myear, mmonth, mday, msec
    integer :: ncid,varid,n
    integer :: idimid,jdimid,kdimid
    character(len=20) :: vname
    character(len=40) :: subname = 'setup_icerestart'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    call nf90_err(nf90_open(trim(fin), nf90_nowrite, ncid), 'open: '//trim(fin))
    ! get the global attributes for time from the input restart file
    call nf90_err(nf90_get_att(ncid, nf90_global, 'istep1', istep1), 'get global attribute istep1 '//trim(fin))
    call nf90_err(nf90_get_att(ncid, nf90_global,  'myear',  myear), 'get global attribute myear '//trim(fin))
    call nf90_err(nf90_get_att(ncid, nf90_global, 'mmonth', mmonth), 'get global attribute mmonth '//trim(fin))
    call nf90_err(nf90_get_att(ncid, nf90_global,   'mday',   mday), 'get global attribute mday '//trim(fin))
    call nf90_err(nf90_get_att(ncid, nf90_global,   'msec',   msec), 'get global attribute msec '//trim(fin))
    call nf90_err(nf90_close(ncid), 'close: '//trim(fin))

    ! create the restart file
    call nf90_err(nf90_create(trim(fout), nf90_clobber, ncid), 'create: '//trim(fout))
    call nf90_err(nf90_def_dim(ncid, 'ni', nxr, idimid), 'define dimension: ni')
    call nf90_err(nf90_def_dim(ncid, 'nj', nyr, jdimid), 'define dimension: nj')
    call nf90_err(nf90_def_dim(ncid, 'ncat',  nlevs, kdimid), 'define dimension: ncat')

    if (allocated(b2d)) then
       do n = 1,nbilin2d
          vname = trim(b2d(n)%var_name)
          call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid/), varid), 'define variable: '// vname)
       enddo
    end if
    !if (allocated(c2d)) then
    !   do n = 1,nconsd2d
    !      vname = trim(c2d(n)%var_name)
    !      call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid/), varid), 'define variable: '// vname)
    !   enddo
    !end if
    if (allocated(b3d)) then
       do n = 1,nbilin3d
          vname = trim(b3d(n)%var_name)
          call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,kdimid/), varid), 'define variable: '// vname)
       enddo
    end if

    call nf90_err(nf90_put_att(ncid, nf90_global, 'istep1', istep1), 'put global attribute istep1')
    call nf90_err(nf90_put_att(ncid, nf90_global,  'myear',  myear), 'put global attribute myear')
    call nf90_err(nf90_put_att(ncid, nf90_global, 'mmonth', mmonth), 'put global attribute mmonth')
    call nf90_err(nf90_put_att(ncid, nf90_global,   'mday',   mday), 'put global attribute mday')
    call nf90_err(nf90_put_att(ncid, nf90_global,   'msec',   msec), 'put global attribute msec')
    call nf90_err(nf90_enddef(ncid), 'enddef: '// trim(fout))
    call nf90_err(nf90_close(ncid), 'close: '// trim(fout))

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine setup_icerestart

  subroutine setup_ocnrestart(fin, fout, fgrid)

    character(len=*), intent(in) :: fin, fout, fgrid

    real(kind=8)              :: timestamp
    character(len= 40)        :: timeunit
    character(len= 20)        :: vname, vunit
    character(len=120)        :: vlong
    real(kind=8), allocatable :: Layer(:)          !< the vertical grid center

    integer :: ncid,varid,n,dims3(3),dims4(4)
    integer :: idimid,jdimid,kdimid,qidimid,qjdimid,timid
    real(kind=8), dimension(nxr,nyr) :: lonCt, latCt, lonCu, latCv

    call nf90_err(nf90_open(trim(fin), nf90_nowrite, ncid), 'open: '//trim(fin))
    ! get the time and layer information from the input restart file
    call nf90_err(nf90_inq_varid(ncid, 'Time', varid), 'get variable Id: Time '//trim(fin))
    call nf90_err(nf90_get_var(ncid, varid, timestamp), 'get variable: timestamp '//trim(fin))
    call nf90_err(nf90_get_att(ncid, varid, 'units', timeunit), 'get variable attribute : units '//trim(fin))
    ! layer
    allocate(Layer(nlevs)) ; Layer = 0.0
    call nf90_err(nf90_inq_varid(ncid, 'Layer', varid), 'get variable Id: Layer '//trim(fin))
    call nf90_err(nf90_get_var(ncid, varid, layer), 'get variable: Layer '//trim(fin))
    call nf90_err(nf90_close(ncid), 'close: '//trim(fin))

    call nf90_err(nf90_create(trim(fout), nf90_clobber, ncid), 'create: '//trim(fout))
    call nf90_err(nf90_def_dim(ncid, 'nx', nxr, idimid), 'define dimension: nx')
    call nf90_err(nf90_def_dim(ncid, 'ny', nyr, jdimid), 'define dimension: ny')
    call nf90_err(nf90_def_dim(ncid, 'Layer',  nlevs, kdimid), 'define dimension: Layer')
    call nf90_err(nf90_def_dim(ncid, 'Time', nf90_unlimited, timid), 'define dimension: Time')
    ! define the time variable
    call nf90_err(nf90_def_var(ncid, 'Time', nf90_double, (/timid/), varid), 'define variable: Time')
    call nf90_err(nf90_put_att(ncid, varid,    'units', trim(timeunit)), 'put variable attribute: units')
    ! vertical grid
    call nf90_err(nf90_def_var(ncid, 'Layer', nf90_double,  (/kdimid/), varid), 'define variable: Layer')
    call nf90_err(nf90_put_att(ncid, varid, 'units', 'm'), 'put variable attribute: units')

    if (allocated(b2d)) then
       do n = 1,nbilin2d
          vname = trim(b2d(n)%var_name)
          vunit = trim(b2d(n)%units)
          vlong = trim(b2d(n)%long_name)
          call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,timid/), varid), 'define variable: '// vname)
          call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units')
          call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name')
       enddo
    end if
!     if (allocated(c2d)) then
!        do n = 1,nconsd2d
!           vname = trim(c2d(n)%var_name)
!           vunit = trim(c2d(n)%units)
!           vlong = trim(c2d(n)%long_name)
!           if (trim(c2d(n)%var_grid) == 'Cu') dims3 = (/qidimid,jdimid,timid/)
!           if (trim(c2d(n)%var_grid) == 'Cv') dims3 = (/idimid,qjdimid,timid/)
!           if (trim(c2d(n)%var_grid) == 'Ct') dims3 = (/idimid,jdimid,timid/)
!           call nf90_err(nf90_def_var(ncid, vname, nf90_double, dims3, varid), 'define variable: '// vname)
!           call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units' )
!           call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name' )
!        enddo
!     end if
!
    if (allocated(b3d)) then
       do n = 1,nbilin3d
          vname = trim(b3d(n)%var_name)
          vunit = trim(b3d(n)%units)
          vlong = trim(b3d(n)%long_name)
          call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,kdimid,timid/), varid), 'define variable: '// vname)
          call nf90_err(nf90_put_att(ncid, varid,      'units', vunit), 'put variable attribute: units' )
          call nf90_err(nf90_put_att(ncid, varid,  'long_name', vlong), 'put variable attribute: long_name' )
       enddo
    end if
    call nf90_err(nf90_enddef(ncid), 'enddef: '// trim(fout))

    ! time
    call nf90_err(nf90_inq_varid(ncid, 'Time', varid), 'get variable Id: Time')
    call nf90_err(nf90_put_var(ncid, varid, timestamp), 'put variable: Time')
    ! vertical
    call nf90_err(nf90_inq_varid(ncid, 'Layer', varid), 'get variable Id: Layer')
    call nf90_err(nf90_put_var(ncid, varid, Layer)    , 'put variable: Layer')
    call nf90_err(nf90_close(ncid), 'close: '//trim(fout))

  end subroutine setup_ocnrestart
end module restarts_mod
