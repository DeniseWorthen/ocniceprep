program ocniceprep

  use ESMF
  use netcdf
  use init_mod   ,     only : nxt, nyt, nlevs, nxr, nyr, outvars, readnml, readcsv
  use init_mod   ,     only : wgtsdir, griddir, ftype, fsrc, fdst, input_file, angvar
  use init_mod   ,     only : do_ocnprep, debug, logunit
  use arrays_mod ,     only : b2d, c2d, b3d, rgb2d, rgb3d, rgc2d, setup_packing
  use arrays_mod ,     only : nbilin2d, nbilin3d, nconsd2d, bilin2d, bilin3d, consd2d
  use arrays_mod ,     only : mask3d, maskspval
  use utils_mod  ,     only : getfield, packarrays, remap, dumpnc, nf90_err
  use utils_esmf_mod , only : createRH, remapRH, ChkErr, rotremap
  use restarts_mod ,   only : setup_icerestart, setup_ocnrestart

  implicit none

  type(ESMF_VM)      :: vm
  character(len=160) :: gridfile
  character(len=160) :: wgtsfile
  character(len=160) :: fout

  real(kind=8), allocatable, dimension(:) :: angsrc         !< the rotation angle at the Ct points for the src grid
  real(kind=8), allocatable, dimension(:) :: angdst         !< the rotation angle at the Ct points for the dst grid

  ! work arrays for output netcdf
  real(kind=8), allocatable, dimension(:,:)   :: out2d !< 2D destination grid output array
  real(kind=8), allocatable, dimension(:,:,:) :: out3d !< 3D destination grid output array

  character(len=120) :: meshfsrc, meshfdst

  integer           :: nvalid
  integer           :: n,nn,rc,ncid,varid
  integer           :: idx1,idx2
  character(len=20) :: vname
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
  call nf90_err(nf90_close(ncid), 'close: '//trim(gridfile))

  ! --------------------------------------------------------
  ! get the 3rd (vertical or ncat) dimension and the masking
  ! variable
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

  ! --------------------------------------------------------
  ! get the masking variable for ocean 3-d remapping
  ! --------------------------------------------------------

  if (do_ocnprep) then
     allocate(mask3d(nxt*nyt,nlevs)); mask3d = 0.0
     call getfield(trim(input_file), 'h', dims=(/nxt,nyt,nlevs/), field=mask3d)

     where(mask3d .le. 1.0e-3)mask3d = maskspval
     where(mask3d .ne. maskspval)mask3d = 1.0
     call dumpnc(trim(ftype)//'.'//trim(fdst)//'.mask3d.nc', 'mask3d', dims=(/nxt,nyt,nlevs/), field=mask3d)
  end if

  ! --------------------------------------------------------
  ! create packed arrays for mapping and remap packed arrays
  ! to the destination grid
  ! --------------------------------------------------------

  !in reorder branch
  ! (nbilin2d,nxt*nyt), (nbilin3d,nlevs,nxt*nyt)
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

  ! 2D conserv
  if (allocated(consd2d)) then
     call packarrays(trim(input_file), trim(wgtsdir)//fsrc(3:5)//'/', cos(angsrc), sin(angsrc),         &
          c2d, dims=(/nxt,nyt/), nflds=nconsd2d, fields=consd2d)
     rgc2d = 0.0
     call remapRH(src_field=consd2d, dst_field=rgc2d)

     if (debug) then
        write(logunit,'(a)')'remap 2D fields conserv with '//trim(wgtsfile)
        write(logunit,'(a)')'packed min/max values, mapped min/max values'
        do n = 1,nconsd2d
           write(logunit,'(i4,a10,3(a2,a6),4g14.4)')n,trim(c2d(n)%var_name),'  ',                       &
                trim(c2d(n)%var_grid),'  ',trim(c2d(n)%var_pair),'  ', trim(c2d(n)%var_pair_grid),      &
		minval(consd2d(:,n)), maxval(consd2d(:,n)), minval(rgc2d(:,n)), maxval(rgc2d(:,n))
        end do
        call dumpnc(trim(ftype)//'.'//trim(fsrc)//'.consd2d.nc', 'consd2d', dims=(/nxt,nyt/),           &
             nflds=nconsd2d, field=consd2d)
	call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgconsd2d.nc', 'rgconsd2d', dims=(/nxr,nyr/),       &
             nflds=nconsd2d, field=rgc2d)
     end if
  end if

  ! 3D bilin
  if (allocated(bilin3d))then
     call packarrays(trim(input_file), trim(wgtsdir)//fsrc(3:5)//'/', cos(angsrc), sin(angsrc),         &
          b3d, dims=(/nxt,nyt,nlevs/), nflds=nbilin3d, fields=bilin3d)
     rgb3d = 0.0
     !call remapRH(src_field=bilin3d, dst_field=rgb3d, srcdims=(/nxt*nyt,nlevs,nbilin3d/),dstdims=(/nxr*nyr,nlevs,nbilin3d/))
     !(nx*ny,nlevs,nfields)
     if (do_ocnprep) then
        do n = 1,nlevs
           call remapRH(src_field=bilin3d(:,n,:), dst_field=rgb3d(:,n,:),hmask=mask3d(:,n))
        end do
     else
        do n = 1,nlevs
           call remapRH(src_field=bilin3d(:,n,:), dst_field=rgb3d(:,n,:))
        end do
     end if

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
     call rotremap(trim(wgtsdir)//fdst(3:5)//'/', b2d, cos(angdst), sin(angdst), dims=(/nxr,nyr/),         &
          nflds=nbilin2d, fields=rgb2d)
     if (debug) then
        call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin2d.ij.nc', 'rgbilin2d', dims=(/nxr,nyr/),       &
             nflds=nbilin2d, field=rgb2d)
     end if
  end if
  if (allocated(consd2d)) then
     call rotremap(trim(wgtsdir)//fdst(3:5)//'/', c2d, cos(angdst), sin(angdst), dims=(/nxr,nyr/),         &
          nflds=nconsd2d, fields=rgc2d)
     if (debug) then
        call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin2d.ij.nc', 'rgbilin2d', dims=(/nxr,nyr/),       &
             nflds=nconsd2d, field=rgc2d)
     end if
  end if
  if (allocated(bilin3d)) then
     call rotremap(trim(wgtsdir)//fdst(3:5)//'/', b3d, cos(angdst), sin(angdst), dims=(/nxr,nyr,nlevs/),   &
          nflds=nbilin3d, fields=rgb3d)
     if (debug) then
        call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin3d.ij.nc', 'rgbilin3d', dims=(/nxr,nyr,nlevs/), &
             nk=nlevs, nflds=nbilin3d, field=rgb3d)
     end if
  end if

  ! --------------------------------------------------------
  ! write the mapped fields
  ! --------------------------------------------------------

  allocate(out2d(nxr,nyr)); out2d = 0.0
  allocate(out3d(nxr,nyr,nlevs)); out3d = 0.0

  fout = trim(ftype)//'.'//trim(fdst)//'.nc'
  if (debug) write(logunit, '(a)')'output file: '//trim(fout)

  gridfile = trim(griddir)//fdst(3:5)//'/'//'tripole.'//trim(fdst)//'.nc'
  if (do_ocnprep) then
     call setup_ocnrestart(trim(input_file),trim(fout),trim(gridfile))
  else
     call setup_icerestart(trim(input_file),trim(fout))
  end if

  call nf90_err(nf90_open(trim(fout), nf90_write, ncid), 'write: '//trim(fout))
  if (allocated(rgb2d)) then
     do n = 1,nbilin2d
        out2d(:,:) = reshape(rgb2d(:,n), (/nxr,nyr/))
        ! temp workaround
        if (b2d(n)%var_grid(1:2) == 'Bu') out2d(:,nyr) = out2d(:,nyr-1)
        vname = trim(b2d(n)%var_name)
        print *,b2d(n)%var_grid(1:2),vname
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
        ! temp workaround
        if (b3d(n)%var_grid(1:2) == 'Cv') out3d(:,nyr,:) = out3d(:,nyr-1,:)
        vname = trim(b3d(n)%var_name)
        call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable Id: '//vname)
        call nf90_err(nf90_put_var(ncid,   varid, out3d), 'put variable: '//vname)
     end do
  end if
  call nf90_err(nf90_close(ncid), 'close: '// trim(fout))
  write(logunit,'(a)')trim(fout)//' done'

  stop

end program ocniceprep
