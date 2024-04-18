module utils_esmf_mod

  use ESMF
  use netcdf
  use init_mod  , only : debug, logunit, vardefs, fsrc, fdst, ftype, nxr, nyr
  use utils_mod , only : dumpnc, remap

  implicit none

  private

  type(ESMF_RouteHandle) :: rh
  type(ESMF_Mesh)        :: meshsrc, meshdst
  type(ESMF_Field)       :: fldsrc, flddst

  interface remapRH
     module procedure remapRH2d
     module procedure remapRH3d
  end interface remapRH

  interface rotremap
     module procedure rotremap2d
     module procedure rotremap3d
  end interface rotremap

  public createRH
  public remapRH
  public rotremap
  public ChkErr

  character(len=*), parameter :: u_FILE_u = &
       __FILE__
contains
  !----------------------------------------------------------
  ! create a RH
  !----------------------------------------------------------
  subroutine createRH(srcmeshfile, dstmeshfile,rc)

    character(len=*), intent(in)    :: srcmeshfile, dstmeshfile
    integer,          intent(inout) :: rc

    ! local variables
    real(kind=8) , pointer  :: srcptr(:), dstptr(:)

    meshsrc = ESMF_MeshCreate(filename=trim(srcmeshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    fldsrc = ESMF_FieldCreate(meshsrc, ESMF_TYPEKIND_R8, name='mshsrc', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    meshdst = ESMF_MeshCreate(filename=trim(dstmeshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    flddst = ESMF_FieldCreate(meshdst, ESMF_TYPEKIND_R8, name='mshdst', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=rh, &
         srcMaskValues=(/0/),                                  &
         dstMaskValues=(/0/),                                  &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR,              &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD,          &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                    &
         ignoreDegenerate=.true.,                              &
         !dstStatusField=dststatusfield,                       &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine createRH

  !----------------------------------------------------------
  ! remap a R8 packed field of nflds via ESMF RH
  !----------------------------------------------------------
  subroutine remapRH2d(src_field,dst_field)

    real(kind=8), intent(in)  :: src_field(:,:)
    real(kind=8), intent(out) :: dst_field(:,:)

    integer               :: rc
    real(kind=8), pointer :: srcptr(:,:), dstptr(:,:)
    character(len=20)     :: subname = 'remap2dRH'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    fldsrc = ESMF_FieldCreate(meshsrc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/size(src_field,2)/),       &
         gridToFieldMap=(/1/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    flddst = ESMF_FieldCreate(meshdst, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/1/), ungriddedUbound=(/size(dst_field,2)/),       &
         gridToFieldMap=(/1/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(fldsrc, farrayptr=srcptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(flddst, farrayptr=dstptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dstptr = 0.0
    srcptr = src_field

    print '(a,2g14.7)','src min/max ',minval(srcptr), maxval(srcptr)
    call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    print '(a,2g14.7)','dst min/max ',minval(dstptr), maxval(dstptr)
    dst_field = dstptr

  end subroutine remapRH2d

  !----------------------------------------------------------
  ! remap a R8 packed field of nlevs,nflds via ESMF RH
  !----------------------------------------------------------
  subroutine remapRH3d(src_field,dst_field)

    real(kind=8), intent(in)  :: src_field(:,:,:)
    real(kind=8), intent(out) :: dst_field(:,:,:)

    integer               :: rc
    real(kind=8), pointer :: srcptr(:,:,:), dstptr(:,:,:)
    character(len=20)     :: subname = 'remap3dRH'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    fldsrc = ESMF_FieldCreate(meshsrc, farrayPtr=srcptr, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    flddst = ESMF_FieldCreate(meshdst, farrayPtr=dstptr, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine remapRH3d

  !----------------------------------------------------------
  ! rotate vectors from EW->IJ and map back to native staggers
  !----------------------------------------------------------
  subroutine rotremap2d(wdir, vars, cosrot, sinrot, dims, nflds, fields)

    character(len=*), intent(in)    :: wdir
    real(kind=8),     intent(in)    :: cosrot(:),sinrot(:)
    type(vardefs),    intent(in)    :: vars(:)
    integer,          intent(in)    :: dims(:)
    integer,          intent(in)    :: nflds
    real(kind=8),     intent(inout) :: fields(:,:)

    integer           :: n, idx1, idx2
    real(kind=8)      :: urot, vrot
    character(len=10) :: vgrid1, vgrid2
    character(len=240) :: wgtsfile
    character(len=20) :: subname = 'rotremap2d'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    idx1 = 0; idx2 = 0
    do n = 1,nflds
       if (len_trim(vars(n)%var_pair) > 0 .and. idx1 .eq. 0) then
          idx1 = n
          idx2 = n+1
       end if
    end do

    do n = 1,dims(1)*dims(2)
       urot = fields(n,idx1)*cosrot(n) - fields(n,idx2)*sinrot(n)
       vrot = fields(n,idx2)*cosrot(n) + fields(n,idx1)*sinrot(n)
       fields(n,idx1) = urot
       fields(n,idx2) = vrot
    end do
    vgrid1 = vars(idx1)%var_grid(1:2)
    vgrid2 = vars(idx2)%var_grid(1:2)

    call dumpnc(trim(ftype)//'.'//trim(fdst)//'.rgbilin2d.Ct.ij.nc', 'rgbilin2d', dims=(/nxr,nyr/),       &
         nflds=nflds, field=fields)

    wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid1)//'.bilinear.nc'
    print '(a)','X0 '//trim(wgtsfile)
    call remap(trim(wgtsfile), fields(:,idx1), fields(:,idx1))
    print '(a)','X1 '//trim(wgtsfile)
    wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid2)//'.bilinear.nc'
    call remap(trim(wgtsfile), fields(:,idx2), fields(:,idx2))

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine rotremap2d

  !----------------------------------------------------------
  ! rotate nlevs vectors from EW->IJ and map back to native staggers
  !----------------------------------------------------------
  subroutine rotremap3d(vars, cosrot, sinrot, dims, nflds, fields)

    real(kind=8),     intent(in)    :: cosrot(:),sinrot(:)
    type(vardefs),    intent(in)    :: vars(:)
    integer,          intent(in)    :: dims(:)
    integer,          intent(in)    :: nflds
    real(kind=8),     intent(inout) :: fields(:,:,:)

    integer           :: n, idx1, idx2
    real(kind=8)      :: urot, vrot
    character(len=20) :: subname = 'rotremap2d'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    ! idx1 = 0; idx2 = 0
    ! do n = 1,nflds
    !    if (len_trim(vars(n)%var_pair) > 0 .and. idx1 .eq. 0) then
    !       idx1 = n
    !       idx2 = n+1
    !    end if
    ! end do
    ! do n = 1,dims(1)*dims(2)
    !    urot = fields(n,idx1)*cosrot(n) - fields(n,idx2)*sinrot(n)
    !    vrot = fields(n,idx2)*cosrot(n) + fields(n,idx1)*sinrot(n)
    !    fields(n,idx1) = urot
    !    fields(n,idx2) = vrot
    ! end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine rotremap3d

  !----------------------------------------------------------
  ! handle ESMF errors
  !----------------------------------------------------------
  logical function ChkErr(rc, line, file)
    integer, intent(in) :: rc            !< return code to check
    integer, intent(in) :: line          !< Integer source line number
    character(len=*), intent(in) :: file !< User-provided source file name
    integer :: lrc
    ChkErr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErr = .true.
    endif
  end function ChkErr
end module utils_esmf_mod
