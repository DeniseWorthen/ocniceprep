module utils_esmf_mod

  use ESMF
  use netcdf
  use init_mod  , only : debug, logunit, vardefs, fsrc, fdst, ftype, nxr, nyr
  use utils_mod , only : dumpnc, remap
  use arrays_mod, only : maskspval, mask3d

  implicit none

  private

  type(ESMF_RouteHandle) :: rh
  type(ESMF_DynamicMask) :: dynamicLevMask
  type(ESMF_Mesh)        :: meshsrc, meshdst
  type(ESMF_Field)       :: fldsrc, flddst

  interface remapRH
     module procedure remapRH2d
     module procedure remapRHdyn
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
    integer :: srcTermProcessing = 0

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
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,         &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                    &
         ignoreDegenerate=.true.,                              &
         srcTermProcessing=srcTermProcessing,                  &
         !dstStatusField=dststatusfield,                       &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create a dynamic mask object
    call ESMF_DynamicMaskSetR8R8R8(dynamicLevMask, dynamicSrcMaskValue=maskspval, &
         dynamicMaskRoutine=DynLevMaskProc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine createRH

  !----------------------------------------------------------
  ! remap a R8 packed field of nlen,nflds via ESMF RH
  !----------------------------------------------------------
  subroutine remapRH2d(src_field,dst_field)

    real(kind=8),    intent(in)           :: src_field(:,:)
    real(kind=8),    intent(out)          :: dst_field(:,:)

    integer               :: rc
    real(kind=8), pointer :: srcptr(:,:), dstptr(:,:)
    character(len=20)     :: subname = 'remapRH2d'

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
    call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dst_field = dstptr

  end subroutine remapRH2d

  !----------------------------------------------------------
  ! remap a R8 packed field of nlen,nflds via ESMF RH
  ! with dyanmic masking
  ! !(nx*ny,nlevs,nfields)
  !----------------------------------------------------------
  subroutine remapRHdyn(src_field,dst_field,hmask)

    !nx*ny,nfields
    real(kind=8), intent(in)  :: src_field(:,:)
    real(kind=8), intent(in)  :: hmask(:)
    real(kind=8), intent(out) :: dst_field(:,:)

    integer               :: n,rc
    real(kind=8), pointer :: srcptr(:,:), dstptr(:,:)
    character(len=20)     :: subname = 'remapRHdyn'

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

    srcptr = src_field
    do n = 1,ubound(src_field,2)
       where(hmask .eq. maskspval)srcptr(:,n) = maskspval
    end do
    dstptr = maskspval

    call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, dynamicMask=dynamicLevMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    dst_field = dstptr

  end subroutine remapRHdyn

  ! !----------------------------------------------------------
  ! ! remap a R8 packed field of nlen,nflds via ESMF RH with
  ! ! dynamic masking on nlev
  ! !----------------------------------------------------------
  ! subroutine remapRH3d(src_field,dst_field,hmask)

  !   real(kind=8), intent(in)  :: src_field(:,:)
  !   real(kind=8), intent(in)  :: hmask(:)
  !   real(kind=8), intent(out) :: dst_field(:,:)

  !   integer               :: rc
  !   real(kind=8), pointer :: srcptr(:,:), dstptr(:,:)
  !   character(len=20)     :: subname = 'remap3dRH'

  !   if (debug)write(logunit,'(a)')'enter '//trim(subname)

  !   fldsrc = ESMF_FieldCreate(meshsrc, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
  !        ungriddedLbound=(/1/), ungriddedUbound=(/size(src_field,2)/),       &
  !        gridToFieldMap=(/1/), rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   flddst = ESMF_FieldCreate(meshdst, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
  !        ungriddedLbound=(/1/), ungriddedUbound=(/size(dst_field,2)/),       &
  !        gridToFieldMap=(/1/), rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !   call ESMF_FieldGet(fldsrc, farrayptr=srcptr, rc=rc)
  !   if (chkerr(rc,__LINE__,u_FILE_u)) return
  !   call ESMF_FieldGet(flddst, farrayptr=dstptr, rc=rc)
  !   if (chkerr(rc,__LINE__,u_FILE_u)) return

  !   dstptr = 0.0
  !   srcptr = src_field
  !   call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   dst_field = dstptr



  !   ! fldsrc = ESMF_FieldCreate(meshsrc, farrayPtr=srcptr, meshloc=ESMF_MESHLOC_ELEMENT, &
  !   !      ungriddedLbound=(/1/), ungriddedUbound=(/size(src_field,2)/),       &
  !   !      gridToFieldMap=(/1/), rc=rc)
  !   ! if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  !   ! flddst = ESMF_FieldCreate(meshdst, farrayPtr=dstptr, meshloc=ESMF_MESHLOC_ELEMENT, &
  !   !      ungriddedLbound=(/1/), ungriddedUbound=(/size(dst_field,2)/),       &
  !   !      gridToFieldMap=(/1/), rc=rc)
  !   ! if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !   ! Do mapping using dynamic masking
  !   ! call ESMF_FieldRegrid(fldsrc, flddst, routehandle=rh, dynamicMask=dynamicLevMask, rc=rc)
  !   ! if (chkerr(rc,__LINE__,u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! end subroutine remapRH3d

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

    integer            :: n, idx1, idx2
    real(kind=8), allocatable, dimension(:) :: urot, vrot
    character(len=10)  :: vgrid1, vgrid2
    character(len=240) :: wgtsfile
    character(len=20)  :: subname = 'rotremap2d'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    idx1 = 0; idx2 = 0
    do n = 1,nflds
       if (len_trim(vars(n)%var_pair) > 0 .and. idx1 .eq. 0) then
          idx1 = n
          idx2 = n+1
       end if
    end do
    print *,idx1,idx2
    if (idx1 .eq. 0)return

    vgrid1 = vars(idx1)%var_grid(1:2)
    vgrid2 = vars(idx1)%var_pair_grid(1:2)

    allocate(urot(1:dims(1)*dims(2))); urot = 0.0
    allocate(vrot(1:dims(1)*dims(2))); vrot = 0.0
    urot(:) = fields(:,idx1)*cosrot(:) - fields(:,idx2)*sinrot(:)
    vrot(:) = fields(:,idx2)*cosrot(:) + fields(:,idx1)*sinrot(:)

    if (debug) write(logunit,'(a)')'restagger from Ct to '//trim(vgrid1)//' and '//trim(vgrid2)

    wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid1)//'.bilinear.nc'
    call remap(trim(wgtsfile), urot, fields(:,idx1))
    wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid2)//'.bilinear.nc'
    call remap(trim(wgtsfile), vrot, fields(:,idx2))

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine rotremap2d

  !----------------------------------------------------------
  ! rotate nlevs vectors from EW->IJ and map back to native staggers
  !----------------------------------------------------------
  subroutine rotremap3d(wdir, vars, cosrot, sinrot, dims, nflds, fields)

    character(len=*), intent(in)    :: wdir
    real(kind=8),     intent(in)    :: cosrot(:),sinrot(:)
    type(vardefs),    intent(in)    :: vars(:)
    integer,          intent(in)    :: dims(:)
    integer,          intent(in)    :: nflds
    real(kind=8),     intent(inout) :: fields(:,:,:)

    integer            :: k, n, idx1, idx2
    real(kind=8), allocatable, dimension(:) :: urot, vrot
    character(len=10)  :: vgrid1, vgrid2
    character(len=240) :: wgtsfile
    character(len=20)  :: subname = 'rotremap3d'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    !  rgbilin3d(nxr*nyr,nlevs,nbilin3d)

    idx1 = 0; idx2 = 0
    do n = 1,nflds
       if (len_trim(vars(n)%var_pair) > 0 .and. idx1 .eq. 0) then
          idx1 = n
          idx2 = n+1
       end if
    end do
    if (idx1 .eq. 0)return

    vgrid1 = vars(idx1)%var_grid(1:2)
    vgrid2 = vars(idx1)%var_pair_grid(1:2)
    if (debug) write(logunit,'(a)')'restagger from Ct to '//trim(vgrid1)//' and '//trim(vgrid2)

    allocate(urot(1:dims(1)*dims(2))); urot = 0.0
    allocate(vrot(1:dims(1)*dims(2))); vrot = 0.0
    do k = 1,dims(3)
       urot(:) = fields(:,k,idx1)*cosrot(:) - fields(:,k,idx2)*sinrot(:)
       vrot(:) = fields(:,k,idx2)*cosrot(:) + fields(:,k,idx1)*sinrot(:)
       wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid1)//'.bilinear.nc'
       call remap(trim(wgtsfile), urot, fields(:,k,idx1))
       wgtsfile = trim(wdir)//'tripole.'//trim(fdst)//'.Ct.to.'//trim(vgrid2)//'.bilinear.nc'
       call remap(trim(wgtsfile), vrot, fields(:,k,idx2))
    end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine rotremap3d

  !----------------------------------------------------------
  !
  !----------------------------------------------------------

  subroutine dynLevMaskProc(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    use ESMF, only : ESMF_RC_ARG_BAD

    ! input/output arguments
    type(ESMF_DynamicMaskElementR8R8R8) , pointer :: dynamicMaskList(:)
    real(ESMF_KIND_R8), intent(in), optional      :: dynamicSrcMaskValue
    real(ESMF_KIND_R8), intent(in), optional      :: dynamicDstMaskValue
    integer           , intent(out)               :: rc

    ! local variables
    integer  :: i, j
    real(ESMF_KIND_R8)  :: renorm
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (associated(dynamicMaskList)) then
       do i=1, size(dynamicMaskList)
          dynamicMaskList(i)%dstElement = 0.0d0 ! set to zero
          renorm = 0.d0 ! reset
          do j = 1, size(dynamicMaskList(i)%factor)
             if (dynamicSrcMaskValue /= dynamicMaskList(i)%srcElement(j)) then
                dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement + &
                     (dynamicMaskList(i)%factor(j) * dynamicMaskList(i)%srcElement(j))
                renorm = renorm + dynamicMaskList(i)%factor(j)
             endif
          enddo
          if (renorm > 0.d0) then
             dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement / renorm
          else if (present(dynamicSrcMaskValue)) then
             dynamicMaskList(i)%dstElement = dynamicSrcMaskValue
          else
             rc = ESMF_RC_ARG_BAD  ! error detected
             return
          endif
       enddo
    endif

  end subroutine DynLevMaskProc

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
