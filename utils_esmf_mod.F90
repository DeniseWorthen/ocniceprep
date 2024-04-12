module utils_esmf_mod

  use ESMF

  implicit none

  private

  type(ESMF_RouteHandle) :: rh
  type(ESMF_Mesh)        :: meshsrc, meshdst
  type(ESMF_Field)       :: fldsrc, flddst

  public createRH
  public remapRH
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

    call ESMF_FieldDestroy(fdsrc)
    call EMSF_FieldDestroy(flddst)

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
