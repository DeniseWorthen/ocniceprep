module utils_mod

  use netcdf
  use init_mod, only : debug, logunit, vardefs, fsrc

  implicit none

  private

  interface getfield
     module procedure getfield2d_R4
     module procedure getfield3d_R4
     module procedure getfield2d_R8
     module procedure getfield3d_R8
  end interface getfield

  interface packarrays
     module procedure packarrays2d_R4
     module procedure packarrays3d_R4
     module procedure packarrays2d_R8
     module procedure packarrays3d_R8
  end interface packarrays

  interface getvecpair
     module procedure getvecpair2d_R4
     module procedure getvecpair3d_R4
     module procedure getvecpair2d_R8
     module procedure getvecpair3d_R8
  end interface getvecpair

  interface remap
     module procedure remap1d_R4
     module procedure remap2d_R4
     module procedure remap3d_R4
     module procedure remap1d_R8
     module procedure remap2d_R8
     module procedure remap3d_R8
  end interface remap

  interface remapRH
     module procedure remapRH2d
     module procedure remapRH3d
  end interface remapRH

  interface dumpnc
     module procedure dumpnc1d_R4
     module procedure dumpnc2d_R4
     module procedure dumpnc3d_R4
     module procedure dumpnc3dk_R4
     module procedure dumpnc1d_R8
     module procedure dumpnc2d_R8
     module procedure dumpnc3d_R8
     module procedure dumpnc3dk_R8
  end interface dumpnc

  public getfield
  public packarrays
  public remap
  public dumpnc
  public nf90_err

contains
  !----------------------------------------------------------
  ! pack 2D R8 fields into arrays by mapping type
  !----------------------------------------------------------
  subroutine packarrays2d_R8(filesrc, wgtsdir, cosrot, sinrot, vars, dims, nflds, fields)

    character(len=*), intent(in)  :: filesrc,wgtsdir
    real(kind=8),     intent(in)  :: cosrot(:),sinrot(:)
    type(vardefs),    intent(in)  :: vars(:)
    integer,          intent(in)  :: dims(:)
    integer,          intent(in)  :: nflds
    real(kind=8),     intent(out) :: fields(:,:)

    ! local variables
    integer                   :: n, nn
    real(kind=8), allocatable :: vecpair(:,:)
    character(len=20)         :: subname = 'packarrays2d_R8'

    fields=0.0
    if (debug)write(logunit,'(a)')'enter '//trim(subname)
    ! obtain vector pairs
    do n = 1,nflds
       if (trim(vars(n)%var_grid) == 'Cu' .or. trim(vars(n)%var_grid) == 'Bu_x') then
          allocate(vecpair(dims(1)*dims(2),2)); vecpair = 0.0
          call getvecpair(trim(filesrc), trim(wgtsdir), cosrot, sinrot,   &
               trim(vars(n)%var_name), trim(vars(n)%var_grid(1:2)), &
               trim(vars(n)%var_pair), trim(vars(n)%var_pair_grid(1:2)),  &
               dims=(/dims(1),dims(2)/), vecpair=vecpair)
       end if
    end do

    ! create packed array
    nn = 0
    do n = 1,nflds
       if (len_trim(vars(n)%var_pair) == 0) then
          nn = nn + 1
          call getfield(trim(filesrc), trim(vars(n)%var_name), dims=(/dims(1),dims(2)/), &
               field=fields(:,nn))
       else ! fill with vector pairs
          nn = nn+1
          ! ocn vectors
          if (trim(vars(n)%var_grid) == 'Cu')fields(:,nn) = vecpair(:,1)
          if (trim(vars(n)%var_grid) == 'Cv')fields(:,nn) = vecpair(:,2)
          ! ice vectors
          if (trim(vars(n)%var_grid) == 'Bu_x')fields(:,nn) = vecpair(:,1)
          if (trim(vars(n)%var_grid) == 'Bu_y')fields(:,nn) = vecpair(:,2)
       end if
    end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine packarrays2d_R8

  !----------------------------------------------------------
  ! pack 3D R8 fields into arrays by mapping type
  !----------------------------------------------------------
  subroutine packarrays3d_R8(filesrc, wgtsdir, cosrot, sinrot, vars, dims, nflds, fields)

    character(len=*), intent(in)  :: filesrc,wgtsdir
    real(kind=8),     intent(in)  :: cosrot(:),sinrot(:)
    type(vardefs),    intent(in)  :: vars(:)
    integer,          intent(in)  :: dims(:)
    integer,          intent(in)  :: nflds
    real(kind=8),     intent(out) :: fields(:,:,:)

    ! local variables
    integer                   :: n, nn
    real(kind=8), allocatable :: vecpair(:,:,:)
    character(len=20)         :: subname = 'packarrays3d_R8'

    fields=0.0

    if (debug)write(logunit,'(a)')'enter '//trim(subname)
    ! obtain vector pairs
    do n = 1,dims(3)
       if (trim(vars(n)%var_grid) == 'Cu') then
          allocate(vecpair(dims(1)*dims(2),dims(3),2)); vecpair = 0.0
          call getvecpair(trim(filesrc), trim(wgtsdir), cosrot, sinrot, &
               trim(vars(n)%var_name), trim(vars(n)%var_grid),    &
               trim(vars(n)%var_pair), trim(vars(n)%var_pair_grid),     &
               dims=(/dims(1),dims(2),dims(3)/), vecpair=vecpair)
       end if
    end do

    ! create packed array
    nn = 0
    do n = 1,nflds
       if (len_trim(vars(n)%var_pair) == 0) then
          nn = nn + 1
          call getfield(trim(filesrc), trim(vars(n)%var_name), dims=(/dims(1),dims(2),dims(3)/), &
               field=fields(:,:,nn))
       else ! fill with vector pairs
          nn = nn+1
          if (trim(vars(n)%var_grid) == 'Cu')fields(:,:,nn) = vecpair(:,:,1)
          if (trim(vars(n)%var_grid) == 'Cv')fields(:,:,nn) = vecpair(:,:,2)
       end if
    end do

    if (debug)write(logunit,'(a)')'exit '//trim(subname)
  end subroutine packarrays3d_R8

  !----------------------------------------------------------
  ! obtain 2D R8 vector pairs mapped to Ct and rotated to EW
  !----------------------------------------------------------
  subroutine getvecpair2d_R8(fname, wdir, cosrot, sinrot, vname1, vgrid1, &
       vname2, vgrid2, dims, vecpair)

    character(len=*), intent(in)  :: fname
    character(len=*), intent(in)  :: wdir
    real(kind=8),     intent(in)  :: cosrot(:), sinrot(:)
    character(len=*), intent(in)  :: vname1, vgrid1, vname2, vgrid2
    integer,          intent(in)  :: dims(:)
    real(kind=8),     intent(out) :: vecpair(:,:)

    ! local variables
    integer :: ii
    real(kind=8), dimension(dims(1)*dims(2)) :: urot, vrot
    character(len=240) :: wgtsfile
    character(len=20) :: subname = 'getvecpair2d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    wgtsfile = trim(wdir)//'tripole.'//trim(fsrc)//'.'//vgrid1//'.to.Ct.bilinear.nc'
    call getfield(fname, vname1, dims=dims, field=vecpair(:,1), wgts=trim(wgtsfile))
    if (debug)write(logunit,'(a)')'wgtsfile for 2d vector '//trim(vname1)//'   '//trim(wgtsfile)
    wgtsfile = trim(wdir)//'tripole.'//trim(fsrc)//'.'//vgrid2//'.to.Ct.bilinear.nc'
    call getfield(fname, vname2, dims=dims, field=vecpair(:,2), wgts=trim(wgtsfile))
    if (debug)write(logunit,'(a)')'wgtsfile for 2d vector '//trim(vname2)//'   '//trim(wgtsfile)

    urot = 0.0; vrot = 0.0
    do ii = 1,dims(1)*dims(2)
       urot(ii) = vecpair(ii,1)*cosrot(ii) + vecpair(ii,2)*sinrot(ii)
       vrot(ii) = vecpair(ii,2)*cosrot(ii) - vecpair(ii,1)*sinrot(ii)
    end do
    vecpair(:,1) = urot(:)
    vecpair(:,2) = vrot(:)

    if (debug) write(logunit,'(a)')'exit '//trim(subname)
  end subroutine getvecpair2d_R8

  !----------------------------------------------------------
  ! obtain 3D R8 vector pairs, mapped to Ct and rotated to EW
  !----------------------------------------------------------
  subroutine getvecpair3d_R8(fname, wdir, cosrot, sinrot, vname1, vgrid1, &
       vname2, vgrid2, dims, vecpair)

    character(len=*), intent(in)  :: fname
    character(len=*), intent(in)  :: wdir
    real(kind=8),     intent(in)  :: cosrot(:), sinrot(:)
    character(len=*), intent(in)  :: vname1, vgrid1, vname2, vgrid2
    integer,          intent(in)  :: dims(:)
    real(kind=8),     intent(out) :: vecpair(:,:,:)

    ! local variables
    integer :: ii,k
    real(kind=8), dimension(dims(1)*dims(2)) :: urot, vrot
    character(len=240) :: wgtsfile
    character(len=20)  :: subname = 'getvecpair3d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    wgtsfile = trim(wdir)//'tripole.'//trim(fsrc)//'.'//vgrid1//'.to.Ct.bilinear.nc'
    call getfield(fname, vname1, dims=dims, field=vecpair(:,:,1), wgts=trim(wgtsfile))
    wgtsfile = trim(wdir)//'tripole.'//trim(fsrc)//'.'//vgrid2//'.to.Ct.bilinear.nc'
    call getfield(fname, vname2, dims=dims, field=vecpair(:,:,2), wgts=trim(wgtsfile))

    do k = 1,dims(3)
       urot = 0.0; vrot = 0.0
       do ii = 1,dims(1)*dims(2)
          urot(ii)= vecpair(ii,k,1)*cosrot(ii) + vecpair(ii,k,2)*sinrot(ii)
          vrot(ii)= vecpair(ii,k,2)*cosrot(ii) - vecpair(ii,k,1)*sinrot(ii)
       end do
       vecpair(:,k,1) = urot(:)
       vecpair(:,k,2) = vrot(:)
    end do

    if (debug) write(logunit,'(a)')'exit '//trim(subname)
  end subroutine getvecpair3d_R8

  !----------------------------------------------------------
  ! obtain a R8 2D field and return a R8 1-D vector array
  !----------------------------------------------------------
  subroutine getfield2d_R8(fname, vname, dims, field, wgts, chkfill)

    character(len=*),           intent(in)  :: fname, vname
    integer,                    intent(in)  :: dims(:)
    real(kind=8),               intent(out) :: field(:)
    character(len=*), optional, intent(in)  :: wgts
    logical,          optional, intent(in)  :: chkfill

    ! local variable
    logical                   :: lchk
    integer                   :: ncid, varid, rc
    real(kind=8)              :: fval
    real(kind=8), allocatable :: a2d(:,:)
    real(kind=8), allocatable :: atmp(:)
    character(len=20)         :: subname = 'getfield2d_R4'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname

    allocate(a2d(dims(1),dims(2))); a2d = 0.0
    allocate(atmp(dims(1)*dims(2))); atmp = 0.0

    lchk = .false.
    if (present(chkfill)) then
       if (chkfill) lchk = chkfill
    end if

    call nf90_err(nf90_open(fname, nf90_nowrite, ncid), 'nf90_open: '//fname)
    call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable ID: '//vname)
    call nf90_err(nf90_get_var(ncid, varid, a2d), 'get variable: '//vname)
    if (lchk) then
       call nf90_err(nf90_get_att(ncid, varid, '_FillValue', fval), 'get attribute FillValue: '//vname)
    end if
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    atmp(:) = reshape(a2d, (/dims(1)*dims(2)/))
    if (lchk) then
       where(atmp .eq. fval)atmp = 0.0
    end if
    if(present(wgts)) then
       call remap(trim(wgts), src_field=atmp, dst_field=field)
    else
       field = atmp
    end if

    if (debug) write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname
  end subroutine getfield2d_R8

  !----------------------------------------------------------
  ! obtain a R8 3D field and return a R8 2-D vector array
  !----------------------------------------------------------
  subroutine getfield3d_R8(fname, vname, dims, field, wgts, chkfill)

    character(len=*),           intent(in)  :: fname, vname
    integer,                    intent(in)  :: dims(:)
    real(kind=8),               intent(out) :: field(:,:)
    character(len=*), optional, intent(in)  :: wgts
    logical,          optional, intent(in)  :: chkfill

    ! local variable
    logical                   :: lchk
    integer                   :: ncid, varid, rc
    real(kind=8)              :: fval
    real(kind=8), allocatable :: a3d(:,:,:)
    real(kind=8), allocatable :: atmp(:,:)
    character(len=20)         :: subname = 'getfield3d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname

    allocate(a3d(dims(1),dims(2),dims(3))); a3d = 0.0
    allocate(atmp(dims(1)*dims(2),dims(3))); atmp = 0.0

    lchk = .false.
    if (present(chkfill)) then
       if (chkfill) lchk = chkfill
    end if

    call nf90_err(nf90_open(fname, nf90_nowrite, ncid), 'nf90_open: '//fname)
    call nf90_err(nf90_inq_varid(ncid, vname, varid), 'get variable ID: '//vname)
    call nf90_err(nf90_get_var(ncid, varid, a3d), 'get variable: '//vname)
    if (lchk) then
       call nf90_err(nf90_get_att(ncid, varid, '_FillValue', fval), 'get attribute FillValue: '//vname)
    end if
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    atmp(:,:) = reshape(a3d, (/dims(1)*dims(2),dims(3)/))
    if (lchk) then
       where(atmp .eq. fval)atmp = 0.0
    end if
    where(atmp .eq. fval)atmp = 0.0
    if(present(wgts)) then
       call remap(trim(wgts), dim2=dims(3), src_field=atmp, dst_field=field)
    else
       field = atmp
    end if

    if (debug) write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname
  end subroutine getfield3d_R8

  !----------------------------------------------------------
  ! remap a R8 1-D vector array
  !----------------------------------------------------------
  subroutine remap1d_R8(fname, src_field, dst_field)

    character(len=*), intent(in)  :: fname
    real(kind=8),     intent(in)  :: src_field(:)
    real(kind=8),     intent(out) :: dst_field(:)

    ! local variables
    integer :: ncid, rc, id
    integer :: i,ii,jj
    integer :: n_a, n_b, n_s
    integer(kind=4), allocatable, dimension(:) :: col, row
    real(kind=8),    allocatable, dimension(:) :: S
    character(len=20) :: subname = 'remap1d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)

    ! retrieve the weights
    call nf90_err(nf90_open(trim(fname), nf90_nowrite, ncid), 'open: '//fname)
    call nf90_err(nf90_inq_dimid(ncid, 'n_s', id), 'get dimension Id: n_s')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_s), 'get dimension: n_s' )
    call nf90_err(nf90_inq_dimid(ncid, 'n_a', id), 'get dimension Id: n_a')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_a), 'get dimension: n_a' )
    call nf90_err(nf90_inq_dimid(ncid, 'n_b', id), 'get dimension Id: n_b')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_b), 'get dimension: n_b' )

    allocate(col(1:n_s)); col = 0
    allocate(row(1:n_s)); row = 0
    allocate(  S(1:n_s)); S = 0.0

    call nf90_err(nf90_inq_varid(ncid, 'col', id),'get variable Id: col')
    call nf90_err(nf90_get_var(ncid,     id, col),'get variable: col')
    call nf90_err(nf90_inq_varid(ncid, 'row', id),'get variable Id: row')
    call nf90_err(nf90_get_var(ncid,     id, row),'get variable: row')
    call nf90_err(nf90_inq_varid(ncid,   'S', id),'get variable Id: S')
    call nf90_err(nf90_get_var(ncid,      id,  S),'get variable: S')
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    dst_field = 0.0
    do i = 1,n_s
       ii = row(i); jj = col(i)
       dst_field(ii) = dst_field(ii) + S(i)*src_field(jj)
    enddo

    if (debug) write(logunit,'(a)')'exit '//trim(subname)
  end subroutine remap1d_R8

  !----------------------------------------------------------
  ! remap a R8 packed field of either nflds or nlevs
  !----------------------------------------------------------
  subroutine remap2d_R8(fname, dim2, src_field, dst_field)

    character(len=*), intent(in)  :: fname
    integer,          intent(in)  :: dim2
    real(kind=8),     intent(in)  :: src_field(:,:)
    real(kind=8),     intent(out) :: dst_field(:,:)

    ! local variables
    integer :: ncid, rc, id
    integer :: i,ii,jj
    integer :: n_a, n_b, n_s
    integer(kind=4), allocatable, dimension(:) :: col, row
    real(kind=8),    allocatable, dimension(:) :: S
    character(len=20) :: subname = 'remap2d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' weights = '//trim(fname)

    ! retrieve the weights
    call nf90_err(nf90_open(trim(fname), nf90_nowrite, ncid), 'open: '//fname)
    call nf90_err(nf90_inq_dimid(ncid, 'n_s', id), 'get dimension Id: n_s')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_s), 'get dimension: n_s')
    call nf90_err(nf90_inq_dimid(ncid, 'n_a', id), 'get dimension Id: n_a')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_a), 'get dimension: n_a')
    call nf90_err(nf90_inq_dimid(ncid, 'n_b', id), 'get dimension Id: n_b')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_b), 'get dimension: n_b')

    allocate(col(1:n_s)); col = 0
    allocate(row(1:n_s)); row = 0
    allocate(  S(1:n_s)); S = 0.0

    call nf90_err(nf90_inq_varid(ncid, 'col', id),'get variable Id: col')
    call nf90_err(nf90_get_var(ncid,     id, col),'get variable: col')
    call nf90_err(nf90_inq_varid(ncid, 'row', id),'get variable Id: row')
    call nf90_err(nf90_get_var(ncid,     id, row),'get variable: row')
    call nf90_err(nf90_inq_varid(ncid,   'S', id),'get variable Id: S')
    call nf90_err(nf90_get_var(ncid,      id,  S),'get variable: S')
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    dst_field = 0.0
    do i = 1,n_s
       ii = row(i); jj = col(i)
       dst_field(ii,:) = dst_field(ii,:) + S(i)*src_field(jj,:)
    enddo

    if (debug) write(logunit,'(a)')'exit '//trim(subname)
  end subroutine remap2d_R8

  !----------------------------------------------------------
  ! remap a R8 field packed array of nk levels and nflds fields
  !----------------------------------------------------------
  subroutine remap3d_R8(fname, nk, nflds, src_field, dst_field)

    character(len=*), intent(in)  :: fname
    integer,          intent(in)  :: nk, nflds
    real(kind=8),     intent(in)  :: src_field(:,:,:)
    real(kind=8),     intent(out) :: dst_field(:,:,:)

    ! local variables
    integer :: ncid, rc, id
    integer :: i,ii,jj
    integer :: n_a, n_b, n_s
    integer(kind=4), allocatable, dimension(:) :: col, row
    real(kind=8),    allocatable, dimension(:) :: S
    character(len=20) :: subname = 'remap3d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' weights = '//trim(fname)

    ! retrieve the weights
    call nf90_err(nf90_open(trim(fname), nf90_nowrite, ncid), 'open: '//fname)
    call nf90_err(nf90_inq_dimid(ncid, 'n_s', id), 'get dimension Id: n_s')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_s), 'get dimension: n_s')
    call nf90_err(nf90_inq_dimid(ncid, 'n_a', id), 'get dimension Id: n_a')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_a), 'get dimension: n_a')
    call nf90_err(nf90_inq_dimid(ncid, 'n_b', id), 'get dimension Id: n_b')
    call nf90_err(nf90_inquire_dimension(ncid, id, len=n_b), 'get dimension: n_b')

    allocate(col(1:n_s)); col = 0
    allocate(row(1:n_s)); row = 0
    allocate(  S(1:n_s)); S = 0.0

    call nf90_err(nf90_inq_varid(ncid, 'col', id),'get variable Id: col')
    call nf90_err(nf90_get_var(ncid,     id, col),'get variable: col')
    call nf90_err(nf90_inq_varid(ncid, 'row', id),'get variable Id: row')
    call nf90_err(nf90_get_var(ncid,     id, row),'get variable: row')
    call nf90_err(nf90_inq_varid(ncid,   'S', id),'get variable Id: S')
    call nf90_err(nf90_get_var(ncid,      id,  S),'get variable: S')
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    dst_field = 0.0
    do i = 1,n_s
       ii = row(i); jj = col(i)
       dst_field(ii,:,:) = dst_field(ii,:,:) + S(i)*src_field(jj,:,:)
    enddo

    if (debug) write(logunit,'(a)')'exit '//trim(subname)
  end subroutine remap3d_R8

  !----------------------------------------------------------
  ! write a bare netcdf file of a 2D packed R8 field
  !----------------------------------------------------------
  subroutine dumpnc2d_R8(fname, vname, dims, nflds, field)

    character(len=*), intent(in) :: fname, vname
    integer,          intent(in) :: dims(:)
    integer,          intent(in) :: nflds
    real(kind=8),     intent(in) :: field(:,:)

    ! local variable
    integer                   :: ncid, varid, rc, idimid, jdimid, fdimid
    real(kind=8), allocatable :: a3d(:,:,:)
    character(len=20)         :: subname = 'dumpnc2d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname
    allocate(a3d(dims(1),dims(2),nflds)); a3d = 0.0

    call nf90_err(nf90_create(trim(fname), nf90_clobber, ncid), 'create: '//fname)
    call nf90_err(nf90_def_dim(ncid, 'nx', dims(1), idimid), 'define dimension: nx')
    call nf90_err(nf90_def_dim(ncid, 'ny', dims(2), jdimid), 'define dimension: ny')
    call nf90_err(nf90_def_dim(ncid, 'nf', nflds,   fdimid), 'define dimension: nf')
    call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,fdimid/), varid), 'define variable: '//vname)
    call nf90_err(nf90_enddef(ncid), 'nf90_enddef: '//fname)

    a3d(:,:,:) =  reshape(field(1:dims(1)*dims(2),1:nflds), (/dims(1),dims(2),nflds/))
    call nf90_err(nf90_put_var(ncid, varid, a3d), 'put variable: '//vname)
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    if (debug)write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname
  end subroutine dumpnc2d_R8

  !----------------------------------------------------------
  ! write a bare netcdf file of a packed 3D R8 field
  !----------------------------------------------------------
  subroutine dumpnc3d_R8(fname, vname, dims, nk, nflds, field)

    character(len=*), intent(in) :: fname, vname
    integer,          intent(in) :: dims(:)
    integer,          intent(in) :: nk, nflds
    real(kind=8),     intent(in) :: field(:,:,:)

    ! local variable
    integer :: n, ncid, varid, rc, idimid, jdimid, kdimid, fdimid
    real(kind=8), allocatable :: a4d(:,:,:,:)
    character(len=20) :: subname = 'dumpnc3d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname
    allocate(a4d(dims(1),dims(2),dims(3),nflds)); a4d = 0.0

    call nf90_err(nf90_create(trim(fname), nf90_clobber, ncid), 'create: '//fname)
    call nf90_err(nf90_def_dim(ncid, 'nx', dims(1), idimid), 'define dimension: nx')
    call nf90_err(nf90_def_dim(ncid, 'ny', dims(2), jdimid), 'define dimension: ny')
    call nf90_err(nf90_def_dim(ncid, 'nk', dims(3), kdimid), 'define dimension: nk')
    call nf90_err(nf90_def_dim(ncid, 'nf', nflds,   fdimid), 'define dimension: nf')
    call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,kdimid,fdimid/), varid), 'define variable: '//vname)
    call nf90_err(nf90_enddef(ncid), 'nf90_enddef: '//fname)

    do n = 1,nflds
       a4d(:,:,:,n) = reshape(field(1:dims(1)*dims(2),1:dims(3),n), (/dims(1),dims(2),dims(3)/))
    end do
    call nf90_err(nf90_put_var(ncid, varid, a4d), 'put variable: '//vname)
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    if (debug)write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname
  end subroutine dumpnc3d_R8

  !----------------------------------------------------------
  ! write a bare netcdf file of an unpacked 3D R8 field
  !----------------------------------------------------------
  subroutine dumpnc3dk_R8(fname, vname, dims, field)

    character(len=*), intent(in) :: fname, vname
    integer,          intent(in) :: dims(:)
    real(kind=8),     intent(in) :: field(:,:)

    ! local variable
    integer                   :: ncid, varid, rc, idimid, jdimid, kdimid
    real(kind=8), allocatable :: a3d(:,:,:)
    character(len=20)         :: subname = 'dumpnc3dk_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname
    allocate(a3d(dims(1),dims(2),dims(3))); a3d = 0.0

    call nf90_err(nf90_create(trim(fname), nf90_clobber, ncid), 'nf90_create: '//fname)
    call nf90_err(nf90_def_dim(ncid, 'nx', dims(1), idimid), 'define dimension: nx')
    call nf90_err(nf90_def_dim(ncid, 'ny', dims(2), jdimid), 'define dimension: ny')
    call nf90_err(nf90_def_dim(ncid, 'nk', dims(3), kdimid), 'define dimension: nk')
    call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid,kdimid/), varid), 'define variable: '//vname)
    call nf90_err(nf90_enddef(ncid), 'nf90_enddef: '//fname)

    a3d(:,:,:) =  reshape(field(1:dims(1)*dims(2),1:dims(3)), (/dims(1),dims(2),dims(3)/))
    call nf90_err(nf90_put_var(ncid, varid, a3d), 'put variable: '//vname)
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    if (debug)write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname

  end subroutine dumpnc3dk_R8

  !----------------------------------------------------------
  ! write a bare netcdf file of an unpacked 2D R8 field
  !----------------------------------------------------------
  subroutine dumpnc1d_R8(fname, vname, dims, field)

    character(len=*), intent(in) :: fname, vname
    integer,          intent(in) :: dims(:)
    real(kind=8),     intent(in) :: field(:)

    ! local variable
    integer                   :: ncid, varid, rc, idimid, jdimid
    real(kind=8), allocatable :: a2d(:,:)
    character(len=20)         :: subname = 'dumpnc1d_R8'

    if (debug)write(logunit,'(a)')'enter '//trim(subname)//' variable '//vname
    allocate(a2d(dims(1),dims(2))); a2d = 0.0

    call nf90_err(nf90_create(trim(fname), nf90_clobber, ncid), 'nf90_create: '//fname)
    call nf90_err(nf90_def_dim(ncid, 'nx', dims(1), idimid), 'define dimension: nx')
    call nf90_err(nf90_def_dim(ncid, 'ny', dims(2), jdimid), 'define dimension: ny')
    call nf90_err(nf90_def_var(ncid, vname, nf90_double, (/idimid,jdimid/), varid), 'define variable: '//vname)
    call nf90_err(nf90_enddef(ncid), 'nf90_enddef: '//fname)

    a2d(:,:) =  reshape(field(1:dims(1)*dims(2)), (/dims(1),dims(2)/))
    call nf90_err(nf90_put_var(ncid, varid, a2d), 'put variable: '//vname)
    call nf90_err(nf90_close(ncid), 'close: '//fname)

    if (debug)write(logunit,'(a)')'exit '//trim(subname)//' variable '//vname

  end subroutine dumpnc1d_R8

  !----------------------------------------------------------
  ! handle netcdf errors
  !----------------------------------------------------------
  subroutine nf90_err(ierr, string)

    integer ,         intent(in) :: ierr
    character(len=*), intent(in) :: string
    if (ierr /= nf90_noerr) then
      write(0, '(a)') 'FATAL ERROR: ' // trim(string)// ' : ' // trim(nf90_strerror(ierr))
      ! This fails on WCOSS2 with Intel 19 compiler. See
      ! https://community.intel.com/t5/Intel-Fortran-Compiler/STOP-and-ERROR-STOP-with-variable-stop-codes/m-p/1182521#M149254
      ! When WCOSS2 moves to Intel 2020+, uncomment the next line and remove stop 99
      !stop ierr
      stop 99
    end if
  end subroutine nf90_err
end module utils_mod
