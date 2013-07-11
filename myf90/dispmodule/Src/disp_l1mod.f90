MODULE DISP_L1MOD

  ! Add-on module to DISPMODULE to display 1-byte logical items
  ! (assuming that these have kind = 1)
  !
  ! This module is obtained by copying the section DEFAULT LOGICAL PROCEDURES from
  ! dispmodule.f90, replacing dlog with log1 and 'default logical' with '1-byte
  ! logical' (only appears in comments), and adding the DECLARATIONS section below.
  !
  ! Copyright (c) 2008, Kristján Jónasson, Dept. of Computer Science, University of
  ! Iceland (jonasson@hi.is). This software is free. For details see the file README.
  
  ! ******************************** DECLARATIONS ********************************************

  use dispmodule_util

  PUBLIC DISP
  PUBLIC TOSTRING

  PRIVATE

  interface disp
    module procedure disp_s_log1, disp_ts_log1, disp_v_log1, disp_tv_log1, disp_m_log1, disp_tm_log1
  end interface

  interface tostring
    module procedure tostring_log1, tostring_f_log1, tostring_s_log1, tostring_sf_log1
  end interface

  integer, parameter :: log1 = 1  ! hopefully logical(1) is byte

CONTAINS

  ! ********************************************** 1-BYTE LOGICAL PROCEDURES *************************************************
  subroutine disp_s_log1(x, fmt, advance, sep, trim, unit)
    ! 1-byte logical scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim
    logical(log1), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_ts_log1('', x, fmt, advance, sep, 'left', trim, unit)
  end subroutine disp_s_log1

  subroutine disp_v_log1(x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! 1-byte logical vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    logical(log1), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tv_log1('', x, fmt, advance, lbound, sep, style, trim, unit, orient)
  end subroutine disp_v_log1

  subroutine disp_m_log1(x, fmt, advance, lbound, sep, style, trim, unit)
    ! 1-byte logical matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim
    logical(log1), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tm_log1('', x, fmt, advance, lbound, sep, style, trim, unit)
  end subroutine disp_m_log1

  subroutine disp_ts_log1(title, x, fmt, advance, sep, style, trim, unit)
    ! 1-byte logical scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim
    logical(log1), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_tm_log1(title, reshape((/x/), (/1, 1/)), fmt, advance, sep=sep, style=style, trim=trim, unit=unit)
  end subroutine disp_ts_log1

  subroutine disp_tv_log1(title, x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! 1-byte logical vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    logical(log1), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient)
    if (SE % row) then
      call disp_log1(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_log1(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_log1

  subroutine disp_tm_log1(title, x, fmt, advance, lbound, sep, style, trim, unit)
    ! 1-byte logical matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    logical(log1),intent(in)           :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Format edit descriptor to use for each matrix element (e.g. 'L1')
    integer,      intent(in), optional :: unit       ! Unit to display on
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no trimming, 
    !                                                ! 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    type(settings) :: SE
    !
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit)
    call disp_log1(title, x, SE)
  end subroutine disp_tm_log1

  subroutine disp_log1(title, x, SE)
    ! Write 1-byte logical to box or unit
    character(*),   intent(in)    :: title
    logical(log1),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer wid(size(x,2)), nbl(size(x,2))
    if (SE % w <= 0 .or. SE % trm) then
      SE % ed = '(L1)'
      if (size(x) == 0) then
        wid = 0
      else
        wid = 1
      endif
      SE % w = 1
      nbl = SE % w - wid
    else
      wid = SE % w
      nbl = 0
    endif
    call tobox_log1(title, x, SE, wid, nbl)
  end subroutine disp_log1

  subroutine tobox_log1(title, x, SE, wid, nbl)
    character(*),   intent(in)    :: title
    logical(log1),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer,        intent(inout) :: wid(:)
    integer,        intent(inout) :: nbl(:)
    character(SE % w)  :: s(size(x,1))
    integer            :: m, n, lin1, i, j, wleft, widp(size(wid))
    character, pointer :: boxp(:,:)
    m = size(x,1)
    n = size(x,2)
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      if (m > 0) write(s, SE % ed) (x(i,j), i=1,m)
      call copytobox(s, lin1, wid(j), widp(j), nbl(j), boxp,  wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_log1

  ! ********** 1-BYTE LOGICAL TOSTRING PROCEDURES *********
  function tostring_s_log1(x) result(st)
    logical(log1), intent(in) :: x
    character(1)            :: st
    st = tostring_f_log1((/x/), 'L1')
  end function tostring_s_log1

  function tostring_sf_log1(x, fmt) result(st)
    logical(log1),intent(in)        :: x
    character(*), intent(in)        :: fmt
    character(len_f_log1((/x/), fmt)) :: st
    st = tostring_f_log1((/x/), fmt)
  end function tostring_sf_log1

  function tostring_log1(x) result(st)
    logical(log1), intent(in)                          :: x(:)
    character(1 + (size(x) - 1)*(1 + tosset % seplen)) :: st
    st = tostring_f_log1(x, 'L1')
  end function tostring_log1

  function tostring_f_log1(x, fmt) result(st)
    logical(log1), intent(in)     :: x(:)
    character(*), intent(in)      :: fmt
    character(len_f_log1(x, fmt)) :: st
    character(widthmax_log1(fmt)) :: sa(size(x))
    integer                       :: w, d
    logical                       :: gedit
    character(nnblk(fmt)+2)       :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) then; st = errormsg; return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES') sa = adjustl(sa)
    call tostring_get(sa, st)
  end function tostring_f_log1

  pure function len_f_log1(x, fmt) result(wtot)
    logical(log1), intent(in)  :: x(:)
    character(*), intent(in)   :: fmt
    integer                    :: wtot, w, d
    logical                    :: gedit
    character(nnblk(fmt)+2)    :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) then; wtot = len(errormsg); return; endif
    if (tosset % trimb == 'YES') wtot = size(x)
    if (tosset % trimb == 'NO' ) wtot = w*size(x)
    wtot = wtot + (size(x) - 1)*(tosset % seplen)
  end function len_f_log1

  pure function widthmax_log1(fmt) result(w)
    character(*), intent(in) :: fmt
    integer w, d
    logical gedit
    character(nnblk(fmt)+5) :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) w = 1
  end function widthmax_log1

END MODULE DISP_L1MOD
