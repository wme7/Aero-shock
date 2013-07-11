MODULE DISP_I1MOD

  ! Add-on module to DISPMODULE to display 1-byte integers
  ! (assuming that these are obtained with selected_int_kind(2))
  ! Version that works with the Absoft Fortran 95 compiler, version 10.1
  !
  ! This module is obtained by copying the section DEFAULT INTEGER PROCEDURES from
  ! dispmodule.f90, replacing dint with byte and 'default integer' with 1-byte
  ! integer (only appears in comments), and adding the DECLARATIONS section below.
  !
  ! Copyright (c) 2008, Kristján Jónasson, Dept. of Computer Science, University of
  ! Iceland (jonasson@hi.is). This software is free. For details see the file README.
  
  ! ******************************** DECLARATIONS ********************************************
  use dispmodule_util

  PUBLIC DISP
  PUBLIC TOSTRING

  PRIVATE

  interface disp
    module procedure disp_s_byte, disp_ts_byte, disp_v_byte, disp_tv_byte, disp_m_byte, disp_tm_byte
  end interface

  interface tostring
    module procedure tostring_byte, tostring_f_byte, tostring_s_byte, tostring_sf_byte
  end interface

  integer, parameter :: byte = selected_int_kind(2)

CONTAINS

  ! ******************************** 1-BYTE INTEGER PROCEDURES *******************************
  subroutine disp_s_byte(x, fmt, advance, sep, trim, unit, zeroas)
    ! 1-byte integer scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim, zeroas
    integer(byte), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_ts_byte('', x, fmt, advance, sep, 'left', trim, unit, zeroas)
  end subroutine disp_s_byte

  subroutine disp_v_byte(x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    ! 1-byte integer vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    integer(byte), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tv_byte('', x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
  end subroutine disp_v_byte

  subroutine disp_m_byte(x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
    ! 1-byte integer matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    integer(byte), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tm_byte('', x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
  end subroutine disp_m_byte

  subroutine disp_ts_byte(title, x, fmt, advance, sep, style, trim, unit, zeroas)
    ! 1-byte integer scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    integer(byte), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_tm_byte(title, reshape((/x/), (/1, 1/)), fmt, advance, sep=sep, style=style, trim=trim, unit=unit, &
         zeroas=zeroas)
  end subroutine disp_ts_byte

  subroutine disp_tv_byte(title, x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    ! 1-byte integer vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    integer(byte), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    if (SE % row) then
      call disp_byte(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_byte(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_byte

  subroutine disp_tm_byte(title, x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
    ! 1-byte integer matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    integer(byte),intent(in)           :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Format edit descriptor to use for each matrix element (e.g.'I4')
    integer,      intent(in), optional :: unit       ! Unit to display on
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: zeroas     ! Zeros are replaced by this string
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no trimming,
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, zeroas=zeroas)
    call disp_byte(title, x, SE)
  end subroutine disp_tm_byte

  subroutine disp_byte(title, x, SE)
    ! 1-byte integer item
    character(*),   intent(in)    :: title
    integer(byte),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer wid(size(x,2)), nbl(size(x,2))
    call find_editdesc_byte(x, SE, wid, nbl) ! determine also SE % w
    call tobox_byte(title, x, SE, wid, nbl)
  end subroutine disp_byte

  subroutine tobox_byte(title, x, SE, wid, nbl)
    ! Write 1-byte integer matrix to box
    character(*),   intent(in)    :: title
    integer(byte),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer,        intent(inout) :: wid(:)
    integer,        intent(inout) :: nbl(:)
    character(SE % w)  :: s(size(x,1))
    integer            :: lin1, j, wleft, m, n, widp(size(wid))
    character, pointer :: boxp(:,:)
    m = size(x,1)
    n = size(x,2)
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      if (m > 0) write(s, SE % ed) x(:,j)
      if (SE % lzas > 0) call replace_zeronaninf(s, SE % zas(1:SE % lzas), x(:,j) == 0)
      call copytobox(s, lin1, wid(j), widp(j), nbl(j), boxp,  wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_byte

  subroutine find_editdesc_byte(x, SE, wid, nbl)
    ! Determine SE % ed, SE % w (unless specified) and wid
    integer(byte),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer,        intent(out)   :: wid(size(x,2)), nbl(size(x,2))
    !
    integer(byte) xmaxv(size(x,2)), xminv(size(x,2)), xp, xm
    logical xzero(size(x,2)), xallz(size(x,2))
    character(22) s
    integer ww
    !
    if (SE % w == 0) then
      xp = maxval(x)
      xm = minval(x)
      write(s, '(SS,I0)') xp; ww = len_trim(s)
      write(s, '(SS,I0)') xm; ww = max(ww, len_trim(s))
      SE % w = max(SE % lzas, ww)
      call replace_w(SE % ed, ww)
    elseif (SE % w < 0) then ! obtain max-width of x
      if (size(x) == 0) then
        SE % ed = '()'
        SE % w = 0
        wid = 0
        return
      endif
      xp = maxval(x)
      xm = minval(x)
      write(s, '(SS,I0)') xp; ww = len_trim(s)
      write(s, '(SS,I0)') xm; ww = max(ww, len_trim(s))
      ww = max(SE % lzas, ww)
      SE % ed = '(SS,Ixx)'
      write(SE % ed(6:7), '(SS,I2)') ww
      SE % w = ww
    endif
    if (SE % trm) then
      xmaxv = maxval(x, 1) ! max in each column
      xminv = minval(x, 1) ! min
      xzero = any(x == 0_byte, 1) ! true where column has some zeros
      xallz = all(x == 0_byte, 1) ! true where column has only zeros
      call getwid_byte(xmaxv, xminv, xzero, xallz, SE,  wid, nbl)
    else
      wid = SE % w
      nbl = 0
    endif
  end subroutine find_editdesc_byte

  subroutine getwid_byte(xmaxv, xminv, xzero, xallz, SE,  wid, nbl)
    integer(byte),  intent(in)  :: xmaxv(:), xminv(:)
    logical,        intent(in)  :: xzero(:), xallz(:) ! True for columns with some/all zeros
    type(settings), intent(in)  :: SE                 ! Settings
    integer,        intent(out) :: wid(:)             ! Widths of columns
    integer,        intent(out) :: nbl(:)             ! n of blanks to peel from left (w-wid)
    character(SE % w) :: stmax(size(xmaxv)), stmin(size(xmaxv))
    integer w
    w = SE % w
    write(stmax, SE % ed) xmaxv
    write(stmin, SE % ed) xminv
    nbl = mod(verify(stmin, ' ') + w, w + 1) ! loc. of first nonblank
    nbl = min(nbl, mod(verify(stmax, ' ') + w, w + 1))
    wid = w - nbl
    if (SE % lzas > 0) then
      wid = merge(SE % lzas, wid, xallz)
      wid = max(wid, merge(SE % lzas, 0, xzero))
      nbl = w - wid
    endif
  end subroutine getwid_byte
  
  ! ********* 1-BYTE INTEGER TOSTRING PROCEDURES *********
  function tostring_s_byte(x) result(st)
    ! Scalar to string
    integer(byte), intent(in)                      :: x
    character(len_f_byte((/x/), 1, tosset % ifmt)) :: st
    st = tostring_f_byte((/x/), tosset % ifmt)
  end function tostring_s_byte

  function tostring_sf_byte(x, fmt) result(st)
    ! Scalar with specified format to string
    integer(byte),intent(in)             :: x
    character(*), intent(in)             :: fmt
    character(len_f_byte((/x/), 1, fmt)) :: st
    st = tostring_f_byte((/x/), fmt)
  end function tostring_sf_byte

  function tostring_byte(x) result(st)
    ! Vector to string
    integer(byte), intent(in) :: x(:)
    character(len_f_byte(x, size(x), tosset % ifmt)) :: st
    st = tostring_f_byte(x, tosset % ifmt)
  end function tostring_byte

  function tostring_f_byte(x, fmt) result(st)
    ! Vector with specified format to string
    integer(byte), intent(in)                 :: x(:)
    character(*), intent(in)                  :: fmt
    character(len_f_byte(x,  size(x), fmt))   :: st
    character(widthmax_byte(x, size(x), fmt)) :: sa(size(x))
    integer                                   :: w, d
    logical                                   :: gedit
    character(nnblk(fmt)+5)                   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; st = errormsg; return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    call tostring_get(sa, st)
  end function tostring_f_byte

  pure function len_f_byte(x, ix, fmt) result(wtot)
    ! Total width of tostring representation of x
    integer, intent(in)                       :: ix     ! This parameter is needed for Absoft's compiler
    integer(byte), intent(in)                 :: x(ix)  ! (with other compilers x can be dimensioned with x(:))
    character(*), intent(in)                  :: fmt
    character(widthmax_byte(x, size(x), fmt)) :: sa(size(x))
    integer                                   :: wtot, w, d
    logical                                   :: gedit
    character(nnblk(fmt)+5)                   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    wtot = sum(len_trim(sa)) + (size(x) - 1)*(tosset % seplen)
  end function len_f_byte

  pure function widthmax_byte(x, ix, fmt) result(w)
    ! Maximum width of string representation of an element in x
    integer, intent(in)        :: ix
    integer(byte), intent(in)  :: x(ix)  ! ix is a parameter to avoid bug in Absoft f95
    character(*), intent(in)   :: fmt
    character(range(x)+2) sx(2)
    integer w, d
    logical gedit
    character(nnblk(fmt)+5) :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w<=0) then
      write(sx, '(SS,I0)') maxval(x), minval(x)
      w = maxval(len_trim(sx))
    end if
  end function widthmax_byte
  ! ************************************* END OF 1-BYTE INTEGER PROCEDURES ******************************************

END MODULE DISP_I1MOD
