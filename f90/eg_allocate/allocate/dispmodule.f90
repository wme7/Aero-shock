! DISPMODULE, A FORTRAN 95 MODULE FOR PRETTY-PRINTING MATRICES.
! Version number 1.02 6-Sept-2008
!
! NOTE: THE MAIN MODULE, DISPMODULE, IS LATER IN THIS FILE. 
!
! The usage documentation for DISPMODULE is in a separate document, that exists
! in several formats:
!
!   dispmodule_userman.doc   Word 2003 doc file
!   dispmodule_userman.pdf   PDF file
!   dispmodule_userman.html  HTML file
!   dispmodule_userman.txt   Text file
!    
! Copyright (c) 2008, Kristján Jónasson, Dept. of Computer Science, University of
! Iceland (jonasson@hi.is). This software is free. For details see the file README.


MODULE PUTSTRMODULE ! DUMMY VERSION
  ! An auxilliary module that accompanies DISPMODULE. This module contains dummy versions of the
  ! subroutines putstr and putnl that do nothing. It is needed to avoid an "undefined symbol" link
  ! error for these. In addition it defines the named constant (or parameter) DEFAULT_UNIT = -3,
  ! which makes the asterisk unit (usually the screen) the default to display on.
  !
  ! The purpose of having this module is to make displaying possible in situations where ordinary
  ! print- and write-statements do not work. Then this module should be replaced by one defining
  ! functional versions of putstr and putnl. An example is given by the commented out PUTSTRMODULE 
  ! for Matlab mex files below.
  !
  integer, parameter :: DEFAULT_UNIT = -3
  !
CONTAINS
  subroutine putstr(s)
    character(*), intent(in) :: s
    integer ldummy, ldummy1  ! these variables exist to avoid unused variable warnings
    ldummy = len(s)
    ldummy1 = ldummy
    ldummy = ldummy1      
  end subroutine putstr
  subroutine putnl()
  end subroutine putnl
END MODULE PUTSTRMODULE



! MODULE PUTSTRMODULE  ! for Matlab mex files.
!   ! This module contains functional versions of subroutines putstr and putnl. It also sets
!   ! DEFAULT_UNIT = -2, which makes putstr/putnl the default to display with. Using this module,
!   ! instead of the dummy module above allows DISPMODULE to be used with Matlab mex files.
!   ! used (commented in) instead of the one above (which should then be commented out), then
!   ! DISPMODULE can be used with Matlab mex files. A shorter version (given in the user manual)
!   ! may be used with g95, but the one below works for both g95 and gfortran.
!   !
!   use, intrinsic :: ISO_C_BINDING
!   integer, parameter :: default_unit = -2
!   interface
!     subroutine mexprintf(s) bind(C, name = 'mexPrintf')
!       import c_char
!       character(c_char) s(*)
!     end subroutine mexprintf
!   end interface
! CONTAINS
!   subroutine putstr(s)
!     character(*), intent(in) :: s
!     call mexprintf(s//char(0))
!   end subroutine putstr
!   subroutine putnl()
!     call mexprintf(char(10)//char(0))
!   end subroutine putnl
! END MODULE PUTSTRMODULE



MODULE DISPMODULE_UTIL
  ! Dispmodule_util contains utilities that are used by Dispmodule, and the add-on modules
  ! disp_i1mod, disp_i2mod,..., disp_l1mod and disp_r16mod. Note that the entities that are
  ! declared public below are not exported to the user. The private statements in dispmodule and
  ! the add-on modules prevent that from happening.

  use putstrmodule
  implicit none

  ! ***************** PUBLIC ENTITIES (ONLY PUBLIC TO DISPMODULE, NOT TO USER PROGRAMS) *****************
  private
  public disp_settings, defset, factory_settings, tosset, tosfac, errormsg, tostring_settings
  public nnblk, upper, readfmt, replace_w, trim_real, get_SE, preparebox, copytobox, boxlist, boxnode
  public copyseptobox, finishbox, tostring_get_complex, disp_errmsg, tostring_get, find_editdesc_real
  public check_settings, tostring_check_settings, replace_zeronaninf, settings, trim_s_real

  ! *********************************** GENERAL DECLARATIONS ********************************************
  type disp_settings
    ! Settings used by subroutine disp and the utility procedures.
    character(6) :: advance     = 'YES'
    character(9) :: matsep      = '   '
    character(3) :: orient      = 'COL'
    character(9) :: sep         = '  '
    character(9) :: style       = 'LEFT'
    character(4) :: trim        = 'AUTO'
    character(9) :: zeroas      = ''
    integer      :: digmax      = 6
    integer      :: matseplen   = 3
    integer      :: seplen      = 2
    integer      :: unit        = DEFAULT_UNIT
    integer      :: zaslen      = 0
  end type disp_settings

  type tostring_settings
    ! Settings used by function tostring. 
    character(10) :: ifmt = 'I0'
    character(16) :: rfmt = '1PG12.5'  ! 'SP,1P,G20.11E3' has length 14 and is about max
    character(9)  :: sep = ', '
    integer       :: seplen = 2
    character(3)  :: trimb = 'YES'
    character(4)  :: trimz = 'G'
  end type tostring_settings

  type settings
    ! Settings used (privately) by disp and the utility procedures, in the variable SE.
    character(22) ed
    character(9) sep, tsty, zas
    character(1) tch
    integer lun, dmx, w, d, lsep, lzas, m1, n1, adv
    logical trm, number, vec, row, gedit
  end type settings

  type(disp_settings), save :: DEFSET, &        ! Current default settings for disp
       &                       FACTORY_SETTINGS ! Original (factory) settings for disp
  type(tostring_settings), save :: tosset, & ! Current settings for tostring
       &                           tosfac    ! Factory settings for tostring

  character(*), parameter :: errormsg = 'Illegal format'

  ! ********************* BOX-PACKAGE DECLARATIONS (SEE EXPLANATION ABOUT BOX-PACKAGE BELOW) *****************
  type boxnode
    ! A box is the character representation of a printed item
    character, pointer     :: box(:,:)
    type(boxnode), pointer :: nextbox => null()
  end type boxnode
  !
  type boxlist
    ! There is one list of boxes associated with each logical unit
    integer :: unit = 1
    type(boxnode), pointer :: firstbox => null()
    type(boxnode), pointer :: lastbox => null()
    type(boxlist), pointer :: nextboxlist => null()
  end type boxlist
  !
  type(boxlist), pointer :: firstboxlist => null()
  ! ************************ END OF BOX-PACKAGE DECLARATIONS ******************************

CONTAINS

  ! ***************************** GENERAL PROCEDURES **************************************
  subroutine check_settings()
    ! Sanity check of display settings
    character(9) :: tsty
    character tch
    logical number, ok, dmxerr, orierr, styerr, adverr
    character(6), parameter :: ADVOK(3) = (/'NO    ', 'YES   ', 'DOUBLE'/)
    type(disp_settings) ds
    ds = DEFSET
    call getstyles(ds % style, tsty, tch, number, ok)
    styerr = .not. ok
    dmxerr = ds % digmax < 1 .or. ds % digmax > 89
    orierr = all(ds % orient /= (/'ROW', 'COL'/))
    adverr = all(ds % advance /= ADVOK)
    if (dmxerr) DEFSET % digmax = 6
    if (orierr) DEFSET % orient = 'COL'
    if (styerr) DEFSET % style = 'LEFT'
    if (adverr) DEFSET % advance = 'YES'
    !
    if (dmxerr) call disp_errmsg('DISP_SET: error, illegal digmax (must be 1-89), set to 6')
    if (orierr) call disp_errmsg('DISP_SET: error, illegal orient: ' // trim(ds % orient) // ', set to "COL"')
    if (styerr) call disp_errmsg('DISP_SET: error, illegal style: ' // trim(ds % style) // ', set to "LEFT"')
    if (adverr) call disp_errmsg('DISP_SET: error, illegal advance: ' // trim(ds % advance) // ', set to "YES"')
  end subroutine check_settings

  function number_rows(SE) result(nbr)
    ! Should rows be numbered?
    type(settings), intent(in) :: SE
    logical nbr
    nbr = .false.
    if (.not. SE % number) return
    if (SE % vec .and. SE % row) return
    nbr = .true.
  end function number_rows

  function number_cols(SE) result(nbr)
    ! Should columns be numbered?
    type(settings), intent(in) :: SE
    logical nbr
    nbr = .false.
    if (.not. SE % number) return
    if (SE % vec .and. .not. SE % row) return
    nbr = .true.
  end function number_cols

  subroutine preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    ! Determine format to use to write matrix to box and row where matrix begins, copy 
    character(*),   intent(in)    :: title     ! The title to use for the matrix
    type(settings), intent(in)    :: SE        ! Settings
    integer,        intent(in)    :: m         ! Row count of matrix
    integer,        intent(in)    :: n         ! Column count of matrix
    integer,        intent(inout) :: wid(:)    ! widths of columns in matrix
    integer,        intent(out)   :: widp(:)   ! widths of columns in box (max(wid, width of col nums))
    integer,        intent(out)   :: lin1      ! Row number where matrix begins (tsty='left' 0, 'pad' 1, 'underline' 2)
    integer,        intent(out)   :: wleft     ! Number of spaces on left of matrix (when tsty is left or title long)
    character, pointer            :: boxp(:,:) ! The box

    integer wt, wa          ! Char count of title, idth of matrix in characters (wbox = lm + wa + rm)
    integer wbox, wrow      ! Width of box in characters, width of row numbers in characters
    integer lm              ! Left margin
    integer h,ws            ! Height of box in characters, length of column separator
    integer m1, n1, i       ! lower bounds (for numbering), index
    character(range(0) + 2) sn(2), row_nums(m), col_nums(n)
    character(10) fmt

    ! ----------wbox---------    -----------wbox----------     -----wbox------
    ! ---lm---                   --wleft-                             --wt-
    ! ----wleft---                lm wrow    wa       rm       wrow    wa
    !    wt   wrow    wa         ----====-----------======     ----===========
    ! --------====-----------    THIS-IS-A-VERY-LONG-TITLE            TITLE
    !               1     2                1     2                   1     2      
    ! MATRIX = 1   4.50  6.80         1   4.50  6.80            1   4.50  6.80
    !          2   6.88  9.22         2   6.88  9.22            2   6.88  9.22
    !          3  19.44  0.08         3  19.44  0.08            3  19.44  0.08
    !          ...                    ...                       ...
    !         10   6.18  4.22        10   6.18  4.22           10   6.18  4.22
    ! rm = 0                     wt = wbox                     lm = rm = 0, wleft = wrow
    m1 = SE % m1
    n1 = SE % n1
    ws = SE % lsep
    wt = len(title)
    wrow = 0
    widp = wid
    if (SE % number) then
      fmt = '(SS,I0)'
      if (number_cols(SE)) then 
        write(col_nums, fmt) (/ (i, i = n1, n1 + n - 1) /)
        widp = max(wid, len_trim(col_nums))
      endif
      if (number_rows(SE)) then
        write(sn, fmt) m1, m1 + m - 1
        wrow = maxval(len_trim(sn)) + ws  ! determine max width of row numbers
        call replace_w(fmt, wrow - ws) ! to create e.g. 'I5' from 'I0'
        write(row_nums, fmt) (/ (i, i = m1, m1 + m - 1) /)
      endif
    endif
    wa = max(0,n-1)*ws + sum(widp)
    select case(upper(SE % tsty))
    case('LEFT');      lin1 = 1; wbox = wt + wrow + wa;     h = max(1,m); lm = wt
    case('PAD');       lin1 = 2; wbox = max(wt, wa + wrow); h = m + 1;    lm = max(0, (wt - wa - wrow)/2)
    case('UNDERLINE'); lin1 = 3; wbox = max(wt, wa + wrow); h = m + 2;    lm = max(0, (wt - wa - wrow)/2)
    case default;      lin1 = 1; wbox = 0; h = 0; lm = 0 ! should not happen
    end select
    wleft = lm
    if (number_cols(SE)) h = h + 1
    call newbox(SE % lun, h, wbox, boxp)
    if (number_cols(SE)) then
      call copycolumnnumberstobox(col_nums, wleft + wrow, wid, widp, ws,  boxp, lin1)
    endif
    if (number_rows(SE)) then
      call copytobox(row_nums, lin1, wrow - ws, wrow - ws, nblj = 0, boxp = boxp, wleft = wleft)
      call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp, wleft)
    endif
  end subroutine preparebox

  subroutine copytobox(s, lin1, widj, widpj, nblj, boxp,  wleft)
    ! Copy strings to column in boxp; update wleft to current char column in boxp
    character(*), intent(in)    :: s(:)        ! the strings to copy
    integer,      intent(in)    :: lin1, widj  ! first line in box to copy to, width of column
    integer,      intent(in)    :: nblj, widpj ! number of blank characters to trim from left of s, offset to next col
    character,    intent(inout) :: boxp(:,:)   ! the box to accept the column
    integer,      intent(inout) :: wleft       ! number of char-columns in box already written to
    integer i, j
    wleft = wleft + widpj - widj
    forall(i = 1:widj, j=1:size(s)) boxp(wleft+i, j+lin1-1) = s(j)(i+nblj:i+nblj)
    wleft = wleft + widj
  end subroutine copytobox

  subroutine copyseptobox(sep, m, lin1, boxp,  wleft)
    ! Copy column separator to boxp; update wleft
    character(*), intent(in)    :: sep
    integer,      intent(in)    :: m, lin1
    character,    intent(inout) :: boxp(:,:)
    integer,      intent(inout) :: wleft
    integer i, j
    forall(i = 1:len(sep), j=1:m) boxp(wleft+i, j+lin1-1) = sep(i:i)
    wleft = wleft + len(sep)
  end subroutine copyseptobox

  subroutine copycolumnnumberstobox(s, wleft, wid, widp, lsep, boxp, lin1)
    character(*), intent(in)    :: s(:)      ! strings with left-adjusted column numbers
    integer,      intent(in)    :: wleft     ! char positions on left of 1st col
    integer,      intent(in)    :: wid(:)    ! widths of columns in matrix
    integer,      intent(in)    :: widp(:)   ! widths of columns in box (max(wid, width of col nums))
    integer,      intent(in)    :: lsep      ! width of column separator
    character,    intent(inout) :: boxp(:,:) ! receives the numbers
    integer,      intent(inout) :: lin1      ! line number in box to copy to
    integer ls(size(s)), rmargmax, k, i, lmargin, j
    !
    ls = len_trim(s)
    rmargmax = (max(0, minval(wid) - maxval(ls)))/2 ! locate according to narrowest column, widest number
    k = wleft
    do i = 1, size(wid)
      lmargin = max(0, widp(i) - ls(i) - rmargmax)
      k = k + lmargin
      forall(j = 1:ls(i)) boxp(k+j, lin1) = s(i)(j:j)
      k = k + widp(i) - lmargin + lsep
    enddo
    lin1 = lin1 + 1
  end subroutine copycolumnnumberstobox

  subroutine finishbox(title, SE, boxp)
    ! Finish creating a box and display it if advancing is turned on
    character(*),   intent(in)    :: title     ! The title to use for the matrix
    type(settings), intent(in)    :: SE        ! Settings
    character,      intent(inout) :: boxp(:,:) ! The box
    !
    integer i, wt, w, wpadright, wpadleft ! index, width of title, width of box and spacing on either side of it
    integer lin1 ! line to put left title
    !
    wt = len(title)
    w = size(boxp,1)
    if (upper(SE % tsty) == 'LEFT') then
      lin1 = 1
      if (number_cols(SE)) lin1 = min(2,size(boxp,2))
      forall(i=1:wt) boxp(i,lin1) = title(i:i)
    else
      wpadright = (w - wt)/2
      wpadleft = w - wpadright - wt
      forall(i=1:wt) boxp(wpadleft+i, 1) = title(i:i)
      if (upper(SE % tsty) == 'PAD') then
        boxp(1:wpadleft, 1) = SE % tch
        boxp(w-wpadright+1:w, 1) = SE % tch
      else ! tsty == 'UNDERLINE'
        boxp(:,2) = SE % tch
      endif
    endif
    if (SE % adv >= 1) call dispboxlist(SE % lun, DEFSET % matsep(1:DEFSET % matseplen))
    if (SE % adv >= 2) call dispnewline(SE % lun)
  end subroutine finishbox

  subroutine find_editdesc_real(exp, expm, dmx,  edesc, flen, ndec, posit) 
    ! Subroutine of find_editdesc_sngl and find_editdesc_dble
    integer,       intent(in)    :: expm, dmx
    integer,       intent(inout) :: exp
    character(14), intent(out)   :: edesc
    integer,       intent(out)   :: flen, ndec
    logical,       intent(in)    :: posit
    integer :: neg, nxp
    exp = max(exp, expm)
    neg = 1
    if (exp < dmx .and. exp >= -1) then
      if (posit .or. exp > max(0, expm)) neg = 0
      edesc = '(SS,Fxx.yy)'
      ndec = max(0, dmx - exp - 1)
      flen = neg + 2 + ndec + max(0,exp) ! -X.YYYYY (2 covers X and .)
      write(edesc(6:10), '(SS,I2,".",I2)') flen, ndec
    else
      if (posit) neg = 0
      if     (abs(exp) > 999) then; nxp = 4
      elseif (abs(exp) >  99) then; nxp = 3
      elseif (abs(exp) >   9) then; nxp = 2
      else                        ; nxp = 1
      endif
      flen = neg + 3 + dmx + nxp
      edesc = '(SS,ESxx.yyEz)'
      write(edesc(7:13), '(SS,I2,".",I2,"E",I1)') flen, dmx - 1, nxp
      ndec = dmx - 1
    endif
  end subroutine find_editdesc_real

  pure subroutine readfmt(fmt, fmt1, w, d, gedit)
    ! Returns w and d when fmt is (Xw.d) or (Xw) (then d = 0), X = edit descriptor letter
    ! (I, F, etc). X can also be ES, DS, 1PG or 1PF. Returns w = -1 for illegal fmt.
    ! Returns gedit = .true. if fmt is Gw.d. How about SS,1PES4.3?
    character(*), intent(in)  :: fmt  ! e.g. fmt = F 8.2
    character(*), intent(out) :: fmt1 ! returns '(SS,F8.2)'
    character ch
    integer, intent(out) :: w, d
    logical, intent(out) :: gedit
    integer :: k0, k1, k2, k3, k4
    call sszipfmt(fmt, fmt1)
    w = -1; d = 0; gedit = .false.
    k1 = verify(fmt1(2:), '0123456789') + 1
    if (k1 == 0) return ! only digits
    k2 = verify(fmt1(k1:), 'ABDEFGILNOPSZabdefgilnopsz,') + k1 - 1 ! , for "1P,G12.3"
    if (k2 <= k1) return ! no letter or only letters
    ch = upper(fmt1(k2-1:k2-1))
    if (ch == ',') then ! deal with SS,1PG13.5
      k0 = k2
      k1 = verify(fmt1(k0:),'0123456789') + k0 - 1
      if (k1==0) return
      k2 = verify(fmt1(k1:),'ABDEFGILNOPSZabdefgilnopsz,') + k1 - 1
      if (k2 <= k1) return
      ch = upper(fmt1(k2-1:k2-1))
    endif
    gedit = ch == 'G' .or. ch == 'g'
    k3 = verify(fmt1(k2:), '0123456789') + k2 - 1
    if (k3 == k2) return ! no digits
    read(fmt1(k2:k3-1), *) w
    if (k3 > len(fmt1)) return
    if (fmt1(k3:k3) /= '.') return ! not . after w
    k4 = verify(fmt1(k3+1:), '0123456789') + k3
    if (k4 == k3+1) return ! no digits
    read(fmt1(k3+1:k4-1), *) d
  end subroutine readfmt

  pure subroutine replace_w(fmt, wnew)
    ! Change e.g. '(F0.3)' to '(F5.3)'. Works also for '(SS,I0)' to '(SS,I5)'. If wnew > 999, set it to 999
    character(*), intent(inout) :: fmt
    integer, intent(in) :: wnew
    integer :: k0, k1, k2, k3
    character(3) rw
    k1 = verify(fmt(2:), '0123456789') + 1
    k2 = verify(fmt(k1:), 'ABDEFGILNOPSZabdefgilnopsz,') + k1 - 1
    if (k2 == k1) return ! no letter
    if (fmt(k2-1:k2-1)==',') then ! Handle (SS,1PF10.3)
      k0 = k2
      k1 = verify(fmt(k0:),'0123456789') + 1
      if (k1==0) return
      k2 = verify(fmt(k1:),'ABDEFGILNOPSZabdefgilnopsz,') + k1 - 1
      if (k2 <= k1) return
    end if
    k3 = verify(fmt(k2:), '0123456789') + k2 - 1
    if (k3 == k2) return ! no digits
    write(rw, '(SS,I0)') min(999,wnew)
    fmt = fmt(1:k2-1) // trim(rw) // fmt(k3:)
  end subroutine replace_w

  subroutine get_SE(SE, title, shapex, fmt, advance, lbound, seperator, style, trim, unit, orient, zeroas, digmax)
    ! Get the settings from the optional parameters fmt...zeroas in to the structure SE.
    ! Replace absent arguments with corresponding values from the structure DEFSET.
    type(settings), intent(out)          :: SE
    character(*),   intent(in)           :: title
    integer,        intent(in)           :: shapex(:)
    character(*),   intent(in), optional :: fmt
    integer,        intent(in), optional :: unit, digmax, lbound(:)
    character(*),   intent(in), optional :: advance, seperator, style, zeroas, trim, orient
    logical ok
    !
    character(22) ed
    character(9) sep, tsty, zas
    character(1) tch
    character(6) advchr
    integer lun, dmx, w, d, lsep, lzas, m1, n1, adv
    logical trm, number, vec, row, is_scalar, gedit
    !
    vec = (size(shapex) == 1)
    is_scalar = size(shapex) == 0
    if (vec .and. present(orient)) then
      select case(upper(orient))
      case('ROW');  row = .true.
      case('COL');  row = .false.
      case default; 
        call disp_errmsg('DISP: error, wrong value of orient: '//orient(1:len_trim(orient))//', using "COL"')
        row = .false.
      end select
    elseif (vec) then
      row = DEFSET % orient == 'ROW'
    else
      row = .false.
    endif
    if (present(fmt)) then
      call readfmt(fmt, ed, w, d, gedit)
    else
      ed = '()'
      w = -1; d = 0; gedit = .false.
    endif
    if (present(unit)) then
      lun = unit
    else
      lun = DEFSET % unit
    endif
    if (.not.present(digmax)) then
      dmx = DEFSET % digmax
    elseif (present(fmt)) then
      call disp_errmsg('DISP: error, both FMT and DIGMAX present, ignoring DIGMAX')
      dmx = 1
    elseif (digmax < 1 .or. digmax > 89) then
      call disp_errmsg('DISP: error, digmax must be >= 1 and < 90, using 6')
      dmx = 6
    else
      dmx = digmax
    endif
    if (present(advance)) then
      advchr = upper(advance)
    else
      advchr = DEFSET % advance
    endif
    select case(trims(advchr))
    case('NO');     adv = 0
    case('YES');    adv = 1
    case('DOUBLE'); adv = 2
    case default
      call disp_errmsg('DISP: error, illegal advance: ' // trims(advance) // ', using "YES"')
      adv = 1
    end select
    if (present(trim)) then
      if (upper(trim) /= 'YES' .and. upper(trim) /= 'NO' .and. upper(trim) /= 'AUTO') then
        call disp_errmsg('DISP: error, illegal trim: ' // trims(trim) // ', using "YES"')
        trm = .true.
      else
        trm = upper(trim) == 'YES' .or. upper(trim) == 'AUTO' .and. .not.present(FMT)
      endif
    elseif (w == 0) then
      trm = .true.
    else
      trm = DEFSET % trim == 'YES' .or. DEFSET % trim == 'AUTO' .and. .not.present(FMT)
    endif
    if (present(seperator)) then
      sep = seperator
      lsep = len(seperator)
    else
      sep = DEFSET % sep
      lsep = DEFSET % seplen
    endif
    if (present(style)) then
      call getstyles(style, tsty, tch, number, ok)
      if (.not. ok) call disp_errmsg('DISP: error, illegal style: '//style//'. Using default instead')
    else
      call getstyles(DEFSET % style, tsty, tch, number, ok)
    endif
    if (title == '') tsty = 'LEFT'
    if (is_scalar) number = .false.
    if (present(zeroas)) then
      zas = zeroas
      lzas = len(zeroas)
    else
      zas = DEFSET % zeroas
      lzas = DEFSET % zaslen
    endif
    if (w > 0) lzas = min(w, lzas)
    zas = zas(1:lzas)
    m1 = 1
    n1 = 1
    if (present(lbound)) then
      number = .true.
      if (size(lbound) == 1) then        
        if (vec .and. row) then
          n1 = lbound(1)
        else
          m1 = lbound(1)
        endif
      elseif (size(lbound) >= 2) then
        m1 = lbound(1)
        n1 = lbound(2)
      endif
    endif
    SE = settings(ed, sep, tsty, zas, tch, lun, dmx, w, d, lsep, lzas, m1, n1, adv, trm, number, vec, row, gedit)
  contains
    function trims(s) result(t)
      character(*), intent(in) :: s
      character(len_trim(s)) :: t
      intrinsic trim
      t = trim(s)
    end function trims
  end subroutine get_SE

  subroutine getstyles(style, tsty, tch, number, ok)
    ! Return tsty = 'LEFT', 'PAD', or 'UNDERLINE', tch = x from xPAD or xUNDERLINE, number = .true. if style includes
    ! NUMBER. If style has ABOVE, return tsty = 'PAD' and tch = ' '. Return tsty = 'LEFT' if error. See NOTE 1 below.
    character(*), intent(in) :: style
    character(9), intent(out) :: tsty
    character(1), intent(out) :: tch
    logical,      intent(out) :: number, ok
    integer kamp, i, nsty
    character(len(style))   :: sty(2)
    character(9), parameter :: LPUA(4) = (/'LEFT     ', 'PAD      ', 'UNDERLINE', 'ABOVE    '/)
    character(9), parameter :: PU(2) = (/'PAD      ', 'UNDERLINE'/)
    kamp = scan(upper(style), '&')
    ok = .true.
    if (kamp > 0) then
      sty(1) = adjustl(upper(style(1:kamp-1)))
      sty(2) = adjustl(upper(style(kamp+1:)))
      nsty = 2
    else
      sty(1) = adjustl(upper(style))
      nsty = 1
    end if
    number = .false.
    tsty = 'LEFT'
    tch = '-'
    do i = 1, nsty
      if (sty(i) == 'NUMBER') then
        number = .true.
      elseif (sty(i) == 'ABOVE') then
        tsty = 'PAD'
        tch = ' '
      elseif (any(sty(i) == LPUA)) then
        tsty = sty(i)
      elseif (any(sty(i)(2:) == PU)) then
        tsty = sty(i)(2:)
        tch = sty(i)(1:1)
      else
        ok = .false.
        return
      endif
    enddo
    ok = .true.
  end subroutine getstyles

  subroutine replace_zeronaninf(s, zas, maskz, masknan, maskminf, maskinf)
    ! replace zeros in s (where maskz is true) with zas (i.e. zero-as string) also replace nans with 'NaN',
    ! infinities with '+Inf' and minus infinities with '-Inf'. Zeros are aligned with . if zas contains .
    ! otherwise right-adjusted. Nans, and infs are right adjusted.
    ! NOTE: There are compiler bugs in current versions of both the Absoft and the Pathscale compilers
    ! so the merge calls (commented out below) had to be replaced with do loops.
    character(*), intent(inout) :: s(:)
    logical     , intent(in)    :: maskz(:), masknan(:), maskinf(:), maskminf(:)
    character(*), intent(in)    :: zas
    optional                    :: masknan, maskminf, maskinf
    character(len(s)) z, nan, minf, inf
    integer w, wz, n, i, k, zasdot
    w = len(s)
    wz = len(zas)
    n = size(maskz)
    if (wz /= 0 .and. wz <= w) then ! zas not empty and not too wide
      zasdot = index(zas, '.')
      z = ''
      if (zasdot > 0) then
        do i=1,n
          if (maskz(i)) exit
        enddo
        if (i<=n) then ! some zeros
          k = index(s(i), '.')
          if (k == 0 .or. zasdot > k .or. wz-zasdot > w-k) then ! cannot align .'s
            z(w-wz+1:) = zas ! align right
          else
            z(k-zasdot+1:k-zasdot+wz) = zas
          endif
        endif
      else
        z(w-wz+1:) = zas
      end if
      ! s = merge(z, s, maskz)
      do i=1,n
        if (maskz(i)) s(i) = z
      enddo
    endif
    if (present(masknan)) then
      if (w >= 4) then
        nan = repeat(' ', w-4) // ' NaN'
        minf = repeat(' ', w-4) // '-Inf'
        inf = repeat(' ', w-4) // '+Inf'
      elseif (w == 3) then
        nan = 'NaN'
        minf = '***'
        inf = 'Inf'
      else
        nan = repeat('*',w)
        minf = nan
        inf = nan
      endif
      ! s = merge(nan, s, masknan)
      ! s = merge(minf, s, maskminf)
      ! s = merge(inf, s, maskinf)
      do i=1,n
        if (masknan(i)) s(i) = nan
        if (maskminf(i)) s(i) = minf
        if (maskinf(i)) s(i) = inf
      enddo
    endif
  end subroutine replace_zeronaninf

  pure function upper(s) result(su) ! Change string to upper case
    character(*), intent(in) :: s
    character(len(s)) su
    character(26), parameter :: ll = 'abcdefghijklmnopqrstuvwxyz', &
         ul = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer i, k
    su = s
    do i = 1,len(s)
      k = index(ll, s(i:i))
      if (k > 0) su(i:i) = ul(k:k)
    end do
  end function upper

  pure subroutine sszipfmt(fmt, fmt1)
    ! Set fmt1 to '(SS,'//removeblanks(fmt)//')'. Caller is responsible that
    ! fmt1 has sufficient length.
    character(*), intent(in) :: fmt
    character(*), intent(out) :: fmt1
    integer i,j
    fmt1 = '(SS,'
    j = 5
    do i = 1,len(fmt)
      if (fmt(i:i) /= ' ') then
        fmt1(j:j) = fmt(i:i)
        j = j+1
      endif
    enddo
    fmt1(j:j) = ')'
  end subroutine sszipfmt

  pure function nnblk(s) result(n) ! count nonblanks in s
    character(*), intent(in) :: s
    integer i, n
    n = 0
    do i = 1,len(s)
      if (s(i:i) /= ' ') n = n+1
    enddo
  end function nnblk

  subroutine disp_errmsg(s)
    character(*), intent(in) :: s
    integer wleft
    character(1), pointer :: boxp(:,:)
    wleft = 0
    call newbox(DEFAULT_UNIT, 1, len(s), boxp)
    call copytobox((/s/), lin1 = 1, widj = len(s), widpj = len(s), nblj = 0, boxp = boxp, wleft = wleft)
    call dispboxlist(DEFAULT_UNIT, sep = '')
  end subroutine disp_errmsg
  ! *********************************** END OF GENERAL PROCEDURES *********************************

  ! ************************************* TOSTRING PROCEDURES *************************************
  subroutine tostring_check_settings
    ! Sanity check of tostring settings
    type(tostring_settings) ts
    integer wi, wr, d
    character(max(len(tosset % rfmt), len(tosset % ifmt)) + 5) fmt1
    logical gedit
    ts = tosset
    if (all(ts % trimb /= (/'YES', 'NO '/)))           tosset % trimb = tosfac % trimb
    if (all(ts % trimz /= (/'NONE', 'ALL ', 'G   '/))) tosset % trimz = tosfac % trimz
    call readfmt(tosset % rfmt, fmt1, wr, d, gedit)
    call readfmt(tosset % ifmt, fmt1, wi, d, gedit)
    if (wr < 0) tosset % rfmt = tosfac % rfmt
    if (wi < 0) tosset % ifmt = tosfac % ifmt
    if (all(ts % trimb /= (/'YES ', 'NO  ', 'AUTO'/))) call disp_errmsg( &
         'TOSTRING_SET: error, illegal trimb: '//trim(ts % trimb)//', set to ' // trim(tosfac % trimb))
    if (all(ts % trimz /= (/'NONE', 'ALL ', 'G   '/))) call disp_errmsg( &
         'TOSTRING_SET: error, illegal trimz: '//trim(ts % trimz)//', set to '//trim(tosfac % trimz))
    if (wr < 0) call disp_errmsg( &
         'TOSTRING_SET: error, illegal rfmt: '//trim(ts % rfmt)//', set to '//trim(tosfac % rfmt))
    if (wi < 0) call disp_errmsg( &
         'TOSTRING_SET: error, illegal ifmt: '//trim(ts % ifmt)//', set to '//trim(tosfac % ifmt))
  end subroutine tostring_check_settings

  pure subroutine trim_s_real(sa, gedit, w)
    ! Trim trailing zeros and possibly decimal point from fractional part.
    ! If sa = '52.2000E12' on entry then it is returned as '52.2E12   '.
    ! Whether trimming is actually done depends on tosset, gedit and w.
    character(*), intent(inout) :: sa
    logical, intent(in) :: gedit
    integer, intent(in) :: w
    integer k, k2, k3
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    if (tosset % trimz == 'ALL' .or. tosset % trimz == 'G' .and. gedit) then
      k = scan(sa, '.')
      if (k > 0) then
        k2 = verify(sa(k+1:), '0123456789') + k
        if (k2 == k) k2 = len(sa) + 1
        k3 = verify(sa(k:k2-1), '0.', back=.true.) + k - 1
        sa(k3+1:) = sa(k2:)
      endif
    endif
  end subroutine trim_s_real

  pure subroutine trim_real(sa, gedit, w)
    ! Trim trailing zeros and possibly decimal point from fractional part.
    ! If sa = '52.2000E12' on entry then it is returned as '52.2E12   '.
    ! Whether trimming is actually done depends on tosset, gedit and w.
    character(*), intent(inout) :: sa(:)
    logical, intent(in) :: gedit
    integer, intent(in) :: w
    integer i
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    if (tosset % trimz == 'ALL' .or. tosset % trimz == 'G' .and. gedit) then
      do i=1,size(sa) ! trim trailing zeros from fractional part
        call trim_s_real(sa(i), gedit, w)
      enddo
    endif
  end subroutine trim_real

  pure subroutine tostring_get(sa, st)
    ! Copy trimmed elements of sa (containing individual elements as strings) to the final
    ! tostring result st, separated by tosset % sep strings.
    character(*), intent(in)  :: sa(:)
    character(*), intent(out) :: st
    integer                   :: i, k, n, sepl
    sepl = tosset % seplen
    k = 0
    do i = 1,size(sa)
      if (k>0) st(k+1:k+sepl) = tosset % sep(1:sepl)
      if (k>0) k = k + sepl
      n = len_trim(sa(i))
      st(k+1:k+n) = trim(sa(i))
      k = k + n
    end do
  end subroutine tostring_get

  pure subroutine tostring_get_complex(sar, sgn, sai, st)
    ! Version of tostring_get for complex numbers
    character(*), intent(in)  :: sar(:), sai(:), sgn(*)
    character(*), intent(out) :: st
    integer                   :: i, k, n, sepl
    sepl = tosset % seplen
    k = 0
    do i = 1,size(sar)
      if (k>0) st(k+1:k+sepl) = tosset % sep(1:sepl)
      if (k>0) k = k + sepl
      n = len_trim(sar(i))
      st(k+1:k+n) = trim(sar(i))
      st(k+n+1:k+n+3) = ' '//sgn(i)//' '     
      k = k + n + 3
      n = len_trim(sai(i))
      st(k+1:k+n) = trim(sai(i))
      st(k+n+1:k+n+1) = 'i'
      k = k + n + 1
    end do
  end subroutine tostring_get_complex

  ! ********************************* END OF TOSTRING PROCEDURES *********************************

  ! *********************************** BOX-PACKAGE **********************************************
  !
  ! A "box" is a variable dimension character matrix that can be created dynamically. There are
  ! linked lists of boxes, one for each logical unit. When disp is called the item to be displayed
  ! is written to a box. If advance = 'no' is in effect, the writing out of the items is delayed
  ! until disp is called on the same unit with advance = 'yes' in effect; then all the boxes in
  ! the relevant list are written to the unit. There are two subroutines that are meant to be
  ! called from outside the Box-package: NEWBOX and DISPBOXLIST:
  !
  ! CALL NEWBOX(UNIT, M, N, BOXP) creates a box on unit UNIT. BOXP returns a pointer to the
  ! created box which is of type CHARACTER and DIMENSION (M,N).
  !
  ! CALL DISPBOXLIST(UNIT, SEP) writes all the boxes in the list associated with UNIT to the file
  ! on UNIT, separated with the string SEP. The following example makes this clear: let SEP = ' : '
  ! and let the first box contain XXX and the second have two rows, both equal to YYYY. Then the
  ! text written will be: XXX : YYYY : YYYY
  !
  ! To obtain tab-separated boxes when using ASCII, let SEP = char(9). After writing the boxes,
  ! the complete list is deallocated. If UNIT = -3 the asterisk unit (usually command window) is
  ! written to. If UNIT = -2 the routine putstr from the disp_where unit is used for writing. If
  ! UNIT = -1 all output will be discarded. With the iso_fortran_env module of Fortran 2003, unit
  ! may also equal OUTPUT_UNIT, unless the compiler sets that to -2.

  function getboxlist(unit) result(p)
    ! Return boxlist associated with specified unit. If this list does not exist a new list is started.
    integer, intent(in) :: unit
    type(boxlist), pointer :: p
    p => firstboxlist
    do while(associated(p))
      if (p % unit == unit) return
      p => p % nextboxlist
    end do
    allocate(p)
    p % nextboxlist => firstboxlist  ! put at head of list
    p % unit = unit
    firstboxlist => p
  end function getboxlist

  subroutine clearboxlist(unit)
    ! Deallocate all boxes associated with unit
    integer, intent(in) :: unit
    type(boxnode), pointer :: p, q
    type(boxlist), pointer :: blp
    blp => firstboxlist
    do while(associated(blp))
      if (blp % unit == unit) exit
      blp => blp % nextboxlist
    end do
    if (.not. associated(blp)) return
    p => blp % firstbox
    do while(associated(p))
      q => p
      p => p % nextbox
      deallocate(q % box)
      deallocate(q)
    enddo
    if (associated(firstboxlist, blp)) then
      firstboxlist => blp % nextboxlist
    endif
    deallocate(blp)
  end subroutine clearboxlist

  subroutine newbox(unit, m, n, boxp)
    ! Create a new box
    character, pointer :: boxp(:,:)
    integer, intent(in) :: unit, m, n
    type(boxnode), pointer :: p
    type(boxlist), pointer :: blp
    allocate(p)
    allocate(p % box(n, m))
    blp => getboxlist(unit)
    if (.not.associated(blp % firstbox)) then
      blp % firstbox => p
    else
      blp % lastbox % nextbox => p
    end if
    blp % lastbox => p
    boxp => p % box
    boxp = ' '
  end subroutine newbox

  function tostr(a) result(s)
    ! Copy char array to string
    character, intent(in) :: a(:)
    character(size(a)) s
    integer i
    do i=1,size(a)
      s(i:i) = a(i)
    enddo
  end function tostr

  subroutine dispboxlist(unit, sep)
    ! Display the list of boxes associated with unit
    integer, intent(in) :: unit
    type(boxnode), pointer :: pfirst, p
    type(boxlist), pointer :: blp
    integer k, nlines, h, w, ns
    character(*), intent(in) :: sep
    blp => getboxlist(unit)
    pfirst => blp % firstbox
    nlines = 0
    p => pfirst
    do while (associated(p))
      nlines = max(nlines, size(p % box, 2))
      p => p % nextbox
    enddo
    do k=1,nlines
      p => pfirst
      ns = 0
      do while (associated(p))
        h = size(p % box, 2)
        w = size(p % box, 1)
        if (k <= h) then
          select case(unit)
          case(-1)
            continue
          case(-2)
            call putstr(sep(1:ns) // tostr(p % box(:,k)))
          case(-3)
            write(*,    '(2A)', advance = 'no') sep(1:ns), tostr(p % box(:,k))
          case default
            write(unit, '(2A)', advance = 'no') sep(1:ns), tostr(p % box(:,k))
          end select
        else
          select case(unit)
          case(-1)
            continue
          case(-2)
            call putstr(sep(1:ns) // repeat(' ', w))
          case(-3)
            write(*,    '(2A)', advance = 'no') sep(1:ns), repeat(' ', w)
          case default
            write(unit, '(2A)', advance = 'no') sep(1:ns), repeat(' ', w)
          end select
        end if
        p => p % nextbox
        ns = len(sep)
      enddo
      call dispnewline(unit)
    enddo
    call clearboxlist(unit)
  end subroutine dispboxlist

  subroutine dispnewline(unit)
    integer, intent(in) :: unit
    select case(unit)
    case(-1); continue
    case(-2); call putnl
    case(-3); write(*,*)
    case default; write(unit,*)
    end select
  end subroutine dispnewline

  !   subroutine print_boxes
  !     ! Print info on all boxes (used for debug purposes)
  !     integer :: k
  !     type(boxlist), pointer :: bl
  !     type(boxnode), pointer :: p
  !     bl => firstboxlist
  !     write(*,'("BOXES:")')
  !     do while (associated(bl))
  !       write(*,'("UNIT=",SS,I0,":")') bl % unit
  !       p => bl % firstbox
  !       k = 1
  !       do while(associated(p))
  !         write(*,'("  box ",SS,I0,", size=(",I0,",",I0,")")') k, shape(p % box)
  !         k = k+1
  !         p => p % nextbox
  !       enddo
  !       bl => bl % nextboxlist
  !     enddo
  !   end subroutine print_boxes

  ! ******************************** END OF BOX-PACKAGE *******************************

END MODULE DISPMODULE_UTIL








MODULE DISPMODULE
  use dispmodule_util
  implicit none

  PRIVATE  ! Everything not explicitly declared public should be private (including entities from dispmodule_util)

  ! ********************************** PUBLIC DECLARATIONS *************************************

  PUBLIC DISP                 ! Main routine of package, "pretty-prints" vectors and matrices
  PUBLIC DISP_SET             ! Subroutine to change default settings for DISP
  PUBLIC DISP_GET             ! Obtain current default settings
  PUBLIC DISP_SET_FACTORY     ! Call (without parameters) to restore original default settings
  PUBLIC TOSTRING             ! Convert numbers to strings
  PUBLIC TOSTRING_SET         ! Change settings for tostring
  PUBLIC TOSTRING_SET_FACTORY ! Restore original default settings for tostring
  !
  PUBLIC DISP_SETTINGS        ! Derived type with settings
  !
  PUBLIC ASTERISK_UNIT        ! Constant to specify displaying on asterisk unit (normally the screen)
  PUBLIC PUTSTR_UNIT          ! Constant to specify the use of subroutines putstr and putnl to display
  PUBLIC NULL_UNIT            ! Constant to specify discarding of all displayed output

  ! ********************************** INTERFACE DECLARATIONS *************************************
  interface disp_set
    module procedure disp_set, disp_set_ds
  end interface

  interface disp
    module procedure disp_s_dint, disp_ts_dint, disp_v_dint, disp_tv_dint, disp_m_dint, disp_tm_dint
    module procedure disp_s_sngl, disp_ts_sngl, disp_v_sngl, disp_tv_sngl, disp_m_sngl, disp_tm_sngl
    module procedure disp_s_dble, disp_ts_dble, disp_v_dble, disp_tv_dble, disp_m_dble, disp_tm_dble
    module procedure disp_s_cplx, disp_ts_cplx, disp_v_cplx, disp_tv_cplx, disp_m_cplx, disp_tm_cplx
    module procedure disp_s_cpld, disp_ts_cpld, disp_v_cpld, disp_tv_cpld, disp_m_cpld, disp_tm_cpld
    module procedure disp_s_dlog, disp_ts_dlog, disp_v_dlog, disp_tv_dlog, disp_m_dlog, disp_tm_dlog
    module procedure              disp_ts_dchr, disp_v_dchr, disp_tv_dchr, disp_m_dchr, disp_tm_dchr
  end interface

  interface tostring
    module procedure tostring_dint, tostring_f_dint, tostring_s_dint, tostring_sf_dint
    module procedure tostring_dlog, tostring_f_dlog, tostring_s_dlog, tostring_sf_dlog
    module procedure tostring_sngl, tostring_f_sngl, tostring_s_sngl, tostring_sf_sngl
    module procedure tostring_dble, tostring_f_dble, tostring_s_dble, tostring_sf_dble
    module procedure tostring_cplx, tostring_f_cplx, tostring_s_cplx, tostring_sf_cplx
    module procedure tostring_cpld, tostring_f_cpld, tostring_s_cpld, tostring_sf_cpld
  end interface

  ! *********************** DEFINITION OF TYPED CONSTANTS: UNITS AND KIND PARAMETERS ********************
  integer, parameter ::    &
       ASTERISK_UNIT = -3  ,&
       PUTSTR_UNIT   = -2  ,&
       NULL_UNIT     = -1

  integer, parameter :: dint = kind(0)       ! default integer
  integer, parameter :: sngl = kind(0.0)     ! single precision (default real)
  integer, parameter :: dble = kind(0d0)     ! double precision
  integer, parameter :: dlog = kind(.false.) ! default logical

  ! The above are also used as specific procedure (i.e. module procedure) name extensions, together 
  ! with the following:
  !        cplx = complex single precision (default complex)
  !        cpld = complex double precision

CONTAINS

  ! ******************************* SETTING AND GETTING PROCEDURES *************************************
  subroutine disp_set(advance, digmax, matsep, orient, sep, style, unit, zeroas)
    ! Change display settings according to individual parameters
    character(*), optional, intent(in) :: advance, sep, matsep, orient, style, zeroas
    integer, optional, intent(in) :: digmax, unit
    if (present(advance))    DEFSET % advance = upper(advance)
    if (present(sep))        DEFSET % sep = sep
    if (present(sep))        DEFSET % seplen = min(9, len(sep))
    if (present(zeroas))     DEFSET % zeroas = zeroas
    if (present(zeroas))     DEFSET % zaslen = min(9, len(zeroas))
    if (present(matsep))     DEFSET % matsep = matsep
    if (present(matsep))     DEFSET % matseplen = min(9, len(matsep))
    if (present(orient))     DEFSET % orient = upper(orient)
    if (present(style))      DEFSET % style = style
    if (present(digmax))     DEFSET % digmax = digmax
    if (present(unit))       DEFSET % unit = unit
    call check_settings
  end subroutine disp_set

  subroutine disp_set_factory()
    ! Change display settings to the original default
    DEFSET = FACTORY_SETTINGS
  end subroutine disp_set_factory

  subroutine avoid_compiler_warnings
    ! Routine that exists only to avoid compiler warnings (due to compiler bugs)
    type(boxlist), pointer :: boxl_dummy1 => null(), boxl_dummy2 => null()
    type(boxnode), pointer :: boxn_dummy1 => null(), boxn_dummy2 => null()
    type(tostring_settings), pointer :: ts1 => null(), ts2 => null()
    ts1 => ts2
    ts2 => ts1
    boxl_dummy2 => boxl_dummy1
    boxl_dummy1 => boxl_dummy2
    boxn_dummy2 => boxn_dummy1
    boxn_dummy1 => boxn_dummy2
  end subroutine avoid_compiler_warnings

  subroutine tostring_set(sep, rfmt, ifmt, trimb, trimz)
    character(*), optional, intent(in) :: sep, rfmt, ifmt, trimb, trimz
    if (present(sep))    tosset % sep    = upper(sep)
    if (present(sep))    tosset % seplen = min(9, len(sep))
    if (present(rfmt))   tosset % rfmt   = upper(rfmt)
    if (present(ifmt))   tosset % ifmt   = upper(ifmt)
    if (present(trimb))  tosset % trimb  = upper(trimb)
    if (present(trimz))  tosset % trimz  = upper(trimz)
    call tostring_check_settings
  end subroutine tostring_set

  subroutine tostring_set_factory()
    logical dummy
    dummy = .false.
    if (dummy) call avoid_compiler_warnings
    tosset = tosfac
  end subroutine tostring_set_factory

  subroutine disp_set_ds(settings)
    ! Change display settings according to the structure "settings"
    type(disp_settings), intent(in) :: settings
    DEFSET = settings
    call check_settings
  end subroutine disp_set_ds

  function disp_get() result(defs)
    ! Return current display settings
    type(disp_settings) :: defs
    defs = DEFSET
  end function disp_get

  ! ********************************* DEFAULT INTEGER PROCEDURES ****************************************
  subroutine disp_s_dint(x, fmt, advance, sep, trim, unit, zeroas)
    ! Default integer scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim, zeroas
    integer(dint), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_ts_dint('', x, fmt, advance, sep, 'left', trim, unit, zeroas)
  end subroutine disp_s_dint

  subroutine disp_v_dint(x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    ! Default integer vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    integer(dint), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tv_dint('', x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
  end subroutine disp_v_dint

  subroutine disp_m_dint(x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
    ! Default integer matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    integer(dint), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tm_dint('', x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
  end subroutine disp_m_dint

  subroutine disp_ts_dint(title, x, fmt, advance, sep, style, trim, unit, zeroas)
    ! Default integer scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    integer(dint), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_tm_dint(title, reshape((/x/), (/1, 1/)), fmt, advance, sep=sep, style=style, trim=trim, unit=unit, &
         zeroas=zeroas)
  end subroutine disp_ts_dint

  subroutine disp_tv_dint(title, x, fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    ! Default integer vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    integer(dint), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, zeroas)
    if (SE % row) then
      call disp_dint(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_dint(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_dint

  subroutine disp_tm_dint(title, x, fmt, advance, lbound, sep, style, trim, unit, zeroas)
    ! Default integer matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    integer(dint),intent(in)           :: x(:,:)     ! The matrix to be written
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
    call disp_dint(title, x, SE)
  end subroutine disp_tm_dint

  subroutine disp_dint(title, x, SE)
    ! Default integer item
    character(*),   intent(in)    :: title
    integer(dint),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer wid(size(x,2)), nbl(size(x,2))
    call find_editdesc_dint(x, SE, wid, nbl) ! determine also SE % w
    call tobox_dint(title, x, SE, wid, nbl)
  end subroutine disp_dint

  subroutine tobox_dint(title, x, SE, wid, nbl)
    ! Write default integer matrix to box
    character(*),   intent(in)    :: title
    integer(dint),  intent(in)    :: x(:,:)
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
  end subroutine tobox_dint

  subroutine find_editdesc_dint(x, SE, wid, nbl)
    ! Determine SE % ed, SE % w (unless specified) and wid
    integer(dint),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer,        intent(out)   :: wid(size(x,2)), nbl(size(x,2))
    !
    integer(dint) xmaxv(size(x,2)), xminv(size(x,2)), xp, xm
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
      xzero = any(x == 0_dint, 1) ! true where column has some zeros
      xallz = all(x == 0_dint, 1) ! true where column has only zeros
      call getwid_dint(xmaxv, xminv, xzero, xallz, SE,  wid, nbl)
    else
      wid = SE % w
      nbl = 0
    endif
  end subroutine find_editdesc_dint

  subroutine getwid_dint(xmaxv, xminv, xzero, xallz, SE,  wid, nbl)
    integer(dint),  intent(in)  :: xmaxv(:), xminv(:)
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
  end subroutine getwid_dint
  
  ! ********* DEFAULT INTEGER TOSTRING PROCEDURES *********
  function tostring_s_dint(x) result(st)
    ! Scalar to string
    integer(dint), intent(in)                   :: x
    character(len_f_dint((/x/), tosset % ifmt)) :: st
    st = tostring_f_dint((/x/), tosset % ifmt)
  end function tostring_s_dint

  function tostring_sf_dint(x, fmt) result(st)
    ! Scalar with specified format to string
    integer(dint),intent(in)        :: x
    character(*), intent(in)        :: fmt
    character(len_f_dint((/x/), fmt)) :: st
    st = tostring_f_dint((/x/), fmt)
  end function tostring_sf_dint

  function tostring_dint(x) result(st)
    ! Vector to string
    integer(dint), intent(in)               :: x(:)
    character(len_f_dint(x, tosset % ifmt)) :: st
    st = tostring_f_dint(x, tosset % ifmt)
  end function tostring_dint

  function tostring_f_dint(x, fmt) result(st)
    ! Vector with specified format to string
    integer(dint), intent(in)        :: x(:)
    character(*), intent(in)         :: fmt
    character(len_f_dint(x, fmt))    :: st
    character(widthmax_dint(x, fmt)) :: sa(size(x))
    integer                          :: w, d
    logical                          :: gedit
    character(nnblk(fmt)+5)          :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; st = errormsg; return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    call tostring_get(sa, st)
  end function tostring_f_dint

  pure function len_f_dint(x, fmt) result(wtot)
    ! Total width of tostring representation of x
    integer(dint), intent(in)        :: x(:)
    character(*), intent(in)         :: fmt
    character(widthmax_dint(x, fmt)) :: sa(size(x))
    integer                          :: wtot, w, d
    logical                          :: gedit
    character(nnblk(fmt)+5)          :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES' .or. w == 0) sa = adjustl(sa)
    wtot = sum(len_trim(sa)) + (size(x) - 1)*(tosset % seplen)
  end function len_f_dint

  pure function widthmax_dint(x, fmt) result(w)
    ! Maximum width of string representation of an element in x
    integer(dint), intent(in)  :: x(:)
    character(*), intent(in) :: fmt
    character(range(x)+2) sx(2)
    integer w, d
    logical gedit
    character(nnblk(fmt)+5) :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w<=0) then
      write(sx, '(SS,I0)') maxval(x), minval(x)
      w = maxval(len_trim(sx))
    end if
  end function widthmax_dint
  ! ************************************* END OF DEFAULT INTEGER PROCEDURES ******************************************

  ! **************************************** SINGLE PRECISION PROCEDURES *********************************************
  subroutine disp_s_sngl(x, fmt, advance, digmax, sep, trim, unit, zeroas)
    ! Single precision scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim, zeroas
    real(sngl), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_ts_sngl('', x, fmt, advance, digmax, sep, 'left', trim, unit, zeroas)
  end subroutine disp_s_sngl

  subroutine disp_v_sngl(x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
    ! Single precision vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    real(sngl), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    call disp_tv_sngl('', x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
  end subroutine disp_v_sngl

  subroutine disp_m_sngl(x, fmt, advance, lbound, sep, style, trim, unit, digmax, zeroas)
    ! Single precision matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    real(sngl), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, digmax, lbound(:)
    call disp_tm_sngl('', x, fmt, advance, digmax, lbound, sep, style, trim, unit, zeroas)
  end subroutine disp_m_sngl

  subroutine disp_ts_sngl(title, x, fmt, advance, digmax, sep, style, trim, unit, zeroas)
    ! Single precision scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    real(sngl), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_tm_sngl(title, reshape((/x/), (/1, 1/)), fmt, advance, digmax, sep=sep, style=style, trim=trim, &
         unit=unit, zeroas=zeroas)
  end subroutine disp_ts_sngl

  subroutine disp_tv_sngl(title, x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
    ! Single precision vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    real(sngl), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, zeroas, digmax)
    if (SE % row) then
      call disp_sngl(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_sngl(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_sngl
  
  subroutine disp_tm_sngl(title, x, fmt, advance, digmax, lbound, sep, style, trim, unit, zeroas)
    ! Single precision matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    real(sngl),   intent(in)           :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Editdit descriptor to use for each matrix element (e.g. 'F5.2')
    integer,      intent(in), optional :: unit       ! Unit to display on
    integer,      intent(in), optional :: digmax     ! Nbr of significant digits for largest abs value in x
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: zeroas     ! Zeros are replaced with this string if it is not empty
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no 
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    type(settings) :: SE
    !
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, zeroas=zeroas, digmax=digmax)
    call disp_sngl(title, x, SE)
  end subroutine disp_tm_sngl

  subroutine disp_sngl(title, x, SE)
    ! Single precision item
    character(*),   intent(in)    :: title
    real(sngl),     intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer wid(size(x,2)), nbl(size(x,2))
    call find_editdesc_sngl(x, SE, wid, nbl) ! determine also SE % w
    call tobox_sngl(title, x, SE, wid, nbl)
  end subroutine disp_sngl

  subroutine tobox_sngl(title, x, SE, wid, nbl)
    ! Write single precision matrix to box
    character(*),   intent(in)    :: title   ! title
    real(sngl),     intent(in)    :: x(:,:)  ! item
    type(settings), intent(inout) :: SE      ! settings
    integer,        intent(inout) :: wid(:)  ! widths of columns
    integer,        intent(inout) :: nbl(:)  ! number of blanks to trim from left
    character(SE % w)  :: s(size(x,1))
    integer            :: lin1, j, wleft, m, n, widp(size(wid))
    character, pointer :: boxp(:,:)
    real(sngl)         :: xj(size(x,1)), h
    m = size(x,1)
    n = size(x,2)
    h = huge(x)
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      xj = x(:, j)
      if (m > 0) write(s, SE % ed) xj
      call replace_zeronaninf(s, SE % zas(1:SE % lzas), xj == 0, xj /= xj, xj < -h, xj > h)
      call copytobox(s, lin1, wid(j), widp(j), nbl(j), boxp,  wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_sngl

  pure function maxw_sngl(x, d) result(w)
    ! Find max field width needed (F0.d editing is specified)
    real(sngl), intent(in) :: x(:)
    integer, intent(in) :: d
    integer expmax, expmin, w
    logical xfinite(size(x))
    real(sngl) xmax, xmin, h
    character(12) :: f1, s(2)
    xmin = 0; xmax = 0; h = huge(h)
    xfinite = x == x .and. x >= -h .and. x <= h ! neither NaN, Inf nor -Inf
    if (.not. any(xfinite)) then
      w = 4
    else
      xmax = maxval(x, mask=xfinite)
      xmin = minval(x, mask=xfinite)
      f1 = '(SS,ES9.0E4)'
      write(s,f1) xmax, xmin
      read(s(:)(5:9),'(I5)') expmax, expmin
      w = max(0, expmax, expmin) + d + 4
    end if
    if (.not. all(xfinite)) w = max(w, 4)
  end function maxw_sngl

  subroutine find_editdesc_sngl(x, SE, wid, nbl)
    ! Determine SE % ed, SE % w (unless specified) and wid.
    ! The if-block (*) is for safety: make f wider in case xm is written ok with the
    ! ES format in fmt but overflows with F format (the feature has been tested through
    ! manual changes to the program).
    real(sngl),     intent(in)    :: x(:,:)         ! Item to be written
    type(settings), intent(inout) :: SE             ! Settings
    integer,        intent(out)   :: wid(size(x,2)) ! Widths of individual columns
    integer,        intent(out)   :: nbl(size(x,2)) ! Blanks to trim from left of individual columns
    integer :: expmax, expmin, ww, dd, dmx
    real(sngl) xmaxv(size(x,2)), xminv(size(x,2)), xp, xm, h
    character(14) :: f1 = '(SS,ESxx.xxE4)'  ! could be ES99.89E4; default is ES14.05E4
    character(99) s
    logical xzero(size(x,2)), xallz(size(x,2)), xfinite(size(x,1),size(x,2)), xnonn(size(x,2)), xalln(size(x,2))
    !
    dmx = SE % dmx
    h = huge(h)
    xfinite = x == x .and. x >= -h .and. x <= h ! neither NaN, Inf nor -Inf
    if (SE % w == 0) then  ! Edit descriptor 'F0.d' specified
      ww = maxw_sngl(reshape(x, (/size(x)/)), SE % d)
      if (SE % lzas > 0 .and. any(x == 0._sngl))  ww = max(ww, SE % lzas)
      call replace_w(SE % ed, ww)
      SE % w = ww
    elseif (SE % w < 0) then ! No edit descriptor specified
      if (size(x) == 0) then
        SE % w = 0
        wid = 0
        nbl = 0
        return
      endif
      if (any(xfinite)) then
        xp = maxval(x, mask=xfinite)
        xm = minval(x, mask=xfinite)
        write(f1(7:11), '(SS,I2,".",I2.2)') dmx + 8, dmx - 1
        write(s,f1) xp; read(s(dmx+4:dmx+8),'(I5)') expmax
        write(s,f1) xm; read(s(dmx+4:dmx+8),'(I5)') expmin
        call find_editdesc_real(expmax, expmin, dmx,  SE % ed, ww, dd, xm >= 0)
        if (.not. all(xfinite))                     ww = max(ww, 4)
        if (SE % lzas > 0 .and. any(x == 0._sngl))  ww = max(ww, SE % lzas)
        if (SE % ed(5:5)=='F') then  ! (*)
          write(s, SE % ed) xp; if (s(1:1) == '*') ww = ww + 1
          write(s, SE % ed) xm; if (s(1:1) == '*') ww = ww + 1
          write(SE % ed(6:10), '(SS,I2,".",I2)') ww, dd
        endif
      else
        ww = 4
        SE % ed = '(F4.0)'
      endif
      SE % w = ww
    endif
    if (SE % trm) then
      xmaxv = maxval(x, 1, mask=xfinite)  ! max in each column
      xminv = minval(x, 1, mask=xfinite)  ! min
      xzero = any(x == 0._sngl, 1) ! true where column has some zeros
      xallz = all(x == 0._sngl, 1) ! true where column has only zeros
      xnonn = any(x > h .or. x < -h .or. x /= x, 1)  ! true where column has some nonnormals (inf, -inf, nan)
      xalln = all(x > h .or. x < -h .or. x /= x, 1)  ! true where column has only nonnormals (inf, -inf, nan)
      call getwid_sngl(xmaxv, xminv, xzero, xallz, xnonn, xalln, SE,  wid, nbl)
    else
      wid = SE % w
      nbl = 0
    endif
  end subroutine find_editdesc_sngl

  subroutine getwid_sngl(xmaxv, xminv, xzero, xallz, xnonn, xalln, SE,  wid, nbl)
    ! determine length of the strings that result when writing with edit descriptor SE%ed a 
    ! vector v where v(i) is xmaxv(i) or xminv(i) depending on which gives longer output
    real(sngl),     intent(in)  :: xmaxv(:), xminv(:) ! max and min values in each column
    logical,        intent(in)  :: xzero(:), xallz(:) ! true for columns with some/all zeros
    logical,        intent(in)  :: xnonn(:), xalln(:) ! true for columns with some/all nonnormals
    type(settings), intent(in)  :: SE                 ! settings
    integer,        intent(out) :: wid(:)             ! widths of columns
    integer,        intent(out) :: nbl(:)             ! number of blanks to peel from left (w-wid)
    character(SE % w) :: stmax(size(xmaxv)), stmin(size(xmaxv))
    integer w
    w = SE % w
    write(stmin, SE % ed) xminv
    write(stmax, SE % ed) xmaxv
    nbl = mod(verify(stmin, ' ') + w, w + 1) ! loc. of first nonblank
    nbl = min(nbl, mod(verify(stmax, ' ') + w, w + 1))
    if (SE % gedit) then
      wid = w
    else
      wid = len_trim(adjustl(stmin))
      wid = max(wid, len_trim(adjustl(stmax)))
    endif
    if (SE % lzas > 0) then
      wid = merge(SE % lzas, wid, xallz)
      wid = max(wid, merge(SE % lzas, 0, xzero))
    endif
    wid = merge(4, wid, xalln)
    wid = max(wid, merge(4, 0, xnonn))
    nbl = w - wid
  end subroutine getwid_sngl

  ! ******** TOSTRING SINGLE PRECISION PROCEDURES ***********
  function tostring_s_sngl(x) result(st)
    ! Scalar to string
    real(sngl), intent(in) :: x
    character(len_f_sngl((/x/), tosset % rfmt)) :: st
    st = tostring_f_sngl((/x/), tosset % rfmt)
  end function tostring_s_sngl

  function tostring_sf_sngl(x, fmt) result(st)
    ! Scalar with specified format to string
    real(sngl),   intent(in) :: x
    character(*), intent(in) :: fmt
    character(len_f_sngl((/x/), fmt)) :: st
    st = tostring_f_sngl((/x/), fmt)
  end function tostring_sf_sngl

  function tostring_sngl(x) result(st)
    ! Vector to string
    real(sngl), intent(in) :: x(:)
    character(len_f_sngl(x, tosset % rfmt)) :: st
    st = tostring_f_sngl(x, tosset % rfmt)
  end function tostring_sngl

  function tostring_f_sngl(x, fmt) result(st)
    ! Vector with specified format to string
    real(sngl)    ,       intent(in) :: x(:)
    character(*),         intent(in) :: fmt
    character(len_f_sngl(x, fmt))    :: st
    character(widthmax_sngl(x, fmt)) :: sa(size(x))
    character(nnblk(fmt)+8)          :: fmt1  !(5 for readfmt and 3 for replace_w)
    integer                          :: w, d, ww
    logical                          :: gedit
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then
      st = errormsg
      return
    elseif (w == 0) then
      ww = maxw_sngl(x, d)    
      call replace_w(fmt1, ww)      
    endif
    write(sa, fmt1) x
    call trim_real(sa, gedit, w)
    call tostring_get(sa, st)
  end function tostring_f_sngl

  pure function len_f_sngl(x, fmt) result(wtot)
    ! Total length of returned string, vector s
    real(sngl), intent(in)           :: x(:)
    character(*), intent(in)         :: fmt
    character(widthmax_sngl(x, fmt)) :: sa(size(x))
    integer                          :: wtot, w, d, ww
    logical                          :: gedit
    character(nnblk(fmt)+8)          :: fmt1  !(5 for readfmt and 3 for replace_w)
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    if (w == 0) then
      ww = maxw_sngl(x, d)
      call replace_w(fmt1, ww)
    endif
    write(sa, fmt1) x
    call trim_real(sa, gedit, w)
    wtot = sum(len_trim(sa)) + (size(x) - 1)*(tosset % seplen)
  end function len_f_sngl

  pure function widthmax_sngl(x, fmt) result(w)
    ! Maximum width of an element of x
    real(sngl), intent(in)   :: x(:)
    character(*), intent(in) :: fmt
    character(nnblk(fmt)+5)  :: fmt1
    integer w, d
    logical gedit
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then ! illegal format, use 1
      w = 1
    elseif (w == 0) then
      w = maxw_sngl(x, d)
    endif
  end function widthmax_sngl

  ! *************************************** END OF SINGLE PRECISION PROCEDURES ***************************************

  ! *************************************** SINGLE PRECISION COMPLEX PROCEDURES **************************************
  subroutine disp_s_cplx(x, fmt, fmt_imag, advance, digmax, sep, trim, unit)
    ! single precision complex scalar without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, trim
    complex(sngl), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_ts_cplx('', x, fmt, fmt_imag, advance, digmax, sep, 'left', trim, unit)
  end subroutine disp_s_cplx

  subroutine disp_v_cplx(x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
    ! single precision complex vector without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim, orient
    complex(sngl), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    call disp_tv_cplx('', x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
  end subroutine disp_v_cplx

  subroutine disp_m_cplx(x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
    ! single precision complex matrix without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim
    complex(sngl), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, digmax, lbound(:)
    call disp_tm_cplx('', x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
  end subroutine disp_m_cplx

  subroutine disp_ts_cplx(title, x, fmt, fmt_imag, advance, digmax, sep, style, trim, unit)
    ! single precision complex scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim
    complex(sngl), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_tm_cplx(title, reshape((/x/), (/1, 1/)), fmt, fmt_imag, advance, digmax, sep=sep, style=style, &
                                                       trim=trim, unit=unit)
  end subroutine disp_ts_cplx

  subroutine disp_tv_cplx(title, x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
    ! single precision complex vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim, orient
    complex(sngl), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    type(settings) SE, SEim
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, digmax=digmax)
    if (present(fmt_imag)) then
      if (.not.present(fmt)) then
        call disp_errmsg('DISP: error, FMT must be present if FMT_IMAG is present'); return; 
      endif
      call get_SE(SEim, title, shape(x), fmt_imag)
    else
      SEim = SE
    end if
    if (SE % row) then
      call disp_cplx(title, reshape(x, (/1, size(x)/)), SE, SEim, n = size(x))
    else
      call disp_cplx(title, reshape(x, (/size(x), 1/)), SE, SEim, n = 1)
    end if
  end subroutine disp_tv_cplx

  subroutine disp_tm_cplx(title, x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
    ! single precision complex matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    complex(sngl),  intent(in)         :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Edit descriptor for each element (real element when fmt_imag &
    !                                                ! is present)
    character(*), intent(in), optional :: fmt_imag   ! Edit descriptor for each imaginary element
    integer,      intent(in), optional :: unit       ! Unit to display on
    integer,      intent(in), optional :: digmax     ! Nbr of significant digits for largest abs value in real(x) &
    !                                                ! and aimag(x)
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no 
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    !
    type(settings) :: SE, SEim
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, digmax=digmax)
    if (present(fmt_imag)) then
      if (.not.present(fmt)) then
        call disp_errmsg('DISP: error, FMT must be present if FMT_IMAG is present'); return
      endif
      call get_SE(SEim, title, shape(x), fmt_imag)
    else
      SEim = SE
    end if
    call disp_cplx(title, x, SE, SEim, n = size(x,2))
  end subroutine disp_tm_cplx

  subroutine disp_cplx(title, x, SE, SEim, n)
    ! Single precision item
    character(*),   intent(in)    :: title
    complex(sngl),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE, SEim
    integer,        intent(in)    :: n
    integer, dimension(n) :: widre(n), widim(n), nblre(n), nblim(n)
    call find_editdesc_sngl(real(x), SE, widre, nblre)         ! determine also SE % w
    call find_editdesc_sngl(abs(aimag(x)), SEim, widim, nblim) ! determine also SEim % w
    call tobox_cplx(title, x, SE, SEim, widre, widim, nblre, nblim, m = size(x,1), n = size(x,2))
  end subroutine disp_cplx

  subroutine tobox_cplx(title, x, SE, SEim, widre, widim, nblre, nblim, m, n)
    ! Write single precision complex matrix to box
    character(*),   intent(in)    :: title
    complex(sngl),  intent(in)    :: x(:,:)
    integer,        intent(in)    :: m, n, widre(:), widim(:), nblre(:), nblim(:)
    type(settings), intent(inout) :: SE, SEim
    character(SE % w)   :: s(m)
    character(SEim % w) :: sim(m)
    character(3)        :: sgn(m)
    integer             :: lin1, i, j, wleft, wid(n), widp(n)
    character, pointer  :: boxp(:,:)
    SE % zas = ''
    SEim % zas = ''
    wid = widre + widim + 4
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      if (m > 0) write(s, SE % ed) (real(x(i,j)), i=1,m)
      call copytobox(s, lin1, widre(j), widp(j) - widim(j) - 4, nblre(j), boxp,  wleft)
      do i=1,m
        if (aimag(x(i,j)) < 0) then; sgn(i) = ' - '; else; sgn(i) = ' + '; endif
        enddo
      call copytobox(sgn, lin1, 3, 3, 0, boxp,  wleft)
      if (m > 0) write(sim, SEim % ed) (abs(aimag(x(i,j))), i=1,m)
      call copytobox(sim, lin1, widim(j), widim(j), nblim(j), boxp,  wleft)
      call copyseptobox('i', m, lin1, boxp, wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_cplx

  ! ******* TOSTRING SINGLE PRECISION COMPLEX PROCEDURES ********

  function tostring_s_cplx(x) result(st)
    complex(sngl), intent(in)                   :: x
    character(len_s_cplx(x, tosset % rfmt)) :: st
    st = tostring_f_cplx((/x/), tosset % rfmt)
  end function tostring_s_cplx

  function tostring_sf_cplx(x, fmt) result(st)
    complex(sngl),  intent(in)        :: x
    character(*), intent(in)          :: fmt
    character(len_s_cplx(x, fmt)) :: st
    st = tostring_f_cplx((/x/), fmt)
  end function tostring_sf_cplx

  function tostring_cplx(x) result(st)
    complex(sngl), intent(in)               :: x(:)
    character(len_f_cplx(x, tosset % rfmt)) :: st
    st = tostring_f_cplx(x, tosset % rfmt)
  end function tostring_cplx

  function tostring_f_cplx(x, fmt) result(st)
    complex(sngl),  intent(in)                    :: x(:)
    character(*),   intent(in)                    :: fmt
    character(len_f_cplx(x, fmt))                 :: st
    character(widthmax_sngl(real(x), fmt))        :: sar(size(x))
    character(widthmax_sngl(abs(x-real(x)), fmt)) :: sai(size(x))  ! x-real(x) instead of aimag(x) to enable the fnction
    character(1)                                  :: sgn(size(x))  ! to pass -stand:f95 switch of the ifort compiler.
    integer                                       :: w, d, wr, wi, i
    logical                                       :: gedit
    character(nnblk(fmt)+8)                       :: fmt1  !(5 for readfmt and 3 for replace_w)
    real(sngl)                                    :: xre(size(x)), xim(size(x)), h
    call readfmt(fmt, fmt1, w, d, gedit)
    xre = real(x)
    xim = aimag(x)
    h = huge(h)
    if (w < 0) then
      st = errormsg
      return
    elseif (w == 0) then
      wr = maxw_sngl(xre, d)
      wi = maxw_sngl(xim, d)
      call replace_w(fmt1, max(wr, wi))
    endif
    write(sar, fmt1) real(x)
    write(sai, fmt1) abs(aimag(x))
    call trim_real(sar, gedit, w)
    call trim_real(sai, gedit, w)
    do i = 1,size(x); if (aimag(x(i)) < 0) then; sgn(i) = '-'; else; sgn(i) = '+'; endif; enddo
    call tostring_get_complex(sar, sgn, sai, st)
  end function tostring_f_cplx

  pure function len_s_cplx(x, fmt) result(wtot)
    complex(sngl), intent(in) :: x
    character(*), intent(in)  :: fmt
    integer                   :: wtot, w, d
    logical                   :: gedit
    character(nnblk(fmt)+8)   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    wtot = len_f_sngl((/real(x)/), fmt) + len_f_sngl((/abs(aimag(x))/), fmt) + 4
  end function len_s_cplx

  pure function len_f_cplx(x, fmt) result(wtot)
    complex(sngl), intent(in) :: x(:)
    character(*), intent(in)  :: fmt
    integer                   :: wtot, w, d
    logical                   :: gedit
    character(nnblk(fmt)+8)   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    wtot = len_f_sngl(real(x), fmt) + len_f_sngl(abs(aimag(x)), fmt) + size(x)*4 - (size(x) - 1)*(tosset % seplen)
    ! subtract seplen because it has been added twice in len_f_sngl
  end function len_f_cplx
  ! *************************************** END OF SINGLE PRECISION COMPLEX PROCEDURES ********************************

  ! ************************************* DOUBLE PRECISION PROCEDURES (SEE NOTE 2 BELOW) ******************************
  subroutine disp_s_dble(x, fmt, advance, digmax, sep, trim, unit, zeroas)
    ! Double precision scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim, zeroas
    real(dble), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_ts_dble('', x, fmt, advance, digmax, sep, 'left', trim, unit, zeroas)
  end subroutine disp_s_dble

  subroutine disp_v_dble(x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
    ! Double precision vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    real(dble), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    call disp_tv_dble('', x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
  end subroutine disp_v_dble

  subroutine disp_m_dble(x, fmt, advance, lbound, sep, style, trim, unit, digmax, zeroas)
    ! Double precision matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    real(dble), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, digmax, lbound(:)
    call disp_tm_dble('', x, fmt, advance, digmax, lbound, sep, style, trim, unit, zeroas)
  end subroutine disp_m_dble

  subroutine disp_ts_dble(title, x, fmt, advance, digmax, sep, style, trim, unit, zeroas)
    ! Double precision scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas
    real(dble), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_tm_dble(title, reshape((/x/), (/1, 1/)), fmt, advance, digmax, sep=sep, style=style, trim=trim, &
         unit=unit, zeroas=zeroas)
  end subroutine disp_ts_dble

  subroutine disp_tv_dble(title, x, fmt, advance, digmax, lbound, sep, style, trim, unit, orient, zeroas)
    ! Double precision vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, zeroas, orient
    real(dble), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, zeroas, digmax)
    if (SE % row) then
      call disp_dble(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_dble(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_dble
  
  subroutine disp_tm_dble(title, x, fmt, advance, digmax, lbound, sep, style, trim, unit, zeroas)
    ! Double precision matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    real(dble),   intent(in)           :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Editdit descriptor to use for each matrix element (e.g. 'F5.2')
    integer,      intent(in), optional :: unit       ! Unit to display on
    integer,      intent(in), optional :: digmax     ! Nbr of significant digits for largest abs value in x
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: zeroas     ! Zeros are replaced with this string if it is not empty
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no 
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    type(settings) :: SE
    !
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, zeroas=zeroas, digmax=digmax)
    call disp_dble(title, x, SE)
  end subroutine disp_tm_dble

  subroutine disp_dble(title, x, SE)
    ! Double precision item
    character(*),   intent(in)    :: title
    real(dble),     intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE
    integer wid(size(x,2)), nbl(size(x,2))
    call find_editdesc_dble(x, SE, wid, nbl) ! determine also SE % w
    call tobox_dble(title, x, SE, wid, nbl)
  end subroutine disp_dble

  subroutine tobox_dble(title, x, SE, wid, nbl)
    ! Write double precision matrix to box
    character(*),   intent(in)    :: title   ! title
    real(dble),     intent(in)    :: x(:,:)  ! item
    type(settings), intent(inout) :: SE      ! settings
    integer,        intent(inout) :: wid(:)  ! widths of columns
    integer,        intent(inout) :: nbl(:)  ! number of blanks to trim from left
    character(SE % w)  :: s(size(x,1))
    integer            :: lin1, j, wleft, m, n, widp(size(wid))
    character, pointer :: boxp(:,:)
    real(dble)         :: xj(size(x,1)), h
    m = size(x,1)
    n = size(x,2)
    h = huge(x)
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      xj = x(:, j)
      if (m > 0) write(s, SE % ed) xj
      call replace_zeronaninf(s, SE % zas(1:SE % lzas), xj == 0, xj /= xj, xj < -h, xj > h)
      call copytobox(s, lin1, wid(j), widp(j), nbl(j), boxp,  wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_dble

  pure function maxw_dble(x, d) result(w)
    ! Find max field width needed (F0.d editing is specified)
    real(dble), intent(in) :: x(:)
    integer, intent(in) :: d
    integer expmax, expmin, w
    logical xfinite(size(x))
    real(dble) xmax, xmin, h
    character(12) :: f1, s(2)
    xmin = 0; xmax = 0; h = huge(h)
    xfinite = x == x .and. x >= -h .and. x <= h ! neither NaN, Inf nor -Inf
    if (.not. any(xfinite)) then
      w = 4
    else
      xmax = maxval(x, mask=xfinite)
      xmin = minval(x, mask=xfinite)
      f1 = '(SS,ES9.0E4)'
      write(s,f1) xmax, xmin
      read(s(:)(5:9),'(I5)') expmax, expmin
      w = max(0, expmax, expmin) + d + 4
    end if
    if (.not. all(xfinite)) w = max(w, 4)
  end function maxw_dble

  subroutine find_editdesc_dble(x, SE, wid, nbl)
    ! Determine SE % ed, SE % w (unless specified) and wid.
    ! The if-block (*) is for safety: make f wider in case xm is written ok with the
    ! ES format in fmt but overflows with F format (the feature has been tested through
    ! manual changes to the program).
    real(dble),     intent(in)    :: x(:,:)         ! Item to be written
    type(settings), intent(inout) :: SE             ! Settings
    integer,        intent(out)   :: wid(size(x,2)) ! Widths of individual columns
    integer,        intent(out)   :: nbl(size(x,2)) ! Blanks to trim from left of individual columns
    integer :: expmax, expmin, ww, dd, dmx
    real(dble) xmaxv(size(x,2)), xminv(size(x,2)), xp, xm, h
    character(14) :: f1 = '(SS,ESxx.xxE4)'  ! could be ES99.89E4; default is ES14.05E4
    character(99) s
    logical xzero(size(x,2)), xallz(size(x,2)), xfinite(size(x,1),size(x,2)), xnonn(size(x,2)), xalln(size(x,2))
    !
    dmx = SE % dmx
    h = huge(h)
    xfinite = x == x .and. x >= -h .and. x <= h ! neither NaN, Inf nor -Inf
    if (SE % w == 0) then  ! Edit descriptor 'F0.d' specified
      ww = maxw_dble(reshape(x, (/size(x)/)), SE % d)
      if (SE % lzas > 0 .and. any(x == 0._dble))  ww = max(ww, SE % lzas)
      call replace_w(SE % ed, ww)
      SE % w = ww
    elseif (SE % w < 0) then ! No edit descriptor specified
      if (size(x) == 0) then
        SE % w = 0
        wid = 0
        nbl = 0
        return
      endif
      if (any(xfinite)) then
        xp = maxval(x, mask=xfinite)
        xm = minval(x, mask=xfinite)
        write(f1(7:11), '(SS,I2,".",I2.2)') dmx + 8, dmx - 1
        write(s,f1) xp; read(s(dmx+4:dmx+8),'(I5)') expmax
        write(s,f1) xm; read(s(dmx+4:dmx+8),'(I5)') expmin
        call find_editdesc_real(expmax, expmin, dmx,  SE % ed, ww, dd, xm >= 0)
        if (.not. all(xfinite))                     ww = max(ww, 4)
        if (SE % lzas > 0 .and. any(x == 0._dble))  ww = max(ww, SE % lzas)
        if (SE % ed(5:5)=='F') then  ! (*)
          write(s, SE % ed) xp; if (s(1:1) == '*') ww = ww + 1
          write(s, SE % ed) xm; if (s(1:1) == '*') ww = ww + 1
          write(SE % ed(6:10), '(SS,I2,".",I2)') ww, dd
        endif
      else
        ww = 4
        SE % ed = '(F4.0)'
      endif
      SE % w = ww
    endif
    if (SE % trm) then
      xmaxv = maxval(x, 1, mask=xfinite)  ! max in each column
      xminv = minval(x, 1, mask=xfinite)  ! min
      xzero = any(x == 0._dble, 1) ! true where column has some zeros
      xallz = all(x == 0._dble, 1) ! true where column has only zeros
      xnonn = any(x > h .or. x < -h .or. x /= x, 1)  ! true where column has some nonnormals (inf, -inf, nan)
      xalln = all(x > h .or. x < -h .or. x /= x, 1)  ! true where column has only nonnormals (inf, -inf, nan)
      call getwid_dble(xmaxv, xminv, xzero, xallz, xnonn, xalln, SE,  wid, nbl)
    else
      wid = SE % w
      nbl = 0
    endif
  end subroutine find_editdesc_dble

  subroutine getwid_dble(xmaxv, xminv, xzero, xallz, xnonn, xalln, SE,  wid, nbl)
    ! determine length of the strings that result when writing with edit descriptor SE%ed a 
    ! vector v where v(i) is xmaxv(i) or xminv(i) depending on which gives longer output
    real(dble),     intent(in)  :: xmaxv(:), xminv(:) ! max and min values in each column
    logical,        intent(in)  :: xzero(:), xallz(:) ! true for columns with some/all zeros
    logical,        intent(in)  :: xnonn(:), xalln(:) ! true for columns with some/all nonnormals
    type(settings), intent(in)  :: SE                 ! settings
    integer,        intent(out) :: wid(:)             ! widths of columns
    integer,        intent(out) :: nbl(:)             ! number of blanks to peel from left (w-wid)
    character(SE % w) :: stmax(size(xmaxv)), stmin(size(xmaxv))
    integer w
    w = SE % w
    write(stmin, SE % ed) xminv
    write(stmax, SE % ed) xmaxv
    nbl = mod(verify(stmin, ' ') + w, w + 1) ! loc. of first nonblank
    nbl = min(nbl, mod(verify(stmax, ' ') + w, w + 1))
    if (SE % gedit) then
      wid = w
    else
      wid = len_trim(adjustl(stmin))
      wid = max(wid, len_trim(adjustl(stmax)))
    endif
    if (SE % lzas > 0) then
      wid = merge(SE % lzas, wid, xallz)
      wid = max(wid, merge(SE % lzas, 0, xzero))
    endif
    wid = merge(4, wid, xalln)
    wid = max(wid, merge(4, 0, xnonn))
    nbl = w - wid
  end subroutine getwid_dble

  ! ******** TOSTRING DOUBLE PRECISION PROCEDURES ***********
  function tostring_s_dble(x) result(st)
    ! Scalar to string
    real(dble), intent(in) :: x
    character(len_f_dble((/x/), tosset % rfmt)) :: st
    st = tostring_f_dble((/x/), tosset % rfmt)
  end function tostring_s_dble

  function tostring_sf_dble(x, fmt) result(st)
    ! Scalar with specified format to string
    real(dble),   intent(in) :: x
    character(*), intent(in) :: fmt
    character(len_f_dble((/x/), fmt)) :: st
    st = tostring_f_dble((/x/), fmt)
  end function tostring_sf_dble

  function tostring_dble(x) result(st)
    ! Vector to string
    real(dble), intent(in) :: x(:)
    character(len_f_dble(x, tosset % rfmt)) :: st
    st = tostring_f_dble(x, tosset % rfmt)
  end function tostring_dble

  function tostring_f_dble(x, fmt) result(st)
    ! Vector with specified format to string
    real(dble)    ,       intent(in) :: x(:)
    character(*),         intent(in) :: fmt
    character(len_f_dble(x, fmt))    :: st
    character(widthmax_dble(x, fmt)) :: sa(size(x))
    character(nnblk(fmt)+8)          :: fmt1  !(5 for readfmt and 3 for replace_w)
    integer                          :: w, d, ww
    logical                          :: gedit
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then
      st = errormsg
      return
    elseif (w == 0) then
      ww = maxw_dble(x, d)    
      call replace_w(fmt1, ww)      
    endif
    write(sa, fmt1) x
    call trim_real(sa, gedit, w)
    call tostring_get(sa, st)
  end function tostring_f_dble

  pure function len_f_dble(x, fmt) result(wtot)
    ! Total length of returned string, vector s
    real(dble), intent(in)           :: x(:)
    character(*), intent(in)         :: fmt
    character(widthmax_dble(x, fmt)) :: sa(size(x))
    integer                          :: wtot, w, d, ww
    logical                          :: gedit
    character(nnblk(fmt)+8)          :: fmt1  !(5 for readfmt and 3 for replace_w)
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    if (w == 0) then
      ww = maxw_dble(x, d)
      call replace_w(fmt1, ww)
    endif
    write(sa, fmt1) x
    call trim_real(sa, gedit, w)
    wtot = sum(len_trim(sa)) + (size(x) - 1)*(tosset % seplen)
  end function len_f_dble

  pure function widthmax_dble(x, fmt) result(w)
    ! Maximum width of an element of x
    real(dble), intent(in)   :: x(:)
    character(*), intent(in) :: fmt
    character(nnblk(fmt)+5)  :: fmt1
    integer w, d
    logical gedit
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then ! illegal format, use 1
      w = 1
    elseif (w == 0) then
      w = maxw_dble(x, d)
    endif
  end function widthmax_dble

  ! *************************************** END OF DOUBLE PRECISION PROCEDURES ***************************************

  ! *************************************** DOUBLE PRECISION COMPLEX PROCEDURES **************************************
  subroutine disp_s_cpld(x, fmt, fmt_imag, advance, digmax, sep, trim, unit)
    ! double precision complex scalar without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, trim
    complex(dble), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_ts_cpld('', x, fmt, fmt_imag, advance, digmax, sep, 'left', trim, unit)
  end subroutine disp_s_cpld

  subroutine disp_v_cpld(x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
    ! double precision complex vector without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim, orient
    complex(dble), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    call disp_tv_cpld('', x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
  end subroutine disp_v_cpld

  subroutine disp_m_cpld(x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
    ! double precision complex matrix without title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim
    complex(dble), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, digmax, lbound(:)
    call disp_tm_cpld('', x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
  end subroutine disp_m_cpld

  subroutine disp_ts_cpld(title, x, fmt, fmt_imag, advance, digmax, sep, style, trim, unit)
    ! double precision complex scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim
    complex(dble), intent(in) :: x
    integer, intent(in), optional :: unit, digmax
    call disp_tm_cpld(title, reshape((/x/), (/1, 1/)), fmt, fmt_imag, advance, digmax, sep=sep, style=style, &
                                                       trim=trim, unit=unit)
  end subroutine disp_ts_cpld

  subroutine disp_tv_cpld(title, x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit, orient)
    ! double precision complex vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, fmt_imag, advance, sep, style, trim, orient
    complex(dble), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:), digmax
    type(settings) SE, SEim
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient, digmax=digmax)
    if (present(fmt_imag)) then
      if (.not.present(fmt)) then
        call disp_errmsg('DISP: error, FMT must be present if FMT_IMAG is present'); return; 
      endif
      call get_SE(SEim, title, shape(x), fmt_imag)
    else
      SEim = SE
    end if
    if (SE % row) then
      call disp_cpld(title, reshape(x, (/1, size(x)/)), SE, SEim, n = size(x))
    else
      call disp_cpld(title, reshape(x, (/size(x), 1/)), SE, SEim, n = 1)
    end if
  end subroutine disp_tv_cpld

  subroutine disp_tm_cpld(title, x, fmt, fmt_imag, advance, digmax, lbound, sep, style, trim, unit)
    ! double precision complex matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    complex(dble),  intent(in)         :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Edit descriptor for each element (real element when fmt_imag &
    !                                                ! is present)
    character(*), intent(in), optional :: fmt_imag   ! Edit descriptor for each imaginary element
    integer,      intent(in), optional :: unit       ! Unit to display on
    integer,      intent(in), optional :: digmax     ! Nbr of significant digits for largest abs value in real(x) &
    !                                                ! and aimag(x)
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: style      ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no 
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    !
    type(settings) :: SE, SEim
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, digmax=digmax)
    if (present(fmt_imag)) then
      if (.not.present(fmt)) then
        call disp_errmsg('DISP: error, FMT must be present if FMT_IMAG is present'); return
      endif
      call get_SE(SEim, title, shape(x), fmt_imag)
    else
      SEim = SE
    end if
    call disp_cpld(title, x, SE, SEim, n = size(x,2))
  end subroutine disp_tm_cpld

  subroutine disp_cpld(title, x, SE, SEim, n)
    ! Double precision item
    character(*),   intent(in)    :: title
    complex(dble),  intent(in)    :: x(:,:)
    type(settings), intent(inout) :: SE, SEim
    integer,        intent(in)    :: n
    integer, dimension(n) :: widre(n), widim(n), nblre(n), nblim(n)
    call find_editdesc_dble(real(x), SE, widre, nblre)         ! determine also SE % w
    call find_editdesc_dble(abs(aimag(x)), SEim, widim, nblim) ! determine also SEim % w
    call tobox_cpld(title, x, SE, SEim, widre, widim, nblre, nblim, m = size(x,1), n = size(x,2))
  end subroutine disp_cpld

  subroutine tobox_cpld(title, x, SE, SEim, widre, widim, nblre, nblim, m, n)
    ! Write double precision complex matrix to box
    character(*),   intent(in)    :: title
    complex(dble),  intent(in)    :: x(:,:)
    integer,        intent(in)    :: m, n, widre(:), widim(:), nblre(:), nblim(:)
    type(settings), intent(inout) :: SE, SEim
    character(SE % w)   :: s(m)
    character(SEim % w) :: sim(m)
    character(3)        :: sgn(m)
    integer             :: lin1, i, j, wleft, wid(n), widp(n)
    character, pointer  :: boxp(:,:)
    SE % zas = ''
    SEim % zas = ''
    wid = widre + widim + 4
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      if (m > 0) write(s, SE % ed) (real(x(i,j)), i=1,m)
      call copytobox(s, lin1, widre(j), widp(j) - widim(j) - 4, nblre(j), boxp,  wleft)
      do i=1,m
        if (aimag(x(i,j)) < 0) then; sgn(i) = ' - '; else; sgn(i) = ' + '; endif
        enddo
      call copytobox(sgn, lin1, 3, 3, 0, boxp,  wleft)
      if (m > 0) write(sim, SEim % ed) (abs(aimag(x(i,j))), i=1,m)
      call copytobox(sim, lin1, widim(j), widim(j), nblim(j), boxp,  wleft)
      call copyseptobox('i', m, lin1, boxp, wleft)
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine tobox_cpld

  ! ******* TOSTRING DOUBLE PRECISION COMPLEX PROCEDURES ********

  function tostring_s_cpld(x) result(st)
    complex(dble), intent(in)                   :: x
    character(len_s_cpld(x, tosset % rfmt)) :: st
    st = tostring_f_cpld((/x/), tosset % rfmt)
  end function tostring_s_cpld

  function tostring_sf_cpld(x, fmt) result(st)
    complex(dble),  intent(in)        :: x
    character(*), intent(in)          :: fmt
    character(len_s_cpld(x, fmt)) :: st
    st = tostring_f_cpld((/x/), fmt)
  end function tostring_sf_cpld

  function tostring_cpld(x) result(st)
    complex(dble), intent(in)               :: x(:)
    character(len_f_cpld(x, tosset % rfmt)) :: st
    st = tostring_f_cpld(x, tosset % rfmt)
  end function tostring_cpld

  function tostring_f_cpld(x, fmt) result(st)
    complex(dble),  intent(in)                    :: x(:)
    character(*),   intent(in)                    :: fmt
    character(len_f_cpld(x, fmt))                 :: st
    character(widthmax_dble(real(x), fmt))        :: sar(size(x))
    character(widthmax_dble(abs(x-real(x)), fmt)) :: sai(size(x))  ! x-real(x) instead of aimag(x) to enable the fnction
    character(1)                                  :: sgn(size(x))  ! to pass -stand:f95 switch of the ifort compiler.
    integer                                       :: w, d, wr, wi, i
    logical                                       :: gedit
    character(nnblk(fmt)+8)                       :: fmt1  !(5 for readfmt and 3 for replace_w)
    real(dble)                                    :: xre(size(x)), xim(size(x)), h
    call readfmt(fmt, fmt1, w, d, gedit)
    xre = real(x)
    xim = aimag(x)
    h = huge(h)
    if (w < 0) then
      st = errormsg
      return
    elseif (w == 0) then
      wr = maxw_dble(xre, d)
      wi = maxw_dble(xim, d)
      call replace_w(fmt1, max(wr, wi))
    endif
    write(sar, fmt1) real(x)
    write(sai, fmt1) abs(aimag(x))
    call trim_real(sar, gedit, w)
    call trim_real(sai, gedit, w)
    do i = 1,size(x); if (aimag(x(i)) < 0) then; sgn(i) = '-'; else; sgn(i) = '+'; endif; enddo
    call tostring_get_complex(sar, sgn, sai, st)
  end function tostring_f_cpld

  pure function len_s_cpld(x, fmt) result(wtot)
    complex(dble), intent(in) :: x
    character(*), intent(in)  :: fmt
    integer                   :: wtot, w, d
    logical                   :: gedit
    character(nnblk(fmt)+8)   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    wtot = len_f_dble((/real(x)/), fmt) + len_f_dble((/abs(aimag(x))/), fmt) + 4
  end function len_s_cpld

  pure function len_f_cpld(x, fmt) result(wtot)
    complex(dble), intent(in) :: x(:)
    character(*), intent(in)  :: fmt
    integer                   :: wtot, w, d
    logical                   :: gedit
    character(nnblk(fmt)+8)   :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w < 0) then; wtot = len(errormsg); return; endif
    wtot = len_f_dble(real(x), fmt) + len_f_dble(abs(aimag(x)), fmt) + size(x)*4 - (size(x) - 1)*(tosset % seplen)
    ! subtract seplen because it has been added twice in len_f_dble
  end function len_f_cpld
  ! *************************************** END OF DOUBLE PRECISION COMPLEX PROCEDURES ********************************

  ! ********************************************** DEFAULT LOGICAL PROCEDURES *****************************************
  subroutine disp_s_dlog(x, fmt, advance, sep, trim, unit)
    ! Default logical scalar without title
    character(*), intent(in), optional :: fmt, advance, sep, trim
    logical(dlog), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_ts_dlog('', x, fmt, advance, sep, 'left', trim, unit)
  end subroutine disp_s_dlog

  subroutine disp_v_dlog(x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! Default logical vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    logical(dlog), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tv_dlog('', x, fmt, advance, lbound, sep, style, trim, unit, orient)
  end subroutine disp_v_dlog

  subroutine disp_m_dlog(x, fmt, advance, lbound, sep, style, trim, unit)
    ! Default logical matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim
    logical(dlog), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tm_dlog('', x, fmt, advance, lbound, sep, style, trim, unit)
  end subroutine disp_m_dlog

  subroutine disp_ts_dlog(title, x, fmt, advance, sep, style, trim, unit)
    ! Default logical scalar with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim
    logical(dlog), intent(in) :: x
    integer, intent(in), optional :: unit
    call disp_tm_dlog(title, reshape((/x/), (/1, 1/)), fmt, advance, sep=sep, style=style, trim=trim, unit=unit)
  end subroutine disp_ts_dlog

  subroutine disp_tv_dlog(title, x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! Default logical vector with title
    character(*), intent(in) :: title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    logical(dlog), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient)
    if (SE % row) then
      call disp_dlog(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_dlog(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_dlog

  subroutine disp_tm_dlog(title, x, fmt, advance, lbound, sep, style, trim, unit)
    ! Default logical matrix with title
    character(*), intent(in)           :: title     ! The title to use for the matrix
    logical(dlog),intent(in)           :: x(:,:)    ! The matrix to be written
    character(*), intent(in), optional :: fmt       ! Format edit descriptor to use for each matrix element (e.g. 'L1')
    integer,      intent(in), optional :: unit      ! Unit to display on
    character(*), intent(in), optional :: advance   ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep       ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: style     ! Style(s): See NOTE 1 below
    character(*), intent(in), optional :: trim      ! 'Auto' (the default) to trim if fmt absent, 'no' for no trimming,
    !                                               ! 'yes' for trimming
    integer,      intent(in), optional :: lbound(:) ! Lower bounds of x
    type(settings) :: SE
    !
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit)
    call disp_dlog(title, x, SE)
  end subroutine disp_tm_dlog

  subroutine disp_dlog(title, x, SE)
    ! Write default logical to box or unit
    character(*),   intent(in)    :: title
    logical(dlog),  intent(in)    :: x(:,:)
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
    call tobox_dlog(title, x, SE, wid, nbl)
  end subroutine disp_dlog

  subroutine tobox_dlog(title, x, SE, wid, nbl)
    character(*),   intent(in)    :: title
    logical(dlog),  intent(in)    :: x(:,:)
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
  end subroutine tobox_dlog

  ! ********** DEFAULT LOGICAL TOSTRING PROCEDURES *********
  function tostring_s_dlog(x) result(st)
    logical(dlog), intent(in) :: x
    character(1)            :: st
    st = tostring_f_dlog((/x/), 'L1')
  end function tostring_s_dlog

  function tostring_sf_dlog(x, fmt) result(st)
    logical(dlog),intent(in)        :: x
    character(*), intent(in)        :: fmt
    character(len_f_dlog((/x/), fmt)) :: st
    st = tostring_f_dlog((/x/), fmt)
  end function tostring_sf_dlog

  function tostring_dlog(x) result(st)
    logical(dlog), intent(in)                          :: x(:)
    character(1 + (size(x) - 1)*(1 + tosset % seplen)) :: st
    st = tostring_f_dlog(x, 'L1')
  end function tostring_dlog

  function tostring_f_dlog(x, fmt) result(st)
    logical(dlog), intent(in)     :: x(:)
    character(*), intent(in)      :: fmt
    character(len_f_dlog(x, fmt)) :: st
    character(widthmax_dlog(fmt)) :: sa(size(x))
    integer                       :: w, d
    logical                       :: gedit
    character(nnblk(fmt)+5)       :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) then; st = errormsg; return; endif
    write(sa, fmt1) x
    if (tosset % trimb == 'YES') sa = adjustl(sa)
    call tostring_get(sa, st)
  end function tostring_f_dlog

  pure function len_f_dlog(x, fmt) result(wtot)
    logical(dlog), intent(in)  :: x(:)
    character(*), intent(in)   :: fmt
    integer                    :: wtot, w, d
    logical                    :: gedit
    character(nnblk(fmt)+5)    :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) then; wtot = len(errormsg); return; endif
    if (tosset % trimb == 'YES') wtot = size(x)
    if (tosset % trimb == 'NO' ) wtot = w*size(x)
    wtot = wtot + (size(x) - 1)*(tosset % seplen)
  end function len_f_dlog

  pure function widthmax_dlog(fmt) result(w)
    character(*), intent(in) :: fmt
    integer w, d
    logical gedit
    character(nnblk(fmt)+5) :: fmt1
    call readfmt(fmt, fmt1, w, d, gedit)
    if (w <= 0) w = 1
  end function widthmax_dlog
  ! ****************************** END OF DEFAULT LOGICAL PROCEDURES *******************************

  ! ******************************* DEFAULT CHARACTER PROCEDURES **********************************
  subroutine disp_v_dchr(x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! Default character vector without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    character(*), intent(in) :: x(:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tv_dchr('', x, fmt, advance, lbound, sep, style, trim, unit, orient)
  end subroutine disp_v_dchr

  subroutine disp_m_dchr(x, fmt, advance, lbound, sep, style, trim, unit)
    ! Default character matrix without title
    character(*), intent(in), optional :: fmt, advance, sep, style, trim
    character(*), intent(in) :: x(:,:)
    integer, intent(in), optional :: unit, lbound(:)
    call disp_tm_dchr('', x, fmt, advance, lbound, sep, style, trim, unit)
  end subroutine disp_m_dchr

  subroutine disp_ts_dchr(title, x, fmt, advance, sep, style, trim, unit)
    ! Default character scalar with title
    character(*), intent(in), optional :: title, x, fmt, advance, sep, style, trim
    character(0) empty(1,0)
    integer, intent(in), optional :: unit
    empty = ''
    if (present(title).and.present(x)) then
      call disp_nonopt_dchr(title, x, fmt, advance, sep=sep, style=style, trim=trim, unit=unit)
    elseif (present(x)) then
      call disp_nonopt_dchr('', x, fmt, advance, sep=sep, style='left', trim=trim, unit=unit)
    elseif (present(title)) then
      call disp_nonopt_dchr('', title, fmt, advance, sep=sep, style='left', trim=trim, unit=unit)
    else
      call disp_tm_dchr('', empty, fmt, advance, sep=sep, style=style, trim=trim, unit=unit)
    end if
  end subroutine disp_ts_dchr

  subroutine disp_nonopt_dchr(title, x, fmt, advance, sep, style, trim, unit)
    ! This routine exists to circumvent bug in gfortran, that made it not possible to change scalar strings
    ! to matrices with reshape in calls of disp_tm_dchr. This intermediate routine provides work-around.
    character(*), intent(in) :: title, x, fmt, advance, sep, style, trim
    optional fmt, advance, sep, style, trim
    integer, intent(in), optional :: unit
    character(len(x)) :: xm(1,1)
    xm(1,1) = x
    call disp_tm_dchr(title, xm, fmt, advance, sep=sep, style=style, trim=trim, unit=unit)
  end subroutine disp_nonopt_dchr

  subroutine disp_tv_dchr(title, x, fmt, advance, lbound, sep, style, trim, unit, orient)
    ! Default character vector with title
    character(*), intent(in) :: title, x(:)
    character(*), intent(in), optional :: fmt, advance, sep, style, trim, orient
    integer, intent(in), optional :: unit, lbound(:)
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit, orient)
    if (SE % row) then
      call disp_dchr(title, reshape(x, (/1, size(x)/)), SE)
    else
      call disp_dchr(title, reshape(x, (/size(x), 1/)), SE)
    end if
  end subroutine disp_tv_dchr

  subroutine disp_tm_dchr(title, x, fmt, advance, lbound, sep, style, trim, unit)
    ! Default character matrix with title
    character(*), intent(in)           :: title      ! The title to use for the matrix
    character(*), intent(in)           :: x(:,:)     ! The matrix to be written
    character(*), intent(in), optional :: fmt        ! Format edit descriptor to use for each matrix element (e.g.'A4')
    integer,      intent(in), optional :: unit       ! Unit to display on
    character(*), intent(in), optional :: advance    ! 'No' to print next matrix to right of current, otherewise 'Yes'
    character(*), intent(in), optional :: sep        ! Separator between matrix columns (e.g. ", ")
    character(*), intent(in), optional :: style      ! Style(s): see NOTE 1 below
    character(*), intent(in), optional :: trim       ! 'Auto' (the default) to trim if fmt absent, 'no' for no
    !                                                ! trimming, 'yes' for trimming
    integer,      intent(in), optional :: lbound(:)  ! Lower bounds of x
    !
    type(settings) :: SE
    call get_SE(SE, title, shape(x), fmt, advance, lbound, sep, style, trim, unit)
    call disp_dchr(title, x, SE)
  end subroutine disp_tm_dchr

  subroutine disp_dchr(title, x, SE)
    ! Default character item to box 
    character(*), intent(in)      :: title, x(:,:)
    type(settings), intent(inout) :: SE
    character(13)                 :: edesc
    character, pointer            :: boxp(:,:)
    integer                       :: m, n, j, lin1, wleft, lx, w
    integer, dimension(size(x,2)) :: wid, nbl, n1, n2, widp
    m = size(x,1)
    n = size(x,2)
    lx = len(x)
    w = SE % w
    if (w <= 0) then
      w = lx
      if (w < 0) then
        edesc = '(A__________)'
        write(edesc(3:12), '(SS,I10)') w
        SE % ed = edesc
      end if
    end if
    if (SE % trm .and. size(x) > 0) then
      n1 = minval(mod(verify(x, ' ') - w - 1, w + 1), 1) + w + 1
      n2 = maxval(verify(x, ' ', back = .true.), 1)
      wid = n2 - n1 + 1
      nbl = w - wid
    else
      n1 = 1
      n2 = w
      wid = w
      nbl = 0
    end if
    if (all(wid == 0)) n = 0
    SE % w = w
    call preparebox(title, SE, m, n, wid, widp, lin1, wleft, boxp)
    do j=1,n
      if (SE % trm) then
        call copytobox(x(:,j)(n1(j):n2(j)), lin1, wid(j), widp(j), nbl(j), boxp,  wleft)
      else
        if (widp(j) > lx) call copyseptobox(repeat(' ', widp(j)-lx), m, lin1, boxp,  wleft)
        call copytobox(x(:,j), lin1, lx, lx, 0, boxp,  wleft)
      end if
      if (j<n) call copyseptobox(SE % sep(1:SE % lsep), m, lin1, boxp,  wleft)
    enddo
    call finishbox(title, SE, boxp)
  end subroutine disp_dchr

  ! ************************* END OF DEFAULT CHARACTER PROCEDURES ********************************

  ! NOTE 1: STYLES
  !   Styles can be LEFT, ABOVE, PAD, UNDERLINE or NUMBER. Padding is by default done with hyphen
  !   characters (e.g. ---title---), but can be changed for example to asterisks with style='*PAD'.
  !   Underlining is also with hypens and can also be changed, e.g. with style='*UNDERLINE'. Lower
  !   or mixed case is acceptable: style='above' or style='Above'. It is also possible to specify
  !   both NUMBER and one of the other styles, with e.g. style='ABOVE & NUMBER'.
  !
  ! NOTE 2: DOUBLE PRECISION
  !   The double precision functions and subroutines above (the sections marked DOUBLE PRECISION
  !   PROECDURES and DOUBLE PRECISION COMPLEX PROECEDURES) are copies of the sections marked SINGLE
  !   PRECISION PROCEDURES and SINGLE PRECISION COMPLEX PROCEDURES, with the kind parameter sngl
  !   changed to dble, the procedure name suffixes _sngl and _cplx changed to _dble and _cpld, and
  !   single changed to double (only appears in comments). The add-on module DISP_R16MOD is another
  !   copy of these procedures (for quad precision).
  
END MODULE DISPMODULE
