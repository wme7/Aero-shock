MODULE TEST_NANINF_MOD

  implicit none
  integer, parameter :: verbose = 0  ! 0 = quiet
  !                                  ! 1 = report only names of test routines
  !                                  ! 2 = report also what the tests display 
  integer, parameter :: srk = kind(0d0)
  integer :: compare_no = 0
  character(30) :: assert_string = ''
  character(3)  :: adv = 'no'
  character(90) :: fmt = '("  Testing display of nan and -inf and inf, real kind=",I0,"...")'

CONTAINS

  subroutine test_ni(nan, inf, minf)
    USE DISPMODULE
    !use, intrinsic :: ieee_arithmetic
    real(srk), intent(in) :: nan, inf, minf
    real(srk) :: &
         zero = 0._srk, &
         d1(3,2) = reshape((/-2.146_srk, 0._srk, 1._srk, 0._srk, 20.33_srk, 0._srk/), (/3,2/)), &
         d2(2)   = (/0._srk, 0._srk/), &
         d3(2)   = (/12e20_srk, 0._srk/), &
         d4      = 0._srk
    logical :: mask1(3,2), mask2(2), mask3(2)
    character(12) :: s1(3)
    character(9)  :: s2(1)
    character(16) :: s3(2)
    character(9)  :: s4(2)
    !
  call open_8
    if (verbose > 0) adv = 'yes'
    write(*, fmt, advance = adv) srk

    mask1 = d1 == zero
    mask2 = d2 == zero
    mask3 = d3 == zero
    d1 = merge(nan, d1, mask1)
    d2 = merge(nan, d2, mask2)
    d3 = merge(nan, d3, mask3)
    d4 = nan
    !
    s1 = (/'-2.15    NaN' ,&
        &  '  NaN  20.33' ,&
        &  ' 1.00    NaN' /)
    !
    s2 = ' NaN, NaN'
    !
    s3(1) = 'A = 1.20000E+21'
    s3(2) = '            NaN'
    !
    s4(1) = 'Longtitle'
    s4(2) = '   NaN   '
    !
    call assert_init('TEST_NAN')
    call disp_set(unit=8, style = 'PAD')
    call disp(d1, digmax=4);                  call compare(s1, 'NaN-1')
    call disp(d1, digmax=4, zeroas='0.');     call compare(s1, 'NaN-1')
    call disp(d2, sep = ',', orient = 'row'); call compare(s2, 'NaN-2')
    call disp('A = ', d3, style = 'left');    call compare(s3, 'NaN-3')
    call disp('Longtitle', d4);               call compare(s4, 'NaN-4')

    d1 = merge(inf, d1, mask1)
    d2 = merge(inf, d2, mask2)
    d3 = merge(inf, d3, mask3)
    d4 = inf
    s1(1) = '-2.15   +Inf'
    s1(2) = ' +Inf  20.33'
    s1(3) = ' 1.00   +Inf'
    s2    = '+Inf,+Inf'
    s3(2) = '           +Inf'
    s4(2) = '  +Inf   '
    call disp(d1, digmax=4);                  call compare(s1, 'Inf-1')
    call disp(d1, digmax=4, zeroas='0.');     call compare(s1, 'Inf-1')
    call disp(d2, sep = ',', orient = 'row'); call compare(s2, 'Inf-2')
    call disp('A = ', d3, style = 'left');    call compare(s3, 'Inf-3')
    call disp('Longtitle', d4);               call compare(s4, 'Inf-4')

    d1(1,:) = merge(minf, d1(1,:), mask1(1,:))
    d1(2,:) = merge( nan, d1(2,:), mask1(2,:))
    d1(3,:) = merge( inf, d1(3,:), mask1(3,:))
    d2 = merge(minf, d2, mask2)
    d3 = merge(minf, d3, mask3)
    d4 = minf
    s1(1) = '-2.15   -Inf'
    s1(2) = '  NaN  20.33'
    s1(3) = ' 1.00   +Inf'
    s2    = '-Inf,-Inf'
    s3(2) = '           -Inf'
    s4(2) = '  -Inf   '
    call disp(d1, digmax=4);                  call compare(s1, '-Inf-1')
    call disp(d2, sep = ',', orient = 'row'); call compare(s2, '-Inf-2')
    call disp('A = ', d3, style = 'left');    call compare(s3, '-Inf-3')
    call disp('Longtitle', d4);               call compare(s4, '-Inf-4')

    s1(1) = ' -2.  -Inf'
    s1(2) = ' NaN   20.'
    s1(3) = '  1.  +Inf'
    call disp(d1, digmax=2);           call compare(s1, 'Mixed-1')
    call disp(d1, 'F4.0');             call compare(s1, 'Mixed-2')
    call disp(d1, 'F4.0', trim='yes'); call compare(s1, 'Mixed-3')
    s1(1) = '-2.  ***'
    s1(2) = 'NaN  20.'
    s1(3) = ' 1.  Inf'
    call disp(d1, 'F3.0'); call compare(s1, 'Field-overflow-1')
    
    call disp(d2, 'F2.0'); call compare((/'**','**'/), 'Field-overflow-2')

    call disp_set_factory
    write(*,'("  OK")')
    call close_8

  end subroutine test_ni

  subroutine msg1(st)
    ! Print st if verbose is >= 1
    character(*) st
    if (verbose >= 1) write(*,'(2x,a)') trim(st)
  end subroutine msg1

  subroutine msg2(st)
    ! Print st if verbose is >= 2
    character(*) st
    if (verbose >= 2) write(*,'(4x,a)') trim(st)
  end subroutine msg2

  subroutine compare(sok, message, sok1)
    ! Utility for all the test routines. Print to the screen what the last disp calls
    ! displayed, and assert that what was displayed matches either sok or sok1.
    character(*)       :: sok(:), sok1(:), message
    optional           :: message, sok1
    character(100)     :: s, mess
    integer ios1, i
    rewind(8)
    compare_no = compare_no + 1
    if (present(message)) then
      mess = message
    else
      write(mess, '(a,"-",i0)') trim(assert_string), compare_no
    endif
    do i = 1,9999
      read(8,1,iostat=ios1) s
      if (ios1 < 0) exit
      call msg2(s)
      call assert(i <= size(sok), message)
      if (present(sok1))       call assert(sok(i) == s .or. sok1(i) == s, mess)
      if (.not. present(sok1)) call assert(sok(i) == s, mess)
    enddo
    call assert(i == size(sok) + 1, mess)
    call reopen_8
    1 format(A)
  end subroutine compare

  subroutine open_8
    open(8, file = 'testtmp.dat', status = 'replace')
  end subroutine open_8
  
  subroutine close_8
    close(8)
  end subroutine close_8

  subroutine reopen_8
    call close_8
    call open_8
  end subroutine reopen_8

  subroutine assert_init(st)
    ! Set assert-string (used by subroutine compare) to st, compare_no to zero and display st
    ! if verbose is true. If st is absent, set assert_string to ''
    character(*), optional :: st
    if (present(st)) then
      assert_string = st
      call msg2('')
      call msg1(st)
    else
      assert_string = ''
    endif
    compare_no = 0
  end subroutine assert_init

  subroutine assert(s, msg)
    ! Assert that s is true. If not print "assertion failed" and msg if it is present
    logical s
    character(*), optional :: msg
    if (.not.s) then
      if (present(msg)) then
        print '(a, ": assertion failed")', trim(msg)
      else
        print '("assertion failed")'
      end if
      stop
    endif
  end subroutine assert

END MODULE TEST_NANINF_MOD
