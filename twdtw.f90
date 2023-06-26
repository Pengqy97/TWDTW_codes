! TWDTW algorithm
! compile with `gfortran ./twdtw.f90 -Ofast -shared -fPIC -o ./twdtw.so`
! Four versions:
! double64: Time series in real(8)，date series in interger(8)
! double32: Time series in real(8)，date series in interger(4)
! single64: Time series in real(4)，date series in interger(8)
! single32: Time series in real(4)，date series in interger(4)
! Original algorithm written in R (https://github.com/vwmaus/dtwSat)
! Converted to MATLAB by Jie Dong
! Converted to Fortran by Ruoque Shen (https://github.com/shenrq)
! VERSION: v2023.06.26

! TWDTW 核心函数
module dtw_core
   implicit none
   public

   interface dist
      module procedure :: s_dist, d_dist
   end interface dist

   interface tw_dist
      module procedure :: s_tw_dist, d_tw_dist
   end interface tw_dist

   interface tw
      module procedure :: tw_single32, tw_single64, tw_double32, tw_double64
   end interface tw

   interface core
      module procedure :: s_core, d_core
   end interface core

   contains

   subroutine s_dist(x, m, std, n, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j))
         end do
      end do
   end subroutine s_dist

   subroutine d_dist(x, m, std, n, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j))
         end do
      end do
   end subroutine d_dist

   subroutine d_tw_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) + tweight(i, j)
         end do
      end do
   end subroutine d_tw_dist

   subroutine s_tw_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) + tweight(i, j)
         end do
      end do
   end subroutine s_tw_dist

   subroutine tw_single32(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(4) :: time_diff, alpha, beta
      real(4), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_single32

   subroutine tw_double32(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(8) :: time_diff, alpha, beta
      real(8), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_double32

   subroutine tw_single64(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(4) :: time_diff, alpha, beta
      real(4), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_single64

   subroutine tw_double64(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(8) :: time_diff, alpha, beta
      real(8), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_double64

   real(4) function s_core(d, m, n)
      implicit none
      integer, intent(in)  :: m, n
      real(4) :: d(m, n)
      integer :: i, j
      do i = 2, m
         d(i, 1) = d(i, 1) + d(i-1, 1)
      end do
      do j = 2, n
         d(1, j) = d(1, j) + d(1, j-1)
         do i = 2, m
            d(i, j) = d(i, j) + min(d(i-1, j), d(i-1, j-1), d(i, j-1))
         end do
      end do
      s_core = d(m, n)
   end function s_core

   real(8) function d_core(d, m, n)
      implicit none
      integer, intent(in)  :: m, n
      real(8) :: d(m, n)
      integer :: i, j
      do i = 2, m
         d(i, 1) = d(i, 1) + d(i-1, 1)
      end do
      do j = 2, n
         d(1, j) = d(1, j) + d(1, j-1)
         do i = 2, m
            d(i, j) = d(i, j) + min(d(i-1, j), d(i-1, j-1), d(i, j-1))
         end do
      end do
      d_core = d(m, n)
   end function d_core

end module dtw_core


! TW and DTW step by step.
module tw_dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8) :: d(m, n)
      call tw_dist(x, m, std, n, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4) :: d(m, n)
      call tw_dist(x, m, std, n, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(8) function double16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(4) function single16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module tw_dtw

! TWDTW.
module twdtw
   use dtw_core
   use tw_dtw
   implicit none
   public
   contains

   real(8) function double64(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      double64 = doublex(x, m, std, n, tweight)
   end function double64

   real(8) function double32(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer, intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      double32 = doublex(x, m, std, n, tweight)
   end function double32


   real(4) function single64(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      single64 = singlex(x, m, std, n, tweight)
   end function single64

   real(4) function single32(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer, intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      single32 = singlex(x, m, std, n, tweight)
   end function single32

end module twdtw

! DTW.
module dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, n)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: d(m, n)
      call dist(x, m, std, n, d)
      doublex = core(d, m, n)
   end function doublex


   real(4) function singlex(x, m, std, n)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: d(m, n)
      call dist(x, m, std, n, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in) :: x(m)
      real(8), intent(in)  :: std(n), factor
      real(8) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      double8 = core(d, m, n)
   end function double8

   real(4) function single8(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in) :: x(m)
      real(4), intent(in)  :: std(n), factor
      real(4) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      single8 = core(d, m, n)
   end function single8

   real(8) function double16(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in) :: x(m)
      real(8), intent(in)  :: std(n), factor
      real(8) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single16(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in) :: x(m)
      real(4), intent(in)  :: std(n), factor
      real(4) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      single16 = core(d, m, n)
   end function single16

end module dtw
