module types
  implicit none

  integer, parameter:: dp=kind(0.d0)
end module



!==============================================================================
! interpolation function
!==============================================================================
module interp
  implicit none

  contains


    subroutine idxsearch(ar1, ar2, k1, k2, idx)
      ! =====================================================
      ! Called by linear_interp
      !   ar1 and ar2  assumed to be monotonically decreasing or increasing
      !   idx(k) means ar2(k) is between ar1(k) and ar1(k+1)
      !   idx(k) = 0 or size(ar1) indicates ar2(k) is outside the range of ar1
      ! =====================================================
      use types
      integer,intent(in)   :: k1, k2
      real(dp), intent(in) :: ar1(k1) !sequence value to search
      real(dp), intent(in) :: ar2(k2) !set of value to serach for
      integer, intent(out) :: idx(k2) !interval locations

      integer  :: i, j
      real(dp) :: val

      if (ar1(1) < ar1(k1)) then
        !monotonically increasing case
        do i = 1, k2
          val = ar2(i)
          if (val <= ar1(1)) then
            idx(i) = 0
            cycle
          else if ( val >= ar1(k1)) then
            idx(i) = k1
            cycle
          end if
          do j = 1, k1-1
            if (val <= ar1(j+1)) then
              idx(i) = j
              exit
            end if
          end do
        end do
      else
        !monotonically decreasing case
        do i = 1, k2
          val = ar2(i)
          if (val >= ar1(1)) then
            idx(i) = 0
            cycle
          else if ( val <= ar1(k1)) then
            idx(i) = k1
            cycle
          end if
          do j = 1, k1-1
            if (val >= ar1(j+1)) then
              idx(i) = j
              exit
            end if
          end do
        end do
      end if
    end subroutine idxsearch


    subroutine linear_interp(in,zin,zout,xnin,znin,znout,out)
      ! =====================================================
      ! 2D vertical linear interpolation.
      ! in   : (xnin, znin) values used for interpolation.
      ! zin  : vertical coordinates for in array.
      ! zout : vertical levels which interpolating applied.
      ! xnin,znin,znout : dimensions for in, zin and zout.
      ! =====================================================
      use types
      integer,  intent(in)  :: xnin, znin, znout
      real(dp), intent(in)  :: in(xnin,znin), zin(xnin,znin), zout(znout)
      real(dp), intent(out) :: out(xnin,znout)

      !f2py intent(in) :: in, zin, zout
      !f2py intent(out) :: out

      integer idx(znout)
      real(dp) :: a, x(znin), y(znin)
      integer  :: i, j, jo

      do i = 1, xnin
        x = zin(i,:)
        y = in(i,:)
        call idxsearch(x, zout, znin, znout, idx)
        do j = 1, znout
          jo = idx(j)
          if (jo == 0) then ! outside of bottom level is constant
            out(i,j) = y(1)
          else if (jo == znin) then !outside of top level is constant
            out(i,j) = y(znin)
          else !linear interpolation
            a = (y(jo+1) - y(jo))/(x(jo+1) - x(jo))
            out(i,j) = a*(zout(j) - x(jo)) + y(jo)
          end if
        end do
      end do
    end subroutine linear_interp


    function binarysearch(length, array, value, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      ! https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90
      use types
      implicit none
      integer, intent(in) :: length
      real(dp), dimension(length), intent(in) :: array
      real(dp), intent(in) :: value
      real(dp), intent(in), optional :: delta

      !f2py intent(in) :: length, array, value
      !f2py intent(in), optional :: delta
      !f2py depend(length) :: array

      integer  :: binarysearch
      integer  :: left, middle, right
      real(dp) :: d

      if (present(delta) .eqv. .true.) then
        d = delta
      else
        d = 1e-9
      endif

      left = 1
      right = length
      do
        if (left > right) then
          exit
        endif
        middle = nint((left+right) / 2.0)
        if ( abs(array(middle) - value) <= d) then
          binarySearch = middle
          return
        else if (array(middle) > value) then
          right = middle - 1
        else
          left = middle + 1
        end if
      end do
      binarysearch = right

    end function binarysearch


    subroutine bilinear_interp(x_array, y_array, f, x, y, z, x_len, y_len, p_len, delta)
      ! This function uses bilinear interpolation to estimate the value
      ! of a function f at points (x,y)
      ! f is assumed to be sampled on a regular grid, with the grid x values specified
      ! by x_array and the grid y values specified by y_array
      ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
      ! https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90
      use types
      implicit none
      integer, intent(in) :: x_len, y_len, p_len
      real(dp), dimension(x_len), intent(in) :: x_array
      real(dp), dimension(y_len), intent(in) :: y_array
      real(dp), dimension(x_len, y_len), intent(in) :: f
      real(dp), dimension(p_len), intent(in) :: x, y
      real(dp), dimension(p_len), intent(out) :: z
      real(dp), intent(in), optional :: delta
      !f2py depend(x_len) :: x_array, f
      !f2py depend(y_len) :: y_array, f
      !f2py depend(p_len) :: x, y, z

      real(dp) :: denom, x1, x2, y1, y2
      integer :: i, j, p

      do p=1, p_len
        i = binarysearch(x_len, x_array, x(p))
        j = binarysearch(y_len, y_array, y(p))

        x1 = x_array(i)
        x2 = x_array(i+1)

        y1 = y_array(j)
        y2 = y_array(j+1)

        denom = (x2 - x1)*(y2 - y1)

        z(p) = (f(i,j)*(x2-x(p))*(y2-y(p)) + f(i+1,j)*(x(p)-x1)*(y2-y(p)) + &
            f(i,j+1)*(x2-x(p))*(y(p)-y1) + f(i+1, j+1)*(x(p)-x1)*(y(p)-y1))/denom
      end do
    end subroutine bilinear_interp

end module interp


!==============================================================================
! grid manipulating function
!==============================================================================
module grid
  implicit none

  contains

    subroutine smooth_area_average(in_field, lon, lat, radius, out_field, n, nlon, nlat)
      ! =====================================================
      ! smoothing grid field with area average
      ! in_field  : 3D field array used to smoothing, the 2, 3 dimension are lat and lon.
      ! lon, lat  : 2D array lon and lat coorindates for the in_field.
      ! radius    : smoothing radius.
      ! out_field : return smoothed field.
      ! =====================================================
      use types
      real(dp), intent(in)  :: n, nlon, nlat
      real(dp), intent(in)  :: in_field(n, nlat, nlon), lon(nlat, nlon), lat(nlat, nlon)
      real(dp), intent(in)  :: radius
      real(dp), intent(out) :: out_field(n, nlat, nlon)

      !f2py intent(in) :: in_field, lon, lat, radius
      !f2py intent(out) :: out_field

      integer :: i, j, k, ii, jj
      real(dp) :: dlat(nlat, nlon), dlon(nlat, nlon), lat1, lat2(nlat, nlon)
      real(dp) :: a(nlat, nlon), dist(nlat, nlon), flag(nlat, nlon), numb
      real(dp), parameter :: earth_radius = 6371000.0
      real(dp), parameter :: deg_to_rad = atan(1.0)/45

      flag = 1.0
      do j=1, nlat
        do i=1, nlon
          ! great circle distance
          dlat = (lat-lat(j,i)) * deg_to_rad
          dlon = (lon-lon(j,i)) * deg_to_rad
          lat1 = lat(j,i) * deg_to_rad
          lat2 = lat * deg_to_rad
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          dist = earth_radius*2.0*asin(sqrt(a))
      
          numb = sum(flag, mask=dist <= radius)
          if (numb > 0) then
            do k=1, n
              out_field(k, j, i) = sum(in_field(k, :, :), mask=dist <= radius) / numb
            end do
          else
            out_field(:, j, i) = in_field(:, j, i)
          end if
        end do
      end do
    end subroutine smooth_area_average


    subroutine smooth_gaussian(in_field, lon, lat, radius, sigma, out_field, n, nlon, nlat)
      ! =====================================================
      ! smoothing grid field with area average
      ! in_field  : 3D field array used to smoothing, the 2, 3 dimension are lat and lon.
      ! lon, lat  : 2D array lon and lat coorindates for the in_field.
      ! radius    : smoothing radius.
      ! sigma     : gaussian sigma parameter.
      ! out_field : return smoothed field.
      ! =====================================================
      use types
      real(dp), intent(in) :: n, nlon, nlat
      real(dp), intent(in) :: in_field(n, nlat, nlon), lon(nlat, nlon), lat(nlat, nlon)
      real(dp), intent(in) :: radius, sigma    ! sigma is an adjustable smoothing length scale.
      real(dp), intent(out) :: out_field(n, nlat, nlon)
      !f2py intent(in) :: in_field, lon, lat, radius, sigma
      !f2py intent(out) :: out_field
  
      integer :: i, j, k, ii, jj
      real(dp) :: dlat(nlat, nlon), dlon(nlat, nlon), lat1, lat2(nlat, nlon)
      real(dp) :: a(nlat, nlon), dist(nlat, nlon), field(nlat, nlon), flag(nlat, nlon), numb
      real(dp), parameter :: earth_radius = 6371000.0
      real(dp), parameter :: deg_to_rad = atan(1.0)/45
      real(dp), parameter :: pi = 3.1415926535897932_8

      flag = 1.0
      do j=1, nlat
        do i=1, nlon
          ! great circle distance
          dlat = (lat-lat(j,i)) * deg_to_rad
          dlon = (lon-lon(j,i)) * deg_to_rad
          lat1 = lat(j,i) * deg_to_rad
          lat2 = lat * deg_to_rad
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          dist = earth_radius*2.0*asin(sqrt(a))
      
          numb = sum(flag, mask=dist <= radius)
          if (numb > 0) then
            do k=1, n
              ! gaussian filter value
              field = exp(-0.5*(dist/sigma)**2)*in_field(k, :, :)/(2*pi*sigma**2)
              out_field(k, j, i) = sum(field, mask=dist <= radius) / numb
            end do
          else
            out_field(:, j, i) = in_field(:, j, i)
          end if
        end do
      end do
    end subroutine smooth_gaussian

end module grid


!==============================================================================
! regridding function
!==============================================================================
module regrid
  implicit none

  contains

  subroutine box_average(field0, lon0, lat0, lon1, lat1, width, field1, nlon0, nlat0, nlon1, nlat1)
    ! =====================================================
    ! regridding high resolution field0 to coarse field1 with box_average.
    ! field0: 2D array (nlat1, nlon0) for input high resolution, with missing value is -9999.0
    ! lon0, lat0: 2D coordinates for field0
    ! field1: 2D array (nlat1, nlon0) for output coarse resolution.
    ! lon1, lat1: 2D coordinates for field1
    ! width: box width
    ! =====================================================
    use types
    real(dp), intent(in) :: width, nlon0, nlat0, nlon1, nlat1
    real(dp), intent(in) :: field0(nlat0, nlon0), lon0(nlat0, nlon0), lat0(nlat0, nlon0)
    real(dp), intent(in) :: lon1(nlat1, nlon1), lat1(nlat1, nlon1)
    real(dp), intent(out) :: field1(nlat1, nlon1)
    !f2py intent(in) :: field0, lon0, lat0, lon1, lat1, width
    !f2py intent(out) :: feild1

    integer :: i, j, k, ii, jj, numb
    logical :: mask(nlat0, nlon0)

    do j=1, nlat1
      do i=1, nlon1
        mask = (lon0 >= (lon1(j,i)-width)) .and. (lon0 <= (lon1(j,i)+width)) .and. &
               (lat0 >= (lat1(j,i)-width)) .and. (lat0 <= (lat1(j,i)+width)) .and. &
               (field0 /= -9999.0)
        numb = count(mask)
        if (numb > 0) then
          field1(j, i) = sum(field0, mask=mask) / numb
        else
          field1(j, i) = -9999.0
        end if
      end do
    end do
  end subroutine box_average
end module regrid

