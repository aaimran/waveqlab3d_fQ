!> This module contains some definitions that are used in every module.

module common

  integer, parameter :: sp = selected_real_kind(9,49) !< Short precision
  integer, parameter :: lp = selected_real_kind(15,99) !< Long precision
  integer, parameter :: wp = lp !< Working precision

end module
