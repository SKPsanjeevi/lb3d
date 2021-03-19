AC_DEFUN([ACX_C_RESTRICT],
[AC_CACHE_CHECK([for C restrict keyword], acx_cv_c_restrict,
  [acx_cv_c_restrict=no
   # Try the official restrict keyword, then gcc's __restrict__, then
   # SGI's __restrict.  __restrict has slightly different semantics than
   # restrict (it's a bit stronger, in that __restrict pointers can't
   # overlap even with non __restrict pointers), but I think it should be
   # okay under the circumstances where restrict is normally used.
   for acx_kw in restrict __restrict__ __restrict; do
     AC_COMPILE_IFELSE([AC_LANG_SOURCE(
      [ float * $acx_kw x;
      ])], [acx_cv_c_restrict=$acx_kw; break])
   done
  ])
 case $acx_cv_c_restrict in
   restrict) ;;
   no) AC_DEFINE(restrict,,
	[Define to equivalent of C99 restrict keyword, or to nothing if this
	is not supported.  Do not define if restrict is supported directly.]) ;;
   *)  AC_DEFINE_UNQUOTED(restrict, $acx_cv_c_restrict) ;;
 esac
])
