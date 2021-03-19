#ifndef T_MATH
#define T_MATH 1

#if defined(T_WITH_MATH_EXPORTC)
extern "C" {
#include <math.h>
}
#else
#include <math.h>
#endif


#undef abs
#undef acos
#undef ceil
#undef cos
#undef exp
#undef fabs
#undef floor
#undef log
#undef log10
#undef sin
#undef sqrt
#undef tan

namespace std 
{

  using ::acos;

#if HAVE_ACOSF
  inline float 
  acos(float __x) { return ::acosf(__x); }
#else
  inline float 
  acos(float __x) { return ::acos(static_cast<double>(__x)); }
#endif

  using ::ceil;
  
#if HAVE_CEILF
  inline float 
  ceil(float __x) { return ::ceilf(__x); }
#else
  inline float 
  ceil(float __x) { return ::ceil(static_cast<double>(__x)); }
#endif

  using ::cos;

#if HAVE_COSF
  inline float
  cos(float __x) { return ::cosf(__x); }
#else
  inline float
  cos(float __x) { return ::cos(static_cast<double>(__x)); }
#endif

  using ::exp;

#if HAVE_EXPF
  inline float 
  exp(float __x) { return ::expf(__x); }
#else
  inline float 
  exp(float __x) { return ::exp(static_cast<double>(__x)); }
#endif

  using ::fabs;

  inline double
  abs(double __x) { return ::fabs(__x); }

#if HAVE_FABSF
  inline float
  fabs(float __x) { return ::fabsf(__x); }
  inline float
  abs(float __x) { return ::fabsf(__x); }
#else
  inline float
  fabs(float __x) { return ::fabs(static_cast<double>(__x)); }
  inline float
  abs(float __x) { return ::fabs(static_cast<double>(__x)); }
#endif

  using ::floor;

#if HAVE_FLOORF
  inline float 
  floor(float __x) { return ::floorf(__x); }
#else
  inline float 
  floor(float __x) { return ::floor(static_cast<double>(__x)); }
#endif

  using ::log;

#if HAVE_LOGF
  inline float 
  log(float __x) { return ::logf(__x); }
#else
  inline float log(float __x) { return ::log(static_cast<double>(__x)); }
#endif

  using ::log10;

#if HAVE_LOG10F
  inline float 
  log10(float __x) { return ::log10f(__x); }
#else
  inline float 
  log10(float __x) { return ::log10(static_cast<double>(__x)); }
#endif

  using ::sin;

#if HAVE_SINF
  inline float
  sin(float __x) { return ::sinf(__x); }
#else
  inline float
  sin(float __x) { return ::sin(static_cast<double>(__x)); }
#endif

  using ::sqrt;

#if HAVE_SQRTF
  inline float
  sqrt(float __x) { return ::sqrtf(__x); }
#else
  inline float
  sqrt(float __x) { return ::sqrt(static_cast<double>(__x)); }
#endif

  using ::tan;

#if HAVE_TANF
  inline float 
  tan(float __x) { return ::tanf(__x); }
#else
  inline float 
  tan(float __x) { return ::tan(static_cast<double>(__x)); }
#endif

} // namespace end
  
#endif
