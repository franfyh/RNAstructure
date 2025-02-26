/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 * 
 * This file is not intended to be easily readable and contains a number of 
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG 
 * interface file instead. 
 * ----------------------------------------------------------------------------- */

#define SWIGJAVA


#ifdef __cplusplus
/* SwigValueWrapper is described in swig.swg */
template<typename T> class SwigValueWrapper {
  struct SwigMovePointer {
    T *ptr;
    SwigMovePointer(T *p) : ptr(p) { }
    ~SwigMovePointer() { delete ptr; }
    SwigMovePointer& operator=(SwigMovePointer& rhs) { T* oldptr = ptr; ptr = 0; delete oldptr; ptr = rhs.ptr; rhs.ptr = 0; return *this; }
  } pointer;
  SwigValueWrapper& operator=(const SwigValueWrapper<T>& rhs);
  SwigValueWrapper(const SwigValueWrapper<T>& rhs);
public:
  SwigValueWrapper() : pointer(0) { }
  SwigValueWrapper& operator=(const T& t) { SwigMovePointer tmp(new T(t)); pointer = tmp; return *this; }
  operator T&() const { return *pointer.ptr; }
  T *operator&() { return pointer.ptr; }
};

template <typename T> T SwigValueInit() {
  return T();
}
#endif

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__)) 
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__)) 
# else
#   define SWIGUNUSED 
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif 
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED 
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* exporting methods */
#if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#  ifndef GCC_HASCLASSVISIBILITY
#    define GCC_HASCLASSVISIBILITY
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif 
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif



/* Fix for jlong on some versions of gcc on Windows */
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  typedef long long __int64;
#endif

/* Fix for jlong on 64-bit x86 Solaris */
#if defined(__x86_64)
# ifdef _LP64
#   undef _LP64
# endif
#endif

#include <jni.h>
#include <stdlib.h>
#include <string.h>


/* Support for throwing Java exceptions */
typedef enum {
  SWIG_JavaOutOfMemoryError = 1, 
  SWIG_JavaIOException, 
  SWIG_JavaRuntimeException, 
  SWIG_JavaIndexOutOfBoundsException,
  SWIG_JavaArithmeticException,
  SWIG_JavaIllegalArgumentException,
  SWIG_JavaNullPointerException,
  SWIG_JavaDirectorPureVirtual,
  SWIG_JavaUnknownError
} SWIG_JavaExceptionCodes;

typedef struct {
  SWIG_JavaExceptionCodes code;
  const char *java_exception;
} SWIG_JavaExceptions_t;


static void SWIGUNUSED SWIG_JavaThrowException(JNIEnv *jenv, SWIG_JavaExceptionCodes code, const char *msg) {
  jclass excep;
  static const SWIG_JavaExceptions_t java_exceptions[] = {
    { SWIG_JavaOutOfMemoryError, "java/lang/OutOfMemoryError" },
    { SWIG_JavaIOException, "java/io/IOException" },
    { SWIG_JavaRuntimeException, "java/lang/RuntimeException" },
    { SWIG_JavaIndexOutOfBoundsException, "java/lang/IndexOutOfBoundsException" },
    { SWIG_JavaArithmeticException, "java/lang/ArithmeticException" },
    { SWIG_JavaIllegalArgumentException, "java/lang/IllegalArgumentException" },
    { SWIG_JavaNullPointerException, "java/lang/NullPointerException" },
    { SWIG_JavaDirectorPureVirtual, "java/lang/RuntimeException" },
    { SWIG_JavaUnknownError,  "java/lang/UnknownError" },
    { (SWIG_JavaExceptionCodes)0,  "java/lang/UnknownError" }
  };
  const SWIG_JavaExceptions_t *except_ptr = java_exceptions;

  while (except_ptr->code != code && except_ptr->code)
    except_ptr++;

  jenv->ExceptionClear();
  excep = jenv->FindClass(except_ptr->java_exception);
  if (excep)
    jenv->ThrowNew(excep, msg);
}


/* Contract support */

#define SWIG_contract_assert(nullreturn, expr, msg) if (!(expr)) {SWIG_JavaThrowException(jenv, SWIG_JavaIllegalArgumentException, msg); return nullreturn; } else


#include "DotPlotBackend.h"


#include <string>


#ifdef __cplusplus
extern "C" {
#endif

SWIGEXPORT jlong JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_new_1DotPlotBackend(JNIEnv *jenv, jclass jcls) {
  jlong jresult = 0 ;
  DotPlotBackend *result = 0 ;
  
  (void)jenv;
  (void)jcls;
  result = (DotPlotBackend *)new DotPlotBackend();
  *(DotPlotBackend **)&jresult = result; 
  return jresult;
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_delete_1DotPlotBackend(JNIEnv *jenv, jclass jcls, jlong jarg1) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  
  (void)jenv;
  (void)jcls;
  arg1 = *(DotPlotBackend **)&jarg1; 
  delete arg1;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getBounds(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (arg1)->getBounds();
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getColors(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (arg1)->getColors();
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jdouble JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getCurrentMax(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jdouble jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  double result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (double)(arg1)->getCurrentMax();
  jresult = (jdouble)result; 
  return jresult;
}


SWIGEXPORT jdouble JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getDefaultMin(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jdouble jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  double result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (double)(arg1)->getDefaultMin();
  jresult = (jdouble)result; 
  return jresult;
}


SWIGEXPORT jdouble JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getDefaultMax(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jdouble jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  double result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (double)(arg1)->getDefaultMax();
  jresult = (jdouble)result; 
  return jresult;
}


SWIGEXPORT jdouble JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getCurrentMin(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jdouble jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  double result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (double)(arg1)->getCurrentMin();
  jresult = (jdouble)result; 
  return jresult;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getDotData(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jint jarg2, jint jarg3) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  int arg2 ;
  int arg3 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  arg2 = (int)jarg2; 
  arg3 = (int)jarg3; 
  result = (arg1)->getDotData(arg2,arg3);
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getDotDataByLocation(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jint jarg2, jint jarg3, jdouble jarg4) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  int arg2 ;
  int arg3 ;
  double arg4 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  arg2 = (int)jarg2; 
  arg3 = (int)jarg3; 
  arg4 = (double)jarg4; 
  result = (arg1)->getDotDataByLocation(arg2,arg3,arg4);
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getGridLine(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jint jarg2) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  int arg2 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  arg2 = (int)jarg2; 
  result = (arg1)->getGridLine(arg2);
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getLegendData(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (arg1)->getLegendData();
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT jint JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1getLength(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_) {
  jint jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  int result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  result = (int)(arg1)->getLength();
  jresult = (jint)result; 
  return jresult;
}


SWIGEXPORT jboolean JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1readData(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jstring jarg2, jint jarg3) {
  jboolean jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string arg2 ;
  int arg3 ;
  bool result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  if(!jarg2) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return 0;
  } 
  const char *arg2_pstr = (const char *)jenv->GetStringUTFChars(jarg2, 0); 
  if (!arg2_pstr) return 0;
  (&arg2)->assign(arg2_pstr);
  jenv->ReleaseStringUTFChars(jarg2, arg2_pstr); 
  arg3 = (int)jarg3; 
  result = (bool)(arg1)->readData(arg2,arg3);
  jresult = (jboolean)result; 
  return jresult;
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1setColors(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jint jarg2) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  int arg2 ;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  arg2 = (int)jarg2; 
  (arg1)->setColors(arg2);
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1setRange(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jdouble jarg2, jdouble jarg3) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  double arg2 ;
  double arg3 ;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  arg2 = (double)jarg2; 
  arg3 = (double)jarg3; 
  (arg1)->setRange(arg2,arg3);
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1writePostscriptFile(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jstring jarg2) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string arg2 ;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  if(!jarg2) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return ;
  } 
  const char *arg2_pstr = (const char *)jenv->GetStringUTFChars(jarg2, 0); 
  if (!arg2_pstr) return ;
  (&arg2)->assign(arg2_pstr);
  jenv->ReleaseStringUTFChars(jarg2, arg2_pstr); 
  (arg1)->writePostscriptFile(arg2);
}


SWIGEXPORT jstring JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1writeStructuresFile(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jstring jarg2, jstring jarg3) {
  jstring jresult = 0 ;
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string arg2 ;
  std::string arg3 ;
  std::string result;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  if(!jarg2) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return 0;
  } 
  const char *arg2_pstr = (const char *)jenv->GetStringUTFChars(jarg2, 0); 
  if (!arg2_pstr) return 0;
  (&arg2)->assign(arg2_pstr);
  jenv->ReleaseStringUTFChars(jarg2, arg2_pstr); 
  if(!jarg3) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return 0;
  } 
  const char *arg3_pstr = (const char *)jenv->GetStringUTFChars(jarg3, 0); 
  if (!arg3_pstr) return 0;
  (&arg3)->assign(arg3_pstr);
  jenv->ReleaseStringUTFChars(jarg3, arg3_pstr); 
  result = (arg1)->writeStructuresFile(arg2,arg3);
  jresult = jenv->NewStringUTF((&result)->c_str()); 
  return jresult;
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1writeSVGFile(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jstring jarg2) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string arg2 ;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  if(!jarg2) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return ;
  } 
  const char *arg2_pstr = (const char *)jenv->GetStringUTFChars(jarg2, 0); 
  if (!arg2_pstr) return ;
  (&arg2)->assign(arg2_pstr);
  jenv->ReleaseStringUTFChars(jarg2, arg2_pstr); 
  (arg1)->writeSVGFile(arg2);
}


SWIGEXPORT void JNICALL Java_RNAstructure_1java_1drawing_source_proxy_DotPlotBackendProxyJNI_DotPlotBackend_1writeTextFile(JNIEnv *jenv, jclass jcls, jlong jarg1, jobject jarg1_, jstring jarg2) {
  DotPlotBackend *arg1 = (DotPlotBackend *) 0 ;
  std::string arg2 ;
  
  (void)jenv;
  (void)jcls;
  (void)jarg1_;
  arg1 = *(DotPlotBackend **)&jarg1; 
  if(!jarg2) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "null std::string");
    return ;
  } 
  const char *arg2_pstr = (const char *)jenv->GetStringUTFChars(jarg2, 0); 
  if (!arg2_pstr) return ;
  (&arg2)->assign(arg2_pstr);
  jenv->ReleaseStringUTFChars(jarg2, arg2_pstr); 
  (arg1)->writeTextFile(arg2);
}


#ifdef __cplusplus
}
#endif

