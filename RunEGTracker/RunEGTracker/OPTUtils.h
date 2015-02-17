/*--------------------------------------------------------------------------*/
/*--------------------------- File OPTUtils.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small classes are provided for:
 * - reading the time of a code;
 * - generating random numbers.
 *
 * The classes can be adapted to different environments setting a
 * compile-time switch in this file.
 *
 * Additionally, a function is provided for safely reading parameters
 * out of a stream.
 *
 * \version 2.63
 *
 * \date 29 - 08 - 2005
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright 1994 - 2004
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OPTUtils
 #define __OPTUtils   /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*----------------------------- MACROS -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTUTILS_MACROS Compile-time switches in OPTUtils.h
    These macros control how the classes OPTTimers and OPTrand are
    implemented; choose the appropriate value for your environment,
    or program a new version if no value suits you.
    Also, namespaces can be eliminated if they create problems.
    @{ */


/*----------------------- OPT_USE_NAMESPACES -------------------------------*/

#define OPT_USE_NAMESPACES 0

/**< Setting OPT_USE_NAMESPACES == 0 should instruct all codes that use
   OPTUtils stuff to avoid using namespaces; to start with, the common
   namespace OPTUtils_di_unipi_it, that contains all the types defined
   herein, is *not* defined. */

/*---------------------------- OPT_TIMERS ----------------------------------*/

#define OPT_TIMERS 5

/**< The class OPTtimers is defined below to give an abstract interface to the
   different timing routines that are used in different platforms. This is
   needed since time-related functions are one of the less standard parts of
   the C[++] library. The value of the OPT_TIMERS constant selects among the
   different timing routines:

   - 1 = Use the Unix times() routine in sys/times.h

   - 2 = As 1 but uses sys/timeb.h (typical for Microsoft(TM) compilers)

   - 3 = Still use times() of sys/times.h, but return wallclock time
         rather than CPU time

   - 4 = As 3 but uses sys/timeb.h (typical for Microsoft(TM) compilers)

   - 5 = return the user time obtained with ANSI C clock() function; this
         may result in more accurate running times w.r.t. but may be limited
         to ~ 72 hours on systems where ints are 32bits.

   - 6 = Use the Unix gettimeofday() routine of sys/time.h.

   Any unsupported value would simply make the class to report constant
   zero as the time.

   The values 1 .. 4 rely on the constant CLK_TCK for converting between
   clock ticks and seconds; for the case where the constant is not defined by
   the compiler -- should not happen, but it does -- or it is defined in a
   wrong way, the constant is re-defined below. */

/*---------------------------- OPT_RANDOM ---------------------------------*/

#define OPT_RANDOM 1

/**< The class OPTrand is defined below to give an abstract interface to the
   different random generators that are used in different platforms. This is
   needed since random generators are one of the less standard parts of the
   C[++] library. The value of the OPT_RANDOM constant selects among the
   different timing routines:

   - 0 = an hand-made implementation of a rather good random number generator
         is used; note that this assumes that long ints >= 32 bits

   - 1 = standard rand() / srand() pair, common to all C libreries but not
         very sophisticated

   - 2 = drand48() / srand48(), common on Unix architectures and pretty good.

   Any unsupported value would simply make the functions to report constant
   zero, which is not nice but useful to quickly fix problems if you don't
   use random numbers at all. */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_RANDOM )
 #include <stdlib.h>

 /* For random routines, see OPTrand() below. */
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( OPT_TIMERS <= 4 )
 #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 3 ) )
  #include <sys/times.h>
 #else
  #include <sys/timeb.h>
 #endif
 #include <limits.h>
#elif( OPT_TIMERS == 5 )
 #include <time.h>
#elif( OPT_TIMERS == 6 )
 #include <sys/time.h>
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#include <iostream>

/* For istream and the >> operator, used in DfltdSfInpt(). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace OPTUtils_di_unipi_it
{
 /** @namespace OPTUtils_di_unipi_it
     The namespace OPTUtils_di_unipi_it is defined to hold all the data
     types, constants, classes and functions defined here. It also
     comprises the namespace std. */
#endif

 using namespace std;  // I know it's not elegant, but ...

/*--------------------------------------------------------------------------*/
/*----------------------------- NULL ---------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NULL
 #define NULL 0
#endif

/* Curiously enough, some compilers do not always define NULL; for instance,
   sometimes it is defined in stdlib.h. */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLK_TCK -------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef CLK_TCK
 #define CLK_TCK 100
#endif

/* CLK_TCK is the constant used to transform (integer) time values in seconds
   for times()-based routines. Its normal value is 100. */

/*--------------------------------------------------------------------------*/
/*--------------------------- OPT_TIMERS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTUtils_CLASSES Classes in OPTUtils.h
    @{ */

#if( OPT_TIMERS )

/** Provides a common interface to the different timing routines that are
    available in different platforms. */

class OPTtimers {

 public:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// constructor of the class
  OPTtimers( void )
  {
   ReSet();
   }

  /// start the timer
  void Start( void )
  {
   if( ! ticking )
   {
    #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
     times( &buff );
     t_u = buff.tms_utime;
     t_s = buff.tms_stime;
    #elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
     t_u = times( &buff );
    #elif( OPT_TIMERS == 5 )
     t_u = clock();
    #elif( OPT_TIMERS == 6 )
     struct timeval t;
     gettimeofday( &t , NULL );
     t_u = double( t.tv_sec + t.tv_usec * 1e-6 );
    #endif

    ticking = 1;
    }
   }

  /// stop the timer
  void Stop( void )
  {
   if( ticking )
   {
    Read( u , s );
    ticking = 0;
    }
   }

  /** Return the elapsed time. If the clock is ticking, return the *total*
    time since the last Start() without stopping the clock; otherwise,
    return the total elapsed time of all the past runs of the clock since
    the last ReSet() [see below]. */

  double Read( void )
  {
   double tu = 0;
   double ts = 0;
   Read( tu , ts );
   return( tu + ts );
   }

  /// As Read( void ) but *adds* user and system time to tu and ts.
  void Read( double &tu , double &ts )
  {
   if( ticking )
   {
    #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
     times( &buff );
     tu += ( double( buff.tms_utime - t_u ) ) / double( CLK_TCK );
     ts += ( double( buff.tms_stime - t_s ) ) / double( CLK_TCK );
    #elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
     tu += ( double( times( &buff ) - t_u ) ) / double( CLK_TCK );
    #elif( OPT_TIMERS == 5 )
     tu += ( double( clock() - t_u ) ) / double( CLOCKS_PER_SEC );
    #elif( OPT_TIMERS == 6 )
     struct timeval t;
     gettimeofday( &t , NULL );
     tu += double( t.tv_sec + t.tv_usec * 1e-6 ) - t_u;
    #endif
    }
   else { tu += u; ts += s; }
   }

  /// reset the timer
  void ReSet( void )
  {
   u = s = 0; ticking = 0;
   }

 private:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double u;      // elapsed *user* time, in seconds
 double s;      // elapsed *system* time, in seconds
 char ticking;  // 1 if the clock is ticking

 #if( ( OPT_TIMERS > 0 ) && ( OPT_TIMERS <= 5 ) )
  clock_t t_u;

  #if( OPT_TIMERS <= 4 )
   struct tms buff;

   #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 3 ) )
    clock_t t_s;
   #endif
  #endif
 #elif( OPT_TIMERS == 6 )
  double t_u;
 #endif

 };  // end( class OPTtimers );

#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ OPTrand() ---------------------------------*/
/*--------------------------------------------------------------------------*/

/** Provide a common interface to the different random generators that are
    available in different platforms. */

class OPTrand {

 public:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// constructor of the class
  OPTrand( void )
  {
   #if( OPT_RANDOM == 0 )
    A[ 0 ] = -1;
    srand( long( 123456789 ) );
   #else
    OPTrand::srand( long( 1 ) );
   #endif
   }

  /** Returns a random number uniformly distributed in [0, 1).
      \note each object of class OPTrand has its own sequence, so that
      multiple OPTrand objects being used within the same program do not
      interfere with each other (as opposed to what C random routines
      would do). */

  double rand( void )
  {
   #if( OPT_RANDOM == 0 )
    long nmbr = *(gb_fptr--);
    if( nmbr < 0 )
     nmbr = gb_flip_cycle();
    return( double( nmbr ) / double( (unsigned long)0x80000000 ) );
   #elif( OPT_RANDOM == 1 )
    ::srand( myseed );
    myseed = ::rand();
    return( double( myseed ) / double( RAND_MAX ) );
   #elif( OPT_RANDOM == 2 )
    return( erand48( myseed ) );
   #else
    return( 0 );  // random, eh?
   #endif
   }

  /// Seeds the random generator for this instance of OPTrand.
  void srand( long seed )
  {
   #if( OPT_RANDOM == 0 )
    long prev = seed , next = 1;
    seed = prev = mod_diff( prev , 0 );
    A[ 55 ] = prev;

    for( long i = 21 ; i ; i = ( i + 21 ) % 55 ) {
     A[ i ] = next;
     next = mod_diff( prev , next );
     if( seed & 1 )
      seed = 0x40000000 + ( seed >> 1 );
     else
      seed >>= 1;

     next = mod_diff( next , seed );
     prev = A[ i ];
     }

    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
   #elif( OPT_RANDOM == 1 )
    myseed = int( seed );
   #elif( OPT_RANDOM == 2 )
    long *sp = (long*)( &myseed );
    *sp = seed;                // copy higher 32 bits
    myseed[ 2 ] = 0x330E;      // initialize lower 16 bits
   #endif
   }

 private:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( OPT_RANDOM == 0 )
  long A[ 56 ];
  long *gb_fptr;

  long mod_diff( long x , long y )
  {
   return( ( ( x ) - ( y ) ) & 0x7fffffff );
   }

  long gb_flip_cycle( void )
  {
   long *ii, *jj;
   for( ii = &A[ 1 ] , jj = &A[ 32 ] ; jj <= &A[ 55 ] ; ii++ , jj++ )
    *ii = mod_diff( *ii , *jj );

   for( jj = &A[ 1 ] ; ii <= &A[ 55 ] ; ii++ , jj++ )
    *ii = mod_diff( *ii , *jj );

   gb_fptr = &A[ 54 ];

   return A[ 55 ];
   }
 #elif( OPT_RANDOM == 1 )
  int myseed;
 #elif( OPT_RANDOM == 2 )
  unsigned short int myseed[ 3 ];
 #endif

 };  // end( class( OPTrand ) )

/* @} end( group( OPTUtils_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------- DfltdSfInpt() --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTUtils_FUNCTIONS Functions in OPTUtils.h
    @{ */

/** Template function for reading parameters from a istream. The function is
   "safe" because it works also if the istream is not given, is not be long
   enough or contains erroneous things.

   Given a &istream (possibly NULL), DfltdSfInpt() attempts to read Param out
   of it, skipping any line that begins with the comment carachter (defaulted
   to '#'), any blank line and any line starting with anything that can not
   be interpreted as a `T'. If, for any reason, the read operation fails,
   then the parameter is given the default value `Dflt'. Otherwise, all the
   rest of the line up to the nearest newline ('\n') carachter is flushed.

   \note lines should not be longer than 255 carachters. */

template<class T>
inline void DfltdSfInpt( istream *iStrm , T &Param , const T Dflt ,
                         const char cmntc = '#' )
{
 static char buf[ 255 ];

 if( iStrm && ( ! ( ! (*iStrm) ) ) )
  // the "! ! stream" trick is there to force the compiler to apply the
  // stream -> bool conversion, in case it is too dumb to do it by itself
  for( char c ;; (*iStrm).getline( buf , 255 ) )
  {
   if( ! ( (*iStrm) >> ws ) )         // skip whitespace
    break;
   
   if( ! ( (*iStrm).get( c ) ) )      // read first non-whitespace
    break;

   if( c != cmntc )                   // that's not a comment
   {
    (*iStrm).seekg( -1 , ios::cur );  // backtrack
    if( ( (*iStrm) >> Param ) )       // try reading it
     (*iStrm).getline( buf , 255 );   // upon success, skip the rest of line

    return;                           // done
    }
   }

 Param = Dflt; 

 }  // end( DfltdSfInpt )

/* @} end( group( OPTUtils_FUNCTIONS ) ) */
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace OPTUtils_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OPTUtils.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OPTUtils.h -------------------------------*/
/*--------------------------------------------------------------------------*/
