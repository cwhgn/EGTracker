/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFClass.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the abstract base class MCFClass, which defines a standard
 * interface for (linear or convex quadratic separable) Min Cost Flow Problem
 * solvers, to be implemented as derived classes.
 *
 * \version 3.01
 *
 * \date 30 - 09 - 2011
 *
 * \author Alessandro Bertolini \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * Copyright &copy 1996 - 2011 by Antonio Frangioni, Claudio Gentile
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFClass
 #define __MCFClass  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFCLASS_MACROS Compile-time switches in MCFClass.h
    These macros control some important details of the class interface.
    Although using macros for activating features of the interface is not
    very C++, switching off some unused features may allow some
    implementation to be more efficient in running time or memory.
    @{ */

/*-------------------------------- USENAME0 --------------------------------*/

#define USENAME0 0

/**< Decides if 0 or 1 is the "name" of the first node.
   If USENAME0 == 1, (warning: it has to be *exactly* 1), then the node
   names go from 0 to n - 1, otherwise from 1 to n. Note that this does not
   affect the position of the deficit in the deficit vectors, i.e., the
   deficit of the i-th node - be its "name" `i' or `i - 1' - is always in
   the i-th position of the vector. */

/*@}  end( group( MCFCLASS_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTUtils.h"

/* OPTtypes.h defines standard interfaces for timing and random routines, as
   well as the namespace OPTtypes_di_unipi_it and the macro
   OPT_USE_NAMESPACES, useful for switching off all namespaces in one blow
   for those strange cases where they create problems. */

#include <iomanip>
#include <sstream>
#include <limits>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
 /** @namespace MCFClass_di_unipi_it
     The namespace MCFClass_di_unipi_it is defined to hold the MCFClass
     class and all the relative stuff. It comprises the namespace
     OPTtypes_di_unipi_it. */

 using namespace OPTtypes_di_unipi_it;
#endif

/*@}  end( group( MCFCLASS_CONSTANTS ) ) */
/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFClass --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFCLASS_CLASSES Classes in MCFClass.h
    @{ */

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** This abstract base class defines a standard interface for (linear or
    convex quadartic separable) Min Cost Flow (MCF) problem solvers.

    The data of the problem consist of a (directed) graph G = ( N , A ) with
    n = |N| nodes and m = |A| (directed) arcs. Each node `i' has a deficit
    b[ i ], i.e., the amount of flow that is produced/consumed by the node:
    source nodes (which produce flow) have negative deficits and sink nodes
    (which consume flow) have positive deficits. Each arc `(i, j)' has an
    upper capacity U[ i , j ], a linear cost coefficient C[ i , j ] and a
    (non negative) quadratic cost coefficient Q[ i , j ]. Flow variables
    X[ i , j ] represents the amount of flow  to be sent on arc (i, j).
    Parallel arcs, i.e., multiple copies of the same arc `(i, j)' (with
    possibily different costs and/or capacities) are in general allowed.
    The formulation of the problem is therefore:
    \f[
     \min \sum_{ (i, j) \in A } C[ i , j ] X[ i, j ] +
                                Q[ i , j ] X[ i, j ]^2 / 2
    \f]
    \f[
     (1) \sum_{ (j, i) \in A } X[ j , i ] -
         \sum_{ (i, j) \in A } X[ i , j ] = b[ i ]
         \hspace{1cm} i \in N
    \f]
    \f[
     (2) 0 \leq X[ i , j ] \leq U[ i , j ]
         \hspace{1cm} (i, j) \in A
    \f]
    The n equations (1) are the flow conservation constraints and the 2m
    inequalities (2) are the flow nonnegativity and capacity constraints.
    At least one of the flow conservation constraints is redundant, as the
    demands must be balanced (\f$\sum_{ i \in N } b[ i ] = 0\f$); indeed,
    exactly n - ConnectedComponents( G ) flow conservation constraints are
    redundant, as demands must be balanced in each connected component of G.
    Let us denote by QA and LA the disjoint subsets of A containing,
    respectively, "quadratic" arcs (with Q[ i , j ] > 0) and "linear" arcs
    (with Q[ i , j ] = 0); the (MCF) problem is linear if QA is empty, and
    nonlinear (convex quadratic) if QA is nonempty.

    The dual of the problem is:
    \f[
     \max \sum_{ i \in N } Pi[ i ] b[ i ] -
          \sum_{ (i, j) \in A } W[ i , j ] U[ i , j ] -
          \sum_{ (i, j) \in AQ } V[ i , j ]^2 / ( 2 * Q[ i , j ] )
    \f]
    \f[
     (3.a) C[ i , j ] - Pi[ j ] + Pi[ i ] + W[ i , j ] - Z[ i , j ] = 0
           \hspace{1cm} (i, j) \in AL
    \f]
    \f[
     (3.b) C[ i , j ] - Pi[ j ] + Pi[ i ] + W[ i , j ] - Z[ i , j ] =
           V[ i , j ]
           \hspace{1cm} (i, j) \in AQ
    \f]
    \f[
     (4.a) W[ i , j ] \geq 0 \hspace{1cm} (i, j) \in A
    \f]
    \f[
     (4.b) Z[ i , j ] \geq 0 \hspace{1cm} (i, j) \in A
    \f]

    Pi[] is said the vector of node potentials for the problem, W[] are
    bound variables and Z[] are slack variables. Given Pi[], the quantities
    \f[
     RC[ i , j ] =  C[ i , j ] + Q[ i , j ] * X[ i , j ] - Pi[ j ] + Pi[ i ]
    \f]
    are said the "reduced costs" of arcs.

    A primal and dual feasible solution pair is optimal if and only if the
    complementary slackness conditions
    \f[
     RC[ i , j ] > 0 \Rightarrow X[ i , j ] = 0
    \f]
    \f[
     RC[ i , j ] < 0 \Rightarrow X[ i , j ] = U[ i , j ]
    \f]
    are satisfied for all arcs (i, j) of A.

    The MCFClass class provides an interface with methods for managing and
    solving problems of this kind. Actually, the class can also be used as
    an interface for more general NonLinear MCF problems, where the cost
    function either nonseparable ( C( X ) ) or arc-separable
    ( \f$\sum_{ (i, j) \in A } C_{i,j}( X[ i, j ] )\f$ ). However, solvers
    for NonLinear MCF problems are typically objective-function-specific,
    and there is no standard way for inputting a nonlinear function different
    from a separable convex quadratic one, so only the simplest form is dealt
    with in the interface, leaving more complex NonLinear parts to the
    interface of derived classes. */

class MCFClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    The MCFClass defines four main public types:

    - Index, the type of arc and node indices;

    - FNumber, the type of flow variables, arc capacities, and node deficits;

    - CNumber, the type of flow costs, node potentials, and arc reduced costs;

    - FONumber, the type of objective function value.

    By re-defining the types in this section, most MCFSolver should be made
    to work with any reasonable choice of data type (= one that is capable of
    properly representing the data of the instances to be solved). This may
    be relevant due to an important property of MCF problems: *if all arc
    capacities and node deficits are integer, then there exists an integral
    optimal primal solution*, and  *if all arc costs are integer, then there
    exists an integral optimal dual solution*. Even more importantly, *many
    solution algorithms will in fact produce an integral primal/dual
    solution for free*, because *every primal/dual solution they generate
    during the solution process is naturally integral*. Therefore, one can
    use integer data types to represent everything connected with flows and/or
    costs if the corresponding data is integer in all instances one needs to
    solve. This directly translates in significant memory savings and/or speed
    improvements.

    *It is the user's responsibility to ensure that these types are set to
    reasonable values*. So, the experienced user may want to experiment with
    setting this data properly if memory footprint and/or speed is a primary
    concern. Note, however, that *not all solution algorithms will happily
    accept integer data*; one example are Interior-Point approaches, which
    require both flow and cost variables to be continuous (float). So, the
    viability of setting integer data (as well as its impact on performances)
    is strictly related to the specific kind of algorithm used. Since these
    types are common to all derived classes, they have to be set taking into
    account the needs of all the solvers that are going to be used, and
    adapting to the "worst case"; of course, FNumber == CNumber == double is
    going to always be an acceptable "worst case" setting. MCFClass may in a
    future be defined as a template class, with these as template parameters,
    but this is currently deemed overkill and avoided.

    Finally, note that the above integrality property only holds for *linear*
    MCF problems. If any arc has a nonzero quadratic cost coefficient, optimal
    flows and potentials may be fractional even if all the data of the problem
    (comprised quadratic cost coefficients) is integer. Hence, for *quadratic*
    MCF solvers, a setting like FNumber == CNumber == double is actually
    *mandatory*, for any reasonable algorithm will typically misbehave
    otherwise.
    @{ */

/*--------------------------------------------------------------------------*/

 typedef unsigned int    Index;           ///< index of a node or arc ( >= 0 )
 typedef Index          *Index_Set;       ///< set (array) of indices
 typedef const Index    cIndex;           ///< a read-only index
 typedef cIndex        *cIndex_Set;       ///< read-only index array

/*--------------------------------------------------------------------------*/

 typedef int             SIndex;           ///< index of a node or arc 
 typedef SIndex         *SIndex_Set;       ///< set (array) of indices
 typedef const SIndex   cSIndex;           ///< a read-only index
 typedef cSIndex       *cSIndex_Set;       ///< read-only index array

/*--------------------------------------------------------------------------*/

 typedef double          FNumber;        ///< type of arc flow
 typedef FNumber        *FRow;           ///< vector of flows
 typedef const FNumber  cFNumber;        ///< a read-only flow
 typedef cFNumber      *cFRow;           ///< read-only flow array

/*--------------------------------------------------------------------------*/

 typedef double          CNumber;        ///< type of arc flow cost
 typedef CNumber        *CRow;           ///< vector of costs
 typedef const CNumber  cCNumber;        ///< a read-only cost
 typedef cCNumber      *cCRow;           ///< read-only cost array

/*--------------------------------------------------------------------------*/

 typedef double          FONumber; 
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

 typedef const FONumber cFONumber;       ///< a read-only o.f. value

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type (FNumber, CNumber, Index); just use Inf<type>(). */

   template <typename T>
   class Inf {
    public:
     Inf() {}
     operator T() { return( std::numeric_limits<T>::max() ); }
    };

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "machine epsilon" for a
    basic type (FNumber, CNumber); just use Eps<type>(). */

   template <typename T>
   class Eps {
    public:
     Eps() {}
     operator T() { return( std::numeric_limits<T>::epsilon() ); }
    };

/*--------------------------------------------------------------------------*/
/** Small class for exceptions. Derives from std::exception implementing the
    virtual method what() -- and since what is virtual, remember to always
    catch it by reference (catch exception &e) if you want the thing to work.
    MCFException class are thought to be of the "fatal" type, i.e., problems
    for which no solutions exists apart from aborting the program. Other kinds
    of exceptions can be properly handled by defining derived classes with
    more information. */

 class MCFException : public exception {
 public:
  MCFException( const char *const msg = 0 ) { errmsg = msg; }

  const char* what( void ) const throw () { return( errmsg ); }
 private:
  const char *errmsg;
  };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible parameters of the MCF solver, to be
    used with the methods SetPar() and GetPar(). */

  enum MCFParam { kMaxTime = 0 ,     ///< max time 
                  kMaxIter ,         ///< max number of iteration
                  kEpsFlw ,          ///< tolerance for flows
                  kEpsDfct ,         ///< tolerance for deficits
                  kEpsCst ,          ///< tolerance for costs
                  kReopt ,           ///< whether or not to reoptimize
                  kLastParam         /**< dummy parameter: this is used to
                                        allow derived classes to "extend"
                                        the set of parameters. */
                  };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible status of the MCF solver. */

  enum MCFStatus { kUnSolved = -1 , ///< no solution available
                   kOK = 0 ,        ///< optimal solution found

                   kStopped ,       ///< optimization stopped
                   kUnfeasible ,    ///< problem is unfeasible
                   kUnbounded ,     ///< problem is unbounded
                   kError           ///< error in the solver
                   };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible reoptimization status of the MCF
    solver. */

  enum MCFAnswer { kNo = 0 ,  ///< no 
                   kYes       ///< yes
                   };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible file formats in WriteMCF(). */

  enum MCFFlFrmt { kDimacs = 0 ,    ///< DIMACS file format for MCF
                   kQDimacs ,       ///< quadratic DIMACS file format for MCF
                   kMPS ,           ///< MPS file format for LP
                   kFWMPS           ///< "Fixed Width" MPS format
                   };

/*--------------------------------------------------------------------------*/
/** Base class for representing the internal state of the MCF algorithm. */

 class MCFState {
 public:
   MCFState( void ) {};
   virtual ~MCFState() {};
 };

 typedef MCFState *MCFStatePtr;  ///< pointer to a MCFState

/*@} -----------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructors
    @{ */

   MCFClass( cIndex nmx = 0 , cIndex mmx = 0 )
   {
    nmax = nmx;
    mmax = mmx;
    n = m = 0;

    status = kUnSolved;
    Senstv = true;

    EpsFlw = Eps<FNumber>() * 100;
    EpsCst = Eps<CNumber>() * 100;
    EpsDfct = EpsFlw * ( nmax ? nmax : 100 );

    MaxTime = 0;
    MaxIter = 0;

    MCFt = NULL;
    }

/**< Constructor of the class.

   nmx and mmx, if provided, are taken to be respectively the maximum number 
   of nodes and arcs in the network. If nonzero values are passed, memory
   allocation can be anticipated in the constructor, which is sometimes
   desirable. The maximum values are stored in the protected fields nmax and
   mmax, and can be changed with LoadNet() [see below]; however, changing
   them typically requires memory allocation/deallocation, which is
   sometimes undesirable outside the constructor.

   After that an object has been constructed, no problem is loaded; this has
   to be done with LoadNet() [see below]. Thus, it is an error to invoke any
   method which requires the presence of a problem (typicall all except those
   in the initializations part). The base class provides two protected fields
   n and m for the current number of nodes and arcs, respectively, that are
   set to 0 in the constructor precisely to indicate that no instance is
   currently loaded. */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   virtual void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
                         cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
                         cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
                         cIndex_Set pEn = NULL ) = 0;

/**< Inputs a new network.
 
   The parameters nmx and mmx are the new max number of nodes and arcs,
   possibly overriding those set in the constructor [see above], altough at
   the likely cost of memory allocation and deallocation. Passing nmx == 
   mmx == 0 is intended as a signal to the solver to deallocate everything
   and wait for new orders; in this case, all the other parameters are ignored.

   Otherwise, in principle all the other parameters have to be provided.
   Actually, some of them may not be needed for special classes of MCF
   problems (e.g., costs in a MaxFlow problem, or start/end nodes in a
   problem defined over a graph with fixed topology, such as a complete
   graph). Also, passing NULL is allowed to set default values.

   The meaning of the parameters is the following:

   - pn     is the current number of nodes of the network (<= nmax).
   - pm     is the number of arcs of the network (<= mmax).

   - pU     is the m-vector of the arc upper capacities; capacities must be
            nonnegative, but can in principle be infinite (== F_INF); passing
            pU == NULL means that all capacities are infinite;

   - pC     is the m-vector of the arc costs; costs must be finite (< C_INF);
            passing pC == NULL means that all costs must be 0.

   - pDfct  is the n-vector of the node deficits; source nodes have negative
            deficits and sink nodes have positive deficits; passing pDfct ==
            NULL means that all deficits must be 0 (a circulation problem);

   - pSn    is the m-vector of the arc starting nodes; pSn == NULL is in
            principle not allowed, unless the topology of the graph is fixed;
   - pEn    is the m-vector of the arc ending nodes; same comments as for pSn.

   Note that node "names" in the arrays pSn and pEn must go from 1 to pn if
   the macro USANAME0 [see above] is set to 0, while they must go from 0 to
   pn - 1 if USANAME0 is set to 1. In both cases, however, the deficit of the
   first node is read from the first (0-th) position of pDfct, that is if
   USANAME0 == 0 then the deficit of the node with name `i' is read from
   pDfct[ i - 1 ].

   The data passed to LoadNet() can be used to specify that the arc `i' must
   not "exist" in the problem. This is done by passing pC[ i ] == C_INF;
   solvers which don't read costs are forced to read them in order to check
   this, unless they provide alternative solver-specific ways to accomplish
   the same tasks. These arcs are "closed", as for the effect of CloseArc()
   [see below]. "invalid" costs (== C_INF) are set to 0 in order to being
   subsequently capable of "opening" them back with OpenArc()  [see below].
   The way in which these non-existent arcs are phisically dealt with is
   solver-specific; in some solvers, for instance, this could be obtained by
   simply putting their capacity to zero. Details about these issues should
   be found in the interface of derived classes.

   Note that the quadratic part of the objective function, if any, is not
   dealt with in LoadNet(); it can only be separately provided with
   ChgQCoef() [see below]. By default, the problem is linear, i.e., all
   coefficients of the second-order terms in the objective function are
   assumed to be zero. */

/*--------------------------------------------------------------------------*/

   virtual inline void LoadDMX( istream &DMXs , bool IsQuad = false );

/**< Read a MCF instance in DIMACS standard format from the istream. The
   format is the following. The first line must be

      p min <number of nodes> <number of arcs>

   Then the node definition lines must be found, in the form

      n <node number> <node supply>

   Not all nodes need have a node definition line; these are given zero
   supply, i.e., they are transhipment nodes (supplies are the inverse of
   deficits, i.e., a node with positive supply is a source node). Finally,
   the arc definition lines must be found, in the form

      a <start node> <end node> <lower bound> <upper bound> <flow cost>

   There must be exactly <number of arcs> arc definition lines in the file.

   This method is *not* pure virtual because an implementation is provided by
   the base class, using the LoadNet() method (which *is* pure virtual).
   However, the method *is* virtual to allow derived classes to implement
   more efficient versions, should they have any reason to do so.

   \note Actually, the file format accepted by LoadDMX (at least in the
         base class implementation) is more general than the DIMACS standard
         format, in that it is allowed to mix node and arc definitions in
         any order, while the DIMACS file requires all node information to
         appear before all arc information.

   \note Other than for the above, this method is assumed to allow for
         *quadratic* Dimacs files, encoding for convex quadratic separable
         Min Cost Flow instances. This is a simple extension where each arc
         descriptor has a sixth field, <quadratic cost>. The provided
         istream is assumed to be quadratic Dimacs file if IsQuad is true,
         and a regular linear Dimacs file otherwise. */

/*--------------------------------------------------------------------------*/

   virtual void PreProcess( void ) {}

/**< Extract a smaller/easier equivalent MCF problem. The data of the instance
   is changed and the easier one is solved instead of the original one. In the
   MCF case, preprocessing may involve reducing bounds, identifying
   disconnected components of the graph etc. However, proprocessing is
   solver-specific.

   This method can be implemented by derived classes in their solver-specific
   way. Preprocessing may reveal unboundedness or unfeasibility of the
   problem; if that happens, PreProcess() should properly set the `status'
   field, that can then be read with MCFGetStatus() [see below].
 
   Note that preprocessing may destroy all the solution information. Also, it
   may be allowed to change the data of the problem, such as costs/capacities
   of the arcs.

   A valid preprocessing is doing nothing, and that's what the default
   implementation of this method (that is *not* pure virtual) does. */

/*--------------------------------------------------------------------------*/

   virtual inline void SetPar( int par , int val );

/**< Set integer parameters of the algorithm.

   @param par   is the parameter to be set; the enum MCFParam can be used, but
                'par' is an int (every enum is an int) so that the method can
                be extended by derived classes for the setting of their
                parameters

   @param value is the value to assign to the parameter.  

   The base class implementation handles these parameters: 

   - kMaxIter: the max number of iterations in which the MCF Solver can find
               an optimal solution (default 0, which means no limit)

   - kReopt:   tells the solver if it has to reoptimize. The implementation in
               the base class sets a flag, the protected \c bool field \c
               Senstv; if true (default) this field instructs the MCF solver to
               to try to exploit the information about the latest optimal
	       solution to speedup the optimization of the current problem,
	       while if the field is false the MCF solver should restart the
	       optimization "from scratch" discarding any previous information.
	       Usually reoptimization speeds up the computation considerably,
	       but this is not always true, especially if the data of the
	       problem changes a lot.*/

/*--------------------------------------------------------------------------*/

   virtual inline void SetPar( int par , double val );

/**< Set float parameters of the algorithm.

   @param par   is the parameter to be set; the enum MCFParam can be used, but
                'par' is an int (every enum is an int) so that the method can
                be extended by derived classes for the setting of their
                parameters

   @param value is the value to assign to the parameter.  

   The base class implementation handles these parameters: 

   - kEpsFlw:  sets the tolerance for controlling if the flow on an arc is zero 
               to val. This also sets the tolerance for controlling if a node
	       deficit is zero (see kEpsDfct) to val * < max number of nodes >;
	       this value should be safe for graphs in which any node has less
	       than < max number of nodes > adjacent nodes, i.e., for all graphs
	       but for very dense ones with "parallel arcs"

   - kEpsDfct: sets the tolerance for controlling if a node deficit is zero to 
               val, in case a better value than that autmatically set by
               kEpsFlw (see above) is available (e.g., val * k would be good
	       if no node has more than k neighbours)

   - kEpsCst:  sets the tolerance for controlling if the reduced cost of an arc
               is zero to val. A feasible solution satisfying eps-complementary
               slackness, i.e., such that
               \f[
                RC[ i , j ] < - eps \Rightarrow X[ i , j ] = U[ ij ]
               \f]
               and
               \f[
                RC[ i , j ] > eps \Rightarrow X[ i , j ] == 0 ,
               \f]
               is known to be ( eps * n )-optimal.

   - kMaxTime: sets the max time (in seconds) in which the MCF Solver can find
               an optimal solution  (default 0, which means no limit). */

/*--------------------------------------------------------------------------*/

   virtual inline void GetPar( int par , int &val );

/**< This method returns one of the integer parameter of the algorithm.

   @param par  is the parameter to return [see SetPar( int ) for comments];

   @param val  upon return, it will contain the value of the parameter.

   The base class implementation handles the parameters kMaxIter and kReopt.
   */

/*--------------------------------------------------------------------------*/

   virtual inline void GetPar( int par , double &val );

/**< This method returns one of the integer parameter of the algorithm.

   @param par  is the parameter to return [see SetPar( double ) for comments];

   @param val  upon return, it will contain the value of the parameter.

   The base class implementation handles the parameters kEpsFlw, kEpsDfct,
   kEpsCst, and kMaxTime. */

/*--------------------------------------------------------------------------*/

   virtual void SetMCFTime( bool TimeIt = true )
   {
    if( TimeIt )
     if( MCFt )
      MCFt->ReSet();
     else
      MCFt = new OPTtimers();
    else {
     delete MCFt;
     MCFt = NULL;
     }
    }

/**< Allocate an OPTtimers object [see OPTtypes.h] to be used for timing the
   methods of the class. The time can be read with TimeMCF() [see below]. By
   default, or if SetMCFTime( false ) is called, no timing is done. Note that,
   since all the relevant methods ot the class are pure virtual, MCFClass can
   only manage the OPTtimers object, but it is due to derived classes to
   actually implement the timing.

   @note time accumulates over the calls: calling SetMCFTime(), however,
         resets the counters, allowing to time specific groups of calls.

   @note of course, setting kMaxTime [see SetPar() above] to any nonzero
         value has no effect unless SetMCFTime( true ) has been called. */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   virtual void SolveMCF( void ) = 0;

/**< Solver of the Min Cost Flow Problem. Attempts to solve the MCF instance
   currently loaded in the object. */

/*--------------------------------------------------------------------------*/

   inline int MCFGetStatus( void )
   {
    return( status );
    }

/**< Returns an int describing the current status of the MCF solver. Possible
   return values are:

   - kUnSolved    SolveMCF() has not been called yet, or the data of the
                  problem has been changed since the last call;

   - kOK          optimization has been carried out succesfully;

   - kStopped     optimization have been stopped before that the stopping
                  conditions of the solver applied, e.g. because of the
                  maximum allowed number of "iterations" [see SetPar( int )]
                  or the maximum allowed time [see SetPar( double )] has been
		  reached; this is not necessarily an error, as it might just
                  be required to re-call SolveMCF() giving it more "resources"
                  in order to solve the problem;

   - kUnfeasible  if the current MCF instance is (primal) unfeasible;

   - kUnbounded   if the current MCF instance is (primal) unbounded (this can
                  only happen if the solver actually allows F_INF capacities,
                  which is nonstandard in the interface);

   - kError       if there was an error during the optimization; this
                  typically indicates that computation cannot be resumed,
                  although solver-dependent ways of dealing with
                  solver-dependent errors may exist.

   MCFClass has a protected \c int \c member \c status \c that can be used by
   derived classes to hold status information and that is returned by the
   standard implementation of this method. Note that \c status \c is an
   \c int \c and not an \c enum \c, and that an \c int \c is returned by this
   method, in order to allow the derived classes to extend the set of return
   values if they need to do so. */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading flow solution
    @{ */

   virtual void MCFGetX( FRow F ,  Index_Set nms = NULL ,
                         cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Write the optimal flow solution in the vector F[]. If nms == NULL, F[]
   will be in "dense" format, i.e., the flow relative to arc `i'
   (i in 0 .. m - 1) is written in F[ i ]. If nms != NULL, F[] will be in 
   "sparse" format, i.e., the indices of the nonzero elements in the flow
   solution are written in nms (that is then Inf<Index>()-terminated) and
   the flow value of arc nms[ i ] is written in F[ i ]. Note that nms is
   *not* guaranteed to be ordered. Also, note that, unlike MCFGetRC() and
   MCFGetPi() [see below], nms is an *output* of the method.

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cFRow MCFGetX( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal data structure containing the
   flow solution in "dense" format. Since this may *not always be available*,
   depending on the implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual bool HaveNewX( void )
   {
    return( false );
    }

/**< Return true if a different (approximately) optimal primal solution is
   available. If the method returns true, then any subsequent call to (any
   form of) MCFGetX() will return a different primal solution w.r.t. the one
   that was being returned *before* the call to HaveNewX(). This solution need
   not be optimal (although, ideally, it has to be "good); this can be checked
   by comparing its objective function value, that will be returned by a call
   to MCFGetFO() [see below].

   Any subsequent call of HaveNewX() that returns true produces a new
   solution, until the first that returns false; from then on, no new
   solutions will be generated until something changes in the problem's
   data.

   Note that a default implementation of HaveNewX() is provided which is
   good for those solvers that only produce one optimal primal solution. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading potentials
    @{ */

   virtual void MCFGetPi(  CRow P ,  cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Writes the optimal node potentials in the vector P[]. If nms == NULL,
   the node potential of node `i' (i in 0 .. n - 1) is written in P[ i ]
   (note that here node names always start from zero, regardless to the value
   of USENAME0). If nms != NULL, it must point to a vector of indices in
   0 .. n - 1 (ordered in increasing sense and Inf<Index>()-terminated), and
   the node potential of nms[ i ] is written in P[ i ]. Note that, unlike
   MCFGetX() above, nms is an *input* of the method.

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the nodes `i' with strt <= i < min( MCFn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to nodes which are *both* in nms[] and whose index is
   in the correct range are returned. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cCRow MCFGetPi( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal data structure containing the
   node potentials. Since this may *not always be available*, depending on
   the implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual bool HaveNewPi( void )
   {
    return( false );
    }

/**< Return true if a different (approximately) optimal dual solution is
   available. If the method returns true, then any subsequent call to (any
   form of) MCFGetPi() will return a different dual solution, and MCFGetRC()
   [see below] will return the corresponding reduced costs. The new solution
   need not be optimal (although, ideally, it has to be "good); this can be
   checked by comparing its objective function value, that will be returned
   by a call to MCFGetDFO() [see below].

   Any subsequent call of HaveNewPi() that returns true produces a new
   solution, until the first that returns false; from then on, no new
   solutions will be generated until something changes in the problem's
   data.

   Note that a default implementation of HaveNewPi() is provided which is
   good for those solvers that only produce one optimal dual solution. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading reduced costs
    @{ */

   virtual void MCFGetRC(  CRow CR ,  cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Write the reduced costs corresponding to the current dual solution in
   RC[]. If nms == NULL, the reduced cost of arc `i' (i in 0 .. m - 1) is
   written in RC[ i ]; if nms != NULL, it must point to a vector of indices
   in 0 .. m - 1 (ordered in increasing sense and Inf<Index>()-terminated),
   and the reduced cost of arc nms[ i ] is written in RC[ i ]. Note that,
   unlike MCFGetX() above, nms is an *input* of the method.

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to arcs which are *both* in nms[] and whose index is
   in the correct range are returned.

   @note the output of MCFGetRC() will change after any call to HaveNewPi()
         [see above] which returns true. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cCRow MCFGetRC( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal data structure containing the
   reduced costs. Since this may *not always be available*, depending on the
   implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method.

   @note the output of MCFGetRC() will change after any call to HaveNewPi()
         [see above] which returns true. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual CNumber MCFGetRC( cIndex i ) = 0;

/**< Return the reduced cost of the i-th arc. This information should be
   cheapily available in most implementations.

   @note the output of MCFGetRC() will change after any call to HaveNewPi()
         [see above] which returns true. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading the objective function value
   @{ */

   virtual FONumber MCFGetFO( void ) = 0;

/**< Return the objective function value of the primal solution currently
   returned by MCFGetX().
   
   If MCFGetStatus() == kOK, this is guaranteed to be the optimal objective
   function value of the problem (to within the optimality tolerances), but
   only prior to any call to HaveNewX() that returns true. MCFGetFO()
   typically returns Inf<FONumber>() if MCFGetStatus() == kUnfeasible and
   - Inf<FONumber>() if MCFGetStatus() == kUnbounded. If MCFGetStatus() ==
   kStopped and MCFGetFO() returns a finite value, it must be an upper bound
   on the optimal objective function value (typically, the objective function
   value of one primal feasible solution). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual FONumber MCFGetDFO( void )
   {
    switch( MCFGetStatus() ) {
     case( kUnSolved ):
     case( kStopped ):
     case( kError ):    return( -Inf<FONumber>() );
     default:           return( MCFGetFO() );
     }
    }

/**< Return the objective function value of the dual solution currently
   returned by MCFGetPi() / MCFGetRC(). This value (possibly) changes after
   any call to HaveNewPi() that returns true. The relations between
   MCFGetStatus() and MCFGetDFO() are analogous to these of MCFGetFO(),
   except that a finite value corresponding to kStopped must be a lower
   bound on the optimal objective function value (typically, the objective
   function value one dual feasible solution).

   A default implementation is provided for MCFGetDFO(), which is good for
   MCF solvers where the primal and dual optimal solution values always
   are identical (except if the problem is unfeasible/unbounded). */

/*@} -----------------------------------------------------------------------*/
/** @name Getting unfeasibility certificate
   @{ */

   virtual FNumber MCFGetUnfCut( Index_Set Cut )
   {
    *Cut = Inf<Index>();
    return( 0 );
    }

/**< Return an unfeasibility certificate. In an unfeasible MCF problem,
   unfeasibility can always be reduced to the existence of a cut (subset of
   nodes of the graph) such as either:

   - the inverse of the deficit of the cut (the sum of all the deficits of the
     nodes in the cut) is larger than the forward capacity of the cut (sum of
     the capacities of forward arcs in the cut); that is, the nodes in the
     cut globally produce more flow than can be routed to sinks outside the
     cut;

   - the deficit of the cut is larger than the backward capacity of the cut
     (sum of the capacities of backward arcs in the cut); that is, the nodes
     in the cut globally require more flow than can be routed to them from
     sources outside the cut.

   When detecting unfeasibility, MCF solvers are typically capable to provide
   one such cut. This information can be useful - typically, the only way to
   make the problem feasible is to increase the capacity of at least one of
   the forward/backward arcs of the cut -, and this method is provided for
   getting it. It can be called only if MCFGetStatus() == kUnfeasible, and
   should write in Cut the set of names of nodes in the unfeasible cut (note
   that node names depend on USENAME0), Inf<Index>()-terminated, returning the
   deficit of the cut (which allows to distinguish which of the two cases
   above hold). In general, no special properties can be expected from the
   returned cut, but solvers should be able to provide e.g. "small" cuts.

   However, not all solvers may be (easily) capable of providing this
   information; thus, returning 0 (no cut) is allowed, as in the base class
   implementation, to signify that this information is not available. */

/*@} -----------------------------------------------------------------------*/
/** @name Getting unboundedness certificate
   @{ */

   virtual Index MCFGetUnbCycl( Index_Set Pred , Index_Set ArcPred )
   {
    return( Inf<Index>() );
    }

/**< Return an unboundedness certificate. In an unbounded MCF problem,
   unboundedness can always be reduced to the existence of a directed cycle
   with negative cost and all arcs having infinite capacity.
   When detecting unboundedness, MCF solvers are typically capable to provide
   one such cycle. This information can be useful, and this method is provided
   for getting it. It can be called only if MCFGetStatus() == kUnbounded, and
   writes in Pred[] and ArcPred[], respectively, the node and arc predecessor
   function of the cycle. That is, if node `i' belongs to the cycle then
   `Pred[ i ]' contains the name of the predecessor of `j' of `i' in the cycle
   (note that node names depend on USENAME0), and `ArcPred[ i ]' contains the
   index of the arc joining the two (note that in general there may be
   multiple copies of each arc). Entries of the vectors for nodes not
   belonging to the cycle are in principle undefined, and the name of one node
   belonging to the cycle is returned by the method. 
   Note that if there are multiple cycles with negative costs this method
   will return just one of them (finding the cycle with most negative cost
   is an NO-hard problem), although solvers should be able to produce cycles
   with "large negative" cost.

   However, not all solvers may be (easily) capable of providing this
   information; thus, returning Inf<Index>() is allowed, as in the base class
   implementation, to signify that this information is not available. */

/*@} -----------------------------------------------------------------------*/
/** @name Saving/restoring the state of the solver
    @{ */

   virtual MCFStatePtr MCFGetState( void )
   {
    return( NULL );
    }

 /**< Save the state of the MCF solver. The MCFClass interface supports the
   notion of saving and restoring the state of the MCF solver, such as the
   current/optimal basis in a simplex solver. The "empty" class MCFState is
   defined as a placeholder for state descriptions.

   MCFGetState() creates and returns a pointer to an object of (a proper
   derived class of) class MCFState which describes the current state of the
   MCF solver. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void MCFPutState( MCFStatePtr S ) {}

/**< Restore the solver to the state in which it was when the state `S' was
   created with MCFGetState() [see above].

   The typical use of this method is the following: a MCF problem is solved
   and the "optimal state" is set aside. Then the problem changes and it is
   re-solved. Then, the problem has to be changed again to a form that is
   close to the original one but rather different from the second one (think
   of a long backtracking in a Branch & Bound) to which the current "state"
   referes. Then, the old optimal state can be expected to provide a better
   starting point for reoptimization [see ReOptimize() below].

   Note, however, that the state is only relative to the optimization process,
   i.e., this operation is meaningless if the data of the problem has changed
   in the meantime. So, if a state has to be used for speeding up
   reoptimization, the following has to be done:

   - first, the data of the solver is brought back to *exactly* the same as
     it was at the moment where the state `S' was created (prior than this
     operation a ReOptimize( false ) call is probably advisable);

   - then, MCFPutState() is called (and ReOptimize( true ) is called);

   - only afterwards the data of the problem is changed to the final state
     and the problem is solved.

   A "put state" operation does not "deplete" the state, which can therefore
   be used more than once. Indeed, a state is constructed inside the solver
   for each call to MCFGetState(), but the solver never deletes statuses;
   this has to be done on the outside when they are no longer needed (the
   solver must be "resistent" to deletion of the state at any moment).

   Since not all the MCF solvers reoptimize (efficiently enough to make these
   operations worth), an "empty" implementation that does nothing is provided
   by the base class. */

/*@} -----------------------------------------------------------------------*/
/** @name Time the code
    @{ */

   void TimeMCF( double &t_us , double &t_ss )
   {
    t_us = t_ss = 0;
    if( MCFt )
     MCFt->Read( t_us , t_ss ); 
    }

/**< Time the code. If called within any of the methods of the class that are
   "actively timed" (this depends on the subclasses), this method returns the
   user and sistem time (in seconds) since the start of that method. If
   methods that are actively timed call other methods that are actively
   timed, TimeMCF() returns the (...) time since the beginning of the
   *outer* actively timed method. If called outside of any actively timed
   method, this method returns the (...) time spent in all the previous
   executions of all the actively timed methods of the class.

   Implementing the proper calls to MCFt->Start() and MCFt->Stop() is due to
   derived classes; these should at least be placed at the beginning and at
   the end, respectively, of SolveMCF() and presumably the Chg***() methods,
   that is, at least these methods should be "actively timed". */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   double TimeMCF( void )
   {
    return( MCFt ? MCFt->Read() : 0 );
    }

/**< Like TimeMCF(double,double) [see above], but returns the total time. */

/*@} -----------------------------------------------------------------------*/
/** @name Check the solutions
   @{ */

   inline void CheckPSol( void );

/**< Check that the primal solution returned by the solver is primal feasible.
   (to within the tolerances set by SetPar(kEps****) [see above], if any). Also,
   check that the objective function value is correct.

   This method is implemented by the base class, using the above methods
   for collecting the solutions and the methods of the next section for
   reading the data of the problem; as such, they will work for any derived
   class that properly implements all these methods. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   inline void CheckDSol( void );

/**< Check that the dual solution returned by the solver is dual feasible.
   (to within the tolerances set by SetPar(kEps****) [see above], if any). Also,
   check that the objective function value is correct.

   This method is implemented by the base class, using the above methods
   for collecting the solutions and the methods of the next section for
   reading the data of the problem; as such, they will work for any derived
   class that properly implements all these methods. */

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading graph size
    @{ */

   inline Index MCFnmax( void )
   {
    return( nmax );
    }

/**< Return the maximum number of nodes for this instance of MCFClass.
   The implementation of the method in the base class returns the protected
   fields \c nmax, which is provided for derived classes to hold this
   information. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   inline Index MCFmmax( void )
   {
    return( mmax );
    }

/**< Return the maximum number of arcs for this instance of MCFClass.
   The implementation of the method in the base class returns the protected
   fields \c mmax, which is provided for derived classes to hold this
   information. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   inline Index MCFn( void )
   {
    return( n );
    }

/**< Return the number of nodes in the current graph.
   The implementation of the method in the base class returns the protected
   fields \c n, which is provided for derived classes to hold this
   information. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   inline Index MCFm( void )
   {
    return( m );
    }

/**< Return the number of arcs in the current graph.
   The implementation of the method in the base class returns the protected
   fields \c m, which is provided for derived classes to hold this
   information. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading graph topology
   @{ */

   virtual void MCFArcs( Index_Set Startv ,  Index_Set Endv ,
                         cIndex_Set nms = NULL , cIndex strt = 0 ,
                         Index stp = Inf<Index>() ) = 0;

/**< Write the starting (tail) and ending (head) nodes of the arcs in Startv[]
   and Endv[]. If nms == NULL, then the information relative to all arcs is
   written into Startv[] and Endv[], otherwise Startv[ i ] and Endv[ i ]
   contain the information relative to arc nms[ i ] (nms[] must be
   Inf<Index>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to arcs which are *both* in nms and whose index is in
   the correct range are returned.

   Startv or Endv can be NULL, meaning that only the other information is
   required.

   @note If USENAME0 == 0 then the returned node names will be in the range
         1 .. n, while if USENAME0 == 1 the returned node names will be in
         the range 0 .. n - 1.

   @note If the graph is "dynamic", be careful to use MCFn() e MCFm() to
         properly choose the dimension of nodes and arcs arrays. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual Index MCFSNde( cIndex i ) = 0;

/**< Return the starting (tail) node of the arc `i'.

   @note If USENAME0 == 0 then the returned node names will be in the range
         1 .. n, while if USENAME0 == 1 the returned node names will be in
         the range 0 .. n - 1. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual Index MCFENde( cIndex i ) = 0;

/**< Return the ending (head) node of the arc `i'.

   @note If USENAME0 == 0 then the returned node names will be in the range
         1 .. n, while if USENAME0 == 1 the returned node names will be in
         the range 0 .. n - 1. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cIndex_Set MCFSNdes( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the starting
   (tail) nodes for each arc. Since this may *not always be available*,
   depending on the implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cIndex_Set MCFENdes( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the ending
   (head) nodes for each arc. Since this may *not always be available*,
   depending on the implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading arc costs
   @{ */

   virtual void MCFCosts( CRow Costv , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Write the arc costs into Costv[]. If nms == NULL, then all the costs are
   written, otherwise Costv[ i ] contains the information relative to arc
   nms[ i ] (nms[] must be Inf<Index>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to arcs which are *both* in nms and whose index is in
   the correct range are returned. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual CNumber MCFCost( cIndex i ) = 0;

/**< Return the cost of the i-th arc. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cCRow MCFCosts( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the arc
   costs. Since this may *not always be available*, depending on the
   implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void MCFQCoef( CRow Qv , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() )

/**< Write the quadratic coefficients of the arc costs into Qv[]. If
   nms == NULL, then all the costs are written, otherwise Costv[ i ] contains
   the information relative to arc nms[ i ] (nms[] must be
   Inf<Index>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to arcs which are *both* in nms and whose index is in
   the correct range are returned.

   Note that the method is *not* pure virtual: an implementation is provided
   for "pure linear" MCF solvers that only work with all zero quadratic
   coefficients. */
   {
    if( nms ) {
     while( *nms < strt )
      nms++;

     for( Index h ; ( h = *(nms++) ) < stp ; )
      *(Qv++) = 0;
     }
    else {
     if( stp > m )
      stp = m;

     while( stp-- > strt )
      *(Qv++) = 0;
     }
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual CNumber MCFQCoef( cIndex i )

/**< Return the quadratic coefficients of the cost of the i-th arc. Note that
   the method is *not* pure virtual: an implementation is provided for "pure
   linear" MCF solvers that only work with all zero quadratic coefficients. */
   {
    return( 0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cCRow MCFQCoef( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the arc
   costs. Since this may *not always be available*, depending on the
   implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready (such as "pure linear" MCF solvers that only
   work with all zero quadratic coefficients) does not need to implement the
   method. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading arc capacities
    @{ */

   virtual void MCFUCaps( FRow UCapv , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Write the arc capacities into UCapv[]. If nms == NULL, then all the
   capacities are written, otherwise UCapv[ i ] contains the information
   relative to arc nms[ i ] (nms[] must be Inf<Index>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the arcs `i' with strt <= i < min( MCFm() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to arcs which are *both* in nms and whose index is in
   the correct range are returned. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual FNumber MCFUCap( cIndex i ) = 0;

/**< Return the capacity of the i-th arc. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cFRow MCFUCaps( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the arc
   capacities. Since this may *not always be available*, depending on the
   implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method. */

/*@} -----------------------------------------------------------------------*/
/** @name Reading node deficits
    @{ */

   virtual void MCFDfcts( FRow Dfctv , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Write the node deficits into Dfctv[]. If nms == NULL, then all the
   defcits are written, otherwise Dfctvv[ i ] contains the information
   relative to node nms[ i ] (nms[] must be Inf<Index>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the nodes `i' with strt <= i < min( MCFn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != NULL then only the
   values corresponding to nodes which are *both* in nms and whose index is in
   the correct range are returned.

   @note Here node "names" (strt and stp, those contained in nms[] or `i'
         in MCFDfct()) go from 0 to n - 1, regardless to the value of
         USENAME0; hence, if USENAME0 == 0 then the first node is "named 1",
         but its deficit is returned by MCFDfcts( Dfctv , NULL , 0 , 1 ). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual FNumber MCFDfct( cIndex i ) = 0;

/**< Return the deficit of the i-th node.

   @note Here node "names" go from 0 to n - 1, regardless to the value of
         USENAME0; hence, if USENAME0 == 0 then the first node is "named 1",
         but its deficit is returned by MCFDfct( 0 ). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual cFRow MCFDfcts( void )
   {
    return( NULL );
    }

/**< Return a read-only pointer to an internal vector containing the node
   deficits.  Since this may *not always be available*, depending on the
   implementation, this method can (uniformly) return NULL.
   This is done by the base class already, so a derived class that does not
   have the information ready does not need to implement the method.

   @note Here node "names" go from 0 to n - 1, regardless to the value of
         USENAME0; hence, if USENAME0 == 0 then the first node is "named 1",
         but its deficit is contained in MCFDfcts()[ 0 ]. */

/*@} -----------------------------------------------------------------------*/
/** @name Write problem to file
    @{ */

   virtual inline void WriteMCF( ostream &oStrm , int frmt = 0 );

/**< Write the current MCF problem to an ostream. This may be useful e.g. for
   debugging purposes.

   The base MCFClass class provides output in two different formats, depending
   on the value of the parameter frmt:

   - kDimacs  the problem is written in DIMACS standard format, read by most
              MCF codes available;

   - kMPS     the problem is written in the "modern version" (tab-separated)
              of the MPS format, read by most LP/MIP solvers;

    - kFWMPS  the problem is written in the "old version" (fixed width
              fields) of the MPS format; this is read by most LP/MIP
              solvers, but some codes still require the old format.

   The implementation of WriteMCF() in the base class uses all the above
   methods for reading the data; as such it will work for any derived class
   that properly implements this part of the interface, but it may not be
   very efficient. Thus, the method is virtual to allow the derived classes
   to either re-implement WriteMCF() for the above two formats in a more
   efficient way, and/or to extend it to support other solver-specific
   formats.

   @note None of the above two formats supports quadratic MCFs, so if
         nonzero quadratic coefficients are present, they are just ignored.
   */

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the costs
    @{ */

   virtual void ChgCosts( cCRow NCost , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Change the arc costs. In particular, change the costs that are:

   - listed in into the vector of indices `nms' (ordered in increasing
     sense and Inf<Index>()-terminated),

   - *and* whose name belongs to the interval [`strt', `stp').

   That is, if strt <= nms[ i ] < stp, then the nms[ i ]-th cost will be
   changed to NCost[ i ]. If nms == NULL (as the default), *all* the
   entries in the given range will be changed; if stp > MCFm(), then the
   smaller bound is used.

   @note changing the costs of arcs that *do not exist* is *not allowed*;
         only arcs which have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChgCost( Index arc , cCNumber NCost ) = 0;

/**< Change the cost of the i-th arc.

   @note changing the costs of arcs that *do not exist* is *not allowed*;
         only arcs which have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChgQCoef( cCRow NQCoef = NULL , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() )

/**< Change the quadratic coefficients of the arc costs. In particular,
   change the coefficients that are:

   - listed in into the vector of indices `nms' (ordered in increasing
     sense and Inf<Index>()-terminated),

   - *and* whose name belongs to the interval [`strt', `stp').

   That is, if strt <= nms[ i ] < stp, then the nms[ i ]-th cost will be
   changed to NCost[ i ]. If nms == NULL (as the default), *all* the
   entries in the given range will be changed; if stp > MCFm(), then the
   smaller bound is used. If NQCoef == NULL, all the specified coefficients
   are set to zero.

   Note that the method is *not* pure virtual: an implementation is provided
   for "pure linear" MCF solvers that only work with all zero quadratic
   coefficients.

   @note changing the costs of arcs that *do not exist* is *not allowed*;
         only arcs which have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */
   {
    if( NQCoef )
     throw( MCFException( "ChgQCoef: nonzero coefficients not allowed" ) );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChgQCoef(  Index arc , cCNumber NQCoef )

/**< Change the quadratic coefficient of the cost of the i-th arc.

   Note that the method is *not* pure virtual: an implementation is provided
   for "pure linear" MCF solvers that only work with all zero quadratic
   coefficients.

   @note changing the costs of arcs that *do not exist* is *not allowed*;
         only arcs which have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */
   {
    if( NQCoef )
     throw( MCFException( "ChgQCoef: nonzero coefficients not allowed" ) );
    }

/*@} -----------------------------------------------------------------------*/
/** @name Changing the capacities
    @{ */

   virtual void ChgUCaps( cFRow NCap , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Change the arc capacities. In particular, change the capacities that are:

   - listed in into the vector of indices `nms' (ordered in increasing sense
     and Inf<Index>()-terminated),

   - *and* whose name belongs to the interval [`strt', `stp').

   That is, if strt <= nms[ i ] < stp, then the nms[ i ]-th capacity will be
   changed to NCap[ i ]. If nms == NULL (as the default), *all* the
   entries in the given range will be changed; if stp > MCFm(), then the
   smaller bound is used.

   @note changing the capacities of arcs that *do not exist* is *not allowed*;
         only arcs that have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChgUCap( Index arc , cFNumber NCap ) = 0;

/**< Change the capacity of the i-th arc.

   @note changing the capacities of arcs that *do not exist* is *not allowed*;
         only arcs that have not been closed/deleted [see CloseArc() /
         DelArc() below and LoadNet() above about C_INF costs] can be
         touched with these methods. */

/*@} -----------------------------------------------------------------------*/
/** @name Changing the deficits
    @{ */

   virtual void ChgDfcts( cFRow NDfct , cIndex_Set nms = NULL ,
                          cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Change the node deficits. In particular, change the deficits that are:

   - listed in into the vector of indices `nms' (ordered in increasing sense
     and Inf<Index>()-terminated),

   - *and* whose name belongs to the interval [`strt', `stp').

   That is, if strt <= nms[ i ] < stp, then the nms[ i ]-th deficit will be
   changed to NDfct[ i ]. If nms == NULL (as the default), *all* the
   entries in the given range will be changed; if stp > MCFn(), then the
   smaller bound is used.

   Note that, in ChgDfcts(), node "names" (strt, stp or those contained
   in nms[]) go from 0 to n - 1, regardless to the value of USENAME0; hence,
   if USENAME0 == 0 then the first node is "named 1", but its deficit can be
   changed e.g. with ChgDfcts( &new_deficit , NULL , 0 , 1 ).

   @note changing the capacities of nodes that *do not exist* is *not
         allowed*; only nodes that have not been deleted [see DelNode()
         below] can be touched with these methods. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChgDfct( Index node , cFNumber NDfct ) = 0;

/**< Change the deficit of the i-th node.

   Note that, in ChgDfct[s](), node "names" (i, strt/ stp or those contained
   in nms[]) go from 0 to n - 1, regardless to the value of USENAME0; hence,
   if USENAME0 == 0 then the first node is "named 1", but its deficit can be
   changed e.g. with ChgDfcts( 0 , new_deficit ).

   @note changing the capacities of nodes that *do not exist* is *not
         allowed*; only nodes that have not been deleted [see DelNode()
         below] can be touched with these methods. */

/*@} -----------------------------------------------------------------------*/
/** @name Changing graph topology
    @{ */

   virtual void CloseArc( cIndex name ) = 0;

/**< "Close" the arc `name'. Although all the associated information (name,
   cost, capacity, end and start node) is kept, the arc is removed from the
   problem until OpenArc( i ) [see below] is called.

   "closed" arcs always have 0 flow, but are otherwise counted as any other
   arc; for instance, MCFm() does *not* decrease as an effect of a call to
   CloseArc(). How this closure is implemented is solver-specific. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual bool IsClosedArc( cIndex name ) = 0;

/**< IsClosedArc() returns true if and only if the arc `name' is closed. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void DelNode( cIndex name ) = 0;

/**< Delete the node `name'.

   For any value of `name', all incident arcs to that node are closed [see
   CloseArc() above] (*not* Deleted, see DelArc() below) and the deficit is
   set to zero.

   Il furthermore `name' is the last node, the number of nodes as reported by
   MCFn() is reduced by at least one, until the n-th node is not a deleted
   one. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void OpenArc( cIndex name ) = 0;

/**< Restore the previously closed arc `name'. It is an error to open an arc
   that has not been previously closed. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual Index AddNode( cFNumber aDfct ) = 0;

/**< Add a new node with deficit aDfct, returning its name. Inf<Index>() is
   returned if there is no room for a new node. Remember that the node names
   are either { 0 .. nmax - 1 } or { 1 .. nmax }, depending on the value of
   USENAME0. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void ChangeArc( cIndex name , cIndex nSN = Inf<Index>() ,
                           cIndex nEN = Inf<Index>() ) = 0;

/**< Change the starting and/or ending node of arc `name' to nSN and nEN.
   Each parameter being Inf<Index>() means to leave the previous starting or
   ending node untouched. When this method is called `name' can be either the
   name of a "normal" arc or that of a "closed" arc [see CloseArc() above]:
   in the latter case, at the end of ChangeArc() the arc is *still closed*,
   and it remains so until OpenArc( name ) [see above] is called. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual void DelArc( cIndex name ) = 0;

/**< Delete the arc `name'. Unlike "closed" arcs, all the information
   associated with a deleted arc is lost and `name' is made available as a
   name for new arcs to be created with AddArc() [see below].

   Il furthermore `name' is the last arc, the number of arcs as reported by
   MCFm() is reduced by at least one, until the m-th arc is not a deleted
   one. Otherwise, the flow on the arc is always ensured to be 0. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual bool IsDeletedArc( cIndex name ) = 0;

/**< Return true if and only if the arc `name' is deleted. It should only
   be called with name < MCFm(), as every other arc is deleted by
   definition. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   virtual Index AddArc( cIndex Start , cIndex End , cFNumber aU ,
                         cCNumber aC ) = 0; 

/**< Add the new arc ( Start , End ) with cost aC and capacity aU,
   returning its name. Inf<Index>() is returned if there is no room for a new
   arc. Remember that arc names go from 0 to mmax - 1. */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~MCFClass()
   {
    delete MCFt;
    }

/**< Destructor of the class. The implementation in the base class only
   deletes the MCFt field. It is virtual, as it should be. */

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following fields are believed to be general enough to make it   --*/
/*--  worth adding them to the abstract base class.                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- MANAGING COMPARISONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Managing comparisons.
    The following methods are provided for making it easier to perform
    comparisons, with and without tolerances.
    @{ */

/*--------------------------------------------------------------------------*/
/** true if flow x is equal to zero (possibly considering tolerances). */
        
   template<class T>
   inline bool ETZ( T x , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x == 0 );
    else 
     return( (x <= eps ) && ( x >= -eps ) );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is greater than zero (possibly considering tolerances). */

   template<class T>
   inline bool GTZ( T x , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x > 0 );
    else
     return( x > eps );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is greater than or equal to zero (possibly considering
    tolerances). */

   template<class T>
   inline bool GEZ( T x , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x >= 0 );
    else
     return( x >= - eps );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is less than zero (possibly considering tolerances). */

   template<class T>
   inline bool LTZ( T x , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x < 0 );
    else
     return( x < - eps );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is less than or equal to zero (possibly considering
    tolerances). */

   template<class T>
   inline bool LEZ( T x , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x <= 0 );
    else
     return( x <= eps );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is greater than flow y (possibly considering tolerances).
 */

   template<class T>
   inline bool GT( T x , T y , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x > y );
    else
     return( x > y + eps );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** true if flow x is less than flow y (possibly considering tolerances). */

   template<class T>
   inline bool LT( T x , T y , const T eps )
   {
    if( numeric_limits<T>::is_integer )
     return( x < y );
    else
     return( x < y - eps );
    }

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED DATA STRUCTURES -------------------------*/
/*--------------------------------------------------------------------------*/

 Index n;      ///< total number of nodes
 Index nmax;   ///< maximum number of nodes

 Index m;      ///< total number of arcs
 Index mmax;   ///< maximum number of arcs

 int status;   ///< return status, see the comments to MCFGetStatus() above.
               ///< Note that the variable is defined int to allow derived
               ///< classes to return their own specialized status codes
 bool Senstv;  ///< true <=> the latest optimal solution should be exploited

 OPTtimers *MCFt;   ///< timer for performances evaluation

 FNumber EpsFlw;   ///< precision for comparing arc flows / capacities
 FNumber EpsDfct;  ///< precision for comparing node deficits
 CNumber EpsCst;   ///< precision for comparing arc costs

 double MaxTime;   ///< max time (in seconds) in which MCF Solver can find 
                   ///< an optimal solution (0 = no limits)
 int MaxIter;      ///< max number of iterations in which MCF Solver can find 
                   ///< an optimal solution (0 = no limits)

/*--------------------------------------------------------------------------*/

 };   // end( class MCFClass )

/*@} -----------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void MCFClass::LoadDMX( istream &DMXs , bool IsQuad )
{
 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 char c;
 for(;;) {
  if( ! ( DMXs >> c ) )
   throw( MCFException( "LoadDMX: error reading the input stream" ) );

  if( c != 'c' )  // if it's not a comment
   break;

  do              // skip the entire line
   if( ! DMXs.get( c ) )
    throw( MCFException( "LoadDMX: error reading the input stream" ) );
  while( c != '\n' );
  }

 if( c != 'p' )
  throw( MCFException( "LoadDMX: format error in the input stream" ) );

 char buf[ 4 ];    // free space
 DMXs >> buf;     // skip (in has to be "min")

 Index tn;
 if( ! ( DMXs >> tn ) )
  throw( MCFException( "LoadDMX: error reading number of nodes" ) );

 Index tm;
 if( ! ( DMXs >> tm ) )
  throw( MCFException( "LoadDMX: error reading number of arcs" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FRow      tU      = new FNumber[ tm ];  // arc upper capacities
 FRow      tDfct   = new FNumber[ tn ];  // node deficits
 Index_Set tStartn = new Index[ tm ];    // arc start nodes
 Index_Set tEndn   = new Index[ tm ];    // arc end nodes
 CRow      tC      = new CNumber[ tm ];  // arc costs
 CRow      tQ      = NULL;
 if( IsQuad )
  tQ = new CNumber[ tm ];                // arc quadratic costs


 for( Index i = 0 ; i < tn ; )           // all deficits are 0
  tDfct[ i++ ] = 0;                      // unless otherwise stated

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i = 0;  // arc counter
 for(;;) {
  if( ! ( DMXs >> c ) )
   break;

  switch( c ) {
   case( 'c' ):  // comment line- - - - - - - - - - - - - - - - - - - - - - -
    do
     if( ! DMXs.get( c ) )
      throw( MCFException( "LoadDMX: error reading the input stream" ) );
    while( c != '\n' );
    break;

   case( 'n' ):  // description of a node - - - - - - - - - - - - - - - - - -
    Index j;
    if( ! ( DMXs >> j ) )
     throw( MCFException( "LoadDMX: error reading node name" ) );

    if( ( j < 1 ) || ( j > tn ) )
     throw( MCFException( "LoadDMX: invalid node name" ) );

    FNumber Dfctj;
    if( ! ( DMXs >> Dfctj ) )
     throw( MCFException( "LoadDMX: error reading deficit" ) );

    tDfct[ j - 1 ] -= Dfctj;
    break;

   case( 'a' ):  // description of an arc - - - - - - - - - - - - - - - - - -
    if( i == tm )
     throw( MCFException( "LoadDMX: too many arc descriptors" ) );

	 //printf("%d\n",i);

    if( ! ( DMXs >> tStartn[ i ] ) )
     throw( MCFException( "LoadDMX: error reading start node" ) );

    if( ( tStartn[ i ] < 1 ) || ( tStartn[ i ] > tn ) )
     throw( MCFException( "LoadDMX: invalid start node" ) );

    if( ! ( DMXs >> tEndn[ i ] ) )
     throw( MCFException( "LoadDMX: error reading end node" ) );

    if( ( tEndn[ i ] < 1 ) || ( tEndn[ i ] > tn ) )
     throw( MCFException( "LoadDMX: invalid end node" ) );

    if( tStartn[ i ] == tEndn[ i ] )
     throw( MCFException( "LoadDMX: self-loops not permitted" ) );

    FNumber LB;
    if( ! ( DMXs >> LB ) )
     throw( MCFException( "LoadDMX: error reading lower bound" ) );

    if( ! ( DMXs >> tU[ i ] ) )
     throw( MCFException( "LoadDMX: error reading upper bound" ) );

    if( ! ( DMXs >> tC[ i ] ) )
     throw( MCFException( "LoadDMX: error reading arc cost" ) );

    if( tQ ) {
     if( ! ( DMXs >> tQ[ i ] ) )
      throw( MCFException( "LoadDMX: error reading arc quadratic cost" ) );

     if( tQ[ i ] < 0 )
      throw( MCFException( "LoadDMX: negative arc quadratic cost" ) );
     }

    if( tU[ i ] < LB )
     throw( MCFException( "LoadDMX: lower bound > upper bound" ) );

    tU[ i ] -= LB;
    tDfct[ tStartn[ i ] - 1 ] += LB;
    tDfct[ tEndn[ i ] - 1 ] -= LB;
    #if( USENAME0 )
     tStartn[ i ]--;  // in the DIMACS format, node names start from 1
     tEndn[ i ]--;
    #endif
    i++;
    break; 

   default:  // invalid code- - - - - - - - - - - - - - - - - - - - - - - - -
    throw( MCFException( "LoadDMX: invalid code" ) );

   }  // end( switch( c ) )
  }  // end( for( ever ) )

 if( i < tm )
  throw( MCFException( "LoadDMX: too few arc descriptors" ) );

 // call LoadNet- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LoadNet( tn , tm , tn , tm , tU , tC , tDfct , tStartn , tEndn );

 // then pass quadratic costs, if any

 if( tQ )
  ChgQCoef( tQ );

 // delete the original data structures - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] tQ;
 delete[] tC;
 delete[] tEndn;
 delete[] tStartn;
 delete[] tDfct;
 delete[] tU;

 }  // end( MCFClass::LoadDMX )

/*--------------------------------------------------------------------------*/

inline void MCFClass::SetPar( int par , int val )
{
 switch( par ) {
 case( kMaxIter ):
  MaxIter = val;
  break;

 case( kReopt ):
  Senstv = (val == kYes);
  break;
  
 default:
  throw( MCFException( "Error using SetPar: unknown parameter" ) );
  }
 }

/*--------------------------------------------------------------------------*/

inline void MCFClass::SetPar( int par, double val )
{
 switch( par ) {
 case( kEpsFlw ):
  if( EpsFlw != FNumber( val ) ) {
   EpsFlw = FNumber( val );
   EpsDfct = EpsFlw * ( nmax ? nmax : 100 );
   status = kUnSolved;
   }
  break;

 case( kEpsDfct ):
  if( EpsDfct != FNumber( val ) ) {
   EpsDfct = FNumber( val );
   status = kUnSolved;
   }
  break;

 case( kEpsCst ):
  if( EpsCst != CNumber( val ) ) {
   EpsCst = CNumber( val );
   status = kUnSolved;
   }
  break;
  
 case( kMaxTime ):
  MaxTime = val;
  break;

 default:
  throw( MCFException( "Error using SetPar: unknown parameter" ) );
  }
 }

/*--------------------------------------------------------------------------*/

inline void MCFClass::GetPar( int par, int &val )
{
 switch( par ) {
 case( kMaxIter ):
  val = MaxIter;
  break;

 case( kReopt ):
  if( Senstv ) val = kYes;
  else         val = kNo;
  break;
  
 default:
  throw( MCFException( "Error using GetPar: unknown parameter" ) );
  }
 }

/*--------------------------------------------------------------------------*/

inline void MCFClass::GetPar( int par , double &val )
{
 switch( par ) {
 case( kEpsFlw ):
  val = double( EpsFlw );
  break;

 case( kEpsDfct ):
  val = double( EpsDfct );
  break;

 case( kEpsCst ):
  val = double( EpsCst );
  break;
  
 case( kMaxTime ):
  val = MaxTime;
  break;

 default:
  throw( MCFException( "Error using GetPar: unknown parameter" ) );
  }
 }

/*--------------------------------------------------------------------------*/

inline void MCFClass::CheckPSol( void )
{
 FRow tB = new FNumber[ MCFn() ];
 MCFDfcts( tB );
 #if( ! USENAME0 )
  tB--;
 #endif

 FRow tX = new FNumber[ MCFm() ];
 MCFGetX( tX );
 FONumber CX = 0;
 for( Index i = 0 ; i < MCFm() ; i++ )
  if( GTZ( tX[ i ] , EpsFlw ) ) {
   if( IsClosedArc( i ) )
    throw( MCFException( "Closed arc with nonzero flow (CheckPSol)" ) );

   if( GT( tX[ i ] , MCFUCap( i ) , EpsFlw ) )
    throw( MCFException( "Arc flow exceeds capacity (CheckPSol)" ) );

   if( LTZ( tX[ i ] , EpsFlw ) )
    throw( MCFException( "Arc flow is negative (CheckPSol)" ) );

   CX += ( FONumber( MCFCost( i ) ) +
           FONumber( MCFQCoef( i ) ) * FONumber( tX[ i ] ) / 2 )
         * FONumber( tX[ i ] );
   tB[ MCFSNde( i ) ] += tX[ i ];
   tB[ MCFENde( i ) ] -= tX[ i ];
   }
 delete[] tX;
 #if( ! USENAME0 )
  tB++;
 #endif

 for( Index i = 0 ; i < MCFn() ; i++ )
  if( ! ETZ( tB[ i ] , EpsDfct ) )
   throw( MCFException( "Node is not balanced (CheckPSol)" ) );

 delete[] tB;
 CX -= MCFGetFO(); 
 if( ( CX >= 0 ? CX : - CX ) > EpsCst * MCFn() )
  throw( MCFException( "Objective function value is wrong (CheckPSol)" ) );

 }  // end( CheckPSol )

/*--------------------------------------------------------------------------*/

inline void MCFClass::CheckDSol( void )
{
 CRow tPi = new CNumber[ MCFn() ];
 MCFGetPi( tPi );

 FONumber BY = 0;
 for( Index i = 0 ; i < MCFn() ; i++ )
  BY += FONumber( tPi[ i ] ) * FONumber( MCFDfct( i ) );

 CRow tRC = new CNumber[ MCFm() ];
 MCFGetRC( tRC );

 FRow tX = new FNumber[ MCFm() ];
 MCFGetX( tX );

 #if( ! USENAME0 )
  tPi--;
 #endif

 FONumber QdrtcCmpnt = 0;  // the quadratic part of the objective
 for( Index i = 0 ; i < MCFm() ; i++ ) {
  if( IsClosedArc( i ) )
   continue;

  cFONumber tXi = FONumber( tX[ i ] );
  cCNumber QCi = CNumber( MCFQCoef( i ) * tXi );
  QdrtcCmpnt += QCi * tXi;
  CNumber RCi = MCFCost( i ) + QCi
                + tPi[ MCFSNde( i ) ] - tPi[ MCFENde( i ) ];
  RCi -= tRC[ i ];

  if( ! ETZ( RCi , EpsCst ) )
   throw( MCFException( "Reduced Cost value is wrong (CheckDSol)" ) );

  if( LTZ( tRC[ i ] , EpsCst ) ) {
   BY += FONumber( tRC[ i ] ) * FONumber( MCFUCap( i ) );

   if( LT( tX[ i ] , MCFUCap( i ) , EpsFlw ) )
    throw( MCFException( "Csomplementary Slackness violated (CheckDSol)" ) );
   }
  else
   if( GTZ( tRC[ i ] , EpsCst ) && GTZ( tX[ i ] , EpsFlw ) )
    throw( MCFException( "Complementary Slackness violated (CheckDSol)" ) );

  }  // end( for( i ) )

 delete[] tX;
 delete[] tRC;
 #if( ! USENAME0 )
  tPi++;
 #endif
 delete[] tPi;

 BY -= ( MCFGetFO() + QdrtcCmpnt / 2 );
 if( ( BY >= 0 ? BY : - BY ) > EpsCst * MCFn() )
  throw( MCFException( "Objective function value is wrong (CheckDSol)" ) );

 }  // end( CheckDSol )

/*--------------------------------------------------------------------------*/

inline void MCFClass::WriteMCF( ostream &oStrm , int frmt )
{
 if( ( ! numeric_limits<FNumber>::is_integer ) ||
     ( ! numeric_limits<CNumber>::is_integer ) )
  oStrm.precision( 12 );

 switch( frmt ) {
  case( kDimacs ):  // DIMACS standard format - - - - - - - - - - - - - - - -
                    //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case( kQDimacs ):  // ... linear or quadratic

   // compute and output preamble and size- - - - - - - - - - - - - - - - - -

   oStrm << "p min " << MCFn() << " ";
   {
    Index tm = MCFm();
    for(  Index i = MCFm() ; i-- ; )
     if( IsClosedArc( i ) || IsDeletedArc( i ) )
      tm--;

    oStrm << tm << endl;
    }

   // output arc information- - - - - - - - - - - - - - - - - - - - - - - - -

   for(  Index i = 0 ; i < MCFm() ; i++ )
    if( ( ! IsClosedArc( i ) ) && ( ! IsDeletedArc( i ) ) ) {
     oStrm << "a\t";
     #if( USENAME0 )
      oStrm << MCFSNde( i ) + 1 << "\t" << MCFENde( i ) + 1 << "\t";
     #else
      oStrm << MCFSNde( i ) << "\t" << MCFENde( i ) << "\t";
     #endif
     oStrm << "0\t" << MCFUCap( i ) << "\t" << MCFCost( i );

     if( frmt == kQDimacs )
      oStrm << "\t" << MCFQCoef( i );

     oStrm << endl;
     }

   // output node information - - - - - - - - - - - - - - - - - - - - - - - -

   for(  Index i = 0 ; i < MCFn() ; ) {
    cFNumber Dfcti = MCFDfct( i++ );
    if( Dfcti )
     oStrm << "n\t" << i << "\t" << - Dfcti << endl;
    }

   break;

  case( kMPS ):     // MPS format(s)- - - - - - - - - - - - - - - - - - - - -
  case( kFWMPS ):   //- - - - - - - - - - - - - - - - - - - - - - - - - - - -

   // writing problem name- - - - - - - - - - - - - - - - - - - - - - - - - -

   oStrm << "NAME      MCF" << endl;

   // writing ROWS section: flow balance constraints- - - - - - - - - - - - -

   oStrm << "ROWS" << endl << " N  obj" << endl;

   for(  Index i = 0 ; i < MCFn() ; )
    oStrm << " E  c" << i++ << endl;

   // writing COLUMNS section - - - - - - - - - - - - - - - - - - - - - - - -

   oStrm << "COLUMNS" << endl;

   if( frmt == kMPS ) {  // "modern" MPS
    for(  Index i = 0 ; i < MCFm() ; i++ )
     if( ( ! IsClosedArc( i ) ) && ( ! IsDeletedArc( i ) ) ) {
      oStrm << " x" << i << "\tobj\t" << MCFCost( i ) << "\tc";
      #if( USENAME0 )
       oStrm << MCFSNde( i ) << "\t-1" << endl << " x" << i << "\tc"
             << MCFENde( i ) << "\t1" << endl;
      #else
       oStrm << MCFSNde( i ) - 1 << "\t-1" << endl << " x" << i << "\tc"
             << MCFENde( i ) - 1 << "\t1" << endl;
      #endif
      }
     }
   else              // "old" MPS
    for(  Index i = 0 ; i < MCFm() ; i++ ) {
     ostringstream myField;
     myField << "x" << i << ends;
     oStrm << "    " << setw( 8 )  << left << myField.str()
           << "  "   << setw( 8 )  << left << "obj"
           << "  "   << setw( 12 ) << left << MCFCost( i );

     myField.seekp( 0 );
     myField << "c" << MCFSNde( i ) - ( 1 - USENAME0 ) << ends;

     oStrm << "   " << setw( 8 )  << left << myField.str()
           << "  "  << setw( 12 ) << left << -1 << endl;

     myField.seekp( 0 );
     myField << "x" << i << ends;

     oStrm << "    " << setw( 8 ) << left << myField.str();

     myField.seekp( 0 );
     myField << "c" << MCFENde( i ) - ( 1 - USENAME0 ) << ends;

     oStrm << "  " << setw( 8 ) << left << myField.str() << endl;
     }

   // writing RHS section: flow balance constraints - - - - - - - - - - - - -

   oStrm << "RHS" << endl;

   if( frmt == kMPS )  // "modern" MPS
    for( Index i = 0 ; i < MCFn() ; i++ )
     oStrm << " rhs\tc" << i << "\t" << MCFDfct( i ) << endl;
   else              // "old" MPS
    for( Index i = 0 ; i < MCFn() ; i++ ) {
     ostringstream myField;
     oStrm << "    " << setw( 8 ) << left << "rhs";
     myField << "c" << i << ends;

     oStrm << "  " << setw( 8 )  << left << myField.str()
           << "  " << setw( 12 ) << left << MCFDfct( i ) << endl;
     }

   // writing BOUNDS section- - - - - - - - - - - - - - - - - - - - - - - - -

   oStrm << "BOUNDS" << endl;

   if( frmt == kMPS ) {  // "modern" MPS
    for( Index i = 0 ; i < MCFm() ; i++ )
     if( ( ! IsClosedArc( i ) ) && ( ! IsDeletedArc( i ) ) )
      oStrm << " UP BOUND\tx" << i << "\t" << MCFUCap( i ) << endl;
    }
   else              // "old" MPS
    for(  Index i = 0 ; i < MCFm() ; i++ )
     if( ( ! IsClosedArc( i ) ) && ( ! IsDeletedArc( i ) ) ) {
      ostringstream myField;
      oStrm << " " << setw( 2 ) << left << "UP"
            << " " << setw( 8 ) << left << "BOUND";

      myField << "x" << i << ends;

      oStrm << "  " << setw( 8 )  << left << myField.str()
            << "  " << setw( 12 ) << left << MCFUCap( i ) << endl;
      }

   oStrm << "ENDATA" << endl;
   break;

  default:          // unknown format - - - - - - - - - - - - - - - - - - - -
                    //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   oStrm << "Error: unknown format " << frmt << endl;
  }
 }  // end( MCFClass::WriteMCF )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFClass.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFClass.h ----------------------------*/
/*--------------------------------------------------------------------------*/
