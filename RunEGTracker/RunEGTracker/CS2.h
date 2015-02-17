/*--------------------------------------------------------------------------*/
/*------------------------- File CS2.h -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Linear Min Cost Flow problems solver, based on the cs2-3.9 code by
 * Andrew Goldberg. Conforms to the standard MCF interface defined in
 * MCFClass.h.
 *
 * The CS2 solver is based on a cost-scaling, push-relabel algorithm.
 * This is a member of the family of primal-dual algorithms which
 * essentially operates as follows: a pseudoflow (a flow vector which
 * satisfies bound and non-negativity constraints but not necessarily
 * flow conservation constraints) is kept which satisfies
 * eps-complementarity slackness conditions with the current vector of
 * potentials; that is, only the flow on arcs whose reduced cost
 * \f[
 *  RC[ i , j ] =  C[ i , j ] - Pi[ j ] + Pi[ i ]
 * \f]
 * is <= eps in absolute value can be chosen to any value between 0 and the
 * capacity, while arcs with reduced cost < - eps are saturated (fixed to
 * their capacity) and arcs with reduced cost > eps are empty (fixed to 0).
 *
 * The algorithm attempts to convert the pseudoflow into a flow (i.e.,
 * to satisfy the flow conservation constraints) by essentially running
 * a max-flow algorithm of the push-relabel type. If the flow is found then
 * this is an eps-optimal solution of the problem; the algorithm is stopped
 * if eps is small enough (e.g., < 1 / n if costs are integer), otherwise
 * eps is decreased, thereby forcing some arcs to be either saturated or
 * emptyed to satisfy eps-complementarity slackness conditions with the new
 * eps. If no feasible flow is found, a saturated cut is identified which
 * separates the origins (nodes not yet producing enough flow) to the
 * destinations (nodes not yet consuming enough flow); this cut is used to
 * modify the potentials, thereby creating new arcs with reduced cost <=
 * eps in absolute value, which can be used to push further flow from the
 * origins to the destinations. If no such arcs can be created the problem
 * is declared unfeasible.
 *
 * \warning The original code has been written with integer data in mind.
 *          Therefore, costs are internally scaled by a factor of n, so to be
 *          able to work with integer potentials and reduced costs. This may
 *          create problems with the dual optimal solution when CNumber is an
 *          integer type, since to obtain a dual optimal solution for the
 *          original problem potentials and reduced costs have to be divided
 *          by n; with an integer data type the rounding of the integer
 *          division may lead to an unfeasible dual solution. If a feasible
 *          dual solution is critical, then CNumber has to be chosen as a
 *          float type even if the costs are integer. \n
 *          Furthermore, the internal scaling of costs has to be taken into
 *          account when selecting (an integer) CNumber; even if the costs
 *          of the original problem would allow using a "small" type (say,
 *          short), the internal scaling may require to use a "larger" type
 *          (say, int) in order to acommodate the larger scaled costs.
 *
 * \warning A private type SIndex is defined which is intended to hold arc
 *          and node indices "with a sign", used to represent orientation.
 *          This has to be "in sync" with Index, in the sense that for every
 *          unsigned index value in Index, the two signed values should be
 *          feasible in SIndex. In other words, either Index is not using at
 *          least half of its feasible values, or SIndex has to be a "bigger"
 *          data type than Index. The default value for SIndex is int.
 *
 * \version 1.40
 *
 * \date 11 - 02 - 2011
 *
 * \author <b>(original C code)</b> \n
 *         Andrew Goldberg \n
 *         http://www.avglab.com/andrew/
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Antonio Manca \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 1997 - 2011 by A. Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CS2
 #define __CS2

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CS2_MACROS Compile-time switches in CS2.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

/*----------------------------- DYNMC_MCF_CS2 ------------------------------*/

#define DYNMC_MCF_CS2 0

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_CS2 > 0, the methods of the public interface of class that
   allow to change the topology of the underlying network are actually
   implemented. Possible values of this macro are:

   - 0 => the topology of the graph cannot be changed;

   - 1 => the methods that change the topology of the graph are
          implemented. */

/*------------------------- CS2_STATISTICS ---------------------------------*/

#define CS2_STATISTICS 0

/**< If CS2_STATISTICS > 0, then statistic information about the behaviour of
   the cost-scaling algorithm is computed. */

/*@} -----------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CS2_CLASSES Classes in CS2.h
    @{ */

/** The CS2 class derives from the abstract base class MCFClass, thus sharing
    its (standard) interface, and implements a push-relabel cost-scaling
    algorithm for solving (Linear) Min Cost Flow problems */

class CS2 : public MCFClass {

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

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

  CS2( cIndex nmx = 0 , cIndex mmx = 0 );

/**< Constructor of the class, as in MCFClass::MCFClass(). */

/*--------------------------------------------------------------------------*/
/*---------------------- OTHER INITIALIZATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
		 cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
		 cIndex_Set pEn = NULL );

/* Inputs a new network, as in MCFClass::LoadNet().

   Passing pC[ i ] == Inf<CNumber>() means that the arc `i' does not exist in
   the problem. If DYNMC_MCF_CS2 > 0, these arcs are just "closed" and their
   cost is set to 0: this is done for being subsequently capable of
   "opening" them back with OpenArc(). */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

  void SolveMCF( void );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   void MCFGetX( FRow F , Index_Set nms = NULL ,
		 cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void MCFGetRC( CRow CR , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() ) ;

   inline CNumber MCFGetRC( cIndex i );

/*--------------------------------------------------------------------------*/

   void MCFGetPi( CRow P , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   CS2::FONumber MCFGetFO( void );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   virtual void MCFArcs( Index_Set Startv , Index_Set Endv ,
			 cIndex_Set nms = NULL , cIndex strt = 0 ,
			 Index stp = Inf<Index>() );

   inline Index MCFSNde( cIndex i );

   inline Index MCFENde( cIndex i );

/*--------------------------------------------------------------------------*/
  
   void MCFCosts( CRow Costv , cIndex_Set nms = NULL ,
                  cIndex strt = 0 , Index stp = Inf<Index>() );

   inline CNumber MCFCost( cIndex i );

/*--------------------------------------------------------------------------*/

   void MCFUCaps( FRow UCapv , cIndex_Set nms = NULL ,
	          cIndex strt = 0 , Index stp = Inf<Index>() ) ;

   inline FNumber MCFUCap( cIndex i );

/*--------------------------------------------------------------------------*/

   void MCFDfcts( FRow Dfctv , cIndex_Set nms = NULL ,
                  cIndex strt = 0 , Index stp = Inf<Index>() );

   inline FNumber MCFDfct( cIndex i );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgCost( Index arc , cCNumber NCost );

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms  ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgDfct( Index nod , cFNumber NDfct );

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap , cIndex_Set nms = NULL ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgUCap( Index arc , cFNumber NCap );

/*--------------------------------------------------------------------------*/
/*--------------- Modifying the structure of the graph ---------------------*/
/*--------------------------------------------------------------------------*/

   void CloseArc( cIndex name );

   void DelNode( cIndex name );

   inline bool IsClosedArc( cIndex name );

   void OpenArc( cIndex name );

   Index AddNode( cFNumber aDfct );

   void ChangeArc( cIndex name ,
		   cIndex nSS = Inf<Index>() , cIndex nEN = Inf<Index>() );

   void DelArc( cIndex name );

   Index AddArc( cIndex Start , cIndex End , cFNumber aU , cCNumber aC );

   inline bool IsDeletedArc( cIndex name );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~CS2();

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

  typedef int SIndex;  // index with a sign

  struct node_st;   // forward declaration of node_st, as arc_st contains one
                    // node_st pointer (the head node)

  struct arc_st {   // structure describing arcs- - - - - - - - - - - - - - -
   FNumber  r_cap;          // residual capacity
   CNumber  cost;           // cost of the arc
   node_st *head;           // head node
   arc_st  *sister;         // opposite arc
   SIndex   position;       // arc position: i + 1 if the arc_st represents
                            // the "original" arc i, - i - 1, if it
                            // it represents the sister of i
   };

  struct node_st {  // structure describing nodes - - - - - - - - - - - - - -
   arc_st  *first;          // first outgoing arc
   arc_st  *current;        // current outgoing arc
   arc_st  *suspended;      // first suspended arc
   #if( DYNMC_MCF_CS2 )
    arc_st *closed;         // first closed arc
   #endif
   FNumber  excess;         // excess of the node
   CNumber  price;          // distance from a sink
   node_st *q_next;         // next node in push queue
   node_st *b_next;         // next node in bucket-list
   node_st *b_prev;         // previous node in bucket-list
   SIndex   rank;           // bucket number
   char     inp;            // temporary number of input arcs
   };

  struct bucket_st {  // structure describing buckets - - - - - - - - - - - -
   node_st *p_first;        // 1st node with positive exces or simply 1st
                            // node in the bucket
   };

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void updtarccst( arc_st *arc , cCNumber NCa );

   // updates the cost of arc to NCa, handling reoptimization

/*--------------------------------------------------------------------------*/

   inline void updtarccap( arc_st *arc , cFNumber NCa );

   // updates the capacity of arc to NCa, handling reoptimization

/*--------------------------------------------------------------------------*/

   bool price_update( void );

/*--------------------------------------------------------------------------*/

   bool relabel( node_st *i );

/*--------------------------------------------------------------------------*/

   void discharge( node_st *i );

/*--------------------------------------------------------------------------*/

   Index price_in( void );

/*--------------------------------------------------------------------------*/

   void refine( void ); 

/*--------------------------------------------------------------------------*/

   bool price_refine( void );

/*--------------------------------------------------------------------------*/

   void compute_prices( void );

/*--------------------------------------------------------------------------*/

   void price_out( void );

/*--------------------------------------------------------------------------*/

   inline bool update_epsilon( void );

/*--------------------------------------------------------------------------*/

   inline CNumber REDUCED_COST( const node_st *i , const node_st *j ,
				const arc_st *a );

   // computes the reduced cost of arc (i, j) == a

/*--------------------------------------------------------------------------*/

   inline bool SUSPENDED( const node_st *i , const arc_st *a );

   // true if the arc a, outgoing from i, is suspended

/*--------------------------------------------------------------------------*/

   inline void EXCHANGE( arc_st* a , arc_st* b );

   // exchange arc position 

/*--------------------------------------------------------------------------*/

   inline void INCREASE_FLOW( node_st* i , node_st* j , arc_st* a ,
			      cFNumber df );

   // increas arc flow of an amount df

/*-------------------------------------------------------------------------*/

   void UPDATE_CUT_OFF( void );

/*-------------------------------------------------------------------------*/

   inline void RESET_EXCESS_Q( void );

   // reset excess nodes list

/*-------------------------------------------------------------------------*/

   inline void INSERT_TO_EXCESS_Q( node_st* i );

/*-------------------------------------------------------------------------*/

   inline void INSERT_TO_FRONT_EXCESS_Q( node_st* i );

   // insert node with positive excess in the excess nodes list

/*-------------------------------------------------------------------------*/

   inline void REMOVE_FROM_EXCESS_Q( node_st* &i );

   // remove node from excess nodes list

/*-------------------------------------------------------------------------*/

   inline void RESET_BUCKET( bucket_st *b );

   // empties bucket b

/*-------------------------------------------------------------------------*/

   inline bool NONEMPTY_BUCKET( bucket_st *b );

   // true if bucket b is nonempty

/*-------------------------------------------------------------------------*/

   inline void INSERT_TO_BUCKET( node_st* i , bucket_st *b );

   // insert node i to bucket b

/*-------------------------------------------------------------------------*/
 
   inline void GET_FROM_BUCKET( node_st* &i , bucket_st* b );

/*-------------------------------------------------------------------------*/

   inline void REMOVE_FROM_BUCKET( node_st* i , bucket_st* b );

/*-------------------------------------------------------------------------*/

   inline void STACKQ_PUSH( node_st* i );

/*--------------------------------------------------------------------------*/

   void MemAlloc();

   void MemDeAlloc();

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  FONumber ObjVal;          // value of the objective function

  bool reopt;               // true if only costs have changed, so the first
                            // part of the algorithm is skipped
  CNumber dn;               // number of nodes in CNumber

  node_st *nodes;           // array of nodes
  node_st *sentinel_node;   // next after last

  arc_st *arcs;             // array of arcs
       
  arc_st **pos;             // pos[ i ] = pointer to the arc_st representing
                            // the original arc i

  bucket_st *buckets;       // array of buckets
  bucket_st *l_bucket;      // last bucket

  #if( CS2_STATISTICS )
   long n_push;             // current number of push operations
   long n_relabel;          // current number of relabel operations
   long n_discharge;        // current number of discharge operations
   long n_refine;           // current number of refine operations
   long n_update;           // current number of update operations
   long n_scan;             // current number of up_node_scan calls
   long n_prscan;           // current number of bucket scans
   long n_prscan1;          // current number of stack_push
   long n_prscan2;          // current number of compute longest distances
   long n_prefine;          // current number of price refine
  #endif

  long n_bad_pricein;       // current number of recalculating excess queue
  long n_bad_relabel;       // current number of bad relabel
  long n_rel;               // number of relabels from last price update 
  long n_ref;               // current number of refines 
  char time_for_price_in;   // contains the value for the call to price_in 

  Index n_src;              // current number of nodes with positive excess 

  node_st *excq_first;      // first node in push-queue
  node_st *excq_last;       // last node in push-queue

  SIndex linf;              // number of buckets + 1

  CNumber m_c;              // max arc cost
  CNumber cut_on;           // the bound for returning suspended arcs
  CNumber cut_off;          // the bound for suspending arcs
  CNumber epsilon;          // optimality bound

  double cut_off_factor;    // multiplier to produce cut_on and cut_off from
                            // n and epsilon

  FNumber total_excess;     // total excess

  bool flag_price;          // true = signal to start price-in ASAP, maybe
                            // there is infeasibility because of suspended
                            // arcs

  long empty_push_bound;    // maximal possible number of zero pushes during
                            // one discharge

  Index Blncd;              // it's 1-value if optimal flow is balanced
                            // 0 otherwise.  

  CNumber PRICE_MIN;

  CNumber LOW_BOUND;

/*--------------------------------------------------------------------------*/

 };  // end( class CS2 )

/* @} end( group( CS2_CLASSES ) ) */
/*-------------------------------------------------------------------------*/
/*-------------------inline methods implementation-------------------------*/
/*-------------------------------------------------------------------------*/

inline MCFClass::Index CS2::MCFSNde( cIndex i )
{
 return( ( ( pos[ i ]->sister )->head - nodes ) - USENAME0 );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::Index CS2::MCFENde( cIndex i )
{
 return( ( pos[ i ]->head - nodes ) - USENAME0 );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::CNumber CS2::MCFCost( cIndex i )
{
 return( ( pos[ i ]->cost ) / dn );
 }

/*-------------------------------------------------------------------------*/

inline MCFClass::FNumber CS2::MCFUCap( cIndex i )
{
 return( pos[ i ]->r_cap + pos[ i ]->sister->r_cap );
 }

/*-------------------------------------------------------------------------*/

inline bool CS2::IsClosedArc( cIndex name )
{
 #if( DYNMC_MCF_CS2 )
  arc_st *arc = pos[ name ];
  return( arc <= (arc->sister)->head->suspended );
 #else
  return( false );
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool CS2::IsDeletedArc( cIndex name )
{
 return( CS2::IsClosedArc( name ) );  // limited implementation, on par with
                                      // the one of DelArc()
 }

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*-------------------------------------------------------------------------*/

#endif  /* CS2.h included */

/*-------------------------------------------------------------------------*/
/*---------------------- End File CS2.h -----------------------------------*/
/*-------------------------------------------------------------------------*/
