/*--------------------------------------------------------------------------*/
/*---------------------------- File CS2.C ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Linear Min Cost Flow problems solver, based on the cs2-3.9 code by  --*/
/*--  A. Goldberg. Conforms to the standard (MCF) interface defined in    --*/
/*--  MCFClass.h.                                                         --*/
/*--                                                                      --*/
/*--                            VERSION 1.40                              --*/
/*--                           11 - 02 - 2011                             --*/
/*--                                                                      --*/
/*--                   Original Idea and Implementation:                  --*/
/*--                                                                      --*/
/*--                            Andrew Goldberg                           --*/
/*--                      http://www.avglab.com/andrew/                   --*/
/*--                                                                      --*/
/*--                      C++ translation and polishing:                  --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                            Antonio Manca                             --*/
/*--                                                                      --*/
/*--                       Operations Research Group                      --*/
/*--                      Dipartimento di Informatica                     --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*--           Copyright (C) 1997 - 2011 by Antonio Frangioni             --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CS2.h"
#include <cmath>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ MACROS ------------------------------------*/
/*--------------------------------------------------------------------------*/

#define COMP_DUALS 1

// at the end of the algorithm, check that the solution is actually an
// epsilon-optimal one

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( ! DYNMC_MCF_CS2 )
#define closed suspended
#endif

// in general, the star of one node i is divided in the following subsets:
//  - i->closed <= a < i->suspended  ==> all the closed arcs
//  - i->suspended <= a < i->first   ==> all the suspended arcs
//  - i->first <= a < (i+1)->closed  ==> all the "active" arcs
//
// if DYNMC_MCF_CS2 == 0, the field "closed" is not defined for the arc_st
// descriptor, as i->closed == i->suspended; this is obtained by the above
// macro

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( CS2_STATISTICS )
#define Increase( a ) a++
#else
#define Increase( a )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- CONSTANTS ------------------------------------*/
/*--------------------------------------------------------------------------*/

static MCFClass::cIndex      UPDT_FREQ_S   = 30;
static MCFClass::cCNumber    SCALE_DEFAULT = 12;
static const MCFClass::Index EMPTY_PUSH_COEF = 1;

static const double  UPDT_FREQ          = 0.4;
static const long    PRICE_OUT_START    = 10;  // may not be less than 1

static const double  CUT_OFF_POWER      = 0.44;
static const double  CUT_OFF_COEF       = 1.5;
static const double  CUT_OFF_POWER2     = 0.75;
static const double  CUT_OFF_COEF2      = 1;
static const double  CUT_OFF_GAP        = 0.8;
static const double  CUT_OFF_MIN        = 12;
static const double  CUT_OFF_INCREASE   = 4;

static const char    TIME_FOR_PRICE_IN1 = 2;
static const char    TIME_FOR_PRICE_IN2 = 4;
static const char    TIME_FOR_PRICE_IN3 = 6;

static const long MAX_CYCLES_CANCELLED  = 0;
static const long START_CYCLE_CANCEL    = 100;

static const char WHITE = 0;
static const char GREY  = 1;
static const char BLACK = 2;

/*--------------------------------------------------------------------------*/
/*---------------------------- FUNCTIONS -----------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void Swap( T &v1 , T &v2 )
{
	const T temp = v1;

	v1 = v2;
	v2 = temp;
}

/*--------------------------------------------------------------------------*/
/*--------------------------- COSTRUCTOR -----------------------------------*/
/*--------------------------------------------------------------------------*/

CS2::CS2( cIndex nmx , cIndex mmx )
:
MCFClass( nmx , mmx )
{
	// allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if( nmax && mmax )
		MemAlloc();
	else
		nmax = mmax = 0;

	if( numeric_limits<CNumber>::is_integer ) {
		PRICE_MIN = - ( Inf<CNumber>() / 2 - 1 );
		LOW_BOUND = 1;
	}
	else {
		PRICE_MIN = - ( Inf<CNumber>() / 2 ) * ( CNumber( 1 ) - Eps<CNumber>() );
		LOW_BOUND = 1.00001;
	}
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void CS2::LoadNet( cIndex nmx , cIndex mmx , cIndex pn , cIndex pm ,
				  cFRow pU , cCRow pC , cFRow pDfct , cIndex_Set pSn ,
				  cIndex_Set pEn )
{
	// allocating and deallocating memory- - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	n = pn;
	m = pm;

	if( ( nmx != nmax ) || ( mmx != mmax ) ) {
		if( nmax && mmax ) {
			MemDeAlloc();
			nmax = mmax = 0;
		}

		if( mmx && nmx ) {
			nmax = nmx;
			mmax = mmx;
			MemAlloc();
		}
	}

	if( ( ! nmax ) || ( ! mmax ) ) {  // just sit down in the corner and wait
		nmax = mmax = 0;
		return;
	}

	// setting up the data structures- - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	sentinel_node = nodes + n + 1;

	dn = CNumber( n );

	// allocating temporary vectors for the ordering - - - - - - - - - - - - - -

	Index_Set arc_tail = new Index[ 2 * m ];

	Index_Set arc_first = new Index[ n + 2 ];
	for( Index i = 0 ; i < n + 2 ; )
		arc_first[ i++ ] = 0;

	// setting up arcs[] - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	m_c = 0;
	arc_st* arc = arcs;
	for( Index i = 0 ; i <  m ; i++ , arc++ ) {
		// read arc description - - - - - - - - - - - - - - - - - - - - - - - - - -

		Index tail = *(pSn++) + USENAME0;
		Index head = *(pEn++) + USENAME0;
		FNumber acap  = ( pU ? *(pU++) : Inf<FNumber>() );         
		CNumber dcost = ( pC ? *(pC++) : 0 );

		// update maximum cost- - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( dcost < Inf<CNumber>() ) {
			if( acap > 0 )
				if( dcost > m_c )
					m_c = dcost;
				else
					if( - dcost > m_c )
						m_c = - dcost;

			dcost *= dn;
		}

		// set the pos[] inverse map- - - - - - - - - - - - - - - - - - - - - - - -

		pos[ i ] = arc;

		// count this arc in the star of its head & tail- - - - - - - - - - - - - -
		// at the end, arc_first[ i + 1 ] = number of arcs outgoing from node "i"

#if( ! DYNMC_MCF_CS2 )
		if( dcost < Inf<CNumber>() )
#endif
		{
			arc_first[ tail + 1 ]++;
			arc_first[ head + 1 ]++;
		}

		// setting up the arc - - - - - - - - - - - - - - - - - - - - - - - - - - -

		arc_tail[ 2 * i ] = tail;
		arc->head = nodes + head;
		arc->r_cap = acap;

		arc->cost = dcost;
		arc->sister = arc + 1; 
		arc->position = i + 1;

		// setting up the sister arc- - - - - - - - - - - - - - - - - - - - - - - -

		( arc + 1 )->sister = arc;

		arc_tail[ 2 * i + 1 ] = head;
		( ++arc )->head = nodes + tail;
		arc->r_cap = 0;

		arc->cost = - dcost;
		arc->position = - ( i + 1 );

	}  // end( for( i ) )

	// reordering arcs[] - linear time algorithm - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// The arcs[] vector is ordered by (Forward) Star, i.e., first all the arcs
	// exiting from the first node, then all the arcs exiting from the second
	// and so on. The position of the first arc exiting one given node is
	// contained into node.first. Note that, since each arc is represented
	// twice (the original arc and its "sister"), the forward and backward star
	// of each node are mixed together.
	//
	// before the loop below, arc_first[ i + 1 ] is the number of arcs outgoing
	// from i; after the loop, arc_first[ i ] is the position of the first arc
	// outgoing from node i when the arcs are ordered; this value is transformed
	// to pointer and written to nodes[ i ].first

	for( Index i = 0 ; i++ <= n ; ) {
		arc_first[ i ] += arc_first[ i - 1 ];
		nodes[ i ].first = arcs + arc_first[ i ];
	}

	sentinel_node->suspended = sentinel_node->first;

#if( ! DYNMC_MCF_CS2 )
	// if DYNMC_MCF_CS2 == 0, arcs with INF cost are not put in the star of
	// their head/tail nodes, but rather put (as suspended arcs) in the star
	// of the dummy node n + 1- - - - - - - - - - - - - - - - - - - - - - - - -

	for( arc_st *arc = arcs ; arc < sentinel_node->suspended ; arc++ )
		while( arc->cost == Inf<CNumber>() )
			EXCHANGE( arc , (sentinel_node->first)++ );
#else
	sentinel_node->closed = sentinel_node->suspended;
#endif

	for( Index i = 1 ; i < n ; i++ ) {  // scanning all the nodes except the last
		cIndex last = ( ( nodes + i + 1 )->first ) - arcs;
		for( Index arc_num = arc_first[ i ] ; arc_num < last ; arc_num++ )
			for( Index tail = arc_tail[ arc_num ] ; tail != i ; ) {
				// the arc arc_num is not in place because the arc in this position must
				// go out from i; we'll put it to its place and continue this process
				// until an arc in this position would go out from i

				cIndex arc_new_num = arc_first[ tail ];
				arc_st *arc = arcs + arc_num;
				arc_st *arc_new = arcs + arc_new_num;

				// arc must be in position arc_new: swapping these arcs

				EXCHANGE( arc , arc_new );

				arc_tail[ arc_num ] = arc_tail[ arc_new_num ];
				arc_tail[ arc_new_num ] = tail;

				arc_first[ tail ]++;  // signal that this arc is in a correct position
				// and need not to be examined again
				tail = arc_tail[ arc_num ];
			}
	}

	// deallocating temporaries- - - - - - - - - - - - - - - - - - - - - - - - -

	delete[] arc_first;
	delete[] arc_tail;

#if( DYNMC_MCF_CS2 )
	// put all closed arcs at the beginning of the star of their nodes- - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	for( node_st *node = nodes ; ++node < sentinel_node ; ) {
		node->closed = node->first;

		for( arc_st *arc = node->first ; arc < (node + 1)->first ; arc++ )
			if( arc->cost == Inf<CNumber>() )
				if( arc == node->first )
					(node->first)++;
				else
					EXCHANGE( arc , ++(node->first) );
	}
#endif

	// setting up nodes[]- - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	node_st *node = nodes;
	node->suspended = NULL;
	node->first = NULL;
	node->b_next = NULL;
	node->b_prev = NULL;
	node->rank = 0;
	node->inp = 0;

	for( ; ++node < sentinel_node ; ) {
		node->suspended = node->first;
		node->excess = - *(pDfct++);
		node->b_next = NULL;
		node->b_prev = NULL;
		node->rank = 0;
		node->inp = 0;
	}

	node->b_next = NULL;
	node->b_prev = NULL;
	node->rank = 0;
	node->inp = 0;

	// setting up buckets[]- - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	for( bucket_st *b = buckets ; b < l_bucket ; b++ )
		RESET_BUCKET( b );

	// final initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( CS2_STATISTICS )
	n_push = n_relabel = n_discharge = n_refine = n_update = n_scan = 
		n_prscan = n_prscan1 = n_prscan2 = n_prefine = 0;
#endif

	empty_push_bound = n * EMPTY_PUSH_COEF;
	status = kUnSolved;

}  // end( CS2::LoadNet )

/*-------------------------------------------------------------------------*/
/*--------------- METHODS FOR SOLVING THE PROBLEM -------------------------*/
/*-------------------------------------------------------------------------*/

void CS2::SolveMCF( void )
{
	ObjVal = Inf<CS2::FONumber>();

	if( MCFt )
		MCFt->Start();

	// initialization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if( status == kUnSolved ) {  // construct an initial "empty" solution- - - -
		// set prices to 0

		for( node_st *i = nodes ; ++i < sentinel_node ; ) {
			i->price = 0;   
			i->first = i->suspended;
			i->current = i->first;
			i->q_next = sentinel_node;
		}

		// compute starting epsilon

		epsilon = m_c * dn; 
		if( epsilon < 1 )
			epsilon = 1;

		Blncd = false;
	}
	else {  // restart with the previous solution- - - - - - - - - - - - - - - -
		// make prices non negative

		CNumber rc = 0;
		for( node_st *i = nodes ; ++i < sentinel_node ; ) {
			if( rc < i->price)
				rc = i->price;

			i->first = i->suspended;
			i->current = i->first; 
			i->q_next = sentinel_node;
		}

		for( node_st *i = nodes ; ++i < sentinel_node ; )
			i->price -= rc;

		// compute starting epsilon

		if( Blncd ) {  // the flow is balanced
			epsilon = dn;

			// compute the maximum (in absolute value) of the reduced cost of arcs
			// violating the complementary slackness condition

			CNumber sum = 0;
			for( node_st *i = nodes ; ++i < sentinel_node ; ) {
				CNumber minc = 0;
				for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ;
					a++ )
					if( GTZ( a->r_cap , EpsFlw ) ) {
						cCNumber rc = REDUCED_COST( i , a->head , a );
						if( LTZ( rc , EpsCst ) )
							minc = max( minc , CNumber( - rc ) ) ;
					}

					sum += minc;
			}

#if( Ctype == REAL_TYPE )
			epsilon = ceil( sum / dn );
#else
			epsilon = sum / dn + ( sum % dn ? 1 : 0 );
#endif

			if( epsilon < 1 )
				epsilon = 1;
		}
		else {          // the flow may not be balanced
			epsilon = 1;

			// the new epsilon value is taken as the minimum reduced cost between
			// residual arcs

			for( node_st *i = nodes ; ++i < sentinel_node ; )
				for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
					if ( GTZ( a->r_cap , EpsFlw ) ) { 
						cCNumber rc = REDUCED_COST( i , a->head , a );
						if ( rc < -epsilon )
							epsilon = -rc ;   
					}

					// every arc with absolute cost less then 2 * dn * epsilon is considered
					// closed into subroutine refine

					if( 2 * epsilon < m_c )
						epsilon = m_c * dn * SCALE_DEFAULT;

					// relabel nodes with positive excess

					for( node_st *i = nodes ; ++i < sentinel_node ; )
						if( i->excess > 0 )
							relabel( i );

		}  // end( the flow may not be balanced )
	}  // end( restart with the previous solution ) - - - - - - - - - - - - - -

	cut_off_factor = CUT_OFF_COEF * pow( dn , CUT_OFF_POWER );
	cut_off_factor = max( cut_off_factor , CUT_OFF_MIN );

	n_ref = 0; 
	n_rel = 0;
	flag_price = false;
	n_bad_pricein = n_bad_relabel = 0;
	empty_push_bound = n * EMPTY_PUSH_COEF;

#if( CS2_STATISTICS )
	n_push = n_relabel = n_discharge = n_refine = n_update = n_scan =
		n_prscan = n_prscan1 = n_prscan2 = n_prefine = 0;
#endif

	excq_first = NULL;

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	status = kOK;
	bool cc = update_epsilon();

	if( Blncd ) {  // this is a reoptimization starting from a balanced flow - -
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// price_refine() can be directly called before refine()

		while( ! cc ) {  // scaling loop- - - - - - - - - - - - - - - - - - - - - -
			for(;;) {
				if( ! price_refine() )
					break;

				if( ( n_ref >= PRICE_OUT_START ) && price_in() )
					break;     

				if( ( cc = update_epsilon() ) )
					break;
			}

			if( cc )
				break;

			refine();

			if( status )  // problem unfeasible or error
				break;     

			if( n_ref >= PRICE_OUT_START )
				price_out();

			if( update_epsilon() )
				break;

			if( MaxIter && ( n_ref > MaxIter ) ) {  // iterations limit
				status = kStopped;
				break;
			}

			if( MCFt && MaxTime && ( MCFt->Read() > MaxTime ) ) {  // time limit
				status = kStopped;
				break;
			}
		}  // end( while( scaling loop ) ) - - - - - - - - - - - - - - - - - - - -
	}
	else  // starting from a possibly non balanced flow- - - - - - - - - - - - -
	{     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// refine() has to be called first

		do {  // scaling loop - - - - - - - - - - - - - - - - - - - - - - - - - - -
			refine();

			if( status )  // problem unfeasible or error
				break;  

			if( n_ref >= PRICE_OUT_START )
				price_out();

			if( update_epsilon() ) 
				break;   

			for(;;) {
				if( ! price_refine() )
					break;

				if( ( n_ref >= PRICE_OUT_START ) && price_in() )
					break;

				if( ( cc = update_epsilon() ) )
					break;
			}

			if( MaxIter && ( n_ref > MaxIter ) ) {  // iterations limit
				status = kStopped;
				break;
			}

			if( MCFt && MaxTime && ( MCFt->Read() > MaxTime ) ) {  // time limit
				status = kStopped;
				break;
			}
		} while( ! cc );
	}  // end( else( Blncd ) )- - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// final things- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// reset node prices in case of unfeasible or error

	if( status > 0 )
		for( node_st *i = nodes ; ++i < sentinel_node ; )
			i->price = 0;

#if( COMP_DUALS )
	{
		// for debugging purposes

		compute_prices();

		//for( node_st *i = nodes ; ++i < sentinel_node ; )
		//	for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ;
		//		a++ )
		//		if( GTZ( a->r_cap , EpsFlw ) )
		//			if( (i->price) + (a->cost) - (a->head)->price < - epsilon )
		//				throw( MCFException( "Error: wrong reduced cost" ) );
	}
#endif

	if( status )
		Blncd = true;

	if( MCFt )
		MCFt->Stop();

} // end( CS2::SolveMCF )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void CS2::MCFGetX( FRow F , Index_Set nms , cIndex strt , Index stp )
{
	if( stp > m )
		stp = m; 

	if( nms )
		for( Index i = strt ; i < stp ; i++ ) {
			cFNumber Xi = (pos[ i ]->sister)->r_cap;
			if( GTZ( Xi , EpsFlw ) ) {
				*(F++) = Xi;
				*(nms++) = i;
			}
		}
	else
		for( Index i = strt ; i < stp ; i++ )
			*(F++) = (pos[ i ]->sister)->r_cap;

}  // end( CS2::MCFGetX() )

/*--------------------------------------------------------------------------*/

void CS2::MCFGetRC( CRow CR , cIndex_Set nms , cIndex strt , Index stp  )
{
	if( nms ) {
		while( *nms < strt )
			nms++;

		for( Index h ; ( h = *(nms++) ) < stp ; )
			*(CR++) = CS2::MCFGetRC( h );
	}
	else {
		if( stp > m )
			stp = m;

		for( Index i = strt ; i < stp ; i++ )
			*(CR++) = CS2::MCFGetRC( i );
	}
}  // end( CS2::MCFGetRC( some ) )   

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber CS2::MCFGetRC( cIndex i )
{
	arc_st* arc = pos[ i ];
	return( REDUCED_COST( arc->sister->head , arc->head , arc ) / dn );

} // end( CS2::MCFGetRC( some ) )

/*--------------------------------------------------------------------------*/

void CS2::MCFGetPi( CRow P , cIndex_Set nms , cIndex strt , Index stp )
{
	node_st *nds = nodes + 1;

	if( nms ) {
		while( *nms < strt )
			nms++;

		for( Index h ; ( h = *(nms++) ) < stp ; )
			*(P++) = nds[ h ].price / dn;
	}
	else {
		if( stp > n )
			stp = n;

		for( Index i = strt ; i < stp ; i++ )
			*(P++) = nds[ i ].price / dn;
	}
}  // end( CS2::MCFGetPi( some ) )

/*--------------------------------------------------------------------------*/

CS2::FONumber CS2::MCFGetFO( void )
{ 
	if( status == kOK ) {
		if( ObjVal == Inf<FONumber>() ) {
			ObjVal = 0;
			arc_st **a = pos;
			for( Index i = m ; i-- ; a++ )
				if( ((*a)->sister)->r_cap )    
					ObjVal += FONumber( (*a)->cost ) * FONumber( ((*a)->sister)->r_cap );  

			ObjVal /=  dn;
		}

		return( ObjVal );
	}
	else
		if( status == kUnbounded )
			return( - Inf<FONumber>() );
		else
			return( Inf<FONumber>() );

}  // end( CS2::MCFGetFO )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void CS2::MCFArcs( Index_Set Startv , Index_Set Endv ,
				  cIndex_Set nms , cIndex strt , Index stp )
{
	if( stp > m )
		stp = m;

	if( nms )
		while( *nms < strt )
			nms++;

	if( Startv )
		if( Endv )
			if( nms )
				for( Index h ; ( h = *(nms++) ) < stp ; ) {
					*(Startv++) = CS2::MCFSNde( h );
					*(Endv++)   = CS2::MCFENde( h );
				}
			else
				for( Index i = strt ; i < stp ; i++ ) {
					*(Startv++) = CS2::MCFSNde( i );
					*(Endv++)   = CS2::MCFENde( i );
				}
		else
			if( nms )
				for( Index h ; ( h = *(nms++) ) < stp ; )
					*(Startv++) = CS2::MCFSNde( h );
			else
				for( Index i = strt ; i < stp ; i++ )
					*(Startv++) = CS2::MCFSNde( i );
	else
		if( Endv )
			if( nms )
				for( Index h ; ( h = *(nms++) ) < stp ; )
					*(Endv++) = CS2::MCFENde( h );
			else
				for( Index i = strt ; i < stp ; i++ )
					*(Endv++) = CS2::MCFENde( i );

}  // end( CS2::MCFArcs )

/*--------------------------------------------------------------------------*/

void CS2::MCFCosts( CRow Costv , cIndex_Set nms , cIndex strt , Index stp  )
{
	if( nms ) {
		while( *nms < strt )
			nms++;

		for( Index h ; ( h = *(nms++) ) < stp ; )
			*(Costv++) = pos[ h ]->cost;
	}
	else {
		if( stp > m )
			stp = m;

		for( Index i = strt ; i < stp ; i++ )
			*(Costv++) = pos[ i ]->cost;           
	}
}  // end( CS2::MCFCosts )

/*--------------------------------------------------------------------------*/

void CS2::MCFUCaps( FRow UCapv , cIndex_Set nms , cIndex strt , Index stp  )
{
	if( nms ) {
		while( *nms < strt )
			nms++;

		for( Index h ; ( h = *(nms++) ) < stp ; )
			*(UCapv++) = CS2::MCFUCap( h );
	}
	else {
		if( stp > m )
			stp = n;

		for( Index i = strt ; i < stp ; i++ )
			*(UCapv++) = CS2::MCFUCap( i );
	}
}  // end( CS2::MCFUCaps )

/*--------------------------------------------------------------------------*/

void CS2::MCFDfcts( FRow Dfctv , cIndex_Set nms , cIndex strt , Index stp  )
{
	if( nms ) {
		while( *nms < strt )
			nms++;

		for( Index h ; ( h = *(nms++) ) < stp ; )
			*(Dfctv++) = CS2::MCFDfct( h );
	}
	else {
		if( stp > n )
			stp = n;

		for( Index i = strt ; i < stp ; i++ )
			*(Dfctv++) = CS2::MCFDfct( i );
	}
}  // end( CS2::MCFDfcts )

/*-------------------------------------------------------------------------*/

inline MCFClass::FNumber CS2::MCFDfct( cIndex i )
{
	node_st *node = nodes + i + 1;
	FNumber Dfcti = node->excess;

	arc_st *a_stop = ( node + 1 )->closed;
	for( arc_st *a = node->first ; a < a_stop ; a++ )
		if( a->position > 0 )  // direct arc leaving from i   
			Dfcti += a->sister->r_cap;
		else                   // inverse arcs entering i
			Dfcti -= a->r_cap;

	return( - Dfcti );
}

/*-------------------------------------------------------------------------*/
/*--------- METHODS FOR ADDING / REMOVING / CHANGING DATA -----------------*/
/*-------------------------------------------------------------------------*/

void CS2::ChgCosts( cCRow NCost , cIndex_Set nms , cIndex strt , Index stp )
{
	if( nms ) {
		while( *nms < strt ) {
			nms++;
			NCost++;
		}

		for( Index h ; ( h = *(nms++) ) < stp ; )
			updtarccst( pos[ h ] , *(NCost++) );
	}
	else {
		if( stp > m )
			stp = m;

		for( Index i = strt ; i < stp ; )
			updtarccst( pos[ i++ ] , *(NCost++) );
	}

	if( ! Senstv )
		status = kUnSolved;

}  // end( CS2::ChgCosts )

/*--------------------------------------------------------------------------*/

void CS2::ChgCost( Index arc , cCNumber NCost )
{
	updtarccst( pos[ arc ] , NCost );

	if( ! Senstv )
		status = kUnSolved;

}  // end( CS2::ChgCost )

/*--------------------------------------------------------------------------*/

void CS2::ChgDfcts( cFRow NDfct ,cIndex_Set nms , cIndex strt , Index stp )
{
	if( stp > m )
		stp = m;

	if( nms ) {
		while( *nms < strt ) {
			nms++;
			NDfct++;
		}

		for( Index h ; ( h = *(nms++) ) < stp ; NDfct++) {
			cFNumber oldfct = CS2::MCFDfct( h );  // remember old deficit

			nodes[ h + 1 ].excess += oldfct - *NDfct; 

			if( ETZ( oldfct , EpsDfct ) && ( ! ETZ( *NDfct , EpsDfct ) ) )
				nodes[ h + 1 ].price = 0;
		}
	}
	else {
		if( stp > n )
			stp = n;

		for( Index i = strt ; i++ < stp ; NDfct++ ) {
			cFNumber oldfct = CS2::MCFDfct( i );  // remember old deficit

			nodes[ i + 1 ].excess += oldfct - *NDfct;

			if( ETZ( oldfct , EpsDfct ) && ( ! ETZ( *NDfct , EpsDfct ) ) )
				nodes[ i + 1 ].price = 0;   
		}
	}

	if( ! Senstv )
		status = kUnSolved;

	Blncd = false;

}  // end( CS2::ChgDfcts )

/*--------------------------------------------------------------------------*/

void CS2::ChgDfct( Index nod , cFNumber NDfct )
{  
	cFNumber oldfct = CS2::MCFDfct( nod );  // remember old deficit

	nodes[ nod + 1 ].excess += oldfct - NDfct;

	if( ETZ( oldfct , EpsDfct ) && ( ! ETZ( NDfct , EpsDfct ) ) )
		nodes[ nod + 1 ].price = 0;   

	if( ! Senstv )
		status = kUnSolved;

	Blncd = false;

}  // end( CS2::ChgDfct )

/*--------------------------------------------------------------------------*/

void CS2::ChgUCaps( cFRow NCap , cIndex_Set nms , cIndex strt , Index stp )
{
	if( nms ) {
		while( *nms < strt ) {
			nms++;
			NCap++;
		}

		for( Index h ; ( h = *(nms++) ) < stp ; )
			updtarccap( pos[ h ] , *(NCap++) );
	}
	else {
		if( stp > m )
			stp = m;

		for( Index i = strt ; i < stp ; )
			updtarccap( pos[ i++ ] , *(NCap++) );
	}

	if( ! Senstv )
		status = kUnSolved;

}  // end( CS2::ChgUCaps )

/*--------------------------------------------------------------------------*/

void CS2::ChgUCap( Index arc , cFNumber NCap )
{
	updtarccap( pos[ arc ] , NCap );

	if( ! Senstv )
		status = kUnSolved;

}  // end( CS2::ChgUCap )

/*--------------------------------------------------------------------------*/

void CS2::CloseArc( cIndex name ) 
{
#if( DYNMC_MCF_CS2 )
	arc_st *arc = pos[ name ];
	arc_st *sis = arc->sister;

	if( arc < sis->head->suspended )   // the arc is already closed
		return;                           // quietly return

	node_st *head = arc->head;
	node_st *tail = sis->head;

	// if there is flow on the arc, set the flow to 0 and signal that something
	// other than costs have changed; otherwise, pretend that the cost of the
	// closed arc has in fact gone to +INF

	if( GTZ( sis->r_cap , EpsFlw ) ) {
		INCREASE_FLOW( tail , head , arc , -(sis->r_cap) );
		Blncd = false;
	}

	// close direct arc

	EXCHANGE( arc , (tail->suspended)++ );

	// close reverse arc

	EXCHANGE( sis , (head->suspended)++ );

	if( ! Senstv )
		status = kUnSolved;
#else
	throw(
		MCFException( "CS2::CloseArc() not implemented if DYNMC_MCF_CS2 == 0"
		) );
#endif

}  // end( CS2::CloseArc )

/*--------------------------------------------------------------------------*/

void CS2::DelNode( cIndex name )
{
#if( DYNMC_MCF_CS2 )
	for( arc_st *a = nodes[ name ].suspended ,
		*a_stop = nodes[ name + 1 ].closed ; a < a_stop ; a++ )
		CloseArc( a - arcs );

	CS2::ChgDfct( name , FNumber( 0 ) );
#else
	throw(
		MCFException( "CS2::DelNode() not implemented if DYNMC_MCF_CS2 == 0"
		) );
#endif

}  // end( CS2::DelNode )

/*--------------------------------------------------------------------------*/

void CS2::OpenArc( cIndex name ) 
{
#if( DYNMC_MCF_CS2 )
	arc_st *arc = pos[ name ];
	arc_st *sis = arc->sister;

	if( arc >= sis->head->suspended )   // the arc is already open
		return;

	node_st *head = arc->head;
	node_st *tail = sis->head;

	// open direct arc

	EXCHANGE( arc , --(tail->suspended) );

	// open reverse arc

	EXCHANGE( sis , --(head->suspended) );

	if( ! Senstv )
		status = kUnSolved;

	Blncd = false;

	if( tail->price == PRICE_MIN )
		tail->price = 0;
#else
	throw(
		MCFException( "CS2::OpenArc() not implemented if DYNMC_MCF_CS2 == 0"
		) );
#endif

}  // end( CS2::OpenArc )

/*--------------------------------------------------------------------------*/

MCFClass::Index CS2::AddNode( cFNumber aDfct )
{
	throw( MCFException( "CS2::AddNode() not implemented yet" ) );

	return( 0 );
}

/*--------------------------------------------------------------------------*/

void CS2::ChangeArc( cIndex name , cIndex nSN , cIndex nEN )
{ 
	throw( MCFException( "CS2::ChangeArc() not implemented yet" ) );
}

/*--------------------------------------------------------------------------*/

void CS2::DelArc( cIndex name )
{
	CloseArc( name );  // limited implementation
}

/*--------------------------------------------------------------------------*/

MCFClass::Index CS2::AddArc( cIndex Start , cIndex End , cFNumber aU ,
							cCNumber aC )
{
	throw( MCFException( "CS2::AddArc() not implemented yet" ) );
	return( Inf<Index>() );
}

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

CS2::~CS2()
{
	if( nmax && mmax )
		MemDeAlloc();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void CS2::updtarccst( arc_st *arc , cCNumber NCa )
{
	if( NCa > m_c )
		m_c = NCa;
	else
		if( - NCa > m_c )
			m_c = - NCa;

	cCNumber SNCa = NCa * dn;

	arc->cost = SNCa;
	( arc->sister )->cost = - SNCa;
}

/*--------------------------------------------------------------------------*/

inline void CS2::updtarccap( arc_st *arc , cFNumber NCa )
{
	arc_st *sis = arc->sister;
	if( NCa < sis->r_cap ) {     // new capacity < current flow on arc
		cFNumber DCap = sis->r_cap - NCa;  // decrease flow on arc of DCap > 0

		sis->head->excess += DCap;  // update arc tail excess
		sis->r_cap = NCa;           // update residual capacity of reverse arc

		arc->head->excess -= DCap;  // update arc head excess
		arc->r_cap = 0;             // update residual capacity of forward arc

		Blncd = false;
	}
	else   
		arc->r_cap = NCa - sis->r_cap;  // update residual capacity of forward arc
}

/*--------------------------------------------------------------------------*/

bool CS2::price_update( void )
{
	// The push-relabel method modifies prices locally, one node at a time.
	// The Price Update heuristic modifies prices in a more global way. In
	// particular, a reverse "sort-of-breadth-first" visit of the residual
	// graph is performed in order to build zero rediced-cost paths between
	// all the sources (nodes with positive excess) and the sinks (nodes with
	// negative excess). The visit is started from the sinks and goes
	// backwards. In the maximum flow context, price updates, implemented using
	// breadth-first search, have been show to significantly improve practical
	// performance of the push-relabel method. In this context, the search is
	// "hybridized" trying to take into account costs: the set of active nodes
	// (those touched by the visit but not yet used to continue the search) is
	// are implemented as a set of buckets based on reduced costs of the arcs
	// used to reach them, each bucket being a queue.
	//
	// The method returns true if the price update failed because some sources
	// are unreachable: either the problem is unfeasible, or you have to return
	// suspended arcs

	Increase( n_update );

	for( node_st *i = nodes ; ++i < sentinel_node ; )
		if( LTZ( i->excess , EpsDfct ) ) {
			INSERT_TO_BUCKET( i , buckets );
			i->rank = 0;
		}
		else
			i->rank = linf;

	if( ! n_src )      // there is no source
		return( false );  // nothing to do

	Index remain = n_src;  // number of sources not touched yet

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	bucket_st *b = buckets;
	for( ; b < l_bucket ; b++ ) {
		while( NONEMPTY_BUCKET( b ) ) {
			node_st *i;
			GET_FROM_BUCKET( i , b );

			// update the rank of node "i" and, tor each arc (i, j) of the forward
			// star of "i", of "j" - - - - - - - - - - - - - - - - - - - - - - - - - -

			Increase( n_scan );

			SIndex i_rank = i->rank;
			for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
			{
				arc_st *ra = a->sister;  // an arc (j, i)

				if( GTZ( ra->r_cap , EpsFlw ) ) {
					node_st *j = a->head;
					SIndex j_rank = j->rank;

					if( j_rank > i_rank ) {
						SIndex j_new_rank;
						cCNumber rc = REDUCED_COST( j , i , ra );
						if( LTZ( rc , EpsCst ) )
							j_new_rank = i_rank;
						else {
							SIndex dr = SIndex( rc / epsilon );
							j_new_rank = ( dr < linf ? i_rank + dr + 1 : linf );
						}

						if( j_rank > j_new_rank ) {  // the rank of j has decreased
							j->rank = j_new_rank;
							j->current = ra;            // update j->current

							if( j_rank < linf )         // update position of j in the bucket
								REMOVE_FROM_BUCKET( j , buckets + j_rank );

							INSERT_TO_BUCKET( j , buckets + j_new_rank );
						}
					}
				}
			}  // end( for( all arcs from i ) )

			// update price and rank of the node i - - - - - - - - - - - - - - - - - -

			i->price -= i_rank * epsilon;
			i->rank = -1;

			// check if all sources have been reached- - - - - - - - - - - - - - - - -

			if( GTZ( i->excess , EpsDfct ) )
				if( ! (--remain) )
					break;

		}  // end( while( scanning the bucket ) )

		if( ! remain )
			break;

	}  // end( for( scanning buckets ) )

	// changing prices for nodes which were not scanned during main loop - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Important note: in certain cases, this procedure is used to early detect
	// infeasibility of the problem by discovering that some sources cannot be
	// reached from some sinks (actually, the converse). However, when flows are
	// floats, the concept of "sink" is a little bit fuzzy due to epsilons. In
	// fact, each balanced node, whose excess is comprised between - EpsDfct and
	// + EpsDfct, can as a sink, since "small" pushes of flow that do not make
	// its excess to go beyond EpsDfct does not make it active, and therefore
	// they are actually "sunk" by the node. Hence, when flows are floats the
	// fact that it is not possible to reach any sink from some source with
	// "small" (positive) excess does not imply that the problem is unfeasible,
	// since each balanced node in the same connected component can "absorb" up
	// to 2 * EpsDfct extra flow.
	//
	// Yet, implementing such a check has never shown to be useful in
	// practice, while it appears to be able to send the algorithm into an
	// infinite loop, so we don't do it.

	cCNumber dp = ( b - buckets ) * epsilon;

	/*!!
	FNumber xcss = 0;  // total excess of non-reached nodes that cannot be
	// "absorbed" by balanced nodes; note that we assume
	// all balanced node and non-reached sources to live
	// in the same connected component
	!!*/

	for( node_st *i = nodes ; ++i < sentinel_node ; )
		if( i->rank >= 0 ) {
			if( i->rank < linf )
				REMOVE_FROM_BUCKET( i , buckets + i->rank );

			if( i->price > PRICE_MIN  )
				i->price -= dp;

			/*!!
			if( i->excess <= EpsDfct )     // a balanced node can absorb
			xcss -= EpsDfct - i->excess;  // this much flow staying balanced
			else
			xcss += i->excess;
			!!*/
		}

		/*!!
		if( LEZ( xcss , EpsDfct ) )
		remain = 0;
		!!*/

		if( remain )
			return( true );
		else
			return( false );

}  // end( price_update )

/*--------------------------------------------------------------------------*/

bool CS2::relabel( node_st *i )
{
	// If i is "active" and for each arc (i, j) in the forward star of i the
	// reduced cost is >= 0, then update the price of i. If some reduced cost
	// in the forward star is negative then the method ends.

	CNumber p_max = PRICE_MIN; 
	CNumber i_price = i->price;

	// scan first half of the arcs - - - - - - - - - - - - - - - - - - - - - - -
	// scan arcs from i->current upwards

	arc_st *a_max;
	for( arc_st *a = i->current , *a_stop = (i + 1)->closed ; ++a < a_stop ; )
		if( GTZ( a->r_cap , EpsFlw ) ) {
			cCNumber dp = a->head->price - a->cost;
			if( GT( dp , p_max , EpsCst ) ) {
				if( i_price < dp ) {
					i->current = a;
					return( true );
				}

				p_max = dp;
				a_max = a;
			}
		}

		// scan second half of the arcs- - - - - - - - - - - - - - - - - - - - - - -
		// if nothing is found from i->current upwards, re-start the search from
		// i->first

		for( arc_st *a = i->first , *a_stop = i->current + 1 ; a < a_stop ; a++ )
			if( GTZ( a->r_cap , EpsFlw ) ) {
				cCNumber dp = a->head->price - a->cost;
				if( GT( dp , p_max , EpsCst ) ) {
					if( i_price < dp ) {
						i->current = a;
						return( true );
					}

					p_max = dp;
					a_max = a;
				}
			}

			// finishup- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			if( p_max > PRICE_MIN ) {            // node can be relabelled
				i->price = p_max - epsilon;
				i->current = a_max;
			}
			else                                 // node can't be relabelled
				if( i->suspended == i->first )      // and it has no supended arcs
					if( ETZ( i->excess , EpsDfct ) )  // it has zero excess
						i->price = PRICE_MIN;             // set its price to - INF
					else {
						if( n_ref == 1 )  
							status = kUnfeasible;            // the problem is unfeasible
						else     
							status = kError;                 // overflow in prices 

						return( false );        
					}
				else                                // node can't be relabelled because of
					flag_price = true;                 // supended arcs

			Increase( n_relabel );
			n_rel++;

			return( false );

}  // end( relabel ) 

/*--------------------------------------------------------------------------*/

void CS2::discharge( node_st *i )
{
	// Node "i" is active, i.e., its excess is > 0: apply push/relabel
	// operations to i until it becomes inactive.

	Increase( n_discharge );

	arc_st *a = i->current;
	node_st *j = a->head;

	if( ! ( GTZ( a->r_cap , EpsFlw ) &&
		LT( i->price + a->cost , j->price , EpsCst ) ) ) {
			relabel( i );
			if( status )
				return;

			a = i->current;
			j = a->head;
	}

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	for( long empty_push = 0 ;; ) {
		cFNumber j_exc = j->excess;
		FNumber df;
		bool bora;

		if( numeric_limits<FNumber>::is_integer )
			bora = ( j_exc >= 0 );
		else {
			// Important note: with float flows, a balanced node may in fact have a
			// a slightly negative excess (>= - EpsDfct), so funny situations may
			// occur where i has deficit larger than EpsDfct but smaller than
			// 2 * EpsDfct, j has a deficit ~= - EpsDfct, after the push both are
			// balanced; alternatively, the arc (i, j) may have such a small
			// residual capacity that pushing on j does not make j active. So, in
			// this case pushing on a balanced node without making it unbalanced has
			// to be treated as pushing on a node with negative excess.

			df = min( i->excess , a->r_cap );
			bora = ( j_exc >= - EpsDfct ) && ( j_exc + df > EpsDfct );
		}

		if( bora ) {
			// pushing on a balanced or active node- - - - - - - - - - - - - - - - - -
			// this is a "bad" situation, unless is it easy to push the flow out of j,
			// too, so check if it is the case

			arc_st *b = j->current;

			if( ( GTZ( b->r_cap , EpsFlw ) &&
				LT( j->price + b->cost , b->head->price , EpsCst ) )
				|| relabel( j ) ) {
					// "good" case: it is possible to push flow out of j

					if( status )  // problem unfeasible or error (in relabel)
						return;

					if( numeric_limits<FNumber>::is_integer )
						df = min( i->excess , a->r_cap );

					if( ETZ( j_exc , EpsDfct ) )  // it was a balanced node
						n_src++;                      // but it is no longer so

					INCREASE_FLOW( i , j , a , df );
					Increase( n_push );

					if( j->q_next == sentinel_node )
						INSERT_TO_EXCESS_Q( j );
			}
			else {
				// "bad" case: it is not possible to push flow out of j, push back
				// flow from j to i instead

				if( status )  // problem unfeasible or error (in relabel)
					return;

				arc_st *ra = a->sister;
				cFNumber rdf = min( j->excess , ra->r_cap );
				if( GTZ( rdf , EpsFlw ) ) {
					INCREASE_FLOW( j , i , ra , rdf );
					Increase( n_push );
					if( ETZ( j->excess , EpsDfct ) )
						n_src--;
				}

				if( empty_push++ >= empty_push_bound ) {
					flag_price = true;  // too many unsuccessful attempts to push flow
					return;             // out of i: time for global update 
				}
			}
		}
		else {  // pushing flow to a node with negative excess- - - - - - - - - - -
			if( numeric_limits<FNumber>::is_integer )
				df = min( i->excess , a->r_cap );

			INCREASE_FLOW( i , j , a , df );
			Increase( n_push );

			cFNumber new_j_exc = j->excess;      // the new excess of j
			if( GTZ( new_j_exc , EpsDfct ) ) {  // ... is positive
				n_src++;
				relabel( j );
				if( status )
					return;

				INSERT_TO_EXCESS_Q( j );
				total_excess += j_exc;
			}
			else
				total_excess -= df;

		}  // end( else( pushing flow to a node with negative excess ) )- - - - - -

		if( LEZ( i->excess ,  EpsDfct ) ) {
			n_src--;
			break;
		}

		if( flag_price )
			break;

		relabel( i );
		if( status )
			return;

		a = i->current;
		j = a->head;

	}  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	i->current = a;

}  // end( discharge )

/*--------------------------------------------------------------------------*/

MCFClass::Index CS2::price_in( void )
{
	// send flow through arcs with positive residual capacity and negative
	// reduced cost, update active nodes queue at the end.

	bool bad_found = false;  // true <=> we are at the second scan
	Index n_in_bad = 0;      // number of priced_in arcs with negative reduced
	// cost
restart:

	for( node_st *i = nodes ; ++i < sentinel_node ; )
		for( arc_st *a = i->first , *a_stop = i->suspended ; --a >= a_stop ; ) {
			cCNumber rc = REDUCED_COST( i , a->head , a );

			if( LTZ( rc , EpsCst ) && GTZ( a->r_cap , EpsFlw ) ) {  // bad case
				if( ! bad_found ) {
					bad_found = true;
					UPDATE_CUT_OFF();
					goto restart;
				}

				INCREASE_FLOW( i , a->head , a , a->r_cap );

				arc_st *ra = a->sister;
				node_st *j = a->head;

				arc_st *b = --( i->first );
				EXCHANGE( a , b );

				if( SUSPENDED( j , ra ) ) {
					arc_st *rb = --( j->first );
					EXCHANGE( ra , rb );
				}

				n_in_bad++;
			}
			else  // good case
				if( ( rc < cut_on ) && ( rc > -cut_on ) ) {
					arc_st *b = --(i->first);
					EXCHANGE( a , b );
				}
		}

		if( n_in_bad ) {  // recalculating excess queue
			n_bad_pricein++;
			total_excess = 0;
			n_src = 0;
			RESET_EXCESS_Q();

			for( node_st *i = nodes ; ++i < sentinel_node ; ) {
				i->current = i->first;
				cFNumber i_exc = i->excess;
				if( GTZ( i_exc , EpsDfct ) ) {  // i is a source
					n_src++;
					total_excess += i_exc;
					INSERT_TO_EXCESS_Q( i );
				}
			}

			INSERT_TO_EXCESS_Q( nodes );  // nodes[ 0 ] is a dummy node
		}

		if( time_for_price_in == TIME_FOR_PRICE_IN2 )
			time_for_price_in = TIME_FOR_PRICE_IN3;

		if( time_for_price_in == TIME_FOR_PRICE_IN1 )
			time_for_price_in = TIME_FOR_PRICE_IN2;

		return( n_in_bad );

}  // end( price_in )

/*--------------------------------------------------------------------------*/

void CS2::refine( void )
{
	// while there exist a push or a relabel operation that applies, select one
	// such operation and apply it - - - - - - - - - - - - - - - - - - - - - - -

	Increase( n_refine );
	n_ref++;
	n_rel = 0;
	n_src = 0;

	RESET_EXCESS_Q();

	time_for_price_in = TIME_FOR_PRICE_IN1;

	// initialize the queue of excess nodes- - - - - - - - - - - - - - - - - - -

	total_excess = 0;

	for( node_st *i = nodes ; ++i < sentinel_node ; ) {
		i->current = i->first;
		cFNumber i_exc = i->excess;
		if( GTZ( i_exc , EpsDfct ) ) {  // i is a source
			n_src++;
			total_excess += i_exc;
			INSERT_TO_EXCESS_Q( i );
		}
	}

	if( ! n_src )                   // the flow is feasible
		return;                        // nothing else to do

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	char pr_in_int = 0;   // current number of update between price_in

	for(;;) {
		if( ! excq_first ) {
			if( n_ref > PRICE_OUT_START )
				price_in();

			if( ! excq_first )
				break;
		}

		node_st *i;
		REMOVE_FROM_EXCESS_Q( i );

		if( GTZ( i->excess , EpsDfct ) ) {  // try to push all excess out of i
			discharge( i );
			if( status )   // problem unfeasible or error
				return;

			// number of relabel is greather than  number of nodes * UPDT_FREQ  +
			// current number of source nodes * UPDT_FREQ_S - - - - - - - - - - -

			if( ( n_rel > n * UPDT_FREQ + n_src * UPDT_FREQ_S ) || flag_price ) {
				if( GTZ( i->excess , EpsDfct ) )
					INSERT_TO_EXCESS_Q( i );

				if( flag_price && ( n_ref > PRICE_OUT_START ) ) {
					pr_in_int = 0;
					price_in();
					flag_price = false;
				}

				while( price_update() )
					if( n_ref == 1 ) {
						status = kUnfeasible;
						return;
					}
					else {
						UPDATE_CUT_OFF();
						n_bad_relabel++;
						pr_in_int = 0;
						price_in();
					}

					n_rel = 0;

					if( ( n_ref > PRICE_OUT_START ) &&
						( pr_in_int++ > time_for_price_in ) ) {
							pr_in_int = 0;
							price_in();
					}
			}  // time for update
		}
	}  // end of main loop
}  // end( refine ) 

/*--------------------------------------------------------------------------*/

bool CS2::price_refine( void )
{
	// this euristic decreases epsilon and do not change the flow f while
	// modifyng the potentials in an attempt to find a potential p such that f
	// is epsilon-optimal with respect to p

	Increase( n_prefine );

	bool cc = true;  // return value: true if flow is epsilon optimal, false if
	// refine is needed
	Index snc = 0;   // total number of negative cycle cancelled
	cIndex snc_max = ( n_ref >= START_CYCLE_CANCEL ? MAX_CYCLES_CANCELLED : 0 );

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// while negative cycle is found or eps-optimal solution is constructed

	for(;;) {
		Index nnc = 0;  // n. of negative cycles cancelled during this iteration

		for( node_st *i = nodes ; ++i < sentinel_node ; ) {
			i->rank = 0;
			i->inp = WHITE;
			i->current = i->first;
		}

		RESET_EXCESS_Q();

		for( node_st *i = nodes ; ++i < sentinel_node ; ) {
			if( i->inp == BLACK )
				continue;

			i->b_next = NULL;

			for(;;) {  // depth first search
				i->inp = GREY;

				// scanning arcs from node i starting from current

				arc_st *a = i->current , *a_stop = (i + 1)->closed;
				for( ; a < a_stop ; a++ )
					if( GTZ( a->r_cap , EpsFlw ) ) {
						node_st *j = a->head;
						if( LTZ( REDUCED_COST( i , j , a ) , EpsCst ) ) {
							if( j->inp == WHITE ) {  // fresh node - step forward
								i->current = a;
								j->b_next  = i;
								i = j;
								a = j->current;
								a_stop = (j + 1)->closed;
								break;
							}

							if( j->inp == GREY ) {  // cycle detected
								cc = false;
								nnc++;

								i->current = a;
								node_st *is = i;
								node_st *ir = i;
								FNumber df = Inf<FNumber>();

								for( ;; ir = ir->b_next ) {
									arc_st * ar = ir->current;
									if( ar->r_cap <= df ) {
										df = ar->r_cap;
										is = ir;
									}

									if( ir == j )
										break;
								}

								for( ir = i ;; ir = ir->b_next ) {
									arc_st * ar = ir->current;
									INCREASE_FLOW( ir , ar->head , ar , df );

									if( ir == j )
										break;
								}

								if( is != i ) {
									for( ir = i ; ir != is ; ir = ir->b_next )
										ir->inp = WHITE;

									i = is;
									a = (is->current) + 1;
									a_stop = (is + 1)->closed;
									break;
								}
							}

							// if j-color is BLACK, continue search from i

						}  // end( if( reduced cost negative )
					}  // end( if( residual capacity positive )

					if( a == a_stop ) {
						Increase( n_prscan1 );
						i->inp = BLACK;  // step back
						node_st *j = i->b_next;
						STACKQ_PUSH( i );

						if( j == NULL )
							break;
						i = j;
						i->current++;
					}
			}  // end( depth first search )
		}  // end( for( all nodes ) )

		// no negative cycle found: computing longest paths with eps-precision- - -

		snc += nnc;
		if( snc < snc_max )
			cc = true;

		if( ! cc )
			break;

		SIndex bmax = 0;       // number of farest nonempty bucket
		while( excq_first ) {  // scan all nodes and compute longest distances
			Increase( n_prscan2 );

			node_st *i;
			REMOVE_FROM_EXCESS_Q( i );

			SIndex i_rank = i->rank;
			for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
				if( GTZ( a->r_cap , EpsFlw ) ) {
					node_st *j = a->head;
					cCNumber rc = REDUCED_COST( i , j , a );

					if( LTZ( rc , EpsCst ) ) {  // admissible arc
						SIndex dr;
						if( numeric_limits<CNumber>::is_integer )
							dr = SIndex( ( - rc - CNumber( 1 ) ) / epsilon );
						else
							dr = SIndex( ( - rc - CNumber( 0.5 )  ) / epsilon );

						const SIndex j_rank = dr + i_rank;
						if( ( j_rank < linf ) && ( j_rank > j->rank ) )
							j->rank = j_rank;
					}
				}

				if( i_rank > 0 ) {
					if( i_rank > bmax )
						bmax = i_rank;

					INSERT_TO_BUCKET( i , buckets + i_rank );
				}
		}  // end( while )

		if( ! bmax )  // preflow is eps-optimal
			break;

		for( bucket_st *b = buckets + bmax ; b > buckets ; b-- ) {
			SIndex i_rank = b - buckets;
			cCNumber dp = CNumber( i_rank * epsilon );

			while( NONEMPTY_BUCKET( b ) ) {
				node_st *i;
				GET_FROM_BUCKET ( i, b );

				Increase( n_prscan );
				for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
					if( GTZ( a->r_cap , EpsFlw ) ) {
						node_st *j = a->head;
						SIndex j_rank = j->rank;
						if( j_rank < i_rank ) {
							SIndex j_new_rank;
							cCNumber rc = REDUCED_COST( i , j , a );
							if( LTZ( rc , EpsCst ) )
								j_new_rank = i_rank;
							else {
								SIndex dr = SIndex( rc / epsilon );
								j_new_rank = ( dr < linf ? i_rank - dr - 1 : 0 );
							}

							if( j_rank < j_new_rank ) {
								if( cc ) {
									j->rank = j_new_rank;

									if( j_rank > 0 )
										REMOVE_FROM_BUCKET( j , buckets + j_rank );

									INSERT_TO_BUCKET( j , buckets + j_new_rank );
								}
								else
									INCREASE_FLOW( i , j , a , a->r_cap );
							}
						}
					}

					i->price -= dp;

			}  // end( while( NONEMPTY_BUCKET( b ) ) )
		}  // end( for( b ) )

		if( ! cc )
			break;

	}  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// finishup: if refine needed, saturate non-epsilon-optimal arcs - - - - - -

	if( ! cc )
		for( node_st *i = nodes ; ++i < sentinel_node ; )
			for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
				if( REDUCED_COST( i , a->head , a ) < - epsilon )
					if( GTZ( a->r_cap , EpsFlw ) )
						INCREASE_FLOW( i , a->head , a , a->r_cap );

	return( cc );

}  // end( price_refine )

/*--------------------------------------------------------------------------*/

#if( COMP_DUALS )

void CS2::compute_prices( void )
{
	Increase( n_prefine );

	bool cc = true;  // true = flow is epsilon optimal, false = refine is needed

	// main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// while negative cycle is found or eps-optimal solution is constructed
	node_st *i = nodes;

	for(;;) {
		for( ; ++i < sentinel_node ; ) {
			i->rank = 0;
			i->inp = WHITE;
			i->current = i->first;
		}

		RESET_EXCESS_Q();

		for( i = nodes ; ++i < sentinel_node ; ) {
			if( i->inp == BLACK )
				continue;

			i->b_next = NULL;

			for(;;) {  // depth first search
				i->inp = GREY;

				// scanning arcs from node i
				arc_st *a = i->first , *a_stop = (i + 1)->closed;

				for(  ; a < a_stop ; a++ )
					if( GTZ( a->r_cap , EpsFlw ) ) {
						node_st *j = a->head;
						if( LTZ( REDUCED_COST( i , j , a ) , EpsCst ) ) {
							if( j->inp == WHITE ) {  // fresh node - step forward
								i->current = a;
								j->b_next  = i;
								i = j;
								a = j->current;
								a_stop = (j + 1)->closed;
								break;
							}

							if( j->inp == GREY )  // cycle detected; should not happen
								cc = false;

							// if j-color is BLACK, continue search from i
						}
					}

					if( a == a_stop ) {  // step back
						i->inp = BLACK;
						Increase( n_prscan1 );

						node_st *j = i->b_next;
						STACKQ_PUSH( i );

						if( j == NULL )
							break;

						i = j;
						i->current++;
					}
			}  // end( depth first search )
		}  // end( for( all nodes ) )

		// no negative cycle found: computing longest paths - - - - - - - - - - - -

		if( cc )
			break;

		SIndex bmax = 0;     // number of farest nonempty bucket
		while( excq_first ) {  // scan all nodes and compute longest distances
			Increase( n_prscan2 );

			REMOVE_FROM_EXCESS_Q( i );

			SIndex i_rank = i->rank;
			for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ )
				if ( GTZ( a->r_cap , EpsFlw ) ) {
					node_st *j = a->head;
					cCNumber rc = REDUCED_COST( i , j , a );

					if( LTZ( rc , EpsCst ) ) {  // admissible arc
						SIndex j_rank = SIndex( - rc ) + i_rank;
						if( j_rank < linf )
							if( j_rank > j->rank )
								j->rank = j_rank;
					}
				}

				if( i_rank > 0 ) {
					if( i_rank > bmax )
						bmax = i_rank;

					INSERT_TO_BUCKET( i , buckets + i_rank );
				}
		}  // end( while )

		if( ! bmax )
			break;

		for( bucket_st *b = buckets + bmax ; b > buckets ; b-- ) {
			SIndex i_rank = b - buckets;
			cCNumber dp = CNumber( i_rank );

			while( NONEMPTY_BUCKET( b ) ) {
				GET_FROM_BUCKET ( i, b );

				Increase( n_prscan );

				for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ;
					a++ )
					if( GTZ( a->r_cap , EpsFlw ) ) {
						node_st *j = a->head;
						SIndex j_rank = j->rank;

						if( j_rank < i_rank ) {
							SIndex j_new_rank;
							cCNumber rc = REDUCED_COST( i , j , a );

							if( LTZ( rc , EpsCst ) )
								j_new_rank = i_rank;
							else {
								SIndex dr = SIndex( rc );
								j_new_rank = ( dr < linf ? i_rank - dr - 1 : 0 );
							}

							if( j_rank < j_new_rank )
								if( cc ) {
									j->rank = j_new_rank;

									if( j_rank > 0 )
										REMOVE_FROM_BUCKET( j , buckets + j_rank );

									INSERT_TO_BUCKET( j, ( buckets + j_new_rank ) );
								}
						}
					}

					i->price -= dp;

			}  // end( while( NONEMPTY_BUCKET( b ) ) )
		}  // end( for( b ) )

		if( ! cc )
			break;

	}  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}  // end( compute_prices )

#endif

/*--------------------------------------------------------------------------*/

void CS2::price_out( void )
{
	// suspend arcs having positive reduced cost and zero flow and arcs having
	// negative reduced costs and closed.  

	cCNumber n_cut_off = - cut_off;

	for( node_st *i = nodes ; ++i < sentinel_node ; )
		for( arc_st *a = i->first , *a_stop = (i + 1)->closed ; a < a_stop ; a++ ) {
			cCNumber rc = REDUCED_COST( i , a->head , a );

			if( ( GT( rc , cut_off , EpsCst ) && LEZ( a->sister->r_cap , EpsFlw ) ) ||
				( LT( rc , n_cut_off , EpsCst ) && LEZ( a->r_cap , EpsFlw ) ) ) {
					arc_st *b = (i->first)++;  // suspend the arc
					EXCHANGE( a , b );
			}
		}

}  // end( price_out )

/*--------------------------------------------------------------------------*/

inline bool CS2::update_epsilon( void )
{
	// decrease epsilon of a value dependent from SCALE_DEFAULT after that an
	// epsilon-optimal flow is constructed

	if( epsilon <= LOW_BOUND )
		return( true );
	else {
		epsilon = ceil( epsilon / SCALE_DEFAULT );

		cut_off = CNumber( cut_off_factor ) * epsilon;
		cut_on = cut_off * CNumber( CUT_OFF_GAP );

		return( false );
	}
}  // end( update_epsilon )

/*--------------------------------------------------------------------------*/

inline MCFClass::CNumber CS2::REDUCED_COST( const node_st *i ,
										   const node_st *j ,
										   const arc_st *a )
{
	return( i->price + a->cost - j->price );
}

/*--------------------------------------------------------------------------*/

inline bool CS2::SUSPENDED( const node_st *i , const arc_st *a )
{
	return( a < i->first && a >= i->closed );
}

/*--------------------------------------------------------------------------*/

inline void CS2::EXCHANGE( arc_st *a , arc_st *b )
{
	if( a != b ) {
		arc_st *sa = a->sister;
		arc_st *sb = b->sister;

		Swap( a->head , b->head );
		Swap( a->cost , b->cost );
		Swap( a->r_cap , b->r_cap );
		Swap( a->position , b->position );

		if( a != sb ) {
			b->sister = sa;
			a->sister = sb;
			sa->sister = b;
			sb->sister = a;
		}

		if( a->position > 0 )
			pos[ a->position - 1 ] = a;

		if( b->position > 0 )
			pos[ b->position - 1 ] = b;
	}
}  // end( EXCHANGE )

/*--------------------------------------------------------------------------*/

inline void CS2::INCREASE_FLOW( node_st* i , node_st* j , arc_st* a ,
							   cFNumber df )
{
	i->excess        -= df;
	j->excess        += df;
	a->r_cap         -= df;
	a->sister->r_cap += df;
}

/*--------------------------------------------------------------------------*/

inline void CS2::UPDATE_CUT_OFF( void )
{
	if( n_bad_pricein + n_bad_relabel )
		cut_off_factor *= CUT_OFF_INCREASE;
	else {
		cut_off_factor = CUT_OFF_COEF2 * pow( dn , CUT_OFF_POWER2 );
		if( cut_off_factor < CUT_OFF_MIN )
			cut_off_factor = CUT_OFF_MIN;
	}

	cut_off = CNumber( cut_off_factor ) * epsilon;
	cut_on = cut_off * CNumber( CUT_OFF_GAP );
}

/*--------------------------------------------------------------------------*/

inline void CS2::RESET_EXCESS_Q( void )
{
	for( ; excq_first != NULL ; excq_first = excq_last ) {
		excq_last = excq_first->q_next;
		excq_first->q_next = sentinel_node;
	}
}

/*--------------------------------------------------------------------------*/

inline void CS2::INSERT_TO_EXCESS_Q( node_st* i )
{
	if( excq_first )
		excq_last->q_next = i;
	else
		excq_first = i;

	i->q_next = NULL;
	excq_last = i;
}

/*--------------------------------------------------------------------------*/

inline void CS2::INSERT_TO_FRONT_EXCESS_Q( node_st* i )
{
	if( ! excq_first )
		excq_last = i;

	i->q_next = excq_first;
	excq_first = i;
}

/*--------------------------------------------------------------------------*/

inline void CS2::REMOVE_FROM_EXCESS_Q( node_st* &i )
{
	i = excq_first;
	excq_first = i->q_next;
	i->q_next = sentinel_node;
}

/*-------------------------------------------------------------------------*/

inline void CS2::RESET_BUCKET( bucket_st *b )
{
	b->p_first = nodes;  // nodes[ 0 ] is a dummy node
}

/*-------------------------------------------------------------------------*/

inline bool CS2::NONEMPTY_BUCKET( bucket_st *b )
{
	return( b->p_first != nodes );  // an empty bucket contains only the dummy
	// node nodes[ 0 ]
}

/*--------------------------------------------------------------------------*/

inline void CS2::INSERT_TO_BUCKET( node_st* i , bucket_st* b )
{
	i->b_next = b->p_first;
	b->p_first->b_prev = i;
	b->p_first = i;
}

/*--------------------------------------------------------------------------*/

inline void CS2::GET_FROM_BUCKET( node_st *&i , bucket_st *b )
{
	i = b->p_first;
	b->p_first = i->b_next;
}

/*--------------------------------------------------------------------------*/

inline void CS2::REMOVE_FROM_BUCKET( node_st *i , bucket_st *b )
{
	if( i == b->p_first )
		b->p_first = i->b_next;
	else {
		i->b_prev->b_next = i->b_next;
		i->b_next->b_prev = i->b_prev;
	}
}

/*--------------------------------------------------------------------------*/

inline void CS2::STACKQ_PUSH( node_st* i )
{
	i->q_next = excq_first;
	excq_first = i;
}

/*--------------------------------------------------------------------------*/

void CS2::MemAlloc( void )
{
	arcs = new arc_st[ 2 * mmax ];
	nodes = new node_st[ nmax + 2 ];

	linf = Index( nmax * SCALE_DEFAULT + 2 );
	buckets = new bucket_st[ linf ];
	l_bucket = buckets + linf;

	typedef arc_st* ptr_arc_st;
	pos =  new ptr_arc_st[ mmax ]; 

}  // end( CS2::MemAlloc )

/*--------------------------------------------------------------------------*/

void CS2::MemDeAlloc( void )
{
	delete[] pos;
	delete[] buckets;
	delete[] nodes;
	delete[] arcs;

}  // end( CS2::MemDeAlloc )

/*-------------------------------------------------------------------------*/
/*---------------------- End File CS2.C -----------------------------------*/
/*-------------------------------------------------------------------------*/
