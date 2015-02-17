#include "DAT_CS2.h"

CW_DATbyCS2::CW_DATbyCS2()
{
	c_maxf=0;
	c_nodes=0;
	c_edges=0;
}

void CW_DATbyCS2::Initial(char* path,int nodes,int edges,int maxf)
{
	memcpy(c_path,path,100*sizeof(char));
	c_maxf=maxf;
	c_nodes=nodes;
	c_edges=edges;
}

void CW_DATbyCS2::LoadGraph()
{
	ifstream iFile(c_path);
	if(!iFile) {
		printf("ERROR: opening input file\n");
		//return(-1);
	}
	mcf = new MCFSOLVER();
	mcf->LoadDMX(iFile);

	iFile.close();
	return;
}

void CW_DATbyCS2::Process()
{
	mcf->SetMCFTime();  // do timing

	// PECULIARITY: SetPar for parameters kEpsFlw and kEpsCst is called only if
	// FNumber/CNumber is float
	if( ! numeric_limits<MCFClass::FNumber>::is_integer ) 
	{
		MCFClass::FNumber eF = 1;
		for( MCFClass::Index i = mcf->MCFm() ; i-- ; )
			eF = max( eF , ABS( mcf->MCFUCap( i ) ) );

		for( MCFClass::Index i = mcf->MCFn() ; i-- ; )
			eF = max( eF , ABS( mcf->MCFDfct( i ) ) );   

		mcf->SetPar(MCFSOLVER::kEpsFlw, 
			(double) numeric_limits<MCFClass::FNumber>::epsilon() 
			* eF * mcf->MCFm() * 10);  // the epsilon for flows
	}

	if( ! numeric_limits<MCFClass::CNumber>::is_integer ) 
	{
		MCFClass::CNumber eC = 1;
		for( MCFClass::Index i = mcf->MCFm() ; i-- ; )
			eC = max( eC , ABS( mcf->MCFCost( i ) ) );

		mcf->SetPar(MCFSOLVER::kEpsCst, 
			(double) numeric_limits<MCFClass::CNumber>::epsilon() 
			* eC * mcf->MCFm() * 10);  // the epsilon for costs
	}

	// set other parameters from configuration file (if any)- - - - - - - - - -
	SetParam( mcf );   

	// solver call- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	mcf->SolveMCF();

	switch(mcf->MCFGetStatus()) 
	{
	case( MCFClass::kOK ):
	   cout << "Optimal Objective Function value = " << mcf->MCFGetFO() << endl;

	   double tu , ts;
	   mcf->TimeMCF( tu , ts );
	   cout << "Solution time (s): user " << tu << ", system " << ts << endl;

	   // check solution
	   //mcf->CheckPSol();
	   //mcf->CheckDSol();
	   break;
    case( MCFClass::kUnfeasible ):
	   cout << "MCF problem unfeasible." << endl;
	   break;
    case( MCFClass::kUnbounded ):
	   cout << "MCF problem unbounded." << endl;
	   break;
    default:
	   cout << "Error in the MCF solver." << endl;
	}

	if( ( ! numeric_limits<MCFClass::CNumber>::is_integer ) || ( ! numeric_limits<MCFClass::FNumber>::is_integer ) ) 
	{
		cout.setf( ios::scientific, ios::floatfield );
		cout.precision( 12 );
	}

	
}

void CW_DATbyCS2::GetResult(int* links,int* connects)
{
	MCFClass::FRow x = new MCFClass::FNumber[mcf->MCFm()];
	mcf->MCFGetX(x);
	for( MCFClass::Index i = 0 ; i < mcf->MCFm() ; i++ )
		if(x[i]!=0)
		{
			//cout << "x[" << i << "] = "  << x[ i ] << endl;
			connects[0]++;
		}
	printf("connects=%d\n",connects[0]);

	int idx=0;
	ifstream DMXs;
	DMXs.open(c_path);
	char c;
	while(1) 
	{
		if(!(DMXs>>c))break;

		int Startn,Endn;
		float LB,tU,tC;
		switch( c ) 
		{
		case( 'a' ):  // comment line- - - - - - - - - - - - - - - - - - - - - - -
			DMXs>>Startn;
			DMXs>>Endn;
			DMXs>>LB;
			DMXs>>tU;
			DMXs>>tC;
			if(idx>=mcf->MCFm())
				printf("Too many edges loaded! BUG!\n");
			if(x[idx]!=0)
			{
				if(links[Startn-1]!=0&&Startn!=1)
					printf("Too many edges from one node! BUG!\n");
				links[Startn-1]=Endn-1;
			}
			idx++;
			break;
		default:
			do DMXs.get( c );
			while(c!='\n');
			break;
		}
	}

	delete[] x;
}

////////////////////////////////////////////PRIVATE////////////////////////////////////////////////////////


// This function tries to read the parameter file; if it finds it, the
// corresponding parameters are set in the MCFClass object
void CW_DATbyCS2::SetParam( MCFClass *mcf )
{
	ifstream iParam( "config.txt" ); 
	if( ! iParam.is_open() )
		return;

	string buf;
	int num;
	SkipComments( iParam , buf );
	str2val( buf.c_str(), num );        // get number of int parameters

	for( int i = 0 ; i < num ; i++ ) {  // read all int parameters
		int param , val;

		SkipComments( iParam , buf );
		str2val( buf.c_str(), param );     // parameter name

		SkipComments( iParam , buf );
		str2val( buf.c_str(), val );       // parameter value

		mcf->SetPar( param , val );

	}  // end( for( i ) )

	SkipComments( iParam , buf );
	str2val( buf.c_str() , num );       // get number of double parameters

	for( int i = 0 ; i < num ; i++ ) {  // read all double parameters
		int param;
		double val;
		SkipComments( iParam , buf );
		str2val( buf.c_str(), param );     // parameter name

		SkipComments( iParam , buf );
		str2val( buf.c_str() , val );      // parameter value

		mcf->SetPar( param , val );

	}  // end( for( i ) )
}  // end( SetParam )

/*--------------------------------------------------------------------------*/
// This function skips comment line in a input stream, where comment line is 
// marked by an initial '#' character
void CW_DATbyCS2::SkipComments( ifstream &iParam , string &buf )
{
	do {
		iParam >> ws;
		getline( iParam , buf );
	}
	while( buf[ 0 ] == '#' );
}