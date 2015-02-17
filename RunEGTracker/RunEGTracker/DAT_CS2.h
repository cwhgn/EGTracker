#include <fstream>
#include <sstream>

#include "CS2.h"

#define MCFSOLVER CS2

/*--------------------------------------------------------------------------*/
template<class T>
inline T ABS( const T x )
{
	return( x >= T( 0 ) ? x : -x );
}

/*--------------------------------------------------------------------------*/
// This function reads the first part of a string (before white spaces) and
// copy T value in the variable sthg (of T type)
template<class T>
static inline void str2val( const char* const str , T &sthg )
{
	istringstream( str ) >> sthg;
}

class CW_DATbyCS2//implement MAP with CS2
{
public:
	int c_maxf;
	int c_nodes;
	int c_edges;

	CW_DATbyCS2();
	void Initial(char* path,int nodes,int edges,int maxf);
	void LoadGraph();
	void Process();
	void GetResult(int* links,int* connects);

private:
	MCFClass *mcf;
	char c_path[100];

	void SetParam( MCFClass *mcf );
	void SkipComments( ifstream &iParam , string &buf );
};