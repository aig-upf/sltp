/*
Copyright (C) 2011 by the PSVN Research Group, University of Alberta
*/

#include "psvn.hpp"
#include <assert.h>
#include <sstream>

ostream &operator<<( ostream &os, const PSVN &psvn )
{
  size_t i, j;

  /* find un-used domains */
  vector<bool> domainUsed( psvn.numVariables(), false );
  for( i = 0; i < psvn.numVariables(); ++i ) {
    domainUsed[ psvn.variableDomain( i ) ] = true;
  }

  bool explicitDomains = false;
  for( i = 0; i < psvn.numDomains(); ++i ) {
    if( !domainUsed[ i ] ) {
      continue;
    }
    if( sscanf( psvn.domainName( i ).c_str(), "%zu", &j ) ) {
      /* don't print out integer domains */
      continue;
    }

    explicitDomains = true;
    os << "DOMAIN " << psvn.domainName( i )
	 << " " << psvn.domainValues( i ).size() << endl;
    os << " ";

    for( j = 0; j < psvn.domainValues( i ).size(); ++j ) {
      os << " " << psvn.domainValues( i )[ j ];
    }
    os << endl;
  }
  if( explicitDomains ) {
    os << endl;
  }

  os << psvn.numVariables() << endl;
  for( i = 0; i < psvn.numVariables(); ++i ) {
    if( i ) {
      os << " ";
    }
    os << psvn.domainName( psvn.variableDomain( i ) );
  }
  os << endl << endl;

  for( i = 0; i < psvn.numRules(); ++i ) {

    for( j = 0; j < psvn.numVariables(); ++j ) {
      if( j ) {
	os << " ";
      }
      os << psvn.ruleTests( i )[ j ];
    }

    os << " =>";

    for( j = 0; j < psvn.numVariables(); ++j ) {
      os << " " << psvn.ruleActions( i )[ j ];
    }

    os << " LABEL " << psvn.ruleLabel( i );
    os << " COST " << psvn.ruleCost( i ) << endl;
  }
  os << endl;

  for( i = 0; i < psvn.numGoalRules(); ++i ) {

    os << "GOAL";
    for( j = 0; j < psvn.numVariables(); ++j ) {
      os << " " << psvn.goalRule( i )[ j ];
    }
    os << endl;
  }

  return os;
}

vector<string> PSVN::getTokens( istream &is )
{
  int c;
  vector<string> tokens;
  string tokenString;

  while( ( c = is.get() ) >= 0 ) {

    if( isspace( c ) ) {

      if( tokenString.length() > 0 ) {
	tokens.push_back( tokenString );
	tokenString.clear();
      }
    } else {

      if( tokenString.size() == 0 && ( c == '#' || c == ';' ) ) {
	/* discard rest of line after a leading # or ; */

	while( ( c = is.get() ) >= 0 && c != '\n' ) {}
      } else {

	tokenString.push_back( toupper( c ) );
      }
    }
  }
  if( tokenString.length() > 0 ) {
    tokens.push_back( tokenString );
  }

  return tokens;
}

static bool nameIsIntegerDomain( const string name,
				 size_t *size, bool *oneBased )
{
  size_t pos;

  pos = name.find_first_not_of( "0123456789" );
  if( pos == name.npos ) {
    assert( sscanf( name.c_str(), "%zu", size ) >= 1 );
    *oneBased = false;
    return true;
  }

  if( pos > 0 && pos == name.size() - 1 && name[ pos ] == 'N' ) {
    assert( sscanf( name.c_str(), "%zu", size ) >= 1 );
    *oneBased = true;
    return true;
  }

  return false;
}

size_t PSVN::addIntegerDomain( const size_t size, const bool oneBased )
{
  char t[ 64 ];
  string name;
  size_t ret, i;

  sprintf( t, "%zu", size );
  name.append( t );
  if( oneBased ) {
    name.append( "N" );
  }

  /* check if we've already constructed the domain */
  ret = getDomain( name );
  if( ret < numDomains() ) {
    return ret;
  }

  /* domain not yet created - make it */
  ret = m_domainNames.size();
  m_domainNames.push_back( name );
  m_domainValues.push_back( vector<string>() );
  for( i = 0; i < size; ++i ) {
    sprintf( t, "%zu", oneBased ? i + 1 : i );
    m_domainValues[ ret ].push_back( string( t ) );
  }

  return ret;
}

PSVN *PSVN::fromInput( istream &is )
{
  vector<string> tokens = PSVN::getTokens( is );
  PSVN *psvn = new PSVN;
  size_t tokenIdx, i;
  bool b;

  tokenIdx = 0;

  /* get any named domains */
  size_t domainSize;
  while( tokens[ tokenIdx ].compare( "DOMAIN" ) == 0 ) {
    ++tokenIdx;

    if( tokenIdx + 2 > tokens.size() ) {
      cerr << "must specify a name and size for domain" << endl;
      delete psvn;
      return NULL;
    }

    /* get the name and size */
    if( nameIsIntegerDomain( tokens[ tokenIdx ], &i, &b ) ) {
      cerr << tokens[ tokenIdx ] << " is a reserved domain name" << endl;
      delete psvn;
      return NULL;
    }
    if( sscanf( tokens[ tokenIdx + 1 ].c_str(), "%zu", &domainSize ) < 1
	|| domainSize < 2 ) {
      cerr << "bad domain size " << domainSize << " for "
	   << tokens[ tokenIdx ] << endl;
      delete psvn;
      return NULL;
    }
    psvn->m_domainNames.push_back( tokens[ tokenIdx ] );
    psvn->m_domainValues.push_back( vector<string>() );
    tokenIdx += 2;

    /* get the value names */
    if( tokenIdx + domainSize > tokens.size() ) {
      cerr << "must specify " << domainSize << " value names for domain "
	   << psvn->m_domainNames.back() << endl;
      delete psvn;
      return NULL;
    }
    for( i = 0; i < domainSize; ++i ) {

      /* avoid duplicates */
      for( size_t j = 0; j < psvn->m_domainValues.back().size(); ++j ) {
	if( psvn->m_domainValues.back()[ j ] == tokens[ tokenIdx + i ] ) {
	  cerr << "duplicate domain value " << tokens[ tokenIdx + i ]
	       << " in domain " << psvn->m_domainNames.back() << endl;
	  delete psvn;
	  return NULL;
	}
      }

      psvn->m_domainValues.back().push_back( tokens[ tokenIdx + i ] );
    }
    tokenIdx += domainSize;
  }

  /* get number of variables */
  size_t numVars;
  if( tokenIdx == tokens.size() ) {
    cerr << "missing number of variables" << endl;
    delete psvn;
    return NULL;
  }
  if( sscanf( tokens[ tokenIdx ].c_str(), "%zu", &numVars ) < 1
      || numVars < 1 ) {
    cerr << "bad number of variables: " << tokens[ tokenIdx ] << endl;
    delete psvn;
    return NULL;
  }
  ++tokenIdx;

  /* get domains for each variable */
  if( tokenIdx + numVars > tokens.size() ) {
    cerr << "insufficient domain names for " << numVars << " variables" << endl;
    delete psvn;
    return NULL;
  }
  do {

    if( nameIsIntegerDomain( tokens[ tokenIdx ], &i, &b ) ) {
      i = psvn->addIntegerDomain( i, b );
    } else {
      i = psvn->getDomain( tokens[ tokenIdx ] );
      if( i >= psvn->numDomains() ) {
	cerr << tokens[ tokenIdx ]
	     << " is an undefined variable domain" << endl;
	delete psvn;
	return NULL;
      }
    }

    psvn->m_variableDomains.push_back( i );

    ++tokenIdx;
  } while( psvn->m_variableDomains.size() < numVars );

  /* get the rules */
  size_t numRulesRead = 0;
  while( tokenIdx < tokens.size() ) {

    if( tokens[ tokenIdx ].compare( "GOAL" ) == 0 ) {
      /* goal rule */
      ++tokenIdx;

      if( tokenIdx + numVars > tokens.size() ) {
	cerr << "too few tokens for a GOAL rule" << endl;
	delete psvn;
	return NULL;
      }

      psvn->m_goalRules.push_back( vector<string>() );
      for( i = 0; i < numVars; ++i ) {
	psvn->m_goalRules.back().push_back( tokens[ tokenIdx + i ] );
      }
      tokenIdx += numVars;
    } else {
      /* normal rule */
      ++numRulesRead;

      /* check for sufficient tokens and test/action separator */
      if( tokenIdx + numVars * 2 + 1 > tokens.size() ) {
	cerr << "too few tokens for a normal rule" << endl;
	delete psvn;
	return NULL;
      }
      if( tokens[ tokenIdx + numVars ].compare( "=>" ) != 0 ) {
	cerr << "could not parse rule - found "
	     << tokens[ tokenIdx + numVars ] << " instead of =>" << endl;
	delete psvn;
	return NULL;
      }

      /* get the tests */
      psvn->m_ruleTests.push_back( vector<string>() );
      for( i = 0; i < numVars; ++i ) {
	psvn->m_ruleTests.back().push_back( tokens[ tokenIdx + i ] );
      }
      tokenIdx += numVars + 1;

      /* get the actions */
      psvn->m_ruleActions.push_back( vector<string>() );
      for( i = 0; i < numVars; ++i ) {
	psvn->m_ruleActions.back().push_back( tokens[ tokenIdx + i ] );
      }
      tokenIdx += numVars;

      /* get a name for the rule */
      if( tokenIdx + 2 <= tokens.size()   /** ROB CHANGE < to <= **/
	  && tokens[ tokenIdx ] == "LABEL" ) {
	psvn->m_ruleLabels.push_back( tokens[ tokenIdx + 1 ] );
	tokenIdx += 2;
      } else {
	char cStr[ 64 ];
	sprintf( cStr, "rule_%zu", numRulesRead );
	psvn->m_ruleLabels.push_back( string( cStr ) );
      }

      /* get a cost for the rule */
      if ( tokenIdx + 2 <= tokens.size()   /** ROB CHANGE < to <= **/
           && tokens[ tokenIdx ] == "COST" ) {
          int cost;
          istringstream is( tokens[ tokenIdx + 1 ] );
          if (!(is >> cost)) {
              cerr << "could not parse cost from '" << tokens[ tokenIdx + 1 ]
                   << "'." << endl;
              delete psvn;
              return NULL;
          }
          if (cost < 0) {
              cerr << "invalid cost (" << cost << ")" << endl;
              delete psvn;
              return NULL;
          }
          psvn->m_ruleCosts.push_back( cost );
          tokenIdx += 2;
      } else {
          psvn->m_ruleCosts.push_back( 1 );
      }
    }
  }

  /* deal with any leftover input */
  if( tokenIdx + 1 <= tokens.size() ) {
    cerr << "ignoring input:";
    while( tokenIdx < tokens.size() ) {
      cerr << " " << tokens[ tokenIdx ];
      ++tokenIdx;
    }
    cerr << endl;
  }

  return psvn;
}


void PSVN::domainAbstraction( const vector<vector<size_t> > &mapping )
{
  size_t var, domain, i, val;

  assert( mapping.size() == numDomains() );

  for( var = 0; var < numVariables(); ++var ) {
    domain = variableDomain( var );
    if( mapping[ domain ].size() == 0 ) {
      continue;
    }

    /* map rule tests */
    for( i = 0; i < m_ruleTests.size(); ++i ) {
      val = getDomainValue( domain, m_ruleTests[ i ][ var ] );
      if( val < domainSize( domain ) && mapping[ domain ][ val ] != val ) {
	m_ruleTests[ i ][ var ]
	  = domainValueName( domain, mapping[ domain ][ val ] );
      }
    }

    /* map rule actions */
    for( i = 0; i < m_ruleActions.size(); ++i ) {
      val = getDomainValue( domain, m_ruleActions[ i ][ var ] );
      if( val < domainSize( domain ) && mapping[ domain ][ val ] != val ) {
	m_ruleActions[ i ][ var ]
	  = domainValueName( domain, mapping[ domain ][ val ] );
      }
    }

    /* map goal tests */
    for( i = 0; i < m_goalRules.size(); ++i ) {
      val = getDomainValue( domain, m_goalRules[ i ][ var ] );
      if( val < domainSize( domain ) && mapping[ domain ][ val ] != val ) {
	m_goalRules[ i ][ var ]
	  = domainValueName( domain, mapping[ domain ][ val ] );
      }
    }
  }
}

void PSVN::projectionAbstraction( const vector<bool> &projection )
{
  size_t var, domain, i;

  assert( projection.size() == numVariables() );

  for( var = 0; var < numVariables(); ++var ) {
    if( !projection[ var ] ) {
      continue;
    }

    domain = addIntegerDomain( 1, false );
    m_variableDomains[ var ] = domain;

    /* remove variable from tests */
    for( i = 0; i < m_ruleTests.size(); ++i ) {
      m_ruleTests[ i ][ var ] = string( "-" );
    }

    /* remove variable from actions */
    for( i = 0; i < m_ruleActions.size(); ++i ) {
      m_ruleActions[ i ][ var ] = string( "-" );
    }

    /* remove variable from goals */
    for( i = 0; i < m_goalRules.size(); ++i ) {
      m_goalRules[ i ][ var ] = string( "-" );
    }
  }
}
