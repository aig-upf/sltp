/*
Copyright (C) 2011 by the PSVN Research Group, University of Alberta
*/

#ifndef _GAMEFILE_HPP
#define _GAMEFILE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

/** Marks a variable as unused. 
    More portable than a #pragma directive. */
template <class T>
inline void UNUSED(const T&)
{
}


class PSVN {
private:
  friend ostream &operator<<( ostream &os, const PSVN &psvn );

  PSVN() {}

protected:
  vector<string> m_domainNames;
  vector<vector<string> > m_domainValues;

  vector<size_t> m_variableDomains;

  vector<vector<string> > m_ruleTests;
  vector<vector<string> > m_ruleActions;
  vector<string> m_ruleLabels;
  vector<size_t> m_ruleCosts;

  vector<vector<string> > m_goalRules;

  /* parse a file into whitespace separated uppercase tokens */
  static vector<string> getTokens( istream &os );

public:

  /* read a PSVN description from a file
     returns a newly allocated PSVN description on success
     returns NULL on failure
     NOTE: THIS WILL READ ALL INPUT UP TO END OF FILE! */
  static PSVN *fromInput( istream &is );

  size_t numDomains() const { return m_domainNames.size(); }
  const vector<string> &domainNames() const { return m_domainNames; }
  const string &domainName( const size_t d ) const {
    return m_domainNames[ d ];
  }

  /* returns the index of a named domain 'name' in list 'domainNames',
     or numDomains() if it does not exist */
  size_t getDomain( const string &name ) {
    size_t i;
    for( i = 0; i < m_domainNames.size(); ++i ) {
      if( m_domainNames[ i ].compare( name ) == 0 ) {
	break;
      }
    }
    return i;
  }

  bool domainIsUsed( const size_t domain ) const {
    for( size_t v = 0; v < m_variableDomains.size(); ++v ) {
      if( m_variableDomains[ v ] == domain ) { return true; }
    }
    return false;
  }

  const vector<vector<string> > &domainValues() const {
    return m_domainValues;
  }
  const vector<string> &domainValues( const size_t domain ) const {
    return m_domainValues[ domain ];
  }
  const string &domainValueName( const size_t domain, const size_t val ) const {
    return m_domainValues[ domain ][ val ];
  }
  size_t domainSize( const size_t domain ) const {
    return m_domainValues[ domain ].size();
  }

  /* returns value if 'name' is a value in 'domain',
     or domainSize( domain ) otherwise */
  size_t getDomainValue( const size_t domain, const string name ) const {
    size_t i;
    for( i = 0; i < m_domainValues[ domain ].size(); ++i ) {
      if( m_domainValues[ domain ][ i ].compare( name ) == 0 ) {
	break;
      }
    }
    return i;
  }

  /* creates and returns the index of an integer domain
     returns existing index if it already exists */
  size_t addIntegerDomain( const size_t size, const bool oneBased );

  const size_t numVariables() const { return m_variableDomains.size(); }
  const vector<size_t> &variableDomains() const { return m_variableDomains; }
  const size_t variableDomain( const size_t var ) const {
    return m_variableDomains[ var ];
  }

  const size_t numRules() const { return m_ruleTests.size(); }
  const vector<vector<string> > &ruleTests() const { return m_ruleTests; }
  const vector<string> &ruleTests( const size_t rule ) const {
    return m_ruleTests[ rule ];
  }
  const vector<vector<string> > &ruleActions() const { return m_ruleActions; }
  const vector<string> &ruleActions( const size_t rule ) const {
    return m_ruleActions[ rule ];
  }
  const vector<string> &ruleLabels() const { return m_ruleLabels; }
  const string ruleLabel( const size_t rule ) const {
    return m_ruleLabels[ rule ];
  }
  const size_t ruleCost( const size_t rule) const {
    return m_ruleCosts[ rule ];
  }

  const size_t numGoalRules() const { return m_goalRules.size(); }
  const vector<vector<string> > &goalRules() const { return m_goalRules; }
  const vector<string> &goalRule( const size_t goal ) const {
    return m_goalRules[ goal ];
  }


  /* Apply a domain abstraction to the game

     mapping[ domain ][ value ] gives the new value' that value will
     be mapped to.  mapping[ domain ] may be empty if the entire
     domain is un-touched.  Any values in mapping[domain][] which are
     out of range for domain will be ignored.

     NOTE: applying multiple domain abstractions produces well
     defined, but possibly unexpected results.  For example, a single
     mapping of {A->B,B->C,C->C} will likely produce different results
     than {A->B,B->B,C->C} followed by {A->A,B->C,C->C}.  If there is
     a goal A B C, the first mapping produces the goal B C C, while
     the second pair of mappings produces the goal C C C */
  void domainAbstraction( const vector<vector<size_t> > &mapping );

  /* Apply a projection to the game

     projection[ var ] is true if the variable is removed, and
     false if the projection keeps the variable */
  void projectionAbstraction( const vector<bool> &projection );
};

#endif
