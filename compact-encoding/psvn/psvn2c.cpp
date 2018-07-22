/*
Copyright (C) 2011, 2012, 2014 by the PSVN Research Group, University of Alberta
*/

#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "psvn2c.hpp"

/*
The original (2011) version of the move pruning code tested two things:
(a) is the newly generated operator sequence redundant with a previously generated sequence?, and
(b) is the some previously generated sequence redundant with the newly generated operator sequence?
In 2012 we discovered that it was unsafe to do both of these, so the current version only does (a).
This is thoroughly explained in Burch and Holte's SoCS 2012 paper, and the correctness of (a) is formally proven.

I have left the original 2011 code here in case anyone wants to reproduce those results and/or
test other ways of making move pruning safe.  It has been disabled by putting it inside #ifdef's.
To reinstate it, just move the following #define outside this comment block.
#define UnsafeVersion
*/

/* Verbosity Levels */
#define V_QUIET 0
#define V_USER  1
#define V_DEBUG0 2
#define V_DEBUG1 3
#define V_DEBUG2 4
#define V_DEBUG_MP 5


/* DEFINE THE DEFAULT SETTINGS AND MAXIMUM VALUES FOR THE COMMAND LINE PARAMETERS */
#define  historyLenMAX 99
#define  fwdHistoryLenDEFAULT  1
#define  bwdHistoryLenDEFAULT  1
#define  bwdMovesDEFAULT       true
#define  ruleGroupSizeDEFAULT  100
#define  ruleGroupSizeMAX  1000
#define  removeDuplicateRulesDEFAULT  false
#define  dynamicAbstractionsDEFAULT  false
#define  AbstractionDEFAULT  true
#define  state_mapDEFAULT  true
#define  verbosityDEFAULT V_QUIET
#define  verbosityMAX 9


#define ASSERT( X )

static bool valueInSet( const Value &val, const vector<Value> &values )
{
  size_t i;

  for( i = 0; i < values.size(); ++i ) {
    if( val == values[ i ] ) {
      return true;
    }
  }

  return false;
}

static bool valuesAreSubsetOf( const vector<Value> values,
			       const vector<Value> set )
{
  size_t i;

  for( i = 0; i < values.size(); ++i ) {

    if( !valueInSet( values[ i ], set ) ) {
      return false;
    }
  }

  return true;
}

bool Rule::operator==( const Rule &rule ) const {
  size_t i;

  if( m_id != rule.m_id ) {
    return false;
  }

  if ( m_cost != rule.m_cost ) {
    return false;
  }

  for( i = 0; i < tests().size(); ++i ) {
    if( !valuesAreSubsetOf( tests()[ i ], rule.tests()[ i ] )
	|| !valuesAreSubsetOf( rule.tests()[ i ], tests()[ i ] ) ) {
      return false;
    }
  }

  for( i = 0; i < actions().size(); ++i ) {
    if( !valuesAreSubsetOf( actions()[ i ], rule.actions()[ i ] )
	|| !valuesAreSubsetOf( rule.actions()[ i ], actions()[ i ] ) ) {
      return false;
    }
  }

  return true;
}

bool Rule::isEmpty() const
{
    for (size_t i = 0; i < tests().size(); ++i)
        if (!tests()[i].empty()) 
            return false;
    for (size_t i = 0; i < actions().size(); ++i)
        if (!actions()[i].empty()) 
            return false;
    return true;
}

void Rule::print() const
{
  cerr << toString() << endl;
}

void Rule::useEquivalence( size_t x, size_t y )
{
  size_t i, j;

  /* remove x tests and add them to y */
  for( i = 0; i < m_tests[ x ].size(); ++i ) {

    /* skip x == x? and x == y? test */
    if( m_tests[ x ][ i ].type() == Label ) {
      if( m_tests[ x ][ i ].value() == x || m_tests[ x ][ i ].value() == y ) {
	continue;
      }
    }

    /* make sure test is not a duplicate */
    for( j = 0; j < m_tests[ y ].size(); ++j ) {
      if( m_tests[ y ][ j ] == m_tests[ x ][ i ] ) {
	goto useEquivalenceDoneWithTest;
      }
    }

    /* move test to variable y */
    m_tests[ y ].push_back( m_tests[ x ][ i ] );

  useEquivalenceDoneWithTest:;
  }
  m_tests[ x ].clear();

  /* go through all tests and replace references to x with y */
  for( i = 0; i < m_tests.size(); ++i ) {
    for( j = 0; j < m_tests[ i ].size(); ++j ) {

      if( m_tests[ i ][ j ].type() == Label
	  && m_tests[ i ][ j ].value() == x ) {
	m_tests[ i ][ j ] = Value( Label, y );
      }
    }
  }

  /* actions don't change - nothing to do */
}

bool Rule::useKnowledge( const Context &context )
{
  size_t i, j, k;
  bool proven;

  for( i = 0; i < m_tests.size(); ++i ) {
    for( j = 0; j < m_tests[ i ].size(); ++j ) {

      if( context.testIsProvable( i, m_tests[ i ][ j ], &proven ) ) {
	if( proven ) {
	  /* proven - remove the test */

	  for( k = j + 1; k < m_tests[ i ].size(); ++k ) {
	    m_tests[ i ][ k - 1 ] = m_tests[ i ][ k ];
	  }
	  m_tests[ i ].pop_back();
	  --j;

	  m_modified = true;
	} else {
	  /* falsified - don't do any more work */

	  return false;
	}
      }
    }
  }

  return true;
}

Rule::Rule( const Rule &rule, const vector<pair<size_t,size_t> > &unboundVars,
	    const vector<size_t> &unboundVals )
  : m_id( rule.id() ), m_name( rule.name() ), m_cost( rule.cost() ),
    m_tests( rule.tests() ), m_actions( rule.actions() ),
    m_modified( rule.m_modified )
{
  size_t i, j;

  for( i = 0; i < m_actions.size(); ++i ) {

    if( m_actions[ i ].size() && m_actions[ i ][ 0 ].type() == Label
	&& m_actions[ i ][ 0 ].value() >= m_actions.size() ) {
      /* variable i in the action is an unbound variable */

      for( j = 0; j < unboundVars.size(); ++j ) {

	if( unboundVars[ j ].first == m_actions[ i ][ 0 ].value() ) {
	  break;
	}
      }

      if( j < unboundVars.size() ) {
	/* we have a replacement value */

	m_actions[ i ][ 0 ] = Value( Constant, unboundVals[ j ] );
      }
    }
  }
}

bool Rule::usesVariable( const size_t var ) const
{
  size_t i, j;

  if( tests()[ var ].size() ) {
    return true;
  }

  for( i = 0; i < tests().size(); ++i ) {
    /* we know tests()[ var ].size() == 0, so we don't need to test i == var */

    for( j = 0; j < tests()[ i ].size(); ++j ) {
      if( tests()[ i ][ j ].type() == Label
	  && tests()[ i ][ j ].value() == var ) {
	return true;
      }
    }
  }

  return false;
}

void Rule::generateActionCode( ostream &os,
			       const string &oldState,
			       const string &newState,
			       const size_t indentation ) const
{
  size_t i;

  for( i = 0; i < actions().size(); ++i ) {

    for( size_t j = 0; j < indentation; ++j ) {
      os << " ";
    }
    os << newState << "->vars[ " << i << " ] = ";

    if( actions()[ i ].size() == 0 ) {

      os << oldState << "->vars[ " << i << " ];\n";
    } else {

      ASSERT( actions()[ i ].size() == 1 );

      if( actions()[ i ][ 0 ].type() == Constant ) {

	os << actions()[ i ][ 0 ].value() << ";\n";
      } else {
	ASSERT( actions()[ i ][ 0 ].type() == Label );
	if ( actions()[ i ][ 0 ].value() < actions().size() ) {
	   os << oldState << "->vars[ " << actions()[ i ][ 0 ].value() << " ];\n";
        } else {
	   os << newState << "->vars[ " << (actions()[ i ][ 0 ].value() - actions().size()) << " ];\n";
        }
      }
    }
  }
}

void Rule::generateActionFunction( ostream &os,
				   const string &treeName,
                                   const vector<size_t>& varDomains,
                                   const bool dynamicAbstractions ) const
{
  os << "\nstatic void " << actionFunctionName( treeName )
     << "( const state_t *state, state_t *child_state )\n{\n";
  generateActionCode( os, string( "state" ), string( "child_state" ), 2 );
  os << "}\n";
}

static void replaceVariable( vector<Value> &values, const size_t label,
			     const Value &val )
{
  size_t i;

  for( i = 0; i < values.size(); ++i ) {

    if( values[ i ].type() == Label
	&& values[ i ].value() == label ) {
      /* value i is the label we're replacing */

      values[ i ] = val;
    }
  }
}

bool Rule::applyTo( vector<Value> &values,
		    vector<Value> &context ) const
{
  size_t i, j, label;

  /* check if the rule can be applied to the state */
  ASSERT( values.size() == tests().size() );
  for( i = 0; i < tests().size(); ++i ) {

    for( j = 0; j < tests()[ i ].size(); ++j ) {

      if( values[ i ].type() == Label ) {
	/* value is currently some unknown */

	if( tests()[ i ][ j ].type() == Constant ) {
	  /* rule specifies a constant - replace with that value */

	  replaceVariable( context, values[ i ].value(), tests()[ i ][ j ] );
	  replaceVariable( values, values[ i ].value(), tests()[ i ][ j ] );
	} else {
	  /* rule specifies a label */
	  ASSERT( tests()[ i ][ j ].type() == Label );

	  label = tests()[ i ][ j ].value();
	  if( !( values[ i ] == values[ label ] ) ) {
	    /* value the label points to is different from value[ i ] */

	    if( values[ label ].type() == Constant ) {

	      replaceVariable( context, values[ i ].value(),
			       values[ label ] );
	      replaceVariable( values, values[ i ].value(),
			       values[ label ] );
	    } else {
	      ASSERT( values[ label ].type() == Label );

	      if( values[ label ].value() > values[ i ].value() ) {

		replaceVariable( context, values[ label ].value(),
				 values[ i ] );
		replaceVariable( values, values[ label ].value(),
				 values[ i ] );
	      } else {

		replaceVariable( context, values[ i ].value(),
				 values[ label ] );
		replaceVariable( values, values[ i ].value(),
				 values[ label ] );
	      }
	    }
	  }
	}
      } else {
	/* value is currently a constant */
	ASSERT( values[ i ].type() == Constant );

	if( tests()[ i ][ j ].type() == Constant ) {
	  /* rule specifies a constant - can't be used if it doesn't match */

	  if( tests()[ i ][ j ].value() != values[ i ].value() ) {

	    return false;
	  }
	} else {
	  /* rule specifies a label */
	  ASSERT( tests()[ i ][ j ].type() == Label );

	  label = tests()[ i ][ j ].value();
	  if( values[ label ].type() == Constant ) {
	    /* referenced value is constant */

	    if( tests()[ i ][ j ].value() != values[ i ].value() ) {

	      return false;
	    }
	  } else {
	    /* referenced value is a label - replace it with values[ i ] */
	    ASSERT( values[ label ].type() == Label );

	    replaceVariable( context, values[ label ].value(), values[ i ] );
	    replaceVariable( values, values[ label ].value(), values[ i ] );
	  }
	}
      }
    }
  }

  /* rule applies - apply it! */
  vector<Value> tValues = values;
  ASSERT( tValues.size() == actions().size() );
  for( i = 0; i < actions().size(); ++i ) {
    if( actions()[ i ].size() == 0 ) {
      continue;
    }
    ASSERT( actions()[ i ].size() == 1 );

    if( actions()[ i ][ 0 ].type() == Constant ) {

      values[ i ] = actions()[ i ][ 0 ];
    } else {
      ASSERT( actions()[ i ][ 0 ].type() == Label );

      values[ i ] = tValues[ actions()[ i ][ 0 ].value() ];
    }
  }

  return true;
}


Context::Context( const Game &game )
  : m_game( game ) {
  m_proven.resize( game.numVars(), Value( Label, 0 ) );
  m_falsified.resize( game.numVars() );
}

Context::Context( const Context &context,
		  const VarTest &test, const size_t result )
  : m_game( context.m_game )
{
  m_proven = context.proven();
  m_falsified = context.falsified();

  if( test.type() == Label ) {

    if( result == 0 ) {

      setTestFalsified( test.var(), Value( Label, test.otherLabel() ) );
    } else {

      setTestProven( test.var(), Value( Label, test.otherLabel() ) );
    }
  } else {

    setTestProven( test.var(), Value( test.type(), result ) );
  }

}

bool Context::operator==( const Context &context ) const {
  size_t i;

  for( i = 0; i < m_game.numVars(); ++i ) {

    if( !( proven()[ i ] == context.proven()[ i ] ) ) {
      return false;
    }

    if( !valuesAreSubsetOf( falsified()[ i ], context.falsified()[ i ] )
	|| !valuesAreSubsetOf( context.falsified()[ i ], falsified()[ i ] ) ) {
      return false;
    }
  }

  return true;
}

void Context::removeDuplicates()
{
  size_t i, j, k;

  for( i = 0; i < m_game.numVars(); ++i ) {
    for( j = 1; j < m_falsified[ i ].size(); ++j ) {
      for( k = 0; k < j; ++k ) {
	if( m_falsified[ i ][ j ] == m_falsified[ i ][ k ] ) {
	  /* m_falsified[ i ][ j ] is a duplicate - remove it */

	  for( k = j + 1; k < m_falsified[ i ].size(); ++k ) {
	    m_falsified[ i ][ k - 1 ] = m_falsified[ i ][ k ];
	  }
	  m_falsified[ i ].pop_back();
	  --j;
	}
      }
    }
  }
}

bool Context::testIsProvable( const size_t var, const Value &test,
			      bool *ret_isProven ) const
{
  size_t i;

  if( test.type() == Constant ) {
    /* var == [0,N-1] ? */

    if( m_proven[ var ].type() == Constant ) {
      /* var = value  -  is value test.value()? */

      *ret_isProven = m_proven[ var ].value() == test.value();
      return true;
    }

    for( i = 0; i < m_falsified[ var ].size(); ++i ) {
      if( m_falsified[ var ][ i ].type() == Constant ) {
	/* var != value  -  is value test.value()? */

	if( m_falsified[ var ][ i ].value() == test.value() ) {
	  *ret_isProven = false;
	  return true;
	}
      }
    }
  } else {
    /* var == y ? */
    ASSERT( test.type() == Label );

    /* catch var == var */
    if( test.value() == var ) {
      *ret_isProven = true;
      return true;
    }

    if( m_proven[ var ].type() == Constant ) {
      /* var = value  -  is y == value? */

      if( testIsProvable( test.value(), m_proven[ var ], ret_isProven ) ) {
	return true;
      }
    }

    for( i = 0; i < m_falsified[ var ].size(); ++i ) {
      if( m_falsified[ var ][ i ].type() == Constant ) {
	/* var != value  -  is y == value? */

	bool b;
	if( testIsProvable( test.value(), m_falsified[ var ][ i ], &b ) && b ) {
	  *ret_isProven = false;
	  return true;
	}
      } else {
	/* var != label  -  is label y? */
	ASSERT( m_falsified[ var ][ i ].type() == Label );

	if( m_falsified[ var ][ i ].value() == test.value() ) {
	  *ret_isProven = false;
	  return true;
	}
      }
    }
  }

  return false;
}

void Context::setTestProven( const size_t var, const Value &test )
{
  size_t i;

  if( test.type() == Constant ) {
    /* var = [0,N-1] */

    /* we shouldn't have already proven var = something else */
    ASSERT( m_proven[ var ].type() == Label /* nothing proven */
	    || m_proven[ var ] == test );

    m_proven[ var ] = test;

    /* no longer need var != [0,N-1] for any value */
    vector<Value> oldFalsified = m_falsified[ var ];
    m_falsified[ var ].clear();
    for( i = 0; i < oldFalsified.size(); ++i ) {
      if( oldFalsified[ i ].type() != Constant ) {
	ASSERT( oldFalsified[ i ].type() == Label );

	m_falsified[ var ].push_back( oldFalsified[ i ] );
      }
    }

    /* for every y such that x != y, add y != test.value,
       unless we already know y = [0,N-1]
       NOTE: may add duplicates because we don't check if it's already there */
    for( i = 0; i < m_falsified[ var ].size(); ++i ) {

      if( m_falsified[ var ][ i ].type() == Label
	  && m_proven[ m_falsified[ var ][ i].value() ].type() != Constant ) {
	m_falsified[ m_falsified[ var ][ i ].value() ].push_back( test );
      }
    }
  } else {
    /* var = y */
    ASSERT( test.type() == Label );
    size_t y = test.value();

    /* rewrite all facts to get rid of variable */
    m_proven[ var ] = Value( Label, 0 );

    for( i = 0; i < m_falsified[ var ].size(); ++i ) {
      m_falsified[ y ].push_back( m_falsified[ var ][ i ] );
    }
    m_falsified[ var ].clear();
  }
}

void Context::setTestFalsified( const size_t var, const Value &test )
{
  if( test.type() == Constant ) {
    /* var != [0,N-1] */

    m_falsified[ var ].push_back( test );
  } else {
    /* var != y */
    ASSERT( test.type() == Label );

    m_falsified[ var ].push_back( test );

    if( m_proven[ test.value() ].type() == Constant ) {
      /* if we know y = (C=[0,N-1]), add x != C */

      m_falsified[ var ].push_back( m_proven[ test.value() ] );
    } else {
      /* y is not a constant, so add y != var */

      m_falsified[ test.value() ].push_back( Value( Label, var ) );
    }
  }
}

void Context::print() const {
  string s;
  size_t i, j;

  cout << "proven" << endl;
  for( i = 0; i < m_game.numVars(); ++i ) {

    cout << " " << i << ": ";
    if( m_proven[ i ].type() != Label ) {
      cout << m_proven[ i ].toString();
    }
    cout << endl;
  }

  cout << "falsified" << endl;
  for( i = 0; i < m_game.numVars(); ++i ) {

    cout << " " << i << ": ";
    for( j = 0; j < m_falsified[ i ].size(); ++j ) {
      if( j ) {
	cout << ",";
      }
      cout << m_falsified[ i ][ j ].toString();
    }
    cout << endl;
  }
}

void RuleTree::forgetVariables() {
  vector<bool> variableNeeded( m_game.numVars(), false );
  vector<size_t> stack;
  size_t i, j;

  /* any variable which is explicitly used is needed */
  for( i = 0; i < m_game.numVars(); ++i ) {

    for( j = 0; j < rules().size(); ++j ) {

      if( rules()[ j ].usesVariable( i ) ) {

	variableNeeded[ i ] = true;
	stack.push_back( i );
	break;
      }
    }
  }

  /* any variable which is connected to a used variable
     through a falsified equality is needed */
  while( stack.size() ) {

    i = stack.back();
    stack.pop_back();
    ASSERT( variableNeeded[ i ] );

    const vector<Value> &f = context().falsified()[ i ];

    for( j = 0; j < f.size(); ++j ) {

      if( f[ j ].type() == Label && !variableNeeded[ f[ j ].value() ] ) {

	variableNeeded[ f[ j ].value() ] = true;
	stack.push_back( f[ j ].value() );
      }
    }
  }

  for( i = 0; i < m_game.numVars(); ++i ) {

    if( !variableNeeded[ i ] ) {

      m_context.clearVar( i );
    }
  }
}

void RuleTree::setTest( const size_t ruleIdx,
			const size_t var,
			const size_t testIdx )
{
  m_testedRuleIdx = ruleIdx;

  const vector<vector<Value> > &tests = m_rules[ m_testedRuleIdx ].tests();
  if( tests[ var ][ testIdx ].type() == Label ) {

    m_test = VarTest( var, tests[ var ][ testIdx ].value() );
  } else {

    ASSERT( tests[ var ][ testIdx ].type() == Constant );
    m_test = VarTest( var );
  }
}

static bool ruleSharesVariables( const size_t numVars,
				 const Rule &rule,
				 const vector<size_t> &varUseCounts )
{
  size_t v;

  for( v = 0; v < numVars; ++v ) {

    if( rule.tests()[ v ].size() && varUseCounts[ v ] ) {

      return true;
    }
  }

  return false;
}

/* sets and vars should be empty when called */
static void getIndependentSets( const size_t numVars,
				const vector<Rule> &rules,
				vector<vector<size_t> > &sets,
				vector<vector<size_t> > &varUses,
				vector<vector<size_t> > &modifiedVarUses,
				vector<bool> &modified )
{
  if( rules.size() == 0 ) {
    return;
  }

  size_t i, j, s, v, firstMatch;
  bool foundMatch;

  for( i = 0; i < rules.size(); ++i ) {

    /* see what sets we're tangled up with */
    foundMatch = false;
    s = 0;
    while( s < sets.size() ) {

      if( ruleSharesVariables( numVars, rules[ i ], varUses[ s ] ) ) {

	if( foundMatch ) {
	  /* already found a match: merge s with first matched set */

	  /* add data from s to firstMatch */
	  for( j = 0; j < sets[ s ].size(); ++j ) {
	    sets[ firstMatch ].push_back( sets[ s ][ j ] );
	  }
	  for( v = 0; v < numVars; ++v ) {
	    varUses[ firstMatch ][ v ] += varUses[ s ][ v ];
	  }
	  for( v = 0; v < numVars; ++v ) {
	    modifiedVarUses[ firstMatch ][ v ] += modifiedVarUses[ s ][ v ];
	  }
	  if( modified[ s ] ) {
	    modified[ firstMatch ] = true;
	  }

	  /* remove set s */
	  if( s + 1 != sets.size() ) {

	    sets[ s ].swap( sets.back() );
	    varUses[ s ].swap( varUses.back() );
	    modifiedVarUses[ s ].swap( modifiedVarUses.back() );
	  }
	  sets.pop_back();
	  varUses.pop_back();
	  modifiedVarUses.pop_back();

	  continue;
	} else {

	  foundMatch = true;
	  firstMatch = s;
	}
      }

      ++s;
    }

    if( foundMatch ) {
      /* found at least one matching set - add rule to it */

      sets[ firstMatch ].push_back( i );

      /* add all the variables of rule i to the set */
      for( v = 0; v < numVars; ++v ) {

	if( rules[ i ].tests()[ v ].size() ) {

	  ++varUses[ firstMatch ][ v ];
	  if( rules[ i ].wasModified() ) {

	    ++modifiedVarUses[ firstMatch ][ v ];
	  }
	}
      }

      /* update the modified status */
      if( rules[ i ].wasModified() ) {

	modified[ firstMatch ] = true;
      }
    } else {
      /* no match found in existing set - create a new one with rule i */

      sets.push_back( vector<size_t>() );
      sets.back().push_back( i );
      varUses.push_back( vector<size_t>( numVars, 0 ) );
      modifiedVarUses.push_back( vector<size_t>( numVars, 0 ) );
      for( v = 0; v < numVars; ++v ) {

	if( rules[ i ].tests()[ v ].size() ) {

	  ++varUses.back()[ v ];
	  if( rules[ i ].wasModified() ) {

	    ++modifiedVarUses.back()[ v ];
	  }
	}
      }
      modified.push_back( rules[ i ].wasModified() );
    }
  }
}

bool RuleTree::pickTest()
{
  size_t i, v, bestS, bestV;

  if( m_rules.size() == 0 ) {
    m_testedRuleIdx = 0;
    m_test = VarTest( 0 );
    return false;
  }

  /* get the independent sets of rules */
  vector<vector<size_t> > isSets;
  vector<vector<size_t> > isVarUses;
  vector<vector<size_t> > isModifiedVarUses;
  vector<bool> isModified;
  getIndependentSets( m_game.numVars(), m_rules, isSets, isVarUses,
		      isModifiedVarUses, isModified );

  /* find the smallest modified set, or smallest unmodified set if no
      modified set exists */
  bestS = 0;
  for( i = 1; i < isSets.size(); ++i ) {

    if( isModified[ i ] ) {
      if( !isModified[ bestS ]
	  || isSets[ i ].size() < isSets[ bestS ].size() ) {

	bestS = i;
      }
    } else {
      if( !isModified[ bestS ]
	  && isSets[ i ].size() < isSets[ bestS ].size() ) {
	bestS = i;
      }
    }
  }

#if 0
  /* really simple check - pick first modified rule, first rule if
     no rules have been modified  */
  size_t bestI = 0;
  for( i = 0; i < isSets[ bestS ].size(); ++i ) {

    if( m_rules[ isSets[ bestS ][ i ] ].wasModified() ) {
      bestI = i;
      break;
    }
  }
  for( v = 0; v < m_game.numVars(); ++v ) {

    if( m_rules[ isSets[ bestS ][ bestI ] ].tests()[ v ].size() ) {

      setTest( isSets[ bestS ][ bestI ], v, 0 );
      return true;
    }
  }
  ASSERT( 0 );
#endif

#if 1
  bestV = 0;
  for( v = 0; v < m_game.numVars(); ++v ) {
    if( !isVarUses[ bestS ][ v ] ) {
      continue;
    }

    if( isModifiedVarUses[ bestS ][ v ] ) {

      if( !isModifiedVarUses[ bestS ][ bestV ]
	  || isVarUses[ bestS ][ v ] > isVarUses[ bestS ][ bestV ] ) {

	bestV = v;
      }
    } else {

      if( !isModifiedVarUses[ bestS ][ bestV ]
	  && isVarUses[ bestS ][ v ] > isVarUses[ bestS ][ bestV ] ) {
	
	bestV = v;
      }
    }
  }
#else
  /* in order of preference, test:
     something which tests all rules
     something which tests as many modified rules as possible (at least 1),
       breaking ties by testing fewer unmodified rules
     something which tests as many unmodified rules as possible */
  bestV = 0;
  for( v = 0; v < m_game.numVars(); ++v ) {
    if( !isVarUses[ bestS ][ v ] ) {
      continue;
    }

    if( isVarUses[ bestS ][ v ] == isSets[ bestS ].size() ) {
      /* v tests all rules - we have a big winner, no need to look further */

      bestV = v;
      break;
    }

    if( isModifiedVarUses[ bestS ][ v ] ) {

      if( isModifiedVarUses[ bestS ][ v ]
	  > isModifiedVarUses[ bestS ][ bestV ] ) {
	bestV = v;
      } else if( isModifiedVarUses[ bestS ][ v ]
		 == isModifiedVarUses[ bestS ][ bestV ] ) {
	if( isVarUses[ bestS ][ v ] > isVarUses[ bestS ][ bestV ] ) {
	  bestV = v;
	}
      }
    } else {

      if( !isModifiedVarUses[ bestS ][ bestV ]
	  && isVarUses[ bestS ][ v ] > isVarUses[ bestS ][ bestV ] ) {
	/* best we could do so far is test unmodified rules, and
	   v tests more rules than bestV */

	bestV = v;
      }
    }
  }
#endif

  /* pick first rule which uses the chosen 'best' variable */
  for( i = 0; i < isSets[ bestS ].size(); ++i ) {

    if( m_rules[ isSets[ bestS ][ i ] ].tests()[ bestV ].size() ) {
      setTest( isSets[ bestS ][ i ], bestV, 0 );
      return true;
    }
  }

  return false;
}

static vector<vector<Rule> > splitRulesIntoPieces( const vector<Rule> &rules,
						   const size_t numPieces )
{
  vector<vector<Rule> > pieces( 1 );
  size_t i;

  for( i = 0; i < rules.size(); ++i ) {

    if( i >= (double)( rules.size() * pieces.size() ) / (double)numPieces ) {

      pieces.push_back( vector<Rule>() );
    }

    pieces.back().push_back( rules[ i ] );
  }

  return pieces;
}

RuleTree::RuleTree( const Game &game,
		    const vector<Rule> &rules,
		    const size_t ruleGroupSize,
                    size_t* numNodesInTree,
		    size_t verbosity)
  : m_game( game ), m_id( 0 ), m_context( game ),
    m_test( VarTest( 0 ) ), m_maxChildren( 0 )
{
  size_t maxNodes, i, max;
  vector<Rule> testRules;

  /* split out any immediately usable rules */
  for( i = 0; i < rules.size(); ++i ) {

    if( rules[ i ].hasNoTests() ) {

      m_rootRules.push_back( rules[ i ] );
    } else {

      testRules.push_back( rules[ i ] );
    }
  }

  /* pick a maximum number of nodes */
  maxNodes = 0;
  for( i = 0; i < testRules.size(); ++i ) {

    for( size_t v = 0; v < m_game.numVars(); ++v ) {

      maxNodes += testRules[ i ].tests()[ v ].size();
    }
  }
  maxNodes *= 2;

  if( testRules.size() == 0 ) {
    /* if there are no rules that need to be tested, we're done */

    m_maxChildren = m_rootRules.size();
    *numNodesInTree = 1;
    return;
  }

  max = 0;
  for( i = 0; i < testRules.size(); ++i ) {
    if( testRules[ i ].id() > max ) {
      max = testRules[ i ].id();
    }
  }

  vector<vector<RuleTree *> > trees( max + 1 );
  size_t numPieces, piece, numNodes;
  bool usedTree;
  RuleTree *leafTree, *pieceTree;
  if( ruleGroupSize > 0 ) {

    numPieces = testRules.size() / ruleGroupSize;
    if( numPieces == 0 ) {

      numPieces = 1;
    }
  } else {
    numPieces = 1;
  }

if ( verbosity >= 2) { cerr << "  (RuleTree) Based on the initial rule group size, the rules with preconditions will be put into " << numPieces << " groups." << endl; }

  for( ; true; ++numPieces ) {

    for( i = 0; i <= max; ++i ) {
      trees[ i ].clear();
    }

    vector<vector<Rule> > pieces = splitRulesIntoPieces( testRules, numPieces );

    numNodes = 1; // count the root (this node)
    leafTree = NULL;
    for( piece = numPieces; piece > 0; --piece ) {

      if( piece > 1 ) {

	pieceTree = new RuleTree( m_game, m_context, pieces[ piece - 1 ] );


	pieceTree->m_id = numNodes;
	++numNodes;
      } else {

	pieceTree = this;
	m_rules = pieces[ piece - 1 ];
      }

      pieceTree->pickTest();

      usedTree = false;

if ( verbosity >= 3) { cerr << "  (RuleTree) Process group #" << piece << "." << endl; }
      if( !pieceTree->buildChildren( trees,
				     &numNodes,
				     leafTree,
				     &usedTree,
				     maxNodes ) ) {
	/* ran out of nodes... */

	break;
      }
      if( leafTree ) {
	pieceTree->m_maxChildren += leafTree->m_maxChildren;
      }

      leafTree = pieceTree;
    }

    if( piece == 0 ) {
      /* if we made it down to 0, we finished all pieces successfully */
      break;
    }
if ( verbosity >= 2) { cerr << "  (RuleTree) Rule Tree too big, re-partitioning the rules into "
	 << (numPieces+1) << " groups." << endl; }

    /* delete all of the tree nodes we created */
    for( i = 0; i <= max; ++i ) {

      for( size_t j = 0; j < trees[ i ].size(); ++j ) {

	delete trees[ i ][ j ];
      }
    }

    /* reset this node */
    m_maxChildren = 0;
    m_children.clear();
    m_childIsDup.clear();
    m_childFinishedRules.clear();

  }

if ( verbosity >= 2) { cerr << "  (RuleTree) Final Rule Tree contains " << numNodes << " nodes." << endl; }
  *numNodesInTree = numNodes;

  m_maxChildren += m_rootRules.size();
}

RuleTree::RuleTree( const Game &game,
		    const Context &context,
		    const vector<Rule> &rules )
  : m_game( game ), m_context( context ),
    m_rules( rules ), m_test( VarTest( 0 ) ), m_maxChildren( 0 )
{
  forgetVariables();
}

void RuleTree::destroyTree()
{
  for( size_t i = 0; i < m_children.size(); ++i ) {
    if( m_childIsDup[ i ] || m_children[ i ] == NULL ) {
      continue;
    }

    m_children[ i ]->destroyTree();
  }

  delete this;
}

bool RuleTree::operator==( const RuleTree &ruleTree ) const
{
  size_t i;

  if( rules().size() != ruleTree.rules().size() ) {
    return false;
  }
  for( i = 0; i < rules().size(); ++i ) {
    if( !( rules()[ i ] == ruleTree.rules()[ i ] ) ) {
      return false;
    }
  }

  if( !( context() == ruleTree.context() ) ) {
    return false;
  }

  return true;
}

static RuleTree* ruleTreeFindDuplicate( const RuleTree &tree,
					const vector<RuleTree *> trees )
{
  size_t i;

  for( i = 0; i < trees.size(); ++i ) {

    if( tree == *( trees[ i ] ) ) {
      return trees[ i ];
    }
  }

  return NULL;
}

bool RuleTree::buildChildren( vector<vector<RuleTree *> > &trees,
			      size_t *numNodes,
			      RuleTree *leafTree,
			      bool *leafTreeUsed,
			      const size_t maxNodes ) {
  size_t max, i, j;
  RuleTree *tree, *dup;

  /* check if we're a leaf node */
  if( m_rules.size() == 0 ) {
#ifdef DEBUGbc
cerr << "RH (buildChildren): exit, leafnode " << "\n";
#endif
    return true;
  }

  /* check all possible test results to create children */
  if( test().type() == Constant ) {

    max = m_game.domainSize( test().var() );
  } else {

    ASSERT( test().type() == Label );
    max = 2; /* false, true */
  }
  for( i = 0; i < max; ++i ) {
    m_childFinishedRules.push_back( vector<Rule>() );

    Context childContext( context(), test(), i );
    childContext.removeDuplicates();

    /* get the child rules */
    vector<Rule> childRules;
    for( j = 0; j < rules().size(); ++j ) {
      Rule rule = rules()[ j ];

      /* equality test can rewrite rules */
      if( test().type() == Label && i ) {
	rule.useEquivalence( test().var(), test().otherLabel() );
      }
      /* remove invalid and triggered rules */
      if( rule.useKnowledge( childContext ) ) {
	/* rule is either plausible or proven */

	if( rule.hasNoTests() ) {
	  /* rule is proven - trigger it for this child */

	  m_childFinishedRules[ i ].push_back( rule ); 
	} else {
	  /* rule is plausible - put it in the child rule tree */

	  childRules.push_back( rule );
	}
      }
    }
    if( childRules.size() == 0 ) {
      /* skip leaves */

      m_children.push_back( leafTree );
      if( leafTree ) {
	if( *leafTreeUsed ) {
	  m_childIsDup.push_back( true );
	} else {
	  *leafTreeUsed = true;
	  m_childIsDup.push_back( false );
	}
      } else {
	m_childIsDup.push_back( false );
      }
    } else {
      /* create the tree node */
      size_t childRule;

      tree = new RuleTree( m_game, childContext, childRules );
      {
          bool pickTestSuccess = tree->pickTest();
          UNUSED( pickTestSuccess );
          ASSERT( pickTestSuccess );
      }

      /* check whether it's a duplicate, and link child in */
      childRule = tree->m_rules[ tree->m_testedRuleIdx ].id();
      dup = ruleTreeFindDuplicate( *tree, trees[ childRule ] );
      if( dup != NULL ) {
	/* it's a duplicate: note this and link previous tree as child */

	delete tree;
	m_children.push_back( dup );
	m_childIsDup.push_back( true );
      } else {
	/* it's a new context/rule set: add it to list, then link as child */

	if( maxNodes && *numNodes > maxNodes ) {
	  /* ran out of nodes! */

	  delete tree;
	  return false;
	}

	
	/* link child in */
	m_children.push_back( tree );
	m_childIsDup.push_back( false );
	tree->m_id = *numNodes;

	trees[ childRule ].push_back( tree );
	++(*numNodes);

	/* expand the child */
	if( !tree->buildChildren( trees,
				  numNodes,
				  leafTree,
				  leafTreeUsed,
				  maxNodes ) ) {
	  return false;
	}
	ASSERT( tree->m_children.size() );
      }
    }

    j = m_childFinishedRules[ i ].size();
    if( m_children[ i ] != leafTree ) {
      j += m_children[ i ]->maxChildren();
    }
    if( j > m_maxChildren ) {
      m_maxChildren = j;
    }
  }

  return true;
}

void RuleTree::print() const {
  string s;
  size_t i, j;

  cout << "Tree node " << id()
       << " (maxChildren = " << maxChildren() << ")" << endl;
  context().print();

  if( id() == 0 ) {
    cout << "ROOT RULES" << endl;
    for( i = 0; i < m_rootRules.size(); ++i ) {
      cout << " " << m_rootRules[ i ].toString();
      cout << " LABEL " << m_rootRules[ i ].id() << endl;
    }
  }

  cout << "TEST RULES" << endl;
  for( i = 0; i < m_rules.size(); ++i ) {
    cout << " " << m_rules[ i ].toString();
    cout << " LABEL " << m_rules[ i ].id();
    if( m_rules[ i ].wasModified() ) {
      cout << " (M)";
    }
    cout << endl;
  }

  if( m_children.size() == 0 ) {

    cout << "leaf" << endl;
  } else {

    cout << "test" << endl;
    test().print();

    for( i = 0; i < m_children.size(); ++i ) {

      cout << endl << "child " << i << " of " << id() << endl;

      cout << "TRIGGERED RULES" << endl;
      for( j = 0; j < m_childFinishedRules[ i ].size(); ++j ) {
	cout << " " << m_childFinishedRules[ i ][ j ].toString();
	cout << " LABEL " << m_childFinishedRules[ i ][ j ].id() << endl;
      }

      if( m_childIsDup[ i ] ) {
	cout << " duplicate of " << m_children[ i ]->id() << endl;
      } else if( m_children[ i ] != NULL ) {
	m_children[ i ]->print();
      } else {
	cout << " leaf" << endl;
      }
    }
  }
}

void RuleTree::resultCall( ostream &os,
			   const string &treeName,
			   const size_t result ) const
{
  string s( "return " );

  if( m_childFinishedRules[ result ].size() ) {
    
    if( m_childFinishedRules[ result ].size() > 1 ) {

      os << "    " << nextFuncCode( actionStubName( treeName, result, 1 ) );
    } else {

      if( m_children[ result ] != NULL ) {

	os << "    " << nextFuncCode( m_children[ result ]
				      ->baseFunctionName( treeName ) );
      } else {

	os << "    " << nextFuncCode( "0" );
      }
    }
    size_t r = m_childFinishedRules[ result ][ 0 ].id() - 1;
    os << "    return " << r << ";\n";
  } else if( m_children[ result ] != NULL ) {

    os << "    return " << m_children[ result ]->baseFunctionName( treeName );
    os << "( state, next_func";
    os << " );\n";
  } else {

    os << string( "    return -1;\n" );
  }
}

void RuleTree::generateRootCode( ostream &os,
				 const string &treeName ) const
{
  string nextFuncName;
  for( size_t i = m_rootRules.size(); i > 0; --i ) {

    if( i < m_rootRules.size() ) {
      /* next entry point is next root rule */

      nextFuncName = rootActionFunctionName( treeName, i );
    } else {
      /* next entry point is root, or finished */

      if( m_children.size() == 0 ) {

	nextFuncName = string( "NULL" );
      } else {

	nextFuncName = baseFunctionName( treeName );
      }
    }

    generateActionStubCode( os,
			    rootActionFunctionName( treeName, i - 1 ),
			    treeName,
			    m_rootRules[ i - 1 ],
			    nextFuncName );
  }
}

void RuleTree::generateActionStub( ostream &os,
				   const string &treeName,
				   const size_t result,
				   const size_t i ) const
{
  string nextFuncName;

  if( i + 1 < m_childFinishedRules[ result ].size() ) {

    nextFuncName = actionStubName( treeName, result, i + 1 );
  } else if( m_children[ result ] != NULL ) {

    nextFuncName = m_children[ result ]->baseFunctionName( treeName );
  } else {

    nextFuncName = string( "NULL" );
  }

  generateActionStubCode( os,
			  actionStubName( treeName, result, i ),
			  treeName,
			  m_childFinishedRules[ result ][ i ],
			  nextFuncName );
}

void RuleTree::generateActionStubCode( ostream &os,
				       const string &funcName,
				       const string &treeName,
				       const Rule &rule,
				       const string nextFuncName ) const
{
  os << "\n" << functionHeader( funcName );
  os << "  " << nextFuncCode( nextFuncName );
  os << "  return " << rule.id() - 1 << ";\n";
  os << "}\n";
}

static bool inResultsList( const vector<Rule> &rules, const Rule &rule ) {
  size_t i;

  for( i = 0; i < rules.size(); ++i ) {
    if( rules[ i ].id() == rule.id() ) {
      return true;
    }
  }

  return false;
}

bool RuleTree::resultsSame( size_t a, size_t b ) const
{
  size_t i;

  if( m_children[ a ] != m_children[ b ] ) {
    return false;
  }

  if( m_childFinishedRules[ a ].size() != m_childFinishedRules[ b ].size() ) {
    return false;
  }

  for( i = 0; i < m_childFinishedRules[ a ].size(); ++i ) {
    if( !inResultsList( m_childFinishedRules[ b ],
			m_childFinishedRules[ a ][ i ] ) ) {
      return false;
    }
  }

  return true;
}

size_t RuleTree::getResultType( size_t *trueValue, size_t *falseValue ) const
{
  size_t i, b, aCount, bCount;

  ASSERT( test().type() == Constant );

  if( m_children.size() == 1 ) {
    return 0;
  }

  b = 0;
  aCount = 1;
  bCount = 0;
  for( i = 1; i < m_children.size(); ++i ) {

    if( resultsSame( 0, i ) ) {
      ++aCount;
    } else if( bCount ) {
      if( resultsSame( b, i ) ) {
	++bCount;
      } else {
	return 2;
      }
    } else {
      b = i;
      ++bCount;
    }
  }

  if( bCount == 0 ) {
    return 0;
  }

  if( aCount == 1 ) {
    *trueValue = 0;
    *falseValue = b;
    return 1;
  }

  if( bCount == 1 ) {
    *trueValue = b;
    *falseValue = 0;
    return 1;
  }

  return 2;
}

void RuleTree::generateFunctionCode( ostream &os,
				     const string &treeName,
                                     size_t numNodesInTree,
                                     int* varEntryId,
                                     const bool dynamicAbstractions ) const
{
  size_t i, j;

  /* Ssssssigh.  Some games have so many rules that a) they need to be
     treated one rule at a time, generating a rule tree which is one
     giant chain, and b) generating code recursively, we run out of stack */
  vector<pair<const RuleTree *,size_t> > nodeStack;
  const RuleTree *node;

  vector<VarTestNode> varTestTable(numNodesInTree);
 
  node = this;
  i = 0;
  while( 1 ) {
  generateFunction_while_continue:

    /* add the children to the stack */
    for( ; i < node->m_children.size(); ++i ) {
      if( node->m_children[ i ] != NULL && !node->m_childIsDup[ i ] ) {

	nodeStack.push_back( pair<const RuleTree *, size_t>( node, i + 1 ) );
	node = node->m_children[ i ];
	i = 0;
	goto generateFunction_while_continue;
      }
    }

    vector<int> childVarTestId(node->m_children.size());

    /* generate the action stub functions of the children */
    for( i = 0; i < node->m_children.size(); ++i ) {

        int childId = (node->m_children[i] != NULL) 
            ? node->m_children[i]->id() : -1;
        if (node->m_childFinishedRules[ i ].size()) {
            childVarTestId[i] = varTestTable.size();
            for (j = 0; j < node->m_childFinishedRules[i].size(); ++j) {
                VarTestNode stub;
                stub.type = 0;
                stub.var = node->m_childFinishedRules[i][j].id() - 1;
                if (j < node->m_childFinishedRules[i].size() - 1)
                    stub.other = varTestTable.size() + 1;
                else
                    stub.other = childId;
                varTestTable.push_back(stub);
            }            
        }
        else
            childVarTestId[i] = childId;
        
        for( j = node->m_childFinishedRules[ i ].size(); j > 1; --j ) {
            node->generateActionStub( os, treeName, i, j - 1 );
        }
    }

    if( node->m_children.size() ) {
      /* generate the function for the node */

      os << "\n";
      os << node->functionHeader( node->baseFunctionName( treeName ) );

      /* generate the test code */
      if( node->test().type() == Constant ) {
	size_t resultType, trueValue, falseValue;

	/* check if we can use an if */
	resultType = node->getResultType( &trueValue, &falseValue );
	if( resultType == 0 ) {

           /* don't even need an if... */
            VarTestNode& vt = varTestTable[ node->id() ];            
            vt.type = 0;
            vt.var = -2;
            vt.other = childVarTestId[0];
            node->resultCall( os, treeName, 0 );

	} else if( resultType == 1 ) {

            VarTestNode& vt = varTestTable[ node->id() ];
            vt.type = 2;
            vt.var = node->test().var();
            vt.other = trueValue;
            vt.rule = node->testedRuleIdx();
            vt.edges.resize(2);
            vt.edges[0] = childVarTestId[ falseValue ];
            vt.edges[1] = childVarTestId[ trueValue ];

            os << "  if( state->vars[ " << node->test().var() << " ] == "
               << trueValue << " ) {\n";
            node->resultCall( os, treeName, trueValue );
            os << "  } else {\n";
            node->resultCall( os, treeName, falseValue );
            os << "  }\n";

	} else {

            VarTestNode& vt = varTestTable[ node->id() ];
            vt.type = 3;
            vt.rule = node->testedRuleIdx();
            vt.edges.resize( node->m_children.size() );
            for( i = 0; i < node->m_children.size(); ++i )
                vt.edges[i] = childVarTestId[i];

	  os << "  switch( state->vars[ " << node->test().var() << " ] ) {\n";

	  for( i = 0; i < node->m_children.size(); ++i ) {

	    if( i + 1 < node->m_children.size() ) {

	      os << "  case " << i << ":\n";
	    } else {

	      os << "  default:\n";
	    }

	    node->resultCall( os, treeName, i );
	  }

	  os << "  }\n";
	}
      } else {
	ASSERT( node->test().type() == Label );
	ASSERT( node->m_children.size() == 2 );

        VarTestNode& vt = varTestTable[ node->id() ];
            vt.type = 1;
            vt.var = -1;
            vt.other = node->test().otherLabel();
            vt.rule = node->testedRuleIdx();
            vt.edges.resize( 2 );
            vt.edges[0] = childVarTestId[0];
            vt.edges[1] = childVarTestId[1];

	os << "  if( state->vars[ " << node->test().otherLabel()
	   << " ] == state->vars[ " << node->test().var() << " ] ) {\n";

	node->resultCall( os, treeName, 1 );

	os << "  } else {\n";

	node->resultCall( os, treeName, 0 );

	os << "  }\n";
      }

      os << "}\n";
    }

    /* generate the extra code for the root node */
    if( node->id() == 0 ) {

        int rootId = node->m_children.size() ? 0 : -1;
        if (node->m_rootRules.size()) {

            *varEntryId = varTestTable.size();

            for (j = 0; j < node->m_rootRules.size(); ++j) {
                VarTestNode stub;
                stub.type = 0;
                stub.var = node->m_rootRules[j].id() - 1;
                if (j < node->m_rootRules.size() - 1)
                    stub.other = varTestTable.size() + 1;
                else
                    stub.other = rootId;
                varTestTable.push_back(stub);
            }
        }
        else
        {
            *varEntryId = 0;  // point directly to root node            
        }

        node->generateRootCode( os, treeName );
    }

    /* get next node from stack */
    if( nodeStack.size() == 0 ) {
      break;
    }
    node = nodeStack.back().first;
    i = nodeStack.back().second;
    nodeStack.pop_back();
  }

}


void Game::generateDomainCode( ostream &os, const bool useAbstraction ) const
{
  size_t i, j;

  os << "#define NUMDOMAINS " << m_domains.size() << "\n";

  os << "static var_t domain_sizes[ NUMDOMAINS ] = { ";
  for( i = 0; i < m_domains.size(); ++i ) {

    if( i ) {
      os << ", ";
    }
    os << m_domains[ i ].size();
  }
  os << " };\n";

  if ( useAbstraction ) {
     os << "static const char *name_of_domain[ NUMDOMAINS ] = { ";
     for (i = 0; i < m_domains.size(); ++i) {
         if ( i ) {
             os << ", ";
         }
         os << '\"' << m_domainNames[i] << '\"';
     }
     os << " };\n";
  }

  os << "static int var_domains[ NUMVARS ] = { ";
  for( i = 0; i <numVars(); ++i ) {

    if( i ) {
      os << ", ";
    }
    os << m_variableDomains[ i ];
  }
  os << " };\n";

  for( i = 0; i < m_domains.size(); ++i ) {

    os << "static const char *domain_" << i
       << "[ " << m_domains[ i ].size() << " ] = {";
    for( j = 0; j < m_domains[ i ].size(); ++j ) {

      if( j ) {
	os << ", ";
      }
      os << "\"" << m_domains[ i ][ j ] << "\"";
    }

    os << " };\n";
  }

  if ( useAbstraction ) {
     os << "static const char **domain_to_domain_names[ NUMDOMAINS ] = { ";
     for( i = 0; i < m_domains.size(); ++i ) {
       if( i ) {
         os << ", ";
       }
       os << "domain_" << i;
     }
     os << " };\n";
  }

  os << "static const char **var_domain_names[ NUMVARS ] = { ";
  for( i = 0; i <numVars(); ++i ) {

    if( i ) {
      os << ", ";
    }
    os << "domain_" << m_variableDomains[ i ];
  }
  os << " };\n";
}

void Game::generateRuleNames( ostream &os, const bool bwdMoves ) const
{
  string s;
  size_t i;

  os << "static const char *fwd_rule_name[ " << rules().size() <<" ] = { ";
  for( i = 0; i < rules().size(); ++i ) {

    if( i ) { 
      os << ", ";
    }
    os << "\"" << rules()[ i ].name() << "\"";
  }
  os << " };\n";
  os << "#define get_fwd_rule_label( ruleid ) (fwd_rule_name[(ruleid)]) \n";

  if( bwdMoves ) {
    os << "static const char *bwd_rule_name[ "
       << backwardsRules().size() << " ] = { ";
    for( i = 0; i < backwardsRules().size(); ++i ) {

      if( i ) { 
	os << ", ";
      }
      os << "\"" << backwardsRules()[ i ].name() << "\"";
    }
    os << " };\n";
    os << "#define get_bwd_rule_label( ruleid ) (bwd_rule_name[(ruleid)]) \n";
  }
}

void Game::generateRuleCosts( ostream &os, const bool bwdMoves ) const
{
  string s;
  size_t i;

  int mincost = INT_MAX;
  int maxcost = -1;
  for( i = 0; i < rules().size(); ++i ) {
    int cost = static_cast<int>(rules()[ i ].cost());
    if (cost < mincost) { mincost = cost; }
    if (cost > maxcost) { maxcost = cost; }
  }

  os << "#define cost_of_cheapest_fwd_rule " << mincost << "\n";
  os << "#define cost_of_costliest_fwd_rule " << maxcost << "\n";
  if ( mincost == maxcost ) {
     os << "#define get_fwd_rule_cost( ruleid ) " << maxcost << "\n";
  } else {
     os << "static const int fwd_rule_cost[ " << rules().size() <<" ] = { ";
     for( i = 0; i < rules().size(); ++i ) {
       if( i ) { 
         os << ", ";
       }
       int cost = static_cast<int>(rules()[ i ].cost());
       os << cost;
     }
     os << " };\n";
     os << "#define get_fwd_rule_cost( ruleid ) (fwd_rule_cost[(ruleid)]) \n";
  }
  

  if( bwdMoves ) {
      int mincost = INT_MAX;
      int maxcost = -1;
      for( i = 0; i < backwardsRules().size(); ++i ) {
          int cost = static_cast<int>(backwardsRules()[ i ].cost());
          if (cost < mincost) { mincost = cost; }
          if (cost > maxcost) { maxcost = cost; }
      }
      os << "#define cost_of_cheapest_bwd_rule "  << mincost << "\n";
      os << "#define cost_of_costliest_bwd_rule "  << maxcost << "\n";
      if ( mincost == maxcost ) {
         os << "#define get_bwd_rule_cost( ruleid ) " << maxcost << "\n";
      } else {
         os << "static const int bwd_rule_cost[ "
            << backwardsRules().size() << " ] = { ";
         for( i = 0; i < backwardsRules().size(); ++i ) {
             if( i ) { 
                 os << ", ";
             }
             int cost = static_cast<int>(backwardsRules()[ i ].cost());
             os << cost;
         }
         os << " };\n";
         os << "#define get_bwd_rule_cost( ruleid ) (bwd_rule_cost[(ruleid)]) \n";
      }
  }
}

/* The following is needed for one of the abstraction-related functions */
void Game::generateRuleLabelSets( ostream &os, const bool bwdMoves ) const
{
    os << "\nstatic int fwd_rule_label_sets["
       << m_rules.size() * numVars() << "] = {";
    for (size_t j = 0; j < m_rules.size(); ++j) {
        const Rule& rule = m_rules[j];
        for (size_t i = 0; i < numVars(); ++i) {
            if (i || j) os << ',';
            if (rule.tests()[i].size() == 1
                && rule.tests()[i][0].type() == Label)
                os << rule.tests()[i][0].value();
            else
                os << i;
        }
    }
    os << "};\n";

    if (bwdMoves) {
        os << "\nstatic int bwd_rule_label_sets["
           << m_bwdRules.size() * numVars() << "] = {";
        for (size_t j = 0; j < m_bwdRules.size(); ++j) {
            const Rule& rule = m_bwdRules[j];
            for (size_t i = 0; i < numVars(); ++i) {
                if (i || j) os << ',';
                if (rule.tests()[i].size() == 1
                    && rule.tests()[i][0].type() == Label)
                    os << rule.tests()[i][0].value();
                else
                    os << i;
            }
        }
        os << "};\n";
    }
}



string Game::generateGoalCode( const bool bwdMoves ) const
{
  string s( "/* returns 1 if state is a goal state, 0 otherwise */\n"
            "static int is_goal( const state_t *state )\n{\n" );
  char cStr[ 64 ];
  size_t i, j, k;
  bool firstTest;

  for( i = 0; i < m_goals.size(); ++i ) {

    s.append( "  if(" );

    firstTest = true;
    for( j = 0; j < numVars(); ++j ) {

      if( m_goals[ i ][ j ].size() ) {
	ASSERT( m_goals[ i ][ j ].size() == 1 );

	if( firstTest ) {
	  firstTest = false;
	} else {
	  s.append( " &&" );
	}

	sprintf( cStr, " state->vars[ %zu ] == ", j );
	s.append( cStr );

	if( m_goals[ i ][ j ][ 0 ].type() == Constant ) {

	  sprintf( cStr, "%zu", m_goals[ i ][ j ][ 0 ].value() );
	  s.append( cStr );
	} else {
	  ASSERT( m_goals[ i ][ j ][ 0 ].type() == Label );

	  sprintf( cStr, "state->vars[ %zu ]", m_goals[ i ][ j ][ 0 ].value() );
	  s.append( cStr );
	}
      }
    }

    ASSERT( firstTest == false );

    s.append( " ) {\n    return 1;\n  }\n" );
  }

  s.append( "  return 0;\n}\n\n" );

  if( bwdMoves ) {
      s.append( "static void init_goal_state( state_t *state, int goal_rule )\n{\n" );
      s.append( "  switch( goal_rule ) {\n" );
      for( i = 0; i < m_goals.size(); ++i ) {

        sprintf( cStr, "  case %zu:\n", i );
        s.append( cStr );
        for( j = 0; j < numVars(); ++j ) {

          if( m_goals[ i ][ j ].size() ) {
	    if( m_goals[ i ][ j ][ 0 ].type() == Constant ) {

	      sprintf( cStr, "    state->vars[ %zu ] = %zu;\n", j,
		       m_goals[ i ][ j ][ 0 ].value() );
	    } else {

	      ASSERT( m_goals[ i ][ j ][ 0 ].type() == Label );
	      sprintf( cStr, "    state->vars[ %zu ] = 0;\n", j );
	    }
          } else {

	    sprintf( cStr, "    state->vars[ %zu ] = 0;\n", j );
          }
          s.append( cStr );
        }
        s.append( "    break;\n" );
      }
      s.append( "  }\n" );
      s.append( "}\n\n" );

      s.append( "/* get the first goal state and initialise iterator */\n" );
      s.append( "#define first_goal_state( state_ptr, int_ptr_goal_iter ) init_goal_state(state_ptr,*(int_ptr_goal_iter)=0)\n\n" );

      s.append( "/* get the next goal state\n" );
      s.append( "   returns 1 if there is another goal state, 0 otherwise */\n" );
      s.append( "static int8_t next_goal_state( state_t *state, int *goal_iter )\n{\n" );
      s.append( "  switch( *goal_iter ) {\n" );
      for( i = 0; i < m_goals.size(); ++i ) {

        sprintf( cStr, "  case %zu:\n", i );
        s.append( cStr );

    /* increment free variables */
        k = numVars();
        for( j = 0; j < numVars(); ++j ) {
          if( m_goals[ i ][ j ].size() ) { continue; }

          if( k < numVars() ) {
	    sprintf( cStr, "    state->vars[ %zu ] = 0;\n", k );
	    s.append( cStr );
          }

          sprintf( cStr, "    if( state->vars[ %zu ] < %zu ) {\n",
	           j, domainSize( j ) - 1 );
          s.append( cStr );
          sprintf( cStr, "      ++state->vars[ %zu ];\n", j );
          s.append( cStr );

      /* fix labeled variables */
          for( k = 0; k < numVars(); ++k ) {

	    if( m_goals[ i ][ k ].size()
	        && m_goals[ i ][ k ][ 0 ].type() == Label
	        && m_goals[ i ][ k ][ 0 ].value() <= j ) {

	      sprintf( cStr, "      state->vars[ %zu ] = state->vars[ %zu ];\n",
		       k, m_goals[ i ][ k ][ 0 ].value() );
	      s.append( cStr );
	    }
          }

          s.append( "      return 1;\n    }\n" );
          k = j;
        }

        if( i + 1 < m_goals.size() ) {
          s.append( "    ++(*goal_iter);\n" );
          s.append( "    init_goal_state( state, *goal_iter );\n" );
          s.append( "    return 1;\n" );
        } else {
          s.append( "    return 0;\n" );
        }
      }
      s.append( "  }\n" );
      s.append( "  return 0;\n" );
      s.append( "}\n" );
  }
  return s;
}



static Rule ruleFromTokens( const string name,
                            const size_t cost, 
			    const vector<string> &tests,
			    const vector<string> &actions,
			    const PSVN &psvn )
{
  vector<vector<string> > labelsInDomain( psvn.numDomains() );
  vector<vector<size_t> > varsForLabels( psvn.numDomains() );
  vector<vector<Value> > ruleTests( psvn.numVariables() );
  vector<vector<Value> > ruleActions( psvn.numVariables() );
  size_t i, domain_i, val;

  /* process the tests */
  for( i = 0; i < psvn.numVariables(); ++i ) {
    if( tests[ i ].compare( "-" ) == 0 
	|| tests[ i ][ 0 ] == '*' ) {
      /* variable is ignored, nothing to do */
      continue;
    }

    string test = tests[ i ];

    domain_i = psvn.variableDomain( i );

    val = psvn.getDomainValue( domain_i, test );
    if( val < psvn.domainSize( domain_i ) ) {
      /* token is a value */

      ruleTests[ i ].push_back( Value( Constant, val ) );
    } else {
      /* token is a label */

      bool foundLabel = false;
      for( val = 0; val < labelsInDomain[ domain_i ].size(); ++val ) {
	if( labelsInDomain[ domain_i ][ val ].compare( test ) == 0 ) {
	  foundLabel = true;
	  break;
	}
      }
      if( foundLabel ) {
	/* token is a label we've already seen */

	ruleTests[ i ]
	  .push_back( Value( Label, varsForLabels[ domain_i ][ val ] ) );
      } else {
	/* token is the first time we've seen this label */

	labelsInDomain[ domain_i ].push_back( test );
	varsForLabels[ domain_i ].push_back( i );
      }
    }
  }

  /* process the actions - ALMOST the same as the tests... */
  for( i = 0; i < psvn.numVariables(); ++i ) {
    if( actions[ i ].compare( "-" ) == 0 ) {
      /* variable untouched by action, so nothing to do */
      continue;
    }

    string action;
    if( actions[ i ][ 0 ] == '*' ) {
      action = actions[ i ].substr( 1 );
    } else {
      action = actions[ i ];
    }

    domain_i = psvn.variableDomain( i );

    val = psvn.getDomainValue( domain_i, action );
    if( val < psvn.domainSize( domain_i ) ) {
      /* token is a value */

      ruleActions[ i ].push_back( Value( Constant, val ) );
    } else {
      /* token is a label */

      bool foundLabel = false;
      for( val = 0; val < labelsInDomain[ domain_i ].size(); ++val ) {
	if( labelsInDomain[ domain_i ][ val ].compare( action ) == 0 ) {
	  foundLabel = true;
	  break;
	}
      }
      if( foundLabel ) {
	/* token references a label we've already seen */

	ruleActions[ i ]
	  .push_back( Value( Label, varsForLabels[ domain_i ][ val ] ) );
      } else {
	/* token is the first time we've seen this label */

	labelsInDomain[ domain_i ].push_back( action );
	varsForLabels[ domain_i ].push_back( i + psvn.numVariables() ); /* ROB SUGGESTS: change i to i + psvn.numVariables() */
	ruleActions[ i ].push_back( Value( Label, i + psvn.numVariables() ) );
      }
    }
  }

  /* create the rule - we'll set the id later... */
  return Rule( 0, name, cost, ruleTests, ruleActions );
}
		    
static vector<vector<Value> > goalFromTokens
( const vector<string> &tokens, const PSVN &psvn )
{
  vector<vector<string> > labelsInDomain( psvn.numDomains() );
  vector<vector<size_t> > varsForLabels( psvn.numDomains() );
  vector<vector<Value> > tests( psvn.numVariables() );
  size_t i, domain_i, val;

  /* process the tests */
  for( i = 0; i < psvn.numVariables(); ++i ) {
    if( tokens[ i ].compare( "-" ) == 0 ) {
      /* token is '-', so nothing to do */
      continue;
    }

    domain_i = psvn.variableDomain( i );

    val = psvn.getDomainValue( domain_i, tokens[ i ] );
    if( val < psvn.domainSize( domain_i ) ) {
      /* token is a value */

      tests[ i ].push_back( Value( Constant, val ) );
    } else {
      /* token is a label */

      bool foundLabel = false;
      for( val = 0; val < labelsInDomain[ domain_i ].size(); ++val ) {
	if( labelsInDomain[ domain_i ][ val ].compare( tokens[ i ] ) == 0 ) {
	  foundLabel = true;
	  break;
	}
      }
      if( foundLabel ) {
	/* token references a label we've already seen */

	tests[ i ]
	  .push_back( Value( Label, varsForLabels[ domain_i ][ val ] ) );
      } else {
	/* token is the first time we've seen this label */

	labelsInDomain[ domain_i ].push_back( tokens[ i ] );
	varsForLabels[ domain_i ].push_back( i );
      }
    }
  }

  return tests;
}

static bool ruleIsDuplicate( const Rule &rule,
			     const vector<Rule> &rules )
{
  for( size_t i = 0; i < rules.size(); ++i ) {
    if( rule == rules[ i ] ) {
      return true;
    }
  }
  return false;
}

static void handleGameRule( const bool removeDuplicateRules,
			    const Rule &rule,
			    const vector<vector<string> > &domains,
			    const vector<size_t> variableDomains,
			    vector<Rule> &rules )
{
  if ( rule.isEmpty() )
      return;
  size_t i;
  vector<pair<size_t,size_t> > unboundVars = rule.unboundActionVariables();
  if( unboundVars.size() == 0 ) {
    /* no unbound action variables, no need to re-write */
    if( !removeDuplicateRules
	|| !ruleIsDuplicate( rule, rules ) ) {
      rules.push_back( rule );
    }
  } else {
    /* re-write as a set of rule to eliminate unbound variables */
    vector<size_t> unboundVals( unboundVars.size() );

    for( i = 0; i < unboundVars.size(); ++i ) {
      unboundVals[ i ]
	= domains[ variableDomains[ unboundVars[ i ].second ] ].size() - 1;
    }
    while( 1 ) {

      /* add the re-written rule */
      Rule constructedRule = Rule( rule, unboundVars, unboundVals );
      if( !removeDuplicateRules
	  || !ruleIsDuplicate( constructedRule, rules ) ) {
	rules.push_back( constructedRule );
      }

      /* get next variable assignments */
      /* find the rightmost non-zero value */
      i = unboundVars.size() - 1;
      while( i > 0 && unboundVals[ i ] == 0 ) {
	unboundVals[ i ]
	  = domains[ variableDomains[ unboundVars[ i ].second ] ].size() - 1;
	--i;
      }
      if( unboundVals[ i ] ) {
	/* decrement the rightmost non-zero value */

	--unboundVals[ i ];
      } else {
	/* no rightmost non-zero value, so we're done */

	break;
      }
    }
  }
}

static string newTempVariableName( size_t *i,
				   const vector<string> &tokensA,
				   const vector<string> &tokensB,
				   const vector<string> &tokensC )
{
  bool okay;
  size_t j;
  string name;
  char cStr[ 64 ];

  for( ; 1; ++(*i) ) {

    /* create a name */
    sprintf( cStr, "T%zu", *i );
    name = string( cStr );

    okay = true;

    /* see if the name is in token list a */
    for( j = 0; j < tokensA.size(); ++j ) {

      if( name.compare( tokensA[ j ] ) == 0 ) {

	okay = false;
	break;
      }
    }
    if( !okay ) {
      /* duplicate name - keep looking */

      continue;
    }

    /* see if the name is in token list b */
    for( j = 0; j < tokensB.size(); ++j ) {

      if( name.compare( tokensB[ j ] ) == 0 ) {

	okay = false;
	break;
      }
    }
    if( !okay ) {
      /* duplicate name - keep looking */

      continue;
    }

    /* see if the name is in token list c */
    for( j = 0; j < tokensC.size(); ++j ) {

      if( name.compare( tokensC[ j ] ) == 0 ) {

	okay = false;
	break;
      }
    }
    if( okay ) {
      /* name is unique - stop looking */
      break;
    }
  }

  return name;
}

Game *Game::fromPSVN( const PSVN &psvn,
		      const bool makeBackwardsMoves,
		      const bool removeDuplicateRules,
		      int verbosity )
{
  Game *game;
  size_t i, j;

  /* create the rules */
  size_t tempIdx;
  vector<Rule> rules, bwdRules;
  vector<string> bwdRuleTests, bwdRuleActions;
  tempIdx = 0;
  for( i = 0; i < psvn.numRules(); ++i ) {
  if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to deal with another forward rule\n"; }
    /* forward rules */
    handleGameRule( removeDuplicateRules,
		    ruleFromTokens( psvn.ruleLabel( i ),
                                    psvn.ruleCost( i ),
				    psvn.ruleTests( i ),
				    psvn.ruleActions( i ),
				    psvn ),
		    psvn.domainValues(),
		    psvn.variableDomains(),
		    rules );
    if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to test makeBackwardsMoves\n"; }
    if( !makeBackwardsMoves ) {
      if ( verbosity >= 4 ) { cerr << "  (fromPSVN) makeBackwardsMoves is false\n"; }
      continue;
    }

    /* backwards rules */
    if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to enter the code for making BackwardsMoves\n"; }
    bwdRuleTests.clear();
    bwdRuleActions.clear();
    for( j = 0; j < psvn.numVariables(); ++j ) {
      if( psvn.ruleActions( i )[ j ].compare( "-" ) == 0 ) {

	if( psvn.ruleTests( i )[ j ][ 0 ] == '*' ) {
	  bwdRuleTests.push_back( psvn.ruleTests( i )[ j ].substr( 1 ) );
	} else {
	  bwdRuleTests.push_back( psvn.ruleTests( i )[ j ] );
	}
	bwdRuleActions.push_back( psvn.ruleActions( i )[ j ] );
      } else {

	bwdRuleTests.push_back( psvn.ruleActions( i )[ j ] );
	if( psvn.ruleTests( i )[ j ].compare( "-" ) == 0 ) {
	  bwdRuleActions.push_back
	    ( newTempVariableName( &tempIdx,
				   psvn.ruleTests( i ),
				   psvn.ruleActions( i ),
				   bwdRuleActions ) );
	} else {
	  bwdRuleActions.push_back( psvn.ruleTests( i )[ j ] );
	}
      }
    }
    if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to call handleGameRule with the BackwardsMoves\n"; }
    handleGameRule( removeDuplicateRules,
		    ruleFromTokens( psvn.ruleLabel( i ),
                                    psvn.ruleCost( i ),
				    bwdRuleTests,
				    bwdRuleActions,
				    psvn ),
		    psvn.domainValues(),
		    psvn.variableDomains(),
		    bwdRules );
  }

  if( rules.size() == 0 ) {
    cerr << "warning: your state space contains no non-trivial operators." << endl;
  }
  if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to create goals\n"; }
  /* create the goals */
  vector<vector<vector<Value> > > goals;
  for( i = 0; i < psvn.numGoalRules(); ++i ) {
    goals.push_back( goalFromTokens( psvn.goalRule( i ), psvn ) );
  }

  /* we've ignored the ids because of unbound variables and duplicate
     detection, so assign ids to the rules now */
  for( i = 0; i < rules.size(); ++i ) {
    rules[ i ].setId( i + 1 );
  }

  if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to assign ids to backwards rules\n"; }
  for( i = 0; i < bwdRules.size(); ++i ) {
    bwdRules[ i ].setId( i + 1 );
  }
  if ( verbosity >= 4 ) { cerr << "  (fromPSVN) about to build the game\n"; }
  /* build the game and return it */
  game = new Game( psvn.domainNames(),
                   psvn.domainValues(),
		   psvn.variableDomains(),
		   rules,
		   bwdRules,
		   goals );
  if ( verbosity >= 4 ) { cerr << "  (fromPSVN) finished building the game\n"; }
  return game;
}


void PruneTree::destroyTree( vector<PruneTree *> &allNodes )
{
  size_t i;

  for( i = m_children.size(); i > 0; --i ) {
    if( m_children[ i - 1 ] == NULL ) {
      continue;
    }

    m_children[ i - 1 ]->destroyTree( allNodes );
  }

  ASSERT( allNodes[ m_id ] == this );
  allNodes[ m_id ] = NULL;

  delete this;
}

/* returns true if all situations in a are covered by b
   That is, the conditions under which a are generated are a subset
   of the conditions under which b are generated, and the conditions
   of a, applied to b, contain a */
bool aCoveredByB( const vector<Value> &a,
		  const vector<Value> &a_c,
		  const vector<Value> &b,
		  const vector<Value> &b_c )
{
  size_t i;

  ASSERT( a.size() == a_c.size() );
  ASSERT( b.size() == b_c.size() );
  if( a.size() != b.size() ) {

    return false;
  }

  /* need to check that context b_c is sufficient to cover a_c
     (if it's not, the states b are generated from a straict subset of
     parent states, and thus does not cover all cases of a) */
  for( i = 0; i < a_c.size(); ++i ) {

    if( a_c[ i ].type() == Constant ) {
      if( b_c[ i ].type() == Constant ) {
        if ( a_c[ i ].value() != b_c[ i ].value() ) {
	/* an unknown would match a constant, but a different constant
	   means the context of b doesn't include the context of a */
	  return false;
        }
      } else {
        ASSERT( b_c[ i ].type() == Label );

        if( b_c[ i ].value() != i ) {
          ASSERT( b_c[ i ].value() < i );
          ASSERT( b_c[ b_c[ i ].value() ].type() == Label );
          /* bi==bj so it must be true that aj==ai */
          if( !( a_c[ b_c[ i ].value() ] == a_c[ i ] ) ) {
            return false;
          }
        }
      }
    } else {
      ASSERT( a_c[ i ].type() == Label );

      if( b_c[ i ].type() == Constant ) {
	/* a_c[ i ] an unknown, b_c[ i ] a constant, so a_c is more general */

	return false;
      }
      ASSERT( b_c[ i ].type() == Label );

      if( a_c[ i ].value() == i ) {
	/* a_c[ i ] is the original unmodified unknown for variable i,
	   so b_c[ i ] must not be equal to some other variable */

	if( b_c[ i ].value() != i ) {

	  return false;
	}
      } else {
	/* a_c[ i ] is equal to some other unknown */
	ASSERT( a_c[ i ].value() < i );
	ASSERT( a_c[ a_c[ i ].value() ].type() == Label );

	if( b_c[ i ].value() != i /* ai==aj, bi free */
	    && b_c[ i ].value() != a_c[ i ].value() /* ai==aj, bi==bj */
	    && ( a_c[ b_c[ i ].value() ].type() != Label
		 || a_c[ b_c[ i ].value() ].value() != a_c[ i ].value() )
	    /* ai==ak==aj, bi==bk */ ) {
	  /* b_c[ i ] is equal to a variable, and it's a different variable */

	  return false;
	}
      }
    }
  }

  /* apply the context of a_c to b */
  vector<Value> mod_b( b );
  for( i = 0; i < b.size(); ++i ) {

    if( mod_b[ i ].type() == Label ) {

      mod_b[ i ]  = a_c[ mod_b[ i ].value() ];
    }
  }

  for( i = 0; i < a.size(); ++i ) {

    if( !( a[ i ] == mod_b[ i ] ) ) {

      return false;
    }
  }

  return true;
}


/* ALL THE REAL WORK FOR CREATING THE OPERATOR-SEQUENCE TREE for move pruning IS DONE HERE */
void PruneTree::expand( const vector<Rule> &rules,
			vector<PruneTree *> &allNodes,
			int verbosity )
{
  size_t i, j, childCost;
  bool addChild;
  vector<Value> childValues;
  vector<Value> childContext;

  m_children.resize( rules.size(), NULL );
  for( i = 0; i < rules.size(); ++i ) {


    if( m_suffix->m_children[ i ] == NULL ) {
      /* we've already pruned rule i */
       if ( verbosity==V_DEBUG_MP ) {
          cerr << "      (PruneTree::expand): no need to extend the sequence with rule i=" << i << "\n";
       }
      continue;
    }

    childCost = m_cost + rules[ i ].cost();
    childValues = m_values;
    childContext = m_context;
    if( rules[ i ].applyTo( childValues, childContext ) ) {
      /* rule can be applied - have a potentially new child */

      /* check if the child can be generated by some other move sequence */
      addChild = true;
      for( j = 0; j < allNodes.size(); ++j ) {
	if( allNodes[ j ] == NULL ) {
	  continue;
	}

	/* CHECK IF THE NEW NODE IS REDUNDANT WITH A PREVIOUSLY GENERATED NODE */
	if( allNodes[ j ]->m_cost <= childCost 
	    && aCoveredByB( childValues,
			    childContext,
			    allNodes[ j ]->values(),
			    allNodes[ j ]->context() ) ) {
	  /* child is redundant with a previously generated rule sequence */
          if ( verbosity==V_DEBUG_MP ) {
             cerr << "      (PruneTree::expand): extending the sequence with rule i=" << i << " creates a redundant sequence\n";
          }

	  addChild = false;
	  break;
#ifdef UnsafeVersion
	/* XXX THE FOLLOWING CODE, FROM THE ORIGINAL (2011) VERSION, IS UNSAFE, IT MUST NOT BE DONE */
	/* CHECK IF A PREVIOUSLY GENERATED NODE IS REDUNDANT WITH THE NEW NODE */
	} else if( childCost <= allNodes[ j ]->m_cost
		   && aCoveredByB( allNodes[ j ]->values(),
				   allNodes[ j ]->context(),
				   childValues,
				   childContext ) ) {
	  /* child is more general than other rule sequence - remove it XXX NOT SAFE */

	  allNodes[ j ]->removeFromTree( allNodes );
	  addChild = true;
#endif
	}
      }

      if( addChild ) {
	/* child is not redundant - add it */

#ifdef UnsafeVersion
	/* XXX THE FOLLOWING TEST CAN ONLY BE TRUE IN THE UNSAFE (2011) VERSION. */
	if( m_suffix->m_children[ i ] == NULL ) {
	  /* we've pruned m_suffix of the child (the path that we will
	     remember after forgetting the first move in the sequence) */

	  PruneTree *root = this;
	  while( root->m_parent ) {
	    root = root->m_parent;
	  }
	  m_children[ i ] = new PruneTree( this,
					   i,
					   root,
					   allNodes.size(),
					   childCost,
					   childValues,
					   childContext );
	  allNodes.push_back( m_children[ i ] );
	  m_children[ i ]->expand( rules, allNodes );
	} else {
#endif
	  /* normal case */
	  assert(m_suffix->m_children[ i ] != NULL);
	  m_children[ i ] = new PruneTree( this,
					   i,
					   m_suffix->m_children[ i ],
					   allNodes.size(),
					   childCost,
					   childValues,
					   childContext );
	  allNodes.push_back( m_children[ i ] );
#ifdef UnsafeVersion
	}
#endif
      } else {   /* child is redundant, so don't add it */

	m_children[ i ] = NULL;
      }
    } else {

      m_children[ i ] = NULL;
    }
  }
}


/* DESCEND TO DEPTH depth (along valid sequences) AND ADD A NEW LAYER TO THE OPERATOR-SEQUENCE TREE */
void PruneTree::expandTo( const vector<Rule> &rules,
			  const size_t depth,
			  vector<PruneTree *> &allNodes,
			  int verbosity )
{
  if( depth == 0 ) {
    ASSERT( m_children.size() == 0 );
    expand( rules, allNodes, verbosity );
    return;
  }

  size_t i;
  for( i = 0; i < m_children.size(); ++i ) {
    if( m_children[ i ] == NULL ) {
      continue;
    }

    ASSERT( m_suffix->m_children[ i ] != NULL );
    if ( verbosity==V_DEBUG_MP ) {
       cerr << "    (PruneTree::expandTo): (d=" << depth << ",i=" << i << ")" << endl;
    }
    m_children[ i ]->expandTo( rules, depth - 1, allNodes, verbosity );
  }
}



void PruneTree::removeFromTree( vector<PruneTree *> &allNodes )
{
  ASSERT( allNodes[ m_id ] == this );
  ASSERT( m_parent != NULL );
  ASSERT( m_parent->m_children[ m_parentRuleUsed ] == this );
  m_parent->m_children[ m_parentRuleUsed ] = NULL;

  for( size_t i = 0; i < allNodes.size(); ++i ) {
    if( allNodes[ i ] == NULL ) {
      continue;
    }

    if( allNodes[ i ]->m_suffix == this ) {
      /* node i has this as a suffix, so it should be removed too */

      allNodes[ i ]->removeFromTree( allNodes );
    }
  }

  allNodes[ m_id ]->destroyTree( allNodes );
}

void PruneTree::print( const size_t depth ) const
{
  size_t i;

  for( i = 0; i < depth; ++i ) {

    cout << " ";
  }
  cout << m_id;
  if( m_parent != NULL ) {

    cout << "(" << m_parentRuleUsed << ")";
  }
  cout << ":";

  for( i = 0; i < m_values.size(); ++i ) {

    cout << " " << m_values[ i ].toString();
  }

  cout << "  context(";
  for( i = 0; i < m_context.size(); ++i ) {

    cout << " " << m_context[ i ].toString();
  }
  cout << ")" << endl;

  cout << endl;

  for( i = 0; i < m_children.size(); ++i ) {
    if( m_children[ i ] == NULL ) {
      continue;
    }

    m_children[ i ]->print( depth + 1 );
  }
}

void PruneTree::fillTable( vector<int> &table ) {
  size_t i;

  if( m_children.size() == 0 ) {
    return;
  }

  for( i = 0; i < m_children.size(); ++i ) {

  /* REMINDER: table[ m_id * m_children.size() + i ] is effectively row m_id, column i */

    if( m_children[ i ] == NULL ) {
      /* child i is pruned */

      table[ m_id * m_children.size() + i ] = 0;  /* 0 means action i is not allowed in context m_id */
    } else if( m_children[ i ]->m_children.size() == 0 ) {
      /* we've traversed the whole history, so use the suffix */

      table[ m_id * m_children.size() + i ]
	= m_children[ i ]->m_suffix->m_id * m_children.size();
      m_children[ i ]->fillTable( table );
    } else {
      /* we haven't made enough moves to truncate the history */

      table[ m_id * m_children.size() + i ]
	= m_children[ i ]->m_id * m_children.size();
      m_children[ i ]->fillTable( table );
    }
  }
}



/* move pruning is done in this constructor */
PruneTree::PruneTree( const Game &game,
		      const vector<Rule> &rules,
		      const size_t depth,
		      int verbosity )
  : m_parent( NULL ), m_parentRuleUsed( 0 ), m_id( 0 ), m_cost( 0 )
{
  size_t d, i;

  m_values.reserve( game.numVars() );
  m_context.reserve( game.numVars() );

/* set each state variable to be a distinct "variable" (Label), representing the generic start state */
  for( i = 0; i < game.numVars(); ++i ) {
    m_values.push_back( Value( Label, i ) );
    m_context.push_back( Value( Label, i ) );
  }

  vector<PruneTree *> allNodes;
  allNodes.push_back( this );

  /* create a dummy node above the root with all moves leading to the
     root in order to make the suffixes easier to deal with */
  m_suffix = new PruneTree( NULL, 0, NULL, 0, 0, m_values, m_context );
  m_suffix->m_children.resize( rules.size(), this );

  if ( verbosity==V_DEBUG_MP ) {
     cerr << "  (PruneTree): about to start building tree of depth=" << depth << "\n";
  }

/* ITERATIVE DEEPENING DEPTH-FIRST GENERATION OF THE OPERATOR-SEQUENCE TREE */
  for( d = 0; d <= depth; ++d ) {
    if ( verbosity==V_DEBUG_MP ) {
       cerr << "  (PruneTree): about to call expandTo with d=" << d << "\n";
    }
    expandTo( rules, d, allNodes, verbosity );
  }
}



vector<int> PruneTree::makeTable()
{
  size_t num;
  vector<int> table;

  /* expansion process leaves holes in the set of IDs
     re-label everything to remove the holes */
  num = 0;
  reLabel( &num );

  ASSERT( m_children.size() != 0 );
  table.resize( num * m_children.size() );  /* effectively it has num rows and m_children.size() columns */
  fillTable( table );

  return table;
}


static void printRuleActions( const string &treeName,
			      const vector<Rule> &rules,
                              const vector<size_t>& varDomains,
                              const bool dynamicAbstractions )
{
  for( size_t i = 0; i < rules.size(); ++i ) {
      rules[ i ].generateActionFunction( cout, treeName, varDomains,
                                         dynamicAbstractions );
  }
  cout << "\nstatic actfuncptr " << treeName << "_rules[ " << rules.size() << " ] = { ";
  for( size_t i = 0; i < rules.size(); ++i ) {

    if( i ) {
      cout << ", ";
    }
    cout << rules[ i ].actionFunctionName( treeName );
  }
  cout << " };\n";
}

static void printPruneTable( const vector<int> &table, const string prefix ) {
  size_t i;

  cout << "static int " << prefix << "_prune_table[ " << table.size() << " ] = {";
  for( i = 0; i < table.size(); ++i ) {
    if( i ) {
      cout << ",";
    }
    cout << " " << table[ i ];
  }
  cout << " };\n";
}

// Returns true if str is a valid C identifier name
bool isValidName(char* str)
{
    size_t len = strlen(str);
    if (len == 0)
        return false;
    if (!isalpha(str[0]) && str[0] != '_')
        return false;
    for (size_t i = 1; i < len; ++i)
        if (!isalpha(str[i]) && !isdigit(str[i]) && str[i]!='_')
            return false;
    return true;
}

int main( int argc, char **argv )
{
  Game *game;
  size_t fwdHistoryLen, bwdHistoryLen, ruleGroupSize, verbosity;
  bool bwdMoves, removeDuplicateRules, dynamicAbstractions, useAbstraction, useStateMap ;
  RuleTree *ruleTree, *bwdRuleTree = NULL;

  ruleGroupSize = ruleGroupSizeDEFAULT;
  fwdHistoryLen = fwdHistoryLenDEFAULT;
  bwdHistoryLen = bwdHistoryLenDEFAULT;
  bwdMoves = bwdMovesDEFAULT;
  removeDuplicateRules = removeDuplicateRulesDEFAULT;
  dynamicAbstractions = dynamicAbstractionsDEFAULT;
  useAbstraction = AbstractionDEFAULT;
  useStateMap = state_mapDEFAULT;
  verbosity = verbosityDEFAULT;

//  string stateSpaceName("psvn_state_space");


/* PROCESS THE COMMAND LINE ARGUMENTS */

  for( int i = 1; i < argc; ++i ) {

    if(  !strcmp( argv[ i ], "-h" ) || !strcmp( argv[ i ], "--help" ) ) {
      cout << "usage: psvn2c [options]" << endl;
      cout << "The PSVN state space definition is read from stdin" << endl;
      cout << "and the generated code is printed on stdout." << endl;
      cout << "Command line options can be one or more of:" << endl;
      cout << "--help or -h: print this message" << endl;
      cout << "--backwards_moves: generate backwards move code";
      if ( bwdMovesDEFAULT==true ) { cout << " (default)"; }
      cout << endl;
      cout << "--no_backwards_moves: don't generate backwards move code";
      if ( bwdMovesDEFAULT==false ) { cout << " (default)"; }
      cout << endl;
      cout << "--history_len=len: use move histories of length len for pruning forward and backwards moves";
      cout << " (0 <= len <= " << historyLenMAX << ")";
      if ( fwdHistoryLenDEFAULT==bwdHistoryLenDEFAULT ) { cout << " (default len=" << fwdHistoryLenDEFAULT << ")"; }
      cout << endl;
      cout << "--fwd_history_len=len: use move histories of length len for pruning forward moves";
      cout << " (0 <= len <= " << historyLenMAX << ")";
      cout << " (default len=" << fwdHistoryLenDEFAULT << ")";
      cout << endl;
      cout << "--bwd_history_len=len: use move histories of length len for pruning backwards moves";
      cout << " (0 <= len <= " << historyLenMAX << ")";
      cout << " (default len=" << bwdHistoryLenDEFAULT << ")";
      cout << endl;
      cout << "--rule_group_size=num: consider num rules at a time";
      cout << " (0 <= num <= " << ruleGroupSizeMAX << "). num=0 considers all rules together.";
      cout << " (default num=" << ruleGroupSizeDEFAULT << ")";
      cout << endl;
      cout << "--keep_duplicate_rules: keep rules with identical preconditions and effects" ;
      if (  removeDuplicateRulesDEFAULT==false ) { cout << " (default)"; }
      cout << endl;
      cout << "--remove_duplicate_rules: remove rules with identical preconditions and effects" ;
      if (  removeDuplicateRulesDEFAULT==true ) { cout << " (default)"; }
      cout << endl;
      cout << "--abstraction: generate the code needed for state space abstraction" ;
      if ( AbstractionDEFAULT==true ) { cout << " (default)"; }
      cout << endl;
      cout << "--no_abstraction: don't generate the code for state space abstraction" ;
      if ( AbstractionDEFAULT==false ) { cout << " (default)"; }
      cout << endl;
      cout << "--state_map: generate the code for the state_map data structure" ;
      if ( state_mapDEFAULT==true ) { cout << " (default)"; }
      cout << endl;
      cout << "--no_state_map: don't generate the state_map code" ;
      if ( state_mapDEFAULT==false ) { cout << " (default)"; }
      cout << endl;
      cout << "--verbosity=val: set verbosity level to val";
      cout << " (0 <= val <= " << verbosityMAX << "). val=0 produces minimal output.";
      cout << " (default val=" << verbosityDEFAULT << ")";
      cout << endl;
      exit( 0 );

    } else if( !strncmp( argv[ i ], "--state_map", 11 ) ) {
              useStateMap = true;
    } else if( !strncmp( argv[ i ], "--no_state_map", 14 ) ) {
              useStateMap = false;
    } else if( !strncmp( argv[ i ], "--abstraction", 13 ) ) {
              useAbstraction = true;
    } else if( !strncmp( argv[ i ], "--no_abstraction", 16 ) ) {
              useAbstraction = false;
    } else if( !strncmp( argv[ i ], "--backwards_moves", 17 ) ) {
              bwdMoves = true;
    } else if( !strncmp( argv[ i ], "--no_backwards_moves", 20 ) ) {
              bwdMoves = false;
    } else if( !strncmp( argv[ i ], "--history_len=", 14 ) ) {
              if( !sscanf( &argv[ i ][ 14 ], "%zu", &fwdHistoryLen )
	          || fwdHistoryLen > historyLenMAX ) {
	        cerr << "bad history length argument:" << argv[ i ] << endl;
	        exit( -1 );
              }
              bwdHistoryLen = fwdHistoryLen;
    } else if ( !strncmp( argv[ i ], "--fwd_history_len=", 18 ) ) {
              if( !sscanf( &argv[ i ][ 18 ], "%zu", &fwdHistoryLen )
	          || fwdHistoryLen > historyLenMAX ) {
	        cerr << "bad forward history length argument:" << argv[ i ] << endl;
	        exit( -1 );
              }
    } else if( !strncmp( argv[ i ], "--bwd_history_len=", 18 ) ) {
              if( !sscanf( &argv[ i ][ 18 ], "%zu", &bwdHistoryLen )
	          || bwdHistoryLen > historyLenMAX ) {
	        cerr << "bad backward history length argument:" << argv[ i ] << endl;
	        exit( -1 );
              }
    } else if( !strncmp( argv[ i ], "--rule_group_size=", 18 ) ) {
              if( !sscanf( &argv[ i ][ 18 ], "%zu", &ruleGroupSize )
	          || ruleGroupSize > ruleGroupSizeMAX ) {
	        cerr << "bad rule group size argument:" << argv[ i ] << endl;
	        exit( -1 );
              }
    } else if( !strncmp( argv[ i ], "--remove_duplicate_rules", 24 ) ) {
              removeDuplicateRules = true;
    } else if( !strncmp( argv[ i ], "--keep_duplicate_rules", 24 ) ) {
              removeDuplicateRules = false;
    } else if( !strncmp( argv[ i ], "--verbosity=", 12 ) ) {
              if( !sscanf( &argv[ i ][ 12 ], "%zu", &verbosity )
	          || verbosity > verbosityMAX ) {
	        cerr << "bad verbosity argument:" << argv[ i ] << endl;
	        exit( -1 );
              }
    } 
  }

cerr << "psvn2c command line options in effect:" << endl ;
cerr << "Backwards rules will " ;
if ( !bwdMoves ) { cerr << "not "; }
cerr << "be generated." << endl;

if( fwdHistoryLen == 0 ) {
   cerr << "No move pruning code will be generated for forward search." << endl;
} else { cerr << "Move pruning for forward search will use a history length of " << fwdHistoryLen << "." << endl ; }

if( bwdMoves ) {
   if( bwdHistoryLen == 0 ) {
      cerr << "No move pruning code will be generated for backwards search." << endl;
   } else { cerr << "Move pruning for backwards search will use a history length of " << bwdHistoryLen << "." << endl ; }
}

cerr << "The initial size of each rule group will be " ;
if ( ruleGroupSize == 0 ) {
   cerr << "all rules." << endl ;
} else {
   cerr << ruleGroupSize << "." << endl ;
}

cerr << "Duplicate rules will be ";
if ( removeDuplicateRules ) {
   cerr << "removed." << endl ;
} else {
   cerr << "kept." << endl ;
}

cerr << "Code supporting abstraction will " ;
if ( !useAbstraction ) { cerr << "not "; }
cerr << "be generated." << endl;

cerr << "Code for the state_map data structure will " ;
if ( !useStateMap ) { cerr << "not "; }
cerr << "be generated." << endl;

cerr << "Verbosity is set to " << verbosity << "." << endl;


/* READ and INTERNALIZE THE PSVN */

if ( verbosity >= 2 ) { cerr << "(main) Read the PSVN file." << endl ; }

  PSVN *psvn = PSVN::fromInput( cin );
  if( psvn == NULL ) { exit( -1 ); }

  if( (*psvn).numGoalRules() == 0 ) {
    cerr << "a PSVN file must have at least one GOAL condition" << endl;
    exit( -1 );
  }


/* GENERATE ALL THE INTERNALS NEEDED FOR THE API */

if ( verbosity >= 2 ) { cerr << "(main) Convert the PSVN to internal Game format." << endl ; }
  game = Game::fromPSVN( *psvn,
			 bwdMoves,
			 removeDuplicateRules,
			 verbosity );
  if( game == NULL ) { exit( -1 ); }
  delete psvn;

if ( verbosity >= 1 ) { cerr << "(main) There are " << game->rules().size() << " forward rules in the internal format." << endl; }


/* GENERATE THE C code for the parts of the API that are not in the standard (non-custom) files. */

if ( verbosity >= 2 ) { cerr << "(main) Begin generating the API code." << endl ; }

  cout << "#include <stdlib.h>\n";
  cout << "#include <stdio.h>\n";
  cout << "#include <stdint.h>\n";
  cout << "#define __STDC_FORMAT_MACROS\n";
  cout << "#define __STDC_LIMIT_MACROS\n";
  cout << "#include <inttypes.h>\n";
  cout << "#include <string.h>\n";
  cout << "#include <assert.h>\n";
  cout << "#include <ctype.h>\n";
  cout << "\n\n";
  
  cout << "#define psvn2c_PSVN_API\n\n";


  if( fwdHistoryLen ) {
    cout << "#define HAVE_FWD_MOVE_PRUNING\n" ;
  }
  if( bwdMoves ) {
    cout << "#define HAVE_BWD_MOVES\n" ;
    if( bwdHistoryLen ) {
      cout << "#define HAVE_BWD_MOVE_PRUNING\n" ;
    }
  }
  if( useAbstraction ) {
     cout << "#define HAVE_ABSTRACTION\n" ;
  }
  cout << "\n\n";

  cout << "/* number of variables in a state */\n";
  cout << "#define NUMVARS " << game->numVars() << "\n\n\n";


/* choose the number of bits needed to represent domain values
   and the corresponding i/o formatting codes.
*/
  size_t maxDomainSize =  game->maxDomainSize();
  if( maxDomainSize >= 2147483648LL ) {
    cout << "typedef int64_t var_t;\n";
    cout << "#define PRI_VAR PRId64\n";
    cout << "#define SCN_VAR SCNd64\n";
  } else if( maxDomainSize >= 32768 ) {
    cout << "typedef int32_t var_t;\n";
    cout << "#define PRI_VAR PRId32\n";
    cout << "#define SCN_VAR SCNd32\n";
  } else if( maxDomainSize >= 128 ) {
    cout << "typedef int16_t var_t;\n";
    cout << "#define PRI_VAR PRId16\n";
    cout << "#define SCN_VAR SCNd16\n";
  } else {
    cout << "typedef int8_t var_t;\n";
    cout << "#define PRI_VAR PRId8\n";
    cout << "#define SCN_VAR SCNd8\n";
  }
  cout << "\n";


if ( verbosity >= 2 ) { cerr << "(main) Generating the Domain code." << endl ; }
  game->generateDomainCode( cout, useAbstraction );
  cout << "\n";

  cout << "typedef struct {\n"
       << " var_t vars[ NUMVARS ];\n"
       << "} state_t;\n"
       << "\n" ;

  cout << "typedef int (*funcptr)( const state_t *, void * );\n";
  cout << "typedef void (*actfuncptr)( const state_t *, state_t * );\n";
  cout << "\n";

  cout << "typedef struct {\n"
       << " const state_t *my_state;\n"
       << " funcptr my_funcptr;\n"
       << "} ruleid_iterator_t;\n"
       << "\n" ;


  cout << "#define num_fwd_rules " << game->rules().size() << '\n';
  if( bwdMoves ) {
     cout << "#define num_bwd_rules " << game->backwardsRules().size() << '\n';
  }

  game->generateRuleNames( cout, bwdMoves );

  game->generateRuleCosts( cout, bwdMoves );

  if ( useAbstraction ) {
     game->generateRuleLabelSets( cout, bwdMoves );  // only needed for abstraction
  }


  if( fwdHistoryLen ) {
/* This is where the move pruning is done for the forward rules */
if ( verbosity >= 1 ) { cerr << "(main) Perform move pruning analysis (forward rules)." << endl ; }
      PruneTree *fwdPruneTree = new PruneTree( *game, game->rules(),
                                               fwdHistoryLen, verbosity );
      vector<int> fwdPruneTable = fwdPruneTree->makeTable();
      fwdPruneTree->destroyTree();
      cout << "\n";
      printPruneTable( fwdPruneTable, string( "fwd" ) );
  }
if ( verbosity >= 2 ) { cerr << "(main) Create the code for the effects of the forward rules." << endl ; }
  printRuleActions( "fwd", game->rules(), game->varDomains(), 
                    dynamicAbstractions );


  int fwdVarEntryId;
  int bwdVarEntryId;
  size_t numFwdNodes;
if ( verbosity >= 2 ) { cerr << "(main) Create the code for the preconditions of the forward rules." << endl ; }
  ruleTree = new RuleTree( *game, game->rules(), ruleGroupSize, &numFwdNodes, verbosity );
  ruleTree->generateFunctionCode( cout, string( "fwd" ), numFwdNodes,
                                  &fwdVarEntryId, dynamicAbstractions );


  if( bwdMoves ) {

if ( verbosity >= 1 ) { cerr << "(main) There are " << game->backwardsRules().size() << " backwards rules in the internal format." << endl; }

    if( bwdHistoryLen ) {
/* This is where the move pruning is done for the backward rules */
if ( verbosity >= 1 ) { cerr << "(main) Perform move pruning analysis (backwards rules)." << endl ; }
        PruneTree *bwdPruneTree = new PruneTree( *game, game->backwardsRules(),
                                                 bwdHistoryLen, verbosity );
        vector<int> bwdPruneTable = bwdPruneTree->makeTable();
        bwdPruneTree->destroyTree();
        cout << "\n";
        printPruneTable( bwdPruneTable, string( "bwd" ) );
    }

if ( verbosity >= 2 ) { cerr << "(main) Create the code for the effects of the backwards rules." << endl ; }
    printRuleActions( "bwd", game->backwardsRules(), game->varDomains(),
                      dynamicAbstractions );

if ( verbosity >= 2 ) { cerr << "(main) Create the code for the preconditions of the backwards rules." << endl ; }
    size_t numBwdNodes;
    bwdRuleTree = new RuleTree( *game, game->backwardsRules(), ruleGroupSize, 
                                &numBwdNodes, verbosity );
    bwdRuleTree->generateFunctionCode( cout, string( "bwd" ), numBwdNodes,
                                       &bwdVarEntryId, dynamicAbstractions);
  }

if ( verbosity >= 2 ) { cerr << "(main) Create the code related to iterators." << endl ; }
  cout << "\n\n";
  cout << "#define init_history 0\n\n";
  cout << "#define max_fwd_children " << ruleTree->maxChildren() << "\n\n";

  cout << "/* NOTE: FOR ALL OF THE MOVE ITERATOR DEFINITIONS funcptr\n";
  cout << "   MUST BE A VARIABLE. */\n\n";
  cout << "/* initialise a forward move iterator */\n";
  cout << "#define init_fwd_iter( ruleid_iter, state ) { \\\n" ;
  cout << "  (*ruleid_iter).my_funcptr = " ;
  if ( game->rules().size() > 0 ) {
       cout << ruleTree->entryName( string( "fwd" ) ) << " ; \\\n";
  } else {
       cout << "NULL ;\\\n";
  }
  cout << "  (*ruleid_iter).my_state = state; \\\n}\n\n";


  cout << "/* use iterator to generate next applicable rule to apply \n";
  cout << "   returns rule to use, -1 if there are no more rules to apply */\n";
  cout << "#define next_ruleid( ruleid_iter ) (((*ruleid_iter).my_funcptr)?((*ruleid_iter).my_funcptr)((*ruleid_iter).my_state,&((*ruleid_iter).my_funcptr)):-1)\n\n";


  cout << "/* apply a rule to a state */\n";
  cout << "#define apply_fwd_rule( rule, state, result ) fwd_rules[(rule)](state,result)\n";

  cout << "/* returns 0 if the rule is pruned, non-zero otherwise */\n";
  if( fwdHistoryLen ) {
      cout << "#define fwd_rule_valid_for_history( history, rule_used ) (fwd_prune_table[(history)+(rule_used)])\n";
      cout << "/* generate the next history from the current history and a rule */\n";
      cout << "#define next_fwd_history( history, rule_used ) (fwd_prune_table[(history)+(rule_used)])\n\n";
      cout << "\n";
  } else {
      cout << "#define fwd_rule_valid_for_history( history, rule_used ) 1 \n";
      cout << "/* generate the next history from the current history and a rule */\n";
      cout << "#define next_fwd_history( history, rule_used ) 0 \n\n";
      cout << "\n";
  }

  if( bwdMoves ) {

    if ( verbosity >= 2 ) { cerr << "(main) about to generate code related to backwards moves.\n"; }

    cout << "static const int bw_max_children = " 
         << bwdRuleTree->maxChildren() << ";\n";
    cout << "#define max_bwd_children " << bwdRuleTree->maxChildren() << "\n\n";

    cout << "/* initialise a backwards move iterator */\n";

    cout << "#define init_bwd_iter( ruleid_iter, state ) { \\\n" ;
    cout << "  (*ruleid_iter).my_funcptr = " ;
    if ( game->backwardsRules().size() > 0 ) {
         cout << bwdRuleTree->entryName( string( "bwd" ) ) << " ; \\\n";
    } else {
         cout << "NULL ;\\\n";
    }
    cout << "  (*ruleid_iter).my_state = state; \\\n}\n\n";

    cout << "/* apply a rule to a state */\n";
    cout << "#define apply_bwd_rule( rule, state, result ) bwd_rules[(rule)](state,result)\n";

    cout << "/* returns 0 if the rule is pruned, non-zero otherwise */\n";
    if( bwdHistoryLen ) {
        cout << "#define bwd_rule_valid_for_history( history, rule_used ) (bwd_prune_table[(history)+(rule_used)])\n";
        cout << "/* generate the next history from the current history and a rule */\n";
        cout << "#define next_bwd_history( history, rule_used ) (bwd_prune_table[(history)+(rule_used)])\n\n";
        cout << "\n";
    } else {
        cout << "#define bwd_rule_valid_for_history( history, rule_used ) 1 \n";
        cout << "/* generate the next history from the current history and a rule */\n";
        cout << "#define next_bwd_history( history, rule_used ) 0 \n\n";
        cout << "\n";
    }
  }


if ( verbosity >= 2 ) { cerr << "(main) Create the code related to the goal conditions." << endl ; }

  cout << game->generateGoalCode( bwdMoves ) << "\n";



/* append a copy of the appropriate files */
if ( verbosity >= 2 ) { cerr << "(main) Append non-custom code." << endl ; }

  fstream fs; char line[ 1024 ];

/* CORE */
  fs.open( "psvn2c_core.c", fstream::in );
  while( fs.getline( line, 1024 ) ) {
    cout << line << endl;
  } 
  fs.close();

/* STATE_MAP */
  if ( useStateMap ) {
     fs.open( "psvn2c_state_map.c", fstream::in );
     while( fs.getline( line, 1024 ) ) {
       cout << line << endl;
     } 
     fs.close();
  } 

/* ABSTRACTION */
  if ( useAbstraction ) {
     fs.open( "psvn2c_abstraction.c", fstream::in );
     while( fs.getline( line, 1024 ) ) {
       cout << line << endl;
     } 
     fs.close();
  } 



/* CLEAN UP */
  ruleTree->destroyTree();
  if (bwdMoves) { bwdRuleTree->destroyTree(); }
  delete game;

if ( verbosity >= 2 ) { cerr << "(main) Done." << endl ; }
}
