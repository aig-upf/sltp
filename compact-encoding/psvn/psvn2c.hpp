/*
Copyright (C) 2011 by the PSVN Research Group, University of Alberta
*/

#ifndef _GAME_H
#define _GAME_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "psvn.hpp"

using namespace std;


class Game;
class Context;
class VarTest;

enum ValueType { Constant, Label };

class Value {
protected:
  ValueType m_type;
  size_t m_value;

public:
  Value( ValueType type, size_t value )
    : m_type( type ), m_value( value ) {}

  bool operator==( const Value &v ) const {
    if( m_type == v.m_type && m_value == v.m_value ) { return true; }
    return false;
  }

  const ValueType type() const { return m_type; }
  const size_t value() const { return m_value; }

  string toString() const {
    string s; char t[ 32 ];
    if( m_type == Label ) { s.append( "[" ); }
    sprintf( t, "%zu", m_value ); s.append( t );
    if( m_type == Label ) { s.append( "]" ); }
    return s;
  }
};


class Rule {
protected:
  size_t m_id;
  string m_name;
  size_t m_cost;
  vector<vector<Value> > m_tests;
  vector<vector<Value> > m_actions;
  bool m_modified;	// used in RuleTrees.  A test taken resulted
			// in the rule being modified (e.g. test removed)

public:
  Rule( const size_t id, const string name, const size_t cost,
	const vector<vector<Value> > &tests,
	const vector<vector<Value> > &actions )
    : m_id( id ),
      m_name( name ),
      m_cost ( cost ), 
      m_tests( tests ),
      m_actions( actions ),
      m_modified( false ) {
    assert( tests.size() == actions.size() );
  }

  bool operator==( const Rule &rule ) const;

  bool isEmpty() const;

  const size_t id() const { return m_id; }
  const string name() const { return m_name; }
  const size_t cost() const { return m_cost; }
  const vector<vector<Value> > &tests() const { return m_tests; }
  const vector<vector<Value> > &actions() const { return m_actions; }

  /* only used in RuleTrees */
  const bool wasModified() const { return m_modified; }

  void setId( const size_t id ) { m_id = id; }

  string toString() const {
    string s; size_t i, j;
    for( i = 0; i < m_tests.size(); ++i ) {
      if( i ) { s.append( " " ); }
      for( j = 0; j < m_tests[ i ].size(); ++j ) {
	if( j ) { s.append( "," ); }
	s.append( m_tests[ i ][ j ].toString() );
      }
      if( m_tests[ i ].size() == 0 ) { s.append( "-" ); }
    }
    s.append( " => " );
    for( i = 0; i < m_actions.size(); ++i ) {
      if( i ) { s.append( " " ); }
      if( m_actions[ i ].size() == 0 ) { s.append( "-" ); }
      else { s.append( m_actions[ i ][ 0 ].toString() ); }
    }
    return s;
  }

  void print() const;

  /* rewrite rule with knowledge that [x] = [y]
     removes references to x */
  void useEquivalence( size_t x, size_t y );

  /* use proven/falsified knowledge in RuleTree, removing any proven tests
     returns false if rule is provably unusable, true otherwise
     (ie true means either proven or plausible)
     NOTE: even if rule is unusable, some true tests may still be removed */
  bool useKnowledge( const Context &context );

  /* true if there are no tests left, false otherwise
     if this is called after useKnowledge, there will be no
     tests if and only if the rule is provable */
  bool hasNoTests() const {
    for( size_t i = 0; i < tests().size(); ++i ) {
      if( tests()[ i ].size() ) { return false; }
    }
    return true;
  }

  /* give a list of action variables which aren't bound by the tests
     each element is <varLabel, varIdx> pair, where varLabel is the
     (integer) label of the unbound variable, and varIdx is an
     index into the actions where it occurs (note that the variable
     may occur at more than one index) */
  vector<pair<size_t,size_t> > unboundActionVariables() const {
    vector<pair<size_t,size_t> > vars;
    size_t i, j;
    bool used[ m_actions.size() ];

    for( i = 0; i < m_actions.size(); ++i ) { used[ i ] = false; }
    for( i = 0; i < m_actions.size(); ++i ) {
      if( m_actions[ i ].size() && m_actions[ i ][ 0 ].type() == Label ) {
	j = m_actions[ i ][ 0 ].value();
	if( j >= m_actions.size() ) {
	  assert( j - m_actions.size() < m_actions.size() );
	  if( !used[ j - m_actions.size() ] ) {
	    used[ j - m_actions.size() ] = true;
	    vars.push_back( pair<size_t,size_t>( j, i ) );
	  }
	}
      }
    }
    return vars;
  }

  /* create a rule with unbound variables replaced with constant values */
  Rule( const Rule &rule, const vector<pair<size_t,size_t> > &unboundVars,
	const vector<size_t> &unboundVals );

  /* returns true if knowledge of variable may affect whether or not
     the rule is usable.  returns false if variable is irrelevant */
  bool usesVariable( const size_t var ) const;

  /* generate code that performs the action of a rule */
  void generateActionCode( ostream &os,
			   const string &oldState,
			   const string &newState,
			   const size_t indentation ) const;

  void generateDynActionCode( ostream &os,
                              const string &oldState,
                              const string &newState,
                              const string &treeName,
                              const size_t indentation,
                              const vector<size_t>& varDomains) const;

  string actionFunctionName( const string &treeName ) const {
    string s = treeName;
    char cStr[ 64 ];
    sprintf( cStr, "rule%zu", id() );
    s.append( cStr );
    return s;
  }

  /* generate a function "treeName""rule"id() */
  void generateActionFunction( ostream &os,
			       const string &treeName,
                               const vector<size_t>& varDomains,
                               const bool dynamicAbstractions ) const;

  /* given a set of values for all variables (as in PruneTree)
     try applying the rule and modifying the values.
     adds any variable changes to satisfy preconditions to contexts 
     Returns false if the rule can't be applied.
     NOTE: VALUES AND CONTEXTS MAY STILL BE MODIFIED
     Returns true if the rule can be applied and applies rule to values */
  bool applyTo( vector<Value> &values,
		vector<Value> &context ) const;
};


class Context {
protected:
  const Game &m_game;

  /* m_proven[ var ].type() will be Label if nothing has been proven
     for var, because  var=y is handled by substitution */
  vector<Value> m_proven;
  vector<vector<Value> > m_falsified;

public:
  Context( const Game &game );
  Context( const Context &context,
	   const VarTest &test, const size_t result );

  bool operator==( const Context &context ) const;

  const vector<Value> &proven() const { return m_proven; }
  const vector<vector<Value> > &falsified() const { return m_falsified; }

  /* go through proven/falsified tests to remove duplicate knowledge */
  void removeDuplicates();

  bool testIsProvable( const size_t var, const Value &test,
		       bool *ret_isProven ) const;

  void setTestProven( const size_t var, const Value &test );
  void setTestFalsified( const size_t var, const Value &test );

  /* throw away all knowledge about variable */
  void clearVar( const size_t var ) {
    m_proven[ var ] = Value( Label, 0 );
    m_falsified[ var ].clear();
  }

  void print() const;
};

class VarTest {
protected:
  size_t m_var;
  ValueType m_type;
  size_t m_otherLabel; /* used if m_type == Label, test is var ?= otherLabel */

public:

  /* simple constant value test */
  VarTest( const size_t var )
    : m_var( var ), m_type( Constant ), m_otherLabel( 0 ) {}

  /* variable equality test */
  VarTest( const size_t var, const size_t otherLabel )
    : m_var( var ), m_type( Label ), m_otherLabel( otherLabel ) {}

  const size_t var() const { return m_var; }
  const ValueType type() const { return m_type; }

  /* only used if type == Labe*/
  const size_t otherLabel() const { return m_otherLabel; }

  void print() const {
    cout << "[" << m_var << "] == ";
    if( m_type == Label ) {
      cout << "[" << m_otherLabel << "]";
    } else {
      assert( m_type == Constant );
      cout << "?";
    }
    cout << endl;
  }
};

class RuleTree {
protected:
  const Game &m_game;

  size_t m_id;

  Context m_context;

  vector<Rule> m_rules; /* rules which we might be able to apply */

  /* rules we now know we can apply at root - only used at root */
  vector<Rule> m_rootRules;

  VarTest m_test;
  size_t m_testedRuleIdx; /* more than rule may be modified by test, but
			     this rule is the source of m_test */
  vector<RuleTree *> m_children;
  vector<bool> m_childIsDup;
  vector<vector<Rule> > m_childFinishedRules;

  size_t m_maxChildren;

  RuleTree( const Game &game,
	    const Context &context,
	    const vector<Rule> &rules );

  /* forget any information about variables that we no longer need */
  void forgetVariables();

  void setTest( const size_t ruleIdx, const size_t var, const size_t testIdx );

  /* pick a new test to make from the rule list
     sets m_testVar and m_testValue
     all provable rules must already have been removed
     returns true/false on succes/failure */
  bool pickTest();

  bool hasRule( const Rule &rule ) const;

  string baseFunctionName( const string &treeName ) const {
    string s( treeName );
    char cStr[ 64 ];
    sprintf( cStr, "fn%zu", id() );
    s.append( cStr );
    return s;
  }
  string rootActionFunctionName( const string &treeName,
				 const size_t i ) const {
    string s = baseFunctionName( treeName );
    char cStr[ 64 ]; sprintf( cStr, "_r%zu", i ); s.append( cStr );
    return s;
  }
  string actionStubName( const string &treeName,
			 const size_t result, const size_t i ) const {
    string s = baseFunctionName( treeName );
    char cStr[ 64 ];
    sprintf( cStr, "_a%zu_%zu", result, i ); s.append( cStr );
    return s;
  }
  string functionHeader( const string &name ) const {
    string s;

    s.append( "static int " );
    s.append( name );
    s.append( "( const state_t *state, void *next_func" );
    s.append( " )\n{\n" );
    return s;
  }
  string nextFuncCode( const string &func ) const {
    string s( "*((funcptr *)next_func) = " );
    s.append( func );
    s.append( ";\n" );
    return s;
  }
  void resultCall( ostream &os,
		   const string &treeName,
		   const size_t result ) const;
  void generateRootCode( ostream &os,
			 const string &treeName ) const;
  void generateActionStub( ostream &os,
			   const string &treeName,
			   const size_t result, const size_t i ) const;

  void generateActionStubCode( ostream &os,
			       const string &funcName,
			       const string &treeName,
			       const Rule &rule,
			       const string nextFuncName ) const;

  /* returns true if a and b lead to the same results */
  bool resultsSame( size_t a, size_t b ) const;

  /* returns 0 if there is no test needed,
     returns 1 if trueValue leads to one node, and all other results
     lead to the same node as falseValue,
     returns 2 in all other cases (ie a switch statement is needed) */
  size_t getResultType( size_t *trueValue, size_t *falseValue ) const;

public:
  /* build tree of rules
     If ruleGroupSize is non-zero, Try grouping rules into disjoint sets
     with ruleGroupSize rules (for speeding up analysis)
     If groupSize is 0, use all rules at once */
  RuleTree( const Game &game,
	    const vector<Rule> &rules,
	    const size_t ruleGroupSize,
            size_t* numNodesInTree,
            size_t verbosity);

  /* destroy a tree, and all children */
  void destroyTree();

  /* NOTE: both 'this' and 'ruleTree' must be sorted according to the
     same total ordering (for example, the order in the original
     list of the original list of rules) */
  bool operator==( const RuleTree &ruleTree ) const;

  bool buildChildren( vector<vector<RuleTree *> > &trees,
		      size_t *numNodes,
		      RuleTree *leafTree,
		      bool *leafTreeUsed,
		      const size_t maxNodes = 0 );

  const size_t id() const { return m_id; }
  const Context &context() const { return m_context; }
  const vector<Rule> &rules() const { return m_rules; }
  const VarTest test() const { return m_test; }
  const size_t testedRuleIdx() const { return m_testedRuleIdx; }
  const size_t maxChildren() const { return m_maxChildren; }

  void print() const;

  void generateFunctionCode( ostream &os,
			     const string &treeName,
                             size_t numNodesInTree,
                             int* varEntryId,
                             const bool dynamicAbstractions ) const;
  string entryName( const string &treeName ) const {
    if( m_rootRules.size() ) {
      return rootActionFunctionName( treeName, 0 );
    } else {
      return baseFunctionName( treeName );
    }
  }
};


class PruneTree {
protected:
  /* pointers for tree structure */
  PruneTree *m_parent;			/* NULL at root */
  const size_t m_parentRuleUsed;	/* 0 at root */
  PruneTree *m_suffix;
  vector<PruneTree *> m_children;

  /* identifier - also used as index into allNodes in expand() */
  size_t m_id;			/* 0 at root */

  /* cost of taking the sequence of moves to reach this node */
  size_t m_cost;

  /* values of state variables - values >= domain size are variables
     with unknown values corresponding to an arbitrary initial state */
  vector<Value> m_values;

  /* the preconditions that had to be satisfied in order to use the
     rules which reached this node in the tree */
  vector<Value> m_context;

  /* create a new tree node with the given parent, suffix, id,
     values, and contexts. */
  PruneTree( PruneTree *parent,
	     const size_t ruleUsed,
	     PruneTree *suffix,
	     const size_t id,
	     const size_t cost,
	     const vector<Value> &values,
	     const vector<Value> &context )
    : m_parent( parent ),
      m_parentRuleUsed( ruleUsed ),
      m_suffix( suffix ),
      m_id( id ),
      m_cost( cost ),
      m_values( values ),
      m_context( context ) {}

  /* destroy a tree we're currently building with expand */
  void destroyTree( vector<PruneTree *> &allNodes );

  /* generate all the children of this node */
  void expand( const vector<Rule> &rules,
	       vector<PruneTree *> &allNodes,
	       int verbosity );

  /* try expanding all children at the specified depth
     NOTE: THIS VALUE SHOULD BE SMALL, LIKE 2 or 3 BECAUSE THE
     RESULTING TREE IS LIKELY GOING TO GROW VERY QUICKLY */
  void expandTo( const vector<Rule> &rules,
		 const size_t depth,
		 vector<PruneTree *> &allNodes,
	         int verbosity );

  void removeFromTree( vector<PruneTree *> &allNodes );

  /* re-label the interior nodes with sequentially increasing IDs */
  void reLabel( size_t *nextLabel ) {
    if( m_children.size() == 0 ) {
      return;
    }
    m_id = *nextLabel;
    ++( *nextLabel );
    for( size_t i = 0; i < m_children.size(); ++i ) {
      if( m_children[ i ] != NULL ) {
	m_children[ i ]->reLabel( nextLabel );
      }
    }
  }

  /* set the entries in a pruning table */
  void fillTable( vector<int> &table );

public:
  /* create a tree root and expand to the specified depth
     sets numInteriorNodes to the number of interior (non-leaf) nodes */
  PruneTree( const Game &game,
	     const vector<Rule> &rules,
	     const size_t depth,
	     int verbosity );

  /* delete node, but not children */
  ~PruneTree() {}

  /* delete node and all children */
  void destroyTree() {
    for( size_t i = m_children.size(); i > 0; --i ) {
      if( m_children[ i - 1 ] == NULL ) {
	continue;
      }

      m_children[ i - 1 ]->destroyTree();
    }

    if( m_parent == NULL ) {
      /* this is the root of a tree, so m_suffix is a specal dummy
	 node which needs to be destroyed */

      delete m_suffix;
    }

    delete this;
  }

  const vector<Value> &values() const { return m_values; }
  const vector<Value> &context() const { return m_context; }

  /* print out the pruning tree */
  void print( const size_t depth ) const;

  /* generate a table from the tree for doing move pruning
     will modify the IDs of the nodes */
  vector<int> makeTable();
};


class Game {
protected:
  vector<string> m_domainNames;
  vector<vector<string> > m_domains;
  vector<size_t> m_variableDomains;
  vector<string> m_ruleNames;
  vector<Rule> m_rules;
  vector<Rule> m_bwdRules;
  vector<vector<vector<Value> > > m_goals;

public:
  Game( const vector<string> domainNames,
        const vector<vector<string> > domains,
	const vector<size_t> variableDomains,
	const vector<Rule> rules,
	const vector<Rule> bwdRules,
	const vector<vector<vector<Value> > > goals )
    : m_domainNames(domainNames),
      m_domains( domains ),
      m_variableDomains( variableDomains ),
      m_rules( rules ),
      m_bwdRules( bwdRules ),
      m_goals( goals )
  {}

  static Game *fromPSVN( const PSVN &psvn,
			 const bool makeBackwardsMoves,
			 const bool removeDuplicateRules,
			 int verbosity );

  const size_t numVars() const { return m_variableDomains.size(); }
  const vector<size_t>& varDomains() const
  { return m_variableDomains; }
  const vector<string> &varDomain( const size_t var ) const
  { return m_domains[ m_variableDomains[ var ] ]; }
  const size_t domainSize( const size_t var ) const
  { return m_domains[ m_variableDomains[ var ] ].size(); }
  const size_t maxDomainSize() const
  {
    size_t m, i;
    for( m = 0, i = 0; i < m_variableDomains.size(); ++i ) {
      if( domainSize( i ) > m ) { m = domainSize( i ); }
    }
    return m;
  }

  void addRule( const Rule &rule ) { m_rules.push_back( rule ); }
  const vector<Rule> &rules() const { return m_rules; }

  void addBackwardsRule( const Rule &rule ) { m_bwdRules.push_back( rule ); }
  const vector<Rule> &backwardsRules() const { return m_bwdRules; }

  void generateDomainCode( ostream &os, const bool useAbstraction ) const;
  void generateRuleNames( ostream &os, const bool bwdMoves ) const;
  void generateRuleCosts( ostream &os, const bool bwdMoves ) const;
  void generateRuleLabelSets( ostream &os, const bool bwdMoves ) const;
  string generateGoalCode( const bool bwdMoves ) const;
  string generateDynGoalCode() const;
};

struct VarTestNode
{
    int type;
    int var; 
    int other;
    int rule;
    vector<int> edges;
    VarTestNode() : type(-1), rule(-1) { };
};

#endif
