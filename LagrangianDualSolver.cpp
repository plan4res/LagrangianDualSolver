/*--------------------------------------------------------------------------*/
/*--------------------- File LagrangianDualSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the LagrangianDualSolver class, which implements the
 * CDASolver interface within the SMS++ framework for a "generic"
 * Lagrangian-based Solver.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2020 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "LagrangianDualSolver.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "RBlockConfig.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Bitwise-coded macro that activate costly checks that should never be done
  * in production, but can be useful during debugging. Currently supported
  * checks are:
  *
  * - bit 0 (+ 1): is_correct() is called on the Lagrangian Dual
  *   AbstractBlock to verify that all Variable and Constraint are properly
  *   linked. */
#else
 #define CHECK_DS 0
 // never change this
#endif


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

using VarValue = Function::FunctionValue;
using c_VarValue = Function::c_FunctionValue;
using Vec_FunctionValue = Function::Vec_FunctionValue;

using Vec_VarValue = Function::Vec_FunctionValue;
using c_Vec_VarValue = Function::c_Vec_FunctionValue;

using coeff_pair = LinearFunction::coeff_pair;
using v_coeff_pair = LinearFunction::v_coeff_pair;
using v_c_coeff_pair = LinearFunction::v_c_coeff_pair;

using LinearCombination = C05Function::LinearCombination;
using c_LinearCombination = C05Function::c_LinearCombination;

using dual_pair = LagBFunction::dual_pair;
using v_dual_pair = std::vector< dual_pair >;
using v_c_dual_pair = const v_dual_pair;

using p_AB = AbstractBlock *;
using p_FRC = FRowConstraint *;
using p_FRO = FRealObjective *;
using p_LF = LinearFunction *;
using p_LBF = LagBFunction *;

using p_BC = BlockConfig *;
using p_BSC = BlockSolverConfig *;

using p_SConf_p_p = SimpleConfiguration< std::pair< Configuration * ,
						    Configuration * > > *;

/*--------------------------------------------------------------------------*/
/* LagrangianDualSolver always need to have an "inner CDASolver" set, but at
 * the beginning it would have none: it will then have a FakeCDASolver one,
 * which does nothing. Clearly, a "real" one will have to be set of
 * LagrangianDualSolver::compute() is to work. */

class FakeCDASolver : public CDASolver {
 public:
 FakeCDASolver( void ) : CDASolver() {}
 virtual ~FakeCDASolver() {}
 int compute( bool changedvars = true ) override final { return( kError ); }
 bool has_var_solution( void ) override final { return( false ); }
 void get_var_solution( Configuration *solc = nullptr ) override final {}
 bool has_dual_solution( void ) override final { return( false ); }
 void get_dual_solution( Configuration *solc = nullptr ) override final {}
 void add_Modification( sp_Mod & mod ) override final {}
 SMSpp_insert_in_factory_h;
 };

SMSpp_insert_in_factory_cpp_0( FakeCDASolver );

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr VarValue NaNshift
                              = std::numeric_limits< VarValue >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

static constexpr VarValue INFshift = Inf< VarValue >();
 ///< convenience constexpr for "Infty"

static constexpr Index InINF = SMSpp_di_unipi_it::Inf< Index >();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register LagrangianDualSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( LagrangianDualSolver );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names

const std::vector< std::string > LagrangianDualSolver::int_pars_str = {
 "int_LDSlv_iBCopy" , "int_LDSlv_NNMult" , "int_LDSlv_CloneCfg"
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of double parameters names

// const std::vector< std::string > LagrangianDualSolver::dbl_pars_str = { };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of string parameters names

const std::vector< std::string > LagrangianDualSolver::str_pars_str = {
 "str_LDSlv_ISName" , "str_LagBF_BCfg" , "str_LagBF_BSCfg" ,
 "str_LDBlck_BCfg" , "str_LDBlck_BSCfg"
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of vector-of-int parameters names

const std::vector< std::string > LagrangianDualSolver::vint_pars_str = {
 "vint_LDSl_WBCfg"  , "vint_LDSl_W2BCfg" ,
 "vint_LDSl_WBSCfg" , "vint_LDSl_W2BSCfg"
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the vector of vector-of-string parameters names

// const std::vector< std::string > LagrangianDualSolver::vstr_pars_str = {
//  "vstr_LDSl_Cfg"
//  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for int parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::int_pars_map = {
 { "int_LDSlv_iBCopy"   , LagrangianDualSolver::int_LDSlv_iBCopy } ,
 { "int_LDSlv_NNMult"   , LagrangianDualSolver::int_LDSlv_NNMult } ,
 { "int_LDSlv_CloneCfg" , LagrangianDualSolver::int_LDSlv_CloneCfg }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for double parameters names

// const std::map< std::string , LagrangianDualSolver::idx_type >
// LagrangianDualSolver::dbl_pars_map = { };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for string parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::str_pars_map = {
 { "str_LDSlv_ISName" , LagrangianDualSolver::str_LDSlv_ISName } ,
 { "str_LagBF_BCfg"   , LagrangianDualSolver::str_LagBF_BCfg } ,
 { "str_LagBF_BSCfg"  , LagrangianDualSolver::str_LagBF_BSCfg } ,
 { "str_LDBlck_BCfg"  , LagrangianDualSolver::str_LDBlck_BCfg } ,
 { "str_LDBlck_BSCfg" , LagrangianDualSolver::str_LDBlck_BSCfg }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for vector-of-int parameters names

const std::map< std::string , LagrangianDualSolver::idx_type >
 LagrangianDualSolver::vint_pars_map = {
 { "vint_LDSl_WBCfg"   , LagrangianDualSolver::vint_LDSl_WBCfg } ,
 { "vint_LDSl_W2BCfg"  , LagrangianDualSolver::vint_LDSl_W2BCfg } ,
 { "vint_LDSl_WBSCfg"  , LagrangianDualSolver::vint_LDSl_WBSCfg } ,
 { "vint_LDSl_W2BSCfg" , LagrangianDualSolver::vint_LDSl_W2BSCfg }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the map for vector-of-string parameters names

// const std::map< std::string , LagrangianDualSolver::idx_type >
//  LagrangianDualSolver::vstr_pars_map = {
//  { "vstr_LDSl_Cfg"  , LagrangianDualSolver::vstr_LDSl_Cfg }
//  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default int parameters

const std::vector< int > LagrangianDualSolver::dflt_int_par = {
 0 ,  // int_LDSlv_iBCopy
 1 ,  // int_LDSlv_NNMult
 0    // int_LDSlv_CloneCfg
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default double parameters

// const std::vector< double > LagrangianDualSolver::dflt_dbl_par = { };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// define and initialize here the default string parameters

const std::vector< std::string > LagrangianDualSolver::dflt_str_par = {
 "FakeCDASolver" ,  // str_LDSlv_ISName
 "" ,               // str_LagBF_BCfg
 "" ,               // str_LagBF_BSCfg
 "" ,               // str_LDBlck_BCfg
 ""                 // str_LDBlck_BSCfg
 };

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF LagrangianDualSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_Block( Block * block )
{
 if( f_Block ) {  // changing from a previous Block- - - - - - - - - - - - - -
                  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cleanup_LagrDual();
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 // lock the Block - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

 // create the Lagrangian Dual and its sub-Block- - - - - - - - - - - - - - - 
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // start right away to construct the Lagrangian Dual Block, even before
 // having fully checked that f_Block satisfies all the necessary conditions;
 // this is made so because checking the conditions requires looking at the
 // abstract representation, which may need to be generated for this very
 // purpose, but the generation of the abstract representation may differ
 // according to the BlockConfig, so ensure that all BlockConfig (but *not*
 // the BlockSolverConfig, see later on for why) that can be apply()-ed to
 // the sub-Block are before doing the checks
 // but at the very least children are required to exist - - - - - - - - - -

 const auto & sb = f_Block->get_nested_Blocks();
 f_nsb = sb.size();
 if( ! f_nsb )
  throw( std::invalid_argument( "LagrangianDualSolver: no sub-Block" ) );

 // create the default BlockConfig
 if( ! LagBF_BCfg.empty() ) {
  auto c = Configuration::deserialize( LagBF_BCfg );
  f_DBCfg = dynamic_cast< BlockConfig * >( c );
  if( ! f_DBCfg ) {
   delete c;
   throw( std::invalid_argument(
		  "LagrangianDualSolver: LagBF_BCfg not a BlockConfig" ) );

   }
  }

 LagrDual = new AbstractBlock;  // create the AbstractBlock

 // resize the sub-Block dictionary and the pointers to the LagBFunction
 blck_to_idx.resize( f_nsb );
 v_LBF.resize( f_nsb );

 // first loop: create the sub-Block and their LagBFunction - - - - - - - - -
 // here comes the crucial decision: copy the sub-Block or "evict" them

 Index iW2BCfg = 0;  // index in W2BCfg
 
 for( Index i = 0 ; i < f_nsb ; ++i ) {
  auto sbi = new AbstractBlock;
  blck_to_idx[ i ].first = sbi;
  blck_to_idx[ i ].second = i;

  Block * csbi;
  if( iBCopy )  // copying the sub-Block
   csbi = sb[ i ]->get_R3_Block( nullptr );  // the copy R3B
  else         // evicting the sub-Block
   csbi = sb[ i ];  // use the original sub-Block; note that its father
                    // will be changed when used in LagBFunction constructor

  // BlockConfig-ure the inner Block
  auto BCi = f_DBCfg;
  if( ! WBCfg.empty() ) {  // individual BlockConfig are provided
   Index h;                // the index in WBCfg

   if( W2BCfg.empty() )    // in dense format
    h = i;
   else                    // in sparse format
    if( ( iW2BCfg < W2BCfg.size() ) && ( W2BCfg[ iW2BCfg ] == i ) )
     h = iW2BCfg++;
    else
     h = WBCfg.size();
 
   if( ( h < WBCfg.size() ) &&
       ( WBCfg[ h ] >= 0 ) && ( WBCfg[ h ] < int( v_Cfg.size() ) ) )
    if( auto c = dynamic_cast< BlockConfig * >( v_Cfg[ WBCfg[ h ] ] ) )
     BCi = c;
   }

  if( BCi )
   BCi->apply( csbi );

  // now construct the LagBFunction; note that doing so may cause the
  // Objective of the inner Block (and therefore the Variable) to be
  // defined because LagBFunction needs it, hence this is only done
  // after the BlockConfig-uration, in particular for the case where the
  // sub-Block is R3B-copied
  auto lbfi = new LagBFunction( csbi );

  // since the Block is locked by f_if, lend the same identity to the
  // LagBFunction in case they need to lock it
  lbfi->set_id( f_id );

  // surely now the Objective is defined: check that all the senses agree
  if( ! i )
   f_max = ( csbi->get_objective()->get_sense() == Objective::eMax );
  else
   if( f_max != ( csbi->get_objective()->get_sense() == Objective::eMax ) )
    throw( std::invalid_argument(
	      "LagrangianDualSolver: mixed min/max sub-Block Objective" ) );

  v_LBF[ i ] = lbfi;
  auto osbi = new FRealObjective( sbi , lbfi );
  osbi->set_sense( f_max ? Objective::eMin : Objective::eMax , eNoMod );
  sbi->set_objective( osbi );

  LagrDual->add_nested_Block( sbi );  // add the sub-Block
  }

 // sort the sub-Block dictionary by Block address
 std::sort( blck_to_idx.begin() , blck_to_idx.end() );

 // BlockConfig-ure the Lagrangian Dual Block as a whole - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if a BlockConfig is not present but a name is, load it from file
 if( ( ! f_BCfg ) && ( ! LDBlck_BCfg.empty() ) ) {
  auto BC = Configuration::deserialize( LDBlck_BCfg );
  f_BCfg = dynamic_cast< BlockConfig * >( BC );
  if( ! f_BCfg ) {
   delete BC;
   throw( std::invalid_argument(
		   "LagrangianDualSolver: invalid BlockConfig from file" ) );

   }
  }

 // if a BlockConfig is present, apply() it
 if( f_BCfg )
  f_BCfg->apply( LagrDual );

 // check conditions on f_Block- - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // first ensure that the abstract representation is available; now all
 // possible BlockConfig-urations have been done
 f_Block->generate_abstract_variables();
 f_Block->generate_abstract_constraints();
 f_Block->generate_objective();

 // the Block must not contain any variable- - - - - - - - - - - - - - - - -
    
 if( ! f_Block->get_static_variables().empty() )
  throw( std::invalid_argument(
		    "LagrangianDualSolver: static Variable not allowed" ) );
    
 if( ! f_Block->get_dynamic_variables().empty() )
  throw( std::invalid_argument(
		   "LagrangianDualSolver: dynamic Variable not allowed" ) );
    
 // there must be no Objective- - - - - - - - - - - - - - - - - - - - - - - -
 // or it must be "empty"

 if( auto obj = f_Block->get_objective() ) {
  if( obj->get_num_active_var() != 0 )
   throw( std::invalid_argument(
		  "LagrangianDualSolver: nonempty Objective not allowed" ) );

  if( f_max != ( obj->get_sense() == Objective::eMax ) )
   throw( std::invalid_argument(
	  "LagrangianDualSolver: Block sense differs form sub-Block one" ) );
  }
 
 // count and check the FRowConstraint in the Block - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = 0;

 // count and check the static FRowConstraint - - - - - - - - - - - - - - - -
 // meanwhile construct the static dictionaries

 // number of groups of static constraints
 auto scn = f_Block->get_static_constraints().size();

 // resize the static constraints<-->Lagrangian-variables dictionaries
 scon_to_idx.resize( scn );
 idx_to_scon.resize( scn );

 {
  Index pos = 0;
  for( const auto & el : f_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] =
			 con_int_int( & var , NumVar , 1 );
			idx_to_scon[ pos++ ] = int_const( NumVar++ , & var );
		        } ) )
    continue;
   // Vectors
   if( un_any_thing_1( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] =
			 con_int_int( var.data() , NumVar , var.size() );
			idx_to_scon[ pos++ ] =
			 int_const( NumVar , var.data() );
			NumVar += var.size();
		        } ) )
    continue;
   // Multiarrays
   if( un_any_thing_K( FRowConstraint , el ,
		       {
			scon_to_idx[ pos ] =
			 con_int_int( var.data() , NumVar ,
				      var.num_elements() );
			idx_to_scon[ pos++ ] =
			 int_const( NumVar , var.data() );
			NumVar += var.num_elements();
		        } ) )
    continue;
   throw( std::invalid_argument(
       "LagrangianDualSolver: static constraint not a FRowConstraint" ) );
   }
  }

 static_cons = NumVar;

 // sort the static constraints-->Lagrangian-variables dictionary
 // this is not necessary for the Lagrangian-variables-->static constraints
 // one since it's surely sorted already
 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 // count and check the dynamic FRowConstraint- - - - - - - - - - - - - - - -
 for( const auto & el : f_Block->get_dynamic_constraints() ) {
  // Singles lists
  if( un_any_thing_0( std::list< FRowConstraint > , el ,
		      { NumVar += var.size(); } ) )
   continue;
  // Vectors of lists
  if( un_any_thing_1( std::list< FRowConstraint > , el ,
                      {
		       for( auto & lel: var )
			NumVar += lel.size();
		       } ) )
    continue;
  // Multiarrays of lists
  if( un_any_thing_K( std::list< FRowConstraint > , el ,
		      {
		       auto it = var.data();
		       for( auto i = var.num_elements() ; i-- ; ++it )
			NumVar += it->size();
		       } ) )
   continue;
  throw( std::invalid_argument(
         "LagrangianDualSolver: dynamic constraint not a FRowConstraint" ) );
  }

 // create the static and dynamic Lagrangian variables- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto Ls = new std::vector< ColVariable >( static_cons );
 auto Ld = new std::list< ColVariable >( NumVar - static_cons );

 // pass the Lagrangian variables to the Lagrangian Dual
 LagrDual->add_static_variable( *Ls , "Lambda_s" );
 LagrDual->add_dynamic_variable( *Ld , "Lambda_d" );

 // scan all FRowConstraints- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // meanwhile construct the linear objective function

 v_coeff_pair objcf( NumVar );
 auto objit = objcf.begin();

 // construct the auxiliary data structure to hold the Lagrangian terms;
 // LagTerms[ i ][ h ] contains the v_coeff_pair corresponding to the
 // Lagrangian term of sub-Block h for the i-th variable

 std::vector< std::vector< v_coeff_pair > > LagTerms( NumVar );
 auto LTit = LagTerms.begin();

 // scan all static FRowConstraints - - - - - - - - - - - - - - - - - - - - -
 {
  auto Lit = Ls->begin();

  // define a lambda that does the job
  auto scan = [ & ]( FRowConstraint & con ) -> void {
   // initialize the value of the Lagrangian variable with the current dual
   // solution of the FRowConstraint, for the odd chance that someone has
   // already put there a meaningful value
   Lit->set_value( con.get_dual() );
   
   // check the LHS/RHS
   auto lhs = con.get_lhs();
   auto rhs = con.get_rhs();

   if( ( ( lhs == -INFshift ) && ( rhs == INFshift ) ) || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( & *(Lit++) , 0 );
    ++LTit;
    return;
    }

   auto coef = constr2val( con , *Lit );

   // write the coefficient in the objective
   *(objit++) = std::make_pair( & *(Lit++) , coef );

   // split the linear constraint among the sub-Block
   split_constraint( con , *(LTit++) );
   };

  // finally apply the lambda to all static constraints
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el , scan , un_any_type< FRowConstraint >() );
  }

 // scan all dynamic FRowConstraints- - - - - - - - - - - - - - - - - - - - -
 // meanwhile construct the dynamic dictionaries

 // resize the dynamic constraints<-->Lagrangian-variables dictionaries
 dcon_to_idx.resize( NumVar - static_cons );
 idx_to_dcon.resize( NumVar - static_cons );

 {
  Index i = static_cons;
  auto Lit = Ld->begin();
  auto dc2iit = dcon_to_idx.begin();
  auto i2dcit = idx_to_dcon.begin();

  // define a lambda that does the job
  auto scan = [ & ]( FRowConstraint & con ) -> void {
   // initialize the value of the Lagrangian variable with the current dual
   // solution of the FRowConstraint, for the odd chance that someone has
   // already put there a meaningful value
   Lit->set_value( con.get_dual() );

   // first write the dictionaries
   *(dc2iit++) = std::make_pair( & con , i++ );
   *(i2dcit++) = & con;

   // then check the LHS/RHS
   auto lhs = con.get_lhs();
   auto rhs = con.get_rhs();

   if( ( ( lhs == -INFshift ) && ( rhs == INFshift ) ) || con.is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( & *(Lit++) , 0 );
    ++LTit;
    return;
    }

   auto coef = constr2val( con , *Lit );

   // write the coefficient in the objective
   *(objit++) = std::make_pair( & *(Lit++) , coef );

   // split the linear constraint among the sub-Block
   split_constraint( con , *(LTit++) );
   };

  // finally apply the lambda to all dynamic constraints
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_dynamic( el , scan , un_any_type< FRowConstraint >() );
  }

 // sort the dynamic constraints-->Lagrangian-variables dictionary
 // this must not be done for the Lagrangian-variables-->dynamic constraints
 // one since the mapping is positional
 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

 // create the Objective of the Lagrangian Dual - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto obj = new FRealObjective( LagrDual ,
				new LinearFunction( std::move( objcf ) ) );
 obj->set_sense( f_max ? Objective::eMin : Objective::eMax , eNoMod );
 LagrDual->set_objective( obj , eNoMod );

 // pass the Lagrangian terms to the corresponding LagBFunction - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index h = 0 ; h < f_nsb ; ++h ) {
  v_dual_pair dp( NumVar );  // construct the dual pairs

  Index i = 0;
  for( ; i < static_cons ; ++i ) {
   dp[ i ].first = & (*Ls)[ i ];
   dp[ i ].second = new LinearFunction( std::move( LagTerms[ i ][ h ] ) );
   }

  auto Lit = Ld->begin();
  for( ; i < NumVar ; ++i ) {
   dp[ i ].first = & (*Lit++);
   dp[ i ].second = new LinearFunction( std::move( LagTerms[ i ][ h ] ) );
   }

  // Since this LagBFunction is already part of an Objective (and, therefore,
  // has an Observer), the method LagBFunction::set_dual_pairs() cannot be
  // invoked. Notice also that the FRealObjective to which this LagBFunction
  // belongs cannot be constructed here, nor that this LagBFunction is set as
  // the Function of its FRealObjective. This is because this FRealObjective
  // needs to have its LagBFunction when the method split_constraint() (and
  // thus the method Block2Index()) is invoked above.

  v_LBF[ h ]->add_dual_pairs( std::move( dp ) );
  }

 // release the Block- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this must be done before the BlockSolverConfig-uration, because
 // Solver::set_Block() is called there, and it can reasonably lock() the
 // sub-Block, which therefore need be un-lock()-ed

 // since the Block is going to be unlock()-ed, retire the identity from the
 // LagBFunctions
 for( auto lbf : v_LBF )
  lbf->set_id();

 if( ! owned )
  f_Block->unlock( f_id );

 // BlockSolverConfig-ure the individual inner Block- - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // create the default BlockSolverConfig
 if( ! LagBF_BSCfg.empty() ) {
  auto c = Configuration::deserialize( LagBF_BSCfg );
  f_DBSCfg = dynamic_cast< BlockSolverConfig * >( c );
  if( ! f_DBSCfg ) {
   delete c;
   throw( std::invalid_argument(
	   "LagrangianDualSolver: LagBF_BCfg not a BlockSolverConfig" ) );
   }
  }

 Index iW2BSCfg = 0;  // index in W2BSCfg
 for( Index i = 0 ; i < f_nsb ; ++i ) {
  auto BSCi = f_DBSCfg;

  if( ! WBSCfg.empty() ) {  // individual BlockSolverConfig provided
   Index h;                 // the index in WBSCfg

   if( W2BSCfg.empty() )    // in dense format
    h = i;
   else                     // in sparse format
    if( ( iW2BSCfg < W2BSCfg.size() ) && ( W2BSCfg[ iW2BSCfg ] == i ) )
     h = iW2BSCfg++;
    else
     h = WBSCfg.size();
 
   if( ( h < WBSCfg.size() ) &&
       ( WBSCfg[ h ] >= 0 ) && ( WBSCfg[ h ] < int( v_Cfg.size() ) ) )
    if( auto c = dynamic_cast< BlockSolverConfig * >( v_Cfg[ WBSCfg[ h ] ] ) )
     BSCi = c;
   }

  if( BSCi ) {
   Block * csbi = v_LBF[ i ]->get_inner_block();
   if( CloneCfg ) {
    auto cBSCi = BSCi->clone();
    cBSCi->apply( csbi );
    delete cBSCi;
    }
   else
    BSCi->apply( csbi );
   }
  }

 // BlockSolverConfig-ure the Lagrangian Dual Block as a whole - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if a BlockSolverConfig is not present but a name is, load it from file
 if( ( ! f_BSCfg ) && ( ! LDBlck_BSCfg.empty() ) ) {
  auto BSC = Configuration::deserialize( LDBlck_BSCfg );
  f_BSCfg = dynamic_cast< BlockSolverConfig * >( BSC );
  if( ! f_BSCfg ) {
   delete BSC;
   throw( std::invalid_argument(
	    "LagrangianDualSolver: invalid BlockSolverConfig from file" ) );
   }
  }

 // if a BlockSolverConfig is present, apply() it
 if( f_BSCfg ) {
  if( CloneCfg ) {
   auto cBSC = f_BSCfg->clone();
   cBSC->apply( LagrDual );
   delete cBSC;
   }
  else
   f_BSCfg->apply( LagrDual );
  }
 
 #if CHECK_DS & 1
  // check that the Lagrangian Dual Block is correct - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LagrDual->is_correct();
 #endif

 // register the inner Solver to the Lagrangian Dual Block - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that the inner Solver could in principle cause some other
 // Block[Solver]Config-uration to the inner Block of the LagBFunction via
 // their ComputeConfig; since there are plenty of other ways to obtain the
 // same result, this should not happen at least for BlockConfig if the
 // corresponding information was needed before

 register_inner_Solver();

 // register the appropriate UpdateSolver- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is done *after* all other Block[Solver]Config-uration to ensure that
 // the UpdateSolver do not mess up with the first Solver registered to the
 // inner Block of the LagBFunction, that is the "crucial" one
 //
 // note that the issue only arises when the original sub-Block is used, as
 // otherwise the UpdateSolver is attached to it but the LagBFunction uses a
 // copy and therefore the UpdateSolver is not registered there
 //
 // HOWEVER, THIS IS CONCEPTUALLY WRONG, AS A Modification COULD ARISE FROM
 // THE UN-lock()-ED Block AND BE MISSED BECAUSE THE UpdateSolver IS NOT
 // PROPERLY REGISTERED YET!!
 //
 // solving this likely requires re-thinking to the whole idea that Solver
 // are "positional", and add the concept of "Solver name" that can be used
 // in BlockSolverConfig to set the Solver and in LagBFunction to have a
 // positional-independent notion of what the inner Solver is. no time to
 // do this properly now
 
 v_US.resize( f_nsb );  // meanwhile handle the UpdateSolver

 if( iBCopy )  // copying the sub-Block 
  // register an UpdateSolver to the original sub-Block so that any
  // Modification to it is immediately forwarded to the copy
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   v_US[ i ] = new UpdateSolver( v_LBF[ i ]->get_inner_block() );
   f_Block->get_nested_Block( i )->register_Solver( v_US[ i ] );
   }
 else          // evicting the sub-Block
  // register an UpdateSolver to the original sub-Block (which now lives
  // in the LagBFunction) so that any Modification to it is immediately
  // forwarded to the former father as it were still its son; note that
  // the UpdateSolver forwards, as opposed to map_*, the Modification
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   v_US[ i ] = new UpdateSolver( f_Block , nullptr , 2 );
   v_LBF[ i ]->get_inner_block()->register_Solver( v_US[ i ] );
   }

 // release the Block- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // and now, finally, all is done

 }  // end( LagrangianDualSolver::set_Block )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( int_LDSlv_iBCopy ):
   if( LagrDual )
    throw( std::logic_error( "changing iBCopy with registered Block" ) );
   iBCopy = bool( value );
   break;
  case( int_LDSlv_NNMult ):
   if( LagrDual )
    throw( std::logic_error( "changing NNMult with registered Block" ) );
   NNMult = bool( value );
   break;
  case( int_LDSlv_CloneCfg ):
   CloneCfg = bool( value );
   break;
  default:
   InnerSolver->set_par(  int_par_lds( par ) , value );
  }
 }

/*--------------------------------------------------------------------------*/
/*!!
void LagrangianDualSolver::set_par( idx_type par , double value )
{
 InnerSolver->set_par( dbl_par_lds( par ) , value );
 }
 !!*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par , std::string && value )
{
 switch( par ) {
  case( str_LDSlv_ISName ): {
   std::string tval( SMSpp_classname_normalise( std::move( value ) ) );
   if( ISName != tval ) {
    ISName = tval;
    unregister_inner_Solver();
    auto ts = new_Solver( ISName );
    InnerSolver = dynamic_cast< CDASolver * >( ts );
    if( ! InnerSolver ) {
     delete ts;
     throw( std::logic_error( ISName + " not a CDASolver" ) );
     }
    register_inner_Solver();
    }
   break;
   }
  case( str_LagBF_BCfg ):
   LagBF_BCfg = value;
   break;
  case( str_LagBF_BSCfg ):
   LagBF_BSCfg = value;
   break;
  case( str_LDBlck_BCfg ):
   LDBlck_BCfg = value;
   break;
  case( str_LDBlck_BSCfg ):
   LDBlck_BSCfg = value;
   break;
  default:
   InnerSolver->set_par( str_par_lds( par ) , std::move( value  ) );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par ,
				    std::vector< int > && value )
{
 switch( par ) {
  case( vint_LDSl_WBCfg ):
   WBCfg = value;
   break;
  case( vint_LDSl_W2BCfg ):
   W2BCfg = value;
   break;
  case( vint_LDSl_WBSCfg ):
   WBSCfg = value;
   break;
  case( vint_LDSl_W2BSCfg ):
   W2BSCfg = value;
   break;
  default:
   InnerSolver->set_par( vint_par_lds( par ) , std::move( value ) );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_par( idx_type par ,
				    std::vector< std::string > && value )
{
 switch( par ) {
  case( vstr_LDSl_Cfg ):
   FCfg = value;
   for( auto el : v_Cfg )
    delete el;
   v_Cfg.resize( FCfg.size() );
   for( Index i = 0 ; i < v_Cfg.size() ; ++i )
    v_Cfg[ i ] = Configuration::deserialize( FCfg[ i ] );
   break;
  default:
   InnerSolver->set_par( vstr_par_lds( par ) , std::move( value ) );
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::set_ComputeConfig( ComputeConfig * scfg )
{
 if( ! scfg ) {  // factory reset
  delete f_BCfg;
  f_BCfg = nullptr;
  delete f_BSCfg;
  f_BSCfg = nullptr;
  ThinComputeInterface::set_ComputeConfig();
  return;
  }

 // first of all check if the inner Solver is changing and act upon it
 for( const auto & pair : scfg->str_pars )
  if( pair.first == "str_LDSlv_ISName" )
   set_par( str_LDSlv_ISName , std::string( pair.second ) );

 // now call the base ThinComputeInterface to do the bulk of work; note
 // that str_LDSlv_ISName is called twice (if ever), but set_par( string )
 // checks if the class remains the same and does nothing
 ThinComputeInterface::set_ComputeConfig( scfg );

 // finally, take care of extra_Configuration (if any)
 if( ! scfg->f_extra_Configuration )
  return;
 
 if( auto scpp = dynamic_cast< p_SConf_p_p >(
			                   scfg->f_extra_Configuration ) ) {
  if( scpp->f_value.first ) {
   clear_LD_BlockSolverConfig();
   f_BSCfg = dynamic_cast< p_BSC >( scpp->f_value.first );
   if( ! f_BSCfg )
    throw( std::invalid_argument(
       "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config.fist"
				 ) );
   scpp->f_value.first = nullptr; // set it to nullptr since it has been
                                  // "extracted"
   }

  if( scpp->f_value.second ) {
   clear_LD_BlockConfig();
   f_BCfg = dynamic_cast< p_BC >( scpp->f_value.second );
   if( ! f_BCfg )
    throw( std::invalid_argument(
       "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config.second"
				 ) );
   scpp->f_value.second = nullptr; // set it to nullptr since it has been
                                   // "extracted"
   }
  return;
  }

 if( auto BSC = dynamic_cast< p_BSC >( scfg->f_extra_Configuration ) ) {
  clear_LD_BlockSolverConfig();
  f_BSCfg = BSC;
  scfg->f_extra_Configuration = nullptr; // set it to nullptr since it has
                                         // been "extracted"
  return;
  }
  
 if( auto BC = dynamic_cast< p_BC >( scfg->f_extra_Configuration ) ) {
  clear_LD_BlockConfig();
  f_BCfg = BC;
  scfg->f_extra_Configuration = nullptr; // set it to nullptr since it has
                                         // been "extracted"
  return;
  }

 throw( std::invalid_argument(
         "LagrangianDualSolver::set_ComputeConfig: invalid extra_Config" ) );
   
 }  // end( LagrangianDualSolver::set_ComputeConfig )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::compute( bool changedvars )
{
 lock();  // lock the mutex

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

 process_outstanding_Modification();

 if( ! owned )
  f_Block->unlock( f_id );

 // if iBCopy == false, inhibit all Modification from the UpdateSolver; these
 // would reach the (disconnected) original Block, but there is no reason for
 // this because these are all "local" changes that will be undone at the end
 // of the solution process, so the Solver attached to the father (and other
 // ancestors) have no reason to act upon them
 if( ! iBCopy )
  for( auto us : v_US )
   us->inhibit_Modification( true );

 auto res = InnerSolver->compute( changedvars );

 // if iBCopy == false, bring back the inner Block to its original objective
 // "like if nothing ever happened", and re-enable the Modification from the
 // UpdateSolver
 if( ! iBCopy ) {
  for( auto lbf : v_LBF )
   lbf->cleanup_inner_objective();
  for( auto us : v_US )
   us->inhibit_Modification( false );
  }

 // because the inner Solver is solving the dual of the original Block,
 // the unbounded an unfeasible return states have to be exchanged
 if( res == kUnbounded )
  res = kInfeasible;
 else
  if( res == kInfeasible )
   res = kUnbounded;

 unlock();  // unlock the mutex

 return( res );

 }  // end( LagrangianDualSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_var_solution( Configuration * solc )
{
 if( ! LagrDual )
  throw( std::logic_error(
    "LagrangianDualSolver::get_var_solution: Lagrangian Dual not formed" ) );

 // first ensure that the optimal convex multipliers are written as the
 // important_linearization_coefficients() of the LagBFunction
 InnerSolver->get_dual_solution();

 // define a lambda that does the solution (computation and) retrieval
 // for a specific sub-Block
 auto getsoli = [ this ] ( Index i ) -> void {
  Index szi = v_LBF[ i ]->get_int_par( C05Function::intGPMaxSz );
  if( ! szi )
   throw( std::invalid_argument(
           "LagrangianDualSolver::get_var_solution: no Solution stored" ) );
  auto & lc = v_LBF[ i ]->get_important_linearization_coefficients();
  if( lc.empty() )
   throw( std::invalid_argument(
       "LagrangianDualSolver::get_var_solution: no coefficients stored" ) );

  Index pos;
  if( lc.size() == 1 )        // the solution is already computed
   pos = lc.front().first;    // this is its position
  else {
   // the solution need be computed: find a free spot where to put it
   for( pos = 0 ; pos < szi ; ++pos )
    if( ! v_LBF[ i ]->is_linearization_there( pos ) )
     break;

   // if no free spot can be found, put it anywhere
   if( pos == szi )
    pos = 0;

   // compute the solution and store it in the global pool
   v_LBF[ i ]->store_combination_of_linearizations( lc , pos );

   // now the important linearization is in the pool, recall this
   LinearCombination nlc( { std::pair( pos ,
				       Function::FunctionValue( 1 ) ) } );
   v_LBF[ i ]->set_important_linearization( std::move( nlc ) );
   }

  // if necessary move the solution back from the global pool to the Block
  v_LBF[ i ]->global_pool_to_block( pos );

  // if sub-Block is a copy, map_back the solution to the original
  if( iBCopy )
    f_Block->get_nested_Block( i )->map_back_solution(
			    v_LBF[ i ]->get_nested_Block( 0 ) , nullptr );
  };

 auto SCvi = dynamic_cast< SimpleConfiguration< std::vector< int > > * >( solc );
 if( SCvi ) {
  assert( std::is_sorted( SCvi->value().begin() , SCvi->value().end() ) );
  if( ( ! SCvi->value().empty() ) && ( ( SCvi->value().front() < 0 ) ||
      ( Index( SCvi->value().back() ) >= f_nsb ) ) )
   throw( std::invalid_argument(
  "LagrangianDualSolver::get_var_solution: wrong indices in solc->value" ) );

  for( auto el : SCvi->value() )
   getsoli( el );
  }
 else
  for( Index i = 0 ; i < f_nsb ; ++i )
   getsoli( i );
 
 }  // end( LagrangianDualSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::get_dual_solution( Configuration * solc )
{
 if( ! LagrDual )
  throw( std::logic_error(
    "LagrangianDualSolver::get_var_solution: Lagrangian Dual not formed" ) );

 // define a Lambda which does the configuration - - - - - - - - - - - - - -
 // Notice that this lambda cannot be static in its current form, since it
 // captures the "this" pointer. If "lcfg" were static, the "this" pointer
 // within this lambda would always point to the same LagrangianDualSolver
 // object (namely, the one associated with the initialization of "lcfg"),
 // making this code wrong if other instances of LagrangianDualSolver call
 // this function.
 auto lcfg = [ this ]( Index b , Configuration * cfg ) {
  auto LSBb = v_LBF[ b ]->get_nested_Block( 0 );
  if( LSBb->get_registered_solvers().empty() )
   return;

  if( auto SBSb = dynamic_cast< CDASolver * >(
				LSBb->get_registered_solvers().front() ) ) {
   SBSb->get_dual_solution( cfg );
   if( iBCopy )  // the sub-Block is a copy
    f_Block->get_nested_Block( b )->map_back_solution( LSBb , nullptr , cfg );
   }
  };
 
 // if solc == nullptr get the dual solutions of all sub-Block with nullptr
 // Configuration- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! solc ) {
  for( Index i = 0 ; i < f_nsb ; ++i )
   lcfg( i , nullptr );

  goto get_duals;  // then go to also get those of the relaxed constraints
  }

 // the case of std::vector< std::pair< int , Configuration * > >- - - - - -
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< std::pair<
                            int  , Configuration * > > > * >( solc ) ) {
  bool get_rel = false;
  for( auto el : c->f_value )
   if( ( el.first < 0 ) || ( el.first >= f_nsb ) )
    get_rel = true;
   else
    lcfg( el.first , el.second );

  if( get_rel )
   goto get_duals;

  return;
  }

 // the case of std::vector< std::pair< int , int > >- - - - - - - - - - - -
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< std::pair<
                            int  , int > > > * >( solc ) ) {
  bool get_rel = false;
  for( auto el : c->f_value )
   if( ( el.first < 0 ) || ( el.first >= f_nsb ) )
    get_rel = true;
   else
    lcfg( el.first ,
	  ( ( el.second >= 0 ) && ( el.second < v_Cfg.size() ) ) ?
	  v_Cfg[ el.second ] : nullptr );

  if( get_rel )
   goto get_duals;

  return;
  }

 // the case of std::vector< Configuration * > - - - - - - - - - - - - - - -
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< Configuration *
                                                 > > * >( solc ) ) {
  for( Index i = 0 ; i < c->f_value.size() ; ++i )
   if( i > f_nsb )
    goto get_duals;
   else
    lcfg( i , c->f_value[ i ] );

  return;
  }

 // the case of std::vector< int > - - - - - - - - - - - - - - - - - - - - -
 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< int > > * >(
								  solc ) ) {
  for( Index i = 0 ; i < c->f_value.size() ; ++i )
   if( i > f_nsb )
    goto get_duals;
   else {
    auto h = c->f_value[ i ];
    lcfg( i , ( ( h >= 0 ) && ( h < v_Cfg.size() ) ) ? v_Cfg[ h ] : nullptr );
    }

  return;  // if none of the above, do nothing
  }

 throw( std::invalid_argument(
    "LagrangianDualSolver::get_dual_solution: wrong Configuration type" ) );

 // get the duals of the relaxed constraints - - - - - - - - - - - - - - - -

 get_duals:

 auto Ls = LagrDual->get_static_variable_v< ColVariable >( "Lambda_s" );
 auto Lsit = Ls->begin();
 auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );
 auto Ldit = Ld->begin();

 if( NNMult ) {
  // get the static part
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el , [ & ]( FRowConstraint & con ) -> void {
     auto val = (Lsit++)->get_value();
     if( to_be_reversed( con ) )
      val = - val;
     con.set_dual( val );
     } , un_any_type< FRowConstraint >() );
  // get the dynamic part
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_static( el , [ & ]( FRowConstraint & con ) -> void {
     auto val = (Ldit++)->get_value();
     if( to_be_reversed( con ) )
      val = - val;
     con.set_dual( val );
     } , un_any_type< FRowConstraint >() );
  }
 else {
  // get the static part
  for( const auto & el : f_Block->get_static_constraints() )
   un_any_const_static( el , [ & ]( FRowConstraint & con ) -> void {
                               con.set_dual( (Lsit++)->get_value() );
                               } , un_any_type< FRowConstraint >() );
  // get the dynamic part
  for( const auto & el : f_Block->get_dynamic_constraints() )
   un_any_const_static( el , [ & ]( FRowConstraint & con ) -> void {
                               con.set_dual( (Ldit++)->get_value() );
                               } , un_any_type< FRowConstraint >() );
  }
 }  // end( LagrangianDualSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/

int LagrangianDualSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( int_LDSlv_iBCopy ):   return( iBCopy );
  case( int_LDSlv_NNMult ):   return( NNMult );
  case( int_LDSlv_CloneCfg ): return( CloneCfg );
  }

 return( InnerSolver->get_int_par( int_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualSolver::get_dbl_par( idx_type par ) const
{
 return( InnerSolver->get_dbl_par( dbl_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & LagrangianDualSolver::get_str_par( idx_type par ) const
{
 switch( par ) {
  case( str_LDSlv_ISName ): return( ISName );
  case( str_LagBF_BCfg ):   return( LagBF_BCfg );
  case( str_LagBF_BSCfg ):  return( LagBF_BSCfg );
  case( str_LDBlck_BCfg ):  return( LDBlck_BCfg );
  case( str_LDBlck_BSCfg ): return( LDBlck_BSCfg );
  }

 return( InnerSolver->get_str_par( str_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & LagrangianDualSolver::get_vint_par( idx_type par )
 const {
 switch( par ) {
  case( vint_LDSl_WBCfg ):   return( WBCfg );
  case( vint_LDSl_W2BCfg ):  return( W2BCfg );
  case( vint_LDSl_WBSCfg ):  return( WBSCfg );
  case( vint_LDSl_W2BSCfg ): return( W2BSCfg );
  }

 return( InnerSolver->get_vint_par( vint_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & LagrangianDualSolver::get_vstr_par(
					              idx_type par ) const {
 if( par == vstr_LDSl_Cfg )
  return( FCfg );

 return( InnerSolver->get_vstr_par( vstr_par_lds( par ) ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::clear_LD_BlockSolverConfig( bool keepcfg )
{
 if( ! f_BSCfg )
  return;

 if( LagrDual ) {
  if( keepcfg ) {
   auto BSC = f_BSCfg->clone();
   BSC->clear();
   BSC->apply( LagrDual );
   delete BSC;
   }
  else {
   f_BSCfg->clear();
   f_BSCfg->apply( LagrDual );
   }
  }

 if( ! keepcfg ) {
  delete f_BSCfg;
  f_BSCfg = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::clear_LD_BlockConfig( bool keepcfg )
{
 if( ! f_BCfg )
  return;

 if( LagrDual ) {
  if( keepcfg ) {
   auto BC = f_BCfg->clone();
   BC->clear();
   BC->apply( LagrDual );
   delete BC;
   }
  else {
   f_BCfg->clear();
   f_BCfg->apply( LagrDual );
   }
  }

 if( ! keepcfg ) {
  delete f_BCfg;
  f_BCfg = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::clear_inner_BlockSolverConfig( bool keepcfg )
{
 if( ( ! f_DBSCfg ) && ( v_Cfg.empty() || WBSCfg.empty() ) )
  return;

 if( LagrDual ) {
  auto DBSC = f_DBSCfg;
  if( DBSC ) {
   if( keepcfg )
    DBSC = f_DBSCfg->clone();
   DBSC->clear();
   }

  Index iW2BSCfg = 0;  // index in W2BSCfg
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   BlockSolverConfig * BSCi = nullptr;

   if( ! WBSCfg.empty() ) {  // individual BlockSolverConfig provided
    Index h;                 // the index in WBSCfg

    if( W2BSCfg.empty() )    // in dense format
     h = i;
    else                     // in sparse format
     if( ( iW2BSCfg < W2BSCfg.size() ) && ( W2BSCfg[ iW2BSCfg ] == i ) )
      h = iW2BSCfg++;
     else
      h = WBSCfg.size();
 
    if( ( h < WBSCfg.size() ) &&
	( WBSCfg[ h ] >= 0 ) && ( WBSCfg[ h ] < int( v_Cfg.size() ) ) )
     if( auto c = dynamic_cast< BlockSolverConfig * >( v_Cfg[ WBSCfg[ h ] ] ) )
      BSCi = c;
    }

   if( BSCi ) {
    if( keepcfg )
     BSCi = BSCi->clone();
    BSCi->clear();
    }
   else
    BSCi = DBSC;

   if( BSCi )
    BSCi->apply( v_LBF[ i ]->get_inner_block() );

   if( ! WBSCfg.empty() ) {
    if( W2BSCfg.empty() ) {
     if( i >= WBSCfg.size() )
      break;
     }
    else
     if( iW2BSCfg >= W2BSCfg.size() )
      break;
    }
   }
  }

 if( ! keepcfg ) {
  delete f_DBSCfg;
  f_DBSCfg = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::clear_inner_BlockConfig( bool keepcfg )
{
 if( ( ! f_DBCfg ) && ( v_Cfg.empty() || WBCfg.empty() ) )
  return;

 if( LagrDual ) {
  auto DBC = f_DBCfg;
  if( DBC ) {
   if( keepcfg )
    DBC = f_DBCfg->clone();
   DBC->clear();
   }

  Index iW2BCfg = 0;  // index in W2BCfg
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   BlockConfig * BCi = nullptr;

   if( ! WBCfg.empty() ) {  // individual BlockConfig provided
    Index h;                // the index in WBCfg

    if( W2BCfg.empty() )    // in dense format
     h = i;
    else                     // in sparse format
     if( ( iW2BCfg < W2BCfg.size() ) && ( W2BCfg[ iW2BCfg ] == i ) )
      h = iW2BCfg++;
     else
      h = WBCfg.size();
 
    if( ( h < WBCfg.size() ) &&
	( WBCfg[ h ] >= 0 ) && ( WBCfg[ h ] < int( v_Cfg.size() ) ) )
     if( auto c = dynamic_cast< BlockConfig * >( v_Cfg[ WBCfg[ h ] ] ) )
      BCi = c;
    }

   if( BCi ) {
    if( keepcfg )
     BCi = BCi->clone();
    BCi->clear();
    }
   else
    BCi = DBC;

   if( BCi )
    BCi->apply( v_LBF[ i ]->get_inner_block() );

   if( ! WBCfg.empty() ) {
    if( W2BCfg.empty() ) {
     if( i >= WBCfg.size() )
      break;
     }
    else
     if( iW2BCfg >= W2BCfg.size() )
      break;
    }
   }
  }

 if( ! keepcfg ) {
  delete f_DBCfg;
  f_DBCfg = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

Index LagrangianDualSolver::index_of_static_constraint(
						 const FRowConstraint * con )
{
 if( scon_to_idx.empty() )
  return( InINF );

 assert( std::is_sorted( scon_to_idx.begin() , scon_to_idx.end() ) );
 auto it = upper_bound( scon_to_idx.begin() , scon_to_idx.end() ,
                        std::make_tuple( con , 0 , 0 ),
                        []( const auto & a , const auto & b ) {
                         return( std::get< 0 >( a ) < std::get< 0 >( b ) );
			 } );

 // it now refers to the first (group of) element(s) greater than i
 if( it == scon_to_idx.begin() )  // all elements are greater
  return( InINF );                // it is not there

 --it;  // the previous group is the one it belongs to

 // compute the distance from the first element of the constraint group
 auto dist = std::distance( static_cast< const FRowConstraint * >(
					       std::get< 0 >( *it ) ) , con );

 if( ( dist >= 0 ) && ( Index( dist ) < std::get< 2 >( *it ) ) )
  return( std::get< 1 >( *it ) + Index( dist ) );  // it belongs to this group

 return( InINF );  // it doesn't exist

 }  // end( LagrangianDualSolver::index_of_static_constraint )

/*--------------------------------------------------------------------------*/

Index LagrangianDualSolver::index_of_dynamic_constraint(
						  const FRowConstraint * con )
{
 assert( std::is_sorted( dcon_to_idx.begin() , dcon_to_idx.end() ) );
 auto it = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
                        std::make_pair( con , 0 ) ,
                        []( const auto & a , const auto & b ) {
                         return( a.first < b.first );
                         } );

 if( ( it != dcon_to_idx.end() ) && ( it->first == con ) )
  return( it->second );

 return( InINF );
 }

/*--------------------------------------------------------------------------*/

FRowConstraint * LagrangianDualSolver::static_constraint_with_index( Index i )
{
 #ifdef NDEBUG
  if( ( i >= static_cons ) || idx_to_scon.empty() )
   throw( std::invalid_argument(
	       "LagrangianDualSolver::invalid index of static constraint" ) );

  assert( std::is_sorted( idx_to_scon.begin() , idx_to_scon.end() ) );
 #endif

 auto it = upper_bound( idx_to_scon.begin() , idx_to_scon.end() ,
                        std::make_pair( i , nullptr ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                         } );

 // it now refers to the first (group of) element(s) greater than i
 #ifdef NDEBUG
  if( it == idx_to_scon.begin() )  // all elements are greater
   throw( std::invalid_argument(
			 "LagrangianDualSolver::inconsistent idx_to_scon" ) );
 #endif

 --it;  // the previous group is the one it belongs to

 return( it->second + ( i - it->first ) );
 }

/*--------------------------------------------------------------------------*/

bool LagrangianDualSolver::to_be_reversed( const FRowConstraint & con )
{
 if( f_max ) {  // maximization problem
  if( con.get_lhs() == -Inf< RowConstraint::RHSValue >() )  // <= constraint
   return( true );
  }
 else           // minimization problem
  if( con.get_rhs() == Inf< RowConstraint::RHSValue >() )    // >= constraint
   return( true );

 return( false );
 }

/*--------------------------------------------------------------------------*/

double LagrangianDualSolver::constr2val( const FRowConstraint & con ,
					 ColVariable & lvar )
{
 // returns the coefficient to be put in the Linear Objective for the
 // Lagrangian variable lvar corresponding to the relaxed constraint con
 // note that con has one of the three forms
 //
 //           variable part == RHS
 //           variable part <= RHS
 //    LHS <= variable part
 //
 // which naturally correspond to Lagrangian terms of the form
 //
 //    lvar * ( variable part - RHS )
 //    lvar * ( variable part - LHS )
 //
 // The natural sign constraint on lvar is:
 //
 //  - if f_max, i.e., the Lagrangian Dual is a minimum because the
 //    original problem is a maximum
 //    = variable part <= RHS has a <= 0 Lagrangian multiplier
 //    = variable part >= LHS has a >= 0 Lagrangian multiplier
 //
 //  - if ! f_max, i.e., the Lagrangian Dual is a maximum because the
 //    original problem is a minimum
 //    = variable part <= RHS has a >= 0 Lagrangian multiplier
 //    = variable part >= LHS has a <= 0 Lagrangian multiplier
 //
 // (equality constraints never have sign-constrained Lagrangian multipliers)
 // if NNMult == true, <= 0 Lagrangian multipliers are made >= 0 by
 // multiplying everything by -1
 //
 // however, note the following: if the Lagrangian term is, say
 //
 //    lvar * ( variable part - RHS )
 //
 // then THE COEFFICIENT TO BE WRITTEN IN THE Linear Objective IS - RHS,
 // the "-" being crucial of course
 
 auto lhs = con.get_lhs();
 auto rhs = con.get_rhs();

 if( NNMult ) {
  if( lhs < rhs ) {                    // an inequality constraint
   lvar.is_positive( true , eNoMod );  // always a >= multiplier
   if( rhs == INFshift )               // a >= constraint
    return( f_max ? - lhs : lhs );

   if( lhs == -INFshift )              // a <= constraint
    return( f_max ? rhs : - rhs );

   throw( std::invalid_argument(
            "LagrangianDualSolver: ranged constraints not supported yet" ) );
   }

  return( - rhs );
  }

 // define the sign constraints on the multiplier (if any)
 if( f_max ) {  // for a max original problem
  if( rhs == INFshift ) {                // a >= constraint 
   lvar.is_positive( true , eNoMod );    // ==> a >= multiplier     
   return( - lhs );
   }

  if( lhs == -INFshift )                // a <= constraint 
   lvar.is_negative( true , eNoMod );   // ==> a <= multiplier
  }
 else {         // for a min original problem
  if( rhs == INFshift ) {               // a >= constraint 
   lvar.is_negative( true , eNoMod );   // ==> a <= multiplier
   return( - lhs );
   }

  if( lhs == -INFshift )                // a <= constraint 
   lvar.is_positive( true , eNoMod );   // ==> a >= multiplier
  }

 return( - rhs );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::split_constraint( const FRowConstraint & con ,
			 std::vector< LinearFunction::v_coeff_pair > & split )
{
 auto lf = dynamic_cast< const p_LF >( con.get_function() );
 if( ! lf )
  throw( std::invalid_argument(
			"LagrangianDualSolver: FRowConstraint not linear" ) );

 if( lf->get_constant_term() != 0 )
  throw( std::invalid_argument(
	   "LagrangianDualSolver: nonzero constant term not handled yet" ) );
 
 auto & vc = lf->get_v_var();

 split.resize( f_nsb );

 if( f_nsb == 1 ) {    // easy case: only one sub-Block, nothing to split
  split.front() = vc;
  if( NNMult && to_be_reversed( con ) )
   for( auto & el : split.front() )
    el.second = - el.second;

  return;
  }

 for( auto & el : split )
  el.clear();

 if( vc.empty() )   // easy case: empty constraint, nothing to split
  return;

 std::vector< Index > blckidx( vc.size() );  // Block to which the var belongs
 std::vector< Index > cntr( f_nsb , 0 );

 // first pass: count the size of each split[ h ]; meanwhile save the
 // Variable-to-sub-Block-index information in blckidx to avoid computing
 // it twice;
 for( Index i = 0 ; i < vc.size() ; ) {
  auto bi = Block2Index( vc[ i ].first->get_Block() );
  blckidx[ i++ ] = bi;
  ++cntr[ bi ];
  }

 // properly size all split[ h ]; meanwhile, reset the counter
 for( Index h = 0 ; h < f_nsb ; ++h )
  if( cntr[ h ] ) {
   split[ h ].resize( cntr[ h ] );
   cntr[ h ] = 0;
   }
 
 // second pass: construct all split[ h ]
 for( Index i = 0 ; i < vc.size() ; ++i )
  split[ blckidx[ i ] ][ cntr[ blckidx[ i ] ]++ ] = vc[ i ];

 // if necessary change the sign
 if( NNMult && to_be_reversed( con ) )
  for( auto & el : split )
   for( auto & lel : el )
    lel.second = - lel.second;

 }  // end( LagrangianDualSolver::split_constraint )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::cleanup_LagrDual( void )
{
 if( ! LagrDual )  // nothing to be cleaned up
  return;          // all done

 // first detach the inner Solver
 unregister_inner_Solver();

 // unregister and delete the UpdateSolver
 for( Index i = 0 ; i < f_nsb ; ++i )
  v_LBF[ i ]->get_inner_block()->unregister_Solver( v_US[ i ] , true );

 v_US.clear();

 clear_LD_BlockSolverConfig( true );
 clear_inner_BlockSolverConfig( true );
 clear_LD_BlockConfig( true );
 clear_inner_BlockConfig( true );

 // if necessary put back the sub_Block
 if( ! iBCopy ) {
  bool owned = f_Block->is_owned_by( f_id );
  if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
   throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

  // remove the sub-Block from the LagBFunction, but do not delete them;
  // rather,  re-attach them Block to their original father
  const auto & sb = f_Block->get_nested_Blocks();
  for( Index i = 0 ; i < f_nsb ; ++i ) {
   v_LBF[ i ]->set_inner_block( nullptr , false );
   sb[ i ]->set_f_Block( f_Block );
   }

  if( ! owned )
   f_Block->unlock( f_id );
  }

 // LagrDual is an AbstractBlock and therefore its destructor deletes
 // everything inside it, comprised the LagBFunction that therefore must
 // not to be deleted here
 delete LagrDual;
 LagrDual = nullptr;

 v_LBF.clear();

 }  // end( LagrangianDualSolver::cleanup_LagrDual )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::guts_of_destructor( void )
{
 unregister_inner_Solver();
 delete InnerSolver;
 InnerSolver = nullptr;
 clear_LD_BlockSolverConfig();
 clear_inner_BlockSolverConfig();
 clear_LD_BlockConfig();
 clear_inner_BlockConfig();
 cleanup_LagrDual();
		 
 }  // end( LagrangianDualSolver:guts_of_destructor )

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::flatten_Modification_list( Lst_sp_Mod & vmt ,
						      sp_Mod mod )
{
 if( const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod ) )
  for( auto submod : tmod->sub_Modifications() )
   flatten_Modification_list( vmt , submod );
 else
  // keep only Modification coming directly from f_Block, i.e., discard all
  // those coming from the sub-Block
  if( mod->get_Block() == f_Block )
   vmt.push_back( mod );
 }

/*--------------------------------------------------------------------------*/

void LagrangianDualSolver::process_outstanding_Modification( void )
{
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 0-th loop: "atomically flatten" v_mod into a temporary list to better
 // handle it, then clear it; meanwhile discard all Modification not
 // coming directly from f_Block (i.e., coming from its sub-Block)

 Lst_sp_Mod v_mod_tmp;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod )
  flatten_Modification_list( v_mod_tmp , mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( v_mod_tmp.empty() )  // no Modification coming directly from f_Block
  return;                 // all done

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 1st loop: only consider addition and deletion of (dynamic) FRowConstraint
 // Only consider addition and removal of constraints; changes to existing
 // constraints will be considered in the second loop, so that changes to
 // constraints that are going to be deleted or changes to added constraints
 // can be ignored. indeed, constraints will only be added at the end, hence
 // directly in their current state
 //
 // IMPORTANT NOTE: since deleted constraints are kept in the Modification
 // up until it is processed and destroyed, a deleted Constraint cannot be
 // found in the added list, and an added Constraint can only be deleted at
 // most once

 std::set< p_FRC > Dltds;  // set of pointers to deleted constraints
 std::vector< p_FRC > Addd;
 // pointers to added constraints in the order they have been added
 std::set< p_FRC > Addds;  // set of pointers to added constraints
 std::set< p_FRC > AddDltd;
 // set of pointers to constraints that have been added and then deleted
 
 bool to_delete;  // should have been defined inside, but there is not
                  // visible by the lambda

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // adding new (dynamic) FRowConstraint == (dynamic) Lagrangian variables
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = std::dynamic_pointer_cast<
                                 BlockModAdd< FRowConstraint > >( *imod ) ) {

   Addd.insert( Addd.end() , tmod->added().begin() , tmod->added().end() );
   Addds.insert( tmod->added().begin() , tmod->added().end() );
   to_delete = true;
   continue;
   }

  // removing new (dynamic) FRowConstraint == (dynamic) Lagrangian variables
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = std::dynamic_pointer_cast<
                                 BlockModRmv< FRowConstraint > >( *imod ) ) {

   if( Addd.empty() )  // nothing added yet, they can only be original constr
    for( auto el : tmod->removed() )
     Dltds.insert( & el );
   else {  // have to check if some was a just added constraints
    for( auto el : tmod->removed() ) {
     auto sit = Addds.find( & el );
     if( sit == Addds.end() )  // one of the original constraints
      Dltds.insert( & el );    // just mark it as removed
     else {                    // a previously added constraint
      AddDltd.insert( & el );  // mark it so
      Addds.erase( sit );      // remove it from the set of added
      // note: the element is *not* removed from the *vector* of added
      // ones since this would be a costly operation, this is done only
      // once at the end
      }
     }
    }
   
   to_delete = true;
   continue;
   }

  // any other Modification: no nothing, it'll be dealt with later
  to_delete = false;

  }  // end( first loop )

 // if necessary, adjust the std::vector of added constraints - - - - - - - -
 if( ( ! Addd.empty() ) && ( ! AddDltd.empty() ) ) {
  auto Awit = Addd.begin();

  // look up first added-then-deleted constraint
  while( AddDltd.find( *Awit ) == AddDltd.end() )
   ++Awit;

  // now copy skipping all the added-then-deleted constraints
  auto Arit = ++Awit;
  for( Index cnt = 1 ; cnt < AddDltd.size() ; ++Arit )
   if( AddDltd.find( *Awit ) == AddDltd.end() )
    *(Awit++) = *Arit;
   else
    ++cnt;

  // finish copying the last part after the last added-then-deleted constraint
  while( Arit != Addd.end() )
   *(Awit++) = *(Arit++);

  // consistency check
  assert( decltype( Addd )::size_type( std::distance( Addd.begin() , Awit ) )
	  == Addd.size() - AddDltd.size() );

  // resize the set of added constraints
  Addd.resize( std::distance( Addd.begin() , Awit ) );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 2nd loop: only consider Modification coming from the LinearFunction
 // inside the FRowConstraint or from changes of the LHS/RHS of the
 // FRowConstraint, ignoring any coming from a FRowConstraint that
 // has been added, deleted or both
 //
 // No particular care is taken into trying to bunch the Modification: each
 // Modification is "immediately reflected" into the LinearFunction inside
 // the LagBFunction, even if multiple Modification may in principle change
 // the same LinearFunction. The idea is that the InnerSolver will do this
 // anyway. The only care that is taken is to bunch all the Modification of
 // this type to the same LagBFunction into a single GroupModification, so
 // that the check if it is necessary to check if the Solutions remained
 // feasible (answer: no, the inner Block is not touched) is only done once
 //
 // note that we don't care to remove acted-upon Modification from the list:
 // everything that remains there at the end of the looop is anyway ignored

 // the chanel opened in each sub-Block (if any)
 std::vector< Block::ChnlName > chnls( f_nsb , 0 );

 std::set< Index > lrhschgd;
 // set of indices of (variables corresponding to) changed LHS/RHS
 Subset nms;  // indices of (variables ...) in the order they have been found
 LinearFunction::Vec_FunctionValue lrhsval;  // values of new LHS/RHS
 
 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ; ++imod ) {

  // the LHS/RHS of a FRowConstraint is changed - - - - - - - - - - - - - - -
  if( auto tmod =
      std::dynamic_pointer_cast< const RowConstraintMod >( *imod ) ) {
   auto cnst = static_cast< const p_FRC >( tmod->constraint() );
   // if the FRowConstraint is deleted or added (or both), do nothing
   if( ( Dltds.find( cnst ) != Dltds.end() ) ||
       ( Addds.find( cnst ) != Addds.end() ) ||
       ( AddDltd.find( cnst ) != AddDltd.end() ) )
    continue;

   auto pos = index_of_constraint( cnst );
   if( lrhschgd.find( pos ) != lrhschgd.end() )  // changed more than once
    continue;                                    // done already

   lrhschgd.insert( pos );
   nms.push_back( pos );
   RowConstraint::RHSValue lrhs;
   switch( tmod->type() ) {
    case( FRowConstraintMod::eChgLHS ):
     lrhs = cnst->get_lhs();
     if( lrhs == -INFshift )
      throw( std::logic_error( "LagrangianDualSolver: changing -INF LHS" ) );
     break;
    case( FRowConstraintMod::eChgRHS ):
     lrhs = cnst->get_rhs();
     if( lrhs == INFshift )
      throw( std::logic_error( "LagrangianDualSolver: changing INF RHS" ) );
     break;
    case( FRowConstraintMod::eChgBTS ):
     #ifndef NDEBUG
     {
      ColVariable * Lpos;
      if( pos < static_cons ) {
       auto Ld = LagrDual->get_static_variable_v< ColVariable >( "Lambda_s" );
       Lpos = & (*Ld)[ pos ];
       }
      else {
       auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );
       auto lvit = std::next( Ld->begin() , pos - static_cons );
       Lpos = & *lvit;
       }
      
      if( Lpos->is_positive() || Lpos->is_negative() )
       throw( std::logic_error(
            "LagrangianDualSolver: changing inequality constraint to equality"
			       ) );
      }
     #endif
     lrhs = cnst->get_rhs();
     break;
    default:
     throw( std::logic_error(
        "LagrangianDualSolver: relaxing/enforcing constraints not handled yet"
			     ) );
    }

   if( NNMult && to_be_reversed( *cnst ) )
    lrhs = - lrhs;
   lrhsval.push_back( - lrhs );  // note the necessary "-"

   }  // end( RowConstraintMod )

  // some new coefficients added to the LinearFunction- - - - - - - - - - - -
  if( auto tmod =
      std::dynamic_pointer_cast< const C05FunctionModVarsAddd >( *imod ) ) {
   auto lf = static_cast< const p_LF >( tmod->function() );
   auto cnst = static_cast< const p_FRC >( lf->get_Observer() );
   // if the FRowConstraint is deleted or added (or both), do nothing
   if( ( Dltds.find( cnst ) != Dltds.end() ) ||
       ( Addds.find( cnst ) != Addds.end() ) ||
       ( AddDltd.find( cnst ) != AddDltd.end() ) )
    continue;

   // have to split the added ColVariable among the sub-Block
   std::vector< v_coeff_pair > split( f_nsb );
   std::vector< Index > blckidx( tmod->vars().size() );
   // Block to which the var belongs
   std::vector< Index > cntr( f_nsb , 0 );

   // first pass: count the size of each split[ h ]; meanwhile save the
   // Variable-to-sub-Block-index information in blckidx to avoid computing
   // it twice;
   for( Index i = 0 ; i < tmod->vars().size() ; ) {
    auto bi = Block2Index( tmod->vars()[ i ]->get_Block() );
    blckidx[ i++ ] = bi;
    ++cntr[ bi ];
    }

   // properly size all split[ h ]; meanwhile, reset the counter
   for( Index h = 0 ; h < f_nsb ; ++h )
    if( cntr[ h ] ) {
     split[ h ].resize( cntr[ h ] );
     cntr[ h ] = 0;
     }

   // second pass: construct all split[ h ]
   for( Index i = 0 ; i < tmod->vars().size() ; ++i )
    split[ blckidx[ i ] ][ cntr[ blckidx[ i ] ]++ ] =
     coeff_pair( static_cast< ColVariable * >( tmod->vars()[ i ] ) ,
		 lf->get_coefficient( lf->is_active( tmod->vars()[ i ] ) )
		 );

   // if necessary change the sign
   if( NNMult && to_be_reversed( *cnst ) )
    for( auto & el : split )
     for( auto & lel : el )
      lel.second = - lel.second;

   // now call add_variables() for all the appropriate LinearFunction
   auto pos = index_of_constraint( cnst );
   for( Index h = 0 ; h < f_nsb ; ++h ) {
    if( ! cntr[ h ] )
     continue;

    if( ! chnls[ h ] )
     chnls[ h ] = LagrDual->get_nested_Block( h )->open_channel();

    static_cast< p_LF >( v_LBF[ h ]->get_Lagrangian_term( pos )
			 )->add_variables( std::move( split[ h ] ) ,
					   Observer::make_par( eModBlck ,
							       chnls[ h ] ) );
    }
   }  // end( C05FunctionModVarsAddd )

  // some coefficients removed from the LinearFunction- - - - - - - - - - - -
  // we capture any FunctionModVars; since the C05FunctionModVarsAddd have
  // been captured yet, all that remains are the C05FunctionModVarsRngd and
  // C05FunctionModVarsSbst, that we don't need to distinguish
  if( auto tmod =
      std::dynamic_pointer_cast< const FunctionModVars >( *imod ) ) {
   auto lf = static_cast< const p_LF >( tmod->function() );
   auto cnst = static_cast< const p_FRC >( lf->get_Observer() );
   // if the FRowConstraint is deleted or added (or both), do nothing
   if( ( Dltds.find( cnst ) != Dltds.end() ) ||
       ( Addds.find( cnst ) != Addds.end() ) ||
       ( AddDltd.find( cnst ) != AddDltd.end() ) )
    continue;

   // have to split the removed ColVariable among the sub-Block
   std::vector< Subset > split( f_nsb );
   std::vector< Index > blckidx( tmod->vars().size() );
   // Block to which the var belongs
   std::vector< Index > cntr( f_nsb , 0 );

   // first pass: count the size of each split[ h ]; meanwhile save the
   // Variable-to-sub-Block-index information in blckidx to avoid computing
   // it twice;
   for( Index i = 0 ; i < tmod->vars().size() ; ) {
    auto bi = Block2Index( tmod->vars()[ i ]->get_Block() );
    blckidx[ i++ ] = bi;
    ++cntr[ bi ];
    }

   // properly size all split[ h ]; meanwhile, reset the counter
   for( Index h = 0 ; h < f_nsb ; ++h )
    if( cntr[ h ] ) {
     split[ h ].resize( cntr[ h ] );
     cntr[ h ] = 0;
     }

   // second pass: construct all split[ h ]
   auto pos = index_of_constraint( cnst );
   for( Index i = 0 ; i < tmod->vars().size() ; ++i ) {
    auto bidx = blckidx[ i ];
    split[ bidx ][ cntr[ bidx ]++ ] =
     static_cast< p_LF >( v_LBF[ bidx ]->get_Lagrangian_term( pos )
			  )->is_active( tmod->vars()[ i ] );
    }

   // now call remove_variables() for all the appropriate LinearFunction
   for( Index h = 0 ; h < f_nsb ; ++h ) {
    if( ! cntr[ h ] )
     continue;

    if( ! chnls[ h ] )
     chnls[ h ] = LagrDual->get_nested_Block( h )->open_channel();

    static_cast< p_LF >( v_LBF[ h ]->get_Lagrangian_term( pos )
			 )->remove_variables( std::move( split[ h ] ) ,
					      false ,
					      Observer::make_par( eModBlck ,
								  chnls[ h ] )
					      );
    }
   }  // end( FunctionModVars )

  // some coefficients changed in the LinearFunction- - - - - - - - - - - - -
  if( auto tmod =
      std::dynamic_pointer_cast< const C05FunctionModLin >( *imod ) ) {
   auto lf = static_cast< const p_LF >( tmod->function() );
   auto cnst = static_cast< const p_FRC >( lf->get_Observer() );
   // if the FRowConstraint is deleted or added (or both), do nothing
   if( ( Dltds.find( cnst ) != Dltds.end() ) ||
       ( Addds.find( cnst ) != Addds.end() ) ||
       ( AddDltd.find( cnst ) != AddDltd.end() ) )
    continue;

   // have to split the changed ColVariable among the sub-Block
   std::vector< Subset > split( f_nsb );
   std::vector< Function::Vec_FunctionValue > delta( f_nsb );
   std::vector< Index > blckidx( tmod->vars().size() );
   // Block to which the var belongs
   std::vector< Index > cntr( f_nsb , 0 );

   // first pass: count the size of each split[ h ]; meanwhile save the
   // Variable-to-sub-Block-index information in blckidx to avoid computing
   // it twice;
   for( Index i = 0 ; i < tmod->vars().size() ; ) {
    auto bi = Block2Index( tmod->vars()[ i ]->get_Block() );
    blckidx[ i++ ] = bi;
    ++cntr[ bi ];
    }

   // properly size all split[ h ] and delta[ h ] and reset the counter
   for( Index h = 0 ; h < f_nsb ; ++h )
    if( cntr[ h ] ) {
     split[ h ].resize( cntr[ h ] );
     delta[ h ].resize( cntr[ h ] );
     cntr[ h ] = 0;
     }

   // second pass: construct all split[ h ] and delta[ h ]
   auto pos = index_of_constraint( cnst );
   for( Index i = 0 ; i < tmod->vars().size() ; ++i ) {
    auto bidx = blckidx[ i ];
    auto tc = cntr[ bidx ]++;
    split[ bidx ][ tc ] =
     static_cast< p_LF >( v_LBF[ bidx ]->get_Lagrangian_term( pos )
			  )->is_active( tmod->vars()[ i ] );
    delta[ bidx ][ tc ] = tmod->delta()[ i ];
    }

   // if necessary change the sign
   if( NNMult && to_be_reversed( *cnst ) )
    for( auto & el : delta )
     for( auto & lel : el )
      lel = - lel;

   // now call modify_coefficients() for all the appropriate LinearFunction
   for( Index h = 0 ; h < f_nsb ; ++h ) {
    if( ! cntr[ h ] )
     continue;

    if( ! chnls[ h ] )
     chnls[ h ] = LagrDual->get_nested_Block( h )->open_channel();

    auto lfh = static_cast< p_LF >( v_LBF[ h ]->get_Lagrangian_term( pos ) );
    auto & vcp = lfh->get_v_var();
    for( Index i = 0 ; i < cntr[ h ] ; ++i )
     delta[ h ][ i ] += vcp[ split[ h ][ i ] ].second;

    if( cntr[ h ] == 1 )
     lfh->modify_coefficient( split[ h ].front() , delta[ h ].front() ,  
			      Observer::make_par( eModBlck , chnls[ h ] ) );
    else
     lfh->modify_coefficients( std::move( delta[ h ] ) ,
			       std::move( split[ h ] ) , false ,
			       Observer::make_par( eModBlck , chnls[ h ] ) );
    }
   }  // end( C05FunctionModLin )

  // if it's anything else, ignore it - - - - - - - - - - - - - - - - - - - -

  }  // end( 2nd loop )

 // close any channel that has actually been opened - - - - - - - - - - - - -
 for( Index h = 0 ; h < f_nsb ; ++h )
  if( chnls[ h ] )
   LagrDual->get_nested_Block( h )->close_channel( chnls[ h ] );

 // if any LHS/RHS has changed, do the change in the Objective of LagrDual
 // note that nms is not ordered
 if( ! lrhschgd.empty() )
  static_cast< p_LF >( static_cast< p_FRO >( LagrDual->get_objective()
					     )->get_function()
		       )->modify_coefficients( std::move( lrhsval ) ,
					       std::move( nms ) , false );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now actually remove all Lagrangian variables from all the LagBFunction
 // and the LinearFunction: ensure that all the corresponding Modification
 // are bunched into a unique GroupModification of the LagrDual

 if( ! Dltds.empty() ) {
  // open a channel where to bunch all the removal Modifications
  const auto chnl = LagrDual->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );

  // the list from which the variable have to be removed
  auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );

  if( Dltds.size() == 1 ) {  // just one constraint
   auto cnst = *(Dltds.begin());
   Index i = index_of_constraint( cnst );

   // remove the variable in the LagBFunction
   for( Index h = 0 ; h < f_nsb ; )
    v_LBF[ h++ ]->remove_variable( i , mp );

   // remove the variable in the Objective
   static_cast< p_LF >( static_cast< p_FRO >( LagrDual->get_objective()
					      )->get_function()
			)->remove_variable( i , mp );

   // now actually remove the dynamic variable
   LagrDual->remove_dynamic_variable( *Ld , std::next( Ld->begin() ,
						       i - static_cons ) );    
   // now adjust the dictionaries
   --NumVar;

   // index to dynamic constraint
   auto i2dit = idx_to_dcon.begin() + i - static_cons;
   std::copy( i2dit + 1 , idx_to_dcon.end() , i2dit );
   idx_to_dcon.resize( idx_to_dcon.size() - 1 );

   // dynamic constraint to index
   auto d2iit = std::lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
				  std::make_pair( cnst , 0 ) ,
				  []( auto & p1 , auto & p2 ) {
				   return( p1.first < p2.first );
				   } );
   std::copy( d2iit + 1 , dcon_to_idx.end() , d2iit );
   dcon_to_idx.resize( dcon_to_idx.size() - 1 );
   for( auto & el : dcon_to_idx )
    if( el.second > i )
     --el.second;
   }
  else {  // multiple constraints
   Subset Dltdn( Dltds.size() );  // set of indices of deleted constraint
   auto Dnit = Dltdn.begin();
   for( auto el : Dltds )
    *(Dnit++) = index_of_constraint( el );

   std::sort( Dltdn.begin() , Dltdn.end() );

   // check if it actually was a range
   bool isrange = true;
   for( auto Dnit = Dltdn.begin() ;  ; ) {
    auto tit = Dnit++;
    if( Dnit == Dltdn.end() )
     break;
    if( *Dnit != (*tit) + 1 ) {
     isrange = false;
     break;
     }
    }

   if( isrange ) {  // it actually was a range
    auto rng = Range( Dltdn.front() , Dltdn.back() );

    // remove the variables in the LagBFunction
    for( Index h = 0 ; h < f_nsb ; )
     v_LBF[ h++ ]->remove_variables( rng , mp );

    // remove the variables in the Objective
    static_cast< p_LF >( static_cast< p_FRO >( LagrDual->get_objective()
					       )->get_function()
			 )->remove_variables( rng , mp );

    // shift range so that it is in [ 0 , n. dynamic constraints )
    rng.first -= static_cons;
    rng.second -= static_cons;

    // adjust the index to dynamic constraint dictionary
    std::copy( idx_to_dcon.begin() + rng.second , idx_to_dcon.end() ,
	       idx_to_dcon.begin() + rng.first );

    // now actually remove the dynamic variable
    LagrDual->remove_dynamic_variables( *Ld , rng );
    }
   else {  // it was a generic subset
    // remove the variables in the LagBFunction (copy the names)
    for( Index h = 0 ; h < f_nsb ; )
     v_LBF[ h++ ]->remove_variables( Subset( Dltdn ) , true , mp );

    // remove the variables in the Objective (give away the names)
    static_cast< p_LF >( static_cast< p_FRO >( LagrDual->get_objective()
					       )->get_function()
			 )->remove_variables( std::move( Dltdn ) , true , mp );

    // adjust the index to dynamic constraint dictionary
    // shift names so that they are in [ 0 , n. dynamic constraints )
    for( auto & el : Dltdn )
     el -= static_cons;

    // nullptr-mark the element to delete
    for( auto & el : Dltdn )
     idx_to_dcon[ el ] = nullptr;

    // now remove() them
    std::remove( idx_to_dcon.begin() , idx_to_dcon.end() , nullptr );

    // now actually remove the dynamic variable
    LagrDual->remove_dynamic_variables( *Ld , std::move( Dltdn ) );
    }

   Index i = static_cons;
   NumVar -= Dltds.size();
   Index nl = NumVar - i;
   idx_to_dcon.resize( nl );

   // now adjust the dynamic constraint to index dictionary, that is,
   // rebuid it anew based on the index to constraint one and re-sort it
   dcon_to_idx.resize( nl );
   for( Index j = 0 ; j < nl ; ++j )
    dcon_to_idx[ j ] = std::pair( idx_to_dcon[ j ] , i++ );

   std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

   }  // end( multiple constraints )
  
  LagrDual->close_channel( chnl );  // close the channel
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now actually add all Lagrangian variables to all the LagBFunction and
 // the LinearFunction: ensure that all the corresponding Modification are
 // bunched into a unique GroupModification of the LagrDual

 if( ! Addd.empty() ) {
  // the list to which the variable have to be added
  auto Ld = LagrDual->get_dynamic_variable< ColVariable >( "Lambda_d" );
  auto NAddd = Addd.size();

  // the list of variables to be added
  std::list< ColVariable > NLd( NAddd );

  // the iterator to the variables to be added (remains valid)
  auto NLDLit = NLd.begin();
  auto Lit = NLDLit;  // working copy

  // actually add the dynamic variable
  LagrDual->add_dynamic_variables( *Ld , NLd );

  // construct and set the new Lagrangian terms and new objective piece
  v_coeff_pair objcf( NAddd );
  auto objit = objcf.begin();

  std::vector< std::vector< v_coeff_pair > > LagTerms( Addd.size() );
  auto LTit = LagTerms.begin();

  Index i = NumVar;
  NumVar += NAddd;
  dcon_to_idx.resize( NumVar - static_cons );
  idx_to_dcon.resize( NumVar - static_cons );
  auto dc2iit = dcon_to_idx.begin() + i;
  auto i2dcit = idx_to_dcon.begin() + i;

  for( auto el : Addd ) {
   // check the LHS/RHS
   auto lhs = el->get_lhs();
   auto rhs = el->get_rhs();

   if( ( ( lhs == -INFshift ) && ( rhs == INFshift ) ) || el->is_relaxed() ) {
    // this constraint is eiter "infinitely loose" or relaxed: its rhs is
    // 0 and the Lagrangian term is empty
    *(objit++) = std::make_pair( & *(Lit++) , 0 );
    ++LTit;
    return;
    }

   auto coef = constr2val( *el , *Lit );

   // write the coefficient in the objective
   *(objit++) = std::make_pair( & *(Lit++) , coef );

   // split the linear constraint among the sub-Block
   split_constraint( *el , *(LTit++) );
   }

  // re-sort the dynamic constraints-->Lagrangian-variables dictionary
  std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );
  
  // open a channel where to bunch all the addition Modifications
  const auto chnl = LagrDual->open_channel();
  const auto mp = Observer::make_par( eModBlck , chnl );

  // add the variables to the objective of the Lagrangian Dual
  static_cast< p_LF >( static_cast< p_FRO >( LagrDual->get_objective()
					     )->get_function()
		       )->add_variables( std::move( objcf ) , mp );

  // add the variables in the LagBFunctions
  for( Index h = 0 ; h < f_nsb ; ++h ) {
   Lit = NLDLit;
   v_dual_pair dp( NAddd );  // construct the dual pairs
   for( Index i = 0 ; i < NAddd ; ++i ) {
    dp[ i ].first = & (*Lit++);
    dp[ i ].second = new LinearFunction( std::move( LagTerms[ i ][ h ] ) );
    }

   v_LBF[ h ]->add_dual_pairs( std::move( dp ) );
   }

  LagrDual->close_channel( chnl );  // close the channel
  }

 // and now, finally, all is done- - - - - - - - - - - - - - - - - - - - - - -
  
 }  // end( LagrangianDualSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------- End File LagrangianDualSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
