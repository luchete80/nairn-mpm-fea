/********************************************************************************
	Expression.hpp
	NairnMPMFEA
 
 	Created by John Nairn on 8/7/18.
 	Copyright (c) 2018 John A. Nairn, All rights reserved.
	
 	A recursive definition of an expression is
 
 		Expression = sign (Expression) op function(Expression,Expression,...) ...
 
 	In evaluations, all Expressions in groupings ( (...) ) are recursively evaluated
	until them become a function or an Atomic. Similar, function arguments
	are recursively evaluated until they are Atomics. In final state, Expression
 	is reduced to
 
 		Expression = sign Atomic op function(Atomic,Atomic...) op ...
 
 		sign is leading + or -
 
 		Atomic is
 			Unsigned Number (23, 2.234, 3e2, 4.56e-3, ...)
 			Variable (varname - start with letter, contain letters or _)
 
		op is + - * / ^ and groupings ( , )
 
 		functions are created in Atomic::FindFunctionCode()
 			and evaulated in DoFunction()
 
 		Expression formed with following operator precedence
 			Variables replaced buy numerical values
 			group evaluation (including function arguments)
 			functions evaluated
 			^ (left to right)
 			sign applied to first term
 			* and / (left to right)
 			+ and - (left to right)
	
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Read_XML/Expression.hpp"
#include "Read_XML/Atomic.hpp"
#include "Exceptions/CommonException.hpp"
#include "System/UnitsController.hpp"

// only here to help XCODE dynamically compile
#include <string>
#include <unordered_map>
#include <cstring>
#include <cmath>
#include "System/MPMPrefix.hpp"

// extra function defined at the bottom
double erfcc(double x);

// To print results during function evaluation
//#define DEBUGEXPR

// globals
Expression *exfxn[12]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

// All variables in faster calls (no use of unordered map)
// Most used for MPM functions need while running (for speed)
// In this mode, only use first letter on variable, must be unique, rest ingnored
// Supported variables are (not hard to add more when needed):
// t:1, 2:x, 3:y, 4:z, 5:dt 6:q
// Create vars[7], set needed values, set var[0] last index+0.5
//                     A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,[,\,],^,_,`,
static short vmap[58]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	                   0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,1,0,0,0,2,3,4};
//                     a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z

#pragma mark CONSTRUCTORS AND DESTRUCTORS

// string with expression
Expression::Expression(const char *s)
{
	exprStr = NULL;
	firstAtom = NULL;
	numAtoms = 0;
	SetString(s);
}

// Destructor
Expression::~Expression()
{
	if(exprStr!=NULL) delete [] exprStr;
	
	Atomic *nextAtom = firstAtom;
	while(nextAtom!=NULL)
	{	Atomic *tempAtom = nextAtom->GetNextAtom();
		delete nextAtom;
		nextAtom = tempAtom;
	}
}

#pragma mark METHODS

// Tokenize expression (but only called once)
void Expression::TokenizeExpr(void)
{
	size_t exprLength = strlen(exprStr);
	bool aNum = true;
	bool hasDecimal = false;
	bool hasExpo = false;
	
	// empty string, return empty tokens array
	if(exprLength==0)
	{	throw CommonException("The expression is an empty string",exprStr);
	}
	
	// last operator
	unsigned int operLoc=0;
	
	// look for initial sign only
	char exprChar=exprStr[0];
	if(exprChar=='+' || exprChar=='-')
	{	AddToken(new Atomic(exprChar));
		operLoc++;
		if(exprLength==1)
			throw CommonException("An expression with just a plus or minus sign",exprStr);
	}
	
	// search for all tokens
	ExprRange expr;
	expr.location = operLoc;
	while(operLoc<exprLength)
	{	// next character
		exprChar = exprStr[operLoc];
		
		// an operator found
		if(exprChar=='+' || exprChar=='-' || exprChar=='*' || exprChar=='/' ||
		   exprChar=='^' || exprChar=='(' || exprChar==',' || exprChar==')')
		{   expr.length = operLoc - expr.location;
			bool openGroup = exprChar=='(' ? true : false ;
			AddSubstringToken(expr,aNum,openGroup);
			
			// add the operator
			AddToken(new Atomic(exprChar));
			
			// reset for next element
			expr.location=operLoc+1;
			aNum = true;
			hasDecimal = false;
			hasExpo = false;
		}
		
		// possible exponential operator (number - e - digit or number - e - sign)
		//   (note 1e- or 1e+ will ignore e- and e+ or be like 1e+0 or 1e-0)
		else if(exprChar=='e' && aNum && !hasExpo &&
						operLoc>expr.location && operLoc<exprLength-1)
		{   exprChar=exprStr[operLoc+1];
			if(exprChar>='0' && exprChar<='9')
			{	// has number after e
				hasExpo = true;
				operLoc++;
			}
			else if(exprChar=='+' || exprChar=='-')
			{	// has + or - after e, must have number next
				if(operLoc+2<exprLength)
				{   exprChar=exprStr[operLoc+2];
					if(exprChar>='0' && exprChar<='9')
					{	hasExpo = true;
						operLoc += 2;
					}
					else
						aNum = false;
				}
				else
					aNum = false;
			}
			else
				aNum=false;
		}
		
		// see if current block is still an unsigned number
		else if(exprChar=='.')
		{   if(hasDecimal || hasExpo)
				aNum = false;
			else
				hasDecimal = true;
		}
		
		// still unsigned number?
		else if(exprChar<'0' || exprChar>'9')
			aNum = false;
		
		// next character
		operLoc++;
	}
	
	// last one
	expr.length = operLoc - expr.location;
	AddSubstringToken(expr,aNum,false);
	
	// check for some errors
	Atomic *prevAtom = firstAtom;
	if(prevAtom==NULL)
		throw CommonException("Expression with no valid tokens",exprStr);
	
	// Check code
	int prevCode = prevAtom->GetCode();
	if(prevCode==OP_CLOSE_GROUP)
		throw CommonException("Expression begins in ')'",exprStr);
	
	// check each subsequent object
	int groupDepth;
	if(prevCode==OP_OPEN_GROUP)
	{	groupDepth = 1;
		exprHasGroups = true;
	}
	else
	{	groupDepth = 0;
		exprHasGroups = false;
	}
	
	numAtoms = 1;
	prevAtom = prevAtom->GetNextAtom();
	while(prevAtom!=NULL)
	{	// check if this code can follow previous code
		int nextCode = prevAtom->GetCode();
		ValidateCodeOrder(prevCode,nextCode);
		
		// track groups
		prevCode = nextCode;
		if(prevCode==OP_OPEN_GROUP)
		{	groupDepth++;
			exprHasGroups = true;
		}
		else if(prevCode==OP_CLOSE_GROUP)
			groupDepth--;
		else if(prevCode==OP_DIVIDE_GROUP)
		{	if(groupDepth<1)
				throw CommonException("Comma not inside function arguments",exprStr);
		}
		if(groupDepth<0)
			throw CommonException("Unmatched parentheses in expression",exprStr);
		
		// next one
		numAtoms++;
		prevAtom = prevAtom->GetNextAtom();
	}
	
	// must end with match groups and cannot end in operator
	if(groupDepth!=0)
		throw CommonException("Unmatched parentheses in expression",exprStr);
	else if(prevCode>=OP_PLUS)
		throw CommonException("Expression ends in an operator",exprStr);
}

// evaluate current expression with set a variables
// When called in parallel, must trap exceptions in that thread
double Expression::EvaluateFunction(double *usevar) const
{
	// must tokenize before call this const method
	if(firstAtom==NULL)
		throw CommonException("Expression has no tokens",exprStr);

	// Copy all atomics and replace variables with their values
	int maxVar = (int)usevar[0];
	Atomic *tmpFirstAtom = NULL,*nextAtom = NULL;
	Atomic *sourceAtom = firstAtom;
	while(sourceAtom!=NULL)
	{	// copy Atomic to new list
		Atomic *copiedAtom = sourceAtom->GetCopy(true);
		if(tmpFirstAtom==NULL)
			tmpFirstAtom = copiedAtom;
		else
			nextAtom->SetNextAtom(copiedAtom);
		nextAtom = copiedAtom;
		
		if(nextAtom->GetCode()==ATOM_VARIABLE)
		{	// get the value (caller better provide them all)
			int mapID = vmap[nextAtom->GetVarID()];
			if(mapID > maxVar || mapID<1)
				throw CommonException("Expression has an undefined variable",nextAtom->GetVarName());
			double varValue = usevar[mapID];
			
			// change atomic to number
			nextAtom->SetCode(ATOM_NUMBER);
			nextAtom->SetValue(varValue);
		}
		
		// next one
		sourceAtom = sourceAtom->GetNextAtom();
	}
	
	// evalute expression and delete the last copy
	double varValue=EvaluateTokens(tmpFirstAtom,exprHasGroups);
	delete tmpFirstAtom;
	return varValue;
}

// evaluate current expression with set a variables
// When called in parallel, must trap exceptions in that thread
double Expression::EvaluateFunction(unordered_map<string, double> usevar) const
{
	// must tokenize before call this const method
	if(firstAtom==NULL)
		throw CommonException("Expression has no tokens",exprStr);
	
	// Copy all atomics and replace variables with their values
	Atomic *tmpFirstAtom = NULL,*nextAtom = NULL;
	Atomic *sourceAtom = firstAtom;
	while(sourceAtom!=NULL)
	{	// copy Atomic to new list
		Atomic *copiedAtom = sourceAtom->GetCopy(true);
		if(tmpFirstAtom==NULL)
			tmpFirstAtom = copiedAtom;
		else
			nextAtom->SetNextAtom(copiedAtom);
		nextAtom = copiedAtom;
		
		if(nextAtom->GetCode()==ATOM_VARIABLE)
		{	// get the value
			string key = nextAtom->GetVarName();
			if(usevar.find(key) == usevar.end())
				throw CommonException("Expression has an undefined variable",nextAtom->GetVarName());
			double varValue = usevar[key];
			
			// change atomic to number
			nextAtom->SetCode(ATOM_NUMBER);
			nextAtom->SetValue(varValue);
		}
		
		// next one
		sourceAtom = sourceAtom->GetNextAtom();
	}
	
	// evalute expression and delete the last copy
	double varValue=EvaluateTokens(tmpFirstAtom,exprHasGroups);
	delete tmpFirstAtom;
	return varValue;
}

// evaluate tokens starting at callFirstAtom
double Expression::EvaluateTokens(Atomic *callFirstAtom,bool hasGroups) const
{
	Atomic *atom,*prevAtom,*startAtom,*nextAtom;
	int groupDepth,opCode;

#ifdef DEBUGEXPR
	cout << "***** Evaluate Tokens with hasGroups = " << hasGroups << ":";
	Describe(callFirstAtom);
#endif
	
	// parce groups recursively
	if(hasGroups)
	{	atom = callFirstAtom;
		while(atom)
		{   if(atom->GetCode()==OP_OPEN_GROUP)
			{   // find matching close group
				prevAtom = atom;
				startAtom = nextAtom = atom->GetNextAtom();
				groupDepth = 1;
				hasGroups = false;
				while(nextAtom)
				{   opCode = nextAtom->GetCode();
					if(opCode==OP_OPEN_GROUP)
					{	groupDepth++;
						hasGroups = true;
					}
					else if(opCode==OP_CLOSE_GROUP)
						groupDepth--;
					if((opCode==OP_CLOSE_GROUP && groupDepth==0) || (opCode==OP_DIVIDE_GROUP && groupDepth==1))
					{	// evaluate subexpression within group
						// group is atom=( or , : subexpression in startAtom to prevAtom : nextAtom=) or ,
						prevAtom->SetNextAtom(NULL);			// mark as end of subexpression
						double exprValue = EvaluateTokens(startAtom,hasGroups);
						
						// replace previous group open token with new value
						atom->SetCode(ATOM_NUMBER);
						atom->SetValue(exprValue);

						// current subexpression from startAtom to prevAtom (the close)
						// and prevAtom has changed its next atom to NULL
						// delete startAtom to prevAtom
						while(startAtom!=NULL)
						{	prevAtom = startAtom->GetNextAtom();
							delete startAtom;
							startAtom = prevAtom;
						}

						// all done if was close group
						if(groupDepth==0)
						{	// atom, which was open group, now the expression value and it
							// should be followed by atom after the close group in nextAtom
							atom->SetNextAtom(nextAtom->GetNextAtom());		// To outside the grounp
							
							// delete close group in nextAtom too
							delete nextAtom;
							
							break;
						}
						
						// nextAtom is a divider or , so go on to next subgroup
						atom->SetNextAtom(nextAtom);			// The initial ( (now a value) links to the , --- which will soon be value
						atom = nextAtom;						// The , becomes the new open group
						startAtom = atom->GetNextAtom();		// First atom in next subexpression is one after ,
					}
					
					// on to next element
					prevAtom = nextAtom;
					nextAtom = prevAtom->GetNextAtom();
				}
			
				// did not find the end group
				if(nextAtom == NULL)
					throw CommonException("Unbalanced parentheses in the expression",exprStr);
			}
			atom = atom->GetNextAtom();
		}
	}

#ifdef DEBUGEXPR
	cout << "***** Finished groups: ";
	Describe(callFirstAtom);
#endif
	
	// evaluate functions
	atom = callFirstAtom;
	while(atom)
	{	if(atom->GetCode()==FUNCTION_NAME)
		{	// Evaluates the function, put in atom, which then points to one after arguments
			DoFunction(atom);
		}
		atom = atom->GetNextAtom();
	}
#ifdef DEBUGEXPR
	cout << "***** Functions Done: ";
	Describe(callFirstAtom);
#endif
	
	// evaluate powers
	opCode = callFirstAtom->GetCode();
	atom = (opCode==OP_PLUS || opCode==OP_MINUS) ? callFirstAtom->GetNextAtom() : callFirstAtom;
	OperateTerms(atom,OP_POW,-1);
#ifdef DEBUGEXPR
	cout << "***** ^ Done: ";
	Describe(callFirstAtom);
#endif
	
	// evalate unary sign (using opCode from above)
	if(opCode==OP_PLUS || opCode==OP_MINUS)
	{	nextAtom = callFirstAtom->GetNextAtom();
		if(nextAtom->GetCode()!=ATOM_NUMBER)
			throw CommonException("Negative sign at start of expression not valid",exprStr);
		double firstValue = nextAtom->GetValue();
		if(opCode==OP_MINUS) firstValue = -firstValue;
		
		// convert first item to number (can't skip or delete because caller needs it to remain)
		callFirstAtom->SetValue(firstValue);
		callFirstAtom->SetCode(ATOM_NUMBER);
		
		// skip over the previoius number
		callFirstAtom->SetNextAtom(nextAtom->GetNextAtom());
		delete nextAtom;

#ifdef DEBUGEXPR
		cout << "***** Uniary Sign Done: ";
		Describe(callFirstAtom);
#endif
	}
	
	
	// evaluate * and /
	OperateTerms(callFirstAtom,OP_TIMES,OP_DIV);
#ifdef DEBUGEXPR
	cout << "***** * and / Done: ";
	Describe(callFirstAtom);
#endif
	
	// evaluate + and -
	OperateTerms(callFirstAtom,OP_PLUS,OP_MINUS);
#ifdef DEBUGEXPR
	cout << "***** + and - Done: ";
	Describe(callFirstAtom);
#endif

	// all done, must be just one left
	if(callFirstAtom->GetNextAtom()!=NULL)
		throw CommonException("Expression ended with unevaluated tokens",exprStr);
	
#ifdef DEBUGEXPR
	cout << "***** Evaluated to Value:";
	cout << "Value = " << callFirstAtom->GetValue() << endl;;	// all done, must be just one left
#endif
	
	return callFirstAtom->GetValue();
}

// The function is in atom and consective ATOM_NUMBER atomics are arguments
// Evaulate function, set atom to value of function
// Delete arguments and point atom to first one after the arguments
void Expression::DoFunction(Atomic *atom) const

{
	// get first argument (need special case if add zero-argument function
	int fxnCode = atom->GetFunctionCode();
	Atomic *arg = atom->GetNextAtom();
	if(arg==NULL)
		throw CommonException("Function does not have a numeric argument",atom->GetFunctionName());
	if(arg->GetCode()!=ATOM_NUMBER)
		throw CommonException("Function does not have a numeric argument",atom->GetFunctionName());
	double arg1 = arg->GetValue();
	
	// If need more arguments, delete current arg, point to next, and get it
	double fxnValue = 0.;			// default assumed in some function (keep as zero)
	switch(fxnCode)
	{	case SIN_FXN:
			fxnValue = sin(arg1);
			break;
		case COS_FXN:
			fxnValue = cos(arg1);
			break;
		case TAN_FXN:
			fxnValue = tan(arg1);
			break;
		case ASIN_FXN:
			fxnValue = asin(arg1);
			break;
		case ACOS_FXN:
			fxnValue = acos(arg1);
			break;
		case ATAN_FXN:
			fxnValue = atan(arg1);
			break;
		case SINH_FXN:
			fxnValue = sinh(arg1);
			break;
		case COSH_FXN:
			fxnValue = cosh(arg1);
			break;
		case TANH_FXN:
			fxnValue = tanh(arg1);
			break;
		case LOG_FXN:
			fxnValue = log(arg1);
			break;
		case LOG10_FXN:
			fxnValue = log10(arg1);
			break;
		case ABS_FXN:
			fxnValue = fabs(arg1);
			break;
		case INT_FXN:
			fxnValue = floor(arg1);
			break;
		case SQRT_FXN:
			fxnValue = sqrt(arg1);
			break;
		case SIGN_FXN:
			if(arg1>0.) fxnValue = 1.;
			break;
		case EXP_FXN:
			fxnValue = exp(arg1);
			break;
		case RAND_FXN:
			fxnValue = arg1*Random();		// i.e., (0,arg1) +/- halfRandomBoxSize
			break;
		case ERF_FXN:
			fxnValue = 1.-erfcc(arg1);
			break;
		case ERFC_FXN:
			fxnValue = erfcc(arg1);
			break;
		case COSRAMP_FXN:
		{	// function = 0.5*arg1*(1-cos(pi*arg2)) when arg2 between 0 and 1
			double arg2 = 0.;
			arg = GetFunctionArg(atom,arg,arg2);
			if(arg2>=1.)
				fxnValue = arg1;
			else if(arg2>0.)
				fxnValue = 0.5*arg1*(1.-cos(PI_CONSTANT*arg2));
			break;
		}
		case RAMP_FXN:
		{	// function = arg1*arg2 when arg2 between 0 and 1
			double arg2 = 0.;
			arg = GetFunctionArg(atom,arg,arg2);
			if(arg2>=1.)
				fxnValue = arg1;
			else if(arg2>0.)
				fxnValue = arg1*arg2;
			break;
		}
		case BOX_FXN:
		{	// function = arg1 when arg2 between 0 and 1
			double arg2 = 0.;
			arg = GetFunctionArg(atom,arg,arg2);
			if(arg2>0. && arg2<1.)
				fxnValue = arg1;
			break;
		}
		case SINBOX_FXN:
		{	// function = arg1*sin(pi*arg2) when arg2 between 0 and 1
			double arg2 = 0.;
			arg = GetFunctionArg(atom,arg,arg2);
			if(arg2>0. && arg2<1.)
				fxnValue = arg1*sin(PI_CONSTANT*arg2);
			break;
		}
        case MOD_FXN:
        {    // function = mod(arg1,arg2)
            double arg2 = 1.;
            arg = GetFunctionArg(atom,arg,arg2);
            fxnValue = fmod(arg1,arg2);
            break;
        }
		case SGN_FXN:
			if(arg1<0.)
				fxnValue = -1.;
			else if(arg1>0.)
				fxnValue = 1.;
			break;
		case TRI_FXN:
		{
			fxnValue = 1. - fabs(arg1);
			fxnValue = fmax(0, fxnValue);
			break;
		}
        case NORMALCDFINVERSE_FXN:
            fxnValue = NormalCDFInverse(arg1);
            break;
		default:
			// will be zero
			break;
	}
	
	// point to one after the argument(s) and delete argument
	Atomic *nextAtom = arg->GetNextAtom();
	delete arg;
	
	if(nextAtom!=NULL)
	{	if(nextAtom->GetCode()==ATOM_NUMBER)
			throw CommonException("Function does has too many numeric arguments",atom->GetFunctionName());
	}
	
	// set function atom to its value
	atom->SetCode(ATOM_NUMBER);
	atom->SetValue(fxnValue);
	
	// atom point to atom after argument
	atom->SetNextAtom(nextAtom);
}

// Group next atom as function argument and find its value
// Also delete current arg atom and return new one
Atomic *Expression::GetFunctionArg(Atomic *atom,Atomic *arg,double &argValue) const
{	// value next argument
	Atomic *tmpArg = arg->GetNextAtom();
	if(tmpArg==NULL)
		throw CommonException("Function does not have a numeric argument",atom->GetFunctionName());
	if(tmpArg->GetCode()!=ATOM_NUMBER)
		throw CommonException("Function does not have a numeric argument",atom->GetFunctionName());
	
	// get value
	argValue = tmpArg->GetValue();

	// delete prior (no longer needed)
	delete arg;
	
	// return new one
	return tmpArg;
}

// operate from firstOpAtom doing one or two types of numeric operations with equal precedence
void Expression::OperateTerms(Atomic *firstOpAtom,int op1,int op2) const
{
	Atomic *prevAtom = firstOpAtom,*nextAtom;
	Atomic *atom = firstOpAtom->GetNextAtom();
	int opCode;
	double prevNumber,nextNumber;
	
		while(atom)
	{	opCode = atom->GetCode();
		if(opCode==op1 || opCode==op2)
		{	// get left argument
			if(prevAtom->GetCode()!=ATOM_NUMBER)
				throw CommonException("Invalid argument to left of an operator","");
			prevNumber = prevAtom->GetValue();
			
			// get right argument
			nextAtom = atom->GetNextAtom();
			if(nextAtom==NULL)
				throw CommonException("Invalid argument to right of an operator","");
			if(nextAtom->GetCode()!=ATOM_NUMBER)
				throw CommonException("Invalid argument to right of an operator","");
			nextNumber = nextAtom->GetValue();
			
			// store result in prevAtom
			prevAtom->SetCode(ATOM_NUMBER);
			switch(opCode)
			{	case OP_TIMES:
					prevAtom->SetValue(prevNumber*nextNumber);
					break;
				case OP_DIV:
					prevAtom->SetValue(prevNumber/nextNumber);
					break;
				case OP_PLUS:
					prevAtom->SetValue(prevNumber+nextNumber);
					break;
				case OP_MINUS:
					prevAtom->SetValue(prevNumber-nextNumber);
					break;
				case OP_POW:
					prevAtom->SetValue(pow(prevNumber,nextNumber));
					break;
				default:
					prevAtom->SetValue(0.);
					break;
			}
			
			// to to past right argument
			prevAtom->SetNextAtom(nextAtom->GetNextAtom());

			// delete operator and right argument
			delete atom;
			delete nextAtom;
			
			// point to prior result
			atom = prevAtom;
		}
		
		// on to next one
		prevAtom = atom;
		atom = prevAtom->GetNextAtom();
			
	}
}

// Add new atom of the list of tokens
void Expression::AddToken(Atomic *nextAtom)
{
	if(firstAtom==NULL)
		firstAtom = nextAtom;
	else
		currentAtom->SetNextAtom(nextAtom);
	currentAtom = nextAtom;
}

// Add current string as token. The only options are
//  1. If before group then must be a valid function
//	2. a valid number
//	3. "pi" for that number
//  4. valid variable name
// Otherwise and error
void Expression::AddSubstringToken(ExprRange expr,bool aNum,bool beforeGroup)
{
	// exit if nothing
	if(expr.length==0) return;
	
	// get substring
	char *subString = new char[expr.length+1];
	for(int i=0;i<expr.length;i++)
		subString[i] = exprStr[expr.location+i];
	subString[expr.length]=0;
	
	// create Atomic
	Atomic *nextAtom = NULL;
	if(beforeGroup)
	{	nextAtom = new Atomic(subString,FUNCTION_NAME);
		delete [] subString;
		if(nextAtom->GetFunctionCode()<0)
			throw CommonException("Expression with invalid function name",subString);
	}
	else if(aNum)
	{	double value;
		sscanf(subString,"%lf",&value);
		nextAtom = new Atomic(value);
		delete [] subString;
	}
	else if(strcmp(subString,"pi")==0 || strcmp(subString,"Pi")==0 || strcmp(subString,"PI")==0)
	{	nextAtom = new Atomic(PI_CONSTANT);
		delete [] subString;
	}
	else
	{	// mush be valid variables name
		for(int i=0;i<expr.length;i++)
		{	char strChar = subString[i];
			
			// uppercase
			if(strChar>='A' && strChar<='Z') continue;
			
			// lowercase
			if(strChar>='a' && strChar<='z') continue;
			
			// underscore
			if(strChar=='_') continue;
			
			// number, but not first
			if(strChar>='0' && strChar<='9' && i>0) continue;
			
			throw CommonException("Expression has invalid variable name",subString);
		}
		
		// create Atomic (which is responsible to delete the substring)
		nextAtom = new Atomic(subString,ATOM_VARIABLE);
	}
	
	// Add to list
	AddToken(nextAtom);
}

// Check token pair to make sure it could be allowed
void Expression::ValidateCodeOrder(int prevCode,int nextCode) const
{
	switch(prevCode)
	{	case ATOM_NUMBER:
		case ATOM_VARIABLE:
			switch(nextCode)
			{	case ATOM_NUMBER:
				case ATOM_VARIABLE:
				case FUNCTION_NAME:
				case OP_OPEN_GROUP:
					throw CommonException("Invalid expression: number, variable, or ')' followed by number, variable, function, or '('",exprStr);
				default:
					// Allowed: +, -, *, /, ^, ), ','
					break;
			}
			break;
			
		case FUNCTION_NAME:
			// only defined before ( and therefore always allowed
			break;
			
		case OP_CLOSE_GROUP:
			switch(nextCode)
			{	case ATOM_NUMBER:
				case ATOM_VARIABLE:
				case FUNCTION_NAME:
				case OP_OPEN_GROUP:
				case OP_DIVIDE_GROUP:
					throw CommonException("Invalid expression: ')' followed by number, variable, function, ',', or '('",exprStr);
				default:
					// Allowed: +, -, *, /, ^, )
					break;
			}
			break;
			
		case OP_PLUS:
		case OP_MINUS:
		case OP_TIMES:
		case OP_DIV:
		case OP_POW:
			switch(nextCode)
			{	case OP_PLUS:
				case OP_MINUS:
				case OP_TIMES:
				case OP_DIV:
				case OP_POW:
				case OP_CLOSE_GROUP:
				case OP_DIVIDE_GROUP:
					throw CommonException("Invalid expression: +, -, *, /, or ^ followed by +, -, *, /, ^, ',', or ')'",exprStr);
				default:
					// Allowed num, var, fxn, (
					break;
			}
			break;
			
		case OP_OPEN_GROUP:
			switch(nextCode)
			{	case OP_TIMES:
				case OP_DIV:
				case OP_POW:
				case OP_CLOSE_GROUP:
				case OP_DIVIDE_GROUP:
					throw CommonException("Invalid expression: '(' followed by *, /, ^, ',', or ')'",exprStr);
				default:
					// Allowed: num, var, fxn, +, -, (,
					break;
			}
			break;
			
		case OP_DIVIDE_GROUP:
			switch(nextCode)
			{	case OP_TIMES:
				case OP_DIV:
				case OP_POW:
				case OP_CLOSE_GROUP:
					throw CommonException("Invalid expression: ',' followed by *, /, or ^",exprStr);
				default:
					// Allowed: num, var, fxn, +, -, (
					break;
			}
			break;
		default:
			break;
	}
}

// Evaluate function at specific values of x,y,z, and time (see Expression vmap)
double Expression::XYZTValue(Vector *vec,double etime) const
{
	double vars[5];
	vars[0] = 4.5;
	vars[1] = etime;		//t
	vars[2] = vec->x;		//x
	vars[3] = vec->y;		//y
	vars[4] = vec->z;		//z
	return EvaluateFunction(vars);
}

// Evaluate function at specific values time (see Expression vmap)
double Expression::TValue(double etime) const
{
	double vars[2];
	vars[0] = 1.5;
	vars[1] = etime;		//t
	return EvaluateFunction(vars);
}

#pragma mark ACCESSORS

// Set string and remove spaces
void Expression::SetString(const char *s)
{
	if(exprStr!=NULL) delete [] exprStr;
	exprStr = new char[strlen(s) + 1];
	strcpy(exprStr, s);
	
	// remove spaces
	char *last = exprStr;
	char *current = exprStr;
	while(*current!=0)
	{	if(*current!=' ')
		{	*last = *current;
			last++;
		}
		current++;
	}
	*last = 0;
}
const char *Expression::GetString(void) const { return exprStr; }

// For dubugging
void Expression::Describe(void) const
{
	Describe(firstAtom);
}

// For dubugging
void Expression::Describe(Atomic *tmpFirstAtom) const
{
	cout << "Expression: '" << exprStr << "'" << endl;
	Atomic *nextAtom = tmpFirstAtom;
	while(nextAtom!=NULL)
	{	nextAtom->Describe();
		nextAtom = nextAtom->GetNextAtom();
	}
}

#pragma mark CLASS METHODS

// Create new expression
// If error, throw SAX exception (so only use while reading XML file)
Expression *Expression::CreateExpression(char *exprString,const char *notValid)
{	Expression *function = NULL;
	try
	{	function = new Expression(exprString);
		function->TokenizeExpr();
	}
	catch(CommonException& err)
	{	char errMsg[200];
		strcpy(errMsg,notValid);
		strcat(errMsg,": ");
		strcat(errMsg,err.Message());
		ThrowSAXException(errMsg);
	}
	catch(...)
	{	ThrowSAXException(notValid);
	}
	return function;
}

// Create extra function of x,y,Z (or r=x and z=y) or D and T for polar coordinates from x and y
// Caller should DeleteFunction(i) when done using it
// These fill list of functions in global variables; cannot be in parallel
bool Expression::CreateFunction(char *&eqn,int i)
{
	// fails if base not there or index bad
	if(i<1 || i>MAX_EXPRESSIONS) return false;
	
	try
	{	exfxn[i-1] = new Expression(eqn);
		exfxn[i-1]->TokenizeExpr();
	}
	catch(CommonException& err)
	{	cout << "**** Expression error (" << err.Message() << ") in '" << eqn << "' ****" << endl;
		return false;
	}
	catch(...)
	{	cout << "**** Expression error in '" << eqn << "' ****" << endl;
		return false;
	}
	return true;
}

// If i<=0, delete all functions, otherwise delete function i
void Expression::DeleteFunction(int i)
{
	if(i>0 && exfxn[i-1]!=NULL)
	{	delete exfxn[i-1];
		exfxn[i-1] = NULL;
	}
	else
	{	for(int j=0;j<MAX_EXPRESSIONS && exfxn[j]!=NULL;j++)
		{	delete exfxn[j];
			exfxn[j] = NULL;
		}
	}
}

// get function value for (x,y,z) or (R,Z) or polar (D,T)
double Expression::FunctionValue(int i,double x,double y,double z,double xorig,double yorig,double zorig)
{
	unordered_map<string,double> vars;
	vars["x"] = x;
	vars["y"] = y;
	vars["z"] = z;
	
	// polar coordinates
	vars["R"] = z;
	vars["Z"] = y;
	x-=xorig;
	y-=yorig;
	double dvalue=sqrt(x*x+y*y),thetavalue;
	if(dvalue==0)
		thetavalue=0.;
	else if(x>=0.)
		thetavalue=asin(y/dvalue);
	else
		thetavalue=PI_CONSTANT-asin(y/dvalue);
	vars["D"] = dvalue;
	vars["T"] = thetavalue;
	
	// returns, but returns zero if error
	double result = 0.;
	if(exfxn[i-1]!=NULL)
	{	try
		{	result = exfxn[i-1]->EvaluateFunction(vars);
		}
		catch(...)
		{	result = 0.;
		}
	}
	return result;
}

// error function complement (Numerical Recipes in C, pg 176)
double erfcc(double x)
{	double z=fabs(x);
	double t=1.0/(1.0+0.5*z);
	double ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
						t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
						t*(-0.82215223+t*0.17087277)))))))));
	return x>=0. ? ans : 2.-ans;
}

