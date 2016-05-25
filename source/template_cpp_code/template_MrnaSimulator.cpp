/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "MrnaSimulator.hpp"

// additionally needed function
void sparmatfill(NRvector<NRsparseCol> &sparmat, MatDoub &fullmat)
{
    Int n,m,nz,nn=fullmat.nrows(),mm=fullmat.ncols();
    if (sparmat.size() != mm) throw("bad sizes");
    for (m=0;m<mm;m++)
    {
        for (nz=n=0;n<nn;n++) if (fullmat[n][m]) nz++;
        sparmat[m].resize(nn,nz);
        for (nz=n=0;n<nn;n++) if (fullmat[n][m])
        {
            sparmat[m].row_ind[nz] = n;
            sparmat[m].val[nz++] = fullmat[n][m];
        }
    }
}


// constructor
MrnaSimulator::MrnaSimulator ( ModelParameters* modelParameters , Int seed )
    : mf_ModelParameters(modelParameters), ran(seed),
    mm(4*mf_ModelParameters->mf_NumGenes), nn(3*mf_ModelParameters->mf_NumGenes),
    a(mm,0.), outchg(mm), depend(mm), pr(mm), t(0.),asum(0.),
    dispatch(new rateptr[mm])
{
    Int i,j,k,d;
    describereactions();
    sparmatfill(outchg,outstate);
    MatDoub dep(mm,mm);
    for (i=0;i<mm;i++) for (j=0;j<mm;j++)
    {
        d = 0;
        for (k=0;k<nn;k++) d = d || (instate[k][i] && outstate[k][j]);
        dep[i][j] = d;
    }
    sparmatfill(depend,dep);
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
    }
}


// preparation
void
MrnaSimulator::prepareForSteps ( CellState* cellState )
{
    t=0;
    Int i;
    asum = 0.;
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
        a[i] = (this->*dispatch[i])(cellState->mf_GeneMrnas);
        asum += a[i];
    }
}

// simulation
Doub
MrnaSimulator::doStep ( CellState* cellState , Doub targetTime )
{
    Int i,n,m,k=0;
    Doub tau,atarg,sum,anew;
    if (asum == 0.) {t = targetTime; return t;}
    tau = -log(ran.doub())/asum;
    if (t+tau>targetTime)
    {
        t=targetTime;
        return t;
    }
    atarg = ran.doub()*asum;
    sum = a[pr[0]];
    while (sum < atarg) sum += a[pr[++k]];
    m = pr[k];
    if (k > 0) SWAP(pr[k],pr[k-1]);
    if (k == mm-1) asum = sum;
    n = outchg[m].nvals;
    for (i=0;i<n;i++)
    {
        k = outchg[m].row_ind[i];
        cellState->mf_GeneMrnas[k] += outchg[m].val[i];
    }
    n = depend[m].nvals;
    for (i=0;i<n;i++)
    {
        k = depend[m].row_ind[i];
        anew = (this->*dispatch[k])(cellState->mf_GeneMrnas);
        asum += (anew - a[k]);
        a[k] = anew;
    }
    if (t*asum < 0.1)
        for (asum=0.,i=0;i<mm;i++) asum += a[i];
    return (t += tau);
}

// sample a cell from steady-state
CellState*
MrnaSimulator::sampleFromOnlyNativeSteadyState (Doub toReachSteadyStateDuration)
{
    CellState* cell = new CellState (mf_ModelParameters) ;
    prepareForSteps (cell) ;
    Doub dt,oldt,newt;
    while (t<toReachSteadyStateDuration)
    {
        oldt = t ;
        newt = doStep ( cell , toReachSteadyStateDuration ) ;
        dt = newt - oldt ;
        for (Int i=0;i<mf_ModelParameters->mf_NumGenes;i++)
        {
            cell->mf_AllProts[i] = mf_ModelParameters->mf_ksps[i]*cell->mf_GeneMrnas[3*i+2]/mf_ModelParameters->mf_rps[i]+(cell->mf_AllProts[i]-mf_ModelParameters->mf_ksps[i]*cell->mf_GeneMrnas[3*i+2]/mf_ModelParameters->mf_rps[i])*exp(-mf_ModelParameters->mf_rps[i]*dt);
        }
    }
    return cell ;
}


// others function, notably for init structure ( called by constructor )
void
MrnaSimulator::describereactions ()
{
	instate = MatDoub ( nn , mm , 0. ) ;
	Int m = 0 ;
	for ( Int g=0 ; g<mf_ModelParameters->mf_NumGenes ; g++ )
	{
		instate [3*g+0][m] = 1. ; m++ ;
		instate [3*g+1][m] = 1. ; m++ ;
		instate [3*g+0][m] = 1. ; m++ ;
		instate [3*g+2][m] = 1. ; m++ ;
	}

	outstate = MatDoub ( nn , mm , 0. ) ;
	m = 0 ;
	for ( Int g=0 ; g<mf_ModelParameters->mf_NumGenes ; g++ )
	{
		outstate [3*g+0][m] = -1. ;
		outstate [3*g+1][m] = 1. ; m++ ;

		outstate [3*g+0][m] = 1. ;
		outstate [3*g+1][m] = -1. ; m++ ;

		outstate [3*g+2][m] = 1. ; m++ ;
		outstate [3*g+2][m] = -1. ; m++ ;
	}

	// for dispatching
placeholder_dispatch
}
