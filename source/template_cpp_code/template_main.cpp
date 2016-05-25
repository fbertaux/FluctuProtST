/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "HybridSimulator.hpp"


int main ()
{
	// loading model parameters
	ModelParameters* modelParameters = new ModelParameters () ;

    // construction of simulator
    HybridSimulator* hybridSimulator = new HybridSimulator ( modelParameters ) ;

	// construction of a cell
	CellState *cell = new CellState (modelParameters) ;
    hybridSimulator->simulate ( cell , 7.*24. ) ;

    // apply a stimulus and simulate the cell response
    // cell->set_NameOfSpecies_Level (0.1) ;
    hybridSimulator->simulate ( cell , 12. ) ;

    // check the result
    // cout << "Level of _NameOfSpecies_ after stimulus = " << cell->get_NameOfSpecies_Level () << endl ;

    // free memory used by the cell if not needed
    delete cell ;

    return 0 ;
}
