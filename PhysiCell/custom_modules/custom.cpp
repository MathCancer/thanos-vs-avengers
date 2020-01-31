/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

Cell_Definition civilian; 
Cell_Definition Thanos; 
Cell_Definition avenger;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 
	
	// add custom data here, if any 
	
	cell_defaults.custom_data.add_variable( "health" , "dimensionless", 1.0 );
	cell_defaults.custom_data.add_variable( "strength" , "dimensionless", 1.0 );
	cell_defaults.custom_data.add_variable( "attack_rate" , "1/min" , 1.0 ); 
	cell_defaults.custom_data.add_variable( "recovery_rate" , "1/min" , 0.01 ); 
	cell_defaults.custom_data.add_variable( "death_threshold" , "dimensionless" , 0.1 ); 
	
	cell_defaults.custom_data.add_variable( "attacking" , "dimensionless" , 0 ); 
	
	cell_defaults.custom_data.add_variable( "infinity_stones" , "dimensionless" , 0 ); 
	
	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	civilian = cell_defaults; 
	civilian.type = 1; 
	civilian.name = "civilian"; 
	
	// make sure the new cell type has its own reference phenotype
	
	civilian.parameters.pReference_live_phenotype = &( civilian.phenotype ); 
	
	// enable random motility 
	civilian.phenotype.motility.is_motile = true; 
	civilian.phenotype.motility.persistence_time = 1.0; 
	civilian.phenotype.motility.migration_speed = parameters.doubles( "civilian_speed" ); // 0.25 micron/minute 
	civilian.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// set birth rate
	civilian.phenotype.cycle.data.transition_rate(0,0) = 
		parameters.doubles( "civilian_birth_rate" ); // 0.0; 
	
	// set death rate 
	civilian.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "civilian_death_rate" ); // 0.0; 
	
	// Set cell-cell adhesion to 5% of other cells 
	civilian.phenotype.mechanics.cell_cell_adhesion_strength = 0;
	
	
	civilian.custom_data["strength"] = parameters.doubles("civilian_strength" ); 
	
	// Thanos setup 
	
	Thanos = civilian; 
	Thanos.type = 2; 
	Thanos.name = "Thanos"; 
	
	Thanos.phenotype.cycle.data.transition_rate(0,0) = 0; 
	
	Thanos.phenotype.motility.migration_speed = parameters.doubles("thanos_speed"); 
	Thanos.phenotype.motility.migration_bias = parameters.doubles("thanos_motility_bias");  
	Thanos.phenotype.motility.persistence_time = 1; 
	
	Thanos.functions.custom_cell_rule = Thanos_function; 
	Thanos.functions.update_migration_bias = avenger_taxis; 
	
	Thanos.custom_data["strength"] = parameters.doubles("thanos_strength" ); 
	
	int thanos_sig_i = microenvironment.find_density_index( "thanos" );
	Thanos.phenotype.secretion.secretion_rates[thanos_sig_i] = 10; 
	Thanos.phenotype.secretion.saturation_densities[thanos_sig_i] = 1; 
	
	// Avenger setup 
	
	avenger = Thanos; 
	avenger.type = 3; 
	avenger.name = "avenger";

	avenger.phenotype.cycle.data.transition_rate(0,0) = 0; 
	
	avenger.phenotype.motility.migration_speed = parameters.doubles("avenger_speed"); 
	avenger.phenotype.motility.migration_bias = parameters.doubles("avenger_motility_bias"); 
	avenger.phenotype.motility.persistence_time = 1; 
	
	avenger.functions.custom_cell_rule = avenger_function; 
	avenger.functions.update_migration_bias = thanos_taxis; 	

	avenger.custom_data["strength"] = parameters.doubles( "avenger_strength" ); 

	int avenger_sig_i = microenvironment.find_density_index( "avenger" );
	avenger.phenotype.secretion.secretion_rates[thanos_sig_i] = 0; 
	avenger.phenotype.secretion.saturation_densities[thanos_sig_i] = 1; 

	avenger.phenotype.secretion.secretion_rates[avenger_sig_i] = 10; 
	avenger.phenotype.secretion.saturation_densities[avenger_sig_i] = 1; 

	// Thanos is a big guy 
	Thanos.phenotype.volume.multiply_by_ratio( 4.0 ); 
	avenger.phenotype.volume.multiply_by_ratio( 2.0 ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	
	// get domain size 
	int X_nodes = microenvironment.mesh.x_coordinates.size(); 
	double X_left  = microenvironment.mesh.x_coordinates[0] - microenvironment.mesh.dx / 2.0 ; 
	double X_right = microenvironment.mesh.x_coordinates[X_nodes-1] + microenvironment.mesh.dx / 2.0 ; 
	double X_length = X_right - X_left; 
	
	int Y_nodes = microenvironment.mesh.y_coordinates.size(); 
	double Y_left  = microenvironment.mesh.y_coordinates[0] - microenvironment.mesh.dy / 2.0 ; 
	double Y_right = microenvironment.mesh.y_coordinates[Y_nodes-1] + microenvironment.mesh.dy / 2.0 ; 
	double Y_length = Y_right - Y_left; 
	
	int number_of_civilians = parameters.ints( "civilian_initial_count" ); 
	for( int n = 0; n < number_of_civilians ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double x = X_left + 0.1*X_length + 0.8*X_length*UniformRandom(); 
		double y = Y_left + 0.1*Y_length + 0.8*Y_length*UniformRandom(); 
		position[0] = x;
		position[1] = y; 
		
		pC = create_cell( civilian ); 
		pC->assign_position( position );
		
	}
	
	// seed avengers 
	
	int number_of_avengers = parameters.ints( "avenger_initial_count" ); 
	for( int n = 0; n < number_of_avengers ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double x = X_left + 0.1*X_length + 0.8*X_length*UniformRandom(); 
		double y = Y_left + 0.1*Y_length + 0.8*Y_length*UniformRandom(); 
		position[0] = x;
		position[1] = y; 
		
		pC = create_cell( avenger ); 
		pC->assign_position( position );
		
		// hide the 6 stones
		// these avengers are less directed towards Thanos 
		if( n < 6 )
		{
			pC->custom_data["infinity_stones"] = 1.0; 
			pC->phenotype.motility.migration_bias = 
				parameters.doubles("avenger_motility_bias_with_stone"); 
		} 
	}
	
	pC = create_cell( Thanos ); 
	std::vector<double> position = {0,0,0}; 
	pC->assign_position( position ); 

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start off with black 
	std::vector<std::string> output( 4 , "black" ); 
	
	// stay black if dead 
	if( pCell->phenotype.death.dead == true )
	{ return output; } 
		
	// Civilian cells are green 
	if( pCell->type == civilian.type )
	{
		output[0] = "limegreen"; 
		output[2] = "limegreen"; 
		output[3] = "limegreen"; 
		return output; 
	}
	
	// Thanos is purple 
	if( pCell->type == Thanos.type )
	{
		std::string color = "mediumpurple"; 
		
		if( pCell->custom_data["attacking"] > 0.5 )
		{ color = "lavender"; }

		output[0] = color;  
		output[2] = color; 
		output[3] = color; 
		
		double stones = pCell->custom_data["infinity_stones"]; 
		
		std::cout << std::endl << "\tThanos has " << (int) (stones) 
				  << " infinity stones!" << std::endl << std::endl; 
		stones /= 6.0; 
		
		char stone_color [1024]; 
		sprintf( stone_color, "rgb(%u,%u,0)",(int)(255.0*stones),(int)(215.0*stones) );
		output[2] = stone_color; 
		
		return output; 
	}	
	
	// avengers are red 
	if( pCell->type == avenger.type )
	{
		std::string color = "red"; 

		if( pCell->custom_data["attacking"] > 0.5 )
		{ color = "mistyrose"; }

		output[0] = color;  
		output[2] = color; 
		output[3] = color; 
		
		if( pCell->custom_data["infinity_stones"] > 0.5 )
		{ output[2] = "gold"; } 
	
		return output; 
	}	
	
	
	
	return output; 
}

void thanos_snap( void )
{
	int number_of_cells = (*all_cells).size(); 
	for( int n=0; n < number_of_cells ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{
			if( pC->type != Thanos.type )
			{
				
				if( UniformRandom() <= 0.5 )
				{
					pC->start_death( 0); 
					pC->functions.custom_cell_rule = sad_blowing_away; 
				}
				
			}
			
		}
		
	}
	
	
	return; 
}

void sad_blowing_away( Cell* pCell, Phenotype& phenotype , double dt )
{
	pCell->velocity[0] += 1.0; 
	return; 
}

void cell1_attacks_cell2( Cell* pCell1 , Cell* pCell2 , double dt )
{
	if( pCell1 == pCell2 )
	{ return; } 
	
	#pragma omp critical 
	{
		std::cout << pCell1->type_name << " attacks " << pCell2->type_name << " ... " << std::endl; 
	
		// force of battle by pCell1 and pCell2 
		double force1 = pCell1->custom_data["health"] * pCell1->custom_data["strength"]; 
		double force2 = pCell2->custom_data["health"] * pCell2->custom_data["strength"]; 
		double force = force1+force2; 
		
		// distribution of the force (damage) on pCell1 and pCell2 
		double rel_damage1 = force2 / force; 
		double rel_damage2 = force1 / force; 
		
		// mess up the health 
		// implicit scheme 
		// dH/dt = -r_attack * rel_damage * H ;  
		
		double r_A = pCell1->custom_data["attack_rate"]; 
		double r_heal1 = pCell1->custom_data["recovery_rate"]; 
		double r_heal2 = pCell2->custom_data["recovery_rate"]; 
		
		pCell1->custom_data["health"] /= (1.0 + dt*r_A*rel_damage1); 
		pCell2->custom_data["health"] /= (1.0 + dt*r_A*rel_damage2); 
		
		// if Cell1 is too damaged, kill it off 
		
		if( pCell1->custom_data["health"] < pCell1->custom_data["death_threshold"] )
		{
			std::cout << "\t" << pCell1->type_name << " died in battle " << std::endl; 
			
			pCell1->start_death( 0 ); 
			pCell1->functions.custom_cell_rule = NULL; 
			pCell1->functions.update_phenotype = NULL; 
			
			// its stones go to #2
			double stones = pCell1->custom_data["infinity_stones"]; 
			pCell1->custom_data["infinity_stones"] = 0; 
			pCell2->custom_data["infinity_stones"] += stones; 			
		}
		
		// if Cell2 is too damaged, kill it off 
		
		if( pCell2->custom_data["health"] < pCell2->custom_data["death_threshold"] )
		{
			std::cout << "\t" << pCell2->type_name << " died in battle " << std::endl; 
			
			pCell2->start_death( 0 ); 
			pCell2->functions.custom_cell_rule = NULL; 
			pCell2->functions.update_phenotype = NULL; 
			
			// its stones go to #1
			double stones = pCell2->custom_data["infinity_stones"]; 
			pCell2->custom_data["infinity_stones"] = 0; 
			pCell1->custom_data["infinity_stones"] += stones; 			
		}
		
	}
	return; 
}

void Thanos_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// am I dead? 
	if( phenotype.death.dead == true )
	{
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// try the snap 
	static bool snap_done = false; 
//	if( PhysiCell_globals.current_time >= parameters.doubles( "thanos_snap_time" ) && snap_done == false )
	if( pCell->custom_data["infinity_stones"] > 5.5 && snap_done == false )
	{
		std::cout << std::endl << std::endl << 
		"======================" << std::endl << 
		"SNAP! I am inevitable!" << std::endl << 
		"======================" << std::endl << std::endl;  
		thanos_snap(); 
		snap_done = true; 
		
		// give some time for the cinematic ending 
		
		PhysiCell_settings.max_time = 
			PhysiCell_globals.current_time + 300.0; 
	}				
	
	// look for nearby things to attack 
	Cell* pC; 
	for( int n =0 ; n < pCell->cells_in_my_container().size() ; n++ )
	{
		pC= pCell->cells_in_my_container()[n];
		pCell->custom_data["attacking"] = 0; 
		if( pC!= pCell && pC->phenotype.death.dead == false ) 
		{
			cell1_attacks_cell2( pCell , pC, dt ); 
			pCell->custom_data["attacking"] = 1; 			
		}
		
	}
	
	return; 
};

void avenger_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// try the snap 
	
	
	// look for nearby things to attack
	// only attack Thanos 
	
	Cell* pC; 
	for( int n =0 ; n < pCell->cells_in_my_container().size() ; n++ )
	{
		pCell->custom_data["attacking"] = 0; 
		pC= pCell->cells_in_my_container()[n];
		if( pC!= pCell && pC->phenotype.death.dead == false && pC->type == Thanos.type ) 
		{
			cell1_attacks_cell2( pCell , pC, dt ); 
			pCell->custom_data["attacking"] = 1; 
		}
	}
	
	return; 
};


void avenger_taxis( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int avenger_i = microenvironment.find_density_index( "avenger" );

	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(avenger_i); 
	normalize( &(phenotype.motility.migration_bias_direction) ); 

	return; 
}

void thanos_taxis( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int thanos_i = microenvironment.find_density_index( "thanos" );

	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(thanos_i); 
	normalize( &(phenotype.motility.migration_bias_direction) ); 

	return; 
}

