#include "../include/solver_structure.hpp"
#include <vector>

CBodyForceModelSolver::CBodyForceModelSolver(void) : CSolver() {

}
CBodyForceModelSolver::CBodyForceModelSolver(CGeometry *geometry, CConfig *config, CSolver *fluidsolver, unsigned short iMesh) : CSolver() {
	unsigned short iZone = config->GetiZone();
	unsigned short BFM_Formulation = config->GetBFM_Formulation();
	SetBFM_Formulation(BFM_Formulation);

	cout << "Initializing body-force model" << endl;
	Rotation_vector[0] = config->GetRotation_Rate_X(iZone);
	Rotation_vector[1] = config->GetRotation_Rate_Y(iZone);
	Rotation_vector[2] = config->GetRotation_Rate_Z(iZone);
	cout << "Stored rotation vector components " << endl;
	Rotation_rate = config->GetBody_Force_Rotation()* 2 * M_PI/60.0;
	unsigned long n_points = geometry->GetnPoint();

	proj_vector_axial  = new su2double*[3];
	proj_vector_tangential  = new su2double*[3];
	proj_vector_radial  = new su2double*[3];
	relative_vel = new su2double*[3];
	cyl_coordinates = new su2double*[2];
	cyl_coordinates[0] = new su2double[n_points];
	cyl_coordinates[1] = new su2double[n_points];
	for(size_t ic=0; ic<3; ic++){
		proj_vector_axial[ic] = new su2double[n_points];
		proj_vector_tangential[ic] = new su2double[n_points];
		proj_vector_radial[ic] = new su2double[n_points];
		relative_vel[ic] = new su2double[n_points];
	}
	cout << "End of BFM constructor" << endl;
	
}

void CBodyForceModelSolver::LoadRestart(CGeometry **geometry, CSolver ***fluidsolver, CConfig *config, int val_iter, bool val_update_geo){
	cout <<"Loading restart from body-force model" << endl;
}

void CBodyForceModelSolver::PreprocessBFMParams(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	cout << "Preprocessing function called" << endl;
    ReadInputFile(config, geometry);
    InterpolateGeometry(geometry, config, fluidSolver);
	if(kind_bfm==1){
    	ComputeBlockageGradient(geometry, config, fluidSolver);
	}
	ComputeCylProjections(geometry, config, fluidSolver);
}
void CBodyForceModelSolver::ComputeRelVelocity(CGeometry *geometry, CConfig *config, CSolver *fluidsolver){
	su2double *Coord, *U_i,*Geometric_Parameters;
	su2double W_ax, W_r, W_th, rotFac;
	for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){
		Coord = geometry->node[iPoint]->GetCoord();
		U_i = fluidsolver->node[iPoint]->GetSolution();
		// Computing relative velocity components in axial, radial and tangential direction
		W_ax = 0;
		W_r = 0;
		W_th = 0;

		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		rotFac = Geometric_Parameters[6];
		// Adding to the relative velocity components through dot product of absolute velocity and respective unit vector
		for(int iDim = 0; iDim < nDim; iDim++){
			W_ax += proj_vector_axial[iDim][iPoint] * U_i[iDim+1]/U_i[0];
			W_th += proj_vector_tangential[iDim][iPoint] * U_i[iDim+1]/U_i[0];
			W_r += proj_vector_radial[iDim][iPoint] * U_i[iDim+1]/U_i[0];
		}
		// Induced velocity due to rotation is subtracted from tangential velocity
		W_th -= rotFac * Rotation_rate * cyl_coordinates[1][iPoint];
		relative_vel[0][iPoint] = W_ax;
		relative_vel[1][iPoint] = W_th;
		relative_vel[2][iPoint] = W_r;
	}
}

void CBodyForceModelSolver::ComputeCylProjections(CGeometry *geometry, CConfig *config, CSolver *fluidsolver){
	cout << "Computing cylindrical unit vectors and coordinates" << endl;
	su2double *Coord;
	su2double rot_dot_x, rot_dot_rot;
	su2double ax, radius;
	su2double axial_vector[nDim], radial_vector[nDim];

	for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){
		Coord = geometry->node[iPoint]->GetCoord();

		rot_dot_x = 0;	// dot product between axial unit vector and coordinate vector
		rot_dot_rot = 0;	// dot product of axial unit vector
		radius = 0;	// Radial coordinate
		ax = 0;		// Axial coordinate
		if (nDim == 3){
			// Computing axial dot products
			for(int iDim = 0; iDim < nDim; iDim ++){
				rot_dot_x += Rotation_vector[iDim]*Coord[iDim];
				rot_dot_rot += Rotation_vector[iDim]* Rotation_vector[iDim];
			}
			//
			for(int iDim = 0; iDim < nDim; iDim ++){
				axial_vector[iDim] = (rot_dot_x/rot_dot_rot)* Rotation_vector[iDim]; // Parallel projection of coordinate vector onto rotation axis
				radial_vector[iDim] = Coord[iDim] - axial_vector[iDim];	// Orthogonal projection of coordinate vector onto rotation axis
				ax += axial_vector[iDim]*axial_vector[iDim];	
				radius += radial_vector[iDim]*radial_vector[iDim];
				proj_vector_axial[iDim][iPoint] = Rotation_vector[iDim];	// Appending axial unit vector component
			}
			radius = sqrt(radius); // Computing radial coordinate
			ax = sqrt(ax);	// Compuring axial coordinate
			if(rot_dot_x < 0.0){ax = -ax;}	// In case the coordinate parallel projection is negative, the axial coordinate is flipped
			cyl_coordinates[0][iPoint] = ax;
			cyl_coordinates[1][iPoint] = radius;
			// Computing tangential unit vector components through cross product of radial and axial vectors
			proj_vector_tangential[0][iPoint] = (proj_vector_axial[1][iPoint]*radial_vector[2] - proj_vector_axial[2][iPoint]*radial_vector[1])/radius;
			proj_vector_tangential[1][iPoint] = (-proj_vector_axial[0][iPoint]*radial_vector[2] + proj_vector_axial[2][iPoint]*radial_vector[0])/radius;
			proj_vector_tangential[2][iPoint] = (proj_vector_axial[0][iPoint]*radial_vector[1] - proj_vector_axial[1][iPoint]*radial_vector[0])/radius;

			// Computing radial unit vector through normilization of radial vector
			for(int iDim = 0; iDim < nDim; iDim ++){proj_vector_radial[iDim][iPoint] = radial_vector[iDim]/radius;}
			
			// Currently, the method is incompatible with 2D simulations
			}else{
				radius = config->GetBody_Force_Radius();
				ax = Coord[0];
				cyl_coordinates[0][iPoint] = ax;
				cyl_coordinates[1][iPoint] = radius;
			}
	}
}
void CBodyForceModelSolver::BFM_Hall(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	unsigned long iPoint;
	su2double *U_i, *V_i, *Coord_i, *Geometric_Parameters;
	su2double pi = M_PI, pitch;
	su2double bfFac, b, Nx, Nt, Nr, x_le, rotFac, BF_blades, radius;
	su2double W_ax, W_th, W_r;
	su2double WdotN, W_nx, W_nr, W_nth, W_px, W_pth, W_pr, W_p, W;
	su2double delta, F_n;
	su2double F_ax, F_th, F_r, e_source;
	su2double BF_source[nDim+2], F[nDim];


	for(iPoint=0; iPoint<geometry->GetnPoint(); iPoint++){
		U_i = fluidsolver->node[iPoint]->GetSolution();
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		bfFac = Geometric_Parameters[0]; 	// Body-force factor. Multiplies all body-forces.
		b = Geometric_Parameters[1];			// Metal blockage factor
		Nx = Geometric_Parameters[2]; 		// Camber normal component in axial direction
		Nt = Geometric_Parameters[3]; 			// Camber normal component in tangential direction
		Nr = Geometric_Parameters[4];			// Camber normal component in radial direction
		x_le = Geometric_Parameters[5]; 		// Axial distance from leading edge
		rotFac = Geometric_Parameters[6];	// Rotation factor. Multiplies the body-force rotation value.
		BF_blades = Geometric_Parameters[7];				// Blade row blade count

		W_ax = relative_vel[0][iPoint];
		W_th = relative_vel[1][iPoint];
		W_r = relative_vel[2][iPoint];

		radius = cyl_coordinates[1][iPoint];
		pitch = 2 * pi * radius / BF_blades;

		WdotN = W_ax * Nx + W_r * Nr + W_th * Nt;		// Dot product of relative velocity and camber normal vector
		W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
		W_px = W_ax - W_nx, W_pr = W_r - W_nr, W_pth = W_th - W_nth;  // Relative velocity components parallel to the blade 
		
		W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
		W = sqrt(W_ax * W_ax + W_r * W_r + W_th * W_th);					// Relative velocity magnitude
		// Calculating the deviation angle
		delta = asin(WdotN / W);

		F_n = -bfFac * pi * delta * (1 / pitch) * (1 / abs(Nt)) * W * W;

		// Transforming the normal and parallel force components to cyllindrical coordinates
			
		F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / W_p));
		F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / W_p));
		F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / W_p));
		e_source = rotFac * Rotation_rate * radius * F_th;
		
		// Appending Cartesial body-forces to body-force vector
		// for(int iDim=0; iDim<nDim; iDim++){
		// 	F[iDim] = U_i[0] * (F_ax*i_ax[iDim] + F_th*i_theta[iDim] + F_r*i_r[iDim]);
		// }
					// Appending Cartesial body-forces to body-force vector
		for(int iDim=0; iDim<nDim; iDim++){
			F[iDim] = U_i[0] * (F_ax*proj_vector_axial[iDim][iPoint] + F_th*proj_vector_tangential[iDim][iPoint] + F_r*proj_vector_radial[iDim][iPoint]);
		}
		// Appending source term values to the body-force source term vector
		BF_source[0] = 0.0;	// empty density component
		// Appending momentum source terms
		for(int iDim = 0; iDim < nDim; iDim ++) {
			BF_source[iDim + 1] = F[iDim];
		}
			// Appending energy source term
		BF_source[nDim+1] = U_i[0] * e_source;
		fluidsolver->node[iPoint]->SetBody_Force_Source(BF_source);
	}
}
void CBodyForceModelSolver::BFM_Thollet(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	/*
	This function computes the body-force field for the Euler solver. 
	*/
    /*--- Load all relevant values from config file ---*/
    unsigned long iPoint;
    su2double *U_i, *V_i, *Coord_i, *Geometric_Parameters;
    su2double gamma = config->GetGamma(), R_gas = config->GetGas_Constant(), BF_radius = config->GetBody_Force_Radius();

    /*--- Initialize common variables ---*/
    su2double pi = M_PI, pitch;
	su2double radius, ax;
	su2double W_ax, W_r, W_th, W;
	su2double bfFac, b, Nx, Nt, Nr, x_le, rotFac, BF_blades;
	su2double WdotN, W_nx, W_nth, W_nr, W_px, W_pth, W_pr, W_p;
	su2double delta, V_sound, M_rel;
	su2double F_n_inc, F_n, F_p;
	su2double F_ax, F_r, F_th, F_x, F_y, F_z, e_source, F[nDim];
	su2double C_f, Re_x, mu=1.716E-5;
	su2double Kprime, K=1;
	su2double BF_source[nDim+2];

    for ( iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
		
		// Extracting flow variables and primitive variables.
        U_i = fluidsolver->node[iPoint]->GetSolution();
        V_i = fluidsolver->node[iPoint]->GetPrimitive();
		
		// Extracting node coordinates 
        Coord_i = geometry->node[iPoint]->GetCoord();
		// // In order to determine the relative velocity components, the local unit vectors in axial
		// // ,radial and tangential direction have to be determined.
		
		// rot_dot_x = 0;	// dot product between axial unit vector and coordinate vector
		// rot_dot_rot = 0;	// dot product of axial unit vector
		// radius = 0;	// Radial coordinate
		// ax = 0;		// Axial coordinate
		// if (nDim == 3){
		// 	// Computing axial dot products
		// 	for(int iDim = 0; iDim < nDim; iDim ++){
		// 		rot_dot_x += Rotation_vector[iDim]*Coord_i[iDim];
		// 		rot_dot_rot += Rotation_vector[iDim]* Rotation_vector[iDim];
		// 	}
		// 	//
		// 	for(int iDim = 0; iDim < nDim; iDim ++){
		// 		axial_vector[iDim] = (rot_dot_x/rot_dot_rot)* Rotation_vector[iDim]; // Parallel projection of coordinate vector onto rotation axis
		// 		radial_vector[iDim] = Coord_i[iDim] - axial_vector[iDim];	// Orthogonal projection of coordinate vector onto rotation axis
		// 		ax += axial_vector[iDim]*axial_vector[iDim];	
		// 		radius += radial_vector[iDim]*radial_vector[iDim];
		// 		i_ax[iDim] = Rotation_vector[iDim];	// Appending axial unit vector component
		// 	}
		// 	radius = sqrt(radius); // Computing radial coordinate
		// 	ax = sqrt(ax);	// Compuring axial coordinate
		// 	if(rot_dot_x < 0.0){ax = -ax;}	// In case the coordinate parallel projection is negative, the axial coordinate is flipped

		// 	// Computing tangential unit vector components through cross product of radial and axial vectors
		// 	i_theta[0] = (i_ax[1]*radial_vector[2] - i_ax[2]*radial_vector[1])/radius;
		// 	i_theta[1] = (-i_ax[0]*radial_vector[2] + i_ax[2]*radial_vector[0])/radius;
		// 	i_theta[2] = (i_ax[0]*radial_vector[1] - i_ax[1]*radial_vector[0])/radius;

		// 	// Computing radial unit vector through normilization of radial vector
		// 	for(int iDim = 0; iDim < nDim; iDim ++){i_r[iDim] = radial_vector[iDim]/radius;}
			
		// 	// Currently, the method is incompatible with 2D simulations
		// 	}else{
		// 	radius = BF_radius;
		// 	ax = Coord_i[0];
		// 	}	
			
			
			Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
			bfFac = Geometric_Parameters[0]; 	// Body-force factor. Multiplies all body-forces.
			b = Geometric_Parameters[1];			// Metal blockage factor
			Nx = Geometric_Parameters[2]; 		// Camber normal component in axial direction
			Nt = Geometric_Parameters[3]; 			// Camber normal component in tangential direction
			Nr = Geometric_Parameters[4];			// Camber normal component in radial direction
			x_le = Geometric_Parameters[5]; 		// Axial distance from leading edge
			rotFac = Geometric_Parameters[6];	// Rotation factor. Multiplies the body-force rotation value.
			BF_blades = Geometric_Parameters[7];				// Blade row blade count
			
			// // Computing relative velocity components in axial, radial and tangential direction
			// W_ax = 0;
			// W_r = 0;
			// W_th = 0;


			W_ax = relative_vel[0][iPoint];
			W_th = relative_vel[1][iPoint];
			W_r = relative_vel[2][iPoint];

			ax = cyl_coordinates[0][iPoint];
			radius = cyl_coordinates[1][iPoint];

			// Adding to the relative velocity components through dot product of absolute velocity and respective unit vector
			// for(int iDim = 0; iDim < nDim; iDim++){
			// 	W_ax += proj_vector_axial[iDim][iPoint] * U_i[iDim+1]/U_i[0];
			// 	W_th += proj_vector_tangential[iDim][iPoint] * U_i[iDim+1]/U_i[0];
			// 	W_r += proj_vector_radial[iDim][iPoint] * U_i[iDim+1]/U_i[0];
			// }
			// Induced velocity due to rotation is subtracted from tangential velocity
			//W_th -= rotFac * omegaR * radius;

			WdotN = W_ax * Nx + W_r * Nr + W_th * Nt;		// Dot product of relative velocity and camber normal vector
			W_nx = WdotN * Nx, W_nr = WdotN * Nr, W_nth = WdotN * Nt;		// Relative velocity components normal to the blade
			W_px = W_ax - W_nx, W_pr = W_r - W_nr, W_pth = W_th - W_nth;  // Relative velocity components parallel to the blade 
			
			W_p = sqrt(W_px * W_px + W_pr * W_pr + W_pth * W_pth);	// Parallel relative velocity magnitude
			W = sqrt(W_ax * W_ax + W_r * W_r + W_th * W_th);					// Relative velocity magnitude
			
			// Calculating the deviation angle
			delta = asin(WdotN / W);
			// Calculating te local, relative Mach number
			V_sound = sqrt(gamma * R_gas * V_i[0]);
			M_rel = W / V_sound;
			
			// Blade pitch 
			pitch = 2 * pi * radius / BF_blades;
		
			// Starting the process of calculating the actual body-force magnitudes
			
			
			// Incompressible normal force magnitude 
			F_n_inc = pi * delta * (1 / pitch) * (1 / abs(Nt)) * (1 / b) * W * W;
			
			// Calculating the friction factor for parallel force calculation
			
			
			// Axial Reynolds number, based on leading edge distance.
			Re_x = (abs(ax - x_le) * W * U_i[0]) / mu;
			if(Re_x == 0.0){
				Re_x = (0.001 * W * U_i[0]) / mu;
			}
			C_f = 0.0592 * pow(Re_x, -0.2) ;
		
			// Calculating the normal force compressibility factor
			if (M_rel < 1) {
				Kprime = 1 / (sqrt(1 - (M_rel * M_rel)));
				if (Kprime <= 3) {
					K = Kprime;
				}
				if (Kprime > 3) {
					K = 3;
				}
			}
			if (M_rel > 1) {
				Kprime = 2 / (pi * sqrt((M_rel * M_rel) - 1));
				if (Kprime <= 3) {
					K = Kprime;
				}
				if (Kprime > 3) {
					K = 3;
				}
			}
		
			
			// Calculating the final values for the normal and parallel force
			F_n = -bfFac * K * F_n_inc;
			F_p = -bfFac * C_f * W * W * (1 / pitch) * (1 / abs(Nt)) * (1 / b);
			
			
			// Transforming the normal and parallel force components to cyllindrical coordinates
			
			F_ax = F_n * (cos(delta) * Nx - sin(delta) * (W_px / W_p)) + F_p * W_ax / W;
			F_r = F_n * (cos(delta) * Nr - sin(delta) * (W_pr / W_p)) + F_p * W_r / W;
			F_th = F_n * (cos(delta) * Nt - sin(delta) * (W_pth / W_p)) + F_p * W_th / W;
			e_source = rotFac * Rotation_rate * radius * F_th;
			
			// Appending Cartesial body-forces to body-force vector
			// for(int iDim=0; iDim<nDim; iDim++){
			// 	F[iDim] = U_i[0] * (F_ax*i_ax[iDim] + F_th*i_theta[iDim] + F_r*i_r[iDim]);
			// }
						// Appending Cartesial body-forces to body-force vector
			for(int iDim=0; iDim<nDim; iDim++){
				F[iDim] = U_i[0] * (F_ax*proj_vector_axial[iDim][iPoint] + F_th*proj_vector_tangential[iDim][iPoint] + F_r*proj_vector_radial[iDim][iPoint]);
			}
			// Appending source term values to the body-force source term vector
			BF_source[0] = 0.0;	// empty density component
			// Appending momentum source terms
			for(int iDim = 0; iDim < nDim; iDim ++) {
				BF_source[iDim + 1] = F[iDim];
			}
				// Appending energy source term
			BF_source[nDim+1] = U_i[0] * e_source;
			fluidsolver->node[iPoint]->SetBody_Force_Source(BF_source);
			//fluidsolver->node[iPoint]->SetBodyForceVector_Turbo(BF_source);
		};
	
};
void CBodyForceModelSolver::AllocateMemory(){
	rows_axial.resize(n_rows);
	rows_axial_chord.resize(n_rows);
	rows_blade_count.resize(n_rows);
	rows_blockage.resize(n_rows);
	rows_Nr.resize(n_rows);
	rows_Nt.resize(n_rows);
	rows_Nx.resize(n_rows);
	rows_radial.resize(n_rows);
	rows_rotation.resize(n_rows);
	rows_X_le.resize(n_rows);
	for(int i_row=0; i_row<n_rows; i_row++){
		rows_axial.at(i_row).resize(n_sec);
		rows_axial_chord.at(i_row).resize(n_sec);
		rows_blockage.at(i_row).resize(n_sec);
		rows_Nr.at(i_row).resize(n_sec);
		rows_Nt.at(i_row).resize(n_sec);
		rows_Nx.at(i_row).resize(n_sec);
		rows_radial.at(i_row).resize(n_sec);
		rows_X_le.at(i_row).resize(n_sec);
		for(int i_sec=0; i_sec<n_sec; i_sec++){
			rows_axial.at(i_row).at(i_sec).resize(n_axial);
			rows_axial_chord.at(i_row).at(i_sec).resize(n_axial);
			rows_blockage.at(i_row).at(i_sec).resize(n_axial);
			rows_Nr.at(i_row).at(i_sec).resize(n_axial);
			rows_Nt.at(i_row).at(i_sec).resize(n_axial);
			rows_Nx.at(i_row).at(i_sec).resize(n_axial);
			rows_radial.at(i_row).at(i_sec).resize(n_axial);
			rows_X_le.at(i_row).at(i_sec).resize(n_axial);
		}
	}
}
void CBodyForceModelSolver::ReadInputFile(CConfig *config, CGeometry *geometry){
    cout << "Reading Body Force Model geometry input file" <<endl;

    int start_line = 0;
	int end_line = 0;
    nDim = geometry->GetnDim();
    
    
	string line;
	// Opening body-force model input file
	std::ifstream inputFile((config->GetBFM_inputName()).c_str());
	
	// Reading number of blade rows, number of radial blade sections and axial points from first line
	inputFile >> n_rows >> n_sec >> n_axial;
    //TODO:check for proper values of n_rows, n_blade and n_points
	AllocateMemory();

	su2double axial_coordinate, radial_coordinate, N_axial, N_tangential, N_radial,
	blockage_factor, axial_c_LE, axial_chord, rotation_factor;
	int blade_count;
    for (int q=0; q < n_rows; q++){
		for (int p=0; p < n_sec; p++){
			getline(inputFile, line);
			start_line = p * n_axial + 1;
			end_line = (p + 1) * n_axial;
			for (int i=0; i<n_axial; i++){
				getline(inputFile, line);
				// Storing values read on the line
				inputFile >> axial_coordinate >> radial_coordinate >> N_axial >> N_tangential >> N_radial >> blockage_factor >> axial_c_LE >> axial_chord >> rotation_factor >> blade_count;
				rows_axial_chord.at(q).at(p).at(i) = axial_chord;
				rows_axial.at(q).at(p).at(i) = axial_coordinate;
				rows_radial.at(q).at(p).at(i) = radial_coordinate;
				rows_Nx.at(q).at(p).at(i) = N_axial;
				rows_Nt.at(q).at(p).at(i) = N_tangential;
				rows_Nr.at(q).at(p).at(i) = N_radial;
				rows_blockage.at(q).at(p).at(i) = blockage_factor;
				rows_X_le.at(q).at(p).at(i) = axial_c_LE;

				// section_axial.push_back(axial_coordinate);
				// section_radial.push_back(radial_coordinate);
				// section_Nx.push_back(N_axial);
				// section_Nt.push_back(N_tangential);
				// section_Nr.push_back(N_radial);
				// section_blockage.push_back(blockage_factor);
				// section_X_le.push_back(axial_c_LE);
				// section_axial_chord.push_back(axial_chord);
				// section_rotation.push_back(rotation_factor);
				if(axial_coordinate <= x_min){
					x_min = axial_coordinate;
				}
					
				
				
			}
			// blade_axial.push_back(section_axial);
			// blade_radial.push_back(section_radial);
			// blade_Nx.push_back(section_Nx);
			// blade_Nt.push_back(section_Nt);
			// blade_Nr.push_back(section_Nr);
			// blade_blockage.push_back(section_blockage);
			// blade_X_le.push_back(section_X_le);
			// blade_axial_chord.push_back(section_axial_chord);
		}
		// rows_axial.push_back(blade_axial);
		// rows_radial.push_back(blade_radial);
		// rows_Nx.push_back(blade_Nx);
		// rows_Nt.push_back(blade_Nt);
		// rows_Nr.push_back(blade_Nr);
		// rows_blockage.push_back(blade_blockage);
		// rows_X_le.push_back(blade_X_le);
		// rows_axial_chord.push_back(blade_axial_chord);
		rows_rotation.at(q) = rotation_factor;
		rows_blade_count.at(q) = blade_count;

	}
}
void CBodyForceModelSolver::InterpolateGeometry(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	cout << "Interpolating geometric parameters to mesh" << endl;
  unsigned long iPoint;
  unsigned short iMesh, iDim;
  unsigned short iZone = config->GetiZone();
	su2double BF_radius = config->GetBody_Force_Radius();	
	// Subtracting 1 meter from minimum axial coordinate to ensure it's located outside of the domain
	x_min = x_min - 1.0;
	
	su2double b=1.0, Nx = 0.0, Nt = 1.0, Nr = 0.0, bfFac = 0.0, x_le = 0.0, rotFac = 0.0, bladeCount = 1, chord = 1.0;
	su2double BodyForceParams[9] = {bfFac, b, Nx, Nt, Nr, x_le, rotFac, bladeCount, chord};
	
	su2double x1, x2, r1, r2;
	int nInt=0;
	bool inside = false;
	su2double dist = 0, deNom=0;
	su2double eNum_b = 0, eNum_Nx = 0, eNum_Nt = 0, eNum_Nr = 0, eNum_x_le = 0, eNum_chord = 0;
	su2double cell_axial_coordinate, cell_radial_coordinate;
	su2double rot_dot_x, rot_dot_rot;
	su2double axial[nDim], radial[nDim];
	su2double x_side [4][2] {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
	//vector<vector<su2double>> x_side;
	vector<su2double> blockage_side, Nx_side, Nt_side, Nr_side, x_le_side, chord_side;
	vector<su2double> edge {0, 0}, v_x {0, 0}, v_opp {0, 0}, proj_x {0, 0}, proj_opp {0, 0};
	su2double orientation;
	// Looping over points in the domain to interpolate body-force parameters onto the respective cells
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		su2double *Coord = geometry->node[iPoint]->GetCoord(); // Coordinate vector
		rot_dot_x = 0;	// dot product between rotation axis and coordinate vector
		rot_dot_rot = 0;	// dot product of rotation axis
		cell_axial_coordinate = 0;
		cell_radial_coordinate = 0;
		if (nDim == 3){
			// Computing dot products
			for(int iDim = 0; iDim < nDim; iDim ++){
				rot_dot_x += Rotation_vector[iDim]*Coord[iDim];
				rot_dot_rot += Rotation_vector[iDim]* Rotation_vector[iDim];
			}
			for(int iDim = 0; iDim < nDim; iDim ++){
				axial[iDim] = (rot_dot_x/rot_dot_rot)* Rotation_vector[iDim]; // computing axial coordinate through parallel projection of coordinate vector onto rotation axis
				radial[iDim] = Coord[iDim] - axial[iDim];	// computing radial vector through normal projection of coordinate vector onto rotation axis
				cell_axial_coordinate += axial[iDim]*axial[iDim];
				cell_radial_coordinate += radial[iDim]*radial[iDim];
			}
			
			cell_axial_coordinate = sqrt(cell_axial_coordinate);	// Computing axial coordinate
			if(rot_dot_x < 0.0){cell_axial_coordinate = -cell_axial_coordinate;}	// In case parallel projection is negative, axial coordinate is flipped
			cell_radial_coordinate = sqrt(cell_radial_coordinate);	// Computing radial coordinatecout << iPoint << " " << Coord[0] << " " << Coord[1] << " " << Coord[2] << " " << cell_axial_coordinate << " " << cell_radial_coordinate << " " << sqrt(Coord[0]*Coord[0] + Coord[1]*Coord[1])<< endl;
		}else{
			cell_radial_coordinate = BF_radius;
			cell_axial_coordinate = Coord[0];
		}
		//su2double b=1.0, Nx = 0.0, Nt = 1.0, Nr = 0.0, bfFac = 0.0, x_le = 0.0, rotFac = 0.0, bladeCount = 1, chord = 1.0;
		//su2double BodyForceParams[9] = {bfFac, b, Nx, Nt, Nr, x_le, rotFac, bladeCount, chord};
		// Initializing body-force parameters with outside values
		b = 1.0;
		Nx = 0.0;
		Nt = 1.0;
		Nr = 0.0;
		bfFac = 0.0;
		x_le = 0.0;
		rotFac = 0.0;
		bladeCount = 1;
		chord = 1.0;
		//su2double BodyForceParams[9] = {bfFac, b, Nx, Nt, Nr, x_le, rotFac, bladeCount, chord};
		BodyForceParams[0] = bfFac;
		BodyForceParams[1] = b;
		BodyForceParams[2] = Nx;
		BodyForceParams[3] = Nt;
		BodyForceParams[4] = Nr;
		BodyForceParams[5] = x_le;
		BodyForceParams[6] = rotFac;
		BodyForceParams[7] = bladeCount;
		BodyForceParams[8] = chord;

		
		
		for (int q=0; q < n_rows; q++){
			for(int j=0; j < n_sec-1; j++){
				for(int n = 0; n < n_axial - 1; n++){
					// Appending interpolation values to sides loop
					x_side[0][0] = rows_axial.at(q).at(j).at(n);
					x_side[0][1] = rows_radial.at(q).at(j).at(n);
					x_side[1][0] = rows_axial.at(q).at(j+1).at(n);
					x_side[1][1] = rows_radial.at(q).at(j+1).at(n);
					x_side[2][0] = rows_axial.at(q).at(j+1).at(n+1);
					x_side[2][1] = rows_radial.at(q).at(j+1).at(n+1);
					x_side[3][0] = rows_axial.at(q).at(j).at(n+1);
					x_side[3][1] = rows_radial.at(q).at(j).at(n+1);
					// x_side = {{rows_axial.at(q).at(j).at(n), rows_radial.at(q).at(j).at(n)},
					// 			{rows_axial.at(q).at(j+1).at(n), rows_radial.at(q).at(j+1).at(n)},
					// 			{rows_axial.at(q).at(j+1).at(n+1),rows_radial.at(q).at(j+1).at(n+1)},
					// 			{rows_axial.at(q).at(j).at(n+1),rows_radial.at(q).at(j).at(n+1)}};
					blockage_side = {rows_blockage.at(q).at(j).at(n), rows_blockage.at(q).at(j+1).at(n), rows_blockage.at(q).at(j+1).at(n+1), rows_blockage.at(q).at(j).at(n+1)};
					Nx_side = {rows_Nx.at(q).at(j).at(n), rows_Nx.at(q).at(j+1).at(n), rows_Nx.at(q).at(j+1).at(n+1), rows_Nx.at(q).at(j).at(n+1)};
					Nt_side = {rows_Nt.at(q).at(j).at(n), rows_Nt.at(q).at(j+1).at(n), rows_Nt.at(q).at(j+1).at(n+1), rows_Nt.at(q).at(j).at(n+1)};
					Nr_side = {rows_Nr.at(q).at(j).at(n), rows_Nr.at(q).at(j+1).at(n), rows_Nr.at(q).at(j+1).at(n+1), rows_Nr.at(q).at(j).at(n+1)};
					x_le_side = {rows_X_le.at(q).at(j).at(n), rows_X_le.at(q).at(j+1).at(n), rows_X_le.at(q).at(j+1).at(n+1), rows_X_le.at(q).at(j).at(n+1)};
					chord_side = {rows_axial_chord.at(q).at(j).at(n), rows_axial_chord.at(q).at(j+1).at(n), rows_axial_chord.at(q).at(j+1).at(n+1), rows_axial_chord.at(q).at(j).at(n+1)};
					nInt = 0; // number of intersections
					inside = true;	// point is outside of the blade region by default
					// Looping over sides to count number of intersections
					
					for(int p=0; p < 4; p++){
						int p_next = (p+1) % 4;
						int p_opp = (p+2) % 4;
						
						edge.at(0) = x_side[p_next][0] - x_side[p][0];
						edge.at(1) = x_side[p_next][1] - x_side[p][1];
						v_opp.at(0) = x_side[p_opp][0] - x_side[p][0];
						v_opp.at(1) = x_side[p_opp][1] - x_side[p][1];
						v_x.at(0) = cell_axial_coordinate - x_side[p][0];
						v_x.at(1) = cell_radial_coordinate - x_side[p][1];
						proj_opp.at(0) = v_opp.at(0) - (vector_dot_product(edge, v_opp)/vector_dot_product(edge, edge))*edge.at(0);
						proj_opp.at(1) = v_opp.at(1) - (vector_dot_product(edge, v_opp)/vector_dot_product(edge, edge))*edge.at(1);
						proj_x.at(0) = v_x.at(0) - (vector_dot_product(edge, v_x)/vector_dot_product(edge, edge))*edge.at(0);
						proj_x.at(1) = v_x.at(1) - (vector_dot_product(edge, v_x)/vector_dot_product(edge, edge))*edge.at(1);
						orientation = vector_dot_product(proj_opp, proj_x);
						if(orientation < 0){
							inside = false;
						}
					}
					// If the point is inside, the body-force parameters are interpolated through distance-weighted average from the corner points
					if (inside){
						dist = 0; deNom=0; eNum_b = 0; eNum_Nx = 0; eNum_Nt = 0; eNum_Nr = 0; eNum_x_le = 0; eNum_chord = 0;
						for(int p = 0; p < 4; p++){
							dist = sqrt((cell_axial_coordinate - x_side[p][0]) * (cell_axial_coordinate - x_side[p][0]) + (cell_radial_coordinate - x_side[p][1]) * (cell_radial_coordinate - x_side[p][1]));
							deNom += 1 / dist;
							eNum_b += blockage_side.at(p) / dist;
							eNum_Nx += Nx_side.at(p) / dist;
							eNum_Nt += Nt_side.at(p) / dist;
							eNum_Nr += Nr_side.at(p) / dist;
							eNum_x_le += x_le_side.at(p) / dist;
							eNum_chord += chord_side.at(p) / dist;
						}
						bfFac = 1.0;
						b = eNum_b / deNom;
						Nx = eNum_Nx / deNom;
						Nt = eNum_Nt / deNom;
						Nr = eNum_Nr / deNom;
						x_le = eNum_x_le / deNom;
						rotFac = rows_rotation.at(q);
						bladeCount = rows_blade_count.at(q);
						chord = eNum_chord / deNom;
					}
					// Storing interpolated geometric parameters into vector
					BodyForceParams[0] = bfFac;
					BodyForceParams[1] = b;
					BodyForceParams[2] = Nx;
					BodyForceParams[3] = Nt;
					BodyForceParams[4] = Nr;
					BodyForceParams[5] = x_le;
					BodyForceParams[6] = rotFac;
					BodyForceParams[7] = bladeCount;
					BodyForceParams[8] = chord;
				}
			}
		}
		// Storing body-force parameter vector onto node
		fluidSolver->node[iPoint]->SetBodyForceParameters(BodyForceParams);
		//cout << x << ", "  << r << ", " << BodyForceParams[0] << ", " << BodyForceParams[1] << ", " << BodyForceParams[2] << ", " << BodyForceParams[3] << ", " << BodyForceParams[4] << endl;
	}
}

void CBodyForceModelSolver::ComputeBlockageGradient(CGeometry *geometry, CConfig *config, CSolver *fluidSolver){
	cout << "Computing metal blockage factor gradient field" << endl;
  unsigned short iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double Blockage_i, Blockage_j, *BFMVec, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2, cvector[nDim], smatrix[nDim][nDim];
  bool singular;
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    //cout << "Coordinate: " << Coord_i[0] << " " << Coord_i[1] << endl;
    /*--- Get primitives from CVariable ---*/
    BFMVec = fluidSolver->node[iPoint]->GetBodyForceParameters();
    Blockage_i = BFMVec[1];
    //cout << "Blockage: " << Blockage_i << endl;
    /*--- Inizialization of variables ---*/
    
    
    for (iDim = 0; iDim < nDim; iDim++){
        cvector[iDim] = 0.0;
    }
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
	/*
    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);
    */
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
	  BFMVec = fluidSolver->node[jPoint]->GetBodyForceParameters();
      Blockage_j = BFMVec[1];
      
	  /*
      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
	  */
	  
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      //cout << "Weight: " << weight << endl;
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        
         for (iDim = 0; iDim < nDim; iDim++){
            cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Blockage_j-Blockage_i)/weight;
		 }
        
      }
      
    }
    //cout << " Cvector[0]: " << cvector[0] << " Cvector[1]: " << cvector[1] << endl;
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    //cout << "r11: " << r11 << " r12: " << r12 << " r22: " << r22 << endl;
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    //cout << "Determinant: " << detR2 << endl;
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        smatrix[0][1] = -r11*r12/detR2;
        smatrix[1][0] = smatrix[0][1];
        smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        smatrix[0][2] = (z13*z33)/detR2;
        smatrix[1][0] = smatrix[0][1];
        smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        smatrix[1][2] = (z23*z33)/detR2;
        smatrix[2][0] = smatrix[0][2];
        smatrix[2][1] = smatrix[1][2];
        smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    //cout << " Smatrix[0][1]: " << smatrix[0][1] <<" Smatrix[1][0]: " << smatrix[1][0] << " Smatrix[1][1]: " << smatrix[1][1] << " Smatrix[0][0]: " << smatrix[0][0] << endl;
    /*--- Computation of the gradient: S*c ---*/
    
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += smatrix[iDim][jDim]*cvector[jDim];
        }
        fluidSolver->node[iPoint]->SetGradient_Blockage(iDim, product);
      }
    
    /*
    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
	*/
  }
}

void CBodyForceModelSolver::ComputeBFMSources(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	ComputeBodyForce_Source(config, geometry, fluidsolver);
	if(kind_bfm==1){
		ComputeBlockage_Source(config, geometry, fluidsolver);
	}
}
void CBodyForceModelSolver::ComputeBodyForce_Source(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	bool thollet=true;
	ComputeRelVelocity(geometry, config, fluidsolver);
	if(kind_bfm==1){
		BFM_Thollet(config, geometry, fluidsolver);
	}
	if(kind_bfm==0){
		BFM_Hall(config, geometry, fluidsolver);
	}
	
}

void CBodyForceModelSolver::ComputeBlockage_Source(CConfig *config, CGeometry *geometry, CSolver *fluidsolver){
	/*
	This function computes the local residual vector due to metal blockage for the body-force Euler solver.
	It uses the blockage factor and its derivatives in x- and r-direction provided by the interpolation function which is 
	performed during the solver initialization.
	*/
	
	// Obtaining the dimension count, flow and primitive variables and the interpolated body-force parameters.
	unsigned long iPoint;
	unsigned long nDim = geometry->GetnDim();
	su2double *U_i, *V_i, enthalpy, *Coord_i, Density, *Geometric_Parameters;
	int iDim;
	
	//  Defining the blockage residual vector.
	su2double Blockage_Vector[nDim + 2];
    //cout << "Blockage function is being called" <<endl;
	// The function loops over all points in the zone, calculating and storing the residual blockage vector for each
	// respective node.
	su2double b = 1.0;
	su2double BGradient[nDim], Blockage_Div = 0.0;
	for ( iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
		// Getting the solution and primitive variables of the current node.
		U_i = fluidsolver->node[iPoint]->GetSolution();
        V_i = fluidsolver->node[iPoint]->GetPrimitive();
		enthalpy = fluidsolver->node[iPoint]->GetEnthalpy();
		// Getting the node coordinates.
		Coord_i = geometry->node[iPoint]-> GetCoord();
		
		// The blockage factor is extracted at the node, resulting from the interpolation during initialization.
		Geometric_Parameters = fluidsolver->node[iPoint]->GetBodyForceParameters();
		b = Geometric_Parameters[1];
		//b = 1.0;
		
		// Obtaining blockage derivative from gradient field, calculated during initialization.
		for(iDim=0; iDim < nDim; iDim++){
			BGradient[iDim] = fluidsolver->node[iPoint]->GetGradient_Blockage(iDim);
			//BGradient[iDim] = 0.0;
		}
		
		
		Blockage_Div = 0.0;
		// Calculating blockage divergence term.
		for(iDim=0; iDim < nDim; iDim ++){
			Blockage_Div += (U_i[iDim + 1] / U_i[0]) * BGradient[iDim];
		}
		
		// Inserting values into the blockage residual vector.
		if(nDim == 2){
		Blockage_Vector[0] = -(1 / b) * U_i[0] * Blockage_Div;
		Blockage_Vector[1] = -(1 / b) * U_i[1] * Blockage_Div;
		Blockage_Vector[2] = -(1 / b) * U_i[2] * Blockage_Div;
		Blockage_Vector[3] = -(1 / b) * U_i[0] * enthalpy * Blockage_Div;
		}else{
		Blockage_Vector[0] = -(1 / b) * U_i[0] * Blockage_Div;
		Blockage_Vector[1] = -(1 / b) * U_i[1] * Blockage_Div;
		Blockage_Vector[2] = -(1 / b) * U_i[2] * Blockage_Div;
		Blockage_Vector[3] = -(1 / b) * U_i[3] * Blockage_Div;
		Blockage_Vector[4] = -(1 / b) * U_i[0] * enthalpy * Blockage_Div;
		}
		
		// Storing the blockage residual vector at the current node.
		fluidsolver->node[iPoint]->SetBlockage_Source(Blockage_Vector);
	};
}

CBodyForceModelSolver::CBodyForceModelSolver::~CBodyForceModelSolver(){
	
	for(size_t ic=0; ic<3; ic++){
		delete [] proj_vector_axial[ic];
		delete [] proj_vector_tangential[ic];
		delete [] proj_vector_radial[ic];
		delete [] relative_vel[ic];
	}
	delete [] cyl_coordinates[0];
	delete [] cyl_coordinates[1];
	delete [] cyl_coordinates;
	delete [] proj_vector_axial;
	delete [] proj_vector_tangential;
	delete [] proj_vector_radial;
	delete [] relative_vel;
	cout << "Deleted CBodyForceModel container" << endl;
}

su2double CBodyForceModelSolver::vector_dot_product(vector<su2double> v_1, vector<su2double> v_2){
	int n = v_1.size();
	su2double dot_product {0};
	for(int i=0; i<n; i++){
		dot_product += v_1.at(i) * v_2.at(i);
	}
	return dot_product;
}