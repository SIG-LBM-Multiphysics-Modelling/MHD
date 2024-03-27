#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = true;
/// Flow quantities
const vector<int> cx = {0, 1, 0, -1,  0, 1, -1, -1,  1},
									cy = {0, 0, 1,  0, -1, 1,  1, -1, -1};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
const int nx = 512, ny = nx, np = 9;
const double U_ref = 2., Reynolds = 10000., Tmax = 100., Scale_time = 1E-4;
const double Scale_length = 2.*M_PI/(nx-1), Scale_velocity = Scale_length/Scale_time, Scale_enstrophy = 1./pow(Scale_time,2), Scale_energy = pow(Scale_velocity,2), Scale_dissipation = Scale_enstrophy*pow(Scale_length,2)/Scale_time, Scale_viscosity = pow(Scale_length,2)/Scale_time;
const double v0 = U_ref/Scale_velocity, ni = U_ref*2.*M_PI/Reynolds/Scale_viscosity, tau = ni*3.+0.5, omega = 1./tau, omega1 = 1.-omega;
const double cs2 = 1./3., cs4 = cs2*cs2, cs6 = cs4*cs2, cs8 = cs6*cs2, rho0 = 1.;
vector<double> f1(nx*ny*np,0.), f2(nx*ny*np,0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), temp_pop(np,0.), vorticity(nx*ny, 0.);
double kin_energy, kin_enstrophy, max_vorticity;
/// MHD quantities
const int npMHD = 5;
const vector<double> wMHD = {1./3., 1./6., 1./6., 1./6., 1./6.};
const double Prm = 1., eta = ni/Prm, tauMHD = eta*3+.5, omegaMHD = 1./tauMHD, omega1MHD = 1.-omegaMHD, b0 = v0;
vector<double> bx(nx*ny, 0.), by(nx*ny, 0.), g1x(nx*ny*npMHD,0.), g2x(nx*ny*npMHD,0.), g1y(nx*ny*npMHD,0.), g2y(nx*ny*npMHD,0.), temp_pop_gx(npMHD,0.), temp_pop_gy(npMHD,0.), current(nx*ny, 0.);
double mag_energy, mag_enstrophy, Jmax;
/// Scaling factor et al
const int nsteps = int(Tmax/Scale_time)+1, n_out = int(0.1/Scale_time);
double A, B, C, R, U, V, ftemp, BX, BY, JZ, geq;
double k0, k1, k2, k3, k4, k5, k6, k7, k8;
double r0, r1, r2, r3, r4, r5, r6, r7, r8;
double CX, CY, U2, V2, UV, U3, V3, BX2, BY2, U2V2, third_order, fourth_order, D, E, F, termMHD, BXY;
double value;
int id, idn, newx, newy, id1, id2;
double gradx_of_v, grady_of_u, current_vorticity;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	/// Write density difference
	output_file << "SCALARS density float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
      output_file << rho[X*ny+Y]-rho0<< "\n";

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
			output_file << u[X*ny+Y]*Scale_velocity << " " << v[X*ny+Y]*Scale_velocity << " 0\n";

	output_file << "VECTORS magnetic_field float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
			output_file << bx[X*ny+Y]*Scale_velocity << " " << by[X*ny+Y]*Scale_velocity << " 0\n";

	output_file << "SCALARS current float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		   output_file << current[X*ny+Y]/Scale_time << "\n";
			
	output_file << "SCALARS vorticity float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
			output_file << vorticity[X*ny+Y]/Scale_time<< "\n";

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double D, E, F, termMHD;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			R = rho[id] = rho0;
			U = u[id] = -sin(2.*M_PI*y/(ny-1))*v0;
      V = v[id] = sin(2.*M_PI*x/(nx-1))*v0;
			BX = bx[id] = -sin(2.*M_PI*y/(ny-1))*b0;
			BY = by[id] = sin(4.*M_PI*x/(nx-1))*b0;
		  C = -1.5*(U*U+V*V);
			D = BX*BX+BY*BY;
			for(int k=0; k<np;k++)
			{
        A = U*cx[k]+V*cy[k];
        B = 4.5*A*A;
				E = pow(sqrt(cx[k]*cx[k]+cy[k]*cy[k]), 2);
				F = pow(cx[k]*BX+cy[k]*BY,2);
				termMHD = wf[k]/(2.*cs2*cs2)*(0.5*D*E-F);
        third_order = 1./(2.*cs6)*((cx[k]*cx[k]-cs2)*cy[k]*U*U*V+(cy[k]*cy[k]-cs2)*cx[k]*U*V*V);
        fourth_order = 1./(4.*cs8)*((cx[k]*cx[k]-cs2)*(cy[k]*cy[k]-cs2)*U*U*V*V);
        f1[id*np+k] = f2[id*np+k] = wf[k]*R*(1.+3.*A+B+C+third_order+fourth_order)+termMHD;
        if(k<npMHD)
        {
					g1x[id*npMHD+k] = g2x[id*npMHD+k] = wMHD[k]*(BX+3.*cy[k]*(V*BX-BY*U));
					g1y[id*npMHD+k] = g2y[id*npMHD+k] = wMHD[k]*(BY+3.*cx[k]*(U*BY-BX*V));
				}
			}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_vorticity()
{
	kin_enstrophy = max_vorticity = 0.;
  for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			gradx_of_v = grady_of_u = 0.;
			for(int k=1; k<np; k++)
			{
				newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
        idn = newx*ny+newy;
				gradx_of_v += 3.*wf[k]*cx[k]*v[idn];
				grady_of_u += 3.*wf[k]*cy[k]*u[idn];
			}
			vorticity[id] = gradx_of_v-grady_of_u;
			kin_enstrophy += vorticity[id]*vorticity[id];
			current_vorticity = fabs(vorticity[id]);
			if(current_vorticity>max_vorticity)
				max_vorticity = current_vorticity;
		}
	kin_enstrophy /= nx*ny;
}
///----------------------------------------------------------------------------------------------------------------------------------
int algoLB()
{
	int check = 0;
	Jmax = kin_energy = mag_energy = mag_enstrophy = 0.;
  for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			U = V = R = BX = BY = JZ = 0.;
			for(int k=0; k<np; k++)
			{
				ftemp = temp_pop[k] = f1[id*np+k];
				R += ftemp;
				U += ftemp*cx[k];
				V += ftemp*cy[k];
				if(k<npMHD)
				{
					temp_pop_gx[k] = g1x[id*npMHD+k];
					temp_pop_gy[k] = g1y[id*npMHD+k];
					BX += temp_pop_gx[k];
					BY += temp_pop_gy[k];
					JZ += cx[k]*temp_pop_gy[k]-cy[k]*temp_pop_gx[k];
				}
			}
			bx[id] = BX;
			by[id] = BY;
			rho[id] = R;
			U /= R;
			V /= R;
			JZ -= 2.*(U*BY-V*BX);
			JZ *= -3*omegaMHD;
			current[id] = JZ;
			JZ = fabs(JZ);
			if(JZ>Jmax)
				Jmax = JZ;
			u[id] = U;
			v[id] = V;
			U2 = U*U;
			V2 = V*V;
			UV = U*V;
			U2V2 = U2*V2;
			BX2 = BX*BX;
			BY2 = BY*BY;
			BXY = BX*BY;

			r4 = temp_pop[1]-temp_pop[2]+temp_pop[3]-temp_pop[4];
			r5 = temp_pop[5]-temp_pop[6]+temp_pop[7]-temp_pop[8];
			k4 = r4-R*(U2-V2);
			k5 = r5-R*UV;

      k0 = R;
			k3 = 2.*R*cs2;
      k4 = omega1*k4+omega*(BY2-BX2);
			k5 = omega1*k5+omega*(-BXY);
			k6 = 0.5*V*(BX2-BY2)+2.*BXY*U;
			k7 = 0.5*U*(BY2-BX2)+2.*BXY*V;
			k8 = R*cs4+0.5*(BX2-BY2)*(U2-V2)-4.*UV*BXY;

      r0 = k0;
			r1 = R*U;
			r2 = R*V;
			r3 = k3+R*(U2+V2);
			r4 = k4+R*(U2-V2);
			r5 = k5+R*UV;
			r6 = k6+2.*U*k5+0.5*V*(k3+k4)+R*U2*V;
			r7 = k7+2.*V*k5+0.5*U*(k3-k4)+R*U*V2;
			r8 = k8+2.*(U*k7+V*k6)+4.*UV*k5-0.5*k4*(U2-V2)+0.5*k3*(U2+V2)+R*U2V2;

      temp_pop[0] = r0-r3+r8;
      temp_pop[1] = 0.5*(r1-r7-r8) + 0.25*(r3+r4);
      temp_pop[2] = 0.5*(r2-r6-r8) + 0.25*(r3-r4);
      temp_pop[3] = 0.5*(-r1+r7-r8) + 0.25*(r3+r4);
      temp_pop[4] = 0.5*(-r2+r6-r8) + 0.25*(r3-r4);
      temp_pop[5] = 0.25*(r5+r6+r7+r8);
      temp_pop[6] = 0.25*(r6-r5-r7+r8);
      temp_pop[7] = 0.25*(r5-r6-r7+r8);
      temp_pop[8] = 0.25*(r7-r6-r5+r8);
			for(int k=0; k<np; k++)
			{
        newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
        idn = newx*ny+newy;
				f2[idn*np+k] = temp_pop[k];
				if(k<npMHD)
				{
					g2x[idn*npMHD+k] = omega1MHD*temp_pop_gx[k] + omegaMHD*wMHD[k]*(BX+3.*cy[k]*(V*BX-BY*U));
					g2y[idn*npMHD+k] = omega1MHD*temp_pop_gy[k] + omegaMHD*wMHD[k]*(BY+3.*cx[k]*(U*BY-BX*V));
				}
			}
			if(fabs(U)>1.)
			  check = 1;
			kin_energy += U2+V2;
			mag_energy += BX*BX+BY*BY;
			mag_enstrophy += JZ*JZ;
		}
	kin_energy *= 0.5/(nx*ny);
	mag_energy *= 0.5/(nx*ny);
	mag_enstrophy /= nx*ny;
  return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  system("mkdir vtk_fluid");
  FILE *data_output = fopen("data.txt","wt");
	int check_mach = 0;
	double dissipation_rate;
	initial_state();
	printf("v0 = %e\n", v0);
	clock_t c_start = std::clock();
  for(int i=0; i<nsteps; i++)
  {
    check_mach = algoLB();
    f1 = f2;
    g1x = g2x;
    g1y = g2y;
		compute_vorticity();
		dissipation_rate = ni*kin_enstrophy+eta*mag_enstrophy;
		fprintf(data_output,"%lf    %e		%e    %e		%e		%e		%e		%e		%e\n", (double)i*Scale_time, Jmax/Scale_time, max_vorticity/Scale_time, kin_energy*Scale_energy, mag_energy*Scale_energy, (mag_energy+kin_energy)*Scale_energy, mag_enstrophy*Scale_enstrophy, kin_enstrophy*Scale_enstrophy, dissipation_rate*Scale_dissipation);
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
    if(i%n_out==0)
      printf("Iteration=%lf of %lf. Jmax=%lf\n", (double)i*Scale_time, (double)nsteps*Scale_time, Jmax/Scale_time);
		if(check_mach==1) // check the Mach number...if too high, it exits!
      goto labelA;
  }
	labelA:
	clock_t c_end = std::clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  fclose(data_output);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
