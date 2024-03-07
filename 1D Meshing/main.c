#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>


int nx, inumber_points, inumber_points_boundary, maxiter, inumber_time_steps;
int neighbor[100][3];
int ib[100], inumber_neighbors[100];
double u[100], v[100], w[100], h[100], p[100], uhat[100], vhat[100], what[100], divg[100], diag[100];
double df1dx[100][7], df1dy[100][7], df1dz[100][7], df2dx[100][7], df2dy[100][7], df2dz[100][7], pcoef[100][7];
double dx, dy, dz, dxsq, dysq, dzsq;
double x[100], y[100], z[100], den[100];
double amu, rho, time_step, omega, lx, sigma;
double ainv[9], a[9], d2phidx2[3], a_coeff[100][3];


void read_point_data();
void neighbor_points();
void local_phi_coefficients(int i);
void local_phi_coefficients_inv();
void d2fdx2(int i);
void find_coeff_A(int i);

int main() {
   int i, j, k;
	inumber_points = 127;
   lx = 1;

   dx = lx/(inumber_points+1);
   sigma = 40*dx;
   inumber_points_boundary = 2;
   nx = inumber_points;
	printf ("dx = %lf\n", dx);
   read_point_data();
   neighbor_points();
   for (i = 1; i <= nx; i++) {
		local_phi_coefficients(i);
		for (k = 0; k <= 8; k++) {
				//printf("pt = %d, a[%d] = %lf\n", i, k, a[k]);
		}
      local_phi_coefficients_inv();
		for (k = 0; k <= 8; k++) {
				//printf("pt = %d, ainv[%d] = %lf\n", i, k, ainv[k]);
		}
      d2fdx2(i);
		for (k = 0; k <= 2; k++) {
				//printf("pt = %d, d2phidx2[%d] = %lf\n", i, k, d2phidx2[k]);
		}
      find_coeff_A(i);
      for (k = 0; k <= 2; k++) {
				printf("a_coeff[%d][%d] = %lf\n", i, k, a_coeff[i][k]);
		}

	}

	return 0;
}

void read_point_data(){
	int i;
	for (i = 1; i <= inumber_points; i++) {
		x[i] = i*dx;
	}
	// boundary points
	x[0] = x[1] - dx;
   x[inumber_points+1] = x[inumber_points] + dx;
   for (i = 0; i <= inumber_points+1; i++) {
		printf ("x(%d) = %lf\n", i, x[i]);
	}
}

void property_data(){

}


void neighbor_points(){
	int i;
	for (i = 1; i<=nx; i++){
		neighbor[i][0] = i-1;
		neighbor[i][1] = i;
		neighbor[i][2] = i+1;
      printf ("neighbor(%d,:) = %d %d %d\n", i, neighbor[i][0], neighbor[i][1], neighbor[i][2]);
	}
}

void local_phi_coefficients(int i){
      a[0] = exp(-pow((x[neighbor[i][0]] - x[neighbor[i][0]])/sigma,2));
		a[1] = exp(-pow((x[neighbor[i][0]] - x[neighbor[i][1]])/sigma,2));
      //printf ("%lf, %lf\n",x[neighbor[i][0]], x[neighbor[i][2]]);
		a[2] = exp(-pow((x[neighbor[i][0]] - x[neighbor[i][2]])/sigma,2));
		a[3] = exp(-pow((x[neighbor[i][1]] - x[neighbor[i][0]])/sigma,2));
		a[4] = exp(-pow((x[neighbor[i][1]] - x[neighbor[i][1]])/sigma,2));
		a[5] = exp(-pow((x[neighbor[i][1]] - x[neighbor[i][2]])/sigma,2));
		a[6] = exp(-pow((x[neighbor[i][2]] - x[neighbor[i][0]])/sigma,2));
		a[7] = exp(-pow((x[neighbor[i][2]] - x[neighbor[i][1]])/sigma,2));
		a[8] = exp(-pow((x[neighbor[i][2]] - x[neighbor[i][2]])/sigma,2));
}

void local_phi_coefficients_inv(){
	double det;
	det = a[0] * a[4] * a[8] + a[3] * a[7] * a[2] + a[6] * a[1] * a[5] - a[0] * a[7] * a[5] - a[6] * a[4] * a[2] - a[3] * a[1] * a[8];
	ainv[0] = (a[4] * a[8] - a[7] * a[5]) / det;
	ainv[1] = (a[2] * a[7] - a[1] * a[8]) / det;
	ainv[2] = (a[1] * a[5] - a[2] * a[4]) / det;
	ainv[3] = (a[5] * a[6] - a[3] * a[8]) / det;
	ainv[4] = (a[0] * a[8] - a[6] * a[2]) / det;
	ainv[5] = (a[2] * a[3] - a[0] * a[5]) / det;
	ainv[6] = (a[3] * a[7] - a[6] * a[4]) / det;
	ainv[7] = (a[6] * a[1] - a[0] * a[7]) / det;
	ainv[8] = (a[0] * a[4] - a[1] * a[3]) / det;
}

void d2fdx2(int i){
	printf ("%lf, %lf,%lf, %lf, %lf\n",x[neighbor[i][1]] , x[neighbor[i][0]],exp(-pow((x[neighbor[i][1]] - x[neighbor[i][0]])/sigma,2)), -2/pow(sigma,2) , 4*pow(x[neighbor[i][1]] - x[neighbor[i][0]],2)/pow(sigma,4));
	d2phidx2[0] = exp(-pow((x[neighbor[i][1]] - x[neighbor[i][0]])/sigma,2))*(-2/pow(sigma,2) + 4*pow(x[neighbor[i][1]] - x[neighbor[i][0]],2)/pow(sigma,4));
   d2phidx2[1] = -2/pow(sigma,2);
	d2phidx2[2]	= exp(-pow((x[neighbor[i][1]] - x[neighbor[i][2]])/sigma,2))*(-2/pow(sigma,2) + 4*pow(x[neighbor[i][1]] - x[neighbor[i][2]],2)/pow(sigma,4));

}


void find_coeff_A(int i){
	a_coeff[i][0] = d2phidx2[0]*ainv[0] + d2phidx2[1]*ainv[3] + d2phidx2[2]*ainv[6];
   a_coeff[i][1] = d2phidx2[0]*ainv[1] + d2phidx2[1]*ainv[4] + d2phidx2[2]*ainv[7];
	a_coeff[i][2] = d2phidx2[0]*ainv[2] + d2phidx2[1]*ainv[5] + d2phidx2[2]*ainv[8];
}
