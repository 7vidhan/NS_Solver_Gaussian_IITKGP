#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>


int nr, ntheta, np, maxiter, inumber_time_steps;
int neighbor[1000000][5];
double T[1000000], Told[1000000], s_t[1000000];
double dr, dtheta, hmin, er_mx;
double r[1000000], theta[1000000];
double omega, r1, r2, sigma, det;
double ainv[5][5], a[5][5], d2phid2[5], a_coeff[1000000][5];


void read_point_data();
void neighbor_points();
void initialise_var();
void source_term();
void coeff_mat();
void boundary_conditions();
double solve();
void local_phi_coefficients(int i, int j, int ij);
void local_phi_coefficients_inv(int np);
double determinant(double a[5][5], int np);
void cofactor(double a[5][5], double np);
void transpose(double a1[5][5], double b1[5][5], double c1);
void d2fd2(int ij);
void find_coeff_A(int i);

int main() {
   int i, j, ij, k, l, iter;
   // no. of internal points in r and theta directions
	nr = 18 ;
   ntheta = 250;
   // domain length in r and theta
   r1 = 10;
	r2 = 20;
	// grid spacing in r and theta
   dr = 0.5;//(r2-r1)/(nr+1);
   dtheta = 0.025;//(2*pi)/(ny+1);
   // find min grid spacing and compute sigma
   hmin = dr;
   FILE *fptr;
   fptr = fopen("GPVidhan.dat","w");

   if (hmin > dtheta) {
		hmin = dtheta;
	}
   sigma = 40*hmin;


	printf ("dr = %lf\n", dr);
	np =5;
   read_point_data();

   neighbor_points();
//   exit(0);
	initialise_var();
   source_term();
   coeff_mat();
   for (iter = 1; iter <= 20000; iter++) {
		boundary_conditions();
		er_mx = solve();
		//printf("%lf %d\n", er_mx, iter);
		if (er_mx > 1e-4) {
			for (i = 1; i <= nr; i++) {
				for (j = 1; j <= ntheta; j++) {
			  		ij  = i-1 + (j-1)*nr;
		      	Told[neighbor[ij][2]] = T[neighbor[ij][2]];
				}
			}
		}
		else {
			break;
		}
	}
	/*for (j = 0; j <= ny+1; j++) {
		i = 0;
		ij = i + j*(nx+2);
      printf ("T[%d] = %lf \n", ij, T[ij]);
	}*/
  	for (j = 0; j <= ntheta+1; j++) {
	   for (i = 0; i <= nr+1; i++) {
            fprintf (fptr,"%lf, %lf, %lf\n", r[ij], theta[ij], T[ij]);
			ij = i + j*(nr+2);
         //printf ("%lf, %lf, %lf\n", x[ij], y[ij], T[ij]);
		}
	}
	return 0;
}

double solve() {
	int i, j, ij, ij1, ij_mx, i_mx, j_mx;
   double  err, er_mx1;
   er_mx1 = 0.;

	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= ntheta; j++) {
		   ij  = i-1 + (j-1)*nr;
         T[neighbor[ij][2]] = (s_t[neighbor[ij][2]]- a_coeff[ij][0]*Told[neighbor[ij][0]] - a_coeff[ij][1]*Told[neighbor[ij][1]] - a_coeff[ij][3]*Told[neighbor[ij][3]] - a_coeff[ij][4]*Told[neighbor[ij][4]])/ a_coeff[ij][2];
			printf("%d %d %d %lf %lf \n", i, j, ij, T[neighbor[ij][2]], Told[neighbor[ij][2]]);
         err = 0;
			err = fabs(T[neighbor[ij][2]] - Told[neighbor[ij][2]]);
         //printf("%lf %lf\n", err, er_mx1);
			if (err > er_mx1) {
				er_mx1 = err;
				ij_mx = neighbor[ij][2];
            i_mx  = i;
				j_mx  = j;
            //printf("wwwwwwwww %lf %d %d %d\n", er_mx1, ij_mx, i, j);
			}

		}
	}
//printf("wwwwwwwww\n");
   //printf("%d %d %d max_err = %lf\n", ij_mx, i_mx, j_mx, er_mx1);
   return (er_mx1);
}


void read_point_data() {
   /* Function to read grid coordinates
   x[0], x[nx+1], y[0], y[ny+1] => boundary point coordinates */
	int i, j, ij;
   for (j = 0; j <= ntheta+1; j++) {
		for (i = 0; i <= nr+1; i++) {
			ij = i + j*(nr+2);
         r[ij] = i*dr;
         theta[ij] = j*dtheta;
         //printf ("%d, x[%d] = %lf, y[%d] = %lf\n", ij, ij, x[ij], ij, y[ij]);
		}
	}
}

double eucdist(int a,int b){
    double xcor,ycor,dist;
    xcor = pow((r[a]*cos(theta[a])-r[b]*cos(theta[b])),2);
    ycor = pow((r[a]*sin(theta[a])-r[b]*sin(theta[b])),2);
    dist = pow((xcor + ycor),0.5);
    //printf(" %lf", dist);
    return(dist);

}

int* sortradist(double* radist){
    int counting,i,m,d,templocation, position;
    counting=(nr+2)*(ntheta+2);
    int locationarray[counting];
    double testarr[(nr+2)*(ntheta+2)],temp;


    for (i=0; i<counting; i++){
        testarr[i] = *(radist + i);
    }
    /*for(m=0;m<counting;m++ ){
        printf(" %lf",testarr[m] );
    }*/
    for (i=0; i<counting; i++){
        locationarray[i]=i;
       // printf(" %d",locationarray[i]);
    }
    for (i=0;i<counting;i++){
        position = i;
        for (d = i + 1; d < counting; d++){
              if (testarr[position] > testarr[d])
                position = d;
            }
    if (position != i)
    {
      temp = testarr[i];
      testarr[i] = testarr[position];
      testarr[position] = temp;
      // Swapping the location of smaller radial distances
      templocation = locationarray[i];
      locationarray[i] = position;
      locationarray[position] = templocation;
    }

    }
    //for(m=0;m<counting;m++ ){
      //  printf(" %d",locationarray[0] );
    //}
   // printf("\n");
   return(locationarray);
}



void neighbor_points() {
	int i, j, ij, ij1,m,counting;
	counting = (nr+2)*(ntheta+2);
	int locationarray[counting];
	int* storeit;
	double radist[counting];
	double *storedist;

	for (j = 1; j<=ntheta; j++) {
		for (i = 1; i<=nr; i++){
		   ij  = i-1 + (j-1)*nr;
		   ij1 = i + j*(nr+2);
         //printf("%d, %d, %d, %d\n", i, j, ij, ij1);
			for (m=0;m<counting;m++){
		   radist[m] = eucdist(ij1,m);
		   }
		   // Sorting the distance from the ij1 for shortest distance
		   storeit = sortradist(&radist);
		   for(m=0; m<counting; m++){
           locationarray[m] = *(storeit+m) ;
		   }



			neighbor[ij][0] = locationarray[1];
		   neighbor[ij][1] = locationarray[2];
		   neighbor[ij][2] = locationarray[0];
		   neighbor[ij][3] = locationarray[3];
		   neighbor[ij][4] = locationarray[4];
		//   printf ("%d %d %d, neighbor(%d,:) =  %d %d %d %d %d\n", i, j, ij1, ij, neighbor[ij][0], neighbor[ij][1], neighbor[ij][2], neighbor[ij][3], neighbor[ij][4]);
		}
	}
}


void initialise_var() {
   /* Function to read grid coordinates
   x[0], x[nx+1], y[0], y[ny+1] => boundary point coordinates */
	int i, j, ij;
   for (j = 0; j <= ntheta+1; j++) {
		for (i = 0; i <= nr+1; i++) {
			ij = i + j*(nr+2);
         T[ij]    = 0;
         Told[ij] = 0;
         //printf ("%d %d %d, T = %lf, Told = %lf\n", i, j,  ij, T[ij], Told[ij]);
		}
	}
}

void source_term() {
	int i, j, ij;
   for (j = 0; j <= ntheta+1; j++) {
		for (i = 0; i <= nr+1; i++) {
			ij = i + j*(nr+2);
         s_t[ij]    = 0;
         //printf ("%d %d %d, S = %lf\n", i, j,  ij, s_t[ij]);
		}
	}
}

void coeff_mat() {
	int i, j, k, l, ij, ij1;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= ntheta; j++) {
			ij = i-1 + (j-1)*nr;
         ij1 = i + j*(nr+2);
         //printf("%d %d %d %lf %lf\n", i, j, ij1, x[ij1], y[ij1]);

			local_phi_coefficients(i, j, ij);

         /*for (k = 0; k <= 4; k++){
				for (l = 0; l <=4; l++) {
			   	printf ("%lf ", a[k][l]);
				}
				printf("\n");
			}*/

         //printf("%d\n", np);
         det = determinant(a, np);

         //printf("***********\n");
			if (det == 0)
   			printf("\nInverse of Entered Matrix is not possible\n");
  			else
  			cofactor(a, np);
      	d2fd2(ij);
         //goto test;
      	find_coeff_A(ij);
         //test :
			for (k = 0; k <= 4; k++) {
         for (l = 0; l <= 4; l++) {
				a[k][l] = 0;
            ainv[k][l] = 0;
			  //printf("%d %d %lf %lf %lf a_coeff[%d][%d] = %lf\n",i, j, y[168], x[i + j*(nx+2)], y[i + j*(nx+2)],  ij, k, a_coeff[ij][k]);
			}
			}
		}
	}

  //exit(0);
}

void local_phi_coefficients(int i, int j, int ij) {
      /*printf("%d %d %d %d %d\n", neighbor[ij][0], neighbor[ij][1], neighbor[ij][2], neighbor[ij][3], neighbor[ij][4]);
      printf("%lf %lf %lf %lf %lf\n", x[neighbor[ij][0]], x[neighbor[ij][1]], x[neighbor[ij][2]], x[neighbor[ij][3]], x[neighbor[ij][4]]);
      printf("%lf %lf %lf %lf %lf\n", y[neighbor[ij][0]], y[neighbor[ij][1]], y[neighbor[ij][2]], y[neighbor[ij][3]], y[neighbor[ij][4]]);
      //exit (0);*/
        a[0][0] = exp(-pow((r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma),2))*exp(-pow((r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma),2));
		a[0][1] = exp(-pow((r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma),2))*exp(-pow((r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma),2));
		a[0][2] = exp(-pow((r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma),2))*exp(-pow((r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma),2));
		a[0][3] = exp(-pow((r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma),2))*exp(-pow((r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma),2));
		a[0][4] = exp(-pow((r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma),2))*exp(-pow((r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma),2));

		a[1][0] = exp(-pow((r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma),2))*exp(-pow((r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma),2));
		a[1][1] = exp(-pow((r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma),2))*exp(-pow((r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma),2));
		a[1][2] = exp(-pow((r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma),2))*exp(-pow((r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma),2));
		a[1][3] = exp(-pow((r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma),2))*exp(-pow((r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma),2));
		a[1][4] = exp(-pow((r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma),2))*exp(-pow((r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma),2));

        a[2][0] = exp(-pow((r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma),2))*exp(-pow((r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma),2));
		a[2][1] = exp(-pow((r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma),2))*exp(-pow((r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma),2));
		a[2][2] = exp(-pow((r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma),2))*exp(-pow((r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma),2));
		a[2][3] = exp(-pow((r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma),2))*exp(-pow((r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma),2));
		a[2][4] = exp(-pow((r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma),2))*exp(-pow((r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma),2));

        a[3][0] = exp(-pow((r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma),2));
		a[3][1] = exp(-pow((r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma),2));
		a[3][2] = exp(-pow((r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma),2));
		a[3][3] = exp(-pow((r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma),2));
		a[3][4] = exp(-pow((r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma),2));

        a[4][0] = exp(-pow((r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma),2));
		a[4][1] = exp(-pow((r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma),2));
		a[4][2] = exp(-pow((r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma),2));
		a[4][3] = exp(-pow((r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma),2));
		a[4][4] = exp(-pow((r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma),2))*exp(-pow((r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma),2));

       int k,l;
       /*for (k=0; k<np;k++){
            for (l=0; l<np; l++){
                printf(" %lf", a[k][l]);
            }
            printf("\n");
		}*/
}

void d2fd2(int ij){
    int l;

	d2phid2[0] = exp(-pow(r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][0]]*cos(theta[neighbor[ij][0]])/sigma,2))*exp(-pow(r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][0]]*sin(theta[neighbor[ij][0]])/sigma,2))*((-4*r[neighbor[ij][2]]/pow(sigma,2)) + (4*r[neighbor[ij][2]]/pow(sigma,4))*(pow(r[neighbor[ij][2]],2) + pow(r[neighbor[ij][0]],2) + r[neighbor[ij][2]]*r[neighbor[ij][0]]*sin(2*(theta[neighbor[ij][2]] - theta[neighbor[ij][0]]))));

	d2phid2[1] = exp(-pow(r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][1]]*cos(theta[neighbor[ij][1]])/sigma,2))*exp(-pow(r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][1]]*sin(theta[neighbor[ij][1]])/sigma,2))*((-4*r[neighbor[ij][2]]/pow(sigma,2)) + (4*r[neighbor[ij][2]]/pow(sigma,4))*(pow(r[neighbor[ij][2]],2) + pow(r[neighbor[ij][1]],2) + r[neighbor[ij][2]]*r[neighbor[ij][1]]*sin(2*(theta[neighbor[ij][2]] - theta[neighbor[ij][1]]))));

	d2phid2[2] = exp(-pow(r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]])/sigma,2))*exp(-pow(r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]])/sigma,2))*((-4*r[neighbor[ij][2]]/pow(sigma,2)) + (4*r[neighbor[ij][2]]/pow(sigma,4))*(pow(r[neighbor[ij][2]],2) + pow(r[neighbor[ij][2]],2) + r[neighbor[ij][2]]*r[neighbor[ij][2]]*sin(2*(theta[neighbor[ij][2]] - theta[neighbor[ij][2]]))));

	d2phid2[3] = exp(-pow(r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][3]]*cos(theta[neighbor[ij][3]])/sigma,2))*exp(-pow(r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][3]]*sin(theta[neighbor[ij][3]])/sigma,2))*((-4*r[neighbor[ij][2]]/pow(sigma,2)) + (4*r[neighbor[ij][2]]/pow(sigma,4))*(pow(r[neighbor[ij][2]],2) + pow(r[neighbor[ij][3]],2) + r[neighbor[ij][2]]*r[neighbor[ij][3]]*sin(2*(theta[neighbor[ij][2]] - theta[neighbor[ij][3]]))));

	d2phid2[4] = exp(-pow(r[neighbor[ij][2]]*cos(theta[neighbor[ij][2]]) - r[neighbor[ij][4]]*cos(theta[neighbor[ij][4]])/sigma,2))*exp(-pow(r[neighbor[ij][2]]*sin(theta[neighbor[ij][2]]) - r[neighbor[ij][4]]*sin(theta[neighbor[ij][4]])/sigma,2))*((-4*r[neighbor[ij][2]]/pow(sigma,2)) + (4*r[neighbor[ij][2]]/pow(sigma,4))*(pow(r[neighbor[ij][2]],2) + pow(r[neighbor[ij][4]],2) + r[neighbor[ij][2]]*r[neighbor[ij][4]]*sin(2*(theta[neighbor[ij][2]] - theta[neighbor[ij][4]]))));

     /*for (l=0; l<5; l++){
                printf(" %lf", d2phid2[l]);
            }
            printf("\n");*/
}


void find_coeff_A(int ij2) {
	//printf("%d \n", ij2);
	int l;
	a_coeff[ij2][0] = d2phid2[0]*ainv[0][0] + d2phid2[1]*ainv[1][0] + d2phid2[2]*ainv[2][0] + d2phid2[3]*ainv[3][0] + d2phid2[4]*ainv[4][0];
	a_coeff[ij2][1] = d2phid2[0]*ainv[0][1] + d2phid2[1]*ainv[1][1] + d2phid2[2]*ainv[2][1] + d2phid2[3]*ainv[3][1] + d2phid2[4]*ainv[4][1];
	a_coeff[ij2][2] = d2phid2[0]*ainv[0][2] + d2phid2[1]*ainv[1][2] + d2phid2[2]*ainv[2][2] + d2phid2[3]*ainv[3][2] + d2phid2[4]*ainv[4][2];
	a_coeff[ij2][3] = d2phid2[0]*ainv[0][3] + d2phid2[1]*ainv[1][3] + d2phid2[2]*ainv[2][3] + d2phid2[3]*ainv[3][3] + d2phid2[4]*ainv[4][3];
	a_coeff[ij2][4] = d2phid2[0]*ainv[0][4] + d2phid2[1]*ainv[1][4] + d2phid2[2]*ainv[2][4] + d2phid2[3]*ainv[3][4] + d2phid2[4]*ainv[4][4];
    /*for (l=0; l<5; l++){
        printf(" %lf", a_coeff[ij2][l]);

    }*/
    printf("\n");
}


void boundary_conditions() {
	int i, j, ij;
	for (j = 0; j <= ntheta+1; j++) {
		i = 0;
		ij = i + j*(nr+2);
		T[ij] = 1;
      Told[ij] = 1;
      //printf ("T[%d] = %lf \n", ij, T[ij]);
		i = nr+1;
		ij = i + j*(nr+2);
		T[ij] = 50;
      Told[ij] = 50;
      //printf ("T[%d] = %lf \n", ij, T[ij]);
	}
	}



/*For calculating Determinant of the Matrix */
double determinant(double a[5][5], int k)
{
  double s = 1, det, b[5][5];
  int i, j, m, n, c;
  det = 0;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 	 //printf ("det = %lf", det);

    return (det);
}

void cofactor(double num[5][5], double f)
{
 double b[5][5], fac[5][5];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}

void transpose(double num[5][5], double fac[5][5], double r)
{
  int i, j;
  double b[5][5],  d;

  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        ainv[i][j] = b[i][j] / d;
        }
    }
   /*printf("\n\n\nThe inverse of matrix is : \n");

   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%f", ainv[i][j]);
        }
    printf("\n");
     }*/
}



