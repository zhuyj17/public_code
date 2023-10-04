#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.1415926
#define nx 64
#define ny 64
#define nz 64
#define N nx*ny*nz
 
double d_mesh=1e-9, itr_limit=20.0;
double voltage_distribute_main[N],epsilon_main[N];
double E_2_main[N], phi_main[N], ita_main[N], norm_cri_main[N];
double H[N];
int kind_main[N];  // if kind=1 means this is a particle; if kind=0 means this is matrix




void generate_electrical_field(double *voltage_distribute, double *E_2)
{


    double *E_x= malloc(nx * ny * nz * sizeof(double));
    double *E_y= malloc(nx * ny * nz * sizeof(double));
    double *E_z= malloc(nx * ny * nz * sizeof(double));
    double *E_x_new= malloc(nx * ny * nz * sizeof(double));
    double *E_y_new= malloc(nx * ny * nz * sizeof(double));
    double *E_z_new= malloc(nx * ny * nz * sizeof(double));
    
    

    for(int i=0; i<nx-1; i++)
    {   
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nz-1; k++)
            {
                *(E_x+(i)*ny*nz+(j)*nz+(k)) = fabs((*(voltage_distribute+(i+1)*ny*nz+j*nz+k)-\
                *(voltage_distribute+i*ny*nz+j*nz+k))/(d_mesh));
                *(E_y+(i)*ny*nz+(j)*nz+(k)) = fabs((*(voltage_distribute+i*ny*nz+(j+1)*nz+k)-\
                *(voltage_distribute+i*ny*nz+j*nz+k))/(d_mesh));
                *(E_z+(i)*ny*nz+(j)*nz+(k)) = fabs((*(voltage_distribute+i*ny*nz+j*nz+(k+1))-\
                *(voltage_distribute+i*ny*nz+j*nz+k))/(d_mesh));

            }         
        }
    } 

    for(int i=0; i<nx-1; i++)
    {   
        for(int j=0; j<ny-1;j++)
        {
            *(E_x+(i)*ny*nz+(j)*nz+(nz-1)) = fabs((*(voltage_distribute+(i+1)*ny*nz+j*nz+nz-1)-\
            *(voltage_distribute+i*ny*nz+j*nz+nz-1))/(d_mesh));
            *(E_y+(i)*ny*nz+(j)*nz+(nz-1)) = fabs((*(voltage_distribute+i*ny*nz+(j+1)*nz+nz-1)-\
            *(voltage_distribute+i*ny*nz+j*nz+nz-1))/(d_mesh));
            *(E_z+(i)*ny*nz+(j)*nz+(nz-1)) = fabs((*(voltage_distribute+i*ny*nz+j*nz+0)-\
            *(voltage_distribute+i*ny*nz+j*nz+nz-1))/(d_mesh));
        }
    }

    for(int i=0; i<nx-1; i++)
    {   
        for(int k=0; k<nz-1;k++)
        {
            *(E_x+(i)*ny*nz+(ny-1)*nz+(k)) = fabs((*(voltage_distribute+(i+1)*ny*nz+(ny-1)*nz+k)-\
            *(voltage_distribute+i*ny*nz+(ny-1)*nz+k))/(d_mesh));
            *(E_y+(i)*ny*nz+(ny-1)*nz+(k)) = fabs((*(voltage_distribute+i*ny*nz+(0)*nz+k)-\
            *(voltage_distribute+i*ny*nz+(ny-1)*nz+k))/(d_mesh));
            *(E_z+(i)*ny*nz+(ny-1)*nz+(k)) = fabs((*(voltage_distribute+i*ny*nz+(ny-1)*nz+(k+1))-\
            *(voltage_distribute+i*ny*nz+(ny-1)*nz+k))/(d_mesh));
        }
    }


    for(int j=0; j<ny-1; j++)
    {   
        for(int k=0; k<nz-1;k++)
        {
            *(E_x+(nx-1)*ny*nz+(j)*nz+(k)) = *(E_x+(nx-2)*ny*nz+(j)*nz+(k)); 
            *(E_y+(nx-1)*ny*nz+(j)*nz+(k)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+(j+1)*nz+k)-\
            *(voltage_distribute+(nx-1)*ny*nz+j*nz+k))/(d_mesh));
            *(E_z+(nx-1)*ny*nz+(j)*nz+(k)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+j*nz+(k+1))-\
            *(voltage_distribute+(nx-1)*ny*nz+j*nz+k))/(d_mesh));
        }
    }



    for(int i=0; i<nx-1; i++)
    {   
        *(E_x+(i)*ny*nz+(ny-1)*nz+(nz-1)) = fabs((*(voltage_distribute+(i+1)*ny*nz+(ny-1)*nz+(nz-1))-\
        *(voltage_distribute+i*ny*nz+(ny-1)*nz+(nz-1)))/(d_mesh));
        *(E_y+(i)*ny*nz+(ny-1)*nz+(nz-1)) = fabs((*(voltage_distribute+i*ny*nz+(0)*nz+(nz-1))-\
        *(voltage_distribute+i*ny*nz+(ny-1)*nz+(nz-1)))/(d_mesh));
        *(E_z+(i)*ny*nz+(ny-1)*nz+(nz-1)) = fabs((*(voltage_distribute+i*ny*nz+(ny-1)*nz+(0))-\
        *(voltage_distribute+i*ny*nz+(ny-1)*nz+(nz-1)))/(d_mesh));
    
    }

    for(int j=0; j<ny-1; j++)
    {   
        *(E_x+(nx-1)*ny*nz+(j)*nz+(nz-1)) = *(E_x+(nx-2)*ny*nz+(j)*nz+(nz-1));
        *(E_y+(nx-1)*ny*nz+(j)*nz+(nz-1)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+(j+1)*nz+(nz-1))-\
        *(voltage_distribute+(nx-1)*ny*nz+j*nz+(nz-1)))/(d_mesh));
        *(E_z+(nx-1)*ny*nz+(j)*nz+(nz-1)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+j*nz+(0))-\
        *(voltage_distribute+(nx-1)*ny*nz+j*nz+(nz-1)))/(d_mesh));
    
    }


    for(int k=0; k<nz-1; k++)
    {   
        *(E_x+(nx-1)*ny*nz+(ny-1)*nz+(k)) = *(E_x+(nx-2)*ny*nz+(ny-1)*nz+(k));
        *(E_y+(nx-1)*ny*nz+(ny-1)*nz+(k)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+(0)*nz+k) -\
        *(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+k))/(d_mesh));
        *(E_z+(nx-1)*ny*nz+(ny-1)*nz+(k)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+(k+1)) -\
        *(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+k))/(d_mesh));
    
    }



    *(E_x+(nx-1)*ny*nz+(ny-1)*nz+(nz-1)) = *(E_x+(nx-2)*ny*nz+(ny-1)*nz+(nz-1)) ;
    *(E_y+(nx-1)*ny*nz+(ny-1)*nz+(nz-1))  = fabs((*(voltage_distribute+(nx-1)*ny*nz+(0)*nz+nz-1)-\
    *(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+nz-1))/(d_mesh));
    *(E_z+(nx-1)*ny*nz+(ny-1)*nz+(nz-1)) = fabs((*(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+0)-\
    *(voltage_distribute+(nx-1)*ny*nz+(ny-1)*nz+nz-1))/(d_mesh));



    for(int i=0; i<nx; i++)
    {   
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                *(E_x_new+(i)*ny*nz+(j)*nz+(k)) = *(E_x+(i)*ny*nz+(j)*nz+(k));
                *(E_y_new+(i)*ny*nz+(j)*nz+(k)) = *(E_y+(i)*ny*nz+(j)*nz+(k));
                *(E_z_new+(i)*ny*nz+(j)*nz+(k)) = *(E_z+(i)*ny*nz+(j)*nz+(k));

            }         
        }
    } 



    for(int i=0; i<nx-1; i++)
    {   
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nz-1; k++)
            {


                *(E_x_new+(i)*ny*nz+(j)*nz+(k)) = 0.25*(*(E_x+(i)*ny*nz+(j+1)*nz+(k+1))+*(E_x+(i)*ny*nz+(j+1)*nz+(k))+\
                *(E_x+(i)*ny*nz+(j)*nz+(k+1))+*(E_x+(i)*ny*nz+(j)*nz+(k)));
                *(E_y_new+(i)*ny*nz+(j)*nz+(k)) = 0.25*(*(E_y+(i+1)*ny*nz+(j)*nz+(k+1))+*(E_y+(i+1)*ny*nz+(j)*nz+(k))+\
                *(E_y+(i)*ny*nz+(j)*nz+(k+1))+*(E_y+(i)*ny*nz+(j)*nz+(k)));
                *(E_z_new+(i)*ny*nz+(j)*nz+(k)) = 0.25*(*(E_z+(i+1)*ny*nz+(j+1)*nz+(k))+*(E_z+(i+1)*ny*nz+(j)*nz+(k))+\
                *(E_z+(i)*ny*nz+(j+1)*nz+(k))+*(E_z+(i)*ny*nz+(j)*nz+(k)));

            }         
        }
    }
    
    /**************************************************************************************/

    for(int i=0; i<nx; i++)
    {   
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                *(E_2+i*ny*nz+j*nz+k) = pow(*(E_x_new+i*ny*nz+j*nz+k),2)+pow(*(E_y_new+i*ny*nz+j*nz+k),2)+\
                pow(*(E_z_new+i*ny*nz+j*nz+k),2);
            }         
        }
    }

    free(E_x);
    free(E_y);
    free(E_z);
    free(E_x_new);
    free(E_y_new);
    free(E_z_new);

}



void electrical_calculate(double *voltage, double *epsilon, int step, int *debug)
{
    double itr=20000.0;
    int N_max=20000, N_step=0;


    while (itr>itr_limit && N_step<N_max)
    {
        itr=0.0;
        double all=0.0;
        for(int i=1;i<nx-1;i++)
        {
            for(int j=1;j<ny-1;j++)
            {
                for(int k=1;k<nz-1;k++)
                {
                    double origin = *(voltage+i*ny*nz+j*nz+k);
                    double a0, a1, a2, a3, a4, a5, a6, aa;
                
                    a1 = 0.25*(*(epsilon+(i)*ny*nz+(j)*nz+(k))+*(epsilon+(i)*ny*nz+(j-1)*nz+(k))+\
                    *(epsilon+(i)*ny*nz+(j)*nz+(k-1))+*(epsilon+(i)*ny*nz+(j-1)*nz+(k-1)));
                    a2 = 0.25*(*(epsilon+(i)*ny*nz+(j)*nz+(k))+*(epsilon+(i-1)*ny*nz+(j)*nz+(k))+\
                    *(epsilon+(i)*ny*nz+(j)*nz+(k-1))+*(epsilon+(i-1)*ny*nz+(j)*nz+(k-1)));
                    a3 = 0.25*(*(epsilon+(i-1)*ny*nz+(j)*nz+(k))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k))+\
                    *(epsilon+(i-1)*ny*nz+(j)*nz+(k-1))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k-1)));
                    a4 = 0.25*(*(epsilon+(i)*ny*nz+(j-1)*nz+(k))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k))+\
                    *(epsilon+(i)*ny*nz+(j-1)*nz+(k-1))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k-1)));
                    a5 = 0.25*(*(epsilon+(i)*ny*nz+(j)*nz+(k))+*(epsilon+(i-1)*ny*nz+(j)*nz+(k))+\
                    *(epsilon+(i)*ny*nz+(j-1)*nz+(k))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k)));
                    a6 = 0.25*(*(epsilon+(i)*ny*nz+(j)*nz+(k-1))+*(epsilon+(i-1)*ny*nz+(j)*nz+(k-1))+\
                    *(epsilon+(i)*ny*nz+(j-1)*nz+(k-1))+*(epsilon+(i-1)*ny*nz+(j-1)*nz+(k-1)));

                    a0 = a1+a2+a3+a4+a5+a6;
                
                    aa = 1/a0*(a1* *(voltage+(i+1)*ny*nz+(j)*nz+(k))+a2* *(voltage+(i)*ny*nz+(j+1)*nz+(k))+\
                    a3* *(voltage+(i-1)*ny*nz+(j)*nz+(k))+a4* *(voltage+(i)*ny*nz+(j-1)*nz+(k))+\
                    a5* *(voltage+(i)*ny*nz+(j)*nz+(k+1))+a6* *(voltage+(i)*ny*nz+(j)*nz+(k-1)));
                
                    *(voltage+(i)*ny*nz+(j)*nz+(k)) += 1*(aa- *(voltage+(i)*ny*nz+(j)*nz+(k)));

                
                    all += *(voltage+(i)*ny*nz+(j)*nz+(k));
                
                    itr += fabs(*(voltage+(i)*ny*nz+(j)*nz+(k))-origin);
                }
            }
        }

        //printf("%f",all);

        itr /= all/(ny-1)/(nx-1)/(nz-1); 

        if (isnan(itr) == 1)
        {
            *debug =0;
            break;
        }



        for(int i=1; i<nx-1; i++)
        {
            for(int j=1; j<ny-1; j++)
            {
                *(voltage+(i)*ny*nz+(j)*nz+(0)) = *(voltage+(i)*ny*nz+(j)*nz+(1));
                *(voltage+(i)*ny*nz+(j)*nz+(nz-1)) = *(voltage+(i)*ny*nz+(j)*nz+(nz-2));
            }
        }

        for(int i=1; i<nx-1; i++)
        {
            for(int k=1; k<nz-1; k++)
            {
                *(voltage+(i)*ny*nz+(0)*nz+(k)) = *(voltage+(i)*ny*nz+(1)*nz+(k));
                *(voltage+(i)*ny*nz+(ny-1)*nz+(k)) = *(voltage+(i)*ny*nz+(ny-2)*nz+(k));
            }
        }

        for(int i=1; i<nx-1; i++)
        {
            *(voltage+(i)*ny*nz+(0)*nz+(0)) = 0.5*(*(voltage+(i)*ny*nz+(1)*nz+(0))+\
            *(voltage+(i)*ny*nz+(0)*nz+(1)));
            *(voltage+(i)*ny*nz+(ny-1)*nz+(0)) = 0.5*(*(voltage+(i)*ny*nz+(ny-2)*nz+(0))+\
            *(voltage+(i)*ny*nz+(ny-1)*nz+(1)));
            *(voltage+(i)*ny*nz+(0)*nz+(nz-1)) = 0.5*(*(voltage+(i)*ny*nz+(1)*nz+(nz-1))+\
            *(voltage+(i)*ny*nz+(0)*nz+(nz-2)));
            *(voltage+(i)*ny*nz+(ny-1)*nz+(nz-1)) = 0.5*(*(voltage+(i)*ny*nz+(ny-2)*nz+(nz-1))+\
            *(voltage+(i)*ny*nz+(ny-1)*nz+(nz-2)));
        }

        
        printf("the itr of step %d error is %5f.\n",step,itr);
        


            
        N_step += 1;
    }
    


}



void phase_cal(double *E_2, double *epsilon, double *norm_cri, double *ita, int *kind, double deltat, \
double epsilonb, double epsilonmatrix, double epsilonpartical, double h, double alpha, double r, double l, double *H, int *debug)
{

    double *f_sep_grad= malloc(nx * ny * nz * sizeof(double));
    double *f_elec_grad= malloc(nx * ny * nz * sizeof(double));


    // initialize
    double epsilon0 = 8.854e-12;

    

    // initialize the gradient matrix
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                if (*(E_2+(i)*ny*nz+(j)*nz+(k)) >= pow(*(norm_cri+(i)*ny*nz+(j)*nz+(k)),2))
                {
                    *(H+(i)*ny*nz+(j)*nz+(k)) = 1.0;
                }
                else
                {
                    *(H+(i)*ny*nz+(j)*nz+(k)) = 0.0;
                }
                *(f_sep_grad+(i)*ny*nz+(j)*nz+(k)) = 2*alpha* *(ita+(i)*ny*nz+(j)*nz+(k))*(1- *(ita+(i)*ny*nz+(j)*nz+(k)))*\
                (1-2* *(ita+(i)*ny*nz+(j)*nz+(k)));
                if (*(kind+(i)*ny*nz+(j)*nz+(k)) == 0)
                {
                    *(f_elec_grad+(i)*ny*nz+(j)*nz+(k)) = 15.0*pow(*(ita+(i)*ny*nz+(j)*nz+(k)),2) *pow((*(ita+(i)*ny*nz+(j)*nz+(k))-1),2) *\
                    epsilon0* *(E_2+(i)*ny*nz+(j)*nz+(k))*(epsilonb-epsilonmatrix);
                }
                if (*(kind+(i)*ny*nz+(j)*nz+(k)) == 1)
                {
                    *(f_elec_grad+(i)*ny*nz+(j)*nz+(k)) = 15.0*pow(*(ita+(i)*ny*nz+(j)*nz+(k)),2) *pow((*(ita+(i)*ny*nz+(j)*nz+(k))-1),2) *\
                    epsilon0* *(E_2+(i)*ny*nz+(j)*nz+(k))*(epsilonb-epsilonpartical);
                }
            }
        }
    }

    

            
    // initialize the Fourier ita matrix 

    fftw_complex *in1, *out1, *in2,*out2,*in3; /* double [2] */
    fftw_plan p,q,pp,qq;
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    

    //int N = nx*ny*nz;
    //fftw_complex *in1, *out1, *in2,*out2,*in3; /* double [2] */
    //fftw_plan p,q,pp,qq;

    //in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    //out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    //in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    //out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    //in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for (int k=0;k < nz; k++)
            {
                in1[(i)*ny*nz+(j)*nz+(k)][0] = *(ita+(i)*ny*nz+(j)*nz+(k));
                in1[(i)*ny*nz+(j)*nz+(k)][1] = 0;
            }   
        }
    }

    /* forward Fourier transform, save the result in 'out' */

    //p = fftw_plan_dft_3d(nx,ny,nz, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
    
    p = fftw_plan_dft_3d(nx,ny,nz, in1, out1, FFTW_FORWARD, FFTW_MEASURE);
    


    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(in1);

    
    
    
    for (int i = 0; i < N; i++) {
        out1[i][0] /= N;
        out1[i][1] /= N;
    }
    //complex double ita_f_rem[nx][ny][nz];

    complex double *ita_f_rem= malloc(nx * ny * nz * sizeof(complex double));




    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for (int k=0;k < nz; k++)
            {
                *(ita_f_rem+(i)*ny*nz+(j)*nz+(k)) = out1[(i)*ny*nz+(j)*nz+(k)][0]+out1[(i)*ny*nz+(j)*nz+(k)][1]*_Complex_I;
            }            
        }
    }


    




    // now the out1 will be the ita_grad_F
    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for (int k=0; k<nz; k++)
            {
                complex double ita_f = out1[(i)*ny*nz+(j)*nz+(k)][0] + out1[(i)*ny*nz+(j)*nz+(k)][1]*_Complex_I, ita_grad_f;

                ita_grad_f = r*ita_f*(cexp(-_Complex_I*2.0*pi*i/nx*1.0)+cexp(-_Complex_I*2.0*pi*j/ny*1.0)+\
                cexp(-_Complex_I*2.0*pi*i/nx*(-1.0))+cexp(-_Complex_I*2.0*pi*j/ny*(-1.0))+\
                cexp(-_Complex_I*2.0*pi*k/nz*(-1.0))+cexp(-_Complex_I*2.0*pi*k/nz*(1.0))-6.0)/pow(h,2);

                out1[(i)*ny*nz+(j)*nz+(k)][0] = creal(ita_grad_f);
                out1[(i)*ny*nz+(j)*nz+(k)][1] = cimag(ita_grad_f);
            }
            
        }
    }

    // the real part of in2 is the contribution of the phase field
    //q = fftw_plan_dft_3d(nx, ny, nz, out1, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    
    q = fftw_plan_dft_3d(nx, ny, nz, out1, in2, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(q);
    fftw_free(out1);

    
    
    
    
    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for ( int k=0;k < nz; k++) 
            {
                double contri = (-in2[(i)*ny*nz+(j)*nz+(k)][0]-*(f_elec_grad+(i)*ny*nz+(j)*nz+(k))+*(f_sep_grad+(i)*ny*nz+(j)*nz+(k)))*\
                (-*(H+(i)*ny*nz+(j)*nz+(k)));
                in2[(i)*ny*nz+(j)*nz+(k)][0] = contri;
                in2[(i)*ny*nz+(j)*nz+(k)][1] = 0;
            }
        }
    }


    
    


    // now the out2 is the whole contribution
    //pp = fftw_plan_dft_3d(nx,ny,nz, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

    pp = fftw_plan_dft_3d(nx,ny,nz, in2, out2, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(pp);
    
    
    
    for (int i = 0; i < N; i++) {
        out2[i][0] /= N;
        out2[i][1] /= N;
    }

    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for (int k=0; k<nz; k++)
            {
                complex double contri_f = out2[(i)*ny*nz+(j)*nz+(k)][0] +\
                out2[(i)*ny*nz+(j)*nz+(k)][1]*_Complex_I, contri_all_f,coeff;

                coeff = (cexp(-_Complex_I*2.0*pi*i/nx*1.0)+cexp(-_Complex_I*2.0*pi*j/ny*1.0)+\
                cexp(-_Complex_I*2.0*pi*i/nx*(-1.0))+cexp(-_Complex_I*2.0*pi*j/ny*(-1.0))+\
                cexp(-_Complex_I*2.0*pi*k/nz*(-1.0))+cexp(-_Complex_I*2.0*pi*k/nz*(1.0))-6.0)/pow(h,2);

                contri_all_f = contri_f*l*deltat/(1.0+0.5*l*r*deltat*coeff)+*(ita_f_rem+(i)*ny*nz+(j)*nz+(k));
                out2[(i)*ny*nz+(j)*nz+(k)][0] = creal(contri_all_f);
                out2[(i)*ny*nz+(j)*nz+(k)][1] = cimag(contri_all_f);
            }
        }
    }

    qq = fftw_plan_dft_3d(nx, ny, nz, out2, in3, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(qq);
    
    
    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for ( int k=0;k < nz; k++) 
            {
                *(ita+(i)*ny*nz+(j)*nz+(k))=in3[(i)*ny*nz+(j)*nz+(k)][0];
                if(*(ita+(i)*ny*nz+(j)*nz+(k))>1.0)
                {
                    *(ita+(i)*ny*nz+(j)*nz+(k))=1;
                }
                if(*(ita+(i)*ny*nz+(j)*nz+(k))<1e-6 && *(H+(i)*ny*nz+(j)*nz+(k))==0)
                {
                    *(ita+(i)*ny*nz+(j)*nz+(k))=0;
                }
                if (isnan(*(ita+(i)*ny*nz+(j)*nz+(k))) == 1)
                {
                    *debug =0;
                    break;
                }

            }
        }
    }

    fftw_destroy_plan(q);
    fftw_destroy_plan(pp);
    fftw_destroy_plan(qq);
    fftw_free(in2); 
    fftw_free(out2);
    fftw_free(in3); 
      
    //update the epsilon distribution
    for ( int i=0;i < nx; i++) 
    {
        for ( int j=0;j < ny; j++) 
        {
            for ( int k=0;k < nz; k++) 
            {
                if(*(kind+(i)*ny*nz+(j)*nz+(k))==0)
                {
                    double ratio = pow(*(ita+(i)*ny*nz+(j)*nz+(k)),3)*\
                    (10.0-15.0**(ita+(i)*ny*nz+(j)*nz+(k))+6.0*pow(*(ita+(i)*ny*nz+(j)*nz+(k)),2));
                    *(epsilon+(i)*ny*nz+(j)*nz+(k)) = epsilonb*ratio + epsilonmatrix*(1-ratio);
                }
                if(*(kind+(i)*ny*nz+(j)*nz+(k))==1)
                {
                    double ratio = pow(*(ita+(i)*ny*nz+(j)*nz+(k)),3)*\
                    (10.0-15.0**(ita+(i)*ny*nz+(j)*nz+(k))+6.0*pow(*(ita+(i)*ny*nz+(j)*nz+(k)),2));
                    *(epsilon+(i)*ny*nz+(j)*nz+(k)) = epsilonb*ratio + epsilonpartical*(1-ratio);
                }
            }
        }
    }



    free(f_sep_grad);
    free(f_elec_grad);
    free(ita_f_rem);
}



 
int main()
{
    


    double h = d_mesh;
    double U_up = 30e3*1e3*nx*h;  // unit: V



    double deltat = 1e-9;
    double alpha = 1e8;   // origin = 1e8
    double r = 1e-10;     // origin = 1e-10
    double l = 1;

    double epsilonb = 1e4, epsilonmatrix = 20, epsilonpartical = 400;
    double norm_cri_matrix = 370*1e3*1e3, norm_cri_particle =50*1e3*1e3;

    int debug = 1;


    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                *(phi_main+i*ny*nz+j*nz+k) = U_up*(nx-1-i)/(nx-1);
            }
        }
    }


    // set initial breakdown region and initial condition
    
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {

                *(kind_main+i*ny*nz+j*nz+k) = 0;


                if (i > nx/2.0-10 && i < nx/2.0+10 && j > ny/2.0-10 && j < ny/2.0+10 && \
                k > nz/2.0-10 && k < nz/2.0+10)
                {
                    *(ita_main+i*ny*nz+j*nz+k) = 1.0;
                }
                else
                {
                    *(ita_main+i*ny*nz+j*nz+k) = 0.0;
                }


                if (*(kind_main+i*ny*nz+j*nz+k) == 0)
                {
                    *(norm_cri_main+i*ny*nz+j*nz+k) = norm_cri_matrix;
                }
                if (*(kind_main+i*ny*nz+j*nz+k) == 1)
                {
                    *(norm_cri_main+i*ny*nz+j*nz+k) = norm_cri_particle;
                }
            }
            
        }
    }

    // update the epsilon distribution
    

    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                if (*(kind_main+i*ny*nz+j*nz+k) == 0)
                {
                    double ratio = pow(*(ita_main+i*ny*nz+j*nz+k),3)*(10-15**(ita_main+i*ny*nz+j*nz+k)+\
                    6*pow(*(ita_main+i*ny*nz+j*nz+k),2));                
                    *(epsilon_main+i*ny*nz+j*nz+k) = epsilonb*ratio + epsilonmatrix*(1-ratio);
                }
                if (*(kind_main+i*ny*nz+j*nz+k) == 1)
                {
                    double ratio = pow(*(ita_main+i*ny*nz+j*nz+k),3)*(10-15**(ita_main+i*ny*nz+j*nz+k)+\
                    6*pow(*(ita_main+i*ny*nz+j*nz+k),2));                
                    *(epsilon_main+i*ny*nz+j*nz+k) = epsilonb*ratio + epsilonpartical*(1-ratio);
                }
            }
        }
    }

    /*********************************************************************************************************/

    // loop the phase field

    clock_t start_all,end_all; 
    start_all = clock();  

    for(int step=0; step<301; step++)
    {
        

        clock_t start,end;  
        start = clock();  
        double U_up_update = U_up+1e3*1e3*nx*h*(50*(step+400)/100);
        for(int i=0; i<nx; i++)
        {
            for(int j=0; j<ny; j++)
            {
                for(int k=0; k<nz; k++)
                {
                    *(phi_main+i*ny*nz+j*nz+k) = U_up_update*(nx-1-i)/(nx-1);
                }
            }
        }
        electrical_calculate(phi_main, epsilon_main, step,&debug);
        
        
        generate_electrical_field(phi_main, E_2_main);
        

        phase_cal(E_2_main, epsilon_main, norm_cri_main, ita_main, kind_main,\
        deltat, epsilonb, epsilonmatrix, epsilonpartical, h, alpha, r, l,H,&debug);

        if(debug==0){
            printf("bad environment!");
            break;
        }

        end = clock();  
        printf("the step %d calculation time is %.2f s.\n",step,(double)(end-start)/1000.0/1000.0);  
        

        if (step%50==0)
        {
            FILE * pFile1;
            char filename1[20];
            sprintf(filename1,"%d_E_2",step);
            pFile1 = fopen (filename1,"w");
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                        fprintf (pFile1, "%f\n", *(E_2_main+i*ny*nz+j*nz+k));
                    }
                }
            }
            fclose (pFile1);


            FILE * pFile2;
            char filename2[20];
            sprintf(filename2,"%d_ita",step);
            pFile2 = fopen (filename2,"w");
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                    fprintf (pFile2, "%f\n", *(ita_main+i*ny*nz+j*nz+k));
                    }
                }
            }
            fclose (pFile2);

            FILE * pFile3;
            char filename3[20];
            sprintf(filename3,"%d_epsilon",step);
            pFile3 = fopen (filename3,"w");
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                    fprintf (pFile3, "%f\n", *(epsilon_main+i*ny*nz+j*nz+k));
                    }
                }
            }
            fclose (pFile3);

            FILE * pFile4;
            char filename4[20];
            sprintf(filename4,"%d_H",step);
            pFile4 = fopen (filename4,"w");
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                    fprintf (pFile3, "%f\n", *(H+i*ny*nz+j*nz+k));
                    }
                }
            }
            fclose (pFile4);
        }


        



    }

    end_all = clock();  
    printf("the total calculation time is %.2f s.\n",(double)(end_all-start_all)/1000.0/1000.0);  



    return 0;
}