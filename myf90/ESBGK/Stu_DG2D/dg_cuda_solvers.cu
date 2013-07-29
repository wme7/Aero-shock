#include <iostream>
#include <cmath>
#include "device_functions.h"
#include "cublas_v2.h"
#include <cuda_runtime.h>  
#define MAX_NGL 5
#define TPB_1 96
#define TPB_2 32


using namespace std;

//------- compute rhs kernel ----------------------------------------------------------



// ----- new version with additional ranks of threads at the nqs level and atomic operations

__global__ void compute_rhs_cuda(float * rhs, float * q, float * u, float * v, 
			    float * ksi_x, float * ksi_y, float * eta_x, 
			    float * eta_y, float * jac, float * psi, 
			    float * psi_ksi, float * psi_eta, int nelem, 
			    int npts, int nqs){

  int l_tid = threadIdx.x;
  int l_block = blockDim.x;
  int tid = threadIdx.x+blockIdx.x*blockDim.x;
  int k = threadIdx.y;

  if(tid<nelem){
    if(k<npts){
      float wq; float e_x; float e_y; float n_x; float n_y; float u_k; float v_k;
      float q_k; float h_k; float h_e; float h_n; float dhdx_k; float dhdy_k;
     
      __shared__ float rhs_s[TPB_1*MAX_NGL*MAX_NGL];
      //cooperatively load rhs for the block into shared memory

      // rhs[nelem*k+tid] = 0.0;
      rhs_s[l_block*k+l_tid] = 0.0;
      __syncthreads();
      
      int dof = k*nelem+tid;
      wq = jac[dof]; e_x = ksi_x[dof]; e_y = ksi_y[dof];
      n_x = eta_x[dof]; n_y = eta_y[dof];
      u_k = 0.0; v_k = 0.0; q_k = 0.0;
      for(int i=0;i<npts;i++){

	h_k = psi[k*npts+i];
	u_k = u_k +h_k*u[i*nelem+tid];
	v_k = v_k+h_k*v[i*nelem+tid];
	q_k = q_k+h_k*q[i*nelem+tid];
      }//if(int i...

      for(int i=0;i<npts;i++){

	h_e = psi_ksi[k*npts+i];
	h_n = psi_eta[k*npts+i];
	dhdx_k=h_e*e_x+h_n*n_x;
	dhdy_k=h_e*e_y+h_n*n_y;

	//atomicAdd(&rhs[i*nelem+tid],wq*q_k*(dhdx_k*u_k+dhdy_k*v_k));
	atomicAdd(&rhs_s[i*l_block+l_tid],wq*q_k*(dhdx_k*u_k+dhdy_k*v_k));
      }//for(int i=...
     
      __syncthreads();
      //cooperatively load rhs_s back to global memory
      rhs[k*nelem+tid]=rhs_s[k*l_block+l_tid];

    }
  }//if(tid<nelem)
}






//-------------- compute flux kernel ---------------------------------------------------
__global__ void compute_flux_cuda(float * rhs, float * q, float * u,
			     float * v, int * psideh, float * nx, 
			     float * ny, float * jac_side, float * psi,
			     int nside, int ngl, int nq, 
			     int * imapl, int * imapr,int nelem){

  //each thread updates a side
  int is = threadIdx.x + blockIdx.x*blockDim.x;
  //int l_is = threadIdx.x;

  if(is<nside){
    
    int iel = psideh[2*nside+is];//element on rhs of side
    //float old;
    if(iel != -6){
     

      //some local variables
      //int nqs = nq*nq;
      //int ilocl,il,jl,kl,ier,ilocr,ir,jr,kr;
      int ilocl, ilocr, ier, i_tmp,j_tmp,k_tmp;
      float wq, nxl,nyl,nxr,nyr,qlq_k, qrq_k,u_k,v_k,ul_k, vl_k, ur_k,vr_k;
      float unl, unr, claml, clamr, clam;// fxl, fyl, fxr, fyr
      //float flux_ql;
      //float flux_qr;
      float flux_q;
      float  diss_q,h_i;

      float ql[MAX_NGL];
      float qr[MAX_NGL];
      float ul[MAX_NGL];
      float vl[MAX_NGL];
     
     

      ilocl = psideh[(1-1)*nside+is]; //just being explicit about the indexing
      for(int l=0;l<ngl;l++){
	//get pointers
	i_tmp = imapl[l*(4*2)+(1-1)*4+(ilocl-1)];
	j_tmp = imapl[l*(4*2)+(2-1)*4+(ilocl-1)];
	k_tmp = (j_tmp-1)*ngl+i_tmp-1; //kl is now zero-based pointer
	//left element
	ql[l] = q[k_tmp*nelem+(iel-1)];
	ul[l] = u[k_tmp*nelem+(iel-1)];
	vl[l] = v[k_tmp*nelem+(iel-1)];
      }//for(int l=0...

      //store right side variables
      ier = psideh[(4-1)*nside + is];
      if(ier != 0){
	ilocr = psideh[(2-1)*nside+is];
	for(int l = 0;l<ngl;l++){
	  i_tmp = imapr[l*(4*2)+(1-1)*4+(ilocr-1)];
	  j_tmp = imapr[l*(4*2)+(2-1)*4+(ilocr-1)];
	  k_tmp = (j_tmp-1)*ngl+i_tmp-1;//kr is now zero-based pointer
	  qr[l] = q[k_tmp*nelem+(ier-1)];
	}//for(int l=0...
      }//if(ier !=0...

      //do gauss-lobatto integration
      for(int l=0;l<nq;l++){
	wq = jac_side[l*nside+is];
	//store normal vectors
	nxl = nx[l*nside+is];
	nyl = ny[l*nside+is];
	nxr = -nxl; 
	nyr = -nyl;

	//interpolate onto quad points
	qlq_k = 0.; qrq_k = 0.; u_k = 0.; v_k = 0.;
	for(int i=0;i<ngl;i++){
	  qlq_k = qlq_k+psi[l*ngl+i]*ql[i];
	  qrq_k = qrq_k+psi[l*ngl+i]*qr[i];
	  u_k = u_k + psi[l*ngl+i]*ul[i];
	  v_k = v_k+psi[l*ngl+i]*vl[i];
	}//for(int i=0...
	ul_k = u_k; vl_k = v_k;
	ur_k = u_k; vr_k = v_k;

	//compute Rusanov flux constant
	unl = nxl*ul_k + nyl*vl_k;
	unr = nxl*ur_k+nyl*vr_k;
	claml=fabs(unl);
	clamr=fabs(unr);
	//clam = max(claml,clamr);
	if(claml > clamr){
	  clam = claml;
	}else{
	  clam = clamr;
	}

	// //flux variables
	// fxl = qlq_k*ul_k;
	// fyl = qlq_k*vl_k;
	// fxr = qrq_k*ur_k;
	// fyr = qrq_k*vr_k;

	// //normal flux component
	// flux_ql = nxl*fxl+nyl*fyl;
	// flux_qr = nxr*fxr+nyr*fyr;
	// flux_q = flux_ql-flux_qr;


	flux_q = nxl*qlq_k*ul_k+nyl*qlq_k*vl_k-(nxr*qrq_k*ur_k+nyr*qrq_k*vr_k);
	//dissipation term
	diss_q = clam*(qrq_k-qlq_k);

	//construct Rusanov flux
	flux_q = 0.5*(flux_q-diss_q);

	// //loop through side interpolation points
	// cout << "side = " << is+1 << endl;

	for(int i=0;i<ngl;i++){
	  h_i = psi[i*ngl+l];
	  //left side
	  i_tmp=imapl[i*(4*2)+(1-1)*4+(ilocl-1)];
	  j_tmp=imapl[i*(4*2)+(2-1)*4+(ilocl-1)];
	  k_tmp = (j_tmp-1)*ngl+i_tmp-1;
	  // cout << "ngl = " << i+1 << ", incr = " << 
	  //   -wq*h_i*flux_q << endl;
	  //rhs[kl*nelem+iel-1]-=wq*h_i*flux_q;
	  atomicAdd(&rhs[k_tmp*nelem+iel-1],-wq*h_i*flux_q);

	  //right side
	  if(ier > 0){
	    i_tmp = imapr[i*(4*2)+(1-1)*4 +(ilocr-1)];
	    j_tmp = imapr[i*(4*2)+(2-1)*4 +(ilocr-1)];
	    k_tmp = (j_tmp-1)*ngl+i_tmp - 1;
	    // cout << "right side, ngl = " << i+1 <<
	    //   ", incr = " << wq*h_i*flux_q << endl;
	    //rhs[kr*nelem+ier-1]+=wq*h_i*flux_q;
	    atomicAdd(&rhs[k_tmp*nelem+ier-1],wq*h_i*flux_q);
	  }//if(ier...
	}//for(int i...
      }//for(int l...



    }//skip sides with parallel BCs
  }//if(is<nside...
}


//-------------------------------------------------------------------------------------

//---- apply inverse mass-matrix kernel -----------------------------------------------

__global__ void apply_invM_cuda(float * rhs, float * invM, int nelem, int npts){

  int l_dof= threadIdx.x; int blk_sz = blockDim.x;
  int ie = threadIdx.x + blockIdx.x*blockDim.x;
  int l = threadIdx.y; //which of npts rows this thread processes

  if(ie<nelem){
    if(l<npts){ //all threads "should" satisfy this...not sure though
      //all threads collaboratively load rhs into shared memory
      __shared__ float rhs_s[TPB_2*MAX_NGL*MAX_NGL];

      rhs_s[blk_sz*l+l_dof] = rhs[l*nelem+ie];
      __syncthreads();


      float tmpRHS = 0.0;
      int invM_dof;
      int rhs_dof;
      //for(int m=0;m<npts;m++){
	//invM_dof = ie*npts*npts+m*npts+l;
	invM_dof = ie*npts*npts+l*npts+l;
	//rhs_dof = m*nelem+ie;
	//rhs_dof = m*blk_sz+l_dof;
	rhs_dof=l*blk_sz+l_dof;
	tmpRHS+=invM[invM_dof]*rhs_s[rhs_dof];
	//}
      __syncthreads();//make sure all l-threads are done using
      //this row of the rhs
      rhs[l*nelem+ie] = tmpRHS;

    }//doing on row of the invM against one row of rhs...
  }//only do the right number of elements
}

// ------------------------------------------------------------------------------------
__global__ void interstage_update(float * q0, float * q1, float *qp, float *rhs,
				  float a0, float a1, float dtt,int nelem, int npts){
  int ie = threadIdx.x+blockIdx.x*blockDim.x;
  int l = threadIdx.y;
  if(ie<nelem){
    if(l<npts){
      int dof = l*nelem+ie;
      qp[dof]=a0*q0[dof]+a1*q1[dof]+dtt*rhs[dof];
    }
  }

}

void print_q(float * q_d, float * q, int nelem, int npts){
  //copy q_d to q
  cudaMemcpy(q,q_d,nelem*npts*sizeof(float),cudaMemcpyDeviceToHost);

  //now print out q
  for(int e = 0;e<nelem;e++){
    cout << endl;
    for(int l=0;l<npts;l++){
      cout << q[l*nelem+e];
      if(l<(npts-1)){
	cout << ", ";
      }
    }
  }
  cout << endl;


}


//-------- main driver function -------------------------------------------------------
extern "C" void dg_advection_2d_ntp_cuda_(int* ntime_ptr, float * dt_ptr, int* kstages_ptr, 
				  float * q0, float * u0, float * v0,
				  float * psi, float * ksi2d_x,
				  float * ksi2d_y, float * eta2d_x,
				  float * eta2d_y, float * jac2d,
				  float * psi2d, float * psi2d_ksi,
				  float * psi2d_eta, int * intma2d,
				  int * psideh, float * jac_side,
				  int * imapl, int * imapr, float * nx,
				  float * ny, int* nelem_ptr, int* npts_ptr,
				  int* nqs_ptr, int* ngl_ptr, int* nq_ptr, 
				  float * Mmatrix_inv, int* nside_ptr){

  int ntime, kstages, nelem, npts, nqs, ngl, nq, nside;
  ntime = *ntime_ptr; kstages = *kstages_ptr; nelem = *nelem_ptr; npts = *npts_ptr;
  nqs = *nqs_ptr; ngl = *ngl_ptr; nq = *nq_ptr; nside = *nside_ptr;
 
  float dt, a0, a1,beta,dtt;
  const int time_step_rep_freq = 100;
  dt = *dt_ptr;
  float * rhs = new float[npts*nelem];
  float * invm = new float[npts*npts*nelem];
  float * q_print = new float[npts*nelem];


  for(int e=0;e<nelem;e++){
    for(int i=0;i<npts;i++){
      for(int j=0;j<npts;j++){
	invm[e*npts*npts+i*npts+j]=Mmatrix_inv[i*nelem*npts+j*nelem+e];
      }
    }
  }

  //declare GPU variables
  float * rhs_d; float * qp_d; float * q1_d; float * q0_d;
  float * u0_d; float * v0_d;
  float * ksi2d_x_d; float * ksi2d_y_d; float * eta2d_x_d;
  float * eta2d_y_d; float * jac2d_d;
  float * psi2d_d; float * psi2d_ksi_d; float * psi2d_eta_d;
  int * psideh_d;
  float * nx_d; float * ny_d; float * jac_side_d;
  float * psi_d;
  int * imapl_d; int * imapr_d;
  float * invm_d;

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);
  if(stat!=CUBLAS_STATUS_SUCCESS){
    cout << "CUBLAS initialization failure!" << endl;
    return;
  }


  //allocate space on GPU
  cudaMalloc((void**)&rhs_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&qp_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&q1_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&q0_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&u0_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&v0_d,(nelem*npts)*sizeof(float));
  cudaMalloc((void**)&ksi2d_x_d,(nelem*nqs)*sizeof(float));
  cudaMalloc((void**)&ksi2d_y_d,(nelem*nqs)*sizeof(float));
  cudaMalloc((void**)&eta2d_x_d,(nelem*nqs)*sizeof(float));
  cudaMalloc((void**)&eta2d_y_d,(nelem*nqs)*sizeof(float));
  cudaMalloc((void**)&jac2d_d,(nelem*nqs)*sizeof(float));
  cudaMalloc((void**)&psi2d_d,(npts*nqs)*sizeof(float));
  cudaMalloc((void**)&psi2d_ksi_d,(npts*nqs)*sizeof(float));
  cudaMalloc((void**)&psi2d_eta_d,(npts*nqs)*sizeof(float));
  cudaMalloc((void**)&psideh_d,(nside*4)*sizeof(int));
  cudaMalloc((void**)&nx_d,(nside*nq)*sizeof(float));
  cudaMalloc((void**)&ny_d,(nside*nq)*sizeof(float));
  cudaMalloc((void**)&jac_side_d,(nside*nq)*sizeof(float));
  cudaMalloc((void**)&psi_d,(ngl*nq)*sizeof(float));
  cudaMalloc((void**)&imapl_d,(4*2*ngl)*sizeof(int));
  cudaMalloc((void**)&imapr_d,(4*2*ngl)*sizeof(int));
  cudaMalloc((void**)&invm_d,(nelem*npts*npts)*sizeof(float));

  //transfer data to the GPU

  cudaMemcpy(rhs_d,rhs,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(qp_d,q0,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(q1_d,q0,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(q0_d,q0,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(u0_d,u0,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(v0_d,v0,nelem*npts*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(ksi2d_x_d,ksi2d_x,nelem*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(ksi2d_y_d,ksi2d_y,nelem*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(eta2d_x_d,eta2d_x,nelem*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(eta2d_y_d,eta2d_y,nelem*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(jac2d_d,jac2d,nelem*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(psi2d_d,psi2d,npts*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(psi2d_ksi_d,psi2d_ksi,npts*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(psi2d_eta_d,psi2d_eta,npts*nqs*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(psideh_d,psideh,nside*4*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(nx_d,nx,nside*nq*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(ny_d,ny,nside*nq*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(jac_side_d,jac_side,nside*nq*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(psi_d,psi,ngl*nq*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(imapl_d,imapl,4*2*ngl*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(imapr_d,imapr,4*2*ngl*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(invm_d,invm,nelem*npts*npts*sizeof(float),cudaMemcpyHostToDevice);

  //establish thread launch configuration
  //dim3 dimBlock1(TPB_1,npts,1);
  dim3 dimBlock1(TPB_1,1,1);
  dim3 dimBlock_rhs(TPB_1,npts,1);
  dim3 dimGrid1((nelem+TPB_1-1)/TPB_1,1,1);//for rhs, inv mass matrix
  dim3 dimGrid2((nside+TPB_1-1)/TPB_1,1,1);//for computing/applying flux
  dim3 dimBlock2(TPB_2,npts,1);//for inverting mass matrix
  dim3 dimGrid3((nelem+TPB_2-1)/TPB_2,1,1);//for inverting mass matrix

  dim3 dimBlock4(TPB_2,npts,1);
  dim3 dimGrid4((nelem+TPB_2-1)/TPB_2,1,1);//for inter-stage updates
 
  cudaFuncSetCacheConfig(compute_flux_cuda,cudaFuncCachePreferL1);
  //commence time-stepping

  // cout << "About to commence time stepping, q0_d = " << endl;
  // print_q(q0_d,q_print,nelem,npts);

  for(int itime = 0; itime<ntime;itime++){

    //outpu progress
    if((itime+1)%time_step_rep_freq == 0){
      cout << "Commencing time step " << itime+1 << endl;
    }

    //Explicit RK stages
    for(int ik=1;ik<=kstages;ik++){
      switch (kstages) {
      case 2:
	switch (ik) {
	case 1:
	  a0 = 1.0; a1 = 0.0; beta = 1.0;
	  break;
	case 2:
	  a0 = 0.5; a1 = 0.5; beta = 0.5;
	}//case (2) switch (ik)
	break;
      case 3:
	switch (ik) {
	case 1:
	  a0 = 1.0; a1 = 0.0; beta = 1.0;
	  break;
	case 2:
	  a0 = 3.0/4.0; a1 = 1.0/4.0; beta = 1.0/4.0;
	  break;
	case 3:
	  a0 = 1.0/3.0; a1 = 2.0/3.0; beta = 2.0/3.0;
	}//switch(ik)
      }//switch (kstages)

      dtt = dt*beta;
      //compute rhs
      compute_rhs_cuda<<<dimGrid1,dimBlock_rhs>>>(rhs_d,qp_d,u0_d,v0_d,ksi2d_x_d,
					  ksi2d_y_d,eta2d_x_d,eta2d_y_d,
					  jac2d_d,psi2d_d,psi2d_ksi_d,
					  psi2d_eta_d,nelem,npts,nqs);

      // cout << "After compute_rhs, rhs_d = "<< endl;
      // print_q(rhs_d,q_print,nelem,npts);
      //compute/apply flux


      compute_flux_cuda<<<dimGrid2,dimBlock1>>>(rhs_d,qp_d,u0_d,v0_d,psideh_d,nx_d,ny_d,
					   jac_side_d,psi_d,nside,ngl,nq,imapl_d,
					   imapr_d,nelem);

      // cout << "After compute_flux, rhs_d = " << endl;
      // print_q(rhs_d,q_print,nelem,npts);


      //apply inverse mass matrix (possibly re-structure inverse mass matrix to allow using cublas/cusparse)
      apply_invM_cuda<<<dimGrid3,dimBlock2>>>(rhs_d,invm_d,nelem,npts);



      // cout << "After applying inverse mass matrix, rhs_d = " << endl;
      // print_q(rhs_d,q_print,nelem,npts);

      // cout << "Before inter-stage update, q0_d = " << endl;
      // print_q(q0_d,q_print,nelem,npts);

      // cout << "Before inter-stage update, q1_d = " << endl;
      // print_q(q1_d,q_print,nelem,npts);

      // cout << "Before inter-stage update, qp_d = " << endl;
      // print_q(qp_d,q_print,nelem,npts);

      // cout << "dtt = " << dtt;
      // cout << "a0 = " << a0;
      // cout << "a1 = " << a1;

      //perform inter-stage updates
      //qp_d = a0*q0_d + a1*q1_d + dtt*rhs_d
      stat=cublasScopy(handle,nelem*npts,rhs_d,1,qp_d,1);//qp_d=rhs_d
      stat=cublasSscal(handle,nelem*npts,&dtt,qp_d,1);//qp_d=dtt*qp_d
      stat=cublasSaxpy(handle,nelem*npts,&a1,q1_d,1,qp_d,1); //qp_d = qp_d+a1*q1_d
      stat=cublasSaxpy(handle,nelem*npts,&a0,q0_d,1,qp_d,1);//qp_d = qp_d+a0*q0_d

      // interstage_update<<<dimGrid4,dimBlock4>>>(q0_d,q1_d,qp_d,rhs_d,a0,a1,dtt,nelem,npts);
      

      // cout << "After inter-stage update, qp_d = " << endl;
      // print_q(qp_d,q_print,nelem,npts);

      //q1_d = qp_d
      stat=cublasScopy(handle,nelem*npts,qp_d,1,q1_d,1);

    }//for(int ik...

    //q0_d = qp_d
    stat=cublasScopy(handle,nelem*npts,qp_d,1,q0_d,1);

    // cout << "After all stages, q0_d = " << endl;
    // print_q(q0_d,q_print,nelem,npts);
     

  }
  //end time-stepping section, transfer needed output data back from GPU
  cudaMemcpy(q0,q0_d,nelem*npts*sizeof(float),cudaMemcpyDeviceToHost);

  //de-allocate memory from GPU

  cublasDestroy(handle);

  cudaFree(rhs_d);
  cudaFree(qp_d);
  cudaFree(q1_d);
  cudaFree(q0_d);
  cudaFree(u0_d);
  cudaFree(v0_d);
  cudaFree(ksi2d_x_d);
  cudaFree(ksi2d_y_d);
  cudaFree(eta2d_x_d);
  cudaFree(eta2d_y_d);
  cudaFree(jac2d_d);
  cudaFree(psi2d_d);
  cudaFree(psi2d_ksi_d);
  cudaFree(psi2d_eta_d);
  cudaFree(psideh_d);
  cudaFree(nx_d);
  cudaFree(ny_d);
  cudaFree(jac_side_d);
  cudaFree(imapl_d);
  cudaFree(imapr_d);
  cudaFree(invm_d);

  //de-allocate any local memory
  delete [] rhs;
  delete [] invm;
  delete [] q_print;

}
