#include <iostream>
#include <cmath>


//contains C/C++ subroutines for performing the DG time-stepping for the
//advection equation.  Since the core routine will only be called from 
//a FORTRAN main program, all array data will be arranged in FORTRAN 
//column-major style.

using namespace std;

//-------------------- subroutine for computing RHS ---------------------

void compute_rhs(float * rhs, float * q, float * u, float * v, 
		 float * ksi_x, float * ksi_y, float * eta_x,
		 float * eta_y, float * jac, float * psi,
		 float * psi_ksi, float * psi_eta, int nelem,
		 int npts, int nqs){

  //local variables
  float wq, e_x, e_y, n_x,n_y,h_k,u_k,v_k,q_k,h_e,h_n,dhdx_k,dhdy_k;

  

    //initialize rhs
    for(int l = 0;l<npts;l++){
      for(int ie=0;ie<nelem;ie++){
	rhs[l*nelem+ie]=0.0;
      }//for(int ie...
    }//for(int l...

  //loop over the elements
  for(int ie=0;ie<nelem;ie++){
    //loop over integration points
    for(int k=0;k<nqs;k++){
      wq = jac[k*nelem+ie];
      e_x = ksi_x[k*nelem+ie];
      e_y = ksi_y[k*nelem+ie];
      n_x = eta_x[k*nelem+ie];
      n_y = eta_y[k*nelem+ie];

      //iterpolate at the integration points
      u_k= 0.; v_k = 0.; q_k = 0.;
      for(int i=0;i<npts;i++){
	h_k = psi[k*npts+i];
	u_k += h_k*u[i*nelem+ie];
	v_k+=h_k*v[i*nelem+ie];
	q_k+=h_k*q[i*nelem+ie];
      }//for(int i=0...

      //interpolate at integration points...more...
      for(int i=0;i<npts;i++){
	h_e = psi_ksi[k*npts+i];
	h_n = psi_eta[k*npts+i];
	dhdx_k = h_e*e_x+h_n*n_x;
	dhdy_k = h_e*e_y+h_n*n_y;
	rhs[i*nelem+ie]+=wq*q_k*(dhdx_k*u_k+dhdy_k*v_k);

      }//for(int i=0...

    }//for(int k=0...

  }//for(int ie...
}
//------------------- end compute_rhs routine -----------------------------

//------------------- routine for computing/applying flux -----------------

void compute_flux(float * rhs, float * q, float * u, float * v, 
		  int * psideh, float * nx, float * ny, float * jac_side,
		  float * psi, int nside, int ngl, int nq, int * imapl,
		  int * imapr, int nelem, int npts){

  //local arrays
  float * ql; float * qr; float * ul; float * vl;
  ql = new float[ngl];
  qr = new float[ngl];
  ul = new float[ngl];
  vl = new float[ngl];

  int nqs = nq*nq;
  //initialize the local arrays
  for(int i=0;i<ngl;i++){
    ql[i] = 0.; qr[i]=0.; ul[i]=0.; vl[i]=0.;
  }

  //some local variables
  int iel, ilocl,il,jl,kl,ier,ilocr,ir,jr,kr;
  float wq, nxl,nyl,nxr,nyr,qlq_k, qrq_k,u_k,v_k,ul_k, vl_k, ur_k,vr_k;
  float unl, unr, claml, clamr, clam, fxl, fyl, fxr, fyr, flux_ql;
  float flux_qr, flux_q, diss_q,h_i;

  // cout << "in compute flux: " << endl;
  // cout << "nside = " << nside << endl;
  // cout << "ngl = " << ngl << endl;


  for(int is = 0; is<nside;is++){
    iel = psideh[(3-1)*nside+is];
    if(iel != -6){ //not periodic
      ilocl = psideh[(1-1)*nside+is]; //just being explicit about the indexing
      for(int l=0;l<ngl;l++){
	//get pointers
	il = imapl[l*(4*2)+(1-1)*4+(ilocl-1)];
	jl = imapl[l*(4*2)+(2-1)*4+(ilocl-1)];
	kl = (jl-1)*ngl+il-1; //kl is now zero-based pointer
	//left element
	ql[l] = q[kl*nelem+(iel-1)];
	ul[l] = u[kl*nelem+(iel-1)];
	vl[l] = v[kl*nelem+(iel-1)];
      }//for(int l=0...

      //store right side variables
      ier = psideh[(4-1)*nside + is];
      if(ier != 0){
	ilocr = psideh[(2-1)*nside+is];
	for(int l = 0;l<ngl;l++){
	  ir = imapr[l*(4*2)+(1-1)*4+(ilocr-1)];
	  jr = imapr[l*(4*2)+(2-1)*4+(ilocr-1)];
	  kr = (jr-1)*ngl+ir-1;//kr is now zero-based pointer
	  qr[l] = q[kr*nelem+(ier-1)];
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

	//flux variables
	fxl = qlq_k*ul_k;
	fyl = qlq_k*vl_k;
	fxr = qrq_k*ur_k;
	fyr = qrq_k*vr_k;

	//normal flux component
	flux_ql = nxl*fxl+nyl*fyl;
	flux_qr = nxr*fxr+nyr*fyr;
	flux_q = flux_ql-flux_qr;

	//dissipation term
	diss_q = clam*(qrq_k-qlq_k);

	//construct Rusanov flux
	flux_q = 0.5*(flux_q-diss_q);

	// //loop through side interpolation points
	// cout << "side = " << is+1 << endl;

	for(int i=0;i<ngl;i++){
	  h_i = psi[i*ngl+l];
	  //left side
	  il=imapl[i*(4*2)+(1-1)*4+(ilocl-1)];
	  jl=imapl[i*(4*2)+(2-1)*4+(ilocl-1)];
	  kl = (jl-1)*ngl+il-1;
	  // cout << "ngl = " << i+1 << ", incr = " << 
	  //   -wq*h_i*flux_q << endl;
	  rhs[kl*nelem+iel-1]-=wq*h_i*flux_q;

	  //right side
	  if(ier > 0){
	    ir = imapr[i*(4*2)+(1-1)*4 +(ilocr-1)];
	    jr = imapr[i*(4*2)+(2-1)*4 +(ilocr-1)];
	    kr = (jr-1)*ngl+ir - 1;
	    // cout << "right side, ngl = " << i+1 <<
	    //   ", incr = " << wq*h_i*flux_q << endl;
	    rhs[kr*nelem+ier-1]+=wq*h_i*flux_q;
	  }//if(ier...
	}//for(int i...
      }//for(int l...
    }//if(iel != -6...

  }//loop over sides


  //deallocate heap memory for local arrays
  delete [] ql; delete [] qr; delete [] ul; delete [] vl;

}

//--------------------- end compute_flux routine ---------------------------

//--------------------- Apply inverse mass matrix routine ------------------

void apply_m_inv(float * rhs, float * Mmatrix_inv, int nelem, int npts){

  float * Mtemp = new float[npts*npts];
  float * rtemp = new float[npts];
  int m_dof, mtemp_dof;
  float rhs_temp;

  
  for(int e = 0;e<nelem;e++){

    //get the inverse mass matrix for the element e
    //M[e,i,j] --> Mtemp(i,j)
    //get Mtemp for this element
    for(int j = 0;j<npts;j++){
      for(int i=0;i<npts;i++){
	m_dof = j*nelem*npts+i*nelem+e;
	mtemp_dof = j*npts + i;
	Mtemp[mtemp_dof] = Mmatrix_inv[m_dof];
      }//for(int i=0...
    }//for(int j=0...

    //rhs(e,:) = Mtemp*rhs(e,:)
    //get rtemp
    for(int i=0;i<npts;i++){
      rtemp[i] = rhs[i*nelem+e];
    }//for(int i=0...

    //do the multiplication
    for(int i=0;i<npts;i++){
      rhs_temp = 0.0;
      for(int j=0;j<npts;j++){
	rhs_temp+=Mtemp[j*npts+i]*rtemp[j];
      }
      rhs[i*nelem+e] = rhs_temp;//assign to rhs
    }//for(int i=0...
  }//for(int e=0...


  delete [] rtemp;

  delete [] Mtemp;
}

//--------------------- end apply inverse mass matrix routine --------------

//--------------------- interstage update of qp ----------------------------

void interstage_update_qp(float * qp, float * q0, float * q1, float * rhs,
			  float a0, float a1, float dtt, int npts,int nelem){

  int dof;
  for(int l=0;l<npts;l++){
    for(int ie=0;ie<nelem;ie++){
      dof = l*nelem+ie;
      qp[dof]=a0*q0[dof]+a1*q1[dof]+dtt*rhs[dof];
    }
  }

}

//---------------------- end interstage update of qp -----------------------

//---------------------- end of time step update of q0 or q1 ----------------
void dof_copy(float * q_a, float * q_b, int npts, int nelem){
  int dof;
  for(int l=0;l<npts;l++){
    for(int ie=0;ie<nelem;ie++){
      dof = l*nelem+ie;
      q_a[dof] = q_b[dof];
    }
  }


}

//----------------------- end of time step update of q0 --------------------

//----------------------- print rhs ----------------------------------------

void dof_print(float * q, int npts, int nelem){

  for(int ie = 0; ie<nelem; ie++){
    cout << endl;
    for(int l=0;l<npts; l++){
      cout << q[l*nelem+ie];
      if(l<(npts-1)){
	cout << ", ";
      }
    }
  }
  cout << endl;
}

//--------------------------------------------------------------------------




extern "C" void dg_advection_2d_ntp_c_(int* ntime_ptr, float* dt_ptr, int* kstages_ptr, 
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

  //local variables
  float dt,a0, a1, beta, dtt;
  const int time_step_rep_freq = 100;
  
  int ntime,kstages,nelem, npts, nqs, ngl, nq, nside;
  ntime = *ntime_ptr;
  kstages = *kstages_ptr;
  nelem = *nelem_ptr;
  npts = *npts_ptr;
  ngl = *ngl_ptr;
  nq = *nq_ptr;
  ngl = *ngl_ptr;
  dt = *dt_ptr;
  nqs = *nqs_ptr;
  nside = *nside_ptr;
 

  //declare local arrays
 
  float * rhs; // nelem x npts
  float * qp; // nelem x npts
  float * q1; // nelem x npts
  // cout << "nelem = " << nelem << ", npts = " << npts << endl;
  // cout << "Allocating memory for local arrays..." << endl;
  // cout << "ntime = " << ntime << endl;
  // cout << "kstages = " << kstages << endl;
  // cout << "nq = " << nq << endl;

 

  //allocate memory for local arrays
 
  rhs = new float[nelem*npts];
  qp = new float[nelem*npts];
  q1 = new float[nelem*npts];

  cout << "Initializing state vector..." << endl;
  //initialize the state vector
  for(int l=0;l<npts;l++){
    for(int ie=0;ie<nelem;ie++){
      rhs[l*nelem+ie] = 0.0;
      q1[l*nelem+ie] = q0[l*nelem+ie];
      qp[l*nelem+ie] = q0[l*nelem+ie];
    }
  }



  cout << "Commencing time-steps " << endl;
  //time-step loop
  for(int itime = 0; itime<ntime;itime++){
    //  cout << "In time-step " << itime << endl;
    // cout << "rhs at start of time-step " << itime+1 << endl;
    // dof_print(rhs,npts,nelem);

    if( (itime+1) % time_step_rep_freq == 0){
      cout << "Commencing time step " << itime + 1 << endl;
    }


    for(int ik = 1; ik<=kstages;ik++){
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
	  a0 = 1.; a1 = 0.; beta = 1.;
	  break;
	case 2:
	  a0 = 3./4.; a1 = 1./4.; beta = 1./4.;
	  break;
	case 3:
	  a0 = 1./3.; a1 = 2./3.; beta = 2./3.;

	}//switch(ik)...

      }//switch (kstages)...
      dtt = dt*beta;

      //cout << "Calling compute_rhs..." << endl;
      //compute RHS
      compute_rhs(rhs,qp,u0,v0,ksi2d_x,ksi2d_y,eta2d_x,eta2d_y,
		  jac2d,psi2d,psi2d_ksi,psi2d_eta,nelem,npts,nqs);

      // // cout << "RHS in timestep " << itime+1 << ", in stage " << ik <<
      // // 	"after calling compute_rhs. ";
      // // dof_print(rhs,npts,nelem);

      // cout << "Calling compute_flux..." << endl;
      //compute/apply flux
      compute_flux(rhs,qp,u0,v0,psideh,nx,ny,jac_side,psi,nside,
		   ngl,nq,imapl,imapr,nelem,npts);


      // cout << "RHS in timestep " << itime+1 << ", in stage " << ik <<
      // 	"after calling compute_flux.";
      // dof_print(rhs,npts,nelem);

      //apply inverse mass matrix
     
      apply_m_inv(rhs,Mmatrix_inv,nelem,npts);
     
      // cout << "RHS in timestep " << itime+1 << ", in stage " << ik <<
      // 	"after applying inverse mass matrix. ";
      // dof_print(rhs,npts,nelem);

      //update inter-stage qp

      interstage_update_qp(qp,q0,q1,rhs,a0,a1,dtt,npts,nelem);

      //update q1
      dof_copy(q1,qp,npts,nelem);

    }//for(int ik...

    //update q0
    dof_copy(q0,qp,npts,nelem);

  }//for(int itime...



  //free memory used for local arrays
 
  delete [] rhs;
  delete [] qp;
  delete [] q1;

}
