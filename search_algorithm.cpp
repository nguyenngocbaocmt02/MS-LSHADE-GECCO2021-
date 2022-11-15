

#include"de.h"
#include"dscg.h"

void searchAlgorithm::initializeParameters() {
  function_number = g_function_number;
  problem_size = g_problem_size;
  max_num_evaluations = g_max_num_evaluations;
  pop_size = g_pop_size;
  initializeFitnessFunctionParameters();
}

void searchAlgorithm::evaluatePopulation(const vector<Individual> &pop, vector<Fitness> &fitness) {
  for (int i = 0; i < pop_size; i++) {
       gl_func(pop[i],  &fitness[i], problem_size, 1, function_number);
  }
}
void searchAlgorithm::evaluateIndividual(const Individual individual, Fitness &fitness) {
   gl_func(individual,  &fitness, problem_size, 1, function_number);
}

void searchAlgorithm::initializeFitnessFunctionParameters() {
  //epsilon is an acceptable error value.
  epsilon = pow(10.0, -8);
  max_region = 100.0;
  min_region = -100.0; 
  optimum=g_optimize_fitness;
}

//set best solution (bsf_solution) and its fitness value (bsf_fitness) in the initial population
void searchAlgorithm::setBestSolution(const vector<Individual> &pop, const vector<Fitness> &fitness, Individual &bsf_solution, Fitness &bsf_fitness) {
  int current_best_individual = 0;

  for (int i = 1; i < pop_size; i++) {
    if (fitness[current_best_individual] > fitness[i]) {
      current_best_individual = i;
    }
  }

  bsf_fitness = fitness[current_best_individual];
  for (int i = 0; i < problem_size; i++) {
    bsf_solution[i] = pop[current_best_individual][i];
  }
}

// make new individual randomly
Individual searchAlgorithm::makeNewIndividual() {
  Individual individual = (variable*)malloc(sizeof(variable) * problem_size);

  for (int i = 0; i < problem_size; i++) {
    individual[i] = ((max_region - min_region) * randDouble()) + min_region;
  }

  return individual;
}

/*
  For each dimension j, if the mutant vector element v_j is outside the boundaries [x_min , x_max], we applied this bound handling method
  If you'd like to know that precisely, please read:
  J. Zhang and A. C. Sanderson, "JADE: Adaptive differential evolution with optional external archive,"
  IEEE Tran. Evol. Comput., vol. 13, no. 5, pp. 945–958, 2009.
 */
void searchAlgorithm::modifySolutionWithParentMedium(Individual child, Individual parent) {
  int l_problem_size = problem_size;
  variable l_min_region = min_region;
  variable l_max_region = max_region;

  for (int j = 0; j < l_problem_size; j++) {
    if (child[j] < l_min_region) {
      child[j]= (l_min_region + parent[j]) / 2.0;
    }
    else if (child[j] > l_max_region) {
      child[j]= (l_max_region + parent[j]) / 2.0;
    }
  }
}


int searchAlgorithm::line_search_inter(Individual x00, Fitness &x00_fitness, Individual result, Fitness &result_fitness,double s0,vector<variable> v) {
	//init
	double anpha=2,beta=0.5,s=s0;
	int eva=0;
	//x
	Individual x=new double[problem_size];
	Fitness x_fitness=x00_fitness;
	for(int i=0;i<problem_size;i++) {
		x[i]=x00[i];
	}
	Individual x0=new double[problem_size];
	Fitness x0_fitness=x00_fitness;
	for(int i=0;i<problem_size;i++) {
		x0[i]=x00[i];
	}
	//
	int ks=1;
	//Step1
	//x=x0+sv
	for(int j=0;j<problem_size;j++) {
		x[j]=x0[j]+s*v[j];
		if(x[j]>max_region) x[j]=max_region;
		if(x[j]<min_region) x[j]=min_region;
	}
	//evaluate fitness of x
	evaluateIndividual(x,x_fitness);
	eva+=1;
	//compare to current x
	if(x_fitness>x0_fitness) {
		//step 2
		ks=0;
		//x=x-2sv
		for(int j=0;j<problem_size;j++) {
		    x[j]=x[j]-2*s*v[j];
	    	if(x[j]>max_region) x[j]=max_region;
	    	if(x[j]<min_region) x[j]=min_region;
	    }
		s=-s;
		evaluateIndividual(x,x_fitness);
		eva+=1;
		if(x_fitness<=x0_fitness) ks=1;
	}
	if(ks==1) {
		while(1) {
		if(eva>=50) break;
		s=anpha*s;
		//x0=x
		for(int i=0;i<problem_size;i++) x0[i]= x[i];
		x0_fitness=x_fitness;
		//x= x0+sv
		for(int j=0;j<problem_size;j++) {
		    x[j]=x0[j]+s*v[j];
	    	if(x[j]>max_region) x[j]=max_region;
	    	if(x[j]<min_region) x[j]=min_region;
	    }
	    evaluateIndividual(x,x_fitness);
	    eva+=1;
	    if (x_fitness>x0_fitness) {
	    	s=s/anpha;
	    	//x=x0+sv
	    	for(int j=0;j<problem_size;j++) x[j]=x0[j]+s*v[j];
	    	evaluateIndividual(x,x_fitness);
			break;
		}
	    }
	}
	
	vector<Individual> interpolytion;
	vector<Fitness> interpolytion_fitness;
	Individual x1= new double[problem_size],x2= new double[problem_size],x3=new double[problem_size],x4= new double[problem_size];
	Fitness x1_fitness,x2_fitness,x3_fitness,x4_fitness;
	for(int j=0;j<problem_size;j++) {
	    x1[j]=x0[j]-s*v[j];
	    x2[j]=x0[j];
	    x3[j]=x0[j]+2*s*v[j];
	    x4[j]=x[0]+s*v[j];
	    if(x1[j]>max_region) x1[j]=max_region;
	    if(x1[j]<min_region) x1[j]=min_region;
	    if(x2[j]>max_region) x2[j]=max_region;
	    if(x2[j]<min_region) x2[j]=min_region;
	    if(x3[j]>max_region) x3[j]=max_region;
	    if(x3[j]<min_region) x3[j]=min_region;
	    if(x4[j]>max_region) x4[j]=max_region;
	    if(x4[j]<min_region) x4[j]=min_region;
	}
	evaluateIndividual(x1,x1_fitness);
	evaluateIndividual(x2,x2_fitness);
	evaluateIndividual(x3,x3_fitness);
	evaluateIndividual(x4,x4_fitness);
	eva+=4;
	interpolytion.push_back(x1);interpolytion.push_back(x2);interpolytion.push_back(x4);interpolytion.push_back(x3);
	interpolytion_fitness.push_back(x1_fitness);interpolytion_fitness.push_back(x2_fitness);interpolytion_fitness.push_back(x4_fitness);interpolytion_fitness.push_back(x3_fitness);
	int reject=0;
	Fitness reject_fitness=interpolytion_fitness[0];
	for(int j=1;j<interpolytion.size();j++) {
		if(interpolytion_fitness[j]>reject_fitness) {
			reject_fitness=interpolytion_fitness[j];
			reject=j;
		}
	}
	interpolytion.erase(interpolytion.begin() + reject);
	interpolytion_fitness.erase(interpolytion_fitness.begin() + reject);
	Individual xt= new double[problem_size];
	Fitness xt_fitness;
	for(int i=0;i<problem_size;i++) {
		xt[i]=(interpolytion[1][i]+s*v[i]*(interpolytion_fitness[0]-interpolytion_fitness[2])/(2.0))/(interpolytion_fitness[0]-2*interpolytion_fitness[1]+interpolytion_fitness[2]);
	    if(xt[i]>max_region) xt[i]=max_region;
	    if(xt[i]<min_region) xt[i]=min_region;
	}
	evaluateIndividual(xt,xt_fitness);
	eva+=1;
	if(xt_fitness<interpolytion_fitness[1]) {
		result_fitness=xt_fitness;
	    for(int i=0;i<problem_size;i++) {
	    	result[i]=xt[i];
	    }
	}
	else {
		result_fitness=interpolytion_fitness[1];
	    for(int i=0;i<problem_size;i++) {
	    	result[i]=interpolytion[1][i];
	    }
	}
	if(x00_fitness<result_fitness) {
		result_fitness=x00_fitness;
	    for(int i=0;i<problem_size;i++) {
	    	result[i]=x00[i];
	    }
	}
	return eva;
}

int searchAlgorithm::DSCG_search(Individual u, Fitness& u_fitness,int dimensions, double s0, double epxilon_DSCG_search, double anpha,int maxEva) {
	//u la ca the can ap dung local search
	//u_fitness la fitness ban dau cua ca the
	//s0 la buoc nhay ban dau cua thuat toan
	//epxilon
	//anpha la he so giam cua s0
	//maxEva la so luong danh gia toi da
	
	//init x,k,i,s,tol
	int n=dimensions;
	vector<Individual> x(n+2);
	vector<Fitness> x_fitness(n+2,0);
	for(int id=0;id<=n+1;id++) {
	    x[id]=new variable[n];
	}
	for(int id=0;id<n;id++) x[0][id]=u[id];
	x_fitness[0]=u_fitness;
	//check x
	//for(int id=0;id<=n+1;id++) {
	//	for(int j=0;j<n;j++) {
	//		cout<<x[id][j]<<" ";
	//	}
	//	cout<<x_fitness[id]<<endl;
	//}
	//end
	int k=1,i=1,eva=0;
	double s=s0,tol=epxilon_DSCG_search;
	//init set of direction
	vector<vector<variable> > direction(n+1);
	for(int id=0;id<n+1;id++) {
		for(int j=0;j<n;j++) {
			if(j==id) direction[id].push_back(1);
			else direction[id].push_back(0);
		}
	}
	//loop
	while(1) {
		//line search
		do {
			int cur=eva;
			eva+=line_search_inter(x[i-1],x_fitness[i-1],x[i],x_fitness[i],s,direction[i-1]);
			if(eva>maxEva) {
				eva=cur;
				break;
			}
			// update the best
			if(x_fitness[i]<u_fitness) {
				u_fitness=x_fitness[i];
            	for(int id=0;id<problem_size;id++) u[id]=x[i][id];
			}
			//consider the evaluations
			if(eva>maxEva) break;
			i=i+1;
		} while (i<=n);
		
		//consider the change
		double sizes=0;
    	Individual z=new double[n];
    	for(int id=0;id<n;id++) {
    		z[id]=x[n][id]-x[0][id];
    		sizes+=pow(z[id],2.0);
		}
		sizes=pow(sizes,0.5);
		
		// if no change
		if(sizes<=0) {
			for(int id=0;id<n;id++) x[n+1][id]=x[n][id];
			x_fitness[n+1]=x_fitness[n];
			s=s*0.1;
		// if s<tol break
		if (s<=tol) break;
		else {
			//next loop
			k=k+1;
			i=1;
			for(int id=0;id<n;id++) x[0][id]=x[n+1][id];
			x_fitness[0]=x_fitness[n+1];
			continue;
		    }
		}
		
		else {
			//if change
			//direction[n]=z/|z|
			for(int id=0;id<n;id++) {
				direction[n][id]=z[id]/sizes;
			}
			i=n+1;
			int cur=eva;
			eva+=line_search_inter(x[n],x_fitness[n],x[n+1],x_fitness[n+1],s,direction[n]);
			if(eva>maxEva) {
				eva=cur;
				break;
			}
			// update the best
			if(x_fitness[n+1]<u_fitness) {
				u_fitness=x_fitness[n+1];
            	for(int id=0;id<problem_size;id++) u[id]=x[n+1][id];
			}
			//consider the evaluations
			//
			
			double sizess=0;
    	    Individual zz=new double[n];
        	for(int id=0;id<n;id++) {
    	    	zz[id]=x[n+1][id]-x[0][id];
    	    	sizess+=pow(zz[id],2.0);
		    }
	    	sizess=pow(sizess,0.5);
	    	if(sizes<s) {
	    		s=s*anpha;
	        	// if s<tol break
	        	if (s<=tol) break;
	        	else {
		    	//next loop
		    	k=k+1;
		    	i=1;
		    	for(int id=0;id<n;id++) x[0][id]=x[n+1][id];
		    	x_fitness[0]=x_fitness[n+1];
		    	continue;
		    }
			}
			else {
				//GramSchmidtOrthogonalization
				GramSchmidtOrthogonalization(direction);
				//next loop
				k=k+1;
			    i=2;
			    //x[0]=x[n]
			    for(int id=0;id<n;id++) x[0][id]=x[n][id];
		    	x_fitness[0]=x_fitness[n];
		    	//x[1]=x[n+1]
		    	for(int id=0;id<n;id++) x[1][id]=x[n+1][id];
		    	x_fitness[1]=x_fitness[n+1];
			    continue;
			}
		}
	}
	//u=x 
	u_fitness=x_fitness[n+1];
	for(int id=0;id<problem_size;id++) u[id]=x[n+1][id];
	return eva;
}

void searchAlgorithm::GramSchmidtOrthogonalization(vector<vector<variable> >& direction) {
	for(int i=0;i<direction.size();i++) {
		for(int j=0;j<i;j++) {// for each previous direction
	    	// calculate proj
			double p1=0,p2=0;
			for(int k=0;k<problem_size;k++) {
				p1+=direction[i][k]*direction[j][k];
				p2+=direction[j][k]*direction[j][k];
			}
		    // update the direction[i]
		    double sizes=0;
		    for(int k=0;k<problem_size;k++) {
				direction[i][k]-=(p1/p2)*direction[j][k];
			}
		}
		double sizes=0;
		for(int k=0;k<problem_size;k++) {
			sizes+=pow(direction[i][k],2.0);
	    }
	    sizes=pow(sizes,0.5);
	    for(int k=0;k<problem_size;k++) {
			direction[i][k]/=sizes;
	    }
	}
}



