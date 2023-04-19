DATA_SECTION
 //Bohai mantis shrimp model
 
 init_int start_mo
 init_int end_mo
 init_int size_n
 init_vector sizes(1,size_n)
  
 init_int surv_dat_mo
 init_ivector surv_months(1,surv_dat_mo)
 init_vector surv_obs_n(1,surv_dat_mo)
 init_matrix surv_size_comp_obs(1,surv_dat_mo,1,size_n)
 
 init_int cat_dat_mo
 init_ivector cat_months(1,cat_dat_mo)
 init_vector cat_obs_n(1,cat_dat_mo)
 init_matrix cat_size_comp_obs(1,cat_dat_mo,1,size_n)

 init_int rec_n
 init_ivector rec_month(1,rec_n)
 
 !!cout<<rec_month<<endl;
 
 init_int rec_n_yr
 init_ivector rec_month_yr(1,rec_n_yr)
 
 init_int growth_n
 init_ivector grow_month(1,growth_n)
 
 init_number catch_cv
 init_number surv_cv
 
 init_number catch_eff_n
 init_number surv_eff_n
 
 !!cout<<surv_eff_n<<endl;
 
 init_number m_mu_prior
 init_number m_sd_prior
 
 init_number growth_var_phase
  !!cout<<m_sd_prior<<endl;
 
   int ipass;
   int proj_n;
  !!proj_n = 2000;
 
PARAMETER_SECTION
  init_bounded_vector log_n_at_l_init(1,size_n,1,30,1)
  init_bounded_vector log_recruits(1,rec_n,1,30,1)
  init_bounded_number nat_m(0.01,4,2)
 
  init_bounded_number surv_q(0.001,1,-3)
  init_bounded_number surv_50(20,100,1)
  init_bounded_number surv_slope(0.01,2,1) 
  
  init_bounded_number fish_50(20,100,1)
  init_bounded_number fish_slope(0.01,2,1)
  init_bounded_vector f_mort(1,cat_dat_mo,0.0001,5,1) 
 
  init_bounded_number growth_alpha(0,10,3)
  init_bounded_number growth_beta(0.0001,5,growth_var_phase)
  init_bounded_number growth_slope(1,1.5,4)
  
  init_bounded_number molt_50(20,200,-1)
  init_bounded_number molt_slope(-2,2,-1)

  //==reference point stuff
  //init_bounded_vector month_f35(1,12)
 
  objective_function_value f
 
  matrix n_at_len(start_mo,end_mo+proj_n,1,size_n)
  vector surv_pred(start_mo,end_mo+proj_n)
  matrix surv_size_comp_pred(start_mo,end_mo+proj_n,1,size_n)
  
  vector catch_pred(start_mo,end_mo+proj_n)
  matrix catch_size_comp_pred(start_mo,end_mo+proj_n,1,size_n)
  matrix pred_c_at_len(start_mo,end_mo+proj_n,1,size_n)
  
  vector fish_sel(1,size_n)
  vector surv_sel(1,size_n)  
  vector molt_prob(1,size_n) 
  vector post_molt_size(1,size_n)
  matrix size_trans_matrix(1,size_n,1,size_n)
  matrix applied_f(start_mo,end_mo+proj_n,1,size_n)
  vector add_recruits(start_mo,end_mo+proj_n)
  vector mean_length(1,size_n)
   
 number catch_like
 number surv_like
 number catch_size_comp_like
 number surv_size_comp_like
  number survey_sel_smooth
  number rec_smooth
  number init_smooth
  number f_smooth
  number IsB0
  number fut_mort
  number fut_rec
  
  number m_prior_like
   
//==============================================================================
PROCEDURE_SECTION

// make fishery selectivity
 for(int size=1;size<=size_n;size++)
   fish_sel(size) = 1/(1+mfexp(-fish_slope*(sizes(size)-fish_50)));

// make survey selectivity
 for(int size=1;size<=size_n;size++)
   surv_sel(size) = 1/(1+mfexp(-surv_slope*(sizes(size)-surv_50)));

// make moltin probability
 for(int size=1;size<=size_n;size++)
   molt_prob(size) = 1/(1+mfexp(-molt_slope*(sizes(size)-molt_50)));

// make fishing mortality
 applied_f.initialize();
  for(int time=1;time<=cat_dat_mo;time++)
	 for(int size=1;size<=size_n;size++)
	  applied_f(cat_months(time),size) = f_mort(time)*fish_sel(size);

// make recruits
  add_recruits.initialize();
   for(int time=1;time<=rec_n;time++)
	   add_recruits(rec_month(time)) = mfexp(log_recruits(time));

// make size transition matrix
  int ilen,il2;
  dvariable devia,alpha,lensum;

  //mean_length is the expected size after molting from a given size bin
  for(ilen=1;ilen<=size_n;ilen++)
      mean_length(ilen)= growth_alpha + growth_slope * sizes(ilen);
   
// using Gamma function for transition matrix
// devia is the bounds of growth bins to evaluate
// the gamma function (x) in prop = integral(i1 to i2) g(x|alpha,beta) dx
// alpha and growth_beta are parameters 
// alpha is the mean growth increment per molt for some premolt length class
// alpha = mean growth increment per molt divided by beta
// beta is the shape parameter - larger beta - more variance 

   for(ilen=1;ilen<=size_n;ilen++)
    {
     // subract the 2.5 from the midpoint of the length bin to get the lower bound
     alpha = (mean_length(ilen)-(sizes(ilen)-2.5))/growth_beta;
    
     lensum = 0;
     for(il2=1;il2<=size_n;il2++)
      if(il2>=ilen)
       {
        devia = sizes(il2)+2.5-sizes(ilen);
        size_trans_matrix(ilen,il2) = mfexp((alpha-1.0)*log(devia)-devia/growth_beta);
        lensum += size_trans_matrix(ilen,il2);
       }  

    //standardize so each row sums to 1.0
    for(il2=1;il2<=size_n;il2++)
      size_trans_matrix(ilen,il2)=size_trans_matrix(ilen,il2)/lensum;
   }

// do pop dy
 get_num_at_len();
  cout<<"obj_fun"<<endl;
 evaluate_the_objective_function();

//========================
FUNCTION get_num_at_len
 
 // initial year numbers at size and selectivity
 for(int size=1;size<=size_n;size++)
   n_at_len(start_mo,size) = exp(log_n_at_l_init(size));

 for(ipass=start_mo;ipass<=end_mo;ipass++) 
    get_num_at_len_yr();

//==========================
FUNCTION get_num_at_len_yr
 int i,j,flag;
 dvar_vector temp_n(1,size_n);
 dvar_vector no_molters(1,size_n);
 dvar_vector molted(1,size_n);
 i = ipass;

 // calculate survey quantities
     surv_size_comp_pred(i) = elem_prod(n_at_len(i),surv_sel);
	 surv_pred(i)=0;
	 for(j=1;j<=size_n;j++)
		surv_pred(i) += surv_size_comp_pred(i,j);
	 for(j=1;j<=size_n;j++)
	 surv_size_comp_pred(i,j) = surv_size_comp_pred(i,j)/surv_pred(i);
    // natural mortality
	temp_n  = n_at_len(i)*mfexp(-nat_m/12);
	
	// fishing mortality
	pred_c_at_len(i) = elem_prod(temp_n,1-mfexp(-applied_f(i)));
	temp_n = elem_prod(temp_n,mfexp(-applied_f(i)));

	// growth
	flag=0;
	for(j=1;j<=growth_n;j++)
	{	 if((i%12+1)==grow_month(j)) flag = 1;}
	if(flag)
	{
	 no_molters = elem_prod(temp_n,1-molt_prob);
	 molted = elem_prod(temp_n,molt_prob);
	 molted = molted * size_trans_matrix;
	 for(j=1;j<=size_n;j++)
      temp_n(j) = no_molters(j) + molted(j);
	}
	// recruitment
	temp_n(1) += add_recruits(i);
	
	// future recruitment
	if(ipass>end_mo)
	{
	flag=0;
	for(j=1;j<=rec_n_yr;j++)
	{	 if((i%12+1)==rec_month_yr(j)) flag = 1;}
	if(flag)
		temp_n(1) += fut_rec;
    }

	n_at_len(i+1) = temp_n;
  
//==============================================================================
FUNCTION evaluate_the_objective_function
  
  f=0;
  
  catch_pred.initialize();
  for(int i=start_mo;i<=end_mo;i++)
   for(int j =1;j<=size_n;j++)
    catch_pred(i) += pred_c_at_len(i,j);
   
  for(int i=start_mo;i<=end_mo;i++)
   for(int j =1;j<=size_n;j++)
	 if(catch_pred(i)>0)
      catch_size_comp_pred(i,j) = pred_c_at_len(i,j)/catch_pred(i);

  //likelihoods
  // catch abundance
  catch_like = 0;
  for(int time=1;time<=cat_dat_mo;time++)
     catch_like += square( log(catch_pred(cat_months(time))) - log(cat_obs_n(time))) / (2.0 * square(catch_cv));

 f += catch_like;

  // catch size composition data
  catch_size_comp_like = 0;
   for(int time=1;time<=cat_dat_mo;time++)
	 for(int size=1;size<=size_n;size++)
      if (cat_size_comp_obs(time,size) >0.001)
       catch_size_comp_like += catch_eff_n*(cat_size_comp_obs(time,size)) * log( catch_size_comp_pred(cat_months(time),size)/ cat_size_comp_obs(time,size));
  catch_size_comp_like = -1*catch_size_comp_like;
 
 f += catch_size_comp_like;
 
  //=================================
 // survey abundance
  //=================================
  surv_like = 0;
  for(int time=1;time<=surv_dat_mo;time++)
     surv_like += square( log(surv_pred(surv_months(time))) - log(surv_obs_n(time))) / (2.0 * square(surv_cv));
 
  f += surv_like;
 
  // survey size composition data
  surv_size_comp_like = 0;
   for(int time=1;time<=surv_dat_mo;time++)
	 for(int size=1;size<=size_n;size++)
      if (surv_size_comp_obs(time,size) >0.001)
       surv_size_comp_like += surv_eff_n*(surv_size_comp_obs(time,size)) * log( surv_size_comp_pred(surv_months(time),size)/ surv_size_comp_obs(time,size));
  surv_size_comp_like = -1*surv_size_comp_like;

  f += surv_size_comp_like;


  init_smooth =0;
  init_smooth = 0.1* (norm2(first_difference(first_difference(log_n_at_l_init))));
     f += init_smooth;

  f_smooth =0;
  f_smooth = 0.1* (norm2(first_difference(first_difference(f_mort))));
     f += f_smooth;
	 
  m_prior_like = 0;
	m_prior_like = square(m_mu_prior-nat_m) / (2.0*square(m_sd_prior));	
	
	f += m_prior_like;
	 
  cout<<surv_size_comp_like<< " " << surv_like << " " << catch_size_comp_like << " " << catch_like << " " << init_smooth << " " << m_prior_like<<endl;
  
 // ==========================================================================

FUNCTION get_fut_mortality
 int i;
  for (i=ipass;i<=end_mo+proj_n;i++)
  {
   if (IsB0 == 0)
    {
   	applied_f(i) = 0; 
	}
   else 
    {
	applied_f(i) = fut_mort; 
	}
    applied_f(i) = elem_prod(applied_f(i),fish_sel); 	
   }
  
//==============================================================================
FUNCTION find_bzero
  int nn; int j;
  // Find B0
  IsB0 = 0;
  fut_rec = 0;
  nn = 0;

 //find average recruitment by month...these could be different
 for(j=1;j<rec_n;j++)
  {
    fut_rec += mfexp(log_recruits(j));
    nn += 1;
  }
  fut_rec = fut_rec/nn;
  cout<<"futrec"<<fut_rec<<endl;
  fut_mort = 0;
  ipass = end_mo+1;
  get_fut_mortality();
   for (ipass=end_mo+1;ipass<end_mo+proj_n;ipass++) get_num_at_len_yr(); 
  
// ========================y==================================================   
REPORT_SECTION
  find_bzero();
  
  report<<"$likelihoods"<<endl;
  report<<surv_size_comp_like<< " " << surv_like << " " << catch_size_comp_like << " " << catch_like <<endl;
  report <<"$recruits" << endl;
  report << mfexp(log_recruits)<<endl;
  report <<"$sizes" << endl;
  report << sizes<<endl;
    report <<"$growth_inc" << endl;
  report << mean_length<<endl;
  
 //===survey inputs and ouputs
  report<<"$survey_months"<<endl;
  report<<surv_months<<endl;
  report<<"$survey_pred_all"<<endl;
  report<<surv_pred<<endl;
  report<<"$survey_obs"<<endl;
  report<<surv_obs_n<<endl;
    report<<"$survey_selectivity"<<endl;
  report<<surv_sel<<endl;
      report<<"$fish_selectivity"<<endl;
  report<<fish_sel<<endl;
  
  report<<"$survey_pred_comp" << endl;
  for(int i=1; i<=surv_dat_mo; i++)
  {
    report << surv_pred(surv_months(i))<<endl;
  }

  report <<"$survey_size_comp_pred_all" << endl;
  for(int i=start_mo; i<=end_mo; i++)
  {
    report << surv_size_comp_pred(i)<<endl;
  }
  report <<"$survey_size_comp_pred_comp" << endl;
  for(int i=1; i<=surv_dat_mo; i++)
  {
    report << surv_size_comp_pred(surv_months(i))<<endl;
  }
  report <<"$survey_size_comp_obs" << endl;
  for(int i=1; i<=surv_dat_mo; i++)
  {
    report << surv_size_comp_obs(i)<<endl;
  }

  //===catch inputs and outputs
  report<<"$catch_months"<<endl;
  report<<cat_months<<endl;
  report<<"$catch_pred_all"<<endl;
  report<<catch_pred<<endl;
  report<<"$cat_obs"<<endl;
  report<<cat_obs_n<<endl;
  report<<"$catch_pred_comp" << endl;
  for(int i=1; i<=cat_dat_mo; i++)
  {
    report << catch_pred(cat_months(i))<<endl;
  }

  report <<"$catch_size_comp_pred_all" << endl;
  for(int i=start_mo; i<=end_mo; i++)
  {
    report << catch_size_comp_pred(i)<<endl;
  }

  report <<"$catch_size_comp_pred_comp" << endl;
  for(int i=1; i<=cat_dat_mo; i++)
  {
    report << catch_size_comp_pred(cat_months(i))<<endl;
  }
    cout<<"2"<<endl;
   report <<"$catch_size_comp_obs" << endl;
  for(int i=1; i<=cat_dat_mo; i++)
  {
    report << cat_size_comp_obs(i)<<endl;
  }

  report <<"$applied_f_matrix" << endl;
  for(int i=start_mo; i<=end_mo+proj_n; i++)
  {
    report << applied_f(i)<<endl;
  }
    report <<"$size_trans_matrix" << endl;
  for(int i=1; i<=size_n; i++)
  {
    report << size_trans_matrix(i)<<endl;
  }
 
   report <<"$n_at_len" << endl;
  for(int i=start_mo; i<=end_mo+proj_n; i++)
  {
    report << n_at_len(i)<<endl;
  }

  save_gradients(gradients);

RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

