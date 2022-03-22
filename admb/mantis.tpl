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
 init_int growth_n
 init_ivector grow_month(1,growth_n)
 
   int ipass;
   int proj_n;
  !!proj_n = 100;
 
PARAMETER_SECTION
  init_bounded_vector log_n_at_l_init(1,size_n,1,30,-1)
  init_bounded_vector log_recruits(1,rec_n,1,30,-1)
  init_bounded_number nat_m(0.01,4,-2)
 
  init_bounded_number surv_q(0.001,1,-3)
  init_bounded_vector survey_select(1,size_n,0.00001,1,-2)
  
  init_bounded_number fish_50(20,80,-1)
  init_bounded_number fish_slope(0.01,2,-1)
  init_bounded_vector f_mort(1,cat_dat_mo,0.0001,5,-1) 
 
  init_bounded_number growth_alpha(0,10,-3)
  init_bounded_number growth_beta(0,2,-3)
  init_bounded_number growth_slope(1,1.5,-4)
  init_number dummy
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
  vector post_molt_size(1,size_n)
  matrix size_trans_matrix(1,size_n,1,size_n)
  matrix applied_f(start_mo,end_mo+proj_n,1,size_n)
  vector add_recruits(start_mo,end_mo+proj_n)
  vector mean_length(1,size_n)
   
  number survey_ind_like
  number size_comp_like
  number survey_sel_smooth
  number rec_smooth

   
//==============================================================================
PROCEDURE_SECTION

// make fishery selectivity
 for(int size=1;size<=size_n;size++)
   fish_sel(size) = 1/(1+mfexp(-fish_slope*(sizes(size)-fish_50)));

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
 i = ipass;

 // calculate survey quantities
     surv_size_comp_pred(i) = elem_prod(n_at_len(i),survey_select);
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
		temp_n = temp_n * size_trans_matrix;

	// recruitment
	temp_n(1) += add_recruits(i);
	
	n_at_len(i+1) = temp_n;
  
//==============================================================================
FUNCTION evaluate_the_objective_function
  
  catch_pred.initialize();
  for(int i=start_mo;i<=end_mo;i++)
   for(int j =1;j<=size_n;j++)
    catch_pred(i) += pred_c_at_len(i,j);
   
  for(int i=start_mo;i<=end_mo;i++)
   for(int j =1;j<=size_n;j++)
	 if(catch_pred(i)>0)
      catch_size_comp_pred(i,j) = pred_c_at_len(i,j)/catch_pred(i);
  
  


  // likelihoods
  // imm_num_like = 0;
  // for (int year=styr;year<=endyr;year++)
    // imm_num_like += square( log(imm_numbers_pred(year)) - log(imm_n_obs(year))) / (2.0 * square(sigma_numbers_imm(year)));

 
 
  // immature numbers at size data
  // imm_like = 0;
  // for (int year=styr;year<=endyr;year++)
   // for (int size=1;size<=size_n;size++)
    // if (imm_n_size_obs(year,size) >0)
     // imm_like += imm_eff_samp*(imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)) * log( (selectivity(year,size)*imm_n_size_pred(year,size)/imm_numbers_pred(year)) / (imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)));
  // imm_like = -1*imm_like;
  
   
  // smooth_q_like = 0;
  // smooth_q_like = smooth_q_weight* (norm2(first_difference(first_difference(q_dev))) +norm2(first_difference(first_difference(q_mat_dev)))) ;
  
 
  f = 0;
  
// ========================y==================================================   
REPORT_SECTION
  report<<"$likelihoods"<<endl;
  report <<"$recruits" << endl;
  report << mfexp(log_recruits)<<endl;
  report <<"$sizes" << endl;
  report << sizes<<endl;
 //===survey inputs and ouputs
  report<<"$survey_months"<<endl;
  report<<surv_months<<endl;
  report<<"$survey_pred_all"<<endl;
  report<<surv_pred<<endl;
  report<<"$survey_obs"<<endl;
  report<<surv_obs_n<<endl;
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
  for(int i=start_mo; i<=end_mo; i++)
  {
    report << applied_f(i)<<endl;
  }
    report <<"$size_trans_matrix" << endl;
  for(int i=1; i<=size_n; i++)
  {
    report << size_trans_matrix(i)<<endl;
  }
 

RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

