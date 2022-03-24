#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <mantis.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  start_mo.allocate("start_mo");
  end_mo.allocate("end_mo");
  size_n.allocate("size_n");
  sizes.allocate(1,size_n,"sizes");
  surv_dat_mo.allocate("surv_dat_mo");
  surv_months.allocate(1,surv_dat_mo,"surv_months");
  surv_obs_n.allocate(1,surv_dat_mo,"surv_obs_n");
  surv_size_comp_obs.allocate(1,surv_dat_mo,1,size_n,"surv_size_comp_obs");
  cat_dat_mo.allocate("cat_dat_mo");
  cat_months.allocate(1,cat_dat_mo,"cat_months");
  cat_obs_n.allocate(1,cat_dat_mo,"cat_obs_n");
  cat_size_comp_obs.allocate(1,cat_dat_mo,1,size_n,"cat_size_comp_obs");
  rec_n.allocate("rec_n");
  rec_month.allocate(1,rec_n,"rec_month");
  growth_n.allocate("growth_n");
  grow_month.allocate(1,growth_n,"grow_month");
  catch_cv.allocate("catch_cv");
  surv_cv.allocate("surv_cv");
  catch_eff_n.allocate("catch_eff_n");
  surv_eff_n.allocate("surv_eff_n");
proj_n = 100;
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_n_at_l_init.allocate(1,size_n,1,30,1,"log_n_at_l_init");
  log_recruits.allocate(1,rec_n,1,30,1,"log_recruits");
  nat_m.allocate(0.01,4,2,"nat_m");
  surv_q.allocate(0.001,1,-3,"surv_q");
  surv_50.allocate(20,80,1,"surv_50");
  surv_slope.allocate(0.01,2,1,"surv_slope");
  fish_50.allocate(20,80,1,"fish_50");
  fish_slope.allocate(0.01,2,1,"fish_slope");
  f_mort.allocate(1,cat_dat_mo,0.0001,5,1,"f_mort");
  growth_alpha.allocate(0,10,3,"growth_alpha");
  growth_beta.allocate(0,5,3,"growth_beta");
  growth_slope.allocate(1,1.5,4,"growth_slope");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  n_at_len.allocate(start_mo,end_mo+proj_n,1,size_n,"n_at_len");
  #ifndef NO_AD_INITIALIZE
    n_at_len.initialize();
  #endif
  surv_pred.allocate(start_mo,end_mo+proj_n,"surv_pred");
  #ifndef NO_AD_INITIALIZE
    surv_pred.initialize();
  #endif
  surv_size_comp_pred.allocate(start_mo,end_mo+proj_n,1,size_n,"surv_size_comp_pred");
  #ifndef NO_AD_INITIALIZE
    surv_size_comp_pred.initialize();
  #endif
  catch_pred.allocate(start_mo,end_mo+proj_n,"catch_pred");
  #ifndef NO_AD_INITIALIZE
    catch_pred.initialize();
  #endif
  catch_size_comp_pred.allocate(start_mo,end_mo+proj_n,1,size_n,"catch_size_comp_pred");
  #ifndef NO_AD_INITIALIZE
    catch_size_comp_pred.initialize();
  #endif
  pred_c_at_len.allocate(start_mo,end_mo+proj_n,1,size_n,"pred_c_at_len");
  #ifndef NO_AD_INITIALIZE
    pred_c_at_len.initialize();
  #endif
  fish_sel.allocate(1,size_n,"fish_sel");
  #ifndef NO_AD_INITIALIZE
    fish_sel.initialize();
  #endif
  surv_sel.allocate(1,size_n,"surv_sel");
  #ifndef NO_AD_INITIALIZE
    surv_sel.initialize();
  #endif
  post_molt_size.allocate(1,size_n,"post_molt_size");
  #ifndef NO_AD_INITIALIZE
    post_molt_size.initialize();
  #endif
  size_trans_matrix.allocate(1,size_n,1,size_n,"size_trans_matrix");
  #ifndef NO_AD_INITIALIZE
    size_trans_matrix.initialize();
  #endif
  applied_f.allocate(start_mo,end_mo+proj_n,1,size_n,"applied_f");
  #ifndef NO_AD_INITIALIZE
    applied_f.initialize();
  #endif
  add_recruits.allocate(start_mo,end_mo+proj_n,"add_recruits");
  #ifndef NO_AD_INITIALIZE
    add_recruits.initialize();
  #endif
  mean_length.allocate(1,size_n,"mean_length");
  #ifndef NO_AD_INITIALIZE
    mean_length.initialize();
  #endif
  catch_like.allocate("catch_like");
  #ifndef NO_AD_INITIALIZE
  catch_like.initialize();
  #endif
  surv_like.allocate("surv_like");
  #ifndef NO_AD_INITIALIZE
  surv_like.initialize();
  #endif
  catch_size_comp_like.allocate("catch_size_comp_like");
  #ifndef NO_AD_INITIALIZE
  catch_size_comp_like.initialize();
  #endif
  surv_size_comp_like.allocate("surv_size_comp_like");
  #ifndef NO_AD_INITIALIZE
  surv_size_comp_like.initialize();
  #endif
  survey_sel_smooth.allocate("survey_sel_smooth");
  #ifndef NO_AD_INITIALIZE
  survey_sel_smooth.initialize();
  #endif
  rec_smooth.allocate("rec_smooth");
  #ifndef NO_AD_INITIALIZE
  rec_smooth.initialize();
  #endif
  init_smooth.allocate("init_smooth");
  #ifndef NO_AD_INITIALIZE
  init_smooth.initialize();
  #endif
  f_smooth.allocate("f_smooth");
  #ifndef NO_AD_INITIALIZE
  f_smooth.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  f =0.0;
 for(int size=1;size<=size_n;size++)
   fish_sel(size) = 1/(1+mfexp(-fish_slope*(sizes(size)-fish_50)));
 for(int size=1;size<=size_n;size++)
   surv_sel(size) = 1/(1+mfexp(-surv_slope*(sizes(size)-surv_50)));
 applied_f.initialize();
  for(int time=1;time<=cat_dat_mo;time++)
	 for(int size=1;size<=size_n;size++)
	  applied_f(cat_months(time),size) = f_mort(time)*fish_sel(size);
  add_recruits.initialize();
   for(int time=1;time<=rec_n;time++)
	   add_recruits(rec_month(time)) = mfexp(log_recruits(time));
  int ilen,il2;
  dvariable devia,alpha,lensum;
  //mean_length is the expected size after molting from a given size bin
  for(ilen=1;ilen<=size_n;ilen++)
      mean_length(ilen)= growth_alpha + growth_slope * sizes(ilen);
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
 get_num_at_len();
  cout<<"obj_fun"<<endl;
 evaluate_the_objective_function();
}

void model_parameters::get_num_at_len(void)
{
 // initial year numbers at size and selectivity
 for(int size=1;size<=size_n;size++)
   n_at_len(start_mo,size) = exp(log_n_at_l_init(size));
 for(ipass=start_mo;ipass<=end_mo;ipass++) 
    get_num_at_len_yr();
}

void model_parameters::get_num_at_len_yr(void)
{
 int i,j,flag;
 dvar_vector temp_n(1,size_n);
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
		temp_n = temp_n * size_trans_matrix;
	// recruitment
	temp_n(1) += add_recruits(i);
	n_at_len(i+1) = temp_n;
}

void model_parameters::evaluate_the_objective_function(void)
{
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
  cout<<surv_size_comp_like<< " " << surv_like << " " << catch_size_comp_like << " " << catch_like << " " << init_smooth << endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report<<"$likelihoods"<<endl;
  report<<surv_size_comp_like<< " " << surv_like << " " << catch_size_comp_like << " " << catch_like <<endl;
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
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
