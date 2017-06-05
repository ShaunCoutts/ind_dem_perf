# Script to run models on severs ect.

## on Zhu machine (first ssh into Zhu machine)

# no space sur model
nohup Rscript ~/ind_perf_pred/scripts/NS_survial_run.R > ~/ind_perf_pred/Rscript_output/NS_sur_cons.out 2>&1 &

# SPP sur model
nohup Rscript ~/ind_perf_pred/scripts/sur_SPP10_stan.R > ~/ind_perf_pred/Rscript_output/SPP_NG_sur_cons.out 2>&1 &

# SPP sur model
nohup Rscript ~/ind_perf_pred/scripts/SPP_NK_sur_stan.R > ~/ind_perf_pred/Rscript_output/SPP_NK_sur_cons.out 2>&1 &



# SA sur model
nohup Rscript ~/ind_perf_pred/scripts/SA_survial_run.R > ~/ind_perf_pred/Rscript_output/SA_sur_cons.out 2>&1 &

# SE sur model
nohup Rscript ~/ind_perf_pred/scripts/sur_SE_stan_run.R > ~/ind_perf_pred/Rscript_output/SE_sur_cons.out 2>&1 &




## set up a run on iceberg
# survival SPP_NK
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_sur_jobcall.sh
  
# survival SPP_25
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_sur25_jobcall.sh

# survival SPP_NK25
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_sur25NN_jobcall.sh

# reproduction SPP_NK
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_rep2_jobcall.sh

# reproduction SPP_25
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_rep25_jobcall.sh

# reproduction SPP_NK
ssh -X bo1src@iceberg.sheffield.ac.uk
qsub /home/bo1src/ind_dem_perf/scripts/ICEBERG_grNK_jobcall.sh
