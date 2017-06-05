# data processing and cleaning script to clean the data up and join the various datasets 
library(plyr) 

# function that takes a row of data frame from 'Hc103burnedb.txt' and turns it into multi-row, one 
# row for each year
detangler = function(row){
  
  years = 2008:2002
  # find the relevent years of data
  sur_ind = which(row[6:12] == 5 | row[6:12] == 1 | row[6:12] == 3)
  sur_years = years[sur_ind]
  death_year = years[which(row[6:12] == 0)]
  life_cycle = c(death_year, sur_years)

  df = data.frame(uID = row$uID, uLoc = row$uLoc, year = life_cycle, sur = NA,
    age = NA, height = NA, height_prev = NA, rep = NA, stems = NA, dup = row$dup)
  df$sur[df$year %in% sur_years] = 1
  df$sur[df$year %in% death_year] = 0
  if(!3 %in% row[6:12]) df$age[df$year %in% sur_years] = (length(sur_years) - 1):0

  # split out the height and reproduction
  height = row[14:20]
  repro = row[22:28]
  stems = row[30:36]
  df$height[df$year %in% sur_years] = as.numeric(height[sur_ind])
  df$rep[df$year %in% sur_years] = as.numeric(repro[sur_ind])
  df$stems[df$year %in% sur_years] = as.numeric(stems[sur_ind])
  
  if(length(life_cycle) > 1){
    h2 = as.numeric(df$height[2:length(life_cycle)])
    df$height_prev[1:(length(life_cycle) - 1)] = h2
  }
  
  return(df)
  
}

setwd('/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data')

# get the post burn vital rate data
vr_post_burn = read.table('Hc103burnedb.txt', header = TRUE)
vr_post_burn$uID = paste0(vr_post_burn$gap, ':', vr_post_burn$tag)

#series of tests to make sure the data is as expected
# check there is no dormancy, which will complicate things
dorm_test = sapply(1:dim(vr_post_burn)[1], FUN = function(x) length(which(vr_post_burn[x, 6:12] == 0)))
sum(dorm_test > 1)
# sum == 0, so no dormancy, which is what we want 

#test to find duplicates
vr_post_burn[duplicated(vr_post_burn$uID), ]
# there are about 30 duplices, but looking at the data they seem like seperate data points, the problem will be 
# assigning a location and kinship to them so mark them as duplicated in the data set, and modify their unique ID
dupes = which(duplicated(vr_post_burn$uID))
row_ind = 1:dim(vr_post_burn)[1]
vr_post_burn$dup = FALSE
vr_post_burn$dup[vr_post_burn$uID %in% vr_post_burn$uID[dupes]] = TRUE
vr_post_burn$uID[vr_post_burn$dup] = paste0(vr_post_burn$uID[vr_post_burn$dup], '.', row_ind[vr_post_burn$dup])
vr_post_burn$uLoc = paste0(vr_post_burn$gap, ':', vr_post_burn$tag)

vr_stack = detangler(vr_post_burn[1, ])
for(i in 2:dim(vr_post_burn)[1]){
  vr_stack = rbind(vr_stack, detangler(vr_post_burn[i, ]))
}

# some basic data exploration
pro_mort_obs = sum(vr_stack$sur == 0) / length(unique(vr_stack$uID))
# looks like we see mortality in about 40% of the individuals
pro_repro = sum(vr_stack$rep, na.rm = TRUE) / dim(vr_stack)[1]
# arround 50% of observations have reproduction

# bring in the space and genotype information
loc_gen = read.csv('location_genotype.csv', header = TRUE, stringsAsFactor = FALSE)
#remove rows where location is NA
loc_gen = loc_gen[!is.na(loc_gen$X), ]

#check the tag IDs match up
prod(loc_gen$tag == loc_gen$tag.1) # = 1, they all do

loc_gen$uLoc = paste0(loc_gen$gap, ':', loc_gen$tag)

# beacuse of the duplicate locations I need to hand roll a dataframe merger 
vr_loc_gen = data.frame(vr_stack, gap = NA, tag = NA, X = NA, Y = NA, MNR = NA, 
  MDH1 = NA, MDH3 = NA, X6PGD = NA, IDH = NA, num_years = NA)

# How many years of data does each individual has, those with only 1 year cannot be used for survival and growth  
year_count = ddply(vr_loc_gen, .(uID), summarise, num_years = length(uID))
hist(year_count$num_years) 
sum(year_count$num_years > 1) / dim(year_count)[1]  # 85% of individuals have more than 1 year of data, median number = 3  

for(i in 1:dim(loc_gen)[1]){
  vr_loc_gen[vr_loc_gen$uLoc %in% loc_gen$uLoc[i], 
    c('gap', 'tag', 'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')] = loc_gen[i, 
    c('gap', 'tag', 'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]
}
  
for(i in 1:dim(year_count)[1]){
    vr_loc_gen$num_years[vr_loc_gen$uID %in% year_count$uID[i]] = year_count$num_years[i]
}
  
write.csv(vr_loc_gen, file = 'vr_loc_gen_postburn.csv')

  
## make survival data  
#need location data so drop the obervations that do not have a location
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X),]

sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]

# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]

# recode year and patch to act as indicies 
sur_dat$year_num = (sur_dat$year - (max(sur_dat$year))) + max(abs(sur_dat$year - (max(sur_dat$year)))) + 1    
sur_dat$gapID = sapply(strsplit(sur_dat$uID, split = ':', fixed = TRUE), FUN = function(x) x[1])
inds = unique(sur_dat$gapID)
sur_dat$gapID_num = NA
for(i in 1:length(inds)) sur_dat$gapID_num[sur_dat$gapID %in% inds[i]] = i

# height: this requires some assumptions because 1) individuals change size over the course of a year, 
# 2) we do not have the final size for indivudals that died. 3) we have indivduals that were only found when they were very large.
# Best we can do here is assue that mortality occurs just after census, so that we can use previous height as a predictor,
# and in the case of indivduals in their first year we say they were 0 
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]
sur_dat$height_mod = sur_dat$height_prev - mean(sur_dat$height_prev)

write.csv(sur_dat, file = 'sur_loc_gen_postburn.csv')

  
## make reproduction data  
#need location data so drop the obervations that do not have a location
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

rep_dat = vr_loc_gen[!is.na(vr_loc_gen$rep), c('uID', 'uLoc', 'rep', 'year', 'height', 'X', 'Y')]
rep_dat = rep_dat[!is.na(rep_dat$height),]
rep_dat = rep_dat[rep_dat$rep <= 1, ]

# recode year and patch to act as indicies 
rep_dat$year_num = (rep_dat$year - (max(rep_dat$year))) + max(abs(rep_dat$year - (max(rep_dat$year)))) + 1    
rep_dat$gapID = sapply(strsplit(rep_dat$uID, split = ':', fixed = TRUE), FUN = function(x) x[1])
inds = unique(rep_dat$gapID)
rep_dat$gapID_num = NA
for(i in 1:length(inds)) rep_dat$gapID_num[rep_dat$gapID %in% inds[i]] = i

rep_dat$height_mod = rep_dat$height - mean(rep_dat$height)

write.csv(rep_dat, file = 'rep_loc_gen_postburn.csv')

# growth data
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

gr_dat = vr_loc_gen[!is.na(vr_loc_gen$height), c('uID', 'uLoc', 'rep', 'year', 'height', 'height_prev', 'X', 'Y')]
gr_dat = gr_dat[!is.na(gr_dat$height_prev), ]

write.csv(gr_dat, file = 'gr_loc_gen_postburn.csv')
