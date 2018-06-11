inventory_wqp <- function(wqp_needs, wqp_variables) {
  
  # identify constituents we need
  constituents <- as.character(unique(wqp_needs$variable))
  
  
  wqp_args <- list(
    characteristicName=NA, # to be updated each time through loop
    siteid=NA #updated for each constituent
  )
  
  # loop over the constituents, getting rows for each
  sample_time <- system.time({
    samples <- bind_rows(lapply(constituents, function(constituent) {
      message(Sys.time(), ': getting inventory for ', constituent)
      wqp_args$characteristicName <- wqp_variables$characteristicName[[constituent]]
      wqp_args$siteid <- wqp_needs$site[wqp_needs$variable %in% constituent]
      tryCatch({
        wqp_wdat <- do.call(whatWQPdata, wqp_args)
        mutate(wqp_wdat, constituent=constituent)
      }, error=function(e) {
        # keep going IFF the only error was that there weren't any matching sites
        if(grepl('arguments imply differing number of rows', e$message)) {
          NULL
        } else {
          stop(e)
        }
      })
    }))
  })
  message(sprintf('sample inventory complete, required %0.0f seconds', sample_time[['elapsed']]))
  
  return(samples)
}

wqp_partition <- function(wqp_inventory, wqp_nrecords_chunk) {
  partitions <- bind_rows(lapply(unique(wqp_inventory$constituent), function(temp_constituent) {
    # an atomic group is a combination of parameters that can't be reasonably
    # split into multiple WQP pulls - in this case we're defining atomic
    # groups as distinct combinations of constituent (a group of
    # characteristicNames) and site ID. 
  
    # right now, just grouping by unique site ID (MonitoringLocationIdentifier), however, if we have crosswalk available at this step,
    # could put in the unique NHDID for each lake to ensure each lake gets in same file
    atomic_groups <-  wqp_inventory %>%
      filter(constituent == temp_constituent) %>%
      group_by(MonitoringLocationIdentifier) %>%
      summarize(NumObs=sum(resultCount)) %>%
      arrange(desc(NumObs))
    
    # split the full pull (combine atomic groups) into right-sized partitions
    # by an old but fairly effective paritioning heuristic: pick the number of
    # partitions desired, sort the atomic groups by descending size, and then
    # go down the list, each time adding the next atomic group to the
    # partition that's currently smallest
    
    # first, pull out sites that by themselves are larger than the cutoff number
    # for these sites, we will "break" the chunk rule, and get in a single pull from WQP
    single_site_partitions <- filter(atomic_groups, NumObs >= wqp_nrecords_chunk)
    n_single_site_partitions <- nrow(single_site_partitions)
    
    multi_site_partitions <- filter(atomic_groups, NumObs < wqp_nrecords_chunk)
    
    num_partitions <- ceiling(sum(multi_site_partitions$NumObs) / wqp_nrecords_chunk)
    
    partition_sizes <- rep(0, num_partitions)
    assignments <- rep(0, nrow(multi_site_partitions))
    for(i in seq_len(nrow(multi_site_partitions))) {
      size_i <- multi_site_partitions[[i,"NumObs"]]
      smallest_partition <- which.min(partition_sizes)
      assignments[i] <- smallest_partition
      partition_sizes[smallest_partition] <- partition_sizes[smallest_partition] + size_i
    }
    
    last_assignment <- max(assignments)
    
    single_site_partitions <- single_site_partitions %>%
      mutate(
        Constituent=constituent,
        PullTask=sprintf('%s_%s_%03d', 'WQP', constituent, seq(last_assignment + 1, last_assignment + nrow(single_site_partitions))))
    
    # create a filename column
   
    multi_site_partitions <- multi_site_partitions %>%
      mutate(
        Constituent=constituent,
        PullTask=sprintf('%s_%s_%03d', 'WQP', constituent, assignments))
    
    partitions <- bind_rows(multi_site_partitions, single_site_partitions)
    
    return(partitions)
    
  }))
}