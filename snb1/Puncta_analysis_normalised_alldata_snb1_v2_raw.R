library(readxl)
library(writexl)
library(ggplot2)
library(pracma)
library(magrittr)
library(tidyverse)
library(gsignal)
library(patchwork)

### NOTE:
# for this script, the sheet name to be analysed has to be manually set in the lapply call at the bottom
# typically this is 'axon' or 'dendrite'

# define functions to detect peaks and plot the outputs
detect_peaks <- function(x,intensities){
  
  target_worm <- intensities %>% 
    dplyr::filter(Worm == x) %>% 
    drop_na()
  
  sd_y <- sd(target_worm$Grey_Value_norm, na.rm = TRUE) # calculate standard deviation
  
  # we need all positive values for peak detection so make a grey_value column normalised to 0
  # min_value <- abs(min(target_worm$Grey_Value_norm))
  # target_worm$Grey_Value_norm_abs <- target_worm$Grey_Value_norm + min_value
  
  # peaks_median <- tibble(findpeaks(target_worm$Grey_Value_norm, minpeakheight = median(target_worm$Grey_Value))) # note we use the median of the un-normalised Grey Value - equivalent to using a threshold of median * 2
  peaks_sd <- as_tibble(pracma::findpeaks(target_worm$Grey_Value_norm, minpeakheight = sd_y))
  
  # colnames(peaks_median) <- c("peak_height","peak_max","peak_start","peak_end")
  colnames(peaks_sd) <- c("peak_height","peak_max","peak_start","peak_end")
  
  # for every identified peak pull the relevant intensity data and calculate FWHM and AUC
  fwhm_vec <- vector(length = nrow(peaks_sd))
  auc_vec <- vector(length = nrow(peaks_sd))
  peak_vec <- vector(length = nrow(peaks_sd))
  start_vec <- vector(length = nrow(peaks_sd))
  end_vec <- vector(length = nrow(peaks_sd))

  for(i in 1:length(fwhm_vec)){
    peak <- target_worm[peaks_sd$peak_start[i]:peaks_sd$peak_end[i],]
    
    fwhm_vec[i] <- fwhm(peak$Distance,peak$Grey_Value_norm)
    auc_vec[i] <- trapz(peak$Distance,peak$Grey_Value_norm)
    peak_vec[i] <- target_worm$Distance[peaks_sd$peak_max[i]]
    start_vec[i] <- target_worm$Distance[peaks_sd$peak_start[i]]
    end_vec[i] <- target_worm$Distance[peaks_sd$peak_end[i]]
  }
  
  peaks_sd$fwhm <- fwhm_vec
  peaks_sd$auc <- auc_vec
  peaks_sd$worm <- x  
  peaks_sd$peak_max_dist <- peak_vec
  peaks_sd$peak_start_dist <- start_vec
  peaks_sd$peak_end_dist <- end_vec
  # calculate the puncta intensity - median ratio
  # note we add median value back to the puncta intensity to undo the normalisation carried out in the master function
  peaks_sd$median_ratio <- (peaks_sd$peak_height + median(target_worm$Grey_Value))/median(target_worm$Grey_Value)
  
  return(peaks_sd)
}

puncta_per_um <- function(x,intensities,peaks){
  
  target_worm <- intensities %>% 
    dplyr::filter(Worm == x) %>% 
    drop_na()
  
  target_peaks <- peaks %>% 
    dplyr::filter(worm == x) %>% 
    drop_na()
  
  peak_density <- nrow(target_peaks)/max(target_worm$Distance)
  
  return(peak_density)
}

plotting_func <- function(x,intensities,peaks){
  
  target_worm <- intensities %>% 
    dplyr::filter(Worm == x) %>% 
    drop_na()
  
  target_peak <- peaks %>% 
    dplyr::filter(worm == x) %>% 
    drop_na()
  
  plot_worm <- ggplot(target_worm, aes(x = Distance, y = Grey_Value_norm)) +
    geom_line() +
    geom_point(data = target_peak, aes(x = peak_max_dist, y = peak_height), color = "black",fill = "firebrick2", size = 2,shape=21,alpha = 0.5) +
    labs(title = paste0(str_to_title(str_extract(x,pattern="[a-z]+"))," ",str_extract(x,pattern="\\d+")),
         x = "Distance (µm)",
         y = "Normalised Intensity (A.U.)") +
    scale_y_continuous(limits = c(-1000,4000)) +
    theme_minimal(base_size = 12)
  
  return(plot_worm)
}

# define master function to apply analysis to every source_data
analyse_intensities <- function(source_data,sheet){
  
  # load data and normalise
  raw_intensities <- read_xlsx(source_data,sheet = sheet)  
  
  source_data_name <- str_sub(source_data,start = 1, end = -6)
  
  num_worms <- ncol(raw_intensities)/2
  
  intensities <- raw_intensities %>%
    pivot_longer(cols = everything(),names_to = c(".value","Distance"),names_sep = "e_")
  
  colnames(intensities) <- c("Worm","Distance","Grey_Value")
  
  # arrange worm names in numeric order to allow handling of any number of non-consecutive worm numbers
  worm_names <- tibble(worm = unique(intensities$Worm),worm_num = as.numeric(str_extract(unique(intensities$Worm),"\\d+")))
  
  worm_names %<>%
    arrange(worm_num)
  
  intensities$Worm <- factor(intensities$Worm,levels = worm_names$worm)
  
  intensities %<>%
    group_by(Worm) %>% 
    mutate(Grey_Value_norm = Grey_Value - median(Grey_Value,na.rm=T))
  
  # # quick visualisation of intensity over neuron length
  # ggplot(intensities,aes(x = Distance, y = Grey_Value_norm)) +
  #   geom_line() +
  #   scale_y_continuous(breaks = c(-4000,0,4000,8000,12000,16000))+
  #   theme_minimal()+
  #   facet_wrap(~Worm)
  # 
  # ggsave(paste0(source_data_name,"_normalised_intensity_plot.png"))
  
  # detect peaks using pracma package
  peaks <- do.call(rbind,lapply(unique(intensities$Worm),detect_peaks,intensities=intensities))
  
  write_xlsx(peaks,paste0("output_data/",source_data_name,"_",sheet,"_peaks_raw.xlsx"))
  
  peaks_summary <- peaks %>% 
    group_by(worm) %>% 
    summarise(n_peaks = n(),mean_peak_intensity = mean(peak_height),mean_peak_auc = mean(auc),mean_peak_median_ratio = mean(median_ratio))
  
  peaks_summary$puncta_per_um <- do.call(c,lapply(unique(intensities$Worm),puncta_per_um,intensities=intensities,peaks=peaks))

  write_xlsx(peaks_summary,paste0("output_data/",source_data_name,"_",sheet,"_peaks_summary.xlsx"))
  
  
  plot_list <- lapply(unique(intensities$Worm),plotting_func,intensities=intensities,peaks=peaks)
  
  wrap_plots(plot_list,nrow=5, ncol=4)
  
    ggsave(paste0("plots/all_peak_plots_",source_data_name,"_",sheet,"_4000cutoff_24h.pdf"),width = 6500, height = 3000, units = "px")
}

filelist_all <- list.files(pattern = paste(".*_snb1.xlsx",sep="|"))

lapply(filelist_all,analyse_intensities,sheet="axon")
