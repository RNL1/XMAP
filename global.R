# Install & Load packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(shiny, dplyr, pracma, ggExtra, plotly, ggplot2, magrittr, data.table, DT, bs4Dash, shinycssloaders, shinyWidgets, thematic, tictoc, umap, tidyr, factoextra, pracma, ggvoronoi)

# global variables
umap_distance_metrics <- c("euclidean","manhattan","cosine")

includance_color <- "#225EA8" #blueish
excludance_color <- "#E31A1C" #redish

# functions
createVoronoiPlots <-
  function(embedding_data,
           plot_labels,
           method,
           method_vals,
           additionalTitle = "",
           logPlot = FALSE,
           base_size = 18,
           calc_first_only = F) {
    suppressWarnings({
      library(gridExtra)
      library(ggvoronoi)
      colorNeg <- "#253494"
      # colorMid <- "white"
      colorMid <- "#F2F3F4"
      colorPos <- "#722F37"
      #colorNa <- "#f768a1"
      colorNa <- "grey"
      colorPointsDefault <- "#161109"
      colorPointsLow <- "#9ebcda"
      colorPointsHigh <- "#6e016b"
      
      plots <- list()
      for (idx in 1:min(ncol(method_vals), 25)) {
        tictoc::tic(msg = paste0("Create Voronoi Plot ", idx))
        plot_dat <-
          bind_cols(embedding_data, z = method_vals[, idx, drop = F])
        colnames(plot_dat) <- c("x", "y", "z")
        
        p1 <- ggplot(plot_dat, aes(x = x, y = y))
        
        
        p1 <-
          p1 + geom_voronoi(aes(fill = z), color=NA) #+ ggtitle(paste0("PC", idx, additionalTitle))
        
        
        p1 <- p1 +
        # p1 <- p1 + stat_voronoi(geom = "path") +
        # p1 <- p1 + stat_voronoi() +
           scale_fill_gradient2(
            low = colorNeg,
            high = colorPos,
            mid=colorMid,
            na.value = colorNa
          ) + labs(fill = "value")
        
        if(!is.null(plot_labels)){
          p1 <- p1 + geom_point(aes(colour = plot_labels))
          if(is.numeric(plot_labels)){
            p1 <- p1 + scale_color_gradient(low=colorPointsLow, high=colorPointsHigh)
          }
        }else{
          p1 <- p1 + geom_point(col=colorPointsDefault)
        }
        
        if (logPlot) {
          plot_dat <- plot_dat %>% mutate(z_log = log10(z))
          p2 <- ggplot(plot_dat, aes(x = x, y = y))
          p2 <-
            p2 + geom_voronoi(aes(fill = z_log)) #+ ggtitle(paste0("log10(", "PC", idx, additionalTitle, ")"))
          p2 <- formatPlot(p2)
          g <- grerid.arrange(p1, p2, nrow = 1)
        } else{
          g <- p1
        }
        plots <- c(plots, list(g))
        tictoc::toc()
        if (calc_first_only == T) {
          return(plots)
        }
      }
      return(plots)
    })
  }


calcluateTContributions <- function(sample, no_pc, pca){

  p_lambda_pt <- (pca$rotation[,1:no_pc, drop=F] %*% inv(diag((pca$sdev[1:no_pc, drop=F]))) %*% t(pca$rotation[,1:no_pc, drop=F]))


  
  if(is.vector(sample)){
    sample_mat <- matrix(sample, nrow=1, dimnames=list(NULL, names(sample)))
  }else{
    sample_mat <- sample
  }
  t_val <- sample_mat %*% p_lambda_pt
  return(t_val)
}

### based on https://wiki.eigenvector.com/index.php?title=T-Squared_Q_residuals_and_Contributions
compareRelTContributions <- function(sample1, sample2, no_pc, pca){
  t_val1 <- calcluateTContributions(sample1, no_pc, pca)
  t_val2 <- calcluateTContributions(sample2, no_pc, pca)
  res <- (t_val1 - t_val2)
  res <- structure(as.vector(res), names=colnames(res))
  return(res)
}


calculateHotelingTs <- function(data, no_pc, pca){
  res <- as.vector(apply(data, 1, function(row, no_pc, pca){
    sample <- calcluateTContributions(row, no_pc, pca)
    sample_mat <- matrix(sample, nrow=1)
    if(ncol(sample_mat) == 1){
      return(sample_mat * sample_mat)
    }
    return(sample_mat %*% t(sample_mat))
  }, no_pc = no_pc, pca=pca))
  return(res)
}

