enableBookmarking(store = "url")

# global_base_size = 40
global_text_size = 18
# global_text_size = 18
global_plot_title_size = 20

shinyServer(function(input, output, session) {
  data_basis <-
    reactiveValues(
      raw_data = NULL,
      # raw data
      raw_data_pivot = NULL,
      # raw data pivoted
      embedded_model = NULL,
      # embedded umap model
      pca_model = NULL,
      # pca model
      no_pcs = NULL,
      # number of pcas
      no_sig_pcs = NULL,
      #  number of significant
      q_vals = NULL,
      # q-val method values
      hotelling_t_vals = NULL,
      # pca scores method values,
      pca_scores = NULL,
      # hotelling method values
      umap_projection_plot = NULL,
      # umap plots
      point_selection_A = NULL,
      # points in selection of cluster A
      point_selection_B = NULL # points in selection of cluster B
    )
  
  observeEvent(input$embedd_model, {
    # return if no data is uploaded
    if (is.null(data_basis$raw_data_pivot)) {
      showNotification(
        "[ERROR] Please upload the raw data or prebuilt embeddings! (icon on top right corner)",
        type = "error"
      )
      return(NULL)
    }
    
    showNotification(
      "[INFO] Started UMAP embedding - Previously upladed embeddings are recalculated with the custom config!",
      type = "message"
    )
    
    withProgress(message = 'Generating UMAP Embedding', value = 0, {
      incProgress(0.2, detail = "Compute PCA")
      
      target_col <- input$target_feature
      
      if (target_col == "None") {
        target_col <- NULL
      }
      
      tictoc::tic("PCA Calculation")
      subset<- data_basis$raw_data_pivot %>% select(-c(target_col)) %>% as.data.frame()
      find.numeric <- sapply(subset, is.numeric)
      subset <- subset[,find.numeric]
      
      pca_input_dat <-
        apply(as.matrix(subset),
        2,
        scale,
        scale = F)
      res_pca <-
        prcomp(pca_input_dat, center = FALSE, scale = FALSE)
      rel_sd <- cumsum(res_pca$sdev ^ 2) / sum(res_pca$sdev ^ 2)
      #umap_data <- res_pca$x[, get_eig(res_pca)$cumulative.variance.percent < 99]
      umap_data <- res_pca$x[, 1:input$n_pca_dim]
      
      data_basis$no_sig_pcs <- NCOL(umap_data)
      
      tictoc::toc()
      
      incProgress(0.1, detail = "Initialize UMAP Parameters")
      
      custom.config = umap.defaults
      custom.config$n_components = 2#input$dims #2
      custom.config$min_dist = input$param_umap_mindist #input$mindist #0.1
      custom.config$n_neighbors = input$param_umap_nneighbours #input$neighbors #15
      custom.config$n_epochs = input$param_umap_epochs
      custom.config$metric = input$param_umap_metric ##input$distance #euclidean
      custom.config$random_state = 42
      
      incProgress(0.2, detail = "Calculate UMAP Embedding")
      
      tictoc::tic("UMAP Embedding")
      
      
      embedding <-
        umap(umap_data,
             custom.config)
      
      tictoc::toc()
      
      incProgress(0.1, detail = "Calculating Q Vals")
      
      tictoc::tic("Q-Vals")
      
      max_comp <- 25
      no_pcs <- min(max_comp, ncol(umap_data)) #res_pca$x
      q_vals <- lapply(1:no_pcs, function(pc) {
        inner <-
          (diag(ncol(res_pca$rotation)) - res_pca$rotation[, pc] %*% t(res_pca$rotation[, pc]))
        qis <- do.call(c, lapply(1:nrow(res_pca$x), function(i) {
          # q_val <- res_pca$x[i,,drop=F] %*% inner %*% t(res_pca$x[i,,drop=F])
          q_val <-
            pca_input_dat[i, , drop = F] %*% inner %*% t(pca_input_dat[i, , drop =
                                                                         F])
          return(q_val)
        }))
        qis <- (qis - min(qis))
        qis <- (qis / max(qis))
        # qis <- abs(qis) * -1
        return(qis)
      })
      names(q_vals) <- sprintf("PC %02d", 1:no_pcs)
      q_vals <- do.call(bind_cols, q_vals)
      
      tictoc::toc()
      
      incProgress(0.1, detail = "Calculating Hotteling's T2")
      
      tictoc::tic("Hotellings T2")
      
      hotelling_t_vals <-
        do.call(cbind, lapply(2:no_pcs, function(no_pc_tmp) {
          print(paste0("Calc Hotellings T2 for PC ", no_pc_tmp))
          hotelingTs <-
            calculateHotelingTs(pca_input_dat, no_pc_tmp, res_pca)
          
          return(hotelingTs)
        }))
      
      colnames(hotelling_t_vals) <- sprintf("#PC %02d", 2:no_pcs)
      
      tictoc::toc()
      
      incProgress(0.1, detail = "Updating UMAP Embedding & Q-Vals & Hotelling T Internally")
      
      data_basis$embedded_model <- embedding
      data_basis$pca_model <- res_pca
      data_basis$q_vals <- q_vals
      data_basis$no_pcs <- no_pcs
      data_basis$hotelling_t_vals <- hotelling_t_vals
      data_basis$pca_scores <- umap_data
      showNotification("[INFO] Finished UMAP embedding", type = "message")
      
    })
  })
  
  observe({
    if (input$select_points_for == "Cluster A") {
      if(!is.null(event_data("plotly_selected")) & length(event_data("plotly_selected"))>0){
        data_basis$point_selection_A = event_data("plotly_selected") %>% select(pointNumber)
      }
    } else if (input$select_points_for == "Cluster B") {
      if(!is.null(event_data("plotly_selected"))  & length(event_data("plotly_selected"))>0){
        data_basis$point_selection_B = event_data("plotly_selected") %>% select(pointNumber)
      }
    }
  })
  
  
  output$pointselectionA <- renderDataTable(data_basis$point_selection_A,extensions = "FixedHeader",
                    style = "bootstrap4",
                    options = list(
                      pageLength = 3,
                      autoWidth = TRUE,
                      paging = TRUE,
                      searching = F,
                      ordering = TRUE,
                      fixedHeader = TRUE
                    ))

  output$pointselectionB <- renderDataTable(data_basis$point_selection_B,extensions = "FixedHeader",
                    style = "bootstrap4",
                    options = list(
                      pageLength = 3,
                      autoWidth = TRUE,
                      paging = TRUE,
                      searching = F,
                      ordering = TRUE,
                      fixedHeader = TRUE
                    ))

  
  
  observeEvent(input$raw_data_upload, {
    file <- input$raw_data_upload
    tictoc::tic("Data Upload")
    raw_data <- fread(file$datapath)
    data_basis$raw_data_pivot <- raw_data
    tictoc::toc()
    
    updateSliderInput(session,
                      "n_pca_dim",
                      max = NCOL(data_basis$raw_data_pivot))
    updateSelectInput(session,
                      "param_umap_targetcol",
                      choices = c("None", colnames(data_basis$raw_data_pivot)))
    updateSelectInput(session,
                      "target_feature",
                      choices = c("None", colnames(data_basis$raw_data_pivot)))
  })
  
  output$explainable_comparison <- renderPlotly({
    if (!is.null(data_basis$point_selection_A) &&
        !is.null(data_basis$point_selection_B)) {
      tictoc::tic("Explainable Comparison")
      # check if all same
      if ((nrow(
        setdiff(
          data_basis$point_selection_A,
          data_basis$point_selection_B
        )
      ) == 0 &
      nrow(
        setdiff(
          data_basis$point_selection_A,
          data_basis$point_selection_B
        )
      ) == 0) == T) {
        return(NULL)
      }
      
      target_col <- input$param_umap_targetcol
      
      if (target_col == "None") {
        target_col <- NULL
      }
      
      subset<- data_basis$raw_data_pivot %>% select(-c(target_col)) %>% as.data.frame()
      find.numeric <- sapply(subset, is.numeric)
      subset <- subset[,find.numeric]
      
      pca_input_dat <-
        apply(as.matrix(subset),
              2,
              scale,
              scale = F)
      
      calcClusterContribution <-
        function(sample_idxs,
                 pca_input_dat,
                 no_pcs,
                 pca_model) {
          as.data.frame(do.call(rbind, (lapply(sample_idxs, function(sample_idx) {
            calcluateTContributions(pca_input_dat[sample_idx, ], no_pcs, pca_model)
          }))))
        }
      
      t_vals_a <-
        calcClusterContribution(
          data_basis$point_selection_A$pointNumber,
          pca_input_dat,
          data_basis$no_sig_pcs,
          data_basis$pca_model
        )
      t_vals_b <-
        calcClusterContribution(
          data_basis$point_selection_B$pointNumber,
          pca_input_dat,
          data_basis$no_sig_pcs,
          data_basis$pca_model
        )
      
      t_vals_a_mean <- colMeans(t_vals_a)
      t_vals_b_mean <- colMeans(t_vals_b)
      
      rel_contr_cluster <- t_vals_a_mean - t_vals_b_mean
      
      plotting_data <- as.data.frame(rel_contr_cluster)
      plotting_data$Feature <- names(rel_contr_cluster)
      colnames(plotting_data) <- c("Rel_Contribution", "Feature")
      
      plotting_data <-
        plotting_data %>% dplyr::filter(
          Rel_Contribution > input$param_contribution_threshold |
            Rel_Contribution < -input$param_contribution_threshold
        )
      
      plt_obj <-
        ggplot(plotting_data, aes(y = Feature, x = Rel_Contribution)) + geom_bar(
          stat = "identity",
          fill = ifelse(
            plotting_data$Rel_Contribution > 0,
            includance_color,
            excludance_color
          )
        ) + scale_x_continuous(limits = c(
          min(plotting_data$Rel_Contribution),
          max(plotting_data$Rel_Contribution)
        )) + ggtitle("Rel. Hotellings T2 Contributions btw. Clusters A and B") + geom_vline(xintercept = 0, alpha = 0.6) +
        theme(
          text = element_text(size = global_text_size),
          plot.title = element_text(size = global_plot_title_size)
        ) +
        # theme_dark(base_size=global_base_size) +
        xlab("Contribution")
      tictoc::toc()
      ggplotly(plt_obj, height=600)
    }
  })
  
  output$pca_graph <- renderPlot({
    if (is.null(data_basis$pca_model)) {
      print("PCA Graph called - no data")
      return (NULL)
    } else{
      tictoc::tic("PCA Graph")
      
      significant_dimensions <-
        get_eig(data_basis$pca_model)$cumulative.variance.percent < 99
      
      
      if (length(significant_dimensions) <= 20) {
        n_plot_dims <- length(significant_dimensions)
      } else{
        n_plot_dims <- 20
      }
      
      fig <-
        fviz_eig(
          data_basis$pca_model,
          ggtheme = NULL,
          ncp = n_plot_dims,
          addlabels = T
        ) +
        # theme_dark(base_size=global_base_size)
        theme(
          text = element_text(size = global_text_size),
          plot.title = element_text(size = global_plot_title_size)
        )
      
      
      fig$layers[[1]]$aes_params$fill <-
        ifelse(significant_dimensions[1:n_plot_dims] == TRUE,
               yes = includance_color,
               no = excludance_color)
      
      tictoc::toc()
      
      return (fig)
    }
  })
  
  output$voronoi_tabs <- renderUI({
    if (!is.null(data_basis$embedded_model)) {
      tictoc::tic("Voronoi Tab Generation")
      
      coloring_values <- data_basis$q_vals
      if (input$param_color_method == "PCA-Scores") {
        coloring_values <- data_basis$pca_scores
      }
      
      target_col <- input$param_umap_targetcol
      embedding_data <-
        as.data.frame(data_basis$embedded_model$layout)
      
      if (target_col %in% colnames(data_basis$raw_data_pivot)) {
        plot_labels <-
          data_basis$raw_data_pivot[[target_col]]
        if(!is.numeric(plot_labels)){
          plot_labels <- as.factor(plot_labels)
        }
      } else{
        plot_labels <- NULL
      }
      
      voronoi_plots <-
        createVoronoiPlots(
          embedding_data,
          plot_labels,
          input$param_color_method,
          coloring_values,
          calc_first_only = F
        )
      nTabs <- length(voronoi_plots)
      
      loading_plots <- lapply(1:nTabs, function(tab) {
        plotting_data <- as.data.frame(data_basis$pca_model$rotation[, tab])
        plotting_data$Feature <-
          names(data_basis$pca_model$rotation[, tab])
        colnames(plotting_data) <- c("Loading", "Feature")
        
        plotting_data <-
          plotting_data %>% dplyr::filter(
            Loading > input$param_loading_threshold  |
              Loading < -input$param_loading_threshold
          )
        #param_loading_threshold
        plt_obj <-
          ggplot(plotting_data,
                 aes(
                   y = Feature,
                   x = Loading,
                   fill = Loading
                 )) +
          geom_bar(
            stat = "identity",
            fill = ifelse(
              plotting_data$Loading > 0,
              includance_color,
              excludance_color
            )
          ) +
          scale_x_continuous(limits = c(-1, 1)) + ggtitle(paste0("Loadings for PC", tab)) +
          geom_vline(xintercept = 0, alpha = 0.6) +
          theme(
            text = element_text(size = global_text_size),
            plot.title = element_text(size = global_plot_title_size)
          )
        
      })
      
      voronoiTabs = lapply(1:nTabs, function(x) {
        tabPanel(paste0("PC", x),
              fluidPage(
                 fluidRow(column(
                   6,
                   renderPlotly(
                     ggplotly(
                       voronoi_plots[[x]] + theme(
                         text = element_text(size = global_text_size),
                         plot.title = element_text(size = global_plot_title_size)
                       ),
                       height = 600,
                       #width = 700
                     ) %>% layout(dragmode = "select")
                   )
                 ),
                 column(6,
                        renderPlotly(
                          ggplotly(loading_plots[[x]], height=600)# height=780
                        )))))
      })
      tictoc::toc()
      do.call(tabsetPanel, voronoiTabs)
    }
  })
})
