library(bs4Dash)
library(shinycssloaders)
library(thematic)
library(dplyr)
library(shinyWidgets)
library(plotly)

source("global.R")

thematic_shiny()

shinyOptions(plot.autocolors = TRUE)
print(getwd())

title <- bs4Dash::dashboardBrand(title = "Embeddings")

header <-
  bs4Dash::dashboardHeader(title = title)


sidebar <- bs4Dash::dashboardSidebar(collapsed = T,
                                     bs4Dash::sidebarMenu(
                                       id = "tabs",
                                       bs4Dash::menuItem(
                                         "UMAP",
                                         icon = icon("dashboard"),
                                         tabName = "umap",
                                         selected = T,
                                         badgeColor = "primary"
                                       )
                                     ))# bookmarkButton(class = "btn-primary"))

body <- bs4Dash::dashboardBody(bs4Dash::tabItems(bs4Dash::tabItem(
  "umap",
  fluidRow(
    column(
      3,
      bs4Dash::box(
        title = "UMAP/PCA Parameters",
        solidHeader = T,
        status = "primary",
        height = "450px",
        width = 12,
        fluidRow(
          column(
            6,
            numericInput(
              "param_umap_nneighbours",
              "# Neighbours",
              min = 2,
              value = 5,
              max = 100
            ),
            selectInput("param_umap_metric",
                        "Metric",
                        choices = umap_distance_metrics)
          ),
          column(
            6,
            numericInput(
              "param_umap_mindist",
              "Min. Dist",
              min = 0.0,
              value = 0.1,
              max = 0.99
            ),
            numericInput(
              "param_umap_epochs",
              "#Epochs",
              min = 1,
              value = 100,
              max = 500
            )
          ),
          column(12,
             sliderInput(
               "n_pca_dim",
               min = 2,
               max = 20,
               value = 3,
               label = "# PCA dimensions"
             )
          )
        ),
        fluidRow(
          fileInput(
            "raw_data_upload",
            "Choose RDS File",
            width = "100%",
            multiple = FALSE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain,.Rds",
              ".csv"
            )
          )
        ),
        fluidRow(
          actionButton(
            "embedd_model",
            "Embedd UMAP Model",
            status = "info",
            width = "100%"
          )
        )
        
      )
    ),
    column(
      4,
      bs4Dash::box(
        title = "Explained Variance (PCA)",
        solidHeader = T,
        height = "450px",
        status = "primary",
        width = 12,
        plotOutput("pca_graph", height = "430px") %>% withSpinner(type = 8)
        # maybe add bitplot ? https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz_pca
      )
    ),
    column(
      5,
      bs4Dash::box(
        height = "450px",
        title = "Selected Data",
        solidHeader = T,
        status = "primary",
        width = 12,
        fluidRow(column(
          6,
          bs4Dash::box(
            width = 12,
            title = "Points in Cluster A",
            DT::dataTableOutput("pointselectionA") %>% withSpinner(type = 8)
          )
        ),
        column(
          6,
          bs4Dash::box(
            width = 12,
            title = "Points in Cluster B",
            DT::dataTableOutput("pointselectionB") %>% withSpinner(type = 8)
          )
        ))
      )
    )
    #,
    # column(
    #   5,
    #   bs4Dash::tabBox(
    #     width = 12,
    #     title = "",
    #     height = "400px",
    #     status = "primary",
    #     solidHeader = T,
    #     id = "raw_data_eda",
    #     tabPanel(
    #       "Linechart",
    #       plotlyOutput("line_graph", height = "380px") %>% withSpinner(type = 8)
    #     ),
    #     tabPanel(
    #       "Correlation Matrix",
    #       plotOutput("correlation_graph") %>% withSpinner(type = 8)
    #     ),
    #     tabPanel(
    #       "Histogram",
    #       plotlyOutput("histogram_graph", height = "380px") %>% withSpinner(type = 8)
    #     )
    #   )
    # )
  ),
  
  fluidRow(column(
    7,
    bs4Dash::box(
      title = "Embedding",
      solidHeader = T,
      height = "1000px",
      status = "primary",
      width = 12,
      maximizable = T,
      fluidRow(
        column(
          3,
          selectInput(
            "select_points_for",
            label = "Select Point(s) for Cluster",
            choices = c("Cluster A", "Cluster B")
          )
        ),
        column(
          3,
          selectInput(
            'param_color_method',
            "Color by (Areas)",
            choice = c("PCA-Scores","Q-Residuals")
          )
        ),
        column(
          3,
          selectInput(
            "param_umap_targetcol",
            "Color by (Points)",
            choices = c("None")
          )
        ),
        column(
          3,
          numericInput(
            "param_loading_threshold",
            "Loadings Threshold",
            value = 0.01,
            min = 0,
            max = 1
          )
        )
      ),
      suppressWarnings({
        uiOutput('voronoi_tabs') %>% withSpinner(type = 8)
      })
    )
  ),
  column(
    5,
    bs4Dash::box(
      height = "1000px",
      title = "Point & Cluster Comparison",
      solidHeader = T,
      status = "primary",
      width = 12,
      fluidRow(column(
        5,
        numericInput(
          "param_contribution_threshold",
          "Feature Threshold",
          value = 0.01,
          min = -10,
          max = 10
        )
      )),
      
      plotlyOutput("explainable_comparison") %>% withSpinner(type = 8)
    )
  ))
)))



bs4Dash::dashboardPage(
  title = "X-MAP GUI",
  header = header,
  sidebar = sidebar,
  body = body,
  dark = T
)
