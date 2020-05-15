
library(shiny)

# Define UI for slider demo application
shinyUI(pageWithSidebar(

  #  Application title
  headerPanel("Random Testing"),

  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
  
    sliderInput("pop_size", h3("Population Size"),
                min = 5000, max = 50000, value = 50000, step=1000),
    sliderInput("tests", h3("Daily tests"),
                min = 200, max = 2000, value = 1000, step=100),  
    sliderInput("ppn_sympt", h3("Proportion of cases symptomatic"),
                min = 0, max = .5, value = .05, step=.05),  
    sliderInput("care_seeking_delay", h3("Days to seeking care"),
                min = 0, max = 5, value = 2, round=TRUE),
    sliderInput("R0", h3("R0"),
                min = 1.5, max = 4, value = 2, step=0.25),  
    sliderInput("ppn_immune", h3("Proportion Immune"),
                min = .5, max = 1, value = .85, step=0.05),
    sliderInput("sensitivity", h3("Testing sensitivity"),
                min = .9, max = 1, value = .99, step=0.05),
    sliderInput("specificity", h3("Testing specificity"),
                min = .9, max = 1, value = .99, step=0.05)
    
	),
  # plot output
  mainPanel(
    tabsetPanel(
		tabPanel("Time to Detection",plotOutput("Time_to_detection_plot")),
		tabPanel("False Positives",plotOutput("False_positives_plot"))
		)
	)	
))