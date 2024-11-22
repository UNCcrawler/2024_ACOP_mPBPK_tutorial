###########################################################
# Project: Shiny app for simulation of pH-dependent antibody
# Author: Dongfen Yuan
# Date: June 11, 2023
###########################################################
library(shiny)
# library(RxODE) # The mdoel was initially built using RxODE, but is compatible with rxode2
library(rxode2)
library(tidyverse)

# mPBPK model with a nested endosome compartment
ode = "
 # mPBPK parameters for a 70 kg human
 Vp = 2.6			# L, plasma volume
 Visf = 15.6			# L, total interstitial fluid volume
 Fisf = 0.8		 	# fraction of interstitial fluid available for IgG1 distribution
 V1 = 10.14		 	# L, interstitial fluid volume of lumped tight tissue
 V2 = 5.46			 # L, interstitial fluid volume of lumped leaky tissue
 Vlymph = 5.2	 		# L, lymph volume
 L = 0.121		 	 # L/h, lymph flow
 L1 = 0.0399		 	# L/h, convection from plasma to lumped tight tissue
 L2 = 0.081	 		# L/h, convection from plasma to lumped leaky tissue
 Rf1 = 0.945		 	# vascular reflection coefficient of lumped tight tissue, average data in PMID 25146360
 Rf2 = 0.697		 	# vascular relection coefficient of lumped leaky tissue, average data in PMID 25146360
 Rflymph = 0.2		 	# lymphatic reflection coefficient
 
 # Endosome parameters
 Ve = 0.0035				# L, endosome volume
 CLup = 0.061667			# L/h, non-specific pinocytosis rate
 CLe = 0.23978			#L/h, manually fitted data for CLup = 0.061667 and CLrec = 0.0189125
 CLrec = 0.693/(8/60)*Ve	# L/h recycling clearance
  
 #k1on = 0.0008676			#1/(pM*h), antibody-FcRn binding on rate
 #k1off = 583.2				#1/h, antibody-FcRn binding off rate
 FcRntotal = 49.8E6			#pM
 
 # Target parameters
 CLp_t = 0.693/Thalf*Vp  # L/h, total target clearance
 CLp_t1 = CLp_t - CLup	# L/h, endosome-independent target clearance
 ksyn = R0*CLp_t  # pmol/h, target synthesis rate

 
 # Total antibody
 Cp_at = Cp_a + Cp_atc

 
 # Total target
  Cp_tt = Cp_t + Cp_atc
  
 # Initial values
 Cp_a(0) = 0
 Cp_t(0) = R0
 Cp_atc(0) = 0
 Ce_a(0) = 0
 Ce_t(0) = 0
 Ce_atc(0) = 0
 FcRnA(0) = 0
 FcRnATC(0) = 0
 FcRn(0) = FcRntotal
 C1_a(0) = 0
 C1_atc(0) = 0
 C2_a(0) = 0
 C2_atc(0) = 0
 Clymph_a(0) = 0
 Clymph_atc(0) = 0
 

 # Plasma
 
 d/dt(Cp_a) = -kon * Cp_a * Cp_t + koff * Cp_atc - (1-Rf1) * L1 * Cp_a/Vp - (1-Rf2) * L2 * Cp_a/Vp  + L * Clymph_a/Vp - CLup*Cp_a/Vp +CLrec*FcRnA/Vp
 
 d/dt(Cp_t) = ksyn/Vp - CLp_t1 * Cp_t/Vp - CLup*Cp_t/Vp - kon * Cp_a * Cp_t + koff * Cp_atc
 
 d/dt(Cp_atc) = kon * Cp_a * Cp_t - koff * Cp_atc - (1-Rf1) * L1 * Cp_atc/Vp - (1-Rf2) * L2 * Cp_atc/Vp + L * Clymph_atc/Vp -CLup*Cp_atc/Vp +CLrec*FcRnATC/Vp
 
 # Endosome
 
 d/dt(Ce_a) = - keon*Ce_a*Ce_t + keoff*Ce_atc - k1on*FcRn*Ce_a + k1off*FcRnA + CLup*Cp_a/Ve - CLe*Ce_a/Ve
 
 d/dt(Ce_t) = - keon*Ce_a*Ce_t + keoff*Ce_atc - keon*FcRnA*Ce_t + keoff*FcRnATC + CLup*Cp_t/Ve - CLe * Ce_t/Ve
 
 d/dt(Ce_atc) = keon*Ce_a*Ce_t - keoff*Ce_atc - k1on*FcRn*Ce_atc + k1off*FcRnATC + CLup*Cp_atc/Ve - CLe*Ce_atc/Ve
 
 d/dt(FcRnA) = k1on*FcRn* Ce_a - k1off*FcRnA - keon*FcRnA*Ce_t + keoff*FcRnATC - CLrec * FcRnA/Ve
 
 d/dt(FcRnATC) = k1on*FcRn*Ce_atc - k1off*FcRnATC + keon*FcRnA*Ce_t - keoff*FcRnATC - CLrec*FcRnATC/Ve 
 
 d/dt(FcRn) = - k1on*FcRn*Ce_a + k1off*FcRnA - k1on*FcRn*Ce_atc + k1off*FcRnATC + CLrec*(FcRnA + FcRnATC)/Ve
 
 
 # Tight tissue
 
 d/dt(C1_a) = (1-Rf1) * L1 * Cp_a / (V1 * Fisf) - (1-Rflymph) * L1 * C1_a / (V1*Fisf) 
 
 d/dt(C1_atc) = (1-Rf1) * L1 * Cp_atc / (V1 * Fisf) -(1-Rflymph) * L1 * C1_atc / (V1 * Fisf)
 
 # Leaky tissue
 
 d/dt(C2_a) = (1-Rf2) * L2 * Cp_a / (V2 * Fisf) - (1-Rflymph) * L2 * C2_a / (V2 * Fisf)
 
 d/dt(C2_atc) = (1-Rf2) * L2 * Cp_atc / (V2 * Fisf) -(1-Rflymph) * L2 * C2_atc / (V2 * Fisf) 
 
 # Lymph
 
 d/dt(Clymph_a) = (1-Rflymph) * L1 * C1_a / Vlymph + (1-Rflymph) * L2 * C2_a / Vlymph - L * Clymph_a / Vlymph
 
 d/dt(Clymph_atc) = (1-Rflymph) * L1 * C1_atc / Vlymph + (1-Rflymph) * L2 * C2_atc / Vlymph - L * Clymph_atc / Vlymph"


mod = RxODE(model = ode)



ui = fluidPage(
  
  # App title ----
  titlePanel("Simulation of pH-dependent antibody in human"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      
      
      # Input: antibody dose ----
      textInput(inputId = "dose",
                label = "Antibody IV dose (mg/kg, multiple doses separated by comma)",
                value = "0.25, 0.5, 1, 3, 5"),
      
      # Input: simulation duration ----
      numericInput(inputId = "duration",
                   label = "Simulation duration (days)",
                   value = 30),
      br(),
      
      
      # Heading
     h4("Antibody-hFcRn binding affinity"),
      
      # Input: K1on ----
      numericInput(inputId = "K1on",
                   label = "K1on (E5/[M*s])",  ##need to be converted to 1/[pM*h] in simulation
                   value = 0.0008676/60/60*10^7), 
      
      # Input: K1off ----
      numericInput(inputId = "K1off",
                   label = "K1off (1/s):",
                   value = 583.2/60/60),  ##need to be converted to 1/h in simulation
      
      
      br(),
      # Heading
      h4("Soluble target kinetics"),
      
      # Input: soluble target baseline ----
      numericInput(inputId = "R0",
                   label = "Baseline (pM):",
                   value = 0.276),
      
      # Input: soluble target half-life
      numericInput(inputId = "Thalf",
                   label = "Soluble target Half-life (min):",
                   value = 13),
      
      br(),
      # Heading
      h4("Antibody-soluble target binding affinity in plasma"),
      
      # Input: Kon ----
      numericInput(inputId = "Kon",
                   label = "Kon (E5/[M*s])", ##need to be converted to 1/[pM*h] in simulation
                   value = 0.002592/60/60*10^7),
      
      # Input: Koff ----
      numericInput(inputId = "Koff",
                   label = "Koff (E-4/s):", 
                   value = 0.468/60/60*10^4),  ##need to be converted to 1/h in simulation
      
      br(),
      # Heading
      h4("Antibody-soluble target binding affinity in acidic endosome"),
      
      # Input: Keon ----
      numericInput(inputId = "Keon",
                   label = "Keon (E5/[M*s]):", ##need to be converted to 1/[pM*h] in simulation
                   value = 0.002592/60/60*10^7),
      
      # Input: Keoff ----
      numericInput(inputId = "Keoff",
                   label = "Keoff (E-4/s):",
                   value = 0.468/60/60*10^4),  ##need to be converted to 1/h in simulation
      
      br(),
      # Heading
      h4("Upload observed data"),
      
      # Input: select observed data ----
      fileInput(inputId = "observed",
                label = "Choose .csv file",
                accept = ".csv")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(
        tabPanel(title = "Simulation",
                 # Heading:
                 h4("Antibody plasma PK"),
                 
                 # Output: Antibody PK ----
                 plotOutput(outputId = "PKPlot"),
                 
                 # Heading:
                 h4("Soluble target plasma PK"),
                 
                 # Output: Soluble target PK ----
                 plotOutput(outputId = "PDPlot")
                 
        ),
        
        tabPanel(title = "About",
                 
                 # Title
                 fixedRow(
                   column(10,
                          h4("A Minimal Physiologically Based Pharmacokinetic Model with a Nested Endosome Compartment for Novel Engineered Antibodies", align = "center"), offset = 0),
                   column(10,
                          h4("Dongfen Yuan, Frederik Rode, Yanguang Cao", align = "center"), offset = 0),
                   column(10,
                          h5("AAPS J. 2018 Mar 14;20(3):48", align = "center"), offset = 0)
                 ),
                 
                 # Section: introduction
                 h4("Introduction"),
                 p("Over the years, several people reached out to me asking questions on this publication. While I am greatful for the interest, I regret that I did not include the model code in the publication. "),
                 p("I took one day and developed a R shiny app. The default setting can be used to reproduce Figure 2b, adalimumab PK in a 70 kg human (digitized data were provided in the data folder). The users can also use it to simulate the PK and PD of pH-dependent antibodies against soluble targets by providing the target baseline/half-life and antibody binding affinities with target in plasma/acidic endosomes."),
                 p("Please refer to the paper for details on model assumptions, methods, and applications. Unfortunately, the model equations for antibody and antibody-target complex in plasma in the paper dropped the return of antibody and antibody-target complex from lymph compartment by mistake. Please refer to the model equations in the R-script."),
                 p("Given that most questions that I received are from industry, I chose to upload the source code to Github so that you can download and adapt to your data without risking uploading sensitive data to an external server."),
                 p("If you find errors or have additional questions, please contact Dongfen Yuan via email dongfen.yuan@gmail.com."),
                 p("Dongfen on 11 June 2023"),
                 
                 # Section: model structure
                 h4("Model structure"),
                 
                 img(src = "Model_structure.png", height = 280, width = 700),
                 
                 # Section: version update
                 h4("Version update"),
                 p("V1-June 11, 2023"),
                 p("V2-November 8, 2024: Include K1on and K1off as input boxes and use more convenient units for all binding rate constants.")
                 
        )
      )
      
      
      
    )
  )
)




server = function(input, output){
  # RxODE simulation
  
  sim = reactive({
    theta = c(Rf1 = input$Rf1, Rf2 = input$Rf2, # antibody vascular reflection coefficient for tight and leaky tissues
              k1on = input$K1on*60*60*10^(-7), k1off = input$K1off*60*60, # antibody-hFcRn binding affinity
              kon = input$Kon*60*60*10^(-7), koff = input$Koff*60*60*10^(-4), # antibody-target binding affinity in plasma
              keon = input$Keon*60*60*10^(-7), keoff = input$Keoff*60*60*10^(-4), # antibody-target binding affinity in acidic endosome
              Thalf = input$Thalf/60, R0 = input$R0)
    
    map_df(strsplit(input$dose, split = ",")%>%unlist%>%as.numeric%>%sort(), 
           function(x){
             ev = eventTable()%>%
               add.dosing(dose = x*70/148000*1000000000/2.6)%>%  # assume mAb 148kDa, mg to pmol/L
               add.sampling(time.units = "h", time = 0:(input$duration*24))
             
             sim = solve(mod, theta, ev)%>%
               select(time, Cp_a, Cp_at, Cp_t, Cp_tt)%>%
               mutate(dose = x)
           }) ->simdat
    
    simdat = simdat%>%
      pivot_longer(cols = starts_with("Cp"))%>%
      mutate(time = time/24, # h to day
             value = if_else(name %in% c("Cp_a", "Cp_at"), 
                             value/1000, # pM to nM for antibody
                             value))%>%
      mutate(dose = factor(dose, levels = unique(dose), 
                           labels = paste(unique(dose), "mg/kg")))
    
    return(simdat)
    
  })
  
  obs = reactive({
    if (!is.null(input$observed)){
      
      read.csv(input$observed$datapath, stringsAsFactors = F)%>%
        pivot_longer(starts_with("Cp"))%>%
        mutate(dose = paste(dose, "mg/kg"))
    }
  })
  
  labelname = list(
    Cp_a  = "Free antibody",
    Cp_at = "Total antibody",
    Cp_t = "Free target",
    Cp_tt = "Total target"
  )
  
  # Output: Antibody PK ----
  output$PKPlot = renderPlot({
    
    ggplot()+
      geom_line(data = filter(sim(), name %in% c("Cp_a", "Cp_at")),
                aes(x = time, y = value, color = dose),
                linewidth = 1)+
      facet_wrap(~name, 
                 labeller = function(variable, value){return(labelname[value])})+
      labs(x = "Days", y = "Plasma conc (nM)")+
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
        )+
      theme_bw()+
      theme(strip.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))+
      annotation_logticks(sides = "l")->p
    
    if(!is.null(obs())){
      if(!is.null(filter(obs(), name %in% c("Cp_a", "Cp_at")))){
        p+geom_point(data = filter(obs(), name %in% c("Cp_a", "Cp_at")),
                     aes(x = time/24, y = value, color = dose),
                     size = 3)
      } else{
        p}
    } else {
      p
    }
    
    
  })
  
  # Output: target PK ----
  output$PDPlot = renderPlot({
    
    ggplot()+
      geom_line(data = filter(sim(), name %in% c("Cp_t", "Cp_tt")),
                aes(x = time, y = value, color = dose),
                size = 1)+
      facet_wrap(~name, 
                 labeller = function(variable, value){return(labelname[value])})+
      labs(x = "Days", y = "Plasma conc (pM)")+
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )+
      theme_bw()+
      theme(strip.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))+
      annotation_logticks(sides = "l")   -> p
    
    if(!is.null(obs())){
      
      if(!is.null(filter(obs(), name %in% c("Cp_t", "Cp_tt")))){
        p+geom_point(data = filter(obs(), name %in% c("Cp_t", "Cp_tt")),
                     aes(x = time/24, y = value, color = dose),
                     size = 3)
      } else {
        p
      }
    } else {
      p
    }
    
  })
  
}

shinyApp(ui = ui, server = server)

# img(src = "rstudio.png", height = 140, width = 400)

## Below is the package information for reproducibility
# sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8 
# [2] LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils    
# [5] datasets  methods   base     
# 
# other attached packages:
#   [1] lubridate_1.9.3 forcats_1.0.0  
# [3] stringr_1.5.1   dplyr_1.1.4    
# [5] purrr_1.0.2     readr_2.1.5    
# [7] tidyr_1.3.1     tibble_3.2.1   
# [9] ggplot2_3.5.1   tidyverse_2.0.0
# [11] rxode2_3.0.1    shiny_1.9.1    
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.5        xfun_0.48          
# [3] bslib_0.8.0         RApiSerialize_0.1.4
# [5] lattice_0.22-6      tzdb_0.4.0         
# [7] vctrs_0.6.5         tools_4.4.1        
# [9] generics_0.1.3      fansi_1.0.6        
# [11] pkgconfig_2.0.3     data.table_1.16.2  
# [13] checkmate_2.3.2     RcppParallel_5.1.9 
# [15] lifecycle_1.0.4     farver_2.1.2       
# [17] compiler_4.4.1      textshaping_0.4.0  
# [19] munsell_0.5.1       qs_0.27.2          
# [21] httpuv_1.6.15       htmltools_0.5.8.1  
# [23] sys_3.4.3           sass_0.4.9         
# [25] later_1.3.2         pillar_1.9.0       
# [27] crayon_1.5.3        jquerylib_0.1.4    
# [29] cachem_1.1.0        nlme_3.1-164       
# [31] mime_0.12           tidyselect_1.2.1   
# [33] digest_0.6.37       lotri_1.0.0        
# [35] stringi_1.8.4       labeling_0.4.3     
# [37] rxode2ll_2.0.11     fastmap_1.2.0      
# [39] grid_4.4.1          colorspace_2.1-1   
# [41] cli_3.6.3           dparser_1.3.1-12   
# [43] magrittr_2.0.3      utf8_1.2.4         
# [45] withr_3.0.1         scales_1.3.0       
# [47] promises_1.3.0      backports_1.5.0    
# [49] timechange_0.3.0    ragg_1.3.3         
# [51] hms_1.1.3           stringfish_0.16.0  
# [53] memoise_2.0.1       knitr_1.48         
# [55] PreciseSums_0.7     rlang_1.1.4        
# [57] Rcpp_1.0.13         xtable_1.8-4       
# [59] glue_1.8.0          rstudioapi_0.16.0  
# [61] jsonlite_1.8.9      R6_2.5.1           
# [63] systemfonts_1.1.0 

