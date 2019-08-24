#Raymond Atta-Fynn
#NYC Data Academy
#Routine for building a model for water molecules for atomistic modeling
#August 2, 2019

library(shiny)
library(rsconnect)
library(ggplot2)
source("shiny_generate_initial_configuration.R")
source("shiny_visualize.R")
source("shiny_rdf.R")
source("shiny_anint.R")
source("shiny_min_vec_dist.R")
source("shiny_atomic_mass.R")
source("shiny_ion_water_12_6_4_parameters.R")

ui <- fluidPage(
  titlePanel(HTML('<p style="color:purple; font-size: 18pt; text-align:center">
  Molecular modeling of water-based systems<br/>\nRaymond Atta-Fynn, 
                  NYC Data Science Academy<br/>\nAugust 3, 2019</p>')),
  tags$br(style="line-height: 6px"),
  tags$style("#selected_var {font-size:18pt; color:red; display:block; text-align:center}"),
  tags$style("#min_max {font-size:18pt; color:blue; display:block; text-align:center}"),
  tags$style("#label {font-size:14pt; color:red; display:block; }"),
  tags$style("#NW {font-size:14pt; color:blue; display:block; }"),
  tags$style("#PERIODIC {font-size:14pt; color:blue; display:block; }"),
  tags$style("#md_tstep {font-size:14pt; color:blue; display:block; }"),
  tags$style("#md_steps {font-size:14pt; color:blue; display:block; }"),
  tags$style("#md_type {font-size:14pt; color:blue; display:block; }"),
  tags$style("#md_temp {font-size:14pt; color:blue; display:block; }"),
  tags$style("#geom_steps {font-size:14pt; color:blue; display:block; }"),
  sidebarLayout(
    sidebarPanel(width=3,
#      helpText(HTML('<p style="color:black; font-size: 18pt">Options:<br/>\n Pure water 
#                    <br/>\ncharged metal ion in water<br/>\n\n</p>')),
#      tags$br(),
      
      tags$style(type='text/css', ".selectize-input { font-size: 14pt; line-height: 14pt;} 
                 .selectize-dropdown { font-size: 14pt; line-height: 14pt; }"),
          selectInput("var", 
                  label = div(style = "font-size: 14pt; color: red",  "Select System (25 choices available)"),
                  choices =c("pure water","Al3+", "Fe3+", "Cr3+", "In3+", "Tl3+", "Y3+", 
                             "La3+", "Ce3+", "Pr3+", "Nd3+", "Sm3+", "Eu3+", "Gd3+", "Tb3+", 
                             "Dy3+", "Er3+", "Tm3+", "Lu3+", "Hf4+", "Zr4+", "Ce4+", "U4+", "Pu4+", "Th4+"),
                  selected = "Al3+"),
      #tags$br(),
      
      numericInput("NW", HTML('<p style="color:red; font-size: 14pt">Number of 
                              water molecules<br/>\n(no more than 256)</p>'), 
                   20, min = 1, max = 256,step=1),
      #tags$br(),

      selectInput("water_model", 
                  label = div(style = "font-size: 14pt; color: red",  "Water Model (3 Choices)"),
                  choices =c("SPCE","TIP3P","MD"),
                  selected = "SPCE"), 
      
      selectInput("PERIODIC", 
                  label = div(style = "font-size: 14pt; color: red",  "Periodic Boundary Condition"),
                  choices =c("Yes","No"),
                  selected = "Yes"), 
      #tags$br(),

      selectInput("run_type", 
                  label = div(style = "font-size: 14pt; color: red",  "Type of run"),
                  choices =c("Only initial structure","Geom. Opt.","MD","MD + Geom. Opt."),
                  selected = "Only initial structure"), 
      #tags$br(),
      
      numericInput("geom_steps", HTML('<p style="color:red; font-size: 14pt">Geom. Opt. Steps</p>'),25),
      #tags$br(),
      
      selectInput("md_type", 
                  label = div(style = "font-size: 14pt; color: red",  "MD Ensemble"),
                  choices =c("NVE","NVT"),
                  selected = "NVE"), 
      #tags$br(),
      
      numericInput("md_steps", HTML('<p style="color:red; font-size: 14pt">Total number of MD steps</p>'),25),
      #tags$br(),
      
      numericInput("md_tstep", HTML('<p style="color:red; font-size: 14pt">MD Time Step in fs</p>'),1.0),
      #tags$br(),
      
      numericInput("md_temp", HTML('<p style="color:red; font-size: 14pt">MD Temperature in K</p>'),300.0),
      #tags$br(),
      
      #submitButton("submit", icon("refresh")) 
      actionButton("SubmitButton", "submit job",style='padding:4px; font-size:14pt'), 
      actionButton("StopButton", "stop job",style='padding:4px; font-size:14pt')
      
    ),
    
    mainPanel(
      textOutput("selected_var"),
      textOutput("min_max"),
      tags$br(),
      tags$br(),
      fluidRow(
        column(5, rglwidgetOutput("myplot",  width = 440, height = 450)),
        column(7, plotOutput("myrdf",width = 570, height = 500))
        )
    )
  )
)

server <- function(input, output) {
  
  output$selected_var <- renderText({ 
    if(input$var=='pure water')paste("System:", input$var)
    if(input$var!='pure water')paste("System:", input$var,' ion in water')
  })
  
  output$min_max <- renderText({ 
    paste("Number of water molecules:",
          input$NW)
  })
  
  ## I would like to evaluate the following block to get the data frame pos
  
  #N <- input$NW
  #system_type <- input$var
  #periodic = TRUE
  #if(input$PERIODIC=="No"){
  #  periodic=FALSE
  #}
  ## Assign the output to the variable tmp
  #tmp <- generate_initial_configuration(N,system_type,periodic)
  
  ## L is the box length; it is required by the RGL visualizer below
  #L=tmp[[1]]
  
  ## pos is a data frame of the atomic positions and atomic symbols
  ## It is also required by the RGL visualizer below
  #pos <- tmp[[2]]
  
  
#  if(input$run_type=="Geom. Opt."){
#    for (i in 1:input$geom_steps){
#      geometry_optimization(pos,input$PERIODIC)
#    }
#  }
  
#  if(input$run_type=="MD"){
#    if(input$md_type=="NVE"){
#      for (i in 1:input$md_steps){
#        md_nve(pos,input$md_tstep,input$PERIODIC)
#      }
#    } else {
#      for (i in 1:input$md_steps){
#        md_nvt(pos,input$md_tstep,input$md_temp,input$PERIODIC)
#      }
#    }
#  }
    
  
  pos <- reactive({
    N          <- input$NW
    model      <- input$water_model
    system_type<- input$var
    periodic   <- TRUE
    if(input$PERIODIC=="No")periodic<-FALSE
    
    if(system_type=='pure water'){
      metal<-"XXX"
    } else {
      metal <- substr(system_type, start =1, stop =nchar(system_type)-2)
    }

    # Assign the output to the variable tmp
    tmp <- generate_initial_configuration(N,system_type,periodic,metal,model)
    
    # L is the box length; it is required by the RGL visualizer below
    L=tmp[[1]]
    
    # pos is a data frame of the atomic positions and atomic symbols
    # It is also required by the RGL visualizer below
    
    pos=tmp[[2]]
    # nr is the total number of atoms (total number of rows)
    # It is required by the RGL visualizer
    nr = nrow(pos)
    # output.xyz is the output written in stanadard materials modeling format
    write.table(nr, file="outfile.xyz", sep="\t", col.names = F, row.names = F,quote = F,append=F)
    write.table(L, file="outfile.xyz", sep="\t", col.names = F, row.names = F,quote = F,append=T)
    write.table(pos, file="outfile.xyz", sep="\t", col.names = F, row.names = F,quote = F,append=T)
    return(list(pos,L,nr,periodic))
  })
  
 
  # First graphical output  
  output$myplot <- renderRglwidget({
    #Set the number of water molecules N
    x=pos()[[1]]$X
    y=pos()[[1]]$Y
    z=pos()[[1]]$Z
    L=pos()[[2]]
    nr=pos()[[3]]
    rcut=1.1
    bond_color="blue"
    bond_thickness=5
    
    #the model will pop up in a small window controled by rgl_init()
    rgl.open(useNULL=T)
    rgl.bg(color = "white" )
    rgl.clear(type = c("shapes", "bboxdeco"))
    rgl.viewpoint(theta = -10, phi = 10, zoom =0.7)
    #rgl_init()
    rgl.spheres(x, y, z, r = 0.3, color = get_colors(pos()[[1]]$ATOM)) 
    if(nr%%3==1)rgl.spheres(x[nr], y[nr], z[nr], r = 0.5, color = "green") 
    rgl_add_axes(x, y, z, L, show.bbox = TRUE)
    rgl_bond(x,y,z,nr,rcut,bond_color,bond_thickness)
    rglwidget()
  })
  
  output$myrdf <- renderPlot({
    tmp <- rdfplot(pos()[[1]],pos()[[2]],pos()[[4]])
    nr<-pos()[[3]]
    # output.xyz is the output written in stanadard materials modeling format
    write.table(nr, file="rdf.dat", sep="\t", col.names = F, row.names = F,quote = F,append=F)
    write.table(tmp, file="rdf.dat", sep="\t", col.names = F, row.names = F,quote = F,append=T)
    ggplot(tmp, aes(x = r0, y =rdfoutput)) + geom_point() + 
      geom_line(size=2)+xlab("r[O-O] in Ang") + ylab("g(r[O-O])") + 
      ggtitle("O-O partial radial distribution function") +  
      theme(
        plot.title = element_text(color="red", size=25, face="bold"),
        axis.title.x = element_text(color="blue", size=25, face="bold"),
        axis.title.y = element_text(color="blue", size=25, face="bold"),
        text = element_text(size=30),
        axis.line = element_line(size = 2, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.ticks.length.y = unit(.5, "cm"),
        axis.ticks.length.x = unit(.5, "cm"),
        axis.text.y = element_text(margin = margin(t =0.3, unit = "cm")),
        axis.text.x = element_text(margin = margin(t =0.3, unit = "cm"))
      )
  })
  
}


shinyApp(ui, server)
