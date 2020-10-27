# 2D Density Generator

library(shiny)
library(mvtnorm)
library(ggplot2)
#The useful density generating functions
bivarnorm.list <- function(arg.list){
    require(mvtnorm)
    
    x.val <- arg.list[1]
    y.val <- arg.list[2]
    mean.x <- arg.list[3]
    mean.y <- arg.list[4]
    sigma.x <- arg.list[5]
    sigma.y <- arg.list[6]
    
    xy <- c(x.val,y.val)
    
    densval <- dmvnorm(x=xy,mean=c(mean.x,mean.y),sigma=diag(c(sigma.x,sigma.y)))
    
    return(densval)
}
calculate.density.2d <- function(x.val,y.val,peak.x,peak.y,spread.x,spread.y){
    xcol <- rep(x.val,length(peak.x))
    ycol <- rep(y.val,length(peak.y))
    
    arg.matrix <- cbind(xcol,ycol,peak.x,peak.y,spread.x,spread.y)
    
    dens.val <- sum(apply(arg.matrix,1,bivarnorm.list))
    
    return(dens.val)
}
calculate.density.2d.list <- function(xylist,peak.x,peak.y,spread.x,spread.y){
    x <- xylist[1]
    y <- xylist[2]
    
    return(calculate.density.2d(x=x,y=y,peak.x=peak.x,peak.y=peak.y,spread.x=spread.x,spread.y=spread.y))
}
grid.density.2d <- function(peak.x,peak.y,spread.x,spread.y,x.lim,y.lim,resolution=1){
    x.grid <- seq(from=x.lim[1],to=x.lim[2],by=resolution)
    y.grid <- seq(from=y.lim[1],to=y.lim[2],by=resolution)
    
    full.grid <- expand.grid(x.grid,y.grid)
    
    density.vals <- apply(full.grid,1,calculate.density.2d.list,peak.x=peak.x,peak.y=peak.y,spread.x=spread.x,spread.y=spread.y)
    
    full.dens <- cbind(full.grid,density.vals)
    
    colnames(full.dens) <- c("x","y","density")
    
    return(full.dens)
}
plot.density.2d <- function(peak.x,peak.y,spread.x,spread.y,x.lim,y.lim,resolution=1){
    require(ggplot2)
    
    dens.grid <- grid.density.2d(peak.x=peak.x,peak.y=peak.y,spread.x=spread.x,spread.y=spread.y,x.lim=x.lim,y.lim=y.lim,resolution=resolution)
    
    dens.plot <- ggplot(data=dens.grid, mapping=aes(x=x, y=y, z = density)) + 
        geom_contour_filled() +
        xlab("X") + ylab("Y")
    
    return(dens.plot)
}
generate.density.func.text <- function(max.x.coord,max.y.coord,spreads.x,spreads.y){
    func.text <- paste0("density.func <- function(xval,yval){return(calculate.density.2d(x = xval,y = yval,peak.x = c(",
                        paste0(max.x.coord,collapse=","),
                        ") ,peak.y = c(",
                        paste0(max.y.coord,collapse=","),
                        "),spread.x = c(",
                        paste0(spreads.x,collapse=","),
                        "),spread.y = c(",
                        paste0(spreads.y,collapse=","),
                        ")))}"
                        ,collapse="")
    
    return(func.text)
}
generate.rand.density.func <- function(x.lim=c(0,100),y.lim=c(0,100),resolution=1,n.peaks=sample(2:5,1)){
    max.x.coord <- sample(seq(from=min(x.lim),to=max(x.lim),by=resolution),n.peaks,replace=T)
    max.y.coord <- sample(seq(from=min(y.lim),to=max(y.lim),by=resolution),n.peaks,replace=T)
    
    x.range <- max(x.lim) - min(x.lim)
    y.range <- max(y.lim) - min(y.lim)
    
    spreads.x <- sample(seq(from=(0.5*x.range),to=(5*x.range),by=resolution),n.peaks,replace=T)
    spreads.y <- sample(seq(from=(0.5*y.range),to=(5*y.range),by=resolution),n.peaks,replace=T)
    
    rand.func <- function(xval,yval){
        calculate.density.2d(x = xval,y = yval,peak.x = max.x.coord,peak.y = max.y.coord,spread.x = spreads.x,spread.y = spreads.y)
    }
    
    return(rand.func)
}
generate.rand.density.func.text <- function(x.lim=c(0,100),y.lim=c(0,100),resolution=1,n.peaks=sample(2:5,1)){
    max.x.coord <- sample(seq(from=min(x.lim),to=max(x.lim),by=resolution),n.peaks,replace=T)
    max.y.coord <- sample(seq(from=min(y.lim),to=max(y.lim),by=resolution),n.peaks,replace=T)
    
    x.range <- max(x.lim) - min(x.lim)
    y.range <- max(y.lim) - min(y.lim)
    
    spreads.x <- sample(seq(from=(0.5*x.range),to=(5*x.range),by=resolution),n.peaks,replace=T)
    spreads.y <- sample(seq(from=(0.5*y.range),to=(5*y.range),by=resolution),n.peaks,replace=T)
    
    func.text <- paste0("rand.func <- function(xval,yval){return(calculate.density.2d(x = xval,y = yval,peak.x = c(",
                        paste0(max.x.coord,collapse=","),
                        ") ,peak.y = c(",
                        paste0(max.y.coord,collapse=","),
                        "),spread.x = c(",
                        paste0(spreads.x,collapse=","),
                        "),spread.y = c(",
                        paste0(spreads.y,collapse=","),
                        ")))}"
                        ,collapse="")
    
    return(func.text)
}
generate.rand.density.params <- function(x.lim=c(0,100),y.lim=c(0,100),resolution=1,n.peaks=sample(2:5,1)){
    max.x.coord <- sample(seq(from=min(x.lim),to=max(x.lim),by=resolution),n.peaks)
    max.y.coord <- sample(seq(from=min(y.lim),to=max(y.lim),by=resolution),n.peaks)
    
    x.range <- max(x.lim) - min(x.lim)
    y.range <- max(y.lim) - min(y.lim)
    
    spreads.x <- sample(seq(from=(0.5*x.range),to=(5*x.range),by=resolution),n.peaks)
    spreads.y <- sample(seq(from=(0.5*y.range),to=(5*y.range),by=resolution),n.peaks)
    
    return(list(max.x.coord,max.y.coord,spreads.x,spreads.y,x.lim,y.lim))
}
generate.rand.density.func.grid <- function(x.lim=c(0,100),y.lim=c(0,100),resolution=1,n.peaks=sample(2:5,1)){
    param.list <- generate.rand.density.params(x.lim=x.lim,y.lim=y.lim,resolution=resolution,n.peaks=n.peaks)
    
    density.func <- function(xval,yval){
        calculate.density.2d(x=xval,y=yval,peak.x=param.list[[1]],peak.y=param.list[[2]],spread.x=param.list[[3]],spread.y=param.list[[4]])
    }
    
    density.func.text <- paste0("density.function <- function(xval,yval){return(calculate.density.2d(x = xval,y = yval,peak.x = c(",
                                paste0(param.list[[1]],collapse=","),
                                ") ,peak.y = c(",
                                paste0(param.list[[2]],collapse=","),
                                "),spread.x = c(",
                                paste0(param.list[[3]],collapse=","),
                                "),spread.y = c(",
                                paste0(param.list[[4]],collapse=","),
                                ")))}"
                                ,collapse="")
    
    density.grid <- grid.density.2d(peak.x=param.list[[1]],peak.y=param.list[[2]],spread.x=param.list[[3]],spread.y=param.list[[4]],x.lim=x.lim,y.lim=y.lim)
    
    density.plot <- plot.density.2d(peak.x=param.list[[1]],peak.y=param.list[[2]],spread.x=param.list[[3]],spread.y=param.list[[4]],resolution=resolution,x.lim=x.lim,y.lim=y.lim)
    
    return(list(density.func.text,density.grid,density.plot))  
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("2D spatial density generator"),

    #
    sidebarLayout(
        sidebarPanel(
            actionButton(inputId = "randParam","Generate Random Parameters")
            ,
            fluidRow(
                column(3,textInput("xmin","Minimum X")),
                column(3,textInput("xmax","Maximum X")),
                column(3,textInput("ymin","Minimum Y")),
                column(3,textInput("ymax","Maximum Y"))
            )
            ,
            fluidRow(
                column(6,textInput("peakx","X Coordinates of Peaks (enter numbers separated by comma)")),
                column(6,textInput("peaky","Y Coordinates of Peaks (enter numbers separated by comma)"))
            )
            ,
            fluidRow(
                column(6,textInput("spreadx","Horizontal Sigmas (enter numbers separated by comma)")),
                column(6,textInput("spready","Vertical Sigmas (enter numbers separated by comma)"))
            )
            ,
            # fluidRow(
            #     column(6,textInput("resolution","Resolution of Plot"))
            # )
            # ,
            fluidRow(
                column(6,actionButton(inputId = "generateDens","Generate Density Function")),
                column(6,actionButton(inputId = "generatePlot","Generate Density Plot"))
            )
            ,
            fluidRow(
                verbatimTextOutput("funcCode",placeholder=T)
            )
        ),

        # Show a plot of the generated density
        mainPanel(
            plotOutput(outputId = "densMap",width="auto",height="auto")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$randParam,{
        newParams <- generate.rand.density.params()
        
        updateTextInput(session,"xmin",value=newParams[[5]][1])
        updateTextInput(session,"xmax",value=newParams[[5]][2])
        updateTextInput(session,"ymin",value=newParams[[6]][1])
        updateTextInput(session,"ymax",value=newParams[[6]][2])
        updateTextInput(session,"peakx",value=newParams[[1]])
        updateTextInput(session,"peaky",value=newParams[[2]])
        updateTextInput(session,"spreadx",value=newParams[[3]])
        updateTextInput(session,"spready",value=newParams[[4]])
        # updateTextInput(session,"resolution",value=1)
    })
    
    observeEvent(input$generateDens,{
        funtext <- generate.density.func.text(input$peakx,input$peaky,input$spreadx,input$spready)
        funstring <- "bivarnorm.list <- function(arg.list){
                require(mvtnorm)
                
                x.val <- arg.list[1]
                y.val <- arg.list[2]
                mean.x <- arg.list[3]
                mean.y <- arg.list[4]
                sigma.x <- arg.list[5]
                sigma.y <- arg.list[6]
                
                xy <- c(x.val,y.val)
                
                densval <- dmvnorm(x=xy,mean=c(mean.x,mean.y),sigma=diag(c(sigma.x,sigma.y)))
                
                return(densval)
            }
            \n
calculate.density.2d <- function(x.val,y.val,peak.x,peak.y,spread.x,spread.y){
                xcol <- rep(x.val,length(peak.x))
                ycol <- rep(y.val,length(peak.y))
                
                arg.matrix <- cbind(xcol,ycol,peak.x,peak.y,spread.x,spread.y)
                
                dens.val <- sum(apply(arg.matrix,1,bivarnorm.list))
                
                return(dens.val)
            }
            \n"
        fulltext <- paste0(funstring,funtext,collapse="")
        output$funcCode <- renderText({fulltext})
    })
    
    observeEvent(input$generatePlot,{
        inputxlim <- c(as.numeric(input$xmin),as.numeric(input$xmax))
        inputylim <- c(as.numeric(input$ymin),as.numeric(input$ymax))
        inputpeakx <- as.numeric(strsplit(input$peakx,split = ",")[[1]])
        inputpeaky <- as.numeric(strsplit(input$peaky,split = ",")[[1]])
        inputspreadx <- as.numeric(strsplit(input$spreadx,split = ",")[[1]])
        inputspready <- as.numeric(strsplit(input$spready,split = ",")[[1]])
        
        plotresolution <- round(mean((inputxlim[[2]]-inputxlim[[1]])/30,(inputylim[[2]]-inputylim[[1]])/30))
        
        output$testpeakx <- renderText({inputpeakx})
        output$testpeaky <- renderText({input$peaky})
        
        output$testspreadx <- renderText({inputspreadx})
        output$testspready <- renderText({inputspready})
        
        widthpix <- round(800*(diff(inputxlim)/diff(inputylim)))
        heightpix <- round(600*(diff(inputylim)/diff(inputxlim)))
        
        output$densMap <- renderPlot(plot.density.2d(peak.x=inputpeakx,peak.y=inputpeaky,
                                    spread.x=inputspreadx,spread.y=inputspready,
                                    x.lim=inputxlim,y.lim=inputylim,
                                    resolution=plotresolution),width=widthpix,height=heightpix)
        
        # output$densMap <- renderPlot(plot.density.2d(peak.x=c(60,40),peak.y=c(50,50),
        #                         spread.x=c(40,50),spread.y=c(100,200),
        #                         x.lim=inputxlim,y.lim=inputylim, resolution=1))
        
        # densGrid <- grid.density.2d(peak.x=inputpeakx,peak.y=inputpeaky,
        #                             spread.x=inputspreadx,spread.y=inputspready,
        #                             x.lim=inputxlim,y.lim=inputylim,
        #                             resolution=as.numeric(input$resolution))
        # 
        # output$densMap <- renderPlot(ggplot(data=densGrid, mapping=aes(x=x, y=y, z = density)) + 
        #                                  geom_contour_filled() +
        #                                  xlab("X") + ylab("Y"))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)