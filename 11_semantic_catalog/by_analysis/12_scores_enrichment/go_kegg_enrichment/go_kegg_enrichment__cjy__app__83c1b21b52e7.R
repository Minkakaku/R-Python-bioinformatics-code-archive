library(stringr)
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(readr)
library(DT)
library(dplyr)
library(ggplot2)
library(shinyjqui)
library(RColorBrewer)
library(plotly)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(shinycssloaders)
library(heatmaply)
library(ggrepel)
library(waiter)
library(pheatmap)
library(shinyFiles)
library(ggsci)
library(RCurl)
library(data.table)
library(visNetwork)
library(rmarkdown)
library(rjson)


mycolors=c("red","yellow","aqua","blue","light-blue","green","navy","teal","olive","lime","orange","fuchsia","purple","maroon","black")
mycoin=c('database','calendar','chart-pie','comments','comment','dna',"comment")
sub.go.list=c("GO:0009507","GO:0005737","GO:0005856","GO:0005783","GO:0005576","GO:0005794","GO:0005739","GO:0016604","GO:0005777","GO:0005886","GO:0005774")
main.data=data.frame(ID=c("Study name: ","Institute: ","Project ID: ","Report date: "),value=c("Quantitative proteomics study of TMT labeled","CWMDA","PJ0001",format(Sys.Date(),'%Y-%m-%d')))

get_uniprot=function(data,api_url = "https://rest.uniprot.org",method = "idmapping/run"){
  uri=paste(api_url,method,sep="/")
  r=postForm(uri,
             ids = paste(data$`Protein accession`,collapse = ","),
             from='UniProtKB_AC-ID',
             to='UniProtKB'
  )
  jobid=fromJSON(json_str = r)$jobId
  while(fromJSON(getURI(paste(api_url,"idmapping/status",jobid,sep="/")))$jobStatus=='RUNNING'){
    Sys.sleep(5)
  }
  if(fromJSON(getURI(paste(api_url,"idmapping/status",jobid,sep="/")))$jobStatus=='FINISHED'){
    result_URL=fromJSON(getURI(paste(api_url,"idmapping/details",jobid,sep="/")))$redirectURL
    qurl=paste0(result_URL,'?format=tsv&size=500&fields=protein_name,gene_names,go_id,xref_kegg,xref_eggnog,xref_pfam,xref_string')
    fread(URLencode(qurl))
  }
  else{
    return(NULL)
  }
}

summary.db=function(){
  list.summary=lapply(dir("DB",pattern = "summary.tsv",recursive = T,full.names = T),function(x) fread(x))
  data=rbindlist(list.summary)
  data$`Create time`=format(data$`Create time`,"%Y-%m-%d %H:%M:%S")
  return(data)
}

enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}
enrich.demo=function(go.data,diff.list){
  go.data.diff=subset(go.data,`Protein accession` %in% diff.list)
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state=subset(go.data.state,mapping>1&`Fold enrich`>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr"))
  return(go.data.state)
}
enrich.combine=function(go.data,diff.list){
  a=split.data.frame(diff.list,diff.list$Regulation)
  enrich.out=lapply(a,function(x) enrich.demo(go.data,x$`Protein accession`))
  enrich.out.data=rbindlist(enrich.out,idcol = "Regulation")
}
classif.go=function(go.data,diff.list){
  go.data=merge(go.data,diff.list,by="Protein accession")
  go.data=subset(go.data,`GO levels`==2,select=c("Category","Terms","Regulation"))
  data=subset(as.data.frame(table(go.data)),Freq>0&nchar(as.character(Terms))<60)
  data=data[order(data$Category,data$Regulation,data$Freq),]
  data$Terms=factor(data$Terms,levels = unique(data$Terms))
  data
}
classif.nog=function(nog.data,diff.list){
  nog.data=merge(nog.data,diff.list,by="Protein accession")
  nog.data=subset(nog.data,select=c("Terms","Regulation"))
  data=subset(as.data.frame(table(nog.data)),Freq>0)
  data=data[order(data$Regulation,data$Freq),]
  data$Terms=factor(data$Terms,levels = unique(data$Terms))
  data
}
classif.subcell=function(go.data,diff.list){
  go.data=merge(go.data,diff.list,by="Protein accession")
  go.data=subset(go.data,ID %in% sub.go.list,select=c("Terms","Regulation"))
  data=subset(as.data.frame(table(go.data)),Freq>0)
  data=data[order(data$Regulation,data$Freq),]
  data$Terms=factor(data$Terms,levels = unique(data$Terms))
  return(data)
}
classif.kegg=function(kegg.data,diff.list){
  kegg.data=merge(kegg.data,diff.list,by="Protein accession")
  kegg.data=subset(kegg.data,select=c("Pathway level 1","Pathway level 2","Regulation"))
  data=subset(data.table(table(kegg.data)),N>0) %>% dplyr::rename(Freq=N)
  data=data[order(data$`Pathway level 1`,data$Regulation,data$Freq),]
  data$`Pathway level 2`=factor(data$`Pathway level 2`,levels = unique(data$`Pathway level 2`))
  return(data)
}
get.ppi=function(annot.data,diff.list){
  data=merge(subset(diff.list,select = c("Protein accession","Regulation")),
             subset(annot.data,!is.na(`Cross-reference (STRING)`),select=c("Entry","Cross-reference (STRING)")),
             by.x="Protein accession",by.y="Entry")
  
  data=separate_rows(data,`Cross-reference (STRING)`)
  data=subset(data,select=c("Cross-reference (STRING)","Regulation")) %>% distinct(`Cross-reference (STRING)`,.keep_all = T)
  colnames(data)=c("Label","Regulation")
  data$ID=c(1:nrow(data))
  string_api_url = "https://string-db.org/api/"
  output_format = "tsv"
  method = "network"
  uri=paste(string_api_url,output_format,method,sep="/")
  species=names(sort(table(sub(data$Label,pattern = "\\..*",replacement = "")),decreasing = T)[1])
  params = list(
    "identifiers"=paste(collapse ="%0d",data$Label),
    "species"=species,
    "required_score"=700
  )
  response=postForm(uri = uri,.params = params)
  out.data=fread(response)
  out.data=merge(out.data,data,by.x="stringId_A",by.y="Label")
  out.data=merge(out.data,data,by.x="stringId_B",by.y="Label")
  edge.data=subset(out.data,select=c("ID.x","ID.y"))
  colnames(edge.data)=c("from","to")
  node.data1=subset(out.data,select=c("ID.x","preferredName_A","Regulation.x"))
  node.data2=subset(out.data,select=c("ID.y","preferredName_B","Regulation.y"))
  colnames(node.data1)=c("id","label","group")
  colnames(node.data2)=c("id","label","group")
  node.data=data.table(rbind(node.data1,node.data2))
  node.data[,value:=length(label),by=id]
  node.data=unique(node.data)
  return(list(edge=edge.data,node=node.data))
}
ppi.top=function(ppi.data,topn){
  node=ppi.data$node
  edge=ppi.data$edge
  node=node[order(-value)]
  topid=head(node,topn)$id
  edge2=subset(edge, from %in% topid | to %in% topid)
  allid=unique(c(edge2$from,edge2$to))
  node=subset(node,id %in% allid)
  edge=subset(edge, from %in% allid & to %in% allid)
  return(list(node=node,edge=edge))
}
ppi.layout=fread(file = "source/ppi_layout.txt")




# -----------------------------------------------------
# -----------------------------------------------------
# -----------------------------------------------------
shinyApp(
  #------------------ ui ------------------------------------
  ui = shinydashboardPlus::dashboardPage(
    skin = "red",
    header = dashboardHeader(title = "PDAP",userOutput("user")),
    sidebar = dashboardSidebar(
      sidebarMenu(
        menuItem("Fileload & QC",icon = icon("hourglass-start"),tabName = "tab1",badgeLabel = "Step 1",badgeColor = "green"),
        menuItem("Classify & Enrich",icon = icon("chart-bar"),tabName = "tab2",badgeLabel = "Step 4",badgeColor = "green"),
        menuItem("Network",icon = icon("project-diagram"),tabName = "tab3",badgeLabel = "Step 6",badgeColor = "green"),
        menuItem("Create report",icon = icon("html5"),tabName = "tab4",badgeLabel = "Step 7",badgeColor = "green")
      )
    ),
    body = dashboardBody(
      autoWaiter(html = spin_whirly(),color = "white"),
      tabItems(
        tabItem("tab1",
                fluidRow(
                  box(width = 12,status = "danger",title = "MaxQuant output file upload",
                      fluidRow(
                        column(width = 6,
                               pickerInput("analysis.db",label = "Select database",options=list(style = "btn-danger",`live-search` = TRUE),choices = summary.db()$`Database name`)
                               ),
                        column(width = 6,
                               fileInput("gene.list",label = "Gene list")
                               )
                        ),
                      dataTableOutput("protein.table")
                  )
                )
        ),
        tabItem("tab2",
                fluidRow(
                  box(status = "danger",headerBorder = F,title = "Functional classification & enrichment",width = 12,
                      sidebar = boxSidebar(id="sidebar2",width = 25,
                                           pickerInput(inputId = 'enrich.color',label = "Enrich plot color",
                                                       choices = c("red","pink","purple","deep-purple","indigo","blue",
                                                                   "light-blue","cyan","teal","green","light-green","lime",
                                                                   "yellow","amber","orange","deep-orange","brown","grey","blue-grey"),
                                                       selected = 'teal'),
                                           sliderInput(inputId = 'enrich.color.a',label = 'color alpha',min = 0.5,max = 1,value = 1),
                                           sliderTextInput(inputId = 'enrich.size.a',label = 'Point size',choices = c(1:12),selected = c(3,8))
                      ),
                      tags$label("Analysis"),tags$br(),
                      actionBttn("do1",label = "Run",icon = icon("satellite-dish"),size = 'sm',color = "danger"),
                      uiOutput('ui1')
                  )
                )
              ),
        tabItem("tab3",
                fluidRow(
                  box(status = "danger",headerBorder = F,title = "Protein-protein interaction network",width = 12,
                      fluidRow(
                        column(width = 3,
                               actionBttn(inputId = 'do2',label = "Get PPI",icon = icon("satellite-dish"),color = 'danger',size = 'sm'),
                               sliderInput(inputId = 'ppi.topn',label = "Maxima network nodes",min = 10,max = 200,value = 50),
                               pickerInput(inputId = "ppi.layout",label = "Network layout",choices = ppi.layout$Layout,multiple = F,selected = "layout_with_kk"),
                               sliderInput(inputId = 'fontsize',label = "Node font size",min = 20,max = 100,value = 30),
                               colorPickr(inputId = 'fontcolor',label = 'Node font color',selected = '#252525',update = "changestop",interaction = list(clear = FALSE,save = FALSE)),
                               pickerInput(inputId = 'colorpal',choices = rownames(brewer.pal.info),selected = 'Dark2',label = "Color palettes"),
                               uiOutput("color.ui")
                        ),
                        column(width = 9,uiOutput('ui3'))
                      )
                  )
                )
              ),
        tabItem("tab4",
                fluidRow(
                  box(status = "danger",headerBorder = F,title = "Create html report doucment",width = 12,collapsible = T,
                      column(width = 9,dataTableOutput('project.info')),
                      column(width=3,
                             selectInput(inputId = 'lang',label = 'Language',choices = c('Chinese','English'),selected = 'Chiniese'),
                             fluidRow(actionBttn(inputId = 'do3',label = "Create report",icon = icon("satellite-dish"),color = 'danger',size = 'sm')),
                             fluidRow(downloadButton(outputId = 'report_zip',label= "Download report",icon = icon("satellite-dish")))
                      )
                  )
                ),
                withSpinner(uiOutput('markdown'))
              )

        
      )
    ),
    title = "Function enrichment"
  ),
  server = function(input,output,session){
    values = reactiveValues(project.info=main.data)
    observeEvent(input$gene.list,{
      data = fread(input$gene.list$datapath)
      if(sum(colnames(data) %in% "Regulation")==0){
        data$Regulation = "All"
      }
      values$protein.data = data
    })
    
    output$protein.table = renderDataTable(values$protein.data,rownames=F,options=list(
      pageLength = 5,scrollX=TRUE,
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
        "}")
    ))
    observeEvent(input$do1,{
      req(values$protein.data)
      go.data=fread(file=paste("DB",input$analysis.db,"GO_annotation.tsv",sep="/"))
      nog.data=fread(file=paste("DB",input$analysis.db,"NOG_annotation.tsv",sep="/"))
      pfam.data=fread(file=paste("DB",input$analysis.db,"Pfam_annotation.tsv",sep="/"))
      kegg.data=fread(file=paste("DB",input$analysis.db,"KEGG_annotation.tsv",sep="/"))
      annot.data=fread(file=paste("DB",input$analysis.db,"annotation.tsv",sep="/"))
      pfam.data$Category="Pfam"
      kegg.data$Category="KEGG"
      gmt.data=rbindlist(
        list(
          subset(go.data,select=c("Protein accession","Category","Terms")),
          subset(pfam.data,select=c("Protein accession","Category","Terms")),
          subset(kegg.data,select=c("Protein accession","Category","Terms"))
        )) %>% distinct_all()
      diff.list = values$protein.data
      diff.list = diff.list %>% group_by(Regulation) %>% distinct(`Protein accession`,.keep_all = T)
      values[["enrich.data"]]=lapply(split(gmt.data,gmt.data$Category),function(x) enrich.combine(x,diff.list))
      values[["goclass.data"]]=classif.go(go.data,diff.list)
      values[["nogclass.data"]]=classif.nog(nog.data,diff.list)
      values[["subcellclass.data"]]=classif.subcell(go.data,diff.list)
      values[["keggclass.data"]]=classif.kegg(kegg.data,diff.list)
      
      values[["Biological process"]]=merge(subset(go.data,Category=='Biological process'),diff.list,by="Protein accession")
      values[["Molecular function"]]=merge(subset(go.data,Category=='Molecular function'),diff.list,by="Protein accession")
      values[["Cellular component"]]=merge(subset(go.data,Category=='Cellular component'),diff.list,by="Protein accession")
      values[["KEGG"]]=merge(kegg.data,diff.list,by="Protein accession")
      values[["Pfam"]]=merge(pfam.data,diff.list,by="Protein accession")
      values[["KOG"]]=merge(nog.data,diff.list,by="Protein accession")
      
      values[["goclass.plot"]]=ggplot(values[["goclass.data"]],aes(x=Terms,y=Freq,fill=Category))+
        geom_bar(stat="identity")+
        scale_fill_manual(values = brewer.pal(3,"Set1"))+theme_test(base_size = 14)+
        facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
        coord_flip()+theme(legend.position = "top",legend.justification = c(1,1))+
        geom_text(aes(label=Freq))+scale_y_continuous(expand = expand_scale(mult = c(0,0.1)))
      values[["nogclass.plot"]]=ggplot(values[["nogclass.data"]],aes(x=Terms,y=Freq,fill=Terms))+
        geom_bar(stat="identity")+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(23))+
        theme_test(base_size = 14)+
        facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
        coord_flip()+theme(legend.position = "none")+
        geom_text(aes(label=Freq))
      values[["subcellclass.plot"]]=ggplot(values[["subcellclass.data"]],aes(x=Terms,y=Freq,fill=Terms))+
        geom_bar(stat="identity")+
        scale_fill_manual(values = brewer.pal(12,"Paired"))+theme_test(base_size = 14)+
        facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
        coord_flip()+theme(legend.position = "none")+
        geom_text(aes(label=Freq))
      values[["keggclass.plot"]]=ggplot(values[["keggclass.data"]],aes(x=`Pathway level 2`,y=Freq,fill=`Pathway level 1`))+
        geom_bar(stat="identity")+
        scale_fill_manual(values = brewer.pal(7,"Set1"))+theme_test(base_size = 14)+
        facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
        coord_flip()+theme(legend.position = "top",legend.justification = c(1,1))+
        guides(fill=guide_legend(nrow = 3))+
        geom_text(aes(label=Freq))+scale_y_continuous(expand = expand_scale(mult = c(0,0.1)))
    }
    )
    observe({
      lapply(names(values[["enrich.data"]]),function(y){
        data.list=values[["enrich.data"]][[y]]
        if(nrow(data.list)<3){values[[paste(y,'enrich.plot',sep=".")]]=NULL}
        else{
          values[[paste(y,'enrich.plot',sep=".")]]=ggplot(subset(data.list,p.value<0.05&nchar(as.character(Terms))<60) %>% group_by(Regulation) %>% slice_min(p.value,n =10), aes(y=reorder(Terms,`Fold enrich`),x=log2(`Fold enrich`),color=p.value))+
            geom_point(aes(size=mapping))+theme_test(base_size = 14)+ylab("Terms")+facet_grid(.~Regulation)+scale_color_material(reverse = T,palette = input$enrich.color,alpha = input$enrich.color.a)+
            scale_size(range = input$enrich.size.a)
        }
      })            
    })
    
    output$ui1 = renderUI({
      output[["goclass.plot"]] = renderPlot({
        req(values[["goclass.plot"]])
        values[["goclass.plot"]]
      })
      output[["nogclass.plot"]] = renderPlot({
        req(values[["nogclass.plot"]])
        values[["nogclass.plot"]]
      })
      output[["subcellclass.plot"]] = renderPlot({
        req(values[["subcellclass.plot"]])
        values[["subcellclass.plot"]]
      })
      output[["keggclass.plot"]] = renderPlot({
        req(values[["keggclass.plot"]])
        values[["keggclass.plot"]]
      })
      
      output[["goclass.data"]]=renderDataTable({
        req(values[["goclass.data"]])
        values[["goclass.data"]]
      },rownames=F,extensions = 'Buttons',server = FALSE,
      options=list(
        dom = 'Blfrtip',scrollX=TRUE,
        buttons = 'excel',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
          "}"))
      )
      
      output[["nogclass.data"]]=renderDataTable({
        req(values[["nogclass.data"]])
        values[["nogclass.data"]]
      },rownames=F,extensions = 'Buttons',server = FALSE,
      options=list(
        dom = 'Blfrtip',scrollX=TRUE,
        buttons = 'excel',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
          "}"))
      )
      output[["subcellclass.data"]]=renderDataTable({
        req(values[["subcellclass.data"]])
        values[["subcellclass.data"]]
      },rownames=F,extensions = 'Buttons',server = FALSE,
      options=list(
        dom = 'Blfrtip',scrollX=TRUE,
        buttons = 'excel',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
          "}"))
      )
      output[["keggclass.data"]]=renderDataTable({
        req(values[["keggclass.data"]])
        values[["keggclass.data"]]
      },rownames=F,extensions = 'Buttons',server = FALSE,
      options=list(
        dom = 'Blfrtip',scrollX=TRUE,
        buttons = 'excel',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
          "}"))
      )
      
      output[["enrich.plot"]] = renderUI({
        req(values$enrich.data)
        data.list = values$enrich.data
        do.call(tabsetPanel,lapply(names(data.list),function(y){
          output[[paste(y,'enrich.plot',sep=".")]]=renderPlot(values[[paste(y,'enrich.plot',sep=".")]])
          output[[paste(y,'enrich.data',sep=".")]]=renderDataTable(data.list[[y]],rownames=F,selection="single",
                                                                     options=list(
                                                                       dom = 'Blfrtip',scrollX=TRUE,
                                                                       buttons = 'excel',
                                                                       initComplete = JS(
                                                                         "function(settings, json) {",
                                                                         "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
                                                                         "}")),extensions = 'Buttons',server = FALSE
          )
          
          output[[paste(y,'enrich.data.detail',sep=".")]]=renderDataTable({
            t1=data.list[[y]]$Terms[input[[paste(y,'enrich.data_rows_selected',sep=".")]]]
            t2=data.list[[y]]$Regulation[input[[paste(y,'enrich.data_rows_selected',sep=".")]]]
            subset(values[[paste(y,sep=".")]],Terms==t1&Regulation==t2)
          },rownames=F,extensions = 'Buttons',server = FALSE,
          options=list(
            dom = 'Blfrtip',scrollX=TRUE,
            buttons = 'excel',
            initComplete = JS(
              "function(settings, json) {",
              "$(this.api().table().header()).css({'background-color': '#555', 'color': '#fff'});",
              "}"))
          )
          
          tabPanel(title = y,width = 12,status = "danger",
                   column(width = 6,
                          withSpinner(jqui_resizable(plotOutput(paste(y,'enrich.plot',sep="."),height = 600,width = 600)))
                   ),
                   column(width = 6,
                          dataTableOutput(paste(y,'enrich.data',sep="."))
                   ),
                   dataTableOutput(paste(y,'enrich.data.detail',sep="."))
          )
          
        }))
      })
      
      
      fluidRow(
        box(width = 12,title = "Functional classfication",solidHeader = T,status = "danger",collapsible = T,
            tabBox(width = 12,
                   tabPanel(title = 'GO',
                            column(width = 6,
                                   withSpinner(jqui_resizable(plotOutput("goclass.plot",height = 800,width = 600)))),
                            column(width = 6,
                                   dataTableOutput("goclass.data"))
                            ),
                   tabPanel(title = 'KEGG',
                            column(width = 6,
                                   withSpinner(jqui_resizable(plotOutput("keggclass.plot",height = 800,width = 800)))),
                            column(width = 6,
                                   dataTableOutput("keggclass.data"))
                            ),
                   tabPanel(title = 'KOG',
                            column(width = 6,
                                   withSpinner(jqui_resizable(plotOutput("nogclass.plot",height = 600,width = 600)))),
                            column(width = 6,
                                   dataTableOutput("nogclass.data"))
                            ),
                   tabPanel(title = 'Subcell',
                            column(width = 6,
                                   withSpinner(jqui_resizable(plotOutput("subcellclass.plot",height = 400,width = 600)))),
                            column(width = 6,
                                   dataTableOutput("subcellclass.data"))
                            )
                   )
            ),
        box(width=12,title = "Functional enrichment",solidHeader = T,status = "danger",collapsible = T,
            uiOutput("enrich.plot")
            )
        )
     })
    
    observeEvent(input$do2,{
      req(values$protein.data)
      annot.data=fread(file=paste("DB",input$analysis.db,"annotation.tsv",sep="/"))
      diff.list=values$protein.data %>% select(`Protein accession`,Regulation)
      diff.list = diff.list %>% group_by(Regulation) %>% distinct(`Protein accession`,.keep_all = T)
      values[["ppi.data"]]=get.ppi(annot.data,diff.list)
    })
    
    output$ui3=renderUI({
      req(values$protein.data)
      if(all(unique(values$protein.data$Regulation) %in% "All") == T){
        output[["ppi.plot"]] <- renderVisNetwork({
          req(values[["ppi.data"]])
          ppi.data=ppi.top(values[["ppi.data"]],input$ppi.topn)
          values[["ppi.plot"]]=visNetwork(ppi.data$node,ppi.data$edge) %>% visInteraction(hover = TRUE) %>%
            visIgraphLayout(layout = input$ppi.layout) %>% visExport(type = "png",float = "left",background = "white") %>%
            visNodes(font=list(size=input$fontsize,color=input$fontcolor)) %>%
            visGroups(groupname = 'All',color=input$nodecolor.up) %>%
            visOptions(highlightNearest = TRUE)
        })
      }else{
        output[["ppi.plot"]] <- renderVisNetwork({
          req(values[["ppi.data"]])
          ppi.data=ppi.top(values[["ppi.data"]],input$ppi.topn)
          values[["ppi.plot"]]=visNetwork(ppi.data$node,ppi.data$edge) %>% visInteraction(hover = TRUE) %>%
            visIgraphLayout(layout = input$ppi.layout) %>% visExport(type = "png",float = "left",background = "white") %>%
            visNodes(font=list(size=input$fontsize,color=input$fontcolor)) %>%
            visGroups(groupname = 'up-regulated',color=input$nodecolor.up) %>%
            visGroups(groupname = 'down-regulated',color=input$nodecolor.down) %>%
            visOptions(highlightNearest = TRUE)
        })
      }
      

        tabPanel(title = "PPI network",
                 jqui_resizable(visNetworkOutput("ppi.plot",height = 800))
        )
    })
    
    
    
    output$color.ui=renderUI({
      div(
        colorPickr(inputId = 'nodecolor.up',label = 'up-regulated',swatches = brewer.pal(brewer.pal.info[input$colorpal,1],input$colorpal),update = "change",interaction = list(clear = FALSE,save = FALSE),selected = "#D95F02"),
        colorPickr(inputId = 'nodecolor.down',label = 'down-regulated',swatches = brewer.pal(brewer.pal.info[input$colorpal,1],input$colorpal),update = "change",interaction = list(clear = FALSE,save = FALSE),selected = "#1B9E77")
      )
    })
    
    
    observeEvent(input$do3,{
      setpn = 5
      withProgress(message = 'Calculation in progress',detail = 'This may take a while...', value = 0, {
        incProgress(amount = 1/setpn,message = "Step1: Create a temporary report folder")
        values$dirname=tempfile()
        dirname=values$dirname
        if(!dir.exists(dirname)){dir.create(path = dirname)}
        
        incProgress(amount = 1/setpn,message = "Step2: DEPs annotation")
        temp.dir=paste(dirname,"1-Annotation_of_DEPs",sep="/")
        if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
        fwrite(values[["Biological process"]],file=paste0(temp.dir,"/Table 1. BP annotation.xls"),sep = "\t")
        fwrite(values[["Cellular component"]],file=paste0(temp.dir,"/Table 1. CC annotation.xls"),sep = "\t")
        fwrite(values[["Molecular function"]],file=paste0(temp.dir,"/Table 1. MF annotation.xls"),sep = "\t")
        fwrite(values[["KEGG"]],file=paste0(temp.dir,"/Table 1. KEGG annotation.xls"),sep = "\t")
        fwrite(values[["Pfam"]],file=paste0(temp.dir,"/Table 1. Pfam annotation.xls"),sep = "\t")
        fwrite(values[["KOG"]],file=paste0(temp.dir,"/Table 1. KOG annotation.xls"),sep = "\t")

        incProgress(amount = 1/setpn,message = "Step3: Functional classification plot")
        temp.dir=paste(dirname,"2-Functional_classification",sep="/")
        if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
        fwrite(values[["goclass.data"]],file=paste0(temp.dir,"/Table 1a. GO classify of.xls"),sep = "\t")
        fwrite(values[["nogclass.data"]],file=paste0(temp.dir,"/Table 1b. KOG classify of.xls"),sep = "\t")
        fwrite(values[["subcellclass.data"]],file=paste0(temp.dir,"/Table 1c. Subcell classify of.xls"),sep = "\t")
        fwrite(values[["keggclass.data"]],file=paste0(temp.dir,"/Table 1d. KEGG classify of.xls"),sep = "\t")
        ggsave(values[["goclass.plot"]],filename=paste0(temp.dir,"/Figure 1a. GO classify of.pdf"),width = 8,height = 8)
        ggsave(values[["nogclass.plot"]],filename=paste0(temp.dir,"/Figure 1b. KOG classify of.pdf"),width = 8,height = 6)
        ggsave(values[["subcellclass.plot"]],filename=paste0(temp.dir,"/Figure 1c. Subcell classify of.pdf"),width = 6,height = 6)
        ggsave(values[["keggclass.plot"]],filename=paste0(temp.dir,"/Figure 1d. KEGG classify of.pdf"),width = 8,height = 8)
        ggsave(values[["goclass.plot"]],filename=paste0(temp.dir,"/Figure 1a. GO classify of.png"),width = 8,height = 8)
        ggsave(values[["nogclass.plot"]],filename=paste0(temp.dir,"/Figure 1b. KOG classify of.png"),width = 8,height = 6)
        ggsave(values[["subcellclass.plot"]],filename=paste0(temp.dir,"/Figure 1c. Subcell classify of.png"),width = 6,height = 6)
        ggsave(values[["keggclass.plot"]],filename=paste0(temp.dir,"/Figure 1d. KEGG classify of.png"),width = 8,height = 8)
        
        incProgress(amount = 1/setpn,message = "Step4: Functional enrichment plot")
        temp.dir=paste(dirname,"3-Functional_enrichment",sep="/")
        if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
        enrich.data=rbindlist(values[["enrich.data"]],idcol = "Fun Category")
        fwrite(enrich.data,file=paste0(temp.dir,"/Table 2. Functional enrichment.xls"),sep = "\t")
        for(j in names(values[["enrich.data"]])){
          ggsave(values[[paste(j,"enrich.plot",sep=".")]],filename=paste0(temp.dir,"/Figure 2. ",j," enrichment.pdf"),width = 8,height = 8)
          ggsave(values[[paste(j,"enrich.plot",sep=".")]],filename=paste0(temp.dir,"/Figure 2. ",j," enrichment.png"),width = 8,height = 8)
        }
        
        incProgress(amount = 1/setpn,message = "Step5: Get PPI")
        temp.dir=paste(dirname,"4-Protein-protein_interaction",sep="/")
        if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
        if(!is.null(values[["ppi.plot"]])){
          visSave(values[["ppi.plot"]], file = paste0(temp.dir,"/Figure 4. PPI.html"),selfcontained = T)
        }
        

        file.copy(from = 'source/report.Rmd',to = paste(dirname,"report.Rmd",sep="/"),overwrite = T)
        file.copy(from = 'source/logo.png',to = paste(dirname,"logo.png",sep="/"),overwrite = T)
        report=paste(dirname,"report.html",sep="/")
        rmarkdown::render(input = paste(dirname,"report.Rmd",sep="/"),output_format = 'html_document',
                          output_file = report)
        unlink(paste(dirname,"report.Rmd",sep="/"))
        unlink(paste(dirname,"logo.png",sep="/"))
        temp_report=paste(gsub(dirname,pattern = ".*/",replacement =""),"report.html",sep="_")
        file.copy(from = report,to = paste("www",temp_report,sep="/"),overwrite = T)
        values$report_ref=temp_report
      })
    })
    
    output$markdown <- renderUI({
      req(values$report_ref)
      tags$iframe(src=values$report_ref,width="100%",height=800)
    })
    
    output$report_zip=downloadHandler(
      filename = "Report.tar.zip",
      content = function(file){
        command=paste("tar czf",file,"-C",tempdir(),gsub(values$dirname,pattern = ".*/",replacement = ""),sep=" ")
        system(command = command)
      }
    )
    
    output$project.info=renderDataTable({
      datatable(values$project.info,selection = 'none',editable=list(target="cell",disable = list(columns = 0)),rownames=F,
                options = list(dom = 't'),style = "auto") %>%
        formatStyle(1,`font-size` = '20px',`font-weight`="bold",color="#377EB8",`text-align`="right") %>%
        formatStyle(2,`font-size` = '20px',`font-weight`="bold",`text-align`="left")
    })
    observeEvent(input$project.info_cell_edit, {
      values$project.info = editData(values$project.info, input$project.info_cell_edit, 'table2',rownames = F)
    })
    
  }
)
                        
