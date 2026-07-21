library(ggplot2)
library(patchwork)
library(RColorBrewer)

gseaplot <- function(gseadata,pathwayID,customtitle = FALSE,showNES = TRUE,showPvalue = TRUE,
                     showFDR = TRUE,showCase = TRUE,showCon = TRUE,ES_linetype = "solid",ES_linecolor = "#41AB5D",
                     ES_linesize = 2,fz_linetype = "solid",fz_linecolor = "grey",fz_linesize = 2,p1ytitlefontface = "bold",
                     p1ytextsize = 12,p1ytextfontface = "plain",p1title = "AAAAAAAAAAA",p1titlesize = 20,p1ytitlesize = 16,
                     p1titlefontface = "bold",p1titlecolor = "black",p1titlehjust = 0.5,p1titlevjust = 0.5,NESx = 0.85,NESy = 0.9,
                     p1bordersize = 2,NESsize = 5,NEScolor = "black",NESfontface = "italic",Pvaluex = 0.85,Pvaluey = 0.8,Pvaluesize = 5,
                     Pvaluecolor = "black",Pvaluefontface = "italic",FDRx = 0.85,FDRy = 0.7,FDRsize = 5,FDRcolor = "black",
                     FDRfontface = "italic",NESdigit = 2,Pvaluedigit = 2,FDRdigit = 2,p2_linetype = 1,p2_linecolor = "black",
                     p2_linesize = 0.8,rec_height_ratio = 0.3,rec_height_alpha = 0.8,rec_col_high = "Reds",rec_col_low = "Blues",
                     p3ytitlesize = 16,p3ytitlefontface = "plain",p3xtitlesize = 16,p3xtitlefontface = "bold",
                     p3xtextsize = 12,p3xtextfontface = "plain",p3bordersize = 1,casename = "Case",
                     conname = "Con",Casex = 0.1,Casey = 0.6,Casesize = 5,Casecolor = "red",Casefontface = "italic",
                     Conx = 0.9,Cony = -0.6,Consize = 5,Concolor = "blue",Confontface = "italic",
                     p1ratio = 1,p2ratio = 0.2,p3ratio=0.5){


    maintitle = gseadata@result[pathwayID, "Description"]
    nes = gseadata@result[pathwayID, "NES"]
    enrichmentScore = gseadata@result[pathwayID, "enrichmentScore"]
    pvalue = gseadata@result[pathwayID, "pvalue"]
    fdr = gseadata@result[pathwayID, "p.adjust"]

    if (!customtitle) {
        Maintitle = maintitle
    } else {
        Maintitle = p1title
    }
    if (showNES) {
        NES = paste0("NES = ",round(nes, NESdigit))
        if (enrichmentScore >= 0) {
        NESx = length(gseadata@geneList) * NESx
        NESy = enrichmentScore * NESy
        } else {
            NESx = length(gseadata@geneList) * (1-NESx)
            NESy = enrichmentScore * (NESy-0.2)
            }
        NESsize = NESsize
        NEScolor = NEScolor
        NESfontface = NESfontface
    } else {
        NES = enrichmentScore
        NESsize = 0
        NESx = length(gseadata@geneList) * 2
        NESy = enrichmentScore * (NESy-0.2)
        NEScolor = "white"
        NESfontface = NESfontface
    }

    if (showPvalue) {
        Pvalue = paste0("P value = ",format(pvalue, scientific = TRUE, digits = Pvaluedigit))
        if (enrichmentScore >= 0) {
            Pvaluex = length(gseadata@geneList) * Pvaluex
            Pvaluey = enrichmentScore * Pvaluey
        } else {
            Pvaluex = length(gseadata@geneList) * (1-Pvaluex)
            Pvaluey = enrichmentScore * Pvaluey
        }
        Pvaluesize = Pvaluesize
        Pvaluecolor = Pvaluecolor
        Pvaluefontface = Pvaluefontface
    } else {
        Pvalue = enrichmentScore
        Pvaluesize = 0
        Pvaluex = length(gseadata@geneList) * 2
        Pvaluey = enrichmentScore * Pvaluey
        Pvaluecolor = "white"
        Pvaluefontface = Pvaluefontface
    }

    if (showFDR) {
        FDR = paste0("FDR = ",format(fdr, scientific = TRUE, digits = FDRdigit))
        if (enrichmentScore >= 0) {
            FDRx = length(gseadata@geneList) * FDRx
            FDRy = enrichmentScore * FDRy
        } else {
            FDRx = length(gseadata@geneList) * (1-FDRx)
            FDRy = enrichmentScore * (FDRy+0.2)
        }
        FDRsize = FDRsize
        FDRcolor = FDRcolor
        FDRfontface = FDRfontface
    } else {
        FDR = enrichmentScore
        FDRsize = 0
        FDRx = length(gseadata@geneList) * 2
        FDRy = enrichmentScore * (FDRy+0.2)
        FDRcolor = "white"
        FDRfontface = FDRfontface
    }
    
    gsdata <- get_gsdata(gseadata, pathwayID)

    p1 <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
        geom_hline(yintercept = 0, color = fz_linecolor, size = fz_linesize, linetype = fz_linetype) +
        geom_line(aes_(y = ~runningScore), color= ES_linecolor, size = ES_linesize, linetype = ES_linetype) +
        theme_bw() +  
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.25,"cm"),
              panel.border = element_rect(size = p1bordersize, fill = NA),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = p1ytextsize, face = p1ytextfontface),
              axis.title.y = element_text(size = p1ytitlesize, face = p1ytitlefontface),
              plot.title=element_text(face=p1titlefontface, #字体
                                color=p1titlecolor, #颜色
                                size=p1titlesize,  #大小
                                hjust=p1titlehjust, #位置
                                vjust=p1titlevjust)
              ) +
        ylab("Enrichment score") + labs(title = Maintitle) +
        scale_x_continuous(expand=c(0,0), limits = c(0,length(gseadata@geneList))) + 
        scale_y_continuous(expand=c(0,0)) +
        annotate(geom = "text", x = NESx, y = NESy, size = NESsize, label = NES, color = NEScolor, fontface = NESfontface) +
        annotate(geom = "text", x = Pvaluex, y = Pvaluey, size = Pvaluesize, label = Pvalue, color = Pvaluecolor, fontface = Pvaluefontface) +
        annotate(geom = "text", x = FDRx, y = FDRy, size = FDRsize, label = FDR, color = FDRcolor, fontface = FDRfontface)

        v <- seq(1, sum(gsdata$position), length.out=9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) inv <- inv + 1
        col <- c(rev(brewer.pal(5, rec_col_low)), brewer.pal(5, rec_col_high))
        
    # part2
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
        geom_linerange(aes_(ymin=~ymin, ymax=~ymax), size = p2_linesize, linetype = as.numeric(p2_linetype), color = p2_linecolor) +
        xlab(NULL) + ylab(NULL) + theme_classic(11) +
        theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.title = element_blank(),
            axis.line.y = element_line(color = "white"),
            plot.margin = margin(t=-.1, b=0,unit="cm"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank()) +
        scale_x_continuous(expand=c(0,0),limits = c(0,length(gseadata@geneList))) +
        scale_y_continuous(expand=c(0,0)) 

    ymin <- min(p2$data$ymin)
    yy <- (max(p2$data$ymax) - min(p2$data$ymin)) * rec_height_ratio + min(p2$data$ymin)
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                            xmin = xmin,
                            xmax = xmax,
                            col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin=~xmin,
                            xmax=~xmax,
                            ymin=~ymin,
                            ymax=~ymax,
                            fill=~I(col)),
                            data=d,
                            alpha=rec_height_alpha,
                            inherit.aes=FALSE) +
            scale_x_continuous(expand=c(0,0),limits = c(0,rev(xmax)[1])) 
        
    if (showCase) {
        Casex = length(gseadata@geneList) * Casex
        Casey = max(abs(gseadata@geneList)) * Casey
        Casesize = Casesize
        Casecolor = Casecolor
        Casefontface = Casefontface
    } else {
        Casesize = 0
        Casex = length(gseadata@geneList) * 2
        Casey = max(abs(gseadata@geneList)) * 0.6
        Casecolor = "white"
        Casefontface = Casefontface
    }

    if (showCon) {
        Conx = length(gseadata@geneList) * Conx
        Cony = max(abs(gseadata@geneList)) * Cony
        Consize = Consize
        Concolor = Concolor
        Confontface = Confontface
    } else {
        Consize = 0
        Conx = length(gseadata@geneList) * 2
        Cony = max(abs(gseadata@geneList)) * -0.6
        Concolor = "white"
        Confontface = Confontface
    }

    df2 <- p1$data 
    df2$y <- p1$data$geneList[df2$x]
    df2$col <- rep(col[unique(inv)], rev(table(inv)))

    p3 <- ggplot() + 
        geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color=df2$col) +
        theme_bw()  + coord_cartesian(expand = 0) +
        theme(axis.ticks.y = element_blank(),
                axis.line = element_line(colour = "black", size = p3bordersize),
                axis.ticks.length = unit(0.25,"cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(size = p3xtextsize, face = p3xtextfontface),
                axis.title.x = element_text(size = p3xtitlesize, face = p3xtitlefontface),
                axis.text.y = element_blank(),
                axis.title.y = element_text(size = p3ytitlesize, face = p3ytitlefontface),
                plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm")
                ) +
        scale_y_continuous(expand=c(0,0),limits = c(max(abs(gseadata@geneList))*-1.2,max(abs(gseadata@geneList))*1.2)) +
        scale_x_continuous(expand=c(0,0),limits = c(0,length(gseadata@geneList))) +
        xlab("Ranked in Ordered Dataset") + ylab("Ranked List Metric") +
        annotate(geom = "text", x = Casex, y = Casey, size = Casesize, label = casename, color = Casecolor, fontface = Casefontface) +
        annotate(geom = "text", x = Conx, y = Cony, size = Consize, label = conname, color = Concolor, fontface = Confontface)
      
    p.final = p1+p2+p3+plot_layout(nrow = 3, height = c(p1ratio,p2ratio,p3ratio))

    return(p.final)
}

get_gsdata <- function(x, geneSetID) {
    if (length(geneSetID) == 1) {
        gsdata <- enrichplot:::gsInfo(x, geneSetID)
        return(gsdata)
    } 
    
    lapply(geneSetID, gsInfo, object = x) |>
        yulab.utils::rbindlist()
}
