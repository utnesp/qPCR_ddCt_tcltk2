#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
tryCatch(expr = library(tcltk2), error = function(e) {biocLite("tcltk2", suppressUpdates = T)})
tryCatch(expr = library(ggplot2), error = function(e) {biocLite("tcltk2", suppressUpdates = T)})
tryCatch(expr = library(ggthemes), error = function(e) {biocLite("tcltk2", suppressUpdates = T)})

## get filename
file <- tclvalue(tkgetOpenFile()) 

if (!nchar(file)) {
    tkmessageBox(message = "No file was selected!")
} else {
    tkmessageBox(message = paste("The file selected was", file))
}

cq_values <- read.csv(file, skip = grep("Well", readLines(file)) - 1)
cq_values <- cq_values[c("Sample.Name", "Detector", "Ct")]
cq_values$Ct <- as.numeric(gsub("Undetermined", 40, cq_values$Ct))

print(cq_values$Sample.Name[!duplicated(cq_values$Sample.Name)])
print(cq_values$Detector[!duplicated(cq_values$Detector)])

## Calculate average CT for reference genes
samples <- cq_values[!duplicated(cq_values$Sample.Name), ]
samples <- samples[grep("NTC", samples$Detector, invert = T), "Sample.Name"]

################################################################
################################################################
win1 <- tktoplevel()
guessed_comparison = (length(samples) / 2)
name <- tclVar(guessed_comparison)
win1$env$entName <-tk2entry(win1, width = "25", textvariable = name)
tkgrid(tk2label(win1, text = "Enter number of comparisons to make:", justify = "left"),
       padx = 10, pady = c(15, 5), sticky = "w")
tkgrid(win1$env$entName, padx = 10, pady = c(0, 15))

onOK <- function() {
    nr_comparisons <<- tclvalue(name)
    tkdestroy(win1)
}
win1$env$butOK <-tk2button(win1, text = "OK", width = -6, command = onOK)
tkgrid(win1$env$butOK, padx = 10, pady = c(5, 15))
tkbind(win1$env$entName, "<Return>", onOK)

tkfocus(win1)




for (i in 1:nr_comparisons) {
assign(paste("comp", i, sep = "_"), tk_select.list(as.character(samples), multiple = T, 
 title = paste("Select comparison", i)), envir = .GlobalEnv)
}

 if (nr_comparisons >= 2) {
 comparisons <- comparisons <- data.frame(rbind(comp_1, comp_2))
        if (nr_comparisons > 2) {
            for (i in 3:nr_comparisons) {
            comparisons <- rbind(comparisons, 
                                get(paste("comp_", i, sep = ""))
                                )
     }}
 
 } else {
     comparisons <- data.frame(Sample1 = comp_1[1], Sample2 = comp_1[2])
     }

row.names(comparisons) <- NULL
colnames(comparisons) <- c("Sample1", "Sample2")

reference_genes <- tk_select.list(as.character(levels(cq_values$Detector)), multiple = T, title = "Select reference genes")


    
    t <- cq_values[cq_values$Detector %in% reference_genes, ]
    
    Ref_CT <- data.frame(Sample.Name = samples)
    for (i in 1:length(samples) ) {
        Ref_CT$Ref_CT[i] <- 2^mean(log2(t[t$Sample.Name == samples[i], "Ct"]))
        Ref_CT$sd[i] <- sd(t[t$Sample.Name == samples[i], "Ct"])
    }
    
    cq_values <- cq_values[grep("NTC", cq_values$Detector, invert = T), ]
    cq_values_ref <- merge(cq_values, Ref_CT, by = "Sample.Name", all = T)
    
    cq_values_ref$dCT <- round(2^(cq_values_ref$Ref_CT - cq_values_ref$Ct), 6)
    
    genes <- levels(cq_values_ref$Detector)
    genes <- genes[!genes %in% reference_genes]
    genes <- genes[grep("NTC", genes, invert = T)]
    
    replicates <- nrow(cq_values_ref[cq_values_ref$Detector == genes[i] & cq_values_ref$Sample.Name == samples[i], ])
    ddct_res <- data.frame(Detector = genes)
    
    for (k in 1:length(genes)) {
        for (j in 1:nrow(comparisons)) {
            for (z in 1:replicates) {
                ddct_res[k, paste(comparisons[j, 2], "_vs_", comparisons[j, 1], "_rep", z, sep = "")] <- cq_values_ref[cq_values_ref$Detector == genes[k] & cq_values_ref$Sample.Name == as.character(comparisons[j, 2]), "dCT"][z] / cq_values_ref[cq_values_ref$Detector == genes[k] & cq_values_ref$Sample.Name == as.character(comparisons[j, 1]), "dCT"][z]
            }
        }
    }
    
    ddct_res[2:ncol(ddct_res)] <- apply(ddct_res[2:ncol(ddct_res)], c(1,2), function(x) ifelse(x < 1, -1/x, x))
    row.names(ddct_res) <- ddct_res$Detector
    ddct_res$Detector <- NULL
    
    t <- t(ddct_res)
    
    row.names(t) <- gsub("_rep1", "", row.names(t))
    row.names(t) <- gsub("_rep2", "", row.names(t))
    row.names(t) <- gsub("_rep3", "", row.names(t))
    row.names(t) <- gsub("_rep4", "", row.names(t))
    
    library(data.table)
    
    t <- data.table::melt(t)
    colnames(t) <- c("X1", "X2", "value")
    
    t <- t[!is.na(t$value), ]
    t$value <- abs(t$value)
    t$X1 <- gsub("IgG_vs_acH3", "acH3_vs_IgG", t$X1)
    ggplot(t, aes(X2, value)) + geom_boxplot() + geom_point(aes(colour = factor(t$X1), shape = factor(t$X1)), size = 4) + 
        labs(x = "", y = "Fold change") + guides(colour = guide_legend(NULL), shape = guide_legend(NULL)) + theme_base() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA, colour = "black", size = 1))
    
    ggsave("ddCT_plot.pdf", device = "pdf", path = dirname(file), dpi = 300, height = 210, width = 297, units = "mm")
    
    # Produce a table of means and sd
    ddct_means <- data.table::dcast(t, X1~X2, mean)
    ddct_sd <- data.table::dcast(t, X1~X2, sd)
    ddct_means_melt <- melt(ddct_means)
    # row.names(ddct_means_melt) <- paste(ddct_means_melt[1], ddct_means_melt[2], sep = "_")
    colnames(ddct_means_melt) <- c("comparison", "mean", "Detector")
    ddct_sd_melt <- melt(ddct_sd)
    # row.names(ddct_sd_melt) <- paste(ddct_sd_melt$X1, ddct_sd_melt$X2, sep = "_")
    colnames(ddct_sd_melt) <- gsub("value", "sd", colnames(ddct_sd_melt))
    ddct_means_melt <- merge(ddct_means_melt, ddct_sd_melt, by = "row.names")
    ddct_means_melt <- ddct_means_melt[c("comparison", "Detector", "mean", "sd")]
    ddct_means_melt <- ddct_means_melt[order(ddct_means_melt$Detector), ]
    
    write.table(ddct_means_melt, file = paste(dirname(file), "ddct_means.txt", sep ="/"), sep = "\t", quote = F)
    write.table(ddct_res, file = paste(dirname(file), "ddct_res.txt", sep ="/"), sep = "\t", quote = F)
    write.table(cq_values_ref, file = paste(dirname(file), "cq_values_dct.txt", sep ="/"), sep = "\t", quote = F)
