# qPCR_ddCt_tcltk2

So you have been strugling with excel to get the qPCR results analysed? It took you one hour in generating the results demanded by your supervisor? I believe you! Then maybe this tool could help you in getting the results ready in 5 minutes!
The program is easy to understand, and works well on OSX (and I guess UNIX as well). Windows user might have problems, but we can sort this out together if you need my help.

To start the software, just double-click on the tktclk_ddct.R file. 

The results will be stored in the folder the CSV file was read from. 

The input csv file looks like this:

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/input.csv.file.png)


If your input csv file is different than the one you see above, then change these two lines of code to your needs:

```R
cq_values <- read.csv(file, skip = grep("Well", readLines(file)) - 1)
cq_values <- cq_values[c("Sample.Name", "Detector", "Ct")]
```


The results should look something like this (different colors correspond to different replicates):
![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/ddCT.plot.png)


The ddCT results will be ready waiting for you in the same folder as your input csv file.


Good luck! :)
