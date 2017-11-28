Note: This a project that is work in progress. 

# Fast analysis of raw qPCR results - qPCR_ddCt_tcltk2

So you have been strugling with excel to get the qPCR results analysed? It took you one hour in generating the results demanded by your supervisor? I believe you! Then maybe this tool could help you getting the results ready in less than 2 minutes!
The program is easy to understand, and works well on OSX (and I guess UNIX as well). Windows user might have problems, but we can sort this out together if you need my help.

To start the software, download and open tcltk2_ddct.Rexec. You can download using the green "Clone or download" button in the upper right corner. If you save as zip, unzip the contents and open the folder. Click on the tcltk2_ddct.Rexec file, and it should start. I tested this on my own machine, and it works for me. If you have any problems, please send me an inquiry so that it can be fixed. You might have to associate the file type .Rexec to open up in a terminal window (command line / cmd in windows), if for example the file opens up in a text editor.

The naming convention of the executable file has it origins from the tlctk R package by [Philippe Grosjean](https://github.com/phgrosjean). For recipes on tcltk2, see [here](http://www.sciviews.org/recipes/tcltk/toc/).

# Input format requirements
The input csv file generated by our qPCR machine looks like this:

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/input.csv.file.png)


If your input csv file is different than the one you see above, then change these two lines of code to your needs:

```R
cq_values <- read.csv(file, skip = grep("Well", readLines(file)) - 1)
cq_values <- cq_values[c("Sample.Name", "Detector", "Ct")]
```
Note: Only the columns with information regarding sample name, detector (gene) and Cq (Ct) are needed. If you have a txt file rather than a csv file, you can change read.csv to read.delim. 

# Step-by-step guide
![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/Step-by-step_guide.png)


# Results generated
The ddCT results will be ready waiting for you in the same folder as your input csv file. 
The results should look like this:

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/ddCT.plot.png)

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/ddct_means.png)

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/ddCT_res.png)

![alt text](https://github.com/utnesp/qPCR_ddCt_tcltk2/blob/master/cq_values_dct.png)



Good luck and hope this helps :)
