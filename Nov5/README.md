# DAWG (Data Analysis Working Group) workshop series 
## 2021 Fall semester, 2nd workshop "Diversity analysis & differential abundance test"

### Getting started

Welcome back! Today, we are going to go over the basic (but important) diversity analysis of the microbiome data (16s rRNA), and differential abundance test using [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and [Metacoder](https://grunwaldlab.github.io/metacoder_documentation/). [Ray](https://plantpath.psu.edu/directory/rog5265) from Bull lab will lead the workshop today. 

#### Ready for the workshop

Before the workshop, I would like to recommend you to either i) install R (>4.0.3) in your local computer or ii) be able to launch R studio session in [Roar](portal.aci.ics.psu.edu). 

#### Launch R studio session in Roar

1. Connect [ICS-ACI portal](portal.aci.ics.psu.edu)
2. Log in with PSU account name and password
3. Click "My Interactive Sessions"
4. Click "R studio server"
5. Set up the desktop setting (Allocation, Number of hours, Node type)
6. Launch the session 

![Screen Shot 2021-11-04 at 1 46 04 PM](https://user-images.githubusercontent.com/77017866/140392070-b94465cc-d17c-48b8-96c1-a239d99ed6e0.png)

  
 * it may take a few minutes to initiate the session

#### Download required file into your directory

In this workshop, you will need two files

 ``` 
  i) .r - R script file
  ii) .rds - R object file (phyloseq object)
```

to download those in your scratch directory, click "terminal" in your R studio (You can find this on your left bottom screen)
then follow the lines below

```
cd gpfs/scratch/[your PSU ID]/
mkdir -p DAWG2021F-2nd
cd DAWG2021F-2nd

wget https://github.com/DAWG-psu/DAWG2021F/raw/main/Nov5/ps_DAWG_2021_F.rds
wget https://github.com/DAWG-psu/DAWG2021F/raw/main/Nov5/Garcia_16S_DAWG_Workshop_November_5.R
```

Now, let's go back to "console", and open the script file (.r file)

Check where you can find the buttons
![Screen Shot 2021-11-04 at 2 04 55 PM](https://user-images.githubusercontent.com/77017866/140395777-e4b38279-8e1d-4b62-9780-d6b769fcf260.png)

Type this - you made the directory from the previous step (/gpfs/scractch/[your PSU ID]/DAWG2021F-2nd
![Screen Shot 2021-11-04 at 2 09 08 PM](https://user-images.githubusercontent.com/77017866/140395696-369e8306-6e2a-4e84-a928-9d99ea85ea6e.png)

Finally, you can load .rds file into your environment by typing this in your console
```
ps.bac <- readRDS("/gpfs/scratch/[your PSU ID]/DAWG2021F-2nd/ps_DAWG_2021_F.rds
```


# **Now we are ready to go!**
