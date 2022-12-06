# scPolyA

#### *A pipeline to repurpose scRNA-seq data for single cell polyadenylation analysis*


## I. Generate Peak Reference

### Step 1: Pre-process BAM

#### Step 1(a) (optional): Subset BAM
* This is an optional step that requires cell barcodes to create a subset BAM file containing reads from only the specified cells
* This is useful when the user wants to perform use scPolyA or APA analysis for a subset (e.g. epithelial subset)

#### Step 1(b): Deduplication
* Cellranger outputs BAM which still contains PCR duplicates and the pipeline removes them using UMI tools (add link)


#### Step 1(c): Filter BAM ##
* oligo dT priming methodss can lead to internal priming events and such reads need to be removed
* This pipeline uses polyAFilter (add github link) to remove these
* Users are encouraged to spend time optimizing parameters for this step

#### Step 1(d): Merge all BAM ##  
* Each operation upto this point was performed on per sample BAM files 
* All BAM files for the given dataset are merged 

#### Step 1(e): Split by strand ##  
* For downstream processes, the merged BAM is split by strand using samtools (add link)

### Step 2: Find peaks and polyA reads

#### Step 2(a): Find Peaks
* MACS2 is used to find peaks

#### Step 2(b): Find polya reads 
* A custom script (add link) from Github is used to find polyA tail containing reads

### Step 3: Intersect peaks with polyA reads 
* Using bedtools (add link), peaks are intersected with polyA reads
* Only those peaks supported by min 3 polyA reads are retained

### Step4: Create final peak reference 

#### Step 4 (a): Assign peaks to genes or Transcription Units (TU) 
* A 5kb extension at the 3' is created for genes to create Transcription Units (TU)
* The peaks are assigned to these TUs to annotate them
* This is done by adapting previously pubished bulk APA analysis code from TingLab (add link to repo)

#### Step 4(b): Clean up peak reference ##
* Remove peaks assigned to >1 TUs
* Filter peaks supported by <3 polya reads


## II. Generate Peak by Cell Matrix #####
*

## III. Perform APA analysis #####
