Module imports can cause problems on HPC clusters witth distributed file systems. Therefor we provide a discussion here so that we are prepared when MOCCaPy goes HPC. 

Potential problems are well explained [here](https://rcs.ucalgary.ca/How_to_work_with_large_number_of_small_files). Some remarks:

- _Network Congestion: In distributed file systems, accessing small files can generate significant network traffic, particularly if the files are scattered across multiple nodes. High metadata traffic, due to frequent file access and updates, can congest the network, leading to performance degradation._ With Python in a HPC context, every process typically will at startup import all its modules which can easily multiply the number of reads by a few thousands ... 

some links:

- [File-access-optimization (HPC admin)](https://www.admin-magazine.com/Archive/2025/87/File-access-optimization-discovery-and-visualization?utm_source=AU)
- [python on LUMI](https://docs.lumi-supercomputer.eu/software/installing/python/)