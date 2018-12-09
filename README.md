# Discrete DP for Solving FAM
This package contains the source code for solving the FAM probelm optimally in 2d described in: S. Zeighami and R. C.-W. Wong, “Finding average regret ratio minimizing set in database,” arXiv preprint: https://arxiv.org/abs/1810.08047

Usage
===========
a. Compilation

	g++ -o run DP.cpp
b. Execution

	./run PointsFile k n
Where PointsFile contains the all the points in the dataset, k is the size of the solution returned and n is the size of the dataset. The algorithm assumes a uniform distribution of linear utility funcitons in a two dimensional database. In the output file, you can see the results of the algorithm. 

About dataset:
The dataset should contain n points in 2 dimensions. Each point must be written in a separate line with its dimensions being separated by tabs.

Example:
The dataset contains 50 points in 2 dimensions. Use

	./run examplePoints.txt 2 50

to select 2 points from the 50 points in the 2-dimensional dataset. The program outputs the results to the file result.txt. 
