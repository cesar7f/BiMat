library(igraph)
source('newman_r_version.r')

#Example 1 - Calculating and plotting the modularity of a real Coda Whale matrix
-------------------------------------------------------------------------------------------------------
matrix <- read.csv('matrix_whales.txt',sep=' ',header = FALSE) # READ A MATRIX FILE
matrix <- as.matrix(matrix) #REQUIRED TO CONVERT TO MATRIX, OTHERWISE THE CODE PRESENT PROBLEMS
data_frame <- read.table("whales_list.txt",header=1)

g <- graph.adjacency(matrix,weighted=TRUE) # Create a graph from the ajacency matrix
g <- as.undirected(graph.data.frame(data_frame)) # Or create the graph from data frame

#leading is an structure with all the required values of modularity
leading <- modularity.leading.eigenvector.commmunities(matrix) #YOU CAN SEND A MATRIX
leading <- modularity.leading.eigenvector.commmunities(g) #OR A IGRAPH OBJECT

#modularity value
print(leading$modularity)

#membership values in the same order that the nodes appear in the matrix
print(leading$membership)

#number of modules either:
print(max(leading$membership))
print(leading$N)

#Plotting the matrix (you need to send the structure you got on line 8 of this code)
#The first block will from top-lef will correspond to module 1, the second to module 2, and so on. In
#this case we will only have two modules.
modularity.plot_modularity(leading)