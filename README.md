# Sales Order Classification #
### A new perspective on the 'Market Basket' problem ###

#### Problem introduction ####

One of the most prevalent sources of data for companies in the retail and logistics industry is data on sales orders. Supermarkets, factories, logistics companies, clothing stores, high-tech companies, car dealerships, all these companies and many more each day receive and send out thousands of product orders. The data generated in this process is a potential treasure trove of information about your customer.

In this blogpost I will discuss the possibility of using network analysis and graphing theory to cluster sales order data. The original setting for which I designed this approach was a high-tech company interested in retrieving the different types of maintenance actions performed by their engineers. The only data I had on the maintenance actions was a list of parts and tools related to each individual maintenance action, similar to typical sales order data. As such the methodology is applicable to pretty much every dataset that relates order numbers to items. 

For the purposes of this blogpost I will make use of the ‘Online Retail’ dataset from the UCI machine learning repository (if you’re not familiar it’s a website listing all sorts of cool free-to-use datasets, and can be found here). The dataset contains a years’ worth of sales transactions from a UK based online retailer. Besides invoices and items the dataset also lists other data like unit price, customer ID and customer country code. 

#### Orders and Items ####
Before we start it is necessary to introduce the concept of orders and items. Comparing it to grocery shopping, the order is your grocery basket, and the items are the products you put in your grocery basket. You can put multiple different items in your basket, and different people can put the same item in their basket. Orders and items have a many to many relationship, one order can contain multiple items, and items can appear in multiple orders.

#### Representation of the problem ####
The objective of this exercise is to group together, or cluster, similar sales orders. The central premise to this approach is that behind each (or a significant amount) sales order there is an unknown intent from the purchaser. Let us again take the example of a person going to the supermarket. Imagine someone puts milk, eggs and flour in their shopping basket. It is than reasonable to assume someone intends to use these to bake pancakes. The same can hold for any set of order data. This perspective is contrary to most research on order data, which focusses on the relationship between individual items, rather than orders. 

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Order_to_pancake.png)

Thus the focus of this methodology is to uncover the hidden intent behind a set of similar orders. And the effectiveness of the algorithm in clustering order sets will strongly depend on the premise that each order is related to some specific intent. The weaker this premise, the less applicable the algorithm will be. In the dataset for which this algorithm was first developed, this premise held quiet strongly. Unfortunately it holds less so for the Online Retail data, but it still serves to show the application of the methodology. 

#### The clustering problem ####
In an ideal situation people would only buy milk, flour and eggs if they want to bake pancakes, and we would relate similar shopping cart by exactly matching their contents. However, in reality we see that some people like to put bacon in their pancakes, whilst others might prefer apples (Yes I’m Dutch….). Still we can argue that both customers are in fact shopping for pancakes.

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Bacon_Apple.PNG)

This is the core of the problem at hand, we want to be able to cluster together highly similar orders, whilst acknowledging that order contents can vary. Furthermore, we would like to know which of the items in the basket are key in characterizing the cluster. From a list of unknown market baskets we would like to end up with a classifier consisting of key items, which we subsequently can use to draw inference on the diffent clusters.

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Clusters.png)

#### Previous work in this field #### 
A well-known field in data science that also aims to extract information from order data is called ‘Market Basket Analysis’ (also known as frequent itemset mining, or association rules mining). This field is aimed at uncovering the relations between individual products that frequently co-occur in an order. If we again take the example of supermarket data, an outcome could be that milk and eggs often occur together in an order. The central premise in Market Basket Analysis is the idea that co-occurrence of items in an order tells us something about the relationship between these two items (For some excellent articles on this topic refer to here, here and here). Ever wondered how Amazon determines their product recommendations based on what you’re looking at currently. Chances are their algorithms are largely based on market basket analysis
The difference between market basket analysis and sales order clustering is that market basket analysis aims at finding the relationship among items, where we are trying to find the relationship between orders. 

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/MarketBasket_vs_Clustering.png)

#### Network analysis and community detection ####
When I first tried to tackle this problem I started out with trying to use more traditional clustering methodologies. Sadly I found that algorithms like hierarchical clustering, DBSCAN, PAM etc. didn’t churn out reasonable results. So I decided to take a slightly different, but related, approach to clustering the sales order data. Namely, using graph theory to construct a network. In such a network each represents an individual sales order, and each edge in the network represents the similarity between two respective orders. 
An important subdomain in graph theory is called community detection. Similar to clustering community detection is concerned with the grouping of similar objects. An excellent blog listing the pros and cons of using graph theory versus traditional clustering, and a great explanation of what community detection is can be found [here](https://blog.insightdatascience.com/graph-based-machine-learning-6e2bd8926a0).

#### Jaccard similarity/distance ####
One of the key determinants of this methodology is the way in which we determine the similarity between individual orders. This similarity (or distance) measure forms the basis of our network and indicates the strength of the relationship between the individual nodes in our network. Key in calculating this network of similarities is the way in which we view our individual orders on a mathematical level. 
I treated each order as an observation on a binary vector with size N (where N = # unique Items). This binary indicates 1 if a certain item is present in the order, and 0 otherwise. Using this representation we can calculate the Jaccard similarity index for each pair of orders. There is a whole range of different similarity/distance metrics for binary vectors, each with their own respective characteristics. The Jaccard similarity is arguable one of the most widely used. Most importantly the Jaccard similarity doesn’t assume any similarity based on non-present features. So an item not belonging to a certain order doesn’t increase the similarity coefficient. More information on similarity measures for binary vectors and an excellent explanation of why to use Jaccard can be found here.

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Jaccard.png)

#### Louvain algorithm ####
I will not go too deeply into the the theory behind the Louvain algorithm, but there are a few things I would like to adress. The Louvain algorithm is a greedy optimization algorithm that optimizes the modularity of community set in a given network. The modularity is a measure for the degree of connectedness within communities versus the degree of connectedness between communities. For the official definition I’d like to refer to here
The algorithm has a time complexity of O(n log(n)) meaning the runtime complexity increases quasi-linearly making it easy to scale up for use on large datasets. Though a major constraint on the performance of the algorithm is size of the distance matrix to cluster on, as this increases exponentially in size as the number of orders/items grows. The problem there is not so much in the runtime of the Louvain algorithm as that it can quiet quickly push the limits of your RAM, especially if your working from a regular office laptop (like me). 

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/CommunityDetection.png)

#### Cutoff score & determining optimum ####
The Louvain algorithm was designed to maximize the modularity, and this does not necessarily mean it will churn out clusters of highly similar orders. Rather, one of the issues I encountered when trying to implement a clustering solution based on this algorithm is that it tends to link together multiple clusters because of some ‘weak connection’ between 2 nodes in clusters. In figure X.X I’ve made a visual representation of this issue. Here the black edges are strong connections (high Jaccard similarity) between nodes, and the red edges represent weak connections (low Jaccard similarity). Because the Louvain algorithm still treats a low similarity edge as a connection between 2 nodes it can happen that multiple clusters are grouped together even though this results in a major decrease in the overall within cluster similarity (i.e. it groups together clusters that don’t belong together). 
 The way I dealt with this issue is by setting all edges below a certain threshold value X to zero. This effectively means removing the low similarity edges from the network. This approach is comparable to setting the cutoff score in hierarchical clustering. 
 
![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Issue1.png)

The issue than becomes how to determine the optimal cutoff score. Though this isn’t an exact science I’ve written a function called ‘optimizer’ to visualize the effect on the clusters at different threshold values. The main parameter for the optimizer function is ‘minMember’, this is the the minimal required number of members for a cluster to be considered ‘relevant’. As the Louvain algorithm start by treating every node in the network as an individual cluster the clustering results would become tedious to interpret if every single small cluster would be incorporated. Taking this ‘minMember’ parameter the optimizer function iterates over a sequence of possible cutoff scores and outputs the % of orders clustered in a ‘relevant’ cluster, the average within cluster Jaccard similarity (a measure for how much the objects in a cluster resemble one another), and the number of relevant clusters. Example of possible output is displayed in figure X.X.

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/results1.png)

Here we can clearly see the tradeoff between within cluster similarity and the amount of orders that end up in relevant clusters. Also striking is that the number of clusters as function of the cutoff score forms a concave set. Up to a certain point the Louvain algorithm will include more observations by assigning them to new ‘relevant' clusters. After the tipping point the algorithm will start merging relevant clusters and assigning new observations to those clusters to increase the overall modularity. (Note that the average within cluster similarity appears to decrease linearly. Which is also to be expected using this approach as by varying the cutoff score we keep adding the less similar observations also in a linear fashion).

#### Results ####
The first thing we have to do after loading the data is determining the cutoff score we want to use to prevent the Louvain algorithm from ‘over clustering’. In this example I’ve chosen to maximize the number of relevant clusters, so a cutoff score of 0.35 with minMember = 7. 
Below an example of one such clusters is given. On the x-axis the individual orders are listed, with on the y-axis the various items. A colored square indicates that a certain item is in a given orders. The green and red indicates whether the item should be considered key for the cluster or not (based on a user-defined value). Besides info on order contents the plot also lists other info like customer country, and date on the x-axis labels. 

![Baking_pankcakes](https://github.com/ysuurmei/Gclust/blob/master/Article/Christmas.png)

What is immediately noticeable about this cluster is that all key items seem to have something to do with Christmas. Looking at the top x-axis this idea is reinforced as we can see that nearly all orders were done in the period September-November, right in time for the holiday season. Looking at other clusters I also noted that most orders contain bulks of different versions of similar items, lunch boxes of various colors, shopping bags, tableware. This leads me to believe that a considerable amount of the orders to this retailer comes from wholesalers, buying items in bulk. 
These insights can be used to optimize marketing strategy, offering discounts on buying items in bulk. Or used in a recommendation system for a web shop. Another extension could be to relate client numbers to the clusters to see if there are any clients that regularly place certain types of orders. Or you could perform seasonality analyses for certain order types, like in the Christmas example to optimize stocking strategies. 

#### R-Code + implementation ####
You can find an R script containing several key functions, visualizations and an example using the “Groceries” dataset in this GitHub repository. The example dataset contains ca. 10.000 individual purchases at a supermarket chain. The implementation focusses on a subset of these orders and succeeds in effectively clustering 30% of orders in significant clusters.  
The code lists several functions that are useful in implementing the Louvain algorithm and visualizing the results. Where possible I’ve used sparse matrices to increase computational efficiency and reduce the required RAM capacity. 
 
 #### Pro’s and con’s ####
The upside to using this integral approach to clustering service orders is the ability to study intent in orders. The original setting for which I designed this approach was a high-tech company interested in retrieving the different types of maintenance actions performed by their engineers. Arguably this was a more ideal setting than grocery shopping, since the premise of specific intent held for most orders. I would therefore recommend this approach, not a replacement for other analysis (like market basket) but rather as a complement to it. 
The major drawback of this approach, like in all clustering techniques, is that the resulting clusters are dependent on the chosen cutoff score. But from several tests it appears that the large prevalent clusters still are rather stable. Equally important this algorithm doesn’t necessarily subdivide each observation into relevant clusters. 
Another limitation to the proposed methodology is the size of the distance matrix, which increases exponentially for large numbers of orders and especially different item sets (large N). I’ve tackled this issue by using sparse matrices and Rcpp optimized functions wherever possible. (More info on the use of  sparse matrices in R can be found here)

#### Possibilities and extensions ####
Do please let me know if you have any comments/suggestions related to the proposed methodology. If you have any issues regarding the R code, please create a XXX via my Git. Furthermore I’m planning to do a series of these types of articles where I treat theory around a Data Science related with a coding example. So if you’re interested in any particular topic please let me know.
