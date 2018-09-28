# Sales Order Classification #
### A new perspective on the 'Market Basket' problem ###

###### Problem introduction ######

One of the most prevalent sources of data for companies in the retail and logistics industry is data on sales orders. Supermarkets, factories, logistics companies, clothing stores, high-tech companies, car dealerships, all these companies and many more each day receive and send out thousands of product orders. The data generated in this process is a potential treasure trove of information about your customer.

In this blogpost I will discuss the possibility of using network analysis and graphing theory to cluster sales order data. The original setting for which I designed this approach was a high-tech company interested in retrieving the different types of maintenance actions performed by their engineers. The only data I had on the maintenance actions was a list of parts and tools related to each individual maintenance action, similar to typical sales order data. As such the methodology is applicable to pretty much every dataset that relates order numbers to items. 

For the purposes of this blogpost I will make use of the ‘Online Retail’ dataset from the UCI machine learning repository (if you’re not familiar it’s a website listing all sorts of cool free-to-use datasets, and can be found here). The dataset contains a years’ worth of sales transactions from a UK based online retailer. Besides invoices and items the dataset also lists other data like unit price, customer ID and customer country code. 

###### Orders and Items ######
Before we start it is necessary to introduce the concept of orders and items. Comparing it to grocery shopping, the order is your grocery basket, and the items are the products you put in your grocery basket. You can put multiple different items in your basket, and different people can put the same item in their basket. Orders and items have a many to many relationship, one order can contain multiple items, and items can appear in multiple orders.

###### Representation of the problem ######
The objective of this exercise is to group together, or cluster, similar sales orders. The central premise to this approach is that behind each (or a significant amount) sales order there is an unknown intent from the purchaser. Let us again take the example of a person going to the supermarket. Imagine someone puts milk, eggs and flour in their shopping basket. It is than reasonable to assume someone intends to use these to bake pancakes. The same can hold for any set of order data. This perspective is contrary to most research on order data, which focusses on the relationship between individual items, rather than orders. 

![Baking pankcakes](./Article/Order_to_pancake.png?raw = true)

Thus the focus of this methodology is to uncover the hidden intent behind a set of similar orders. And the effectiveness of the algorithm in clustering order sets will strongly depend on the premise that each order is related to some specific intent. The weaker this premise, the less applicable the algorithm will be. In the dataset for which this algorithm was first developed, this premise held quiet strongly. Unfortunately it holds less so for the Online Retail data, but it still serves to show the application of the methodology. 

###### The clustering problem ######
In an ideal situation people would only buy milk, flour and eggs if they want to bake pancakes, and we would relate similar shopping cart by exactly matching their contents. However, in reality we see that some people like to put bacon in their pancakes, whilst others might prefer apples (Yes I’m Dutch….). Still we can argue that both customers are in fact shopping for pancakes.

![Baking pankcakes](./Article/Bacon_Apple.png?raw = true)

This is the core of the problem at hand, we want to be able to cluster together highly similar orders, whilst acknowledging that order contents can vary. Furthermore, we would like to know which of the items in the basket are key in characterizing the cluster. From a list of unknown market baskets we would like to end up with a classifier consisting of key items, which we subsequently can use to draw inference on the diffent clusters.

![Baking pankcakes](./Article/Clusters.png?raw = true)

####### Previous work in this field ###### 
A well-known field in data science that also aims to extract information from order data is called ‘Market Basket Analysis’ (also known as frequent itemset mining, or association rules mining). This field is aimed at uncovering the relations between individual products that frequently co-occur in an order. If we again take the example of supermarket data, an outcome could be that milk and eggs often occur together in an order. The central premise in Market Basket Analysis is the idea that co-occurrence of items in an order tells us something about the relationship between these two items (For some excellent articles on this topic refer to here, here and here). Ever wondered how Amazon determines their product recommendations based on what you’re looking at currently. Chances are their algorithms are largely based on market basket analysis
The difference between market basket analysis and sales order clustering is that market basket analysis aims at finding the relationship among items, where we are trying to find the relationship between orders. 

![Baking pankcakes](./Article/MarketBasket_vs_Clustering.png?raw = true)

###### Network analysis and community detection ######
When I first tried to tackle this problem I started out with trying to use more traditional clustering methodologies. Sadly I found that algorithms like hierarchical clustering, DBSCAN, PAM etc. didn’t churn out reasonable results. So I decided to take a slightly different, but related, approach to clustering the sales order data. Namely, using graph theory to construct a network. In such a network each represents an individual sales order, and each edge in the network represents the similarity between two respective orders. 
An important subdomain in graph theory is called community detection. Similar to clustering community detection is concerned with the grouping of similar objects. An excellent blog listing the pros and cons of using graph theory versus traditional clustering, and a great explanation of what community detection is can be found [here](https://blog.insightdatascience.com/graph-based-machine-learning-6e2bd8926a0).

###### Jaccard similarity/distance ######
One of the key determinants of this methodology is the way in which we determine the similarity between individual orders. This similarity (or distance) measure forms the basis of our network and indicates the strength of the relationship between the individual nodes in our network. Key in calculating this network of similarities is the way in which we view our individual orders on a mathematical level. 
I treated each order as an observation on a binary vector with size N (where N = # unique Items). This binary indicates 1 if a certain item is present in the order, and 0 otherwise. Using this representation we can calculate the Jaccard similarity index for each pair of orders. There is a whole range of different similarity/distance metrics for binary vectors, each with their own respective characteristics. The Jaccard similarity is arguable one of the most widely used. Most importantly the Jaccard similarity doesn’t assume any similarity based on non-present features. So an item not belonging to a certain order doesn’t increase the similarity coefficient. More information on similarity measures for binary vectors and an excellent explanation of why to use Jaccard can be found here.

![Alt text](./Article/Jaccard.png?raw = true)

###### Louvain algorithm ######
I will not go too deeply into the the theory behind the Louvain algorithm, but there are a few things I would like to adress. The Louvain algorithm is a greedy optimization algorithm that optimizes the modularity of community set in a given network. The modularity is a measure for the degree of connectedness within communities versus the degree of connectedness between communities. For the official definition I’d like to refer to here
The algorithm has a time complexity of O(n log(n)) meaning the runtime complexity increases quasi-linearly making it easy to scale up for use on large datasets. Though a major constraint on the performance of the algorithm is size of the distance matrix to cluster on, as this increases exponentially in size as the number of orders/items grows. The problem there is not so much in the runtime of the Louvain algorithm as that it can quiet quickly push the limits of your RAM, especially if your working from a regular office laptop (like me). 
