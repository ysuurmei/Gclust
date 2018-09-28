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

Markup : ![picture alt](https://github.com/ysuurmei/Gclust/tree/master/Article/Order_to_pancake.png "Title is optional")

Thus the focus of this methodology is to uncover the hidden intent behind a set of similar orders. And the effectiveness of the algorithm in clustering order sets will strongly depend on the premise that each order is related to some specific intent. The weaker this premise, the less applicable the algorithm will be. In the dataset for which this algorithm was first developed, this premise held quiet strongly. Unfortunately it holds less so for the Online Retail data, but it still serves to show the application of the methodology. 

###### The clustering problem ######
In an ideal situation people would only buy milk, flour and eggs if they want to bake pancakes, and we would relate similar shopping cart by exactly matching their contents. However, in reality we see that some people like to put bacon in their pancakes, whilst others might prefer apples (Yes I’m Dutch….). Still we can argue that both customers are in fact shopping for pancakes.

![alt text](https://github.com/ysuurmei/Gclust/Article/Bacon_Apple.png)

This is the core of the problem at hand, we want to be able to cluster together highly similar orders, whilst acknowledging that order contents can vary. Furthermore, we would like to know which of the items in the basket are key in characterizing the cluster. From a list of unknown market baskets we would like to end up with a classifier consisting of key items, which we subsequently can use to draw inference on the diffent clusters.

