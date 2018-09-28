# Sales Order Classification
###### A new perspective on the 'Market Basket' problem

**Problem introduction**
One of the most prevalent sources of data for companies in the retail and logistics industry is data on sales orders. Supermarkets, factories, logistics companies, clothing stores, high-tech companies, car dealerships, all these companies and many more each day receive and send out thousands of product orders. The data generated in this process is a potential treasure trove of information about your customer.
In this blogpost I will discuss the possibility of using network analysis and graphing theory to cluster sales order data. The original setting for which I designed this approach was a high-tech company interested in retrieving the different types of maintenance actions performed by their engineers. The only data I had on the maintenance actions was a list of parts and tools related to each individual maintenance action, similar to typical sales order data. As such the methodology is applicable to pretty much every dataset that relates order numbers to items. 
For the purposes of this blogpost I will make use of the ‘Online Retail’ dataset from the UCI machine learning repository (if you’re not familiar it’s a website listing all sorts of cool free-to-use datasets, and can be found here). The dataset contains a years’ worth of sales transactions from a UK based online retailer. Besides invoices and items the dataset also lists other data like unit price, customer ID and customer country code. 


