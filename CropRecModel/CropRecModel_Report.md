# Crop Recommendation Model Report

Group Project: Applying Mathematical Concepts to Real-World Issues in Data Science for an academic module 

Aims:
Our model will calculate the suitability of different crops at a theoretical Farm A, focusing on differences in soil nutritional values, analysing the three-dimensional impact of Nitrogen (N), Phosphorus (P), Potassium (K) to provide personalised crop recommendation. This is calculated by Euclidean distance as a measure of suitability for different crops relative to soil data. After recommending the most suitable crops, the second part of our model calculates the optimum method through which crop yields can be effectively increased through the alteration of soil states and nutrition using additional fertilisers. Our final model element integrates cost as a factor in crop selection. Helping farmers to choose the correct crop for their soil helps to maximise expected yield, whilst finding the best way of optimising soil nutrient ratios using gradient descent helps farmers to make economically advantageous decisions.

Both involve the model and the statistical analysis done to better understand the data.


# Descriptive Statistics of the Crop Data

![image](https://github.com/user-attachments/assets/315dd818-007d-4219-975f-e60de3f62fc6)



We also considered the mean and s.d. required for each crop. The mean and
standard deviation of the environmental conditions were calculated, and an “optimal
growth range” was defined for each crop, laying the foundation for crop
recommendation. Five example crops are shown below:


![image](https://github.com/user-attachments/assets/885bcd81-9c92-461d-989a-e4418cc05e25)

Statistical analysis of crop growing conditions allows us to identify similarities
and differences between crops. By comparing crop requirements in terms of N, P , K and
other environmental variables, we can combine the model with deeper knowledge of
different crops, helping the model to recommend the most suitable crop.


# Mathematical Model

Our model will focus primarily on Nitrogen, Phosphorus, and Potassium as the
NPK ratio is a standardised figure which may serve as a reference point in precision
agriculture, with possible variation based on crop type, crop stage, and soil types.
Farmer A’s theoretical soil values were chosen as N = 120, P = 40, and K = 75.
From the descriptive statistics, this places the Nitrogen value at almost the upper limit
of the data set, the Phosphorus value below the mean but within one standard deviation,
and the Potassium value above the mean but within one standard deviation. This selection 
therefore examines an unbalanced and nutrient-deficient soil, identifying the
ideal crop that might thrive under challenging conditions.


# Model Assumptions

Our model assumes that soil nutrition and state and climate variables can all
directly affect crop growth, as well as that soil data provided by farmers is
representative and accurate. As there is only one data point for each individual crop, we
also assume that the nutritional, state, and climate for each crop is recorded at a
consistent point in the planting/harvest process.


# Crop Recommendation Model

The Euclidean formula measures the distances between soil nutrient values of
farmers and average soil nutrition for each crop, allowing the variables to interact in
three dimensions using the equation:
![image](https://github.com/user-attachments/assets/9e883d38-757c-4c0a-952d-682388b8e790)


The results of this calculation show the soil nutrient of the farmer closest to the ideal type of crop. The calculated distances for each crop indicate that watermelon is the closest to the input value, and apple is the furthest.

![image](https://github.com/user-attachments/assets/a342ea1b-83fb-4ffd-91a0-af7d1b87b165)
Three-dimensional scatter plot illustrating NPK distances between Farm A
and crop averages, providing insight into optimal crop selection for unfavorable soil conditions.

Running iterations of gradient descent is a simple and efficient way to
handle multi-variable optimisation problems. Its iterative approach makes it
computationally efficient, and it has scalability, meaning that if we chose to expand the
model we would be able to include more environmental variables. It allows us to find
the closest minimums, i.e., the nearest distance between each of N, P , and K, and crop
average values, presuming that each variable may be optimised individually.

Gradient descent optimisation helps farmers select the most suitable crops based
on soil nutrients (N, P , K) and environmental factors like temperature, humidity, pH, and
precipitation. First, global averages of these variables for all crops are calculated to
understand their growth needs. Then, the nutrient requirements for each crop are
compared with the farmer’s soil conditions by calculating the Euclidean distance; the
smaller the distance, the better the match. Based on this, the model recommends the
most suitable crop. For example, for Farm A, rice may be the best choice, improving
crop efficiency and sustainability.


# Further Integrating Fertiliser Costs as a Model Factor
Gradient descent optimization is particularly important due to the cost
implications of improving soil quality; fertilisers may improve one soil nutritional factor ,
all three, or even affect pH, but any changes have an economic cost . In order to enhance
the utility of our model, integrating the cost of fertiliser based on the model’s
recommendations helps to provide farmers with cost effective options.
![image](https://github.com/user-attachments/assets/7aa63db6-4f36-4e8d-8179-6d9213f088a4)

We created a bar chart of optimised planting recommendations for Farm A. The
different crop suitability scores are positive and the higher the score, the better the crop
is suited to the current soil conditions. By looking at the charts, farmers can get a better
idea of which crops are better adapted and cost-effective in their soil conditions.

![image](https://github.com/user-attachments/assets/13010702-4684-4791-a8c9-d3c6c0f5be2f)



# Model Assessment

The linear reduction in Euclidean
distance over iterations is
steady,demonstrating the model's
performance.
![image](https://github.com/user-attachments/assets/18b051be-266e-4a1c-8c50-345c6ced21a6)


Distance reduction over
time for both training and
testing data during gradient
descent optimisation. The graph
illustrates the distance between
the farmer's soil and the ideal
crop values at each iteration.
![image](https://github.com/user-attachments/assets/9513ea56-8752-4722-b2ef-93fd8c488e85)
Overall, this suggests that the optimisation is
effectively working with the provided crop data. However , there could be potential
underfitting or ineffective learning due to a lack of complexity in the data or the model.
As a result, the model may fail to generalise efficiently to new, unseen data beyond what
it has memorised from the training set.

The suggested future steps for improving and applying the model to other data
include working with a more complex model that has sufficient data points to better
train and optimise the model. Incorporating a wider variety of crops with NPK values
could reveal interesting trends and help address the underfitting issue. In addition,
steps can be taken to ensure the optimisation models learn rather than memorise, such
as integrating machine learning components.



# References

- Harvard.edu, 2023. Crop_recommendation.tab - Harvard Dataverse. Accessed 26 Nov. 2024.
- www.india.gov.in. (n.d.). Website of Ministry of Agriculture & Farmers Welfare| National
Portal of India. [online] Available at:
https://www.india.gov.in/website-ministry-agriculture-farmers-welfare
- Government of India, Department of Fertilizers, 2024. Nutrient Based Subsidy (NBS)
Notification for Rabi 2024-25. [online] Available at:
https://www.fert.nic.in/sites/default/files/What-is-new/NBS%20Notification%20Rabi
%202024-25.pdf [Accessed 12 December 2024].


