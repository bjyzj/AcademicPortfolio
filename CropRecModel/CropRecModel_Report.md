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


























References:

- Harvard.edu, 2023. Crop_recommendation.tab - Harvard Dataverse. Accessed 26 Nov. 2024.
- www.india.gov.in. (n.d.). Website of Ministry of Agriculture & Farmers Welfare| National
Portal of India. [online] Available at:
https://www.india.gov.in/website-ministry-agriculture-farmers-welfare
- Government of India, Department of Fertilizers, 2024. Nutrient Based Subsidy (NBS)
Notification for Rabi 2024-25. [online] Available at:
https://www.fert.nic.in/sites/default/files/What-is-new/NBS%20Notification%20Rabi
%202024-25.pdf [Accessed 12 December 2024].


