""" here i wil take the csv file with the fitnesses and just calculate the averages and standard deviations across the
iterations and and make it into a graph"""

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt

# reading the csv to pandas df
fitness_df = pd.read_csv('fitness.csv', header=None)
print(fitness_df)

# mean
mean_fit = fitness_df.mean()
print(mean_fit)
"""
# Plotting the DataFrame
plt.figure(figsize=(12, 8))  # Set the figure size
plt.scatter(mean_fit.index[::10], mean_fit.values[::10], label="Data", color="blue", s=1)

# Adding titles and labels
plt.title("mean average fitness", fontsize=16)
plt.xlim(0, mean_fit.index.max())
# plt.xticks(range(0, 4050, 50), fontsize=3)
plt.ylim(0, mean_fit.values.max())
# plt.yticks(range(0, 100, 1), fontsize=3)
plt.xlabel("Index", fontsize=12)
plt.ylabel("Values", fontsize=12)

# Add a grid for better readability
plt.grid(True, linestyle='--', alpha=0.6)

# Add a legend
plt.legend() """

# Display the plot
plt.show()
# plt.savefig("mean average fitness.png")

# # std
std_fit = fitness_df.std()
# Plotting the DataFrame
plt.figure(figsize=(12, 8))  # Set the figure size
plt.scatter(std_fit.index[::10], std_fit.values[::10], label="Data", color="blue", s=1)

# Adding titles and labels
plt.title("standard deviation of the fitness", fontsize=16)
plt.xlim(0, std_fit.index.max())
# plt.xticks(range(0, 4050, 50), fontsize=3)
plt.ylim(0, std_fit.values.max())
# plt.yticks(range(0, 100, 1), fontsize=3)
plt.xlabel("Index", fontsize=12)
plt.ylabel("Values", fontsize=12)

# Add a grid for better readability
plt.grid(True, linestyle='--', alpha=0.6)

# Add a legend
plt.legend()

# Display the plot
plt.show()

# plt.figure(figsize=(30,20))
# plt.plot(x, std_fit)
# plt.xlabel("iterations")
# plt.ylabel("standard dev of fitness score")
# plt.xticks(range(0, 4000, 200))
# plt.savefig("std dev average fitness.png")
