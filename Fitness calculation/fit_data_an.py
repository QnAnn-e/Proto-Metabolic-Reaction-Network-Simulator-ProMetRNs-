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

# # std
std_fit = fitness_df.std()
# Plotting the DataFrame
plt.figure(figsize=(12, 8))  # Set the figure size
plt.scatter(std_fit.index[::10], std_fit.values[::10], label="Data", color="blue", s=1)

# Adding titles and labels
plt.title("standard deviation of the fitness", fontsize=16)
plt.xlim(0, (std_fit.index.max()+ 100))
# plt.xticks(range(0, 8000, 100), fontsize=3)
plt.ylim(0, (std_fit.values.max()+ 10))
# plt.yticks(range(0, 140, 2), fontsize=3)
plt.xlabel("Iterations", fontsize=12)
plt.ylabel("Standard Deviation", fontsize=12)

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
plt.savefig("std dev average fitness.png")
