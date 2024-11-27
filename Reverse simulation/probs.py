import random
import math

total_possible_reactions = 100000
jcvi_reactions = 338
amount = 25
target_jcvi = 212 # reactions making up JCVI syn3A central, nucleotide and lipid metabolism

def simulate_draws_to_target(total_possible_reactions, jcvi_reactions, amount,  target_jcvi):
    captured_jcvi = set()
    draws_count = 0

    while len(captured_jcvi) < target_jcvi:
        draw = random.sample(range(total_possible_reactions), amount)
        captured_jcvi.update(num for num in draw if num < jcvi_reactions)
        draws_count += 1

    return draws_count

simulations = 15000
results = [simulate_draws_to_target(total_possible_reactions, jcvi_reactions, amount, target_jcvi) for _ in range(simulations)]
average_draws = sum(results) / len(results)
sets_of_10 = math.ceil(average_draws)
print(sets_of_10)
squared_diffs = [(x - average_draws) ** 2 for x in results]
population_variance = sum(squared_diffs) / len(results)
N = len(results)  # Number of observations
standard_error = (population_variance / N) ** 0.5
print(standard_error)