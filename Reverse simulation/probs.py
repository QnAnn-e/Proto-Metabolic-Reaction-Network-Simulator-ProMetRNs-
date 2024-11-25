import random
import math

# Simulation parameters
total_numbers = 100000
desired_numbers = 338
draw_size = 10
target_desired = 101

def simulate_draws_to_target(total_numbers, desired_numbers, draw_size, target_desired):
    # Track captured unique desired numbers
    captured_desired = set()
    draws_count = 0

    while len(captured_desired) < target_desired:
        # Perform a draw of `draw_size` numbers
        draw = random.sample(range(total_numbers), draw_size)
        # Add desired numbers from this draw to the captured set
        captured_desired.update(num for num in draw if num < desired_numbers)
        # Increment the number of draws performed
        draws_count += 1

    return draws_count

# Run the simulation multiple times for reliability
simulations = 1000
results = [simulate_draws_to_target(total_numbers, desired_numbers, draw_size, target_desired) for _ in range(simulations)]

# Average number of draws needed
average_draws = sum(results) / len(results)
sets_of_10 = math.ceil(average_draws)

print(sets_of_10)
