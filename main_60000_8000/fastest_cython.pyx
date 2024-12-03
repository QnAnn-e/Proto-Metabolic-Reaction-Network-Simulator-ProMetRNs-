import numpy as np

cpdef double fitness_calc(list m_nodes, dict sie_dict):
    cdef list system_sie = []
    cdef int item
    cdef double total_sie = 0.0
    cdef double system_fitness = 0.0

    for item in m_nodes:
        if item in sie_dict:
            system_sie.append(sie_dict[item])

    # Use numpy for faster summation
    total_sie = np.sum(system_sie)

    system_fitness = total_sie - len(m_nodes)

    return system_fitness
