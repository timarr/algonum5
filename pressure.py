import numpy as np

time = 1.0
pressure_air = 1.225

#Computes the speed
def speed(length):
    return length/time

#Computes the dynamic pressure of a given length
def dynamic_pressure(length):
    return (1/2)*pressure_air*(speed(length)**2)

#Computes the total pressure
def calculate_pressure(length_function, slice):
    length = length_function(slice[0], slice[slice.size - 1])
    return dynamic_pressure(length)
