import pickle
from AttenuationFunctionsTesting import max_and_normalized_max_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

with open('near_field_df.pkl', 'rb') as pkl_file:
            df = pickle.load(pkl_file)
            print("Completed: Loaded data from the near_field_df.pkl file") 

max_and_normalized_max_plots(df)