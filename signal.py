import numpy as np
import matplotlib.pyplot as plt
from . import visualizer

class Signal:
	""" Overall Signal definition"""
	
	def __init__(self, samples=np.empty(0), center_frequency=0, sample_rate=1, bits=np.empty(0), symbols=np.empty(0), symbol_rate=1, modulation_format='', modulation_order=0):
		self.samples = samples
		self.center_frequency = center_frequency
		self.sample_rate = sample_rate
		self.bits = bits
		self.symbols = symbols
		self.symbol_rate = symbol_rate
		self.modulation_format = modulation_format
		self.modulation_order = modulation_order


	def plot_spectrum(self):
		visualizer.plot_spectrum(self.samples, self.sample_rate)
