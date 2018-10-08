import os as os
import time as time
import numpy as np
from pandas import DataFrame, HDFStore


# class SafeHDFStore(HDFStore):
# 	# credit: Pietro Battiston (https://stackoverflow.com/questions/22522551/pandas-hdf5-as-a-database/29014295#29014295)
#     def __init__(self, *args, **kwargs):
#         probe_interval = kwargs.pop("probe_interval", 1)
#         self._lock = "%s.lock" % args[0]
#         while True:
#             try:
#                 self._flock = os.open(self._lock, os.O_CREAT |
#                                                   os.O_EXCL |
#                                                   os.O_WRONLY)
#                 break
#             except FileExistsError:
#                 time.sleep(probe_interval)

#         HDFStore.__init__(self, *args, **kwargs)

#     def __exit__(self, *args, **kwargs):
#         HDFStore.__exit__(self, *args, **kwargs)
#         os.close(self._flock)
#         os.remove(self._lock)

class Database:
	def __init__(self, filename):
		self.keys = []
		with HDFStore(filename) as h:
			self._store = h
		self.n_records = {}
		return

	def put(self, key, df, format='t'):
		self._store.open()
		try:
			self._store.put(key, df, format=format)
			if key not in self.keys:
				self.keys.append(key)
			self.n_records[key] = len(df)
		except:
			self._store.close()
			raise Exception('Failed to put %s in %s' % (key, self._store.filename))
		self._store.close()
		return

	def get(self, key):
		if key in self.keys:
			self._store.open()
			try:
				df = self._store.get(key)
			except:
				df = DataFrame()
			self._store.close()
		else:
			df = pd.DataFrame()
		return df

	def append(self, key, records, **opts):
		if key in self.keys:
			self._store.open()
			try:
				assert isinstance(records, DataFrame)
				n_records = len(records)
				self.n_records[key] =  self.n_records[key] + n_records
				self._store.append(key, records, **opts)
			except:
				self._store.close()
				raise Exception('Failed to append %s to %s in %s' % (str(type(records)), key, self._store.filename))
			self._store.close()
		else:
			raise Exception('Key does not exist in %s' % self._store.filename)
		return

	def remove(self, key):
		if key in self.keys:
			self._store.open()
			try:
				self._store.remove(key)
				self.keys.remove(key)
			except:
				self._store.close()
				raise Exception('Failed to remove %s from %s' % (key, self._store.filename))
			self._store.close()
		return

	def reindex(self, key, format='t'):
		df = self.get(key)
		n_records = len(df)
		df.index = np.arange(n_records)
		self.put(key, df, format=format)
		return
