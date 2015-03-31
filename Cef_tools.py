class CEF_obj(object):
	def __init__(self):
		self.headers = 0
		self.row_attr = 0
		self.col_attr = 0
		self.rows = 0
		self.cols = 0
		self.flags = 0
		self.tab_fields = 7

		self.header_names = []
		self.header_values = []

		self.row_attr_names = []
		self.row_attr_values = []

		self.col_attr_names = []
		self.col_attr_values = []

		self.matrix = []

	# def __add__(self, other):
	# 	assert type(other) == CEF_obj, 'Cannot perform addition between CEF_obj and a different data type'
	# 	result = CEF_obj()

	# 	result.row_attr_names = self.row_attr_names +\
	# 	[i for i in other.row_attr_names if i not in self.row_attr_names]
	# 	result.row_attr_values = self.row_attr_values +\
	# 	[other.row_attr_values[i] for i, v in enumerate(other.row_attr_names) if n not in self.row_attr_names]
	# 	result.row_attr = len(result.row_attr_names)


	# 	self.col_attr = 0
	# 	self.rows = 0
	# 	self.cols = 0
	# 	self.flags = 0
	# 	self.tab_fields = 7

	# 	self.header_names = []
	# 	self.header_values = []

	# 	self.row_attr_names = []
	# 	self.row_attr_values = []

	# 	self.col_attr_names = []
	# 	self.col_attr_values = []

	# def __radd__(self,other):
	# 	self.__add__(self, other)

	def update(self):
		self.headers = len( self.header_names)
		self.row_attr = len(self.row_attr_names)
		self.col_attr = len(self.col_attr_names)
		self.rows = len(self.matrix)
		self.cols = len(self.matrix[0])
		self.tab_fields = max(7, self.cols + self.row_attr + 1)
		self.linestr = u'%s' + u'\t%s' * (self.tab_fields-1) + u'\n'

	def add_header(self, name, value):
		self.header_names.append(name)
		self.header_values.append(value)

	def add_row_attr(self, name, value):
		self.row_attr_names.append(name)
		self.row_attr_values.append(list(value))

	def add_col_attr(self, name, value):
		self.col_attr_names.append(name)
		self.col_attr_values.append(list(value))

	def set_matrix(self, matrix):
		for row in matrix:
			self.matrix.append(list(row))

	def readCEF(self, filepath, matrix_dtype = int):
		#Delete all the stored information
		self.__init__()
		#Start parsing
		with open(filepath, 'rb') as fin:
			# Read cef file first line
			self.header, self.row_attr, self.col_attr, self.rows,\
			self.cols, self.flags = fin.readline().rstrip('\n').rstrip('\r').split('\t')[1:7]
			self.header = int(self.header)
			self.row_attr = int( self.row_attr )
			self.col_attr = int(self.col_attr)
			self.rows = int(self.rows)
			self.cols = int(self.cols)
			self.flags = int(self.flags)
			self.row_attr_values = [[] for _ in xrange(self.row_attr)]
			# Read header
			for i in range(self.header):
				name, value = fin.readline().rstrip('\n').rstrip('\r').split('\t')[:2]
				self.header_names.append(name)
				self.header_values.append(value)
			# Read col attr
			for i in range(self.col_attr):
				line_col_attr = fin.readline().rstrip('\n').rstrip('\r').split('\t')[self.row_attr:]
				self.col_attr_names.append( line_col_attr[0] )
				self.col_attr_values.append( line_col_attr[1:] ) 
			#Read row attr and matrix
			self.row_attr_names += fin.readline().rstrip('\n').rstrip('\r').split('\t')[:self.row_attr]
			for _ in xrange(self.rows):
				linelist = fin.readline().rstrip('\n').rstrip('\r').split('\t')
				for n, entry in enumerate( linelist[:self.row_attr] ):
					self.row_attr_values[n].append( entry )
				self.matrix.append( [matrix_dtype(el) for el in linelist[self.row_attr+1:] ])


	def writeCEF(self, filepath):
		self.update()
		with open(filepath, 'wb') as fout:
			#Write cef file first line
			fout.write( self.linestr % ( ('CEF', unicode(self.headers), unicode(self.row_attr),\
				unicode(self.col_attr) , unicode(self.rows), unicode(self.cols), unicode(self.flags) ) +\
				 ('',) * (self.tab_fields - 7) ) )
			#Write header
			for i in range(self.headers):
				fout.write(self.linestr % ( (unicode( self.header_names[i]), unicode( self.header_values[i]) ) + ('',) * (self.tab_fields - 2) ))
			#Write col attributes
			for i in range( self.col_attr ):
				fout.write( self.linestr % ( ('',) * (self.row_attr) + (unicode( self.col_attr_names[i] ),) + tuple( unicode(el) for el in self.col_attr_values[i] ) ))
			#Write headers of row attributes
			fout.write( self.linestr % ( tuple(unicode(el) for el in self.row_attr_names) + ('',)*(self.tab_fields-self.row_attr) ) )
			#Write rows
			for i in range(self.rows):
				for j in range(self.row_attr):
					fout.write( unicode(self.row_attr_values[j][i]) + u'\t')
				fout.write(u'\t')
				fout.write(u'\t'.join( [unicode(el) for el in self.matrix[i]] ) )
				fout.write(u'\n')



