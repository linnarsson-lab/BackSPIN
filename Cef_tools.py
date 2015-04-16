# Copyright (c) 2015 Gioele La Manno and Sten Linnarsson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# This code contains a simple parser reader for Cef files
# typical usage to write a cef_file:
#
# cef = CEF_obj()
# cef.add_header('Mynote', 'This is my message')
# cef.set_matrix(your_array)
# cef.add_col_attr(attribute_name, attribute_values_list)
# cef.add_row_attr(attribute_name, attribute_values_list)
# cef.writeCEF('path/to/your/file.cef')
#
# to read:
# cef = CEF_obj()
# cef.readCEF('path/to/your/file.cef')
# your_array = cef.matrix
# attribute_values_list2 = cef.col_attr_values[2]
# ...

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

	def readCEF(self, filepath, matrix_dtype = 'auto'):
		#Delete all the stored information
		self.__init__()
		#Start parsing
		with open(filepath, 'rbU') as fin:
			# Read cef file first line
			self.header, self.row_attr, self.col_attr, self.rows,\
			self.cols, self.flags = fin.readline().rstrip('\n').split('\t')[1:7]
			self.header = int(self.header)
			self.row_attr = int( self.row_attr )
			self.col_attr = int(self.col_attr)
			self.rows = int(self.rows)
			self.cols = int(self.cols)
			self.flags = int(self.flags)
			self.row_attr_values = [[] for _ in xrange(self.row_attr)]
			# Read header
			for i in range(self.header):
				name, value = fin.readline().rstrip('\n').split('\t')[:2]
				self.header_names.append(name)
				self.header_values.append(value)
			# Read col attr
			for i in range(self.col_attr):
				line_col_attr = fin.readline().rstrip('\n').split('\t')[self.row_attr:]
				self.col_attr_names.append( line_col_attr[0] )
				self.col_attr_values.append( line_col_attr[1:] ) 
			#Read row attr and matrix
			self.row_attr_names += fin.readline().rstrip('\n').split('\t')[:self.row_attr]
			for _ in xrange(self.rows):
				linelist = fin.readline().rstrip('\n').split('\t')
				for n, entry in enumerate( linelist[:self.row_attr] ):
					self.row_attr_values[n].append( entry )
				if matrix_dtype == 'auto':
					if sum(('.' in i) or ('e' in i) for i in linelist[self.row_attr+1:]) != 0:
						matrix_dtype = float
					else:
						matrix_dtype = int
				try:
					self.matrix.append( [matrix_dtype(el) for el in linelist[self.row_attr+1:] ])
				except ValueError:
					print repr(el), ' is invalid'



	def writeCEF(self, filepath, matrix_str_fmt = '%i'):
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
				fout.write(u'\t'.join( [unicode(matrix_str_fmt % el) for el in self.matrix[i]] ) )
				fout.write(u'\n')



