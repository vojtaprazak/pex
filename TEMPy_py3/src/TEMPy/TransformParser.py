##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright 2010-2014 TEMPy Inventors and Birkbeck College University of London.
#                          The TEMPy Inventors are: Maya Topf, Daven Vasishtan, 
#                           Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph,
#                          Harpal Sahota
# 
# 
#     TEMPy is available under Public Licence.
#     
#     Please cite your use of TEMPy in published work:
#     
#     Vasishtan D, Topf M. (2011) J Struct Biol 174:333-343. Scoring functions for cryoEM density fitting.
#
#===============================================================================

from numpy import load, save,savez
class TransformParser:
	"""A class to read and save transformation matrices """
	def __init__(self):
		pass

	def load_matrix(self,matrixname,mmap_mode=None):
		"""
		Load an array(s) from .npy, .npz 
		
		Arguments:
			
			*matrixname*:
				.npy  matrix
				If the filename extension is .gz, the file is first decompressed
				(see numpy.load for more information)
			
			*mmap_mode*:
				default None (memory-map the file)
				It can be set with different mode: 
				'r','r+','w+','c' accordingly with numpy.load (see numpy.memmap for a detailed description of the modes)
				The file is opened in this mode:
					'r'	Open existing file for reading only.
					'r+'	Open existing file for reading and writing.
					'w+'	Create or overwrite existing file for reading and writing.
					'c'	Copy-on-write: assignments affect data in memory, but changes are not saved to disk. The file on disk is read-only.
				A memory-mapped array is kept on disk. However, it can be accessed and sliced like any ndarray. 
				Memory mapping is especially useful for accessing small fragments of large files without reading the entire file into memory.
				
		"""
		return load(matrixname,mmap_mode)

	def save_npy_matrix(self,file,arr):
		"""
		Save an array to a binary file in NumPy .npy format.
		
		Arguments:
			*file* 
				File or filename to which the data is saved. If file is a file-object, then the filename is unchanged. If file is a string, a .npy extension will be appended to the file name if it does not already have one.
			*arr*
				array_like. Array data to be saved.
		
		"""
		return save(file, arr)
	
	def save_npz_matrix(self,file):
		"""
		Save several arrays into a single file in uncompressed .npz format. (See numpy.savez for more information)
		"""
		return savez(file)
	
