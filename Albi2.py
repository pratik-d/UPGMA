from bioservices.services import RESTService
import numpy as np
from Bio import Phylo
from copy import deepcopy
import Bio.Cluster
from random import random, choice
import sys
from itertools import groupby
from Bio import SeqIO
#import random
from Bio import AlignIO
from StringIO import StringIO
#__all__ = ["MUSCLE"]

#MUSCLE webservice as taken from the Bioservices website

class MUSCLE(RESTService):
  

  _url = "http://www.ebi.ac.uk/Tools/services/rest/muscle"
  def __init__(self, verbose=True):
    super(MUSCLE, self).__init__(name='MUSCLE', url=MUSCLE._url, verbose=verbose )
    self._parameters = None
    self._parametersDetails = {}

  def getParameters(self):
    
    request = self.url + "/parameters/"
    res = self.request(request)
    return res

  def _get_parameters(self):
    if self._parameters:
      return self._parameters
    else:
      res = self.getParameters()
      parameters = [x.text for x in res.getchildren()]
      self._parameters = parameters
    return self._parameters
  parameters = property(_get_parameters, doc=r"""Read-only attribute that returns a list of parameters. See :meth:`getParameters`.""")

  def parametersDetails(self, parameterId):
    
    if parameterId not in self.parameters:
      raise ValueError("Invalid parameterId provided(%s). See parameters attribute" % parameterId)

    if parameterId not in self._parametersDetails.keys():
      request = self.url + "/parameterdetails/" + parameterId
      res = self.request(request)
      self._parametersDetails[parameterId] = res

    try:
      # try to interpret the content to return a list of values instead of the XML
      res = [x for x in self._parametersDetails[parameterId].findAll("value")]
      res = [y[0].text for y in [x.findAll('label') for x in res] if len(y)]
    except:
      pass

    return res
    
  def run(self, format=None, sequence=None, tree="none", email=None):
    # There are compulsary arguments:
    if format==None or sequence==None  or email==None:
      raise ValueError("format, sequence and email must be provided")

    from easydev import checkParam

    # Here, we will check the arguments values (not the type)
    # Arguments will be checked by the service itself but if we can
    # catch some before, it is better

    # FIXME: return parameters from server are not valid
    checkParam(format, ['fasta','clw','clwstrict','html','msf','phyi','phys'])

    checkParam(tree, ['none','tree1','tree2'])

    # parameter structure
    params = {
            'format': format,
            'sequence': sequence,
            'email': email}
      
    request = self.url + "/run/"
    res = self.requestPost(request, params)
    return res

  def getStatus(self, jobid):
    requestUrl = self.url + '/status/' + jobid
    res = self.request(requestUrl, format="txt")
    return res

  def getResultTypes(self, jobid):
    if self.getStatus(jobid)!='FINISHED':
      self.logging.warning("waiting for the job to be finished. May take a while")
      self.wait(jobid, verbose=False)
    requestUrl = self.url + '/resulttypes/' + jobid
    res = self.request(requestUrl, format="xml")

    output = {}
    
    def myf(x):
      if len(x)==0: return ""
      else: return x[0].text

    descriptions = [myf(x.findall("description")) for x in res.getchildren()]
    identifiers = [myf(x.findall("identifier")) for x in res.getchildren()]
    mediaTypes = [myf(x.findall("mediaType")) for x in res.getchildren()]
    labels = [myf(x.findall("label")) for x in res.getchildren()]
    suffixes = [myf(x.findall("fileSuffix")) for x in res.getchildren()]

    for i, ident in enumerate(identifiers):
      output[ident] = {'label':labels[i], 'mediaType': mediaTypes[i],
                       'description':descriptions[i], 'fileSuffix':suffixes[i]}
    return output

  def getResult(self, jobid, resultType):
    if self.getStatus(jobid)!='FINISHED':
      self.logging.warning("waiting for the job to be finished. May take a while")
      self.wait(jobid, verbose=False)
    if self.getStatus(jobid) != "FINISHED":
      raise ValueError("job is not finished")
    requestUrl = self.url + '/result/' + jobid + '/' + resultType
    res = self.request(requestUrl, format=resultType)

    return res

  def wait(self, jobId, checkInterval=5, verbose=True):
    import sys
    import time

    if checkInterval<1:
      raise ValueError("checkInterval must be positive and less than minute")
    result = 'PENDING'
    while result == 'RUNNING' or result == 'PENDING':
      result = self.getStatus(jobId)
      if verbose:
        print >> sys.stderr, jobId, " is ", result
      if result == 'RUNNING' or result == 'PENDING':
        time.sleep(checkInterval)
    return result

#generates tree/tree file for alignments without bootstrap values
def generate_tree(deleted_positions, species,):
	for_tree = list(range(0,10))	#list with numbers 0 to 10 which will later be edited to produce a newick tree
	#print for_tree
	for i in range(0, len(deleted_positions)):
		#print deleted_positions[i][0]
		#print deleted_positions[i][1]
		
		#print for_tree
		pos1 = for_tree[deleted_positions[i][1]]
		pos2 = for_tree[deleted_positions[i][0]]
		
		for_tree.pop(deleted_positions[i][1])
		for_tree.pop(deleted_positions[i][0])
		
		#print pos1
		#print '%s\t%s\n' %(pos1, pos2)
		#del for_tree[for_tree.index(deleted_positions[i][0])]
		#del for_tree[for_tree.index(deleted_positions[i][1])]
		
		a = '(' + str(pos1) + ',' + str(pos2) + ')'	#abspeichern der beiden position in newick format
		#a = str(deleted_positions[i][0]) + ',' + str(deleted_positions[i][1])
		#print a
		tree_list.append(a)
		for_tree.append(a)

	tree = str(for_tree).strip("['").strip("']")
	#print tree
	
	
	#print species
	tree = tree.replace('0', species_name[0][0]).replace('1', species_name[0][1]).replace('2', species_name[0][2]).replace('3', species_name[0][3]).replace('4', species_name[0][4]).replace('5', species_name[0][5]).replace('6', species_name[0][6]).replace('7', species_name[0][7]).replace('8', species_name[0][8]).replace('9', species_name[0][9])
	
	#print tree

	aaa = open('D:\\UNI\\Semester5\\AlBi\\Praktikum\\Praktikum2\\make_tree' + str(counter) + '.dnd', 'w')
	aaa.write(tree)
	aaa.close()
	tree_file = Phylo.read("D:\\UNI\\Semester5\\AlBi\\Praktikum\\Praktikum2\\make_tree" + str(counter) + ".dnd", "newick")
	#print tree_file

	tree_file.rooted = True
	Phylo.draw_ascii(tree_file)

#generates tree/tree file for alignments with bootstrap values
def bootstrap_tree(deleted_positions, species):
	print tree_list
	
	for_tree = list(range(0,10))	#list with numbers 0 to 10 which will later be edited to produce a newick tree
	#print for_tree
	#print deleted_positions
	for i in range(0, len(deleted_positions)):
		#print deleted_positions[i][0]
		#print deleted_positions[i][1]
		
		pos1 = for_tree[deleted_positions[i][1]]
		pos2 = for_tree[deleted_positions[i][0]]
		
		for_tree.pop(deleted_positions[i][1])	#every single number from deleted_positions is taken out and made into a node, which at the end produces a tree format
		for_tree.pop(deleted_positions[i][0])
		
		#print pos1
		#print '%s\t%s\n' %(pos1, pos2)
		
		a = '(' + str(pos1) + ',' + str(pos2) + ')'
		#print '%s\n' %(a)
		
		for jj in range(0, len(tree_list)):
			if a == tree_list[jj]:
				bootstrap_clades[jj] += 1	#counter, which determines how often a node has been bootstrapped
					
		#a = str(deleted_positions[i][0]) + ',' + str(deleted_positions[i][1])
		#print a
		for_tree.append(a)
		#print for_tree
	tree = str(tree_list[-1])
	#print tree
	#print tree_new
	#print tree
	

	
def dist_matrix(pim_mat):
	lines = pim_mat.split('\n')[6:] #because the reply from MUSCLE webservice consist of 5 blank lines at the beginning we build our matrix without including those 5 lines
	#print pim_matrix
	#the name of the species for the respective sequences are extracted
	species = [line.split()[1] for line in lines if len(line) >= 10]
	species = [item.split('|')[2] for item in species]
	

	distance_matrix = []	#an empty list which will be filled with values to form our distance matrix
	for i in range(0, len(lines)-1):
		a = lines[i].split(' ')
		#print a
		aa = []	#pim matrix can be seperated/parsed by splitting the elements whenever a blank space, i.e ' ' occurs
		for line in a:
			if line != '':
				aa.append(line) #all elements that are not ' ' will be written into the aa array/list
				#print line


		aa_array = np.array(aa[2:len(aa)]) #the aa list starting from the second position is written into the aa_array
		#print aa_array
		#print float(aa[3])
		
		row = []	#forming an empty list row, where values will be added after conversion into distance matrix values, these values will later be added to the distance matrix empty list formed in line 156 of this code
		for element in aa_array:
			els = float(element)
			dist = (1- els/100) #calculating distance matrix from pim matrix, where each value in divided by 100 and then subtracted by 1
			#print dist
			row.append(dist)	#values for distance matrix are added to the empty list row
			#row.append(element)
		#print row
		distance_matrix.append(row[0:len(row)])	#values from row are added to the empty distance matrix list, this gives us a list of lists which acts like a matrix in python
	
	
	print distance_matrix
		
	#formation of an upper-triangular matrix for the implementation of the UPGMA algorithm, the values on the lower diagonal of the matrix will be replaced by a '-' character which acts as a string
	for i in range(0, len(distance_matrix)):
		for j in range(0, len(distance_matrix[i])):
			if j <= i:
				distance_matrix[i][j] = '-'
	

	#print distance_matrix
	
	original_matrix = deepcopy(distance_matrix) #saving a copy of the upper-triangular matrix as orginal_matrix
	
	
	cluster_references = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]	#represents the number of sequences in a cluster, initial values are all 1s, because every single sequence is a single cluster at the beginning
	#deleted_positions = []
	for kk in range(0, len(original_matrix)):
		#the program enters the if loop when there is only a single value left in the matrix with no minimum to be defined. It will be easier to understand the program by first going through the else loop (line 206) and then looking at the if loop (line 199)
		if (int(len(distance_matrix)) == 2):
			b = distance_matrix[0][1]
			pos_min = [(i, pos.index(b)) for i, pos in enumerate(distance_matrix) if b in pos]
			#print b
			#print pos_min	# the last remaining value is the new minimum whose index and the value itself will be printed
			deleted_positions.append(pos_min[0])
			break

		else:
			#obtaining the minimum value in the distance matrix und pos_min gives the position of the minimum value in the matrix
			b = 100
			for i in range(0, len(distance_matrix)):
				a = min(distance_matrix[i])
				if a < b:
					b = a
			#print b
			pos_min = [(i, pos.index(b)) for i, pos in enumerate(distance_matrix) if b in pos]
			print pos_min
			deleted_positions.append(pos_min[0])
			#print distance_matrix
			
			#deleting row in a distance matrix as required in the UPGMA algorithm
			if pos_min[0][0] > pos_min[0][1]:
				del distance_matrix[int(pos_min[0][0])]
				del distance_matrix[int(pos_min[0][1])]
			else:
				del distance_matrix[int(pos_min[0][1])]
				del distance_matrix[int(pos_min[0][0])]
			
			#deleting a columns from the distance matrix as required in the UPGMA algorithm
			if pos_min[0][0] > pos_min[0][1]:
				[i.remove(i[pos_min[0][0]]) for i in distance_matrix]
				[i.remove(i[pos_min[0][1]]) for i in distance_matrix]
				
			else:
				[i.remove(i[pos_min[0][1]]) for i in distance_matrix]
				[i.remove(i[pos_min[0][0]]) for i in distance_matrix]
			#print a
			
			
			#dzws is an empty list, in which values will be added, the values are the ones that will be added as a new column in the distance matrix as required by UPGMA
			dzws = []
			for i in range(0, len(distance_matrix)):
				Nx = cluster_references[(pos_min[0][0])]
				Ny = cluster_references[(pos_min[0][1])]
				
				#obtaining Dxw value, to use the value in the formula for calculating values of the new column in the matrix
				if original_matrix[i][(pos_min[0][0])] != '-':
					Dxw = original_matrix[i][(pos_min[0][0])]
					#print Dxw
				elif original_matrix[(pos_min[0][0])][i] != '-':
					Dxw = original_matrix[(pos_min[0][0])][i]
					#print Dxw
				else:
					Dxw = 0
				
				#obtaining Dyw value, to use the value in the formula for calculating values of the new column in the matrix
				if original_matrix[i][pos_min[0][1]] != '-':
					Dyw = original_matrix[i][pos_min[0][1]]
					#print Dyw
				elif original_matrix[pos_min[0][1]][i] != '-':
					Dyw = original_matrix[pos_min[0][1]][i]
					#print Dyw
				else:
					Dyw = 0
				
				#calculating dzw as required by UPGMA
				dzw = ((Nx * Dxw) + (Ny * Dyw))/(Nx + Ny)
				#every time a dzw value is calculated, it will be added to the list dzws
				dzws.append(dzw)
				
			#deleting values in the cluster list, because sequences combine to form a cluster and the values must be adjusted accordingly
			if pos_min[0][0] > pos_min[0][1]:
				del cluster_references[pos_min[0][0]]
				del cluster_references[pos_min[0][1]]
			else:
				del cluster_references[pos_min[0][1]]
				del cluster_references[pos_min[0][0]]
			
			#updating new values to the cluster, which will be used in further processing of the algorithm
			cluster_references.append(Nx+Ny)	
			#print dzws
			
			#adding a new column to the matrix with appropriate values
			for i in range(0, len(distance_matrix)):
				distance_matrix[i].append(dzws[i])
			
			#print len(distance_matrix[0])
			
			#adding a new row to the matrix, albeit an empty row. This is because we always work with a square matrix and the square property is to be conserved
			d = ['-']*len(distance_matrix[0])
			distance_matrix.append(d)
			#print distance_matrix
			
			
	#print cluster_references
		
	#print distance_matrix
	#print species
	species_name.append(species)
	#print species[0], species[1]
	#species.append(speciess)
	
	#print deleted_positions


	
#runs MUSCLE for original/mutated sequences without the bootstrap values, also bootstrap alignments are created which later help in bootstrapping of sequences, such alignments are stored in p2.fasta file	
def Muscle_anwendung(sequences):
	k = MUSCLE()	#runs the class MUSCLE which has the code for the muscle webservice
	jobID = k.run(format='fasta',sequence=sequences.read(),email='pratikdhakal@hotmail.com')
	k.wait(jobID)
	pim = k.getResult(jobID, 'pim')
	pim_aln = k.getResult(jobID, 'aln-fasta')
	alignment_file = open(Pfad + 'p2.fasta', 'w')
	alignment_file.write(pim_aln)
	alignment_file.close()
	#print (k.getResultTypes(jobID))
	#pim_matrix = k.getResult(jobID,"pim",) #pim stands for "Percent Identity Matrix"
	dist_matrix(pim)	#funktionsaufruf dist_matrix

#runs MUSCLE for bootstrapped alignments, i.e extraction of pim matrix from alignments in order to generate a distance matrix
def Muscle_bootstrapped(bootstrapped_sequences):
	k = MUSCLE()	#runs the class MUSCLE which has the code for the muscle webservice
	jobID = k.run(format='fasta',sequence=bootstrapped_sequences.read(),email='pratikdhakal@hotmail.com')
	k.wait(jobID)
	pim1 = k.getResult(jobID, 'pim')
	#alignment_file = open(Pfad + 'p2.fasta', 'w')
	#alignment_file.write(pim1)
	#alignment_file.close()
	#print pim1
	dist_matrix(pim1)



	
def bootstrapping(Datei):
	handle = open(Datei, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	#print records[0].id
	#print len(records[0].seq)  #first record
	a = Pfad + 'bootstrap.fasta'		#bootstrapped sequences are written to the file bootstrap.fasta
	boot_file = open(a, 'w')
	
	
	listed = []
	for ppp in range(0, len(records[0].seq)):
		import random
		ppp = random.randrange(0, len(records[0].seq))	#generating random numbers for a list named listed
		listed.append(ppp)
	#print len(listed)
	#print listed
	
	#arranging the columns in the MSA according to the random list. This ultimately leads to a bootstrapped MSA
	for me in range(0,len(records)):
		id_name = records[me].id
		z = ''
		#boot_file.write('>' + id_name + '\n')
		for i in range(0, len(listed)):
			kk = records[me].seq[listed[i]]
			z = z + kk
		boot_file.write('>' + id_name + '\n')
		boot_file.write(z + '\n')
			
	boot_file.close()
	



#for mutation purpose the code was taken from the internet and modified slightly to adjust with the running code
def to_mutate(fasta_name, mutation_freq):
	def fasta_iter(fasta_name):
		"""
		given a fasta file. yield tuples of header, sequence
		"""
		fh = open(fasta_name)
		# ditch the boolean (x[0]) and just keep the header or sequence since
		# we know they alternate.
		faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
		for header in faiter:
			# drop the ">"
			header = header.next()[1:].strip()
			# join all sequence lines to one.
			seq = "".join(s.strip() for s in faiter.next())
			yield header, seq

	def main(fasta_name, mutation_freq):

		for header, seq in fasta_iter(fasta_name):
			seq = list(seq)
			for i, s in enumerate(seq):
				val = random()
				if val < mutation_freq:
					# choose a random nucleotide that's different.
					seq[i] = choice([x for x in "ACTG" if x != s.upper()])
			a.write(">%s\n%s\n" % (header, "".join(seq)))
			
	if __name__ == "__main__":
		main((fasta_name), fasta_name)

	
Pfad = 'D:\\UNI\\Semester5\\AlBi\\Praktikum\\Praktikum2\\'

sequences = open(Pfad + 'praktikum2_sequences.fasta','r')

counter = 0
deleted_positions = []
species_name = []
tree_list = []

Muscle_anwendung(sequences)
#print species_name
print '%s' %(deleted_positions)
generate_tree(deleted_positions, species_name)
#print deleted_positions
print tree_list
original_tree = tree_list[-1]

bootstrap_clades = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#for i in range(0,10):
to_open = Pfad + 'p2.fasta'
for i in range(0,100):
	deleted_positions = []
	#tree_list
	bootstrapping(to_open)
	boot_strap = open(Pfad + 'bootstrap.fasta', 'r')
	Muscle_bootstrapped(boot_strap)
	bootstrap_tree(deleted_positions, species_name)

tree_list = [tr.replace('0', species_name[0][0]).replace('1', species_name[0][1]).replace('2', species_name[0][2]).replace('3', species_name[0][3]).replace('4', species_name[0][4]).replace('5', species_name[0][5]).replace('6', species_name[0][6]).replace('7', species_name[0][7]).replace('8', species_name[0][8]).replace('9', species_name[0][9]) for tr in tree_list]
	
tree = str(tree_list[-1])
for i in range(0, len(tree_list)):
	if tree_list[i] in tree:
		#a = tree_list[i].strip(")")
		tree = tree.replace(tree_list[i], (str(tree_list[i] + ":" + str(bootstrap_clades[i]))))
		tree_list = [t.replace(str(tree_list[i]), (str(tree_list[i] + ":" + str(bootstrap_clades[i])))) for t in tree_list]
#print bootstrap_clades
#print tree_list
#print tree

tree_file1 = open(Pfad + 'boot_tree.dnd', 'w')
tree_file1.write(tree)
tree_file1.close()




#a for statement in order to get 5 different mutated sequences with varying rates of mutation, the new sequences are saved as mutated1.fasta, mutated2.fasta, etc. Variable j provides us with the rate of mutation
for pp in range(1,6):		
#Mutations:
	j = 0.05
	a = open(Pfad + 'mutated' + str(pp) + '.fasta', 'w')
	to_mutate((Pfad + 'praktikum2_sequences.fasta'), j)
	j += 0.05
	a.close()


deleted_positions = []
tree_list = []

#executing the function for mutated sequences
for pp in range(1,6):
	deleted_positions = []
	species_name = []
	tree_list = []
	counter = 1
	a = open(Pfad + 'mutated' + str(pp) + '.fasta', 'r')
	Muscle_anwendung(a)
	#print species_name
	#print '%s' %(deleted_positions)
	generate_tree(deleted_positions, species_name)
	#print deleted_positions
	#print tree_list
	original_tree = tree_list[-1]
	bootstrap_clades = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	to_open = Pfad + 'p2.fasta'
	
	for i in range(0, 100):
		deleted_positions = []
		bootstrapping(to_open)
		boot_strap = open(Pfad + 'bootstrap.fasta', 'r')
		Muscle_bootstrapped(boot_strap)
		bootstrap_tree(deleted_positions, species_name)
	
	#the numbers in the newick tree will be replaced through the name of the respective sequences
	tree_list = [tr.replace('0', species_name[0][0]).replace('1', species_name[0][1]).replace('2', species_name[0][2]).replace('3', species_name[0][3]).replace('4', species_name[0][4]).replace('5', species_name[0][5]).replace('6', species_name[0][6]).replace('7', species_name[0][7]).replace('8', species_name[0][8]).replace('9', species_name[0][9]) for tr in tree_list]
	
	tree = str(tree_list[-1])
	for i in range(0, len(tree_list)):
		if tree_list[i] in tree:
			#a = tree_list[i].strip(")")
			tree = tree.replace(tree_list[i], (str(tree_list[i] + ":" + str(bootstrap_clades[i]))))		#editing the tree to add bootstrap values to the branches
			tree_list = [t.replace(str(tree_list[i]), (str(tree_list[i] + ":" + str(bootstrap_clades[i])))) for t in tree_list]
			#print tree
			

	#print bootstrap_clades
	#print tree_list
	#print tree
	
	tree_file = open(Pfad + 'mutated_bootstrap' + str(pp) + '.dnd', 'w')	#bootstrapped mutated trees 
	tree_file.write(tree)
	tree_file.close()
	
		
	counter += 1

