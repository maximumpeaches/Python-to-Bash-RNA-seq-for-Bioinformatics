import random
import re
import os
import subprocess
import time

#You might want to edit these variables (without necessarily understanding the rest of the file)
note = 'setting threads to 5. should fail and send an email to me.' #insert anything you'd like to output to the note file to this variable
work_dir = '/work/m/maxwell9/'
input_file = '/home/m/maxwell9/RNAseq/input.txt'
email_to_send_output_to = 'maxwell.pietsch@gmail.com'

def development_available():
	"""
	If the development partition of the Research Computing cluster has more than a certain number of CPUs idle
	, usually more than about 10%, then this returns the #SBATCH flag needed to run
	the job on the development partition, otherwise it returns an empty string
	"""
	min_percent_idle = 10
	sinfo_command = 'sinfo -O partition,cpusstate'
	p = subprocess.Popen(sinfo_command.split(), stdout=subprocess.PIPE)
	out, err = p.communicate()
	print(sinfo_command + '\n' + out)
	o = out.split()
	dev_flag = ''
	for i in range(len(o)): # iterate through the output looking for the development partition
		if 'development' in o[i]: # the aiot info immediately follows the word 'development' and some spaces
			dev_aiot = o[i+1].split('/')
			dev_t = dev_aiot[3]
			print('total number of cpus in development partition: ' + dev_t)
			dev_i = dev_aiot[1]
			print('number of cpus idle in development partition: ' + dev_i)
			percent_idle = float(dev_i) / float(dev_t) * 100
			print('percentage of cpus idle in development partition: ' + str(percent_idle) + '%')
			if percent_idle >= min_percent_idle:
				# TODO need to consider whether the time limit will be less than 2 hours to submit it to the development partition. We'll need to add (unless this has already been added by the time you're reading this) a bit of code that determines how long the job should take given the number of FASTQ files we want to process and their size and the size of the reference genome.
				# We might also measure the availability of Circe and submit it to Circe if it is more available than the development partition, if we feel we need to.
				print('Submitting script to development partition.')
				return '\n#SBATCH --partition=development\n'
			else:
				print('Submitting script to Circe partition.')
				return ''

class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg
    def __str__(self):
    	return repr(self.expr + ' ' + self.msg)

def check_pathname(variablename, pathname):
	"""
	TODO this should really be a class that validates upon initalization
	Checks that a given pathname is close to what a pathname should look like and that it exists
	>>> check_pathname('genomefile', 'workmmaxwell9genome')
	Traceback (most recent call last):
		...
	InputError: "genomefile must contain '/'"
	>>> check_pathname('variablex', '/work/m/maxwell9/ genome.fasta/')
	Traceback (most recent call last):
		...
	InputError: 'variablex must not contain empty spaces'
	>>> check_pathname('x', '  /work/m/maxwell9/genome.fasta  ')
	'/work/m/maxwell9/genome.fasta'
	"""
	pathname = pathname.strip()
	if '/' not in pathname:
		raise InputError(variablename, 'must contain \'/\'')
	if len(pathname.split()) != 1:
		raise InputError(variablename, 'must not contain empty spaces and be nonempty')
	if not os.access(pathname, os.F_OK):
		raise InputError(variablename, 'file or directory does not exist')
	return pathname

class SlurmScript(object):
	
	echo_info = '\necho "slurm job id: $SLURM_JOB_ID\nslurm cluster name: $SLURM_CLUSTER_NAME\nslurm nodelist: $SLURM_JOB_NODELIST"\n'

	def __init__(self, obj, objtype):
		self.unique_id = obj.unique_id
		self.dev_flag = objtype.dev_flag
		self._output_file = None
		self._array_size = None
		self._num_nodes = None
		self._dependency_condition = None
		self._filename = None
	
	@property
	def job_name(self):
		return self._job_name
	
	@job_name.setter
	def job_name(self, value):
		self._job_name = job_name

	@property
	def job_time(self):
		return self._job_time
	
	@job_time.setter
	def job_time(self, value):
		self._job_time = value

	@property
	def output_file(self):
		if self._output_file: #this way if someone wants to set a custom output_file that is cool or if they don't want to bother then we'll just output the standard one
			return self._output_file
		else:
			return work_dir + self.job_name + self.unique_id

	@output_file.setter
	def output_file(self, val):
		val = check_pathname(val)
		self._output_file = val


	@property
	def text(self):
		assert(job_time and job_name)
		self._text = ('#!/bin/bash\n\n#SBATCH --workdir=' + work_dir + 
			'\n#SBATCH --job-name=' + self.job_name + self.unique_id
			+ self.dev_flag
			+ '\n#SBATCH --time=' + self.job_time)
		if self._output_file: # remember that we need to not define any output file for when we use the slurm array
			self.text += '\n#SBATCH --output=' + self.output_file
		self._text += self.echo_info
		if self._array_size:
			#TODO test that this works if there's only two reads in the directory. I mean, does the syntax `#SBATCH --array=1` make sense to SLURM?
			# if not I can rewrite it so that the array is just not used at all in the script
			if int(self.array_size) == 2:
				self.text += '\n#SBATCH --array=1'
			if int(self.array_size) > 2:
				self.text += '\n#SBATCH --array=1-' + str(self.array_size)
		if self._num_nodes:
			self.text += '\n#SBATCH -N ' + self.num_nodes
		if self._threads:
			self.text += '\n#SBATCH -n ' + self.threads
		if self._dependency_condition:
			self.text += '\n#SBATCH --dependency=' + self.dependency_condition + ':' + self.dependent_on_id
		return self._text

	@property
	def dependency_condition(self):
		return self._dependency_condition
	
	@dependency_condition.setter
	def dependency_condition(self, val):
		self._dependency_condition = val

	@property
	def threads(self):
		return self._threads
	
	@threads.setter
	def threads(self, val):

	@property
	def num_nodes(self):
		return self._num_nodes
	
	@num_nodes.setter
	def num_nodes(self, val):
		assert(type(int(val)) == int)
		self._num_nodes = val

	@property
	def array_size(self):
		return self._array_size
	
	@array_size.setter
	def array_size(self, val):
		assert(type(int(val)) == int)
		assert(val >= 2) # must have at least two reads
		self._array_size = val

	@text.setter
	def text(self, value):
		self._text += value

	@property
	def slurm_id(self):
		return str(self._slurm_id)
	
	@slurm_id.setter
	def slurm_id(self, val):
		val = val.strip()
		assert(type(int(val)) == int)
		self._slurm_id = str(val)

	@property
	def filename(self):
		return self._filename
	
	@filename.setter
	def filename(self, val):
		self._filename = check_pathname(val)

	def outputAndSubmit():
		if not self.filename:
			filename = work_dir + self.job_name + self.unique_id + '.sh'
		with open(filename, 'w') as bt_script:
			bt_script.write(self.text)
		p = subprocess.Popen(['sbatch', filename], stdout=subprocess.PIPE)
		out, err = p.communicate()	
		print(out.strip())
		print(err.strip())
		for i, c in enumerate(out): #find location of first digit in out. This should be the start of the SLURM ID. Note how we output the SLURM IDs to the notes file below.
			if c.isdigit():
				self.slurm_id = str(out[i:len(out)]).strip()
				break

class SeqJob(object):

	dev_flag = development_available()
	btbuild_script = SlurmScript()
	seq_script = SlurmScript()
	notes_script = SlurmScript()

	def __init__(self, unique_id):
		assert(type(int(unique_id)) == int)
		self.unique_id = str(unique_id)

	@property
	def threads(self):
		return self._threads
	
	@threads.setter
	def threads(self, threads = '4')
		if (type(int(threads)) != int):
			raise InputError('threads', 'must be an integer')
		try:
			if (int(threads) > 32 or int(threads) < 1):
				raise InputError('threads', 'should be less than 32 and greater than 0. 32 was chosen somewhat arbitrarily. feel free to change it. it is probably a bit high.')
		except InputError:
			print('Setting threads to 4 and continuing execution')
			threads = '4'
		self._threads = threads

	@property
	def gff(self):
		return self._gff

	@gff.setter
	def gff(self, value):
		value = check_pathname('gff', value)
		if '.gff' not in value:
			raise InputError('gff','must contain .gff')
		print('Gff file is ' + value)
		self._gff = value

	@property
	def genome_index_path(self):
		return self._genome_index_path

	@genome_index_path.setter
	def genome_index_path(self, genome_index_path):
		genome_index_path = check_pathname('genome_index_path', genome_index_path)
		if not ('.fa' in genome_index_path or '.fasta' in genome_index_path):
			raise InputError('genome_index_path', 'must contain .fa or .fasta')
		print('Genome index path is ' + genome_index_path)
		self._genome_index_path = genome_index_path
		if genome_index_path[len(genome_index_path) - 1] == '/': # if it ends with a forward slash, get rid of the last character
			genome_index_path = genome_index_path[:len(genome_index_path) - 1]
		index = genome_index_path.rfind('/')
		index += 1
		genome_index = ''
		while index < len(genome_index_path):
			genome_index += genome_index_path[index]
			index += 1
		if not ('.fa' in genome_index or '.fasta' in genome_index):
			raise InputError('genome_index', 'must contain .fa or .fasta')
		
		# now we have something that looks like genome.fasta and want to convert it to something that looks like genome
		gilist = genome_index.split('.')
		if len(gilist) != 2:
			raise InputError('genome_index', 'must be in the form \'oneword.anotherword\'')
		genome_index_base = gilist[0]
		self.genome_index_base = genome_index_base

	@property
	def sample_dir(self):
		return self._sample_dir
	
	@sample_dir.setter
	def sample_dir(self, sample_dir):
		sample_dir = check_pathname('sample_dir', sample_dir)
		#if the last character of the sample_dir is not a forward slash, then append a forward slash to the sample_dir
		if sample_dir[len(sample_dir) - 1] != '/':
			sample_dir = sample_dir + '/'
		print('Sample directory is ' + sample_dir)
		self._sample_dir = sample_dir

	def createAndSubmitBtScript(self):
		self.btbuild_script.job_name = 'bt'
		self.btbuild_script.job_time = '00:08:00'
		self.btbuild_script.text += '\ncp ' + self.genome_index_path + ' ' + work_dir + self.genome_index_base + self.unique_id + '.fa'
				+ '\n\nbowtie2-build ' + self.genome_index_path + ' ' + work_dir + self.genome_index_base + self.unique_id
		self.btbuild_script.outputAndSubmit()
		

	def createAndSubmitSeqScript(self):
		"""
		Writes and submits a file which contains the meat of the pipeline.
		"""
		self.seq_script.job_name = 'seq'
		self.seq_script.job_time = '01:55:00'
		print('Scanning ' + self.sample_dir + ' for filenames containing \'.fastq\'')
		fastqs = []
		for file in os.listdir(self.sample_dir):
			if '.fastq' in file:
				fastqs.append(file)
		if len(fastqs) == 0:
			raise InputError(sample_dir, 'contains no files containing .fastq')
		
		reads = find_matching_reads(fastqs)
		
		reads_arr = ''
		sarray_count = -2
		print('Creating SLURM array of pathnames to FASTQ files')
		for key in reads:
			sarray_count += 2
			reads_arr += '\nreads[' + str(sarray_count) + ']="' + self.sample_dir + key + '"'
			reads_arr += '\nreads[' + str(sarray_count + 1) + ']="' + self.sample_dir + reads[key] + '"'

		print('Outputting primary pipeline script and calling sbatch')
		self.seq_script.nodes = '1'
		self.seq_script.array_size = len(reads)
		self.seq_script.threads = threads
		self.seq_script.dependency_condition = 'afterok'
		self.seq_script.dependent_on_id = self.bt_script.slurm_id
		# note how we don't set the output file here with the -o flag, because doing that when
		# using a slurm array will output all the data from the different members to the same
		# file which makes it confusing to read. using $SLURM_ARRAY_TASK_ID with #SBATCH -o didn't work

		output_directory_suffix = dir_id.zfill(2) + '\'_\'$SLURM_ARRAY_TASK_ID'
		tophat_out = 'thout_' + output_directory_suffix
		cufflinks_out = 'clout_' + output_directory_suffix
		cuffquant_out = 'cqout_' + output_directory_suffix	
		
		self.seq_script.text += ('\ntophat2 -G ' + gff.strip()
			+ ' -o ' + tophat_out 
			+ ' -p ' + threads 
			+ ' ' + genome_prefix 
			+ ' ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 2))]} ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 1))]}\n')
		self.seq_script.text += ('\ncufflinks' 
			+ ' -p ' + threads 
			+ ' -o ' + cufflinks_out 
			+ ' ' + tophat_out + '/accepted_hits.bam\n')
		self.seq_script.text += ('\ncuffquant '
			+'-p ' + threads 
			+ ' -o ' + cuffquant_out + ' ' 
			+ gff + ' ' 
			+ tophat_out + '/accepted_hits.bam\n')
		self.seq_script.outputAndSubmit()
		
	def createAndSubmitNotesScript(self):
		"""
		Creates and submits a file with some notes that might be helpful for purposes such as debugging, writing new scripts, etc.
		"""
		# find the time at which the script is being submitted. This information will be output to the notes file.
		curr_time = str(time.localtime().tm_hour).zfill(2) + '.' + str(time.localtime().tm_min).zfill(2)
		curr_day = str(time.localtime().tm_mon).zfill(2) + '.' + str(time.localtime().tm_mday).zfill(2) + '.' + str(time.localtime().tm_year)
		self.notes_script.filename = work_dir + 'notes' + self.unique_id + 'at' + curr_time + 'on' + curr_day
		self.notes_script.text = ('\nJob with custom ID ' + self.unique_id 
				+ ' and slurm IDs ' + self.bt_script.slurm_id + ', ' + self.seq_script.slurm_id 
				+ ' was submitted at ' + curr_time + ' on ' + curr_day)
		self.notes_script.outputAndSubmit()
		
	def createAndSubmitPostScript(self):
		"""
		Creates and submits a script which does some cleanup, and notifies people via email about how a program with any errors encountered
		and emails them relevant information from the job.
		"""
		temp_file = work_dir + 'tmp.txt
		job_notok = ('#!/bin/bash\n\n#SBATCH --workdir=' + work_dir + 
					'\n#SBATCH --job-name=post' + dir_id + 
				'\n#SBATCH --output=' + work_dir + 'slurm-postjob' + dir_id 
				+ dev_flag 
				+ '\n#SBATCH --dependency=afternotok:' + slurm_rnaseq_id
				+ '\n#SBATCH --time=00:10:00\n' 
				+ echo_info 
				+ '\n# use grep to find if job exceeded allocated memory limit. if so take some actions.'
				+ '\ndate >> ' + temp_file
				+ '\nif grep -Fxq "Exceeded job memory limit" my_list.txt'
				+ '\nthen'
				+ '\necho "Job memory limit exceeded with -n ' + str(threads) + ' threads" >> tmp.txt'
				+ '\nelse'
				+ '\n# code if not found'
				+ '\nfi'
				+ '\nmail -s "Email from Rays\' automated RNA-seq pipeline $(date)" ' + email_to_send_output_to + ' < ' + temp_file
				+ '\nrm ' + temp_file)
		job_notok_filename = work_dir + 'postjob' + dir_id
		with open(post_job_filename, 'w') as cleaner_file:
			cleaner_file.write(post_job)
		p = subprocess.Popen(['sbatch', post_job_filename], stdout=subprocess.PIPE)
		out, err = p.communicate()	
		print(out.strip())

#bt2-build files end in .btw
#if there are files ending in .btw in reference folders, don't run. otherwise, run bt2-build.
#after tophat, what percentage of alignment, and how many reads. can find in alignments_summary.txt. email info to swamy.
#need to add in cufflinks, cuffquant, etc.
#from tophat most important for analysis afterwards is accepted_hits.bam. delete all tophat files except that one after finishing cufflinks, cuffdiff, etc.
#do like Tophat6.sh, up to at least cuffnorm, except skip cuffmerge

def find_matching_reads(fastqs):
	"""
	Iterate through a list of filenames looking for .fastq files that are identical
	except one contains 'R1' and the other contains 'R2'
	return a dictionary where the keys are R1s and the values are the matching R2s
	>>> find_matching_reads(['/shares/RF_293_R1_m.fastq','/shares/RF_293_R1_m.fastq','/shares/RF_754_R1_m.fastq','/shares/RF_754_R2_m.fastq'])
	{'/shares/RF_293_R1_m.fastq': '/shares/RF_293_R2_m.fastq', '/shares/RF_754_R1_m.fastq': '/shares/RF_754_R2_m.fastq'}
	>>>
	"""
	print('Matching .fastq R1 and R2 files that have otherwise identical filenames')
	reads = {} # dictionary to hold name of read files where key is read1 and value is read2.
	for f in fastqs:
		r1pos = f.find('R1')
		if r1pos != -1: # if we found R1 in the file name then
			prefix = f[:r1pos] 
			x = ''
			for j in fastqs:
				r2pos = j.find('R2')
				if r2pos != -1 and j[:r2pos] == prefix: # if we found R2 in a filename and it has the same prefix as filename containing R1, then
					if x != '':
						raise InputError(f, 'had more than one matching R2 file')
					x = j
			if x == '':
				raise InputError(f, 'had no matching R2 file')
			else:
				reads[f] = x
	return reads

if __name__ == '__main__':
	doctest.mod()
	work_dir = check_pathname('work_dir',work_dir)
	if work_dir[len(work_dir) - 1] != '/': #work_dir should end with '/' for later parts of script to function properly
		work_dir += '/'
	#create random number to add to filenames we output so we don't replace files from previous runs.
	rand_id = str(random.randint(1,9999))
	dev_flag = development_available()
	input_file = check_pathname('input_file', input_file)
	#Open a file containing the pathnames of the FASTQ files we need to input to Tophat.
	with open(input_file, 'r') as f:
		comment = f.readline()
		num_of_dirs = int(f.readline().strip())
		if num_of_dirs < 1:
			raise InputError('num_of_dirs', 'must be greater than or equal to 1')
		for dir_num in range(num_of_dirs):
			j = SeqJob(str(dir_num) + rand_id)
			print('Reading from file \'' + input_file)
			comment = f.readline()
			j.setSampleDir(f.readline())
			comment = f.readline()
			j.setGff(f.readline())
			comment = f.readline()
			j.setGenomeFasta(f.readline())
			comment = f.readline()
			j.setThreads(f.readline())
			j.createAndSubmitBtScript()
			j.createAndSubmitSeqScript()
			j.createAndSubmitNotesScript()
			j.createAndSubmitPostScript()
