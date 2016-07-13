import random
import re
import os
import subprocess
import time
import sys

#This file was written by Maxwell Pietsch for Dr. Rays Jiang's lab
#If you need any support please contact me at maxwell.pietsch@gmail.com
#Written for use with Python 2.7

#You might want to edit these variables (without necessarily understanding the rest of the file)
note_content = '' #insert anything you'd like to output to the note file to this variable
curr_time = str(time.localtime().tm_hour).zfill(2) + '.' + str(time.localtime().tm_min).zfill(2)
curr_day = str(time.localtime().tm_mon).zfill(2) + '.' + str(time.localtime().tm_mday).zfill(2) + '.' + str(time.localtime().tm_year)


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

class BtScript(object):
	def __init__(self):
		self.slurm = Slurm()	

	def output(self):
		"""
		Writes and submits a Bash file to SLURM. The file tells bowtie2-build to create the files which Tophat2 then takes as input.
		"""
		self.slurm.job_name = 'bt'
		self.slurm.job_time = '00:08:00'
		self.slurm.text = ('\ncp ' + genome_index_path + ' ' + output_dir + genome_index_base + '.fa' 
			+ '\n\n/usr/bin/time bowtie2-build ' + genome_index_path + ' ' + output_dir + genome_index_base)
		self.slurm.outputAndSubmit()

class SeqScript(object):
	def __init__(self):
		self.slurm = Slurm()

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
	def sample_dir(self):
		return self._sample_dir

	@sample_dir.setter
	def sample_dir(self, sample_dir):
		sample_dir = check_pathname('sample directory', sample_dir)
		#if the last character of the sample_dir is not a forward slash, then append a forward slash to the sample_dir
		if sample_dir[len(sample_dir) - 1] != '/':
			sample_dir = sample_dir + '/'
		print('Sample directory is ' + sample_dir)
		self._sample_dir = sample_dir

	def output(self, bt_script_slurm_id):
		"""
		Writes and submits a file which contains the meat of the pipeline.
		"""
		self.slurm.job_name = 'seq'
		self.slurm.dependency = bt_script
		if not self.slurm.hours or not self.slurm.minutes:
			print('\n\n Hours and minutes are not set for main portion of pipeline! Setting to default, but this may fail! \n\n')
			self.slurm.job_time = '08:00:00'
		else:
			seq_script.slurm.job_time = self.slurm.hours + ':' + self.slurm.minutes + ':00'
		print('Scanning ' + self.sample_dir + ' for filenames containing \'.fastq\'')
		fastqs = []
		for file in os.listdir(self.sample_dir):
			if '.fastq' in file:
				fastqs.append(file)
		if len(fastqs) == 0:
			raise InputError(self.sample_dir, 'contains no files containing .fastq')
		
		reads = find_matching_reads(fastqs)
		
		self.slurm.slurm_array = ''
		sarray_count = -2
		print('Creating SLURM array of pathnames to FASTQ files')
		for key in reads:
			sarray_count += 2
			seq_script.slurm.slurm_array += '\nreads[' + str(sarray_count) + ']="' + self.sample_dir + key + '"'
			seq_script.slurm.slurm_array += '\nreads[' + str(sarray_count + 1) + ']="' + self.sample_dir + reads[key] + '"'

		print('Outputting primary pipeline script and calling sbatch')
		self.slurm.nodes = '1'
		self.slurm.array_size = len(reads)
		#self.slurm.threads is defined in main
		self.slurm.dependency_condition = 'afterok'
		self.slurm.dependent_on_id = bt_script_slurm_id
		# note how we don't set the output file here with the -o flag, because doing that when
		# using a slurm array will output all the data from the different members to the same
		# file which makes it confusing to read. using $SLURM_ARRAY_TASK_ID with #SBATCH -o didn't work

		tophat_out = output_dir + 'thout_' + '$SLURM_ARRAY_TASK_ID'
		cufflinks_out = output_dir + 'clout_' + '$SLURM_ARRAY_TASK_ID'
		cuffquant_out = output_dir + 'cqout_' + '$SLURM_ARRAY_TASK_ID'
		
		self.slurm.text += ('\n/usr/bin/time tophat2 -G ' + self.gff
			+ ' -o ' + tophat_out 
			+ ' -p ' + self.slurm.threads 
			+ ' ' + genome_index_base
			+ ' ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 2))]} ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 1))]}\n')
		self.slurm.text += ('\n/usr/bin/time cufflinks' 
			+ ' -p ' + self.slurm.threads 
			+ ' -o ' + cufflinks_out 
			+ ' ' + tophat_out + '/accepted_hits.bam\n')
		self.slurm.text += ('\n/usr/bin/time cuffquant '
			+'-p ' + self.slurm.threads 
			+ ' -o ' + cuffquant_out + ' ' 
			+ self.gff + ' '
			+ tophat_out + '/accepted_hits.bam\n')
		self.slurm.outputAndSubmit()

class Note(object):
	def output(self, bt_script_slurm_id, seq_script, input_contents):
		"""
		Creates and submits a file with some notes that might be helpful for purposes such as debugging, writing new scripts, etc.
		"""
		# find the time at which the script is being submitted. This information will be output to the notes file.
		text = ('\nJob with custom ID '
				+ ' and slurm IDs ' + bt_script_slurm_id + ', ' + seq_script.slurm.slurm_id 
				+ ' was submitted at ' + curr_time + ' on ' + curr_day + ' It contained ' + seq_script.slurm.threads + ' threads and ' + seq_script.slurm.mem + ' megabytes of memory. Input file:\n' + input_contents)
		with open(output_dir + 'note', 'w') as file:
			file.write(text)

class EmailUser(object):
	"""
	see also https://docs.python.org/2.7/library/email-examples.html
	
	import smtplib
	from email.mime.text import MIMEText

	fp = open(textfile, 'rb')
	# Create a text/plain message
	msg = MIMEText(text)
	fp.close()

	# me == the sender's email address
	me = 'automated_rna_seq_Rays_Jiang'
	# you == the recipient's email address
	#you = 'swamyrakesh@mail.usf.edu'
	you = 'maxwell.pietsch@gmail.com'
	msg['Subject'] = 'A message'
	msg['From'] = me 
	msg['To'] = you

	# Send the message via our own SMTP server, but don't include the envelope header.
	s = smtplib.SMTP('localhost')
	s.sendmail(me, [you], msg.as_string())
	s.quit()
	"""
	def __init__(self):
		self.subject = None
		self.sender = 'autoRNAseqRaysJiang'
		#self.receiver = email_to_send_output_to
		self.filename = None

	def output(self):
		pass
		
class JobNotOk(object):
	def __init__(self):
		self.slurm = Slurm()
		self.email_user = EmailUser()

	def output(self):
		"""
		I apologize for the clutter, but I'm pasting some text here in case it becomes useful in the future: 
		'\n# use grep to find if job exceeded allocated memory limit. if so take some actions.'
				+ '\ndate > ' + temp_file
				+ '\nif grep -Fxq "Exceeded job memory limit" ~/RNAseq/test_file.txt'
				+ '\nthen'
				+ '\necho "Job memory limit exceeded with -n ' + str(seq_script_threads) + ' threads" >> tmp.txt'
				+ '\nelse'
				+ '\n# code if not found'
				+ '\nfi'
				+ '\nmail -s "Email from Rays\' automated RNA-seq pipeline $(date)" ' + email_to_send_output_to + ' < ' + temp_file
				+ '\nrm ' + temp_file
		"""
		self.email_user.filename = 'to_email_if_notok.py'
		self.email_user.subject = 'RNA-seq job failed'
		self.email_user.output()
	
		self.slurm.job_name = 'notok'
		self.slurm.job_time = '00:10:00'
		self.slurm.dependency_condition = 'afternotok'
		self.slurm.dependent_on_id = seq_script.slurm.slurm_id
		self.slurm.text = ('python ' + output_dir + self.email_user.filename)
		self.slurm.outputAndSubmit()

class JobOk(object):
	def __init__(self):
		self.slurm = Slurm()
		self.email_user = EmailUser()

	def output(self, seq_script):
		self.email_user.filename = 'to_email_if_ok.py'
		self.email_user.subject = 'RNA-seq job succeeded'
		self.email_user.output()
	
		self.slurm.job_name = 'ok'
		self.slurm.job_time = '00:10:00'
		self.slurm.dependency_condition = 'afterok'
		self.slurm.dependent_on_id = seq_script.slurm.slurm_id
		self.slurm.text = ('python ' + output_dir + self.email_user.filename)
		self.slurm.outputAndSubmit()

class JobAny(object):
	def __init__(self):
		self.slurm = Slurm()

	def output(self, seq_script):
		self.slurm.job_name = 'any'
		self.slurm.job_time = '00:10:00'
		self.slurm.dependency_condition = 'afterany'
		self.slurm.dependent_on_id = seq_script.slurm.slurm_id
		self.slurm.text = ('\nrm_arr=( "a" "b" )\nfor element in rm_arr\ndo\n	rm $element\ndone')
		self.slurm.outputAndSubmit()

class Slurm(object):
	
	echo_info = '\n\necho "slurm job id: $SLURM_JOB_ID\nslurm cluster name: $SLURM_CLUSTER_NAME\nslurm nodelist: $SLURM_JOB_NODELIST"\n'

	def __init__(self):
		self._output_file = None
		self._array_size = None
		self._num_nodes = None
		self._dependency_condition = None
		self._filename = None
		self._threads = None
		self._mem = None
		self._minutes = None
		self._hours = None
		self._text = ''
		self.slurm_array = None
		self.dev_flag = dev_flag
	
	@property
	def job_name(self):
		return self._job_name
	
	@job_name.setter
	def job_name(self, value):
		self._job_name = value

	@property
	def job_time(self):
		return self._job_time
	
	@job_time.setter
	def job_time(self, value):
		self._job_time = value

	@property
	def output_file(self):
		return self._output_file
		
	@output_file.setter
	def output_file(self, val):
		self._output_file = check_pathname(val)	

	@property
	def threads(self):
		return self._threads
	
	@threads.setter
	def threads(self, threads):
		threads = threads.strip()
		if (type(int(threads)) != int):
			raise InputError('threads', 'must be an integer')
		try:
			if (int(threads) > 8 or int(threads) < 1):
				raise InputError('threads', 'should be less than 8 and greater than 0. There (probably!) aren\'t significant performance increases for more than 8 threads.')
		except InputError:
			print('Setting threads to 4 and continuing execution')
			threads = '4'
		self._threads = threads

	@property
	def text(self):
		return self._text

	@text.setter
	def text(self, value):
		self._text = value

	@property
	def dependency_condition(self):
		return self._dependency_condition
	
	@dependency_condition.setter
	def dependency_condition(self, val):
		self._dependency_condition = val

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

	@property
	def mem(self):
		return self._mem

	@mem.setter
	def mem(self, val):
		val = val.strip()
		if (type(int(val)) != int):
			raise InputError('memory requested', 'must be an integer')
		self._mem = val
	
	@property
	def hours(self):
		return self._hours

	@hours.setter
	def hours(self, val):
		val = val.strip()
		if (type(int(val)) != int):
			raise InputError('hours requested', 'must be an integer')
		self._hours = val

	@property
	def minutes(self):
		return self._minutes

	@minutes.setter
	def minutes(self, val):
		val = val.strip()
		if (type(int(val)) != int):
			raise InputError('minutes requested', 'must be an integer')
		if (int(val) > 59 or int(val) < 0):
			raise InputError('minutes requested', 'must be greater than 0 and less than 59')
		self._minutes = val

	@property
	def slurm_id(self):
		return self._slurm_id
	
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
		self._filename = check_pathname('self._filename', val)

	def get_script_content(self):
		"""
		Returns the content we want to write to file, to then submit to SLURM.
		The idea behind the if statements, such as 'if self._threads:' is that in the init() function _threads was set to None. So it's like we're saying if we are specifying a number of threads, append that number to the script to be output. If we never defined the number of threads then skip appending anything.
		"""
		assert(self.job_time and self.job_name)
		self.script_content = ('#!/bin/bash\n\n#SBATCH --workdir=' + output_dir + 
			'\n#SBATCH --job-name=' + self.job_name
			+ '\n#SBATCH --time=' + self.job_time + self.dev_flag )
		if self._output_file: # remember that we need to leave the output_file variable as None when we include a slurm array
			self.script_content += '\n#SBATCH --output=' + self.output_file
		if self._array_size:
			#TODO test that this works if there's only two reads in the directory. I mean, does the syntax `#SBATCH --array=1` make sense to SLURM?
			# if not I can rewrite it so that the array is just not used at all in the script
			if int(self.array_size) == 2:
				self.script_content += '\n#SBATCH --array=1'
			if int(self.array_size) > 2:
				self.script_content += '\n#SBATCH --array=1-' + str(self.array_size)
		self.script_content += '\n#SBATCH --nodes=1'
		if self._threads:
			self.script_content += '\n#SBATCH --ntasks=' + self.threads
		if self._mem:
			self.script_content += '\n#SBATCH --mem=' + self.mem
		if self._dependency_condition:
			self.script_content += '\n#SBATCH --dependency=' + self.dependency_condition + ':' + self.dependent_on_id
		if self.slurm_array:
			self.script_content += '\n' + self.slurm_array
		self.script_content += self.echo_info
		self.script_content += self.text
		return self.script_content

	def outputAndSubmit(self):
		"""
		Submits content to SLURM. Returns SLURM id.
		"""
		if not self.filename:
			self.filename = output_dir + self.job_name + '.sh'
		print('Outputting file ' + self.filename)
		with open(self.filename, 'w') as script:
			script.write(self.get_script_content())
		p = subprocess.Popen(['sbatch', self.filename], stdout=subprocess.PIPE)
		out, err = p.communicate()	
		if out:
			print(out.strip())
			print(self.filename + ' submitted successfully')
			if err:
				# In the past, when there was an error, it printed to console without getting placed in this err variable. Not sure why.
				print(err.strip())
			for i, c in enumerate(out): #find location of first digit in out. This should be the start of the SLURM ID. Note how we output the SLURM IDs to the notes file below.
				if c.isdigit():
					self.slurm_id = str(out[i:len(out)]).strip()
					break
		else:
			print('Something went wrong with submission to SLURM')		
			if self.dev_flag != '': # don't want to create an endless loop. using dev_flag to tell us if we've already tried submitting to Circe
				print('Did you violate the `#SBATCH --qos=development` requirements? Submitting to Circe now in case that was the source of error.')
				self.dev_flag = ''
				self.outputAndSubmit()
			else:
				print('\n\n' + self.filename + ' was not submitted!\n\n')

def development_available():
	"""
	If the development partition of the Research Computing cluster has more than a certain number of CPUs idle
	, usually more than about 10%, then this returns the #SBATCH flag needed to run
	the job on the development partition, otherwise it returns an empty string
	"""
	min_percent_idle = 10
	sinfo_command = 'sinfo -O partition,cpusstate'
	try:
		p = subprocess.Popen(sinfo_command.split(), stdout=subprocess.PIPE)
	except Exception, e:
		#if I'm testing it on my computer at home, an exception will be raised
		#because my computer isn't the research computing cluster
		#but I'd like to be able to test it on my computer, so I'm returning an empty string in that case
		#also this practice seems fine and is perhaps an improvement on the cluster
		print(e)
		return ''
	else:
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
					return '\n#SBATCH --partition=development\n#SBATCH --qos=development'
				else:
					print('Submitting script to Circe partition.')
					return ''

def check_pathname(variablename, pathname):
	"""
	TODO maybe this should be a class that validates upon initalization
	Checks that a given pathname is close to what a pathname should look like and that it exists
	 check_pathname('genomefile', 'workmmaxwell9genome')
	Traceback (most recent call last):
		...
	InputError: "genomefile must contain '/'"
	>>> check_pathname('variablex', '/work/m/maxwell9/ genome.fasta/')
	Traceback (most recent call last):
		...
	InputError: 'variablex must not contain empty spaces and be nonempty'
	>>> check_pathname('x', '  /work/m/maxwell9/genome.fasta  ')
	'/work/m/maxwell9/genome.fasta'
	"""
	pathname = pathname.strip()
	#if '/' not in pathname:
	#	raise InputError(variablename, 'must contain \'/\'')
	if len(pathname.split()) != 1:
		raise InputError(variablename, 'must not contain empty spaces and be nonempty')
	#TODO I was getting some false positives when checking if a file or directory exists, so I've disabled it for now
	#if not os.access(pathname, os.F_OK):
	#	raise InputError(variablename, 'file or directory does not exist')
	return pathname



def extract_genome_index_base(genome_index_path):
	"""
	>>> f.process_genome_index_path('/home/makes/genomeindex.asta')
	Traceback (most recent call last):
		...
	InputError: 'genome_index_path must contain .fa or .fasta'
	>>> f.process_genome_index_path('/home/makes/genomeindex.fa')
	Genome index path is /home/makes/genomeindex.fa
	'genomeindex'
	>>> f.process_genome_index_path('/home/makes/genomeindex.fa/')
	Genome index path is /home/makes/genomeindex.fa/
	'genomeindex'
	>>> f.process_genome_index_path('/home/makes/genomeindex.fasta')
	Genome index path is /home/makes/genomeindex.fasta
	'genomeindex'
	>>> f.process_genome_index_path('/home/makes/.fasta')
	Traceback (most recent call last):
		...
	InputError: 'genome_index_path cannot be a file that consists solely of an extension name'
	"""
	# if it ends with a forward slash, get rid of the last character
	original = genome_index_path
	if genome_index_path[len(genome_index_path) - 1] == '/': 
		genome_index_path = genome_index_path[:len(genome_index_path) - 1]
	# find everything past the last forward slash
	index = genome_index_path.rfind('/')
	index += 1
	genome_index = ''
	while index < len(genome_index_path):
		genome_index += genome_index_path[index]
		index += 1
	if not ('.fa' in genome_index or '.fasta' in genome_index):
		raise InputError('genome_index', 'must contain .fa or .fasta')
	# now we have something that looks like 'genome.fasta' and want to convert it to something that looks like 'genome'
	gilist = genome_index.split('.')
	genome_index_base = gilist[0]
	if genome_index_base == '':
		raise InputError('genome_index_path', 'cannot be a file that consists solely of an extension name')
	assert original == genome_index_path # assert that genome_index_path has not been changed during the course of the function
	return genome_index_base

def process_genome_index_path(genome_index_path):
	genome_index_path = check_pathname('genome_index_path', genome_index_path)
	if not ('.fa' in genome_index_path or '.fasta' in genome_index_path):
		raise InputError('genome_index_path', 'must contain .fa or .fasta')
	print('Genome index path is ' + genome_index_path)
	return genome_index_path

def create_output_dir(work_dir):
		output_dir = work_dir + 'at' + curr_time + 'on' + curr_day
		#add a random id to the folder so we don't overwrite a folder if we run this script twice within the same minute. 
		#if by an unlikely chance the same random number is generated, neverfear, since the script will abort if a folder of the same name already exists.
		output_dir += '_' + str(random.randint(1,99)) + '/' 
		# use subprocess.check_call instead of subprocess.Popen() b/c Popen creates a separate process which might cause a race condition where the folder doesn't exist when it is trying to be written to.
		subprocess.check_call(['mkdir', output_dir], stdout=subprocess.PIPE)
		return output_dir

def confirmIsComment(line):
	if not '#' in line:
		raise InputError('comment lines','must contain #')


def process_work_dir(work_dir):	
	work_dir = check_pathname('work_dir', work_dir)
	#work_dir should end with '/' for later parts of script to function properly
	if work_dir[len(work_dir) - 1] != '/':
		work_dir += '/'
	print('Working directory is ' + work_dir)
	return work_dir

def process_email(email):
	if not '@' in email:
		raise InputError('email','must contain \'@\'')
	return email

def find_matching_reads(fastqs):
	"""
	Iterate through a list of filenames looking for .fastq files that are identical
	except one contains 'R1' and the other contains 'R2'
	return a dictionary where the keys are R1s and the values are the matching R2s
	>>> find_matching_reads(['/shares/RF_293_R2_m.fastq','/shares/RF_293_R1_m.fastq','/shares/RF_754_R1_m.fastq','/shares/RF_754_R2_m.fastq'])
	Matching .fastq R1 and R2 files that have otherwise identical filenames
	{'/shares/RF_754_R1_m.fastq': '/shares/RF_754_R2_m.fastq', '/shares/RF_293_R1_m.fastq': '/shares/RF_293_R2_m.fastq'}
	>>> find_matching_reads(['/shares/RF_293_R2_m.fastq','/shares/RF_293_R1_m.fastq','/shares/RF_75_R1_m.fastq','/shares/RF_754_R2_m.fastq'])
	Traceback (most recent call last):
		...
	InputError: '/shares/RF_75_R1_m.fastq had no matching R2 file'
	"""
	print('Matching .fastq R1 and R2 files that have otherwise identical filenames')
	reads = {} # dictionary to hold name of read files where key is read1 and value is read2.
	for fastq in fastqs:
		if fastq.find('R1') != -1: # if we found R1 in the file name then
			prefix = fastq[:fastq.find('R1')] 
			postfix = fastq[fastq.find('R1')+2:]
			alreadyFoundAMatch = False
			for potential_match in fastqs:
				# if we found R2 in a filename and it has the same prefix and postfix as filename containing R1, then
				if potential_match.find('R2') != -1 and potential_match[:potential_match.find('R2')] == prefix and potential_match[potential_match.find('R2')+2:] == postfix: 
					if alreadyFoundAMatch:
						#I don't think it would ever be possible to have more than one matching R2 file
						#because the filenames come from scanning a directory
						#and the operating system wouldn't allow a directory to contain two files of the same name
						#but oh well might as well leave it here
						raise InputError(fastq, 'had more than one matching R2 file')
					match = potential_match
					alreadyFoundAMatch = True
			if not alreadyFoundAMatch:
				raise InputError(fastq, 'had no matching R2 file')
			else:
				reads[fastq] = match
	return reads

if __name__ == '__main__':
	#import doctest
	#doctest.testmod()
	if len(sys.argv) != 2:
		sys.exit("You have to enter the name of the input file, with info like the location of the FASTQ files, working directory in Circe, etc. after the script name. So for instance if input.txt is the text file you'd enter `python slurm_generator.py input.txt` on the command line.")
	input_file = sys.argv[1]
	print('Reading from file ' + input_file)
	
	dev_flag = development_available()

	with open(input_file, 'r') as f:

		bt_script = BtScript()
		seq_script = SeqScript()
		note = Note()
		job_notok = JobNotOk()
		job_ok = JobOk()
		job_any = JobAny()


		# calling .readline() without assigning the return value gets rid of a comment line
		confirmIsComment(f.readline()) 
		confirmIsComment(f.readline())
			
		output_dir = create_output_dir(process_work_dir(f.readline()))
		
		confirmIsComment(f.readline())
		email_to_send_output_to = process_email(f.readline())
		
		confirmIsComment(f.readline())
		seq_script.sample_dir = f.readline()
		
		confirmIsComment(f.readline())
		seq_script.gff = f.readline()
		
		confirmIsComment(f.readline())
		genome_index_path = process_genome_index_path(f.readline())
		genome_index_base = extract_genome_index_base(genome_index_path)

		confirmIsComment(f.readline())
		seq_script.slurm.threads = f.readline()

		confirmIsComment(f.readline())
		seq_script.slurm.mem = f.readline()
		
		confirmIsComment(f.readline())
		seq_script.slurm.hours = f.readline()
		
		confirmIsComment(f.readline())
		seq_script.slurm.minutes = f.readline()

		bt_script.output()
		seq_script.output(bt_script.slurm.slurm_id)

		#job_notok.output(seq_script)
		#job_ok.output(seq_script)
		#job_any.output(seq_script)
	with open(input_file, 'r') as f:
		input_contents = f.read()	
	note.output(bt_script.slurm.slurm_id, seq_script, input_contents)
