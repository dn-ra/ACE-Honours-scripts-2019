import subprocess
import os

completion_dictionary = {}


read_files_dir = '/srv/home/s4204666/abisko_readfiles'
read_files_2 = '/srv/home/s4204666/abisko/data/flat20180301_2013-2017'
assembly_files_dir = '/srv/home/s4204666/Abisko_Assemblies/filtered_2000bp'
out_dir = '/srv/home/s4204666/Abisko_Assemblies/filtered_2000bp/bam_files'

assemblies = [f for f in os.listdir(assembly_files_dir) if f.startswith('filtered.2000bp')]
read_files = os.listdir(read_files_dir)
read_files2 = os.listdir(read_files_2)
for file in assemblies:
	read_source = read_files_dir
	orig_assembly = file.split(".")[2]
	read_couple = [f for f in read_files if orig_assembly in f]
	if len(read_couple) == 0:
		read_couple = [f for f in read_files2 if orig_assembly in f]
		read_source = read_files_2
	cmd = 'bamm make -d {} -c {} {} --quiet -o {} -t 24'.format("/".join([assembly_files_dir, file]), "/".join([read_source, read_couple[0]]), "/".join([read_source, read_couple[1]]), out_dir)
	print(cmd)
	completion_dictionary[file] = subprocess.call(['bash', '-c', cmd])
print("all assemblies processed")
print(completion_dictionary)
