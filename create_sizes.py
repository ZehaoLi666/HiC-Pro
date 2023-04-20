from sys import argv

if len(argv) != 3:
	print("Run with 2 arguments: input_filename.txt output_filename.wig")
	exit(0)
else:
	pass

script, input_filename, output_filename 

current_file = open(input_filename)