import subprocess
import artiq

"""ARTIQ installation diagnostics.
"""



def main():
	
	print("path to python3 should be ~/anaconda/bin/python")
	bashCommand = "which python3"
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	output = process.communicate()[0]
	print("your path to python3: {}".format(output))
	print()

	print("path to artiq")
	print("your artiq path: {}".format(artiq.__path__[0]))
if __name__ == "__main__":
    main()