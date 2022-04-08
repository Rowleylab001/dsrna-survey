# Note: Don't download more than 100 sequences during peak usage times
# Make sure you have installed biopython: https://biopython.org/wiki/Download
# Make sure you are using Python 3

# import modules
from Bio import Entrez
from Bio import SeqIO
import os

####### Gather user input #######
# propt user for email (tells NCBI who user is)
Entrez.email = input("Enter your email: \n")
# list of genbank IDs to be downloaded
seq_list = ["CP004103","CP004788","CP004789","CP004790","CP004791","CP004792","CP004793","CP004794","CP004795","CP004796","CP004797","CP004798","CP004800","CP004801","CP004802","CP004803","CP004804","CP004805","CP004806","CP004807","CP004808","CP004809","CP004810","CP004811","CP004812","CP004813","CP004814","CP004815","CP004816","CP004817","CP004818","CP004819","CP004820","CP004821","CP004822","CP004823","CP004824","CP004825","CP004826","CP004827","CP004828","CP004829","CP004830","CP004831","CP004832","CP004833","CP004834","CP004835","CP004836","CP004837","CP004838","CP004839","CP004840","CP004841","CP004842","CP004843","CP004844","CP004845","CP004799","CP004754","CP004755","CP004756","CP004757","CP004758","CP004759","CP004760","CP004761","CP004762","CP004763","CP004764","CP004765","CP004766","CP004767","CP004768","CP004769","CP004770","CP004771","CP004772","CP004773","CP004774","CP004775","CP004776","CP004777","CP004778","CP004779","CP004780","CP004781","CP004782","CP004783","CP004784","CP004785","CP004786","CP004787","CP004107","CP005016","CP005017","CP005018","CP005019","CP005020","CP005021","CP005022","CP005023","CP005024","CP005025","CP005026","CP005028","CP005029","CP005030","CP005031","CP005032","CP005033","CP005034","CP005035","CP005036","CP005037","CP005038","CP005039","CP005040","CP005041","CP005042","CP005043","CP005044","CP005045","CP005046","CP005047","CP005048","CP005049","CP005050","CP005051","CP005052","CP005053","CP005054","CP005055","CP005056","CP005057","CP005058","CP005059","CP005060","CP005061","CP005062","CP005063","CP005064","CP005065","CP005066","CP005067","CP005068","CP005069","CP005070","CP005071","CP005072","CP005073","CP005027","CP004982","CP004983","CP004984","CP004985","CP004986","CP004987","CP004988","CP004989","CP004990","CP004991","CP004992","CP004993","CP004994","CP004995","CP004996","CP004997","CP004998","CP004999","CP005000","CP005001","CP005002","CP005003","CP005004","CP005005","CP005006","CP005007","CP005008","CP005009","CP005010","CP005011","CP005012","CP005013","CP005014","CP005015"]

####### MAKE DIRECTORY ############
# define the name of the directory to be created
current_path = os.getcwd()
new_directory = "FGSC_chr_files/"
new_directory_full_path = str(current_path) + "/" + str(new_directory)
# make new sub-directory
try:
    access_rights = 0o755 # define the access rights
    os.mkdir(new_directory, access_rights)
except OSError:
    print ("Creation of the directory %s failed" % new_directory)
else:
    print ("Successfully created the directory %s " % new_directory)

############ DOWNLOAD FILES #############
# reference for "db": https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# reference for "rettype": https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
i = 0
for seq_acc in seq_list:
    handle = Entrez.efetch(db="nucleotide", id=seq_acc, rettype="fasta", retmode="text")
    myfilename = seq_acc + '.fasta'
    myfile = open(os.path.join(new_directory, myfilename), 'w')
    myfile.write(handle.read())
    myfile.close()
    i += 1

# print the number of files downloaded and the directory they're downloaded to
print ("Successfully downloaded", i, "sequence files in", new_directory_full_path)
