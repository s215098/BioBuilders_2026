#!/usr/bin/python

#############################################
###    Technical University of Denmark    ###
#############################################
### When using the CUPP programs or results, please cite:
### Barrett K, Lange L. Peptide-Based Functional Annotation of Carbohydrate-Active Enzymes by Conserved Unique Peptide Patterns (CUPP). 
### Biotechnol Biofuels. BioMed Central; 2019;1–21. 
### Available from: https://doi.org/10.1186/s13068-019-1436-5
#############################################

import os,sys,json,argparse,time,re,psutil,random,math,gzip,itertools,platform
from scipy.cluster.hierarchy import dendrogram,linkage,cut_tree, to_tree
from multiprocessing import Pool, Queue, Manager, Process, Pipe
from scipy.spatial.distance import squareform 
from scipy import ndimage
import numpy as np
random.seed(1)


def cupp_logo():
    '''
    Print out the CUPP logo
    '''
    print ("###############################################################\n"+
           "###############################################################\n"+
           "###                                                         ###\n"+
           "###        #####     ###    ###   ########    ########      ###\n"+
           "###      ###   ###   ###    ###   ###   ###   ###   ###     ###\n"+
           "###     ###          ###    ###   ###   ###   ###   ###     ###\n"+
           "###     ###          ###    ###   ########    ########      ###\n"+
           "###     ###          ###    ###   ###         ###           ###\n"+
           "###      ###   ###    ###  ####   ###         ###           ###\n"+
           "###        #####       ##### ##   ###         ###           ###\n"+
           "###                                                         ###\n"+
           "###############################################################\n"+
           "###       Powered by Technical University of Denmark        ###\n"+
           "###############################################################")

def arg_parser_predict(args):
    '''
    Process command-line user arguments and prepare for further directories
    '''

    ####################################
    #Folders and files:
    parser = argparse.ArgumentParser(description='The input parameters for the CUPP identification technology')
    parser.add_argument('-working_dir',     default=os.getcwd()+"/CUPP",type=str,        help="Change working directory for both clustering and prediction")
    parser.add_argument('-fasta_extension', default=".faa",    type=str,                 help="Define the fasta extension to not be included in the family abbreviation")
    
    ####################################
    #CUPPclustering parameters:
    parser.add_argument('-ambiguous',       default=2,         type=int,                 help="Introduce a number of X's into each peptides as ambiguous amino acids")
    parser.add_argument('-n_mer',           default=8,         type=int,                 help="Peptide length for both clustering and prediction")
    parser.add_argument('-cdhit',           default=0.90,      type=float,               help="Use only representative proteins but meta data from all entries in CD-HIT cluster")
    parser.add_argument('-domain_off',     default=False,     action="store_true",      help='Use full protein fasta only instead of both')
    parser.add_argument('-cc',              default=9,         type=float,               help="Amplification of clustering c_clust")
    
    ####################################
    #CUPPprediction parameters:
#    parser.add_argument('-compiled_json',   default="/work3/kbaka/8x2_90_v1.1.0_transferred_CUPPlibrary.json",        type=str,                 help="Specify full path to CUPPlibrary")
    parser.add_argument('-compiled_json',   default="",        type=str,                 help="Specify full path to CUPPlibrary")
    parser.add_argument('-query',           default="",        type=str,                 help="Select a fasta file for CUPP clustering or prediction")
    parser.add_argument('-dir_query',       default="",        type=str,                 help="Process all files in folder")
    parser.add_argument('-best_hit',        default=0.5,       type=float,               help="When to consider a secondary hit equally relevant as the first")
    parser.add_argument('-occurence',       default=False,     action="store_true",      help="Retrun only the EC number which the highest occurence of the CUPPgroup")
    parser.add_argument('-overlap',         default=0.6,       type=float,               help="How much can a domain overlap with an already assigned domain")
    parser.add_argument('-minimum_cup',     default=5,         type=float,               help="Required frequency for positive hit")
    parser.add_argument('-evidence',        default=0.01,      type=float,               help="Minimum percentage of frequency of CUPPgroup for hit")
    parser.add_argument('-position_cut',    default=0.2,       type=float,               help="Frequency for a position to be included in the domain")
    parser.add_argument('-domain_min',      default=20,        type=int,                 help="Minimum number of positions covered be a domain")
    parser.add_argument('-output_path',     default="",        type=str,                 help="Specify a full path and filename of the output fasta file (including fasta extension)")
    parser.add_argument('-genomic',         default=0,         type=int,                 help="Treat the query as DNA and locate ORFS and predict them individually, specify minimum number of bases")
    parser.add_argument('-max_longst',      default=50000,     type=int,                 help="Maximum length of protein")
    parser.add_argument('-cupp_minimum_score',default=5,       type=int,                 help="The minimum score for a assignment of CUPP group and its EC function to a query protein")
    
    ####################################
    #Parallel computing and documentation:
    parser.add_argument('-table',           default=False,     action="store_true",      help="Show the CUPP library prior to prediction!")
    parser.add_argument('-version',         default="v1.0.14", type=str,                 help="The folder of fastafiles from which clusters are made")
    parser.add_argument('-beta_options',    default=[],        nargs='+',                help="New features under development",choices=["validate","debug","mix_type","no_fam_ec","float","experimental_families","keep","catalytic","relatedness_alone","prediction_fails","complex"])
    parser.add_argument('-exclude_family',  default=[],        nargs='+',                help="Select families which will be ignored from the prediction")
    parser.add_argument('-type',            default="both",    type=str,                 help='Select output formats',choices=["none","quick","fasta","raw","both"])
    parser.add_argument('-sep',             default="|",       type=str,                 help='Select the separator in the header for both input and output')
    parser.add_argument('-live',            default=False,     action="store_true",      help='Keep the CUPP library loaded and predicted in serial FASTA files as query')
    parser.add_argument('-sleep',           default=0,         type=int,                 help="Number of seconds to sleep")
    parser.add_argument('-keep_only',       default=[],        nargs='+',                help="Select few families which are the only ones predicted")
    args = parser.parse_args(args)

    
    ####################################
    #Prepare output paths:
    args.exclude_family = set(args.exclude_family)
    args.working_dir = os.path.normpath(args.working_dir)
    os.makedirs(args.working_dir+"/predicted", exist_ok=True)
    os.makedirs("%s/output" % args.working_dir, exist_ok=True)    
    args.dir_query = dict((os.path.join(args.dir_query,fil),1) for fil in os.listdir(args.dir_query)) if os.path.exists(args.dir_query) else {}
    args.raw_path = "%s/predicted/%s_CUPP.log" % (args.working_dir,os.path.basename(args.query).split(".")[0]) if not args.output_path else args.output_path+".log"
    args.fasta_path = "%s/predicted/%s_CUPP%s" % (args.working_dir,os.path.basename(args.query).split(".")[0],args.fasta_extension) if not args.output_path else args.output_path
    if "." not in os.path.basename(args.fasta_path):
        args.fasta_path += args.fasta_extension
    
    ####################################
    #Avoid overwriting the CUPlibrary accidentally:
    args.setting = "%s%s%sx%s_%s" % ("f" if args.domain_off else "d","c" if args.cdhit else "a",args.n_mer,args.ambiguous,int(args.cc*10))
    args.compiled_json = "%s/CUPPlibrary/%s_%s_CUPPlibrary.json" % (args.working_dir,args.setting[2:],args.version) if not args.compiled_json else args.compiled_json
    args.beta_options = set(args.beta_options)
    time.sleep(args.sleep)
    
    return (args)
    
class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

def merge_fasta_dir(liste=[],output="CUPP/merged.faa"): 
    '''
    Merge several fasta files into one file with merged meta data and non-redundant sequences
    '''
    collection = {}; header = ""
    for fil in liste:
        handle = open(fil)
        for line in handle: 
            line = line.strip()
            
            if line[0] == ">":
                header = line
            elif header: 
                if line not in collection:
                    collection[line] = []
                if header.split("|")[1][:3] != "CBM" and header.split("|")[1][2] != "0":
                    collection[line].append(header)
                header = ""
        handle.close()

    cat = ["Accession","Family","Subfam","raw_function","Uniprot","PDB","Name","Taxid"]
    handle = open(output,"w")
    for seq,v in collection.items():
        meta = [set() for n in range(len(cat))]
        for l in v:
            parts = l[1:].split("|")
            for n in range(len(cat)):
                if n in [2,3] and parts[n]:
                    parts[n] = parts[1] + ":" + parts[n] if "CBM" not in parts[1] else ""
                if parts[n]:
                    meta[n].add(parts[n])
                    
        header = "|".join(["+".join(sorted(e)) for e in meta])
        
        #if 1 < len(meta[1]) or 1 < len(meta[2]):
        #    for acc in meta[0]:
        #        print ("%s\t#41CB3A\t%s" % (acc,"|".join(["+".join(sorted(e)) for e in meta[1:]])))
            
        if meta[0]:
            handle.write(">%s\n%s\n" % (header,seq))
    handle.close()

def export_CUPP_members(CUPPpool_path,original_list,cdhit_list,top_members=5):
    '''
    BETA FUNCTION
    '''

    with open(CUPPpool_path) as json_data:
        lib = json.load(json_data)
        
    redundant_proteins = load_redundancy(cdhit_list,"NEW1") if [f for f in cdhit_list if f] else {}
    collection, fasta_meta = obtain_collection(original_list)

    members = {}; group = {}; limit = {}
    for cupp,v in lib[list(lib.keys())[0]]["meta"].items():
        for acc,score in sorted(v[0].items(),key=lambda x: x[1],reverse=True):
            if not top_members or score:
                real = acc.split(":")[0]
                if real in members:
                    print ("ERROR2 %s" % acc)
                else:
                    if cupp not in limit:
                        limit[cupp] = {}
                    if not top_members or len(limit[cupp]) < top_members + 1:
                        limit[cupp][real] = score
                        members[real] = cupp

    for acc in members:
        if acc in redundant_proteins:
            for a in redundant_proteins[acc]:
                if a[0] not in group:
                    group[a[0]] = set()
                group[a[0]].add(acc)
    
    temp = {}
    for acc,g in members.items():
        temp[acc] = g
        
        if acc in redundant_proteins and not top_members:
            gr = redundant_proteins[acc][0]
            if gr in group:
                for a in group[gr]:
                    temp[acc] = g
    members = temp
    for seq,acc in sorted(collection.items(),key=lambda x: family_sort(members[x[1]]) if x[1] in members else [family_sort(members[a]) for a in fasta_meta[x[1]]["Accession"] if a in members][0] if [family_sort(members[a]) for a in fasta_meta[x[1]]["Accession"] if a in members] else tuple()):
        for a in fasta_meta[acc]["Accession"]:
            if not top_members or a in members:
                print (">%s|%s|%s|%s|%s|%s" % (a,members[a] if a in members else "","|".join(sorted(fasta_meta[acc]["Subfam"])),"|".join(sorted(fasta_meta[acc]["PDB"])),"|".join(sorted(fasta_meta[acc]["raw_function"])),limit[members[a]][a] if a in members and members[a] in limit else ""))#,"|".join(sorted(fasta_meta[acc]["Name"]))))
                print (seq)

def use_guide():
    '''
    BETA FUNCTION
    '''
    ####################################
    #Restrict the clustering to only include peptides from positions shared by a guide domain:
    if args.guide:
        
        ####################################
        #Generate a pool of domain relevant peptides:
        print ("### Restriction of peptides to domain is ON")
        domain_peptides = set()
        guide_collection, guide_fasta_meta = obtain_collection(args.guide, sep = args.sep, common = args.general_name, meta_cat = args.header, beta = args.beta_options)
        for seq,acc in guide_collection.items():
            dom_pep = obtain_peptides({seq:acc},n_mer=args.n_mer,ambiguous=args.ambiguous)
            domain_peptides.update(set(dom_pep[acc]))
        print ("### Number of peptides relevant in guide domains: %s" % len(domain_peptides))
            
        for seq, acc in collection.items():
            original_seq_pep_short = obtain_peptides({seq:acc},use_positions=False,n_mer=args.n_mer,ambiguous=args.ambiguous,beta=args.beta_options,limit=domain_peptides)
            if 20 < len(original_seq_pep_short[acc]):
                collection = dict((seq,v) for seq,v in collection.items() if original_seq_pep_short[collection[seq]])
                original_seq_pep.update(obtain_peptides({seq:acc},use_positions=False,n_mer=args.n_mer,ambiguous=args.ambiguous,beta=args.beta_options))
                original_seq_pep2 = dict((acc,dict((pep,s) for pep,s in v.items() if pep in domain_peptides)) for acc,v in original_seq_pep.items())
        print ("--> Number of proteins left in collection: %s [%s]" % (len(collection),time.time() - start))
        
        for acc, v in original_seq_pep2.items():
            target_pool = set()
            for pep,n in v.items():
                for a in range(args.n_mer):
                    target_pool.add(a+n)
            for pep,n in original_seq_pep[acc].items():
                if n in target_pool:
                    original_seq_pep2[acc][pep] = n
       
        for acc,v in sorted(original_seq_pep.items(),key=lambda x: len(original_seq_pep2[x[0]])/len(original_seq_pep[x[0]])):
            if original_seq_pep[acc]:
                print (acc,len(original_seq_pep[acc]),len(original_seq_pep2[acc]),fasta_meta[acc]["Subfam"])
        original_seq_pep = original_seq_pep2
        
    ####################################
    #Output the peptide-related proteins:
    if collection and (args.guide or "export" in args.beta_options):
        handle = open("output/%s_CUPP_BLAST.faa" % args.general_name,"w")
        for seq,acc in collection.items():
            handle.write(">%s|%s\n%s\n" % (acc,"|".join(["+".join([str(i) for i in sorted(v)]) for f,v in sorted(fasta_meta[acc].items(),key=lambda x:args.header.index(x[0])) if f in args.header]),seq))
        handle.close()
    

def submit_hpc(run="", ID="NEW1", mail="", hours=24, ram=20, jobs = 1):
    '''
    Create a batch file for submission to the hpc server
    '''
    #Run with more than 128GB RAM:
    jobs = int(ram/128)+1 if int(ram/128)+1 < jobs else jobs
    
    ####################################
    #Create the batch file for submission:
    os.makedirs("status", exist_ok=True)    
    batch_file = "status/submit_%s.sh" % ID
    handle = open(batch_file,"w")
    handle.write("#!/bin/sh\n")
    ### -- specify queue -- 
    handle.write("#BSUB -q %s\n" % ("hpc" if ram != 300 else "epyc"))
    ### -- set the job Name -- 
    handle.write("#BSUB -J %s\n" % ID)
    ### -- ask for number of cores (default: 1) -- 
    handle.write("#BSUB -n %s\n" % jobs)
    ### -- specify that the cores must be on the same host -- 
    handle.write('#BSUB -R "span[hosts=1]"\n')
    ### -- specify that we need 2GB of memory per core/slot -- 
    handle.write('#BSUB -R "rusage[mem=%sGB]"\n' % (int(ram/jobs)+1)) 
    ### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
    #handle.write("#BSUB -M %sGB " % ram*2 if ram*2 < 128 else 128)
    ### -- set walltime limit: hh:mm -- 
    handle.write("#BSUB -W %s:00\n" % hours)
    ### -- set the email address -- 
    handle.write("#BSUB -u kbaka@dtu.dk\n")
    ### -- send notification at start -- 
    #BSUB -B 
    ### -- send notification at completion -- 
    #BSUB -N 
    ### -- Specify the output and error file. %J is the job-id -- 
    ### -- -o and -e mean append, -oo and -eo mean overwrite -- 
    handle.write("#BSUB -oo status/"+ID+"_Output_%J.out\n") 
    handle.write("#BSUB -eo status/"+ID+"_Error_%J.err\n")
    handle.write(run)
    handle.close()
    time.sleep(0.1)

    ####################################
    #Submit the batch script to the hpc (running on the hpc):
    os.system("bsub < %s" % batch_file)
    return run,"["+str(ram)+";"+str(hours)+"]"   


    
def submit_hpc_old(run="", ID="NEW1", mail="", hours=48, ram=0, jobs = 1):
    '''
    Create a batch file for submission to the hpc server
    '''
    
    ####################################
    #Create the batch file for submission:
    os.makedirs("status", exist_ok=True)    
    batch_file = "status/submit_%s.txt" % ID
    handle = open(batch_file,"w")
    handle.write("#!/bin/bash\n")
    handle.write("#PBS -o status/%s_$PBS_JOBID.out\n" % ID)
    handle.write("#PBS -e status/%s_$PBS_JOBID.err\n" % ID)
    handle.write("#PBS -l mem=%sgb\n" % ram) if ram else ""
    handle.write("#PBS -l nodes=1:ppn=%s\n" % jobs)
    handle.write("#PBS -l walltime=%s:00:00\n" % hours)
    handle.write("#PBS -M %s\n" % mail) if mail else ""
    handle.write("#PBS -m abe\n\n") if mail else ""
    handle.write("if test X$PBS_ENVIRONMENT = XPBS_BATCH; then cd $PBS_O_WORKDIR; fi\n\n")
    handle.write(run)
    handle.close()
    
    ####################################
    #Submit the batch script to the hpc (running on the hpc):
    time.sleep(0.1)
    os.system("qsub %s" % batch_file)
    return run,"["+str(ram)+";"+str(hours)+"]"
    
    
def show_arguments(args,cluster=False,recom=False,predict=False):
    '''
    Display the currently applied settings
    '''
    if cluster:
        print ("#"*60+"\n###   CUPPclustering by Technical University of Denmark  ###\n"+"#"*60)
    elif predict:
        print ("#"*60+"\n###   CUPPprediction by Technical University of Denmark  ###\n"+"#"*60)
        
    arguments = dict((k,v) for k,v in vars(args).items() if k != "verbose")
    general_included = {"common","final_cut","ambiguous","n_mer","jobs",
                        "rel_cut","setting","version","working_dir","seed",
                        "database_version","ignored","domain_off","beta_options"}
    cluster_included = {"start_positions","deposit","dbcan_folder","cut","outlier"
                        "cup_cut","cc","cdhit","cdhit_folder","user_ncbi",
                        "fasta_extension","gradient","minimum_group_size",
                        "header","expr","domain_query","redundancy_path"}
    predict_included = {"evidence","domain_min","best_hit","output",
                        "minimum_cup","overlap","query","option","profile","genomic"}
    arguments["common"] = " ".join(arguments["common"])
    arguments["gradient"] = " ".join([str(f) for f in arguments["gradient"]])
    
    included_arguments = dict((k,v) for k,v in arguments.items() if k in general_included)
    included_arguments.update(dict((k,v) for k,v in arguments.items() if arguments["clustering"] and k in cluster_included))
    included_arguments.update(dict((k,v) for k,v in arguments.items() if arguments["predict"] and k in predict_included))
    for arg,value in sorted(included_arguments.items()):
        print ("### %s:%s%s" % (arg," "*(20-len(arg)),value if arg != "ignored" else len(value))) if value else None 
        
def family_sort(x,start=2):
    '''
    Sort a list of string according to natural sorting by inital letters and/or the numbers
    '''
    match = re.findall(r'\d+',x)
    return(tuple([x[:start]]+[int(f) for f in match]))

def common_abb(expr,path,extension=".faa"):
    '''
    Use regular expression to identify relevant file names in folder
    '''
    expr = re.compile(expr); targets = {}
    for fil in os.listdir(path):
        if re.search(expr, fil):
            targets[fil.replace("%s.dbcan" % (extension),"").replace(extension,"")] = 1
    return [fil for fil in sorted(targets,key=lambda x: family_sort(x))]
    
def ram():
    '''
    Determine the current RAM usage by the process
    '''
    py = psutil.Process(os.getpid())
    return round(py.memory_info()[0]/2.**30,2)
    
def taxonomy_database(nodes_file="resources/nodes.dmp",names_file="resources/names.dmp",old_file="resources/merged.dmp",delete_file="resources/delnodes.dmp"):
    '''
    Loading the NCBI taxonomy database names and nodes
    '''
    
    ####################################
    #Download the NCBI nodes and names if missing:
    if not os.path.exists(nodes_file) or not os.path.exists(names_file):
        print("### Download the NCBI nodes and names, when extract the content into the folder resources!")
        sys.exit("--> ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    print ("### Loading ncbi taxonomy...")
    ncbi_names, ncbi_nodes = {"1":["root"]},{"1":["1","root"]}
    
    ####################################
    #Handling old taxids:
    for line in open(old_file):
        line = line.strip().split("\t")
        ncbi_nodes[line[0]] = [line[2],"unknown*"]
        ncbi_names[line[0]] = ["unknown*"]
    
    ####################################
    #Handling deleted taxids:
    for line in open(delete_file):
        line = line.strip().split("\t")
        ncbi_nodes[line[0]] = ["1","unknown*"]
        ncbi_names[line[0]] = ["unknown*"]
        
    ####################################
    #Reading the NCBI nodes into dict:
    for line in open(nodes_file):
        line = line.strip().replace("\t","").split("|")
        ncbi_nodes[line[0]] = [line[1],line[2]]
        
    ####################################
    #Reading the NCBI names into dict:
    for line in open(names_file):
        taxID, name, comment, nameClass = line.strip().replace("\t","").split("|")[:-1]
        taxID = taxID
        if taxID in ncbi_names:
            if nameClass == "scientific name":
                ncbi_names[taxID] = [name]
        else: ncbi_names[taxID] = [name]
        
    ####################################
    #Create a framework for the nodes and names as one dict:
    ncbi = {"nodes":ncbi_nodes,"names":ncbi_names,"partial":0,"banned":{}}

    return (ncbi)

def kingdom(tax,ncbi,temp={},verbose=False,text=True,preloaded=False):
    '''
    Returns a ordered list of taxids and their scientific names
    '''
    ####################################
    #Definig the positions of different taxonomical levels:
    kingdom_ranks = ["species", "genus", "family", "order", "subclass", "class", "subphylum", "phylum", "subkingdom", "kingdom", "superkingdom"]
    ranks = [["" for f in range(len(kingdom_ranks))],["" for f in range(len(kingdom_ranks))]]
    
    ####################################
    #Preparing input tax:
    if not str(tax).isdigit():
        sys.stderr.write("?") if verbose else ""
        return (ranks[text],ncbi,temp)
    tax_list = [str(tax)]
    
    ####################################
    #Moving backwards in the lineage until root ("1") is reached or invalid taxid is found: 
    while tax_list[-1] != "1":
    
        ####################################
        #If the big one is loaded, the missing taxids are added to the partial dict:
        if preloaded:
            if tax_list[-1] in ncbi["nodes"]:
                temp["nodes"][tax_list[-1]] = ncbi["nodes"][tax_list[-1]]
                temp["names"][tax_list[-1]] = ncbi["names"][tax_list[-1]]
                
        ####################################
        #Create the tax-line of the id:
        if tax_list[-1] and tax_list[-1] in ncbi["nodes"]:
            tax_list.append(ncbi["nodes"][tax_list[-1]][0])
        else:
            sys.stderr.write("!") if verbose else ""
            return (ranks[text],ncbi,temp)
            
    ####################################
    #Create a list of the correct taxonomical level in specific positions:
    for taxID in tax_list:
        for no in range(len(kingdom_ranks)):
            if ncbi["nodes"][taxID][1] == kingdom_ranks[no]: 
                ranks[0][no] = taxID
                ranks[1][no] = ncbi["names"][taxID][0]
                
    return (ranks[text],ncbi,temp)
    
def load_redundancy(files,common=[],sep="|"):
    '''
    Load protein redundency clustering from CD-HIT.
    '''
        
    ####################################
    #Load the redundancy for the current collection:
    redundancy = {}; no = 0
    
    for path in files:
        acc = ""; cluster = False
        handle = open(path)
        for line in handle:
            line = line.strip()
            if line:
                if acc and cluster:
                    if acc in redundancy:
                        if star:
                            if cluster != redundancy[acc][0]:
                                sys.stderr.write("!!! WARNING: The redundancy file have the same accession number twice [%s]\n" % acc)
                            redundancy[acc] = [cluster,star]
                    else:
                        redundancy[acc] = [cluster,star]
                    acc = ""
                if line[0] == ">":
                    cluster = common[no]+"_"+line.split()[-1]
                elif 3 <= len(line.split()):
                    acc = line.split()[2].strip(">").split(sep)[0]
                    #acc = line.split()[2].strip(">").split("|")[0]
                    star = True if line.split()[3] == "*" else False
                else:
                    acc = ""; cluster = False
        if acc and cluster:
            if acc in redundancy:
                if star:
                    sys.stderr.write("!!! WARNING: The redundancy file have the same accession number twice [%s]\n" % acc)
                    redundancy[acc] = [cluster,star]
            else:
                redundancy[acc] = [cluster,star]
        handle.close(); no += 1
    
    return (redundancy)

def obtain_collection(files,cdhit={},sep="|",ignore={},common="NEW1",meta_cat = ["Accession","Family","Subfam","raw_function","Uniprot","PDB","Name","Taxid"],beta = {}):
    '''
    Take a query folder containing fastafiles for families/groups of sequences
    '''

    ####################################
    #Create set of accession numbers to ignore:
    ignore = set((cdhit[k][0]) for k in ignore) if cdhit else ignore
    ignore = set((k) for k,v in cdhit.items() if v[0] in ignore) if cdhit else ignore
    use_common_family = False if "Family" in meta_cat else True
    meta_cat.append("Family") if use_common_family else None
    #Future Improvement: Take full header as accession number if Accession not in args.header
    
    ####################################
    #Read a fasta files and obtain a dict of sequences (keys) and header (values):
    sequence = ""; cd_seq = {}; collection = {}; meta_data_acc = {}; cd_merge = {}
    
    for path in files:
    
        ####################################
        #Read the lines of a collection:
        temp_collection = {}; sequence = ""; header = ""; count = 0; seq_non = -1; flag = False
        import codecs
        with codecs.open(path, 'r', encoding='utf-8', errors='ignore') as handle_in:        

            for line in handle_in:
            
                ####################################
                #Ignore emtry lines:
                line = line.strip()
                if not line:
                    continue
                    
                ####################################
                #When an entry has a sequence and a header, add it to collection:
                if line[0] == ">" and sequence:
                    
                    #header = header.split(" ")[0] + (header.split(sep,1)[-1] if sep in header else "")
                    acc = header.split(sep)[0].strip(">")
                    if sequence not in temp_collection:
                        temp_collection[sequence] = {}
                    header = header.strip() + sep*(len(meta_cat)-header.count(sep))
                    temp_collection[sequence][header] = cdhit[acc] if acc in cdhit else ["ALONE"+str(len(temp_collection))+str(len(temp_collection[sequence])),False]
                    if cdhit and acc not in cdhit:
                        count += 1
                        cdhit[acc] = ["ALONE"+str(len(temp_collection))+str(len(temp_collection[sequence])),False]
                    sequence = ""; header = ""
                if line[0] == ">":
                    header = line
                else:
                    sequence += re.sub('[^CSTPAGNDEQHRKMILVFYW]', '_', line.upper().strip())

            if header and sequence:
            
                #header = header.split(" ")[0] + (header.split(sep,1)[-1] if sep in header else "")
                acc = header.split(sep)[0].strip(">")
                if sequence not in temp_collection:
                    temp_collection[sequence] = {}
                header = header.strip() + sep*(len(meta_cat)-header.count(sep))
                temp_collection[sequence][header] = cdhit[acc] if acc in cdhit else ["ALONE"+str(len(temp_collection))+str(len(temp_collection[sequence])),False]
                if cdhit and acc not in cdhit:
                    count += 1
                    cdhit[acc] = ["ALONE"+str(len(temp_collection))+str(len(temp_collection[sequence])),False]
                sequence = ""; header = ""
        
        ####################################
        #Check if new proteins not included in the cdhit clusters are found in collection:
        if count:
            sys.stderr.write("!!! The cdhit clusters does not cover all accession numbers, errors: %s" % count)
        
        for sequence,v in temp_collection.items():
        
            ####################################
            #Check if sequence is too short:
            if len(sequence.replace("X","").replace("_","")) < 20:
                continue
                
            ####################################
            #Ignore fragments:
            #if sum(["(fragment)" in h for h in v]) and not "fragments" in beta:
            #    continue
            
            ####################################
            #Loop over each of the headers:
            merge_head = [{} for f in range(len(meta_cat))]; cd_cluster = set()
            for header,cd in v.items():
                acc_meta = header.replace('>','').replace("&","_").split(sep)
                cd_cluster.add(cd[0]) if cdhit else None
                
                ####################################
                #Loop over each element of current header:
                for n in range(len(meta_cat)):
                
                    ####################################
                    #Ignore emtry strings:
                    if acc_meta[n]:
                                    
                        ####################################
                        #Exclude ignored entries:
                        if acc_meta[0] not in ignore:
                        
                            ####################################
                            #Merge headers of identical sequences:
                            for key in acc_meta[n].split("+"):    
                                merge_head[n][key] = 1
                                
            ####################################
            #Add family to subfamily and function:
            if use_common_family:
                family = common
            else:
                family = "+".join(sorted(merge_head[meta_cat.index("Family")]))
            if "Subfam" in meta_cat and merge_head[meta_cat.index("Subfam")]:
                merge_head[meta_cat.index("Subfam")] = {family+":"+"&".join([f for f in sorted(merge_head[meta_cat.index("Subfam")],key=lambda x:family_sort(x,start=0))]):1}

            ####################################
            #Remove subfamily to include family:
            if "raw_function" in meta_cat and merge_head[meta_cat.index("raw_function")]:
                merge_head[meta_cat.index("raw_function")] = {family+":"+"&".join([f.split(":",1)[-1] for f in sorted(pool_pool(dict((family+":"+k,1) for k in merge_head[meta_cat.index("raw_function")])),key=lambda x:family_sort(x,start=0))]):1}
            
            ####################################
            #Choose representative sequences:
            if cdhit:
                c = "+".join(sorted(cd_cluster))
                
                for acc in sorted(merge_head[meta_cat.index("Accession")]):
                    override_flag = False
                    if acc not in cdhit:
                        sys.stderr.write(" !  The cdhit cluster files are missing an accession number: %s\n" % acc)
                    elif cdhit[acc][1]:
                        cd_seq[c] = [sequence,acc,cdhit[acc][1],bool("Subfam" in meta_cat and merge_head[meta_cat.index("Subfam")]),bool("raw_function" in meta_cat and merge_head[meta_cat.index("raw_function")])]
                    
                if merge_head[meta_cat.index("Accession")]:
                    if c not in cd_merge:
                        cd_merge[c] = [{} for f in range(len(meta_cat))]
                    for n in range(len(merge_head)):
                    
                        ####################################
                        #Mark the accession numbers being representatives:
                        if n == meta_cat.index("Accession"):
                            for k in merge_head[n]:
                                cd_merge[c][n][k] = True if k in cdhit and cdhit[k][1] else False #Error recently fixed!
                        else:
                            for k in merge_head[n]:
                                if k not in cd_merge[c][n]:
                                    cd_merge[c][n][k] = 0
                                cd_merge[c][n][k] += 1
                        
            else:
                for acc in sorted(merge_head[0]):
                    if sequence not in collection:
                        collection[sequence] = acc
                    if acc not in meta_data_acc:
                        meta_data_acc[acc] = dict(zip(meta_cat,merge_head))
                    else:
                        for k,v in dict(zip(meta_cat,merge_head)).items():
                            for i,j in v.items():
                                if i not in meta_data_acc[acc][k]:
                                    meta_data_acc[acc][k][i] = 0
                                meta_data_acc[acc][k][i] += 1
                            
    ####################################
    #Create the final collection:
    if cdhit:
        for c,v in sorted(cd_seq.items()):
            if v[0] in collection:
                print ("ERROR",v[0])
            acc = v[1]
            collection[v[0]] = acc
            if acc not in meta_data_acc:
                meta_data_acc[acc] = dict(zip(meta_cat,cd_merge[c]))
            else:
                for k,v in dict(zip(meta_cat,cd_merge[c])).items():  #Error recently fixed!
                    for i,j in v.items():
                        if i not in meta_data_acc[acc][k]:
                            meta_data_acc[acc][k][i] = 0
                        meta_data_acc[acc][k][i] += 1
            
    ####################################
    #Return the meta data along with collection:
    return (collection,meta_data_acc)

def format_meta(collection, fasta_meta, ncbi = {}, meta_cat = []):
    '''
    Identify additional meta data relevant information for the proteins of the collection
    '''
    ####################################
    #For each path load the fasta file:
    kingdom_ranks = ["species", "genus", "family", "order", "subclass", "class", "subphylum", "phylum", "subkingdom", "kingdom", "superkingdom"]
    find_relevant = {"Aspergillaceae":"Asper","Archaea":"Archa","Fungi":"Fungi","Eukaryota":"Eukar","Bacteria":"Bacte","Other":"Other"}
    
    ####################################
    #Loop over each accession to index it and check for redundancy:
    for seq,acc in collection.items():
        
        ####################################
        #Obtaining taxonomical relationship, by taking the entry tax and look it up:
        fasta_meta[acc]["Taxonomy"] = {}
        for item in list(find_relevant.values()):
            fasta_meta[acc]["Taxonomy"][item] = 0
        species_pool = {}
        
        if "Taxid" in meta_cat:
            for tax in fasta_meta[acc]["Taxid"]:
                valid_rank = ""; tax_flag = True
                species,ncbi,temp = kingdom(tax,ncbi,preloaded=False)

                ####################################
                #Find the class/target tax of the entry:
                for r in range(2,11):
                    if species[r] in find_relevant:
                        fasta_meta[acc]["Taxonomy"][find_relevant[species[r]]] += 1
                        tax_flag = False
                if tax_flag:
                    fasta_meta[acc]["Taxonomy"]["Other"] += 1
                for r in range(5,11):
                    if not valid_rank and species[r]:
                        valid_rank = species[r]+"_%s" % kingdom_ranks[r]
                        species_pool[valid_rank] = 1
            
        fasta_meta[acc]["Classes"] = species_pool
        fasta_meta[acc]["Add"] = {}
        #fasta_meta[acc]["Add"]["Fragment"] = sum(["(fragment)" in f.lower() for f in fasta_meta[acc]["Name"]])
        fasta_meta[acc]["Add"]["Patent"] = sum(["_patent_" in f.lower() for f in fasta_meta[acc]["Name"]]) if "Name" in meta_cat else 0
        
        ####################################
        #Make dual meta data:
        for meta_item in ["raw_function","Subfam"]:
            if meta_item in meta_cat and fasta_meta[acc][meta_item]:
                dual_meta = {}
                fam_span = set((k.split(":")[0]) for k in fasta_meta[acc][meta_item])
                unique_ones = pool_pool(fasta_meta[acc][meta_item])
                for fam in fam_span:
                    fam_meta = dict((k.split(":",1)[-1],v) for k,v in unique_ones.items() if k.split(":")[0] == fam)
                    dual_meta.update({fam+":"+"+".join([f for f in sorted(fam_meta,key=lambda x:family_sort(x,start=0))]):min(fam_meta.values())})
                fasta_meta[acc][meta_item] = dual_meta
            
    return (collection, fasta_meta)

def obtain_peptides(collection,n_mer=8,ambiguous=2,use_tandem=False,single=False,use_positions=False,limit={},beta={}):
    '''
    Additional features are ambi in peptides, ignore peptides having "_" and numbering of peptides instead of strings!
    '''
    ####################################
    all_peptides = {}; singles = {}; include_all = True if len(collection) == 1 else False
    single = True if include_all else single
    seq_peptides = dict((acc,{}) for seq,acc in collection.items()) if not single else {}; tandem = {}
    combinations = list(list(k) for k in itertools.combinations(range(n_mer), ambiguous))
    
    for seq,acc in sorted(collection.items()):
            
        ####################################
        #Loop over each peptides of the protein:
        seq = [a for a in ("_"*ambiguous)+seq+("_"*ambiguous)]
        for n in range(len(seq)-n_mer+1):
            
            ####################################
            #Create each valid peptide:
            pep_original = seq[n:n+n_mer]
                
            ####################################
            #Create all combinations of X-variants of the peptide:
            for x_list in combinations:
                pep_list = pep_original[:]
                for x in x_list:
                    pep_list[x] = "X"
                xed = "".join(pep_list)
                
                ####################################
                #Continue if no invalid positions in peptide:
                if "_" not in xed and (not limit or xed in limit):
                
                    ####################################
                    #Obtain all peptides if a single protein is in the collection:
                    if single:
                        if xed not in seq_peptides: 
                            seq_peptides[xed] = n
                        elif use_tandem:
                            if xed not in tandem:
                                tandem[xed] = [n]
                            else:
                                tandem[xed].append(n)   
                    else:
                    
                        ####################################
                        #To save RAM, the peptides can be represented as ints instead:
                        if use_positions:
                            if xed not in all_peptides:
                                all_peptides[xed] = len(all_peptides)
                            xed = all_peptides[xed]
                            
                        ####################################
                        #Disregard peptides only found once among proteins in group:
                        if xed not in singles:
                            singles[xed] = [acc,n]
                        else:
                            if not singles[xed] or singles[xed][0] != acc:
                                if singles[xed]:
                                    seq_peptides[singles[xed][0]][xed] = singles[xed][1]
                                    singles[xed] = False
                                seq_peptides[acc][xed] = n

    ####################################
    if single:
        return seq_peptides,tandem
        
    return seq_peptides

    
def add_kingdom(meta,ncbi):
    '''
    Obtain the taxonomy of a protein based in NCBI Taxid
    '''

    kingdom_ranks = ["species", "genus", "family", "order", "subclass", "class", "subphylum", "phylum", "subkingdom", "kingdom", "superkingdom"]
    find_relevant = {"Aspergillaceae":"Asper","Archaea":"Archa","Fungi":"Fungi","Eukaryota":"Eukar","Bacteria":"Bacte","Other":"Other"}

    for acc,current in sorted(meta.items()):
    
        ####################################
        #Obtaining taxonomical relationship, by taking the entry tax and look it up:
        species_pool = {}
        current["Taxonomy"] = dict((v,0) for k,v in find_relevant.items())
        if "Taxid" in current:
            for tax in current["Taxid"]:
                valid_rank = ""; tax_flag = True
                species,ncbi,temp = kingdom(tax,ncbi,preloaded=False)

                ####################################
                #Find the class/target tax of the entry:
                for r in range(2,11):
                    if species[r] in find_relevant:
                        current["Taxonomy"][find_relevant[species[r]]] += 1
                        tax_flag = False
                if tax_flag:
                    current["Taxonomy"]["Other"] += 1
                for r in range(5,11):
                    if not valid_rank and species[r]:
                        valid_rank = species[r]+"_%s" % kingdom_ranks[r]
                        species_pool[valid_rank] = 1
            
        current["Classes"] = species_pool
        current["Add"] = {}
        current["Add"]["Fragment"] = sum(["(fragment)" in f.lower() for f in current["Name"]]) if "Name" in current else 0
        current["Add"]["Patent"] = sum(["_patent_" in f.lower() for f in current["Name"]]) if "Name" in current else 0
        
        meta[acc] = current
        
    return meta
    
def parallel_predict(content):
    """
    Predict in parallel.
    """
    collection, opts = tuple(content); results = {}
    for seq,acc in collection.items():
        results[acc] = predict(seq,opts)
        
    return results
    
def obtain_groups(seq_peptides,n_jobs=1,cc=9):
    '''
    Determine if the protein groups should be obtained using one or several processors.
    '''

    #Additional Improvement of RAM usage for large families in several cores: Make the peptides into n_jobs dict having all accessions but only peptides in relevant accession (for each q): peptides[q][acc][pep] = n
    print ("### Current usage of RAM:\t%s MB" % (round(ram()*1000)))
    
    if 0 and n_jobs != 1:
        ####################################
        #BETA: Reduce the size of the dict send to each processor (limited to 2 GB-pickle size):
        temp = dict((a,{}) for a in seq_peptides); pos = {}
        for acc,v in seq_peptides.items():
            for pep,n in v.items():
                if pep not in pos:
                    pos[pep] = len(pos)
                temp[acc][pos[pep]] = n
        seq_peptides = temp; temp = {}; pos = {}
    
    ####################################
    #Normal run if only a single core is specified:
    if n_jobs == 1:
        links = CUPP_clustering(seq_peptides,n_jobs=n_jobs,cc=cc)    
    else:
        
        ####################################
        #Divide the task and assemble the part:
        with Pool(processes=n_jobs) as pool:
            condensed_matrix_parts = pool.map(CUPP_clustering, ({"seq_peptides":seq_peptides,"q":q,"n_jobs":n_jobs,"cc":cc} for q in range(n_jobs)))
        condensed_matrix = np.concatenate(condensed_matrix_parts)

        ####################################
        #Simulate pool as control:
        #condensed_matrix = []
        #for q in range(n_jobs):
        #    condensed_matrix += CUPP_clustering(seq_peptides,q=q,n_jobs=n_jobs,cc=cc)    
        links = linkage(condensed_matrix,method="ward", optimal_ordering = True)
            
    return links

def CUPP_clustering(seq_peptides,q=0,n_jobs=1,mat=[],cc=9,debug=False):
    '''
    Cluster proteins of a full or part of a protein family.
    '''
    ####################################
    #Multi-processing:
    if "n_jobs" in seq_peptides:
        q = seq_peptides["q"]
        n_jobs = seq_peptides["n_jobs"]
        seq_peptides = seq_peptides["seq_peptides"]
    
    ####################################
    #Specify orders and sets of accessions:
    n = len(seq_peptides)
    ran = range(n)
    combi = int((n*n + n - n*(n+1)/2 - n))
    accession_order = dict(zip(ran,sorted(seq_peptides)))
    sets = dict(zip(ran,list(set(v) for k,v in sorted(seq_peptides.items()))))
    max_positions = dict((k,len(set(v.values()))) for k,v in seq_peptides.items())
    length = int(combi/n_jobs)
    dif = int(round(((combi/n_jobs - int(combi/n_jobs))*n_jobs if q and q+1 == n_jobs else 0)))
                
    ####################################
    #Create the distance matrix and potential for multiprocessing:
    condensed_matrix = np.zeros(length+dif)
    if q == 0 and n_jobs != 1:
        print ("### RAM usage during clustering: %s" % (ram()))
        
    for i in ran:
        for j in ran:
            if i < j:
                no = int((i*n + j - i*(i+1)/2 - i - 1)-q*length)
                if no < 0 or length+dif <= no:
                    continue
                common_peptides = sets[i]&sets[j]
                max_conserved_length = max(max_positions[accession_order[i]], max_positions[accession_order[j]])
                positions_seq_i = {seq_peptides[accession_order[i]][pep] for pep in common_peptides}
                positions_seq_j = {seq_peptides[accession_order[j]][pep] for pep in common_peptides}
                condensed_matrix[no] = (1-( (len(positions_seq_i)+len(positions_seq_j))/(2*max_conserved_length) ))**cc

    ####################################
    #Return the part of the condensed distance matrix for assembling:
    if n_jobs != 1:
        return (condensed_matrix)
    links = linkage(condensed_matrix,method="ward")
            
    return(links)

def arrange_meta(groups,meta,header=[]):
    '''
    Auto inspect the groups by merging odd functional strings in pool_pool
    '''
    
    ####################################
    CUPP_meta = {}
    for g,members in groups.items():
        if g not in CUPP_meta:
            CUPP_meta[g] = [{} for f in range(len(header))]
        for acc in members:
            for category,string in meta[acc].items():
                if string:
                    string = "+".join(sorted(string.split("+"),key=lambda x: family_sort(x)))
                    string = "&".join(sorted(string.split("&"),key=lambda x: family_sort(x)))
                    if string not in CUPP_meta[g][header.index(category)]:
                        CUPP_meta[g][header.index(category)][string] = 0
                    CUPP_meta[g][header.index(category)][string] += 1
    return CUPP_meta
    
def scoring(set_one,set_two,pos_one,pos_two,all_one,all_two,cc=3):
    '''
    Calculating the score for the relationship between the sets of peptides
    '''
    ####################################
    #Locate the position of each common peptide (only relevant if ambiguous is not 0):
    i_pos,j_pos = set(),set()
    for p in set_one & set_two:
        i_pos.add(pos_one[p])
        j_pos.add(pos_two[p])

    ####################################
    #Calculation og inspection of the scores:
    if all_one and all_two:
        score = (1-( (len(i_pos)+len(j_pos))/(2*max(all_one,all_two)) ))**cc
    else:
        score = 1

    return (score)
   
def check_existence(powder_json,setting):
    '''
    Check if desired model-settings exists for the current protein family
    '''
    common = os.path.basename(powder_json).replace("_powder.json","")
    existing = set(); related = False
    if os.path.exists(powder_json):
        try:
            with open(powder_json) as json_data:
                powder = json.load(json_data)
                existing = set((s) for s in powder if check_settings(s,setting))
                for setting in existing:
                    if "related" in powder[setting]:
                        related = True
                powder = {}
        except:
            sys.stderr.write(" !  The saved data was corrupted and was deleted for %s (setting:%s)!\n" % (common,setting))
            handle_out = open(powder_json,'w')
            json.dump({},handle_out,indent=3)
            handle_out.close()
    if existing:
        print("### Already processed %s %s: %s" % (common,"including relatedness!" if related else "",existing))
    return(existing)
        
def dendro(links,labels,algo="single",cc=0.2,name="",fil="",show=False,note={},horisontal=False):
    '''
    Make a dendrogram of distance matrix
    '''
    import matplotlib.pyplot as plt
    if 3000 < len(links):
        sys.setrecursionlimit(len(links))
        
    #from matplotlib import rcParams
    #    rcParams["lines.linewidth"] = 0.1
        
    ####################################    
    #Save the data for later:
    if show or not fil:
    
        ####################################
        #Set the size of the plot window in inches and also resolution and size-ratio:
        print('### Dendrogram with "%s" and cc at %s' % (algo,cc))
        fig = plt.figure()#figsize=(20,7))
        ax = fig.add_subplot(1, 1, 1)
        fig.canvas.set_window_title("Dendrogram %s" % (name.replace(".json","")))
        
        ####################################
        #The font size of the labels:
        ax.tick_params(axis='y', which='major', labelsize=8)
        
        ####################################
        #Decide the orientation of the plot:
        if horisontal:
            plt.subplots_adjust(left=0.03, right=0.99, top=0.99, bottom=0.2)
            ddata = dendrogram(links,ax=ax,labels=labels,color_threshold=cc)
        else:
            plt.subplots_adjust(left=0.01, right=0.8, top=0.99, bottom=0.05)
            ddata = dendrogram(links,ax=ax,orientation="left",labels=labels,color_threshold=cc)
            
        ####################################
        #If all nodes are below one, set to one:
        xmin, xmax = plt.xlim()         
        if xmin < 1:
            plt.xlim(xmin=1)
        
        ####################################
        #Special note saved with plot:
        if 0:#note:
            if isinstance(note,dict):
                print ("### Ignored accessions:")
                print("\n".join(["### %s" % (i) for i,j in (sorted(note.items()))]))
            else:
                print (note)

    ####################################
    #Save for later inspection:
    if fil:
        handle_out = open(fil,'w')
        json.dump({"links":list([list(e) for e in links]),"labels":labels,"algo":algo,"cc":cc,"horisontal":horisontal,"name":name,"note":note},handle_out,indent=3)
        handle_out.close()
        
    ####################################
    #Or plot it now and display later:
    if show:
        plt.show()    

def hex_code_colors():
    n1, n2, n3 = random.randrange(0,256),random.randrange(0,256),random.randrange(0,256)
    
    #Avoid yellow colors:
    while (200 < n1 and 200 < n2):
        n1, n2, n3 = random.randrange(0,256),random.randrange(0,256),random.randrange(0,256)
        
    a = hex(n1)[2:]
    b = hex(n2)[2:]
    c = hex(n3)[2:]
    if len(a)<2:
        a = "0" + a
    if len(b)<2:
        b = "0" + b
    if len(c)<2:
        c = "0" + c
    z = a + b + c

    return "#" + z.upper()
    
def gen_hex_color(group,custom_colors):

    group = str(group).strip("*")
    if group not in custom_colors:
        custom_colors[group] = hex_code_colors()
        
    return (custom_colors[group])

def gen_color(group,custom_colors):
    '''
    Generation of a random color if group do not have a customized color
    '''
    import matplotlib.pyplot as plt
    if group not in custom_colors:
        if len(custom_colors) <= 20:
            custom_colors[group] = plt.cm.tab20b(len(custom_colors)/20)
        else:
            custom_colors[group] = plt.cm.gist_ncar(np.random.random())
        
    return (custom_colors[group])
    
def getNewick(node, newick, parentdist, leaf_names):
    '''
    Convert dendrogram to Newick tree format for iTOL interaction.
    '''
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
    
def itol(types,meta_merge,general_name,header,meta,newick_path="",labels=[],links=[],custom_colors = {},path = "CUPP"):
    
    fams = {}
    for cupp,v in sorted(meta_merge.items(),key=lambda x: family_sort(x[0])):
        fam = cupp.split(":")[0]
        fams[fam] = 1

    ####################################
    #Create all the colors for cupp groups up front:
    for fam in fams:
        for i in range(1,500):
            for j in range(1,20):
                cupp = fam+":"+str(i)+"."+str(j)
                col = gen_hex_color(cupp,custom_colors)

    ####################################
    #Check if multiple domains of the same family is found in one protein:
    already = {}
    for cupp,v in sorted(meta_merge.items(),key=lambda x: family_sort(x[0])):
        for acc,s in v[header.index("Accession")].items():
            real_acc = acc.split(":")[0]
            if real_acc not in already:
                already[real_acc] = 0
            else:
                already[real_acc] += 1
            
    if "label" in types:
        label_handle = open("%s/itol/%s_label.txt" % (path,general_name),"w")
        label_handle.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t%s_label\nCOLOR\t#ff0000\nDATA" % general_name)
        for cupp,v in sorted(meta_merge.items(),key=lambda x: family_sort(x[0])):
            col = gen_hex_color(cupp,custom_colors)
            for acc,s in v[header.index("Accession")].items():
                real_acc = acc.split(":")[0]
                if real_acc in meta:
                    for a in meta[real_acc]["Accession"]:
                        a = a.split(":")[0]
                        label_handle.write("\n%s\t%s\t%s" % (real_acc,col,cupp+"_"+a+"_"+str(s)+("_DUAL_DOMAIN %s" % (already[real_acc]) if already[real_acc] else "")))
                label_handle.write("\n%s\t%s\t%s" % (real_acc,col,cupp+"_"+real_acc+"_"+str(s)+("_DUAL_DOMAIN %s" % (already[real_acc]) if already[real_acc] else "")))
        label_handle.close()
        
    if "text" in types:
        text_handle = open("%s/itol/%s_text.txt" % (path,general_name),"w")
        text_handle.write("DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,%s_text\nCOLOR,#ff0000\nDATA" % general_name)
        for cupp,v in sorted(meta_merge.items(),key=lambda x: family_sort(x[0])):
            col = gen_hex_color(cupp,custom_colors)
            if cupp.split(":")[1][0] != "0":
                for acc,s in v[header.index("Accession")].items():
                    real_acc = acc.split(":")[0]
                    if real_acc in meta:
                        for a in meta[real_acc]["Accession"]:
                            a = a.split(":")[0]
                            text_handle.write('\n%s,%s,-1,%s,bold,3,0' % (real_acc,cupp,col))
                    text_handle.write('\n%s,%s,-1,%s,bold,3,0' % (real_acc,cupp,col))
                    
        text_handle.close()
        
    if "other" in types:
        for word in ["raw_function","Taxonomy","Name","Subfam"]:
            if word in types or 1:
                text_handle = open("%s/itol/%s_%s.txt" % (path,general_name,word),"w")
                text_handle.write("DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,%s_%s\nCOLOR,#ff0000\nDATA" % (general_name,word))
                for acc,v in sorted(meta.items()):
                    if word in v:
                        text = "+".join([f for f,s in sorted(v[word].items()) if s])
                        text = text.replace(":","_") if word == "Subfam" else text
                        if text:
                            col = gen_hex_color(text,custom_colors)
                            text_handle.write('\n%s,%s,-1,%s,bold,3,0' % (acc.split(":")[0],text,col))
                text_handle.close()

    if newick_path:
        ####################################
        #Save den dendrogram as Newick format for handling e.g. in iTOL:
        tree = to_tree(links,False)
        handle = open(newick_path,"w")
        handle.write(getNewick(tree, "", tree.dist, labels))
        handle.close()

def transfer_ref(old,new):

    #!!!The cdhit clusters of the old models most be included in the CUPPlibrary to trace group after massiv rearrangements!!!

    ####################################
    #Identify the most central member of each group:
    new_pool = {}; fam = set()
    for cupp,v in sorted(new.items()):
        maxi = max(v[0].values())
        new_pool[[acc for acc,s in sorted(v[0].items()) if s == maxi][0]] = [cupp,v[0]]
        fam.add(cupp.split(":")[0].strip("x"))
        
    old = dict((cupp,v) for cupp,v in old.items() if cupp.split(":")[0].strip("x") in fam)
        
    ####################################
    #Find direct transfers of CUPP groups from one clustering to another:
    taken = {}; cdhit = {}; index = {}; transfer = {}; emtry = set(); used_acc = set()
    for cupp,v in sorted(old.items()):
    
        old_type_pool = v[0]; maxi_global = max(old_type_pool.values())
        missing = dict((acc,s) for acc,s in old_type_pool.items() if s < maxi_global*0.5)
        while old_type_pool and cupp not in taken:
            old_type_pool = dict((acc,s) for acc,s in v[0].items() if acc not in missing)
            
            if not old_type_pool:
                continue
                
            maxi = max(old_type_pool.values())
            old_sum = sum(v[0].values())
            old_type = [f for f,s in sorted(old_type_pool.items()) if s == maxi][0]
            
            old_type = [acc for acc in sorted(cdhit[acc]) if acc in new_pool][0] if cdhit and acc in cdhit else old_type
            if old_type in new_pool:
                taken[cupp] = new_pool[old_type][0]
                
                ####################################
                #Check 50%:
                shared = set(v[0]) & set(new_pool[old_type][1])
                other_in_new = set(new_pool[old_type][1]) - shared
                old_sum_transferred = sum([s for acc,s in new_pool[old_type][1].items() if acc in shared])
                other_sum_in_new = sum([s for acc,s in new_pool[old_type][1].items() if acc in other_in_new])
                
                if old_sum*0.5 < old_sum_transferred and other_sum_in_new < old_sum*0.5:
                    print ("###",old_type,cupp,"-->",new_pool[old_type][0],old_sum,old_sum_transferred,other_sum_in_new)
                else:
                    if maxi_global == 0:
                        print ("# #",old_type,cupp,"-->",new_pool[old_type][0],old_sum,old_sum_transferred,other_sum_in_new)
                    else:
                        print ("! !",old_type,cupp,"-->",new_pool[old_type][0],old_sum,old_sum_transferred,other_sum_in_new)
                transfer[new_pool[old_type][0]] = cupp
                used_acc.update(shared)
            else:
                missing[old_type] = 0
                
    ####################################
    #Identify the most central member of each group:
    old_pool = {}
    for cupp,v in sorted(old.items()):
        old_members = dict((acc,s) for acc,s in v[0].items() if acc not in used_acc and s)
        if 3 <= len(old_members) and cupp.strip("x") in fam:
            maxi = max(old_members.values())
            for acc,s in old_members.items():
                if maxi*0.5 <= s:
                    old_pool[acc] = [cupp,old_members]
            
    ####################################
    #Check for any new CUPP groups which have not been placed in the existing system:
    new_ones = {}
    for cupp in [c for c in sorted((set(new) - (set(new) & set(transfer))),key=lambda x: int(x.split(":")[1].split(".")[0]))]:
    
        ####################################
        #Identify members and their rank:
        new_members = new[cupp][0]
        
        ####################################
        #Identify the highest value:
        maxi_global = max(new_members.values())
        
        ####################################
        #Consider members with less than 50% removed:
        missing = dict((acc,s) for acc,s in new_members.items() if s < maxi_global*0.5)
        
        ####################################
        #While a representative protein has not been tested continue to test relationship:
        while new_members and cupp not in transfer:
        
            ####################################
            #Remove members which have already had a try:
            new_members = dict((acc,s) for acc,s in new_members.items() if acc not in missing)
            
            ####################################
            #If no member of initiate the retation, there is none:
            if not new_members:
                continue
                
            ####################################
            #Identify the score of the current highest scoring member:
            maxi = max(new_members.values())
            
            ####################################
            #Identify the member with the highest score:
            new_type = [f for f,s in sorted(new_members.items()) if s == maxi][0]
            
            ####################################
            #Determine the sum of all members in group:
            new_sum = sum(new[cupp][0].values())
            
            ####################################
            #Future implementation using cdhit to extend members in case of major rearrangements:
            new_type = [acc for acc in sorted(cdhit[acc]) if acc in old_pool][0] if cdhit and acc in cdhit else new_type
            
            ####################################
            #Check if the current representative in found in the highest ranking remaining member of an old group:
            if new_type in old_pool:
                taken[cupp] = old_pool[new_type][0]
                
                ####################################
                #Identify accession numbers shared between the old and the new group:
                shared = set(new[cupp][0]) & set(old_pool[new_type][1])
                
                ####################################
                #Identify how many new members also found in the new group:
                other_in_new = set(old_pool[new_type][1]) - shared
                
                ####################################
                #How much of the score is transferred from the old group into the new:
                old_sum_transferred = sum([new[cupp][0][acc] for acc,s in old_pool[new_type][1].items() if acc in shared])
                
                ####################################
                #How much new score is also found in the new group:
                other_sum_in_new = sum([new[cupp][0][acc] for acc,s in old_pool[new_type][1].items() if acc in other_in_new])
                
                ####################################
                #Evaluate how good the relationship is between the old and the new group:
                if new_sum*0.5 < old_sum_transferred and other_sum_in_new < new_sum*0.5:
                    print ("   ",new_type,cupp,"-->",old_pool[new_type][0],new_sum,old_sum_transferred,other_sum_in_new)
                else:
                    if maxi_global == 0:
                        print (" # ",new_type,cupp,"-->",old_pool[new_type][0],new_sum,old_sum_transferred,other_sum_in_new)
                    else:
                        print (" ! ",new_type,cupp,"-->",old_pool[new_type][0],new_sum,old_sum_transferred,other_sum_in_new)

                ####################################
                #Rename the CUPP group to indicate a recent split:
                secondary = old_pool[new_type][0]; no = 1
                while secondary in taken or cupp not in transfer:
                    secondary = old_pool[new_type][0].split(".")[0]+"."+str(int(old_pool[new_type][0].split(".")[1])+no)
                    if secondary not in taken and cupp not in transfer:
                        transfer[cupp] = secondary
                    no += 1
            missing[new_type] = 0
    
        ####################################
        #If a group is found, do not continue
        if cupp in transfer:
            continue
    
        ####################################
        #Identify avaiable number:
        fam = cupp.split(":")[0].strip("x")
        max_old = [int(c.split(":")[1].split(".")[0]) for c in sorted(old,key=lambda x: int(x.split(":")[1].split(".")[0]),reverse=True) if c.split(":")[0].strip("x") == fam]
        max_old = max_old[0] if max_old else 0
        
        ####################################
        #Assign a new number of group did not fit in existing groups:
        if fam not in new_ones:
            new_ones[fam] = set()
        new_ones[fam].add(max(max_old+1,(max(new_ones[fam]) if new_ones[fam] else 0)+1))
        new1 = fam+":"+str(max(new_ones[fam]))+".1"
        
        transfer[cupp] = new1
        print ("NEW",cupp,"-->",new1,max_old)
                
    return transfer

    
def check_settings(current_setting,target_setting,full=False):
    '''
    Check if the setting is consistent, could be modified to merge powder json of differnent settings:
    '''
    ####################################
    if full:
        if current_setting[2:].split("_")[0] != target_setting[2:].split("_")[0]:
            return (0)
        if current_setting[0] != target_setting[0]:
            return (0)
        if current_setting[1] != target_setting[1]:
            return (0)
        return (1)
        
    else:
        if current_setting[2:].split("_")[0] != target_setting[2:].split("_")[0]:
            return (0)
        return (1)
        
def determine_paths(common_list,deposit,query="",fasta_extension=".faa",dbcan_folder="",cdhit_folder="",redundancy_path="",domain_off=False,domains_only=False,cdhit=0.9, beta_options={}):

    ####################################
    #Create CUPP Library or predict several collections (to merge meta data of identical sequences from different families):
    domain_list = []; original_list = []; cdhit_list = []; used_list = []
    for common in common_list:

        ####################################
        #Check if the collections exists:
        original = "%s/%s%s" % (deposit,common,fasta_extension) if not query else query
        protein_domains = "%s/%s%s.dbcan" % (dbcan_folder,common,fasta_extension)
        if not domain_off:
            if not os.path.exists(original):
                sys.stderr.write("!0! The path for full length proteins does not exist: %s\n" % original)
                original = protein_domains
            if not os.path.exists(protein_domains):
                sys.stderr.write("!1! The path for domains does not exist: %s\n" % protein_domains)
                if os.path.exists(original):
                    protein_domains = original
                else:
                    sys.stderr.write("!2! The path for full length proteins does not exist: %s\n" % original)
                    
            if not os.path.exists(original) and not os.path.exists(protein_domains):
                sys.stderr.write("!2! No paths could be found for: %s\n" % common)
            
        original_list.append(original if not domains_only else protein_domains)
        domain_list.append(protein_domains)
        used_list.append(original if domain_off else protein_domains)
            
        ####################################
        #In case the path does not exist, proceeds without clusters:
        cd_dom = "" if domain_off else "_domain"
        cd_dom = "_domain" if domains_only else cd_dom
        cdhit_fil = "%s/%s%s%s%s.clstr" % (cdhit_folder,common,fasta_extension,cd_dom,("_%s" % cdhit))
        if not os.path.exists(cdhit_fil):
            cdhit_fil = "%s/%s%s%s.clstr" % (cdhit_folder,common,cd_dom,("_%s" % cdhit))

        cdhit_fil = cdhit_fil if not redundancy_path else redundancy_path
        if not os.path.exists(cdhit_fil):
            if cdhit != 0:
                sys.stderr.write("!3! The path for redundency clusters does not exist: %s\n" % cdhit_fil)
            cdhit_fil = ""        
        cdhit_list.append(cdhit_fil)  
        
    return (domain_list, original_list, used_list, cdhit_list)

def origin(original_list,cdhit_list,all_meta,common_list,header,sep="|"):
    '''
    Create one collection of multiple fastafiles to merge metadata of identical seqeunces in case of dual domains in one entry
    Also the CUPP group of an entry is recovered for validation of prediction
    '''
    meta_data_acc={}
    locate_accession = {}
    for no in range(len(original_list)):

        ####################################
        #Load the cdhit clustering:
        cdhit = {}
        if os.path.exists(cdhit_list[no]):
            cdhit = load_redundancy([cdhit_list[no]],common=common_list[no])
            
        ####################################
        #Find the representative accession of a cdhit cluster:
        cluster_num = dict((v[0],k) for k,v in cdhit.items() if v[1])
        
        ####################################
        #Loop over the accession and assign the relevant representative accession number:
        assign_cup = {}
        for cup in all_meta:
            fam,number = cup.split(":")
            agg = fam+":"+number.split(".")[0]
            if cup != fam:
                for acc in all_meta[cup]["Accession"]:
                    acc = acc.split(":")[0]
                    if acc not in assign_cup:
                        assign_cup[acc] = [set(),set()]
                    assign_cup[acc][0].add(cup)
                    assign_cup[acc][1].add(agg)
                    
                    if not cdhit:
                        locate_accession[acc] = [assign_cup[acc][0],assign_cup[acc][1]]
        
        for acc in cdhit:
            acc = acc.split(":")[0]
            if acc in assign_cup:
                locate_accession[acc] = [assign_cup[acc][0],assign_cup[acc][1]]
                
    return (locate_accession)

def create_table(all_meta, targets=[], meta_cat=False, silent=False):
    '''
    Display meta data overview of CUPP groups and format the meta data
    from lists into dictionaries
    '''

    if meta_cat:
        temp_meta = {}
        for cup, liste in sorted(all_meta.items(), key=lambda x: family_sort(x[0])):
            temp_meta[cup] = {}
            fam = cup.split(":")[0]
            if fam == cup:
                continue
            for no in range(len(liste)):
                temp_meta[cup][meta_cat[no]] = liste[no]
        all_meta = temp_meta

    # ---- Compact summary table ----
    if not silent and "None" not in targets:

        # KEEP ALL GROUPS, INCLUDING ONES WITH ZERO SCORE VALUES
        filtered_meta = dict(all_meta)

        # Determine valid columns
        valid_cols = []
        for col in meta_cat:
            for _, row in filtered_meta.items():
                val = row.get(col, "")
                if isinstance(val, dict):
                    meaningful = any(True for _ in val.keys())  # FIX: dictionary keys represent presence
                else:
                    meaningful = val not in ("", None, {}, 0)
                if meaningful:
                    valid_cols.append(col)
                    break

        if not filtered_meta:
            print("### CUPPgroup")
        else:
            # Compute dynamic widths
            max_len = max(len(name) for name in filtered_meta.keys())
            acc_width = 6  # digits right-aligned under this

            # HEADER — perfectly aligned
            header = (
                "### "
                + "CUPPgroup".ljust(max_len + 2)
                + "Accession".rjust(acc_width + 1)
                + "  Freq_sum"
            )
            print(header)

            # ROWS
            for cupp, row in sorted(filtered_meta.items(), key=lambda x: family_sort(x[0])):
                gpad = cupp.ljust(max_len + 2)

                # FIX: count members by NUMBER OF KEYS, not score!=0
                acc_dict = row.get("Accession", {})
                if isinstance(acc_dict, dict):
                    acc_count = len(acc_dict.keys())     # CORRECT FIX!
                else:
                    acc_count = 0

                apad = str(acc_count).rjust(acc_width)

                # Freq_sum
                fval = row.get("Freq_sum", "")
                if isinstance(fval, dict):
                    fval = f"{fval.get('Average','')}/{fval.get('sum','')}"

                print(f"### {gpad}{apad}  {fval}")

    return (all_meta)    
    
def rename(CUPPlist,meta_merge,meta_cat,fam,opportunists,assign_cupp,ignored_groups,start_positions,last=False):
    '''
    Determine the number of CUPP group based on how related the peptides of the groups are along with meta data:
    '''
    ####################################
    #Identify groups in which all proteins have been removed in the final clustering:
    ups = set()#; last = False
    if last:
        for cupp in CUPPlist:
            for acc in meta_merge[cupp][meta_cat.index("Accession")]:
            
                if acc not in ignored_groups:
                    ups.add(cupp)
    
    ####################################
    #In the case of two clusters gets merged, handle the renaming of groups:
    long_label = []; long_assign = []; new_name = {}; sub_name = {}
    for cupp in sorted(CUPPlist,key=lambda x:(rename_order(meta_merge[x],meta_cat,fam=True),rename_order(meta_merge[x],meta_cat))):
        
        ####################################
        #Determine name of subfamily:
        emtry = False if not last or cupp in ups else True
        new_cupp = assign_cupp[cupp][0]
        if new_cupp not in sub_name:
            sub_name[new_cupp] = [len(sub_name)+1,0]
            
        ####################################
        #Determine new cupp name:
        if last:
            sub_name[new_cupp][1] += 1
        new_name[new_cupp] = fam+":"+str(sub_name[new_cupp][0])+"."+str(sub_name[new_cupp][1])
        new_cupp = new_name[new_cupp]

        ####################################
        #Prepare remerging of meta data according to new groups:
        for acc in meta_merge[cupp][meta_cat.index("Accession")]:
            if acc not in ignored_groups:
                if start_positions <= opportunists[acc]:
                    long_assign.append(new_cupp)
                    long_label.append(acc)
            elif emtry:
                long_assign.append(new_cupp)
                long_label.append(acc)
                
    ####################################
    #Print out names of groups which almost made it:
    nuked = set(CUPPlist) - ups
    if ups and nuked:
        print ("--> Some CUPP groups (%s) does not have any members yet: %s" %(len(nuked),";".join([new_name[assign_cupp[cupp][0]] for cupp in sorted(nuked)])) )
    
    return long_assign, long_label

def rename_order(m,meta_cat,fam=False):
    '''
    During powderize, the CUPPgroups are being renamed in a user defined manner
    '''
    Subfam = 10**10; Members = 10*10
    if "Subfam" in meta_cat and m[meta_cat.index("Subfam")]:
        sub = [k for k,v in sorted(m[meta_cat.index("Subfam")].items(),key=lambda x:x[1],reverse=True)][0].split("+")[0].split(":")[-1]
        Subfam = int(sub) if sub.isdigit() else int(sub.split(":",1)[-1]) if sub.split(":",1)[-1].isdigit() else 0#sub if sub else Subfam
        Members = int([v for k,v in sorted(m[meta_cat.index("Subfam")].items(),key=lambda x:x[1],reverse=True)][0]) * -1
    Family = "+".join([f for f,v in sorted(m[meta_cat.index("Family")].items(),key=lambda x:x[1],reverse=True)]) if "Family" in meta_cat else ""
    PDB = len(m[meta_cat.index("PDB")]) * -1 if "PDB" in meta_cat else 0
    Function = sum(m[meta_cat.index("raw_function")].values()) * -1 if "raw_function" in meta_cat else 0
    Classes = len(m[meta_cat.index("Classes")]) * -1
    Accessions = len(m[meta_cat.index("Accession")]) * -1
    if fam:
        return str(Family)
    else:
        return Subfam,PDB,Function,Members,Classes,Accessions

def pool_groups(groups,labels,meta,meta_cat):
    '''
    Reduce the complexity of string in meta data of a CUPP group, without losing information
    '''
    
    ####################################
    #Debugging:
    if len(labels) != len(groups):
        print ("!!! Major error in pool_groups function")
    
    ####################################
    #Loop over the labels and create group lookups:
    pair = dict(zip(labels,groups))
    lookup_acc_group = {}; meta_merge = {}
    for acc,group in pair.items():
    
        ####################################
        #
        if acc in meta:
            
            ####################################
            #Lookup the group of an accession number:
            for accession in meta[acc]["Accession"]:
                if accession in lookup_acc_group and lookup_acc_group[accession][0] != group:
                    sys.stderr.write(" !  Accession number found in multiple groups!\n",accession,lookup_acc_group[accession],group)
                lookup_acc_group[accession] = [group,""]
            
            ####################################
            #Merge the meta data by groups:
            current = meta[acc]
            if group not in meta_merge:
                meta_merge[group] = [{} for f in range(len(current))]
            #for key,item in meta[labels[n]].items():
            for no in range(len(meta_cat)):
                
                ####################################
                #Ignore emtry strings and add all other to dict:
                if current[meta_cat[no]]:
                    for key,value in current[meta_cat[no]].items():
                        if key not in meta_merge[group][no]:
                            meta_merge[group][no][key] = 0
                        meta_merge[group][no][key] += value
            
    return (lookup_acc_group,meta_merge)
        
def relatedness(conserved,cup_cut=0.7):
    '''
    Determine how related CUPP groups are based on the peptides they share.
    '''

    ####################################
    #Calculate the summed square-frequencies within each group and the shared parts between groups:  
    connectivity = {}
    for pep in conserved:
        for group,score in conserved[pep].items():
            if group not in connectivity:
                connectivity[group] = {}
            for counter,s in conserved[pep].items():
                if counter not in connectivity[group]:
                    connectivity[group][counter] = 0
                connectivity[group][counter] += s*score
    
    ####################################
    #Create a matrix of the ratios of frequencies of shared peptides between CUPPgroups: 
    order = dict((no,sorted(connectivity)[no]) for no in range(len(connectivity)))
    pos = dict((v,k) for k,v in order.items())
    CUPPlist = [v for k,v in sorted(order.items(),key=lambda x:x[0])]
    
    ####################################
    #Create matrix of  similarities if any between CUPPgroups: 
    if 1 < len(connectivity):
        matrix = np.ones((len(connectivity),len(connectivity)))
        for group in connectivity:
            for counter,ratio in connectivity[group].items():
                matrix[pos[group]][pos[counter]] -= ratio/max(connectivity[group][group],connectivity[counter][counter])
                
        ####################################
        #Calculate distances:
        matrix = squareform(matrix)
        links = linkage(matrix,method="complete")
                
        ####################################
        #Locate clusters with shared peptide frequencies:
        assign_cup = dict(zip(CUPPlist,cut_tree(links,height=cup_cut)))
        #assign_sub = dict(zip(CUPPlist,cut_tree(links,height=cup_cut)))
    else:
        assign_cup = {sorted(CUPPlist)[0]:"Lone_Wolf"}
        #assign_sub = assign_cup.copy()
        
    return assign_cup,CUPPlist 
    
def powderize(collection,seq_pep,lookup_acc_group,min_size,n_mer,ambiguous,conserved_cut,last=False):
    '''
    Determine the conservation of peptides in the protein groups.
    '''
    
    ####################################
    #Connect the cluster of accessions to sequences:
    bean_clusters = {};
    for seq,acc in collection.items():
    
        ####################################
        #Create a lookup for the CUPPgroup of each accession number
        if acc in lookup_acc_group:
            group = lookup_acc_group[acc][0]
            if group not in bean_clusters:
                bean_clusters[group] = {}
            bean_clusters[group][seq] = acc
        
    ####################################
    #Ignore proteins of groups with less than a minimum of members:
    good_clusters = dict((group,seqs) for group,seqs in bean_clusters.items() if min_size <= len(seqs))
    
    ####################################
    #Looping over sequences of each group: 
    conserved = {}; opportunists = {}
    for group,collection in sorted(good_clusters.items()):
       
        ####################################
        #Obtain frequencies of peptides in current collection:
        cupp_peptides = {}
        for seq,acc in collection.items():
            for pep in seq_pep[acc]:
                if pep not in cupp_peptides:
                    cupp_peptides[pep] = 0
                elif not cupp_peptides[pep]:
                    cupp_peptides[pep] += 2/len(collection)
                else:
                    cupp_peptides[pep] += 1/len(collection)
        
        kept_peptides = set()
        for pep,score in cupp_peptides.items():
            if conserved_cut <= score:
                if pep not in conserved:
                    conserved[pep] = {}
                conserved[pep][group] = score
                kept_peptides.add(pep)
        
        ####################################
        #Identify proteins which have lost all peptides due to increased conservedness of peptides:
        for seq,acc in collection.items():
            opportunists[acc] = len(set((n) for pep,n in seq_pep[acc].items() if pep in kept_peptides))
                
    return opportunists, conserved, bean_clusters, good_clusters    

    
#####################################
### Experimental scripts/fuctions ###
#####################################

def cupp_align(query,target,n_mer = 8):

    '''
    Evaluation of protein comparison based on peptides
    '''

    #Ottain peptides of query protein:
    q_pep = obtain_peptides({query:"NEW1"})["NEW1"]
    
    #Obtain peptides of target protein:
    t_pep = obtain_peptides({target:"NEW1"})["NEW1"]
    
    #Identify peptides shared between the two proteins:
    shared = set(q_pep) & set(t_pep)
    
    #Create list of peptide pools:
    l_pep = [q_pep,t_pep]
    l_seq = [query,target]
    
    #Identify shared peptides in each position of both proteins:
    q_pos = {}; t_pos = {}; l_pos = [q_pos,t_pos]
    for no in range(2):
        for p,n in l_pep[no].items():
            n -= 2
            for o in range(2,8):
                if p not in shared or len(l_seq[no]) < n+o or p[o] == "X":
                    continue
                if n+o not in l_pos[no]:
                    l_pos[no][n+o] = set()
                l_pos[no][n+o].add(p)
                
    alternative = {}; taken = set()
    for i,q_pos_pep in sorted(q_pos.items()):
    
        ####################################
        #Skip position if no shared peptides:
        if not q_pos_pep:
            continue
            
        #Determine peptides shared between query position and target peptide pool.
        pep_share = set(q_pos_pep) & set(t_pep)
        
        #Determine the position of which the peptides shared with the position originates from:
        pos_share = set(n for pep,n in t_pep.items() if pep in pep_share and n not in taken)
        
        print (i,"-->",min(pos_share) if pos_share else "?",len(pep_share),pos_share,set(n for pep,n in t_pep.items() if pep in pep_share))
        taken.add(min(pos_share)) if pos_share else None
            
def pool_pool(raw_pool,keep_fam=False,reduce=True):    
    '''
    Locate common denominator in "+" a pool of functions:
    '''
    ####################################
    #Take one family at a time:
    fam_span = set((k.split(":")[0]) for k in raw_pool)
    all_pool = {}
    for fam in fam_span:
        fam_meta = dict((k.split(":",1)[-1],v) for k,v in raw_pool.items() if k.split(":")[0] == fam)

        ####################################
        #Count the total count of all non star functions in the pool:
        pool = {}
        for functions,score in fam_meta.items():
            for f in set(functions.replace("+","&").split("&")):
                if f not in pool:
                    pool[f] = 0
                pool[f] += score
                
        ####################################
        #If the same EC numbers are found as a dual domain keep only that:
        dual_pool = {}
        for k,v in fam_meta.items():
            dual_pool[len(dual_pool)] = set(f for f in k.replace("+","&").split("&") if f in pool)
        within = {}
        for i,n in dual_pool.items():
            for j,m in dual_pool.items():
                if n != m:
                    if len(n & m) == len(n):
                        within[i] = 1
        dual_pool = list((f) for n,f in dual_pool.items() if n not in within)
        
        ####################################
        #Handling the resulting pool if dual function:
        final_pool = {}
        if dual_pool and reduce:
            for functions,count in fam_meta.items():
                for d in dual_pool:
                    current = set(f for f in functions.replace("+","&").split("&") if f in pool)
                    if len(d & current) == len(current):
                        merge_pool = "&".join(sorted(d,key=lambda x:family_sort(x,start=0)))
                        #print (d,merge_pool)
                        if merge_pool not in final_pool:
                            final_pool[merge_pool] = 0
                        if len(d & current) == len(d):
                            final_pool[merge_pool] += count
                        
        ####################################
        #Handling resulting pool without dual function:
        if not final_pool:
            for functions,score in fam_meta.items():
                for f in set(functions.split("+")):
                    if f not in final_pool:
                        final_pool[f] = 0
                    final_pool[f] += score/len(functions.split("+"))
            
        ####################################
        #Make whole float numbers into integers:
        for k,v in final_pool.items():
            if int(v) == v:
                final_pool[k] = int(v)
            else:
                final_pool[k] = round(v,1)
                
        ####################################
        #Put back the family on the function:
        all_pool.update(dict((fam+":"+k,v) for k,v in final_pool.items()))
    
    return (all_pool)
    
def load_powder(json_fil,all_peptides,meta,meta_cat,current_setting,experimentals,mimimum_group=5,rel_cut=0.9,final_cut=0.2,power_factor=2,limit=0,mono_family=True,aa=6,special={}):
    '''
    Rename the CUPP groups according to meta data and exclude proteins which have too few conserved peptides of the group
    '''
    
    ####################################
    #Loop over each json file in a list:
    with open(json_fil) as json_data:
        powder = json.load(json_data)

    if not current_setting and len(powder) == 1:
        current_setting = sorted(powder.keys())[0]
        
    ####################################
    #Controlling the setting are suitable:
    accepted = set((s) for s in sorted(powder) if check_settings(s,current_setting))
    fam = os.path.basename(json_fil).replace("_CUPPpool.json","")
    
    ####################################
    #Check for programming errors/debugging:
    if len(accepted) != 1:
        print ("; ".join(powder.keys()))
        sys.exit(' !  There is different from one accepted powder json (%s)! "%s" \n%s' % (len(accepted),json_fil," & ".join(accepted)))
        
    ####################################
    #Add family the experimental families:
    if not set((s) for s in sorted(powder) if check_settings(s,current_setting,full=True)) and accepted:
        experimentals[fam] = "|".join(sorted(accepted))
    
    ####################################
    #Loop over accepted settings:
    CUPPsize = {}; all_meta = {}
    catalytic = {"H":16.5,"D":14,"E":12.5,"R":9.25,"K":8.75,"C":6.75,"Y":6.25,"S":4.75,"N":4.5,"G":4,"T":3.75,"Q":3,"F":1.5,"W":1.25,"L":1,"A":1,"P":0.25,"M":0.25,"V":0.25,"I":0.25}
    for setting,patterns in powder.items():
        if setting in accepted:
            
            ####################################
            #Recover position of the collections:
            meta_cat += patterns["meta_categories"][len(meta_cat):] if "meta_categories" in patterns else []
            
            ####################################
            #For each groups and its peptide frequency dict load it into all_peptides:
            print ("### Loading powder for (%s): %s" % (setting,os.path.basename(json_fil).replace("_powder.json","")))
            
            ####################################
            #Limit to peptides to be the top conserved peptides only:
            if limit:
                cupp = {}; new = {}
                for pep,freqs in patterns["peptides"].items():
                    for group,score in freqs.items():
                        if group not in cupp:
                            cupp[group] = {}
                        cupp[group][pep] = score
                for group,p in cupp.items():
                    inside = 0
                    for pep,score in sorted(p.items(),key=lambda x: x[1],reverse=True):
                        if inside <= limit:
                            if pep not in new:
                                new[pep] = {}
                            new[pep][group] = score
                            inside += score
                patterns["peptides"] = new
            
            ####################################
            #Construct the all peptide dictionary for fast prediction:
            for pep,freqs in patterns["peptides"].items():
            
                for group,score in freqs.items():
                    if final_cut < score:
                    
                        ####################################
                        #Decrease influence of linker and membrane bound regions
                        if "catalytic" in special:
                            ####################################
                            #Alternative amino acid bias reduction
                            pep_score = 0
                            for p in pep:
                                pep_score += catalytic[p] if p != "X" else 0
                            if pep_score <= aa:
                                continue
                        else:
                            if (pep.count("G") + bool(pep.count("A")) + bool(pep.count("P")) + bool(pep.count("V")) + bool(pep.count("L"))) == aa:
                                continue
                            
                        ####################################
                        #Save all peptides with the square-score for each CUPP group it is conserved in:   
                        if pep not in all_peptides:
                            all_peptides[pep] = {} 
                        if group in all_peptides[pep]:
                            sys.stderr.write("!!! ERROR IN LOAD: %s %s %s\n" % (group,all_peptides[pep][group],score))
                        all_peptides[pep][group] = score**power_factor
                        #NOTE: Using single tuple instead of dict of dicts for each peptides uses half the RAM but double the load time
                        
                        ####################################
                        #Determine summed square-score/frequency within the individual CUPPgroups:
                        if group not in CUPPsize:
                            CUPPsize[group] = {"sum":0,"peps":0}
                        CUPPsize[group]["sum"] += score**power_factor
                        CUPPsize[group]["peps"] += 1
                
            ####################################
            #Determine the sum of frequencies of the CUPPgroup:
            recover = [0,0,0,0]; related_functions = {}
            for cup,index in patterns["meta"].items():
                function = pool_pool(index[meta_cat.index("raw_function")]) if "raw_function" in meta_cat else {}
                if "related" in patterns and cup in patterns["related"] and (not function or patterns["related"][cup][cup] == 1) and "raw_function" in meta_cat:
                
                    ####################################
                    #In case multiple friends are found with the same relatedness:
                    friend_pool = {}
                    for friend,relatedness in patterns["related"][cup].items():
                        for f,s in patterns["meta"][friend][meta_cat.index("raw_function")].items():
                            if relatedness <= rel_cut:
                                if f not in friend_pool:
                                    friend_pool[f] = 0
                                friend_pool[f] += s
                    function = pool_pool(friend_pool)
                    function = dict((k,0) for k in function)
                
                ####################################
                #Determine the function of the family:
                for f,s in function.items():
                    if f not in related_functions:
                        related_functions[f] = 0
                    related_functions[f] += s
                if cup in CUPPsize:
                    function = {} if "CBM" in cup else function
                    index[-1].update(dict((k,v) for k,v in sorted(CUPPsize[cup].items())))
                    all_meta[cup] = index + [function]
                #else:
                #    sys.stderr.write(" !  All peptides in CUPPgroup %s has been removed: %s %s %s\n" % (cup,index[1],index[2],index[3]))
                if len([acc for acc,v in index[meta_cat.index("Accession")].items() if v]) < mimimum_group:
                    experimentals[cup] = len([acc for acc,v in index[meta_cat.index("Accession")].items() if v])
                    
            ####################################
            #If only a single function is known to the family assign that one:
            for cup,index in all_meta.items():
                if mono_family:
                    if not index[-1] and len(pool_pool(related_functions)) == 1 and "CBM" not in cup:
                        all_meta[cup][-1] = dict((k,0) for k,v in pool_pool(related_functions).items())
                recover[0] += 1 if "raw_function" in meta_cat and index[meta_cat.index("raw_function")] else 0
                recover[1] += 1 if index[-1] else 0
                recover[2] += 1
                
            meta.update(all_meta)
            if recover[2]:
                print ("### %s CUPP groups; %s%s with function (recovered %s%s)" % (fam,round(100*recover[0]/recover[2],1),"%",round(100*recover[1]/recover[2],1),"%"))
            else:
                sys.stderr.write("!!! No CUPP groups found in CUPP library\n")
                
def validate_meta(predicted,target,fam_overview,function_overview,cupp_overview,sub_overview,dbcan_domain,string="",lost_meta={},verbose=False,unknown=False,focus={}):
    '''
    Check if the corect EC is exactly the same, partially the same or among the correct
    '''
    bad_flag = False
    for tag in ["function","cup","sub","family","agg"]:

        ####################################
        #Ignore entries without the type of meta:
        if tag not in target or not target[tag]:
            continue
            
        ####################################
        #Determine if some hit are found in families which are not in target:
        predictions_outside_target_families = set()
        holding_families = set()
        for k,v in predicted[tag].items():
            if k not in target["Family"]:
                for f in v:
                    predictions_outside_target_families.add(f)
            elif v:
                holding_families.add(k)

        ####################################
        #Loop over each known family:
        for fam in target["Family"]:
        
            ####################################
            #Preparing set for sensitivity and precision determination:
            cupp = list((k) for k,v in target["cup"].items() if k.split(":")[0] == fam)
            cupp = "+".join(sorted(cupp)) if cupp else ""
            sub = "+".join(sorted(target["sub"]))
            function_string = "+".join(sorted(target["function"]))
            target_for_family = set(f.split(":",1)[-1] for f in target[tag] if f.split(":")[0].split("_")[0] == fam)
            predict_pool = set(); target_pool = set(); better = False
            precision = 0; sensitivity = 0
            
            ####################################
            #Handle special metadata separators ("&","+","-") for predict:
            if fam in predicted[tag]:
                for i in predicted[tag][fam]:
                    for f in i.split("&"):
                        predict_pool.add(f)
                predict_pool_two = set()
                for i in predict_pool:
                    for f in i.split("-"):
                        predict_pool_two.add(f)
                predict_pool = set()
                for i in predict_pool_two:
                    for f in i.split("+"):
                        predict_pool.add(f)
                predict_pool = set((k.split(":",1)[-1]) for k in predict_pool if k.split(":",1)[-1] != "Unknown")
                
            ####################################
            #Handle special metadata separators ("&","+","-") for target:
            if target_for_family:
                for i in target_for_family:
                    for f in i.split("&"):
                        target_pool.add(f.split(":",1)[-1])
                target_pool_two = set()
                for i in target_pool:
                    for f in i.split("-"):
                        target_pool_two.add(f.split(":",1)[-1])
                target_pool = set()
                for i in target_pool_two:
                    for f in i.split("+"):
                        target_pool.add(f.split(":",1)[-1])
                target_pool = set((k.split(":",1)[-1]) for k in target_pool)
          
            ####################################
            #Do not count if there is no relevant target:
            if not target_pool:
                continue
                
            ####################################
            #If the function attempted to be predicted is not in the dataset, continue if still unknown:
            known_unknown = False
            if tag in lost_meta:
                completely_lost = True
                for f in target_pool:
                    if lost_meta[tag][f][1]:
                        completely_lost = False
                if completely_lost and not predict_pool:
                    known_unknown = True
                    
            ####################################
            #Determine if there is a dual function in the target or in the predicted:
            if target_pool and predict_pool:
                no_overlap = True
                dual_predict = set((k) for k in predicted[tag][fam] if "&" in k)
                dual_target =  set((k) for k in target_for_family if "&" in k)
                lookup_dual = {}
                for dual_function in dual_predict:
                    for f in dual_function.split(":",1)[-1].split("&"):
                        if f in lookup_dual and lookup_dual[f] != dual_function:
                            no_overlap = False
                        lookup_dual[f] = dual_function
                for dual_function in dual_target:
                    for f in dual_function.split(":",1)[-1].split("&"):
                        if f in lookup_dual and lookup_dual[f] != dual_function:
                            no_overlap = False
                        lookup_dual[f] = dual_function
                    
                ####################################
                #Check if the individual functions are only found within one dual function and use the dual function instead:
                if no_overlap:
                    predict_pool = set(lookup_dual[f] if f in lookup_dual else f for f in predict_pool)
                    target_pool  = set(lookup_dual[f] if f in lookup_dual else f for f in target_pool)
            
            ####################################
            #Determine sensitivity and precision:
            shared = predict_pool & target_pool
            add_predicted = predict_pool - (predict_pool & target_pool)
            wrong = len(predict_pool) + (len(predictions_outside_target_families)/len(holding_families) if holding_families else 0)
            precision = (len(predict_pool)-len(add_predicted))/wrong if wrong else 1
            sensitivity = len(shared)/len(target_pool) if target_pool else 1
            sensitivity = 1 if known_unknown else sensitivity
                
            ####################################
            #Save the results for later:
            if fam not in fam_overview:
                fam_overview[fam] = {"function":[0,0,0],"sub":[0,0,0],"cup":[0,0,0],"family":[0,0,0,{}],"agg":[0,0,0]}
            fam_overview[fam][tag][0] += precision
            fam_overview[fam][tag][1] += sensitivity
            fam_overview[fam][tag][2] += 1
            
            ####################################
            #Remember the additional families:
            if tag == "family":
                for f in predictions_outside_target_families:
                    if f not in fam_overview[fam][tag][3]:
                        fam_overview[fam][tag][3][f] = 0
                    fam_overview[fam][tag][3][f] += 1

            ####################################
            #Save the results for later of cupp:
            if cupp:
                if cupp not in cupp_overview:
                    cupp_overview[cupp] = {"function":[0,0,0],"sub":[0,0,0],"cup":[0,0,0],"family":[0,0,0],"agg":[0,0,0]}
                cupp_overview[cupp][tag][0] += precision
                cupp_overview[cupp][tag][1] += sensitivity
                cupp_overview[cupp][tag][2] += 1
            
            ####################################
            #Save the results for later of function:
            if function_string:
                if function_string not in function_overview:
                    function_overview[function_string] = {"function":[0,0,0],"sub":[0,0,0],"cup":[0,0,0],"family":[0,0,0],"agg":[0,0,0]}
                function_overview[function_string][tag][0] += precision
                function_overview[function_string][tag][1] += sensitivity
                function_overview[function_string][tag][2] += 1            
           
            ####################################
            #Save the results for later of subfamily:
            if sub:
                if sub not in sub_overview:
                    sub_overview[sub] = {"function":[0,0,0],"sub":[0,0,0],"cup":[0,0,0],"family":[0,0,0],"agg":[0,0,0]}
                sub_overview[sub][tag][0] += precision
                sub_overview[sub][tag][1] += sensitivity
                sub_overview[sub][tag][2] += 1
           
            ####################################
            #Output the precision and sensitivity for selected accession numbers:
            if string in focus:
                print (precision,sensitivity,target_pool,"-->",predict_pool)
   
            ####################################
            #Output the precision and sensitivity for selected accession numbers:
            if precision != 1 and tag == "family":
                bad_flag = string+" "+"+".join(target["Family"])+(" %s" % round(precision,2))+" "+str(round(sensitivity,2))+" "+"+".join(predictions_outside_target_families)+" --> "+"+".join(predicted["family"])
    return bad_flag
   
def determine_relatedness(predictions,all_meta,CUPPpool_path,dendrogram_path,rel_cut=0.7):
    '''
    Determine the frequency of two CUPP groups being located in the same positions when they are found together in the same proteins
    '''
    
    ####################################
    #Prepare the prediction overview:
    connectivity = {}
    for acc,v in predictions.items():
        for one in v["rel"]:
            if one not in connectivity:
                connectivity[one] = {}
            for two in v["rel"][one]:
                if two not in connectivity[one]:
                    connectivity[one][two] = []
                connectivity[one][two] += v["rel"][one][two]

    
    ####################################
    #Identify potential overlap in the position of the individual CUPPgroups:
    mono_cup = set((k) for k,v in connectivity.items() if len([l for c,l in v.items() if sum(l)]) == 1)
    shar_cup = set((k) for k,v in connectivity.items() if v)
    cut = 0; assignments = {}

    if len(mono_cup) != len(shar_cup) and shar_cup:
        order = dict((no,sorted(connectivity)[no]) for no in range(len(connectivity)))
        CUPPlist = [v for k,v in sorted(order.items(),key=lambda x:x[0])]
        pos = dict((v,k) for k,v in order.items())
        matrix = np.ones((len(order),len(order)))

        ####################################
        #Scoring of connectivity of CUPPgroups:
        for cup_one in connectivity:
            for cup_two,ratio in connectivity[cup_one].items():
                matrix[pos[cup_one]][pos[cup_two]] -= sum(ratio)/len(ratio)
                
        ####################################
        #Clustering of CUPP connectivity:
        matrix = squareform(matrix)
        links = linkage(matrix,method="complete")
        
        ####################################
        #Locate best friend with an EC number:
        while cut < rel_cut:
            assign = [f for f in cut_tree(links,height=cut)]
            assign = dict(zip(CUPPlist,assign))
            for cup,i in assign.items():
                if cup not in assignments:
                    assignments[cup] = {cup:cut} 
                for other,j in assign.items():
                    if i == j and other not in assignments[cup]:
                        assignments[cup][other] = round(cut,4) 
            cut += 0.01
            
    if cut:
        labels_dict = dict((v,"") for k,v in order.items())
        for k,cup in order.items():
            labels_dict[cup] += cup+"_"
            labels_dict[cup] += "-".join(["%s#%s" % (k.split(":",1)[-1],v) for k,v in sorted(all_meta[cup]["Subfam"].items())])+"_"
            
            function = pool_pool(all_meta[cup]["raw_function"])
            if not function and assignments:
            
                ####################################
                #In case multiple friends are found with the same relatedness:
                friend_pool = {}
                for friend,related in assignments[cup].items():
                    for f,s in pool_pool(all_meta[friend]["raw_function"]).items():
                        if related <= rel_cut:
                            if f not in friend_pool:
                                friend_pool[f] = 0
                            friend_pool[f] += s
                function = pool_pool(friend_pool)
                function = dict((k,0) for k in function)
            labels_dict[cup] += "-".join(function)
            
        dendro(links,[labels_dict[l] for n,l in sorted(order.items(),key=lambda x:x[0])],algo="complete",cc=rel_cut,name=dendrogram_path,fil=dendrogram_path,show=False,note=";".join(sorted(mono_cup)))
        
    ####################################
    #Adding the relatedness to the powder files prior to recompiling:
    print ("#"*60+"\n### The connectivity of CUPPgroups are being determined",flush=True, end="")
    with open(CUPPpool_path) as json_data:
        powder = json.load(json_data)
    for s,v in powder.items():
        powder[s].update({"related":assignments})
    handle_out = open(CUPPpool_path,'w')
    json.dump(powder,handle_out,indent=3)
    handle_out.close()
    print (" DONE!")
    
    ####################################
    #Output any strong relationships between groups:
    for k,v in sorted(assignments.items(),key=lambda x:family_sort(x[0])):
        if k not in v:
            print ("### %s <-- %s" % (k,"-".join(["%s#%s" % (i,round(j,3)) for i,j in sorted(v.items(),key=lambda x:x[1],reverse=True)])))

def association(seq,residues,domains,predicted,n_mer,ambiguous,final_prediction):
    '''
    Store the relation between CUPP groups during prediction and used it to determine associations between CUPP groups
    '''

    ####################################
    #Keep relevant domains:
    residues = dict((cupp,v) for cupp,v in residues.items() if cupp in domains)
    
    ####################################
    #Inspect the ratios between X and non-X positions of peptides:
    for cupp,v in residues.items():
        liste = [0 for f in range(len(seq)+ambiguous*2)]
        for n in range(len(liste)):
            if v[1][n] and final_prediction[cupp][n]:
                liste[n] = v[0][n]/v[1][n]-(1-ambiguous/n_mer)
                if 0 < liste[n] and 1 < v[1][n]:
                    liste[n] = liste[n]/(ambiguous/n_mer)
                else:
                    liste[n] = 0
                    
        ####################################
        #Normalize the X-ratios:            
        normalized_liste = [0 for f in range(len(seq)+ambiguous*2)]
        for n in range(len(liste)):
            frame = [f for f in liste[0 if n-4 < 0 else n-4:-1 if len(liste) < n+4 else n+4] if 0.1 < f]
            average = sum(frame)/len(frame) if len(frame) else 0
            normalized_liste[n] = -1*liste[n]/average if average else 0
        residues[cupp] = normalized_liste

    ####################################
    #Filter to keep only CUPPgroup of the target family:
    connectivity = {}
    for fam,v in predicted["cupp"].items():
        for cup_one,score in v.items():
            if cup_one not in connectivity:
                connectivity[cup_one] = {}
            pos_one = set(no for no in range(len(residues[cup_one])) if residues[cup_one][no] < 0)
            for cup_two in v:
                pos_two = set(no for no in range(len(residues[cup_two])) if residues[cup_two][no] < 0)
                if max(len(pos_one),len(pos_two)) == 0:
                    continue
                if cup_two not in connectivity[cup_one]:
                    connectivity[cup_one][cup_two] = []
                shared = pos_one & pos_two
                ratio_score = len(shared)/(len(pos_one - shared)+len(pos_two - shared)+len(shared))
                connectivity[cup_one][cup_two].append(ratio_score)
    
    return connectivity

def orf(acc,seq,codontable,minimum=60*3):
    """
    Locate open reading frames within a dna fragment and return a collection of the amino acid sequences
    """
    ####################################
    #Compile regular expression of start and stop codons:
    init = re.compile("ATG")
    stop = re.compile("TAG|TAA|TGA")
    
    ####################################
    #Identify all stop codons in config:
    stop  = {f.end() for f in stop.finditer(seq)}
    collection = {}
    
    ####################################
    #Identify all start codons in contig:
    last_end = set()
    for start in init.finditer(seq):
    
        ####################################
        #Identify start and stop position of each protein:
        start = start.start()
        end = start
        while end < len(seq) and end not in stop:
            end += 3
            end = len(seq) if len(seq) < end else end
            
        ####################################
        #Extract the relevant DNA of config:
        cds = seq[start:end]
        if minimum <= len(cds):
        
            ####################################
            #Keep only the longer variants of each proteins (sharing stop codon):
            if end not in last_end:
                collection["%s:%s..%s" % (acc,start,end)] = "".join([codontable[cds[n:n+3]] if cds[n:n+3] in codontable else "X" if len(cds[n:n+3]) == 3 else "" for n in range(0,len(cds),3)])
            last_end.add(end)
            
    ####################################
    #Identify all stop codons in config in the reverse direction:
    init = re.compile("CAT")
    stop = re.compile("CTA|TTA|TCA")
    stop = {f.start() for f in stop.finditer(seq)}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    
    ####################################
    #Identify all start codons in contig in the reverse direction:
    last_end = set()
    
    ####################################
    #Iterate over each start codon from the end of the config:
    for start in reversed([start.end() for start in init.finditer(seq)]): 
        end = start
        while 0 < end and end not in stop:
            end -= 3
            end = 0 if end < 0 else end
            
        ####################################
        #Extract the relevant DNA of config, reversed and complemented:
        cds = "".join([complement[b] if b in complement else "N" for b in seq[end:start][::-1]])
        
        ####################################
        #Keep only the longer variants of each proteins (sharing stop codon):
        if minimum <= len(cds):
            if end not in last_end:
                collection["%s:%s..%s" % (acc,start,end)] = "".join([codontable[cds[n:n+3]] if cds[n:n+3] in codontable else "X" if len(cds[n:n+3]) == 3 else "" for n in range(0,len(cds),3)])
            last_end.add(end)
        
    return (collection)
       
def determine_function(domains,all_meta,precision_ratio,limit_domains,occurence=False,beta={}):
    '''
    Determine the family,subfamily and EC function of query protein.
    '''
    
    ####################################
    #Domains contain the predicted CUPPgroups:
    predicted = {"best_function":{},"raw_function":{},"Family":{},"cupp":{},"Subfam":{},"agg":{},"additional":{}}
    
    ####################################
    #Prediction of family:
    for cupp,score in domains.items():
        fam = cupp.split(":")[0]
        if fam not in predicted["Family"]:
            predicted["Family"][fam] = {}
            predicted["Family"][fam][fam] = 0
        predicted["Family"][fam][fam] += score

    ####################################
    #Sorting the predicted CUPPgroups:
    for cupp,score in domains.items():
        fam = cupp.split(":")[0]
        if fam not in predicted["cupp"]:
            predicted["cupp"][fam] = {}
        predicted["cupp"][fam][cupp] = score

    ####################################
    #Prediction of subfamily:
    for fam in set(f.split(":")[0] for f in domains):
        if fam not in predicted:
            predicted["Subfam"][fam] = {}
        for cupp,score in predicted["cupp"][fam].items():
        
            ####################################
            #zero groups can not have an EC:
            if cupp not in all_meta:
                continue
        
            sub_string = "+".join([sub for sub,occur in sorted(all_meta[cupp]["Subfam"].items()) if 3 <= occur]) if "Subfam" in all_meta[cupp] else ""
            if sub_string not in predicted["Subfam"][fam]:
                predicted["Subfam"][fam][sub_string] = score
            else:
                predicted["Subfam"][fam][sub_string] = max(score,predicted["Subfam"][fam][sub_string])
            
        ####################################
        #In case the best scoring subfamily is unknown (emtry string), no subfamily will be predicted:
        predicted["Subfam"][fam] = dict((sub,score) for sub,score in predicted["Subfam"][fam].items() if max(predicted["Subfam"][fam].values()) == score and sub)
        
    ####################################
    #Prediction of the function:
    for cupp,score in sorted(limit_domains.items(),key=lambda x:x[1],reverse=True):
        fam = cupp.split(":")[0]
        
        ####################################
        #Remove EC extrapolation:
        if "no_fam_ec" in beta and not all_meta[cupp]["raw_function"]:
            all_meta[cupp]["Function"] = {}

        ####################################
        #Consider all EC function of a CUPPgroup equally:
        target = all_meta[cupp]["Function"]
        for f in target:
            if fam not in predicted["raw_function"]:
                predicted["raw_function"][fam] = {}
            if f not in predicted["raw_function"][fam]:
                predicted["raw_function"][fam][f] = 0
            predicted["raw_function"][fam][f] += score/len(target)
            
        ####################################
        #Consider only the must occuring EC function of the CUPPgroup:
        highest_occur = pool_pool(all_meta[cupp]["raw_function"],reduce=False) if "raw_function" in all_meta[cupp] else {}
        highest_occur = pool_pool(all_meta[cupp]["Function"],reduce=False) if not highest_occur else highest_occur
        highest_occur = dict((k,v) for k,v in highest_occur.items() if v == max(highest_occur.values()))
        for f in highest_occur:
            if fam not in predicted["best_function"]:
                predicted["best_function"][fam] = {}
            if f not in predicted["best_function"][fam]:
                predicted["best_function"][fam][f] = 0
            predicted["best_function"][fam][f] += score/len(highest_occur)
            
        ####################################
        #Only consider the must occuring EC function:
        if occurence:
            predicted["function"] = predicted["best_function"]
            
        ####################################
        #Function is unknown if the group with the highest score in of Unknown function:
        if not all_meta[cupp]["Function"]:
            for item in ["raw_function","best_function"]:
                if fam not in predicted[item]:
                    predicted[item][fam] = {}
                if "%s:Unknown" % fam not in predicted[item][fam]:
                    predicted[item][fam]["%s:Unknown" % fam] = 0
                predicted[item][fam]["%s:Unknown" % fam] += score
    
    ####################################
    #Consider unknown function equal to any other EC function:
    for fam in predicted["raw_function"]:
        for item in ["raw_function","best_function"]:
            if "%s:Unknown" % fam in predicted[item][fam] and predicted[item][fam]["%s:Unknown" % fam] == max(predicted[item][fam].values()):
                predicted[item][fam] = dict((k,v) for k,v in predicted[item][fam].items() if k == "%s:Unknown" % fam)
            else:
                predicted[item][fam] = dict((k,v) for k,v in predicted[item][fam].items() if k != "%s:Unknown" % fam)

    '''    
    ####################################
    #Prediction of CUPPagg:
    for group,score in domains.items():
        fam,number = group.split(":")
        agg = fam+":"+number.split(".")[0]
        if agg:
            if fam not in predicted["agg"]:
                predicted["agg"][fam] = {}
            if agg not in predicted["agg"][fam]:
                predicted["agg"][fam][agg] = 0
            predicted["agg"][fam][agg] += score
    '''
    
    ####################################
    #Keep only CUPPgroup with 50% or greater than the best hit of the family:
    for tag in predicted:
        if tag != "cupp":
            for fam,index in predicted[tag].items():
                predicted[tag][fam] = dict((k,v) for k,v in index.items() if max(index.values())*precision_ratio <= v)
            
    return (predicted)
                
def rank_domains(cupp,short_list,graph,freq_sum,domain_min,cup_cut,evidence,line_cut):
    '''
    Determine order of domains during filtering.
    '''
    
    if domain_min <= len([a for a in graph[cupp] if line_cut < a]):
    
        if freq_sum[cupp]*evidence < short_list[cupp]:
            return short_list[cupp]/freq_sum[cupp]
    return 0

def predict(seq,opts):
    '''
    Prediction Families, Subfamilies, Functions and CUPP groups within each protein
    '''
    
    ####################################
    #Define the prediction variables:
    cup_cut = opts.minimum_cup/opts.evidence; freq_sum = {}
    final_prediction, domains, short_list, amplify = {},{},{},{}
    graph = {}; residues = {}; opts.minimum_quick = opts.domain_min-2*(opts.n_mer-opts.ambiguous)
    predicted = dict((tag,{}) for tag in ["best_function","cupp","raw_function","Family","Subfam","rel","additional"])
    
    ####################################
    #Prefiltering of domains:
    peptides,tandem = obtain_peptides({seq:"NEW1"},n_mer=opts.n_mer,ambiguous=opts.ambiguous,limit=opts.all_peptides,single=True)

    ####################################
    #Loop over each peptide:
    for pep in peptides:
    
        ####################################
        #Loop over each of the group of a single peptide:
        for cupp,score in opts.all_peptides[pep].items():
        
            ####################################
            #Add the weights of peptides for the individual CUPP groups
            if cupp not in short_list:
                short_list[cupp] = 0
            short_list[cupp] += score
            
    ####################################
    #Determine the current sum of weights:
    freq_sum = dict((cupp,(opts.meta[cupp]["Freq_sum"]["sum"] if cup_cut < opts.meta[cupp]["Freq_sum"]["sum"] else cup_cut)*opts.evidence) for cupp in short_list)
                
    ####################################
    #Prefilter the domain hits to increase speed:
    short_list = dict((cupp,v) for cupp,v in short_list.items() if freq_sum[cupp] <= v and cupp.split(":")[0] not in opts.exclude_family and cupp not in opts.exclude_family and (not opts.keep_only or (cupp.split(":")[0] in opts.keep_only or cupp in opts.keep_only)))
    
    ####################################
    #Filter the predictions based on covered positions:
    if short_list:
        
        ####################################
        #Apply filtering of domains based on positions and overlap:
        if opts.type != "quick":
            
            ####################################
            #Prefiltering of domains:
            for pep,pos in peptides.items():
                
                ####################################
                #Loop over each of the group of a single peptide:
                for cupp,score in opts.all_peptides[pep].items():
                    
                    ####################################
                    #Jump prefiltered domains:
                    if cupp not in short_list:
                        continue
                        
                    ####################################
                    #Make a diagram of the domains:
                    if cupp not in graph:
                        graph[cupp] = np.zeros(len(seq)+opts.ambiguous*2)
                    for n in range(opts.n_mer):
                        if pep[n] != "X":
                            graph[cupp][pos+n] += score
                            
                    ####################################
                    #Determine relatedness:
                    if opts.relatedness:
                        if cupp not in residues:
                            residues[cupp] = [np.zeros(len(seq)+opts.ambiguous*2),np.zeros(len(seq)+opts.ambiguous*2)]
                        for n in range(opts.n_mer):
                            if pep[n] != "X":
                                residues[cupp][0][pos+n] += score
                            residues[cupp][1][pos+n] += score
                    
                    ####################################
                    #If one peptide in present multiple places, it may count in range:
                    if pep in tandem:
                    
                        for p in tandem[pep]:
                        
                            ####################################
                            #Make a diagram of the domains:
                            if cupp not in graph:
                                graph[cupp] = np.zeros(len(seq)+opts.ambiguous*2)
                            for n in range(opts.n_mer):
                                if pep[n] != "X":
                                    graph[cupp][p+n] += score
                                    
                            ####################################
                            #Determine relatedness:
                            if opts.relatedness:
                                if cupp not in residues:
                                    residues[cupp] = [np.zeros(len(seq)+opts.ambiguous*2),np.zeros(len(seq)+opts.ambiguous*2)]
                                for n in range(opts.n_mer):
                                    if pep[n] != "X":
                                        residues[cupp][0][p+n] += score
                                    residues[cupp][1][p+n] += score
                            
            ####################################
            #Exclude too short domains:
            graph = dict((cupp,v) for cupp,v in graph.items() if opts.domain_min <= len([a for a in v if opts.line_cut < a]))
            
            ####################################
            #Filter the domain based on the covered positions:
            occupied = {}; freq_domains = {}
            for cupp,profile in sorted(graph.items(),key=lambda x:short_list[x[0]]/freq_sum[x[0]],reverse=True):
            
                ####################################
                #Test how many positions are greatest:
                different_sub = [0,0]
                for pos in range(len(profile)):
                    if pos not in occupied or occupied[pos] < profile[pos]:
                        different_sub[0] += profile[pos]
                        if pos in occupied and occupied[pos] < profile[pos]:
                            different_sub[0] -= occupied[pos]
                    different_sub[1] += profile[pos]

                ####################################
                #The region of the new domain may not already be occupied by another domain:
                if 1-opts.overlap < different_sub[0]/different_sub[1]:
                
                    ####################################
                    #Remember the highest score of the occupied positions:
                    for pos in range(len(profile)):
                        if pos not in occupied or occupied[pos] < profile[pos]:
                            occupied[pos] = profile[pos]
                    
                    ####################################
                    final_prediction[cupp] = [f if opts.line_cut < f else 0 for f in profile]
                    domains[cupp] = short_list[cupp]/freq_sum[cupp]*sum([1 if f else 0 for f in final_prediction[cupp]])/opts.meta[cupp]["Freq_sum"]["Average"]
                    freq_domains[cupp] = short_list[cupp]/freq_sum[cupp]
                    
        ####################################
        #Filter the domain based on the covered positions:
        #freq_domains = {}; 
        last = 0; additional = {}
        if "mix_type" in opts.beta_options or opts.type == "quick":
            for cupp,s in sorted(short_list.items(),key=lambda x:short_list[x[0]]/freq_sum[x[0]],reverse=True):
            
                ####################################
                #Filter based on the initial score and not domain length:
                last = s if last < s else last

                ####################################
                #The region of the new domain may not already be occupied by another domain:
                if 1-opts.overlap < s/last:
                
                    #Reduce number of false positives in the quick prediction:
                    if opts.minimum_quick < len(set((n) for p,n in peptides.items() if cupp in opts.all_peptides[p])):
                        additional[cupp] = short_list[cupp]/freq_sum[cupp]
                    
        #Advanced output format:
        limit_domains = domains; overlapping = {}
        if opts.type in ["both","raw"] and not opts.relatedness:
        
            ####################################
            #Identify range of the domain:
            included  = {}; ranges = {}
            for cupp in domains:
                included[cupp] = {n+1-opts.ambiguous for n in range(len(final_prediction[cupp])) if final_prediction[cupp][n]}
                ranges[cupp] = []
                
                start = -1; end = -1; 
                for n in sorted(included[cupp]):
                
                    ####################################
                    #Continue if the position is in the domain:
                    if n < end:
                        continue
                        
                    ####################################
                    #The first start position:
                    if start == -1:
                        start = n
                
                    ####################################
                    #Is the next gap included in the domain:
                    flag = False
                    for i in range(1,opts.meta[cupp]["Freq_sum"]["Average"]):
                        if n+i in included[cupp]:
                            end = n+i
                            flag = True
                            
                    ####################################
                    #To the gap was too large to be a single domain:
                    if not flag:
                        ranges[cupp].append("%s..%s" % (start,end))
                        start = -1
                        
                ####################################
                #Evaluate if each range can stand alone:
                if not "keep_ranges" in opts.beta_options and 1 < len(ranges[cupp]):
                    new_ranges = []
                    for r in ranges[cupp]:
                        s,e = r.split("..")[0],r.split("..")[1]
                        if opts.domain_min <= len(included[cupp] & set(range(int(s),int(e)))):
                            new_ranges.append("%s..%s" % (s,e))
                    ranges[cupp] = new_ranges
                    
            ####################################
            #Check if ranges overlap:
            ban_range = {}
            if "complex" not in opts.beta_options:
                
                ####################################
                #Handle multiple domains in the same range:
                if 1 < len(ranges):
            
                    ####################################
                    #Pool all ranges:
                    relevant = {}; trace = {}
                    for cupp,ran in sorted(ranges.items()):
                        for r in sorted(ran):
                            trace[len(relevant)] = cupp
                            relevant[len(relevant)] = [int(r.split("..")[0]),int(r.split("..")[1]),cupp,r]
                        
                    ####################################
                    #Delete CUPP groups of overlapping ranges:
                    for i,ran_i in relevant.items():
                        for j,ran_j in relevant.items():
                            if i < j:
                                set_i  = set(range(ran_i[0],ran_i[1]))
                                set_j  = set(range(ran_j[0],ran_j[1]))
                                common = len(set_i & set_j)
                                uni_i  = len(set_i - set_j)
                                uni_j  = len(set_j - set_i)
                                if uni_i/(common+uni_i) < opts.overlap or uni_j/(common+uni_j) < opts.overlap:
                                    #Make the 50% cupp group filtering here instead of before range filtering
                                    if domains[trace[i]] < domains[trace[j]]*opts.precision_ratio:
                                        if relevant[i][2] not in ban_range:
                                            ban_range[relevant[i][2]] = []
                                        ban_range[relevant[i][2]].append(relevant[i][3])
                                    elif domains[trace[j]] < domains[trace[i]]*opts.precision_ratio:
                                        if relevant[j][2] not in ban_range:
                                            ban_range[relevant[j][2]] = []
                                        ban_range[relevant[j][2]].append(relevant[j][3])
                                    else:
                                        overlapping[relevant[j][2]] = relevant[i][2]
                                        overlapping[relevant[i][2]] = relevant[j][2]
                    
                    ####################################
                    #Keep only non-overlapping cupp groups:
                    ranges = dict((cupp,[r for r in v if cupp not in ban_range or r not in ban_range[cupp]]) for cupp,v in ranges.items())
                    
                ####################################
                #Family prediction is fine with freqs around 1, however the correct CUPP group, sub and EC may be difficult thus require a higher freq:
                limit_domains = dict((cupp,v) for cupp,v in domains.items() if opts.cupp_minimum_score <= short_list[cupp]/freq_sum[cupp])
                
            ####################################
            #Determine the function and subfamily of protein:
            predicted = determine_function(domains,opts.meta,opts.precision_ratio,limit_domains,occurence=opts.occurence,beta=opts.beta_options)
            
            ####################################
            #Keep only CUPP groups with ranges, if any have a range:
            if "complex" not in opts.beta_options:
                predicted["cupp"] = dict((fam,dict((cupp,s) for cupp,s in v.items() if not sum([1 for c in v if c in ranges and ranges[c]]) or cupp in ranges and ranges[cupp])) for fam,v in predicted["cupp"].items())
            
            ####################################
            temp = {"cupp":{},"additional":{},"raw_function":predicted["raw_function"],"best_function":predicted["best_function"],"Subfam":predicted["Subfam"],"Family":predicted["Family"]}
            for fam,index in predicted["cupp"].items():
                temp["cupp"][fam] = {}
                for cupp,s in sorted(index.items(),key=lambda x: domains[x[0]],reverse=True):
                
                    ####################################
                    #If a single CUPP group cannot be determined, the CUPP string will be modified to "fam:0.1":
                    new_cupp = cupp if cupp in limit_domains else fam+":0.1"
                    new_cupp = new_cupp if cupp not in overlapping or 1 == sum([1 for f in overlapping[cupp] if f.split(":")[0] == fam]) else fam+":0.2"
                    
                    ####################################
                    #Recall relevant information from prediction for output:
                    temp["cupp"][fam][new_cupp] = s
                    temp["additional"][new_cupp] = [str(round(freq_domains[cupp],1))]
                    temp["additional"][new_cupp].append(str(round(domains[cupp],2)))
                    temp["additional"][new_cupp].append(str(int([f for f in ndimage.center_of_mass(graph[cupp])][0])))
                    temp["additional"][new_cupp].append("&".join(ranges[cupp] if cupp in ranges else []))
                    temp["additional"][new_cupp].append(str(len(included[cupp])))
                    
                    ####################################
                    #If a single CUPP group cannot be determined, the CUPP string will be modified and kept once:
                    if new_cupp != cupp:
                        break
                        
            ####################################
            #Include the weaker prediction of quick, not considered for the strick prediction:
            if "mix_type" in opts.beta_options:
                taken = set()
                if not temp["cupp"]:
                    for cupp,s in sorted(additional.items(),key=lambda x: x[1],reverse=True):
                        fam = cupp.split(":")[0]
                        
                        if fam in taken:
                            continue
                            
                        new_cupp = fam+":0.0" if "debug" not in opts.beta_options else cupp
                        temp["cupp"][fam] = {new_cupp:s}
                        temp["additional"][new_cupp] = [str(round(s,1))]
                        temp["additional"][new_cupp].append(str(round(s,1)))
                        temp["additional"][new_cupp].append("N/A")
                        temp["additional"][new_cupp].append("N/A")
                        temp["additional"][new_cupp].append("N/A")
                        temp["Family"][fam] = {fam:s}
                        taken.add(fam)
                    
            ####################################
            #Replace the prediction with more reduced complexity:
            predicted = temp
            
        else:
            
            ####################################
            #Reduce the information level of the quick prediction to only contain family info:
            if not opts.relatedness:
                predicted = determine_function(additional,opts.meta,opts.precision_ratio,limit_domains,occurence=opts.occurence,beta=opts.beta_options)
                temp = {"cupp":{},"additional":{},"raw_function":{},"best_function":{},"Subfam":predicted["Subfam"],"Family":predicted["Family"]}
                for fam,index in predicted["cupp"].items():
                    for cupp,s in sorted(index.items(),key=lambda x: additional[x[0]],reverse=True):
                        
                        temp["cupp"][fam] = {}
                        temp["cupp"][fam][fam+":0.0"] = s
                        temp["additional"][fam+":0.0"] = s
                        break                        
                    
                predicted = temp
            else:
                predicted = determine_function(domains,opts.meta,opts.precision_ratio,limit_domains,occurence=opts.occurence,beta=opts.beta_options)
            
        ####################################
        #Apply selected output type:
        if opts.type == "quick":
            predicted["additional"]["seq"] = seq
        elif opts.type in ["both","fasta"]:
            if "float" in opts.beta_options:
                predicted["additional"]["seq"] = "\t".join([str(round(max([p[n-opts.ambiguous] for cupp,p in final_prediction.items() if cupp in domains]),3) if final_prediction else 0) for n in range(len(seq))])
            else:
                predicted["additional"]["seq"] = "".join([seq[n].upper() if sum([p[n-opts.ambiguous] for cupp,p in final_prediction.items() if cupp in domains]) else seq[n].lower() for n in range(len(seq))])
                
        if opts.relatedness:
            predicted["rel"] = association(seq,residues,domains,predicted,opts.n_mer,opts.ambiguous,final_prediction)
            
    elif "keep" in opts.beta_options:
        predicted["additional"]["seq"] = seq.lower()
            
    return predicted
    
def out(handle_fas,handle_raw,acc,q,seq,beta,meta,sep="|"):
    '''
    Format the output of CUPP prediction
    '''

    ####################################
    #Output the FASTA file for prediction:
    hit = bool(q["Family"])
    if handle_fas:
        string = ">%s" % (acc+sep+"&"+sep)
        if hit:
            string += "+".join(["-".join(sorted(v)) for k,v in sorted(q["raw_function"].items()) if v])+sep
            string += "+".join(["-".join(sorted(v)) for k,v in sorted(q["best_function"].items()) if v])+sep
            string += "+".join(["-".join([("%s(%s,%s)" % (i,q["additional"][i][0] if i in q["additional"] else "",q["additional"][i][3] if i in q["additional"] else "")) if handle_raw else i for i,j in sorted(v.items(),key=lambda x: x[1],reverse=True)]) for k,v in sorted(q["cupp"].items(),key=lambda x: family_sort(x[0]))])+sep
            string += "+".join(["-".join(sorted(v)) for k,v in sorted(q["Subfam"].items()) if v]).replace(":","_")+"\n"
            string += q["additional"]["seq"].replace("_","x")+"\n"
            handle_fas.write(string)
        elif "keep" in beta and seq:
            string += "|"*3+"\n"+q["additional"]["seq"].replace("_","x").lower()+"\n"
            handle_fas.write(string)
        
    ####################################
    #Output raw file for prediction:
    if handle_raw and hit:
        for fam,v in q["cupp"].items():
            for cupp,content in sorted(v.items(),key=lambda x: x[1],reverse=True):
                string = "%s\t%s\t%s\t%s" % (acc,fam,cupp,"\t".join(q["additional"][cupp]))
                if cupp.split(":")[1][0] == "0":
                    ec1 = sorted(q["raw_function"][fam]) if fam in q["raw_function"] else []
                    ec2 = sorted(q["best_function"][fam]) if fam in q["best_function"] else []
                    sub = sorted(q["Subfam"][fam]) if fam in q["Subfam"] else []
                else:
                    ec1 = [f for f in sorted(meta[cupp]["Function"]) if fam in q["raw_function"] and f in q["raw_function"][fam]]
                    ec2 = [f for f,s in sorted(meta[cupp]["Function"].items()) if fam in q["best_function"] and f in q["best_function"][fam]]
                    sub = [f.replace(":","_") for f,s in sorted(meta[cupp]["Subfam"].items()) if fam in q["Subfam"] and f in q["Subfam"][fam]]
                    
                string += "\t" + ("-".join(ec1))
                string += "\t" + ("-".join(ec2))
                string += "\t" + ("-".join(sub)) + "\n"
                    
                handle_raw.write(string)
                    
    return hit
    
if __name__ == "__main__":

    '''
    Run the CUPPprediction.py script alone with the using CUPPclustering.py script
    '''

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    
    ####################################
    #Load user arguments:
    args = arg_parser_predict(sys.argv[1:])
    
    ####################################
    #Start with one query of a folder:
    if not args.query and args.dir_query:
        args.query = sorted(args.dir_query)[0]
        print ("--> %s" % args.query)
        
    ####################################
    #Check if any of the FASTA files are valid:
    if (args.query and os.path.exists(args.query)) or args.live:
    
        ####################################
        #Load the conserved peptides:
        all_peptides = {}
        if os.path.exists(args.compiled_json):
            print ("### The precompiled CUPPlibrary is being loaded: %s" % args.setting,flush=True,end="")
            with open(args.compiled_json) as json_data:
                powder_merged = json.load(json_data)
            print (" DONE!")
            all_peptides = powder_merged["peptides"]; CUPPmeta = powder_merged["meta"]; meta_cat = powder_merged["meta_categories"]
        if not all_peptides or not os.path.exists(args.compiled_json):
            sys.exit("!!! The CUPP library does not exists or is emtry: %s" % args.compiled_json)
        CUPPmeta = create_table(CUPPmeta,meta_cat=meta_cat,silent=True if not args.table else False); already_processed ={}
        
        ####################################
        #To be included during recompiling instead: 
        for cupp,v in CUPPmeta.items():
            if "Average" not in CUPPmeta[cupp]["Freq_sum"]:
                CUPPmeta[cupp]["Freq_sum"]["Average"] = 201
            else:
                CUPPmeta[cupp]["Freq_sum"]["Average"] = 201 if 200 < v["Freq_sum"]["Average"] else v["Freq_sum"]["Average"]+1 if 40 < v["Freq_sum"]["Average"] else 41

        ####################################
        #Keep the library for another query FASTA file:
        next = True; first = True
        while next:
            
            ####################################
            #Start the live server, without a query initially:
            if args.query and os.path.exists(args.query) and (not args.dir_query or args.query not in already_processed or os.path.getmtime(args.query) != already_processed[args.query]):
                
                ####################################
                #Determine if some family or CUPP groups could be left out:
                if "Experimentals" not in powder_merged["stats"]:
                    args.exclude_family.update({"GH0","CE0","GT0","AA0","PL0","PL28","GH145","GH146","GH147","GH148","GH149","GH150","GH151","GH152","GH153","AA14","AA15","GT105"})
                args.exclude_family = args.exclude_family if "experimental_families" not in args.beta_options else {}
                
                ####################################
                #Predict all family collections in a folder, but excluding the correct family:
                if "validate" in args.beta_options:
                    args.exclude_family = {os.path.basename(args.query).split(".")[0]}
                    
                ####################################
                #Define namespace for prediction:
                opts = argparse.Namespace(all_peptides=all_peptides, meta=CUPPmeta, overlap=args.overlap, n_mer=args.n_mer, ambiguous=args.ambiguous,
                    line_cut=args.position_cut, precision_ratio = args.best_hit,occurence = args.occurence, type = args.type, cupp_minimum_score = args.cupp_minimum_score,
                    minimum_cup = args.minimum_cup, evidence = args.evidence, domain_min = args.domain_min, relatedness = False,
                    exclude_family = args.exclude_family, beta_options = args.beta_options, keep_only = args.keep_only)
                
                ####################################
                #Open handles:
                seq = ""; acc = ""; count = [0,0]; start_time = time.time()
                if not args.output_path:
                    args.raw_path = "%s/predicted/%s_CUPP.log" % (args.working_dir,".".join(os.path.basename(args.query).replace(".gz","").split(".")[:-1]))
                    args.fasta_path = "%s/predicted/%s_CUPP%s" % (args.working_dir,".".join(os.path.basename(args.query).replace(".gz","").split(".")[:-1]),args.fasta_extension)
                else:
                    args.raw_path = "%s.log" % (args.output_path)
                    args.fasta_path = args.output_path
                handle_fas = open(args.fasta_path,"w" if first else "a") if args.type in ["both","fasta","quick"] else False
                handle_raw = open(args.raw_path,"w" if first else "a") if args.type in ["both","raw"] else False
                handle = open(args.query) if args.query.split(".")[-1] != "gz" else gzip.open(args.query,"rb")
                
                ####################################
                #The query can be genomic or protein and optionally as gz compressed:
                for line in handle:
                    line = line.strip() if args.query.split(".")[-1] != "gz" else line.decode("utf-8").strip() 
                    if not line:
                        continue
                    
                    ####################################
                    #Loop over each entry individually:
                    if line[0] == ">":
                        if seq:
                            if args.genomic:
                                for a,s in orf(acc,seq,codontable,minimum=args.genomic).items():
                                    count[0] += out(handle_fas,handle_raw,a,predict(s,opts),s,args.beta_options,CUPPmeta,sep=args.sep); count[1] += 1
                            else:
                                if args.max_longst < len(seq):
                                    sys.stderr.write(" !  Did not expect a single protein longer than %s aa, yours is %s: %s\n" % (args.max_longst,len(seq),acc.split(args.sep)[0]))
                                else:
                                    count[0] += out(handle_fas,handle_raw,acc,predict(seq,opts),seq,args.beta_options,CUPPmeta,sep=args.sep); count[1] += 1
                            seq = ""
                        acc = line[1:]
                    elif acc:
                        seq += line.upper()
                        
                ####################################
                #Process the last entry:
                if acc and seq:
                    if args.genomic:
                        for a,s in orf(acc,seq,codontable,minimum=args.genomic).items():
                            count[0] += out(handle_fas,handle_raw,a,predict(s,opts),s,args.beta_options,CUPPmeta,sep=args.sep); count[1] += 1
                    else:
                        ####################################
                        #Check if the FASTA file is occupied in a meaningful way leaving one very large protein:
                        if args.max_longst < len(seq):
                            sys.stderr.write(" !  Did not expect a single protein longer than %s aa, yours is %s: %s\n" % (args.max_longst,len(seq),acc.split(args.sep)[0]))
                        else:
                            count[0] += out(handle_fas,handle_raw,acc,predict(seq,opts),seq,args.beta_options,CUPPmeta,sep=args.sep); count[1] += 1
                    acc = ""; seq = ""
                    
                handle_fas.close() if handle_fas else None
                handle_raw.close() if handle_raw else None
                already_processed[args.query] = os.path.getmtime(args.query)
                print ("### Prediction COMPLETED on %s proteins, found: %s (%s [s]): %s" % (count[1],count[0],round(time.time() - start_time,2), args.query))
                count = [0,0]
                
            ####################################
            #Format the new arguments, if any:
            if args.live or args.dir_query:
                new_query_coming = True
                while new_query_coming and new_query_coming != "EXIT":
                    try:
                    
                        ####################################
                        #Check if unprocessed collection exists in folder:
                        if sum(args.dir_query.values()):
                            args.query = sorted([path for path,waiting in args.dir_query.items() if waiting])[0]
                            
                            ####################################
                            #Merge the output of several collections into one output file:
                            if args.output_path:
                                first = False
                            
                            ####################################
                            #Flag the current collection:
                            args.dir_query[args.query] = False
                            
                            ####################################
                            #No new collection needed, already have one:
                            new_query_coming = False
                            
                            ####################################
                            #Output status line:
                            print ("--> %s" % args.query) if args.query not in already_processed or os.path.getmtime(args.query) != already_processed[args.query] else None
                            
                        ####################################
                        #Terminate the script or ask for new user arguments:
                        elif args.live:
                            arguments = input("### Specify new arguments using the loaded CUPP library:\n")
                            if "EXIT" == arguments.upper()[0:4]:
                                new_query_coming = "EXIT"
                            else:
                                args = arg_parser_predict(sys.argv[1:] + arguments.split())
                                if (not args.query or not os.path.exists(args.query)) and not sum(args.dir_query.values()):
                                    sys.stderr.write("!!! The qeury path does not exist, try again!\n")
                                else:
                                    new_query_coming = False
                        else:
                            new_query_coming = "EXIT"

                    except:
                        sys.stderr.write("!!! The arguments were not valid, try again or type EXIT!\n")
                if new_query_coming == "EXIT":
                    cupp_logo()
                    sys.exit("### See you soon, CUPP'er...")
            else:
                next = False
    else:
        sys.stderr.write(" !  Please specify a valid query fasta file for CUPP prediction\n")
        
