#!/usr/bin/python
import os,sys,time,json,psutil,argparse
from CUPPprediction_DIRECT import *
from scipy.cluster.hierarchy import fcluster
from multiprocessing import Pool
from statistics import median
import numpy as np
sys.stdout = Unbuffered(sys.stdout) #Try to make the output better!
sys.stderr = Unbuffered(sys.stderr) #Try to make the output better!


def arg_parser(args): 
    '''
    Process command-line user arguments and prepare for further directories
    '''

    ####################################
    #Check user argument range:
    def lower_float(x):
        try:
            if float(x) < 0.0 or 1.0 < float(x):
                raise argparse.ArgumentTypeError("\n!!! The specified float is not between zero and one: %s" % (x))
        except:
            raise argparse.ArgumentTypeError("\n!!! The specified value can't be converted into a float: %s" % (x))

        return float(x)
        
    ####################################
    #Folders and files:
    parser = argparse.ArgumentParser(description='The input parameters for the CUPP identification technology')
    parser.add_argument('-working_dir',       default=os.getcwd()+"/CUPP",type=str,       help="Change working directory for both clustering and prediction")
    parser.add_argument('-database_version',  default="",        type=str,                help="Suffix on clusterings result folders")
    parser.add_argument('-fasta_extension',   default=".faa",    type=str,                help="Define the fasta extension to not be included in the family abbreviation")
    parser.add_argument('-cdhit_folder',      default="cdhit_2019_5c",type=str,           help='Folder for cdhit clsr files, NOTE: Accession can only be 17 chars long ending with "|"')
    parser.add_argument('-dbcan_folder',      default="dbcan_2019_5c",type=str,           help="Folder for domain using same accession numbers with prefix dot and a number")
    parser.add_argument('-deposit',           default="cazy_2019_5_ec",type=str,          help="Folder of the original full length proteins")
    
    ####################################
    #CUPPclustering parameters:
    parser.add_argument('-clustering',        default=False,     action="store_true",      help="Perform CUPPclustering on the proteins or domains")
    parser.add_argument('-user_ncbi',         default=False,     action="store_true",      help="Apply user defined NCBI Taxonomy Database")
    parser.add_argument('-header',            default=["Accession","Family","Subfam","raw_function","Uniprot","PDB","Name","Taxid"], nargs='+', help='The order of the FASTA header, with sep = "|" as defualt')
    parser.add_argument('-common',            default=[],        nargs='+',                help="Select families by name followed by .faa in deposit folder, (need to be identical to position two of the header for validation)")
    parser.add_argument('-regular_exp',       default="",        type=str,                 help="Specify a search string to locate relevante FASTA files in a folder")
    parser.add_argument('-domain_query',      default="",        type=str,                 help="Specify the domain FASTA file for clustering and a full length protein for prediction")
    parser.add_argument('-ambiguous',         default=2,         type=int,                 help="Introduce a number of X's into each peptides as ambiguous amino acids")
    parser.add_argument('-n_mer',             default=8,         type=int,                 help="Peptide length for both clustering and prediction")
    parser.add_argument('-redo',              default=False,     action="store_true",      help='Ignore existing CUPPclusterings and override them')
    parser.add_argument('-gradient',          default=[0.1,0.2,0.3,0.4], nargs='+',        help="Conservedness cutoff in each iteration")
    parser.add_argument('-cdhit',             default=0.90,      type=float,               help="Use only representative proteins but meta data from all entries in CD-HIT cluster")
    parser.add_argument('-redundancy_path',   default="",        type=str,                 help="Specify a direct path to a CDHIT cluster file (overrides)")
    parser.add_argument('-domain_off',        default=False,     action="store_true",      help='Use full protein fasta only instead of both')
    parser.add_argument('-domains_only',      default=False,     action="store_true",      help='Use full protein fasta only instead of both')
    parser.add_argument('-cc',                default=9,         type=float,               help="Amplification of clustering c_clust")
    parser.add_argument('-outlier',           default=10,        type=float,               help="Ignore proteins with less than n percent of the median conserved peptides or less than n peptides")
    parser.add_argument('-start_positions',   default=20,        type=float,               help="Ignore proteins with less than n peptides from different positions")
    parser.add_argument('-minimum_group_size',default=5,         type=float,               help="Ignore groups with this number of protein members")
    parser.add_argument('-cup_cut',           default=0.7,       type=lower_float,         help="Define the threshold for when two CUPPgroups should be merged rather than stay as two based on ratio of shared peptide frequencies")
    parser.add_argument('-variable_cut',      default=False,     action="store_true",      help='Evaluate the best suited threshold for protein group formation and replaces the step using cup_cut (the parameter "final_threshold" was also become unused')
    parser.add_argument('-guide',             default="",        type=str,                 help="Restrict the clustering to only include peptides covered by peptides of a specified domain")
    parser.add_argument('-mix',               default="",        type=str,                 help="A name of the model e.g. if several families are included")
    
    ####################################
    #Relatedness parameters:
    parser.add_argument('-rel_cut',           default=0.2,       type=lower_float,         help="Threshold for making a connection between a unknown group and the groups having meta data based on shared position coverage")
    parser.add_argument('-cupp_cut_last',     default=0,         type=lower_float,         help="Threshold for naming to give an indication of related CUPP groups 47.1 and 48.1 vs 47.1 and 47.2")
    parser.add_argument('-limit',             default=0,         type=int,                 help="Reduce the size of the CUPPlibrary by only keeping peptides with a combined sum below a set limit")
    
    ####################################
    #CUPPlibrary parameters:
    parser.add_argument('-recompile',         default=False,     action="store_true",      help="Compile CUPPgroups of selected families into one CUPPlibrary")
    parser.add_argument('-final_cut',         default=0.2,       type=lower_float,         help="Frequency cutoff for peptides obtained after clustering")
    parser.add_argument('-final_group_size',  default=5,         type=int,                 help="Keep the small groups after clustering")
    parser.add_argument('-final_threshold',   default=0,         type=int,                 help="Alter the default threshold for the linkage of the last round of clustering")
    parser.add_argument('-mono_family',       default=True,      action="store_false",     help="If only a single function is known to a family, assume all its CUPP groups hold that EC function")
    
    ####################################
    #CUPPprediction parameters:
    parser.add_argument('-compiled_json',     default="",        type=str,                 help="Specify full path to CUPPlibrary")
    parser.add_argument('-predict',           default=False,     action="store_true",      help="Perform CUPPprediction")
    parser.add_argument('-query',             default="",        type=str,                 help="Select a fasta file for CUPP clustering or prediction")
    parser.add_argument('-dir_query',         default="",        type=str,                 help="Select a folder CUPP compile")
    parser.add_argument('-best_hit',          default=0.5,       type=lower_float,         help="When to consider a secondary hit equally relevant as the first")
    parser.add_argument('-occurence',         default=False,     action="store_true",      help="Retrun only the EC number which the highest occurence of the CUPPgroup")
    parser.add_argument('-overlap',           default=0.6,       type=lower_float,         help="How much can a domain overlap with an already assigned domain")
    parser.add_argument('-minimum_cup',       default=5,         type=float,               help="Required frequency for positive hit")
    parser.add_argument('-evidence',          default=0.01,      type=lower_float,         help="Minimum percentage of frequency of CUPPgroup for hit")
    parser.add_argument('-position_cut',      default=0.2,       type=lower_float,         help="Frequency for a position to be included in the domain")
    parser.add_argument('-domain_min',        default=20,        type=int,                 help="Minimum number of positions covered be a domain")
    parser.add_argument('-output_path',       default="",        type=str,                 help="Specify a full path and filename of the output fasta file (without extension)")
    parser.add_argument('-genomic',           default=0,         type=int,                 help="Treat the query as DNA and locate ORFS and predict them individually")
    parser.add_argument('-cupp_minimum_score',default=5,         type=int,                 help="The minimum score for a assignment of CUPP group and its EC function to a query protein")
    
    ####################################
    #Parallel computing and documentation:
    parser.add_argument('-pipe',              default=False,     action="store_true",      help='Avoid user questions (auto reply "Y") and do not display plots')
    parser.add_argument('-silent',            default=False,     action="store_true",      help="Do not show warnings!")
    parser.add_argument('-time',              default=int(psutil.Process(os.getpid()).create_time())-1530000000, type=int, help="Time of Creation")
    parser.add_argument('-version',           default="v1.1.1",  type=str,                 help="The folder of fastafiles from which clusters are made")
    parser.add_argument('-ignored',           default=[],        nargs='+',                help="Manually ignore listed accession numbers and all members of their CDHIT clusters")
    parser.add_argument('-jobs',              default=4,         type=int,                 help='Number of processor for use for clustering')
    parser.add_argument('-dendrogram',        default=False,     action="store_true",      help='Print the dendrograms in the end of the run or load previously saved plots')
    parser.add_argument('-beta_options',      default=[],        nargs='+',                help="New features under development",choices=["merge","debug","mix_type","no_fam_ec","keep_ranges","complex","fragments","float","experimental_families","keep","catalytic","relatedness_alone","prediction_fails","export"])
    parser.add_argument('-hpc',               default=[],        nargs='+',                help="State the mail for which to receive notifications, second element is the RAM, third is the run-time in hours")
    parser.add_argument('-type',              default="both",    type=str,                 help='Select output formats',choices=["none","quick","fasta","raw","both"])
    parser.add_argument('-sep',               default="|",       type=str,                 help='Select the separator in the header for both input and output')
    parser.add_argument('-initiated',         default=False,     action="store_true",      help='Used internally for serial or parallel clustering')
    parser.add_argument('-live',              default=False,     action="store_true",      help='Keep the CUPP library loaded and predicted in serial FASTA files as query')
    parser.add_argument('-exclude_family',    default=[],        nargs='+',                help="Select families which will be ignored from the prediction")
    parser.add_argument('-keep_only',         default=[],        nargs='+',                help="Select few families which are the only ones predicted")
    parser.add_argument('-export',            default="",        type=str,                 help="Number of proteins to output from each CUPP group")
    args = parser.parse_args(args)

    ####################################
    #Logic consequences of an user options:
    args.export = int(args.export) if str(args.export).isdigit() else ""
    args.ignored = dict((k,"leave out") for k in args.ignored)
    args.working_dir = os.path.normpath(args.working_dir)
    args.common = [os.path.basename(args.query).split(".")[0]] if args.query else args.common
    args.clustering = False if args.recompile else args.clustering
    try:
        args.gradient = [float(f) for f in args.gradient]+[float(args.gradient[-1])]
    except:
        sys.exit("!!! The specified gradient value can not be converted into a float: %s" % args.gradient)
    args.database_version = "_"+args.database_version if args.database_version and args.database_version[0] != "_" else args.database_version 
    args.predict = True if args.query and not args.clustering and not args.recompile and not args.dendrogram else args.predict
    args.predict = False if args.clustering else args.predict
    args.setting = "%s%s%sx%s_%s" % ("f" if args.domain_off else "d","c" if args.cdhit else "a",args.n_mer,args.ambiguous,int(args.cc*10))
    args.final_group_size = args.minimum_group_size if args.final_group_size == 0 else args.final_group_size
    args.beta_options = set(args.beta_options)
    args.itol_out = []
    args.ref_old = ""#CUPP/CUPPlibrary/8x2_90_AA9_OLD2_CUPPlibrary.json"
    args.fasta_extension = "."+args.fasta_extension if "." not in args.fasta_extension else args.fasta_extension

    ####################################
    #Error handling:
    if args.domain_off and args.domains_only:
        sys.exit("!!! You can not have both domain_off and domain_only at the same time")
    
    ####################################
    #During prediction, some families can be left out [only for cazy_6_2018]:
    args.exclude_family = set(args.exclude_family)
            
    ####################################
    #Accession must be in header:
    if "Accession" not in args.header:
        sys.exit('!!! The specified header order must contain the column "Accession"')
            
    ####################################
    #Hierarchy of options, recompile is higher, clustering is next and predict is last: 
    if sum([1 < f for f in args.gradient]):
        sys.exit(" !  You can not have a cluster coefficient (cc) above 1!")
    
    ####################################
    #Preperation: of folders:
    os.makedirs(args.working_dir+"/predicted", exist_ok=True)
    os.makedirs("%s/output" % args.working_dir, exist_ok=True)    
    os.makedirs("%s/CUPPpools%s" % (args.working_dir,args.database_version), exist_ok=True)
    os.makedirs("%s/graphs%s" % (args.working_dir,args.database_version), exist_ok=True)
    os.makedirs("%s/CUPPlibrary" % args.working_dir, exist_ok=True)
    os.makedirs("%s/itol%s" % (args.working_dir,args.database_version), exist_ok=True)                  
    
    ####################################
    #Process several families for clustering:
    if args.regular_exp and not args.initiated:
        args.common = common_abb(args.regular_exp,args.deposit,extension=args.fasta_extension)
        print ("### The regular expression found common names: %s" % len(args.common))
        if not args.common:
            sys.exit('!!! BAD regular expression for deposit(%s): "%s"' % (args.deposit,args.regular_exp))
    if not args.common and (not args.mix or not args.query):
        sys.exit("!!! Please specify a common name for the current collection")
    
    ####################################
    #Name of the collections will be set to common if only one family:
    args.general_name = args.common[0] if not args.mix else args.mix
        
    ####################################
    #Process each family individually as serial or on hpc for parallel:
    if args.clustering and ((len(args.common) != 1 and not args.mix) or args.hpc) and not args.initiated:
        if args.common:
            if 1 < len(args.common) and not args.pipe:
                serial_answer = input(" ?  ARE YOU SURE YOU WANT TO PROCESS SEVERAL PROTEIN COLLECTIONS (%s) IN A SERIAL (N/Y) ?\n" % len(args.common))
            else:
                serial_answer = "Y"
            if serial_answer and serial_answer.upper()[0] == "Y":
                for common in args.common:
                    if not args.hpc:
                        os.system("python %s -common %s -initiated" % (" ".join(sys.argv),common))
                        print ("--> python %s -common %s" % (" ".join(sys.argv),common))
                    else:
                        submit_hpc("python %s -common %s -initiated" % (" ".join(sys.argv),common),hours=48, ram = 50, mail=args.hpc, jobs = args.jobs, ID = common)
                print ("### Full run-time:\t%s [s]" % str(round(time.time() - psutil.Process(os.getpid()).create_time(),1)))
                sys.exit("### The submission of different families have been successful!")
        sys.exit("### Please exactly one family name as the common (%s)" % len(args.common))
        
    ####################################
    #Prepare output paths:
    args.raw_path = "%s/predicted/%s_CUPP.log" % (args.working_dir,args.general_name) if not args.output_path else args.output_path+".log"
    args.fasta_path = "%s/predicted/%s_CUPP%s" % (args.working_dir,args.general_name,args.fasta_extension) if not args.output_path else args.output_path
    if "." not in os.path.basename(args.fasta_path):
        args.fasta_path += args.fasta_extension
    
    ####################################
    #Avoid overwriting the CUPlibrary accidentally:
    args.compiled_json = "%s/CUPPlibrary/%s_%s_CUPPlibrary.json" % (args.working_dir,args.setting[2:],args.version) if not args.compiled_json else args.compiled_json
    if args.clustering:
        args.compiled_json = "%s/output/%s_CUPPlibrary.json" % (args.working_dir,tuple(args.common)[0])
        args.recompile = True
        args.predict = True
    elif args.recompile:
        if os.path.exists(args.compiled_json) and not args.pipe:
            recompile_answer = input(" P  %s\n ?  ARE YOU SURE YOU WANT TO OVERRIDE THE COMPILED CUPP LIBRARY (N/Y) ?\n" % args.compiled_json)
            if not recompile_answer or recompile_answer[0].upper() != "Y":
                args.recompile = False
                
    return (args)


if __name__ == "__main__":

    ####################################
    #User inputs:
    args = arg_parser(sys.argv[1:])

    ####################################
    #Determine the relevant paths for collections to include:
    domain_list, original_list, used_list, cdhit_list = determine_paths(args.common, args.deposit,
            query=args.query, 
            fasta_extension=args.fasta_extension, 
            dbcan_folder=args.dbcan_folder, 
            cdhit_folder=args.cdhit_folder,
            redundancy_path = args.redundancy_path, 
            domain_off=args.domain_off,
            domains_only=args.domains_only,
            cdhit=args.cdhit, 
            beta_options = args.beta_options)
    if not args.domain_query:
        domain_list, original_list, used_list, cdhit_list = determine_paths(args.common, args.deposit,
            query=args.query, 
            fasta_extension=args.fasta_extension, 
            dbcan_folder=args.dbcan_folder, 
            cdhit_folder=args.cdhit_folder,
            redundancy_path = args.redundancy_path, 
            domain_off=args.domain_off,
            domains_only=args.domains_only,
            cdhit=args.cdhit, 
            beta_options = args.beta_options)
    else:
        domain_list, original_list, used_list, cdhit_list = [args.domain_query], [args.domain_query], [args.domain_query], ['']
        
    ####################################
    #For current clustering:
    CUPPpool_path_list = []
    CUPPpool_path = os.path.join(args.working_dir,"CUPPpools%s" % args.database_version,args.general_name+"_CUPPpool.json")
    
    ####################################
    #Export the top ranking members of each CUPP group within an already processed family:
    if args.export != "":
        export_CUPP_members(CUPPpool_path,domain_list,cdhit_list,args.export)
        print ("### The domains of the high ranking members of each CUPP group have been exported")
        sys.exit()
        
    ####################################
    #Merge several family collections into one common collection of non-redundant proteins:
    if "merge" in args.beta_options:
        merge_fasta_dir(original_list)
        print ("### The %s collections have been merged" % len(original_list))
        sys.exit()
        
    ####################################
    #Change the setting according to changes:
    check = "c" if sum([1 if cdhit_list[n] else 0 for n in range(len(cdhit_list))]) == len(cdhit_list) else "a" if sum([1 if cdhit_list[n] else 0 for n in range(len(cdhit_list))]) == 0 else "m"
    args.setting = "%s%s%sx%s_%s" % ("f" if args.domain_off else "d",check,args.n_mer,args.ambiguous,int(args.cc*10))
    #print ("!!! The setting for the individual common families are not handled yet, which may cause problems!!!")
        
    #####################################
    ### CUPPclustering of collections ###
    #####################################
    start = time.time()
    if args.clustering and (not check_existence(CUPPpool_path,args.setting) or args.redo):

        ####################################
        #Printing current settings:
        show_arguments(args,cluster=True)
        
        ####################################
        #Load NCBI taxonomy nodes and names:
        if args.user_ncbi:
            os.makedirs("resources", exist_ok=True)
            ncbi = taxonomy_database(nodes_file="resources/nodes.dmp",names_file="resources/names.dmp",old_file="resources/merged.dmp",delete_file="resources/delnodes.dmp")
        elif not os.path.exists("resources/ncbi_short.json"):
            ncbi = {"nodes":{"1":["1","root"]},"names":{"1":["root"]},"partial":0,"banned":{}}
        else:
            with open("resources/ncbi_short.json") as json_data:
                ncbi = json.load(json_data)
        
        ####################################
        #Obtain a collection of sequences from a fasta file and high similarity cluster-info (CDHIT):
        redundant_proteins = load_redundancy(cdhit_list,args.common) if cdhit_list[0] else {}
        
        ####################################
        #Create a dict of only the representative protein of each CD-HIT cluster:
        group_rep = {}
        for acc,v in redundant_proteins.items():
            if v[1]:
                if v[0] not in group_rep:
                    group_rep[v[0]] = acc
                else:
                    print ("!!! ERROR IN CDHIT GROUPS %s" % v[0])
                
        collection, fasta_meta = obtain_collection(used_list, cdhit = redundant_proteins, sep = args.sep, common = args.general_name, meta_cat = args.header, beta = args.beta_options)
        print ("### Proteins in collection: %s" % len(collection))
        
        if len(collection) == 1:
            sys.exit("!!! The family does only hold a single member")
        elif len(collection) == 0:
            sys.exit("!!! The family does not contain any proteins")
        
        ####################################
        #Obtain peptides of each protein as a bag-of-words (peptide pools):
        original_seq_pep = obtain_peptides(collection,use_positions=False,n_mer=args.n_mer,ambiguous=args.ambiguous,beta=args.beta_options)
        print ("### Preparation of clustering:\t%s [s]" %(round(time.time() - start,1)))
        
        ####################################
        #Format the meta data of proteins and add taxonomy info:
        collection, meta = format_meta(collection, fasta_meta, ncbi=ncbi, meta_cat=args.header)
        meta_cat = args.header+([] if "Family" in args.header else ["Family"])+["Classes","Taxonomy","Add"]
        
        ####################################
        #Save conserved peptides between the iterations of the incremental clustering:
        all_peptides = {}
        
        ####################################
        #Incremental CUPP clustering:
        for current_iteration in range(len(args.gradient)):
                    
            ####################################
            #Round specific parameters:
            final_round = current_iteration == len(args.gradient)-1; round_start = time.time()
            min_size = args.minimum_group_size + ((2+current_iteration-len(args.gradient))) if not final_round else args.final_group_size
            print ("#"*28+"  "+str(current_iteration+1)+" "+"#"*28)
            print ("### The minimum group size:\t%s" % min_size)
            print ("### Peptide conservation:\t%s" % args.gradient[current_iteration])

            ####################################
            #Remove peptides from the peptide pools if they are not considered important: 
            seq_pep = dict((acc,dict((p,n) for p,n in peptides.items() if p in all_peptides)) for acc,peptides in original_seq_pep.items()) if current_iteration else original_seq_pep.copy()
            
            ####################################
            #Determine the number of positions in which a conserved peptide originates:
            seq_length = dict((acc,len(set((j) for i,j in v.items()))) for acc,v in seq_pep.items())
            
            ####################################
            #Ignore protein sequences which do not belong in the collection:
            args.cut_by_median = median(seq_length.values()) if seq_length else 0
            coverage_outliers = dict((k,v) for k,v in seq_length.items() if v < args.start_positions)
            for seq,acc in collection.items():
                if acc in coverage_outliers:
                    covered_positions = set()
                    for pep,pos in seq_pep[acc].items():
                        for n in range(args.n_mer):
                            if pep[n] != "X":
                                covered_positions.add(pos+n)
                    coverage_outliers[acc] = len(covered_positions)
            args.ignored.update(dict((k,"%s of %s coverage" % (v,args.start_positions)) for k,v in coverage_outliers.items() if v < args.start_positions))
            args.ignored.update(dict((k,"%s of %s median" % (v,args.cut_by_median)) for k,v in seq_length.items() if v < args.cut_by_median*args.outlier/100 ))
            seq_pep = dict((k,v) for k,v in seq_pep.items() if k not in args.ignored)
            print ("### Outlier cut-off:\t\t%s" % (args.cut_by_median*args.outlier/100))
            print ("### Proteins in collection:\t%s" % len(seq_pep))
            
            ####################################
            if len(seq_pep) <= 1:
                sys.stderr.write("!!! Too few proteins in family to form protein groups: %s\n" % ("%s\n" % used_list[0]) if len(used_list) == 1 else "in any of the collection (%s)\n" % len(used_list))
                sys.exit()
            
            ####################################
            #Obtain the protein groups based on the conserved peptides:
            labels = dict(zip(sorted(seq_pep),range(len(seq_pep))))
            links = obtain_groups(seq_pep,n_jobs=args.jobs,cc=args.cc)
            if not final_round or not args.variable_cut:
                clusters = fcluster(links, 1 if not final_round or not args.final_threshold else args.final_threshold, criterion='distance')
            else:
                #clusters = variable_groups(links,sorted(seq_pep),seq_pep,cup_cut=args.gradient[current_iteration],minimum=args.start_positions,minimum_group=min_size)
                clusters = variable_groups(links,sorted(seq_pep),seq_pep,cup_cut=args.gradient[current_iteration],minimum=args.start_positions,minimum_group=min_size)
            fam_groups = ["%s:%s" % (args.general_name,k) for k in clusters]
            lookup_acc_group,meta_merge = pool_groups(fam_groups,sorted(seq_pep),meta,meta_cat)
            
            ####################################
            #Keep only entries which are not removed as outliers:
            collection = dict((seq,acc) for seq,acc in collection.items() if acc in seq_pep)
            meta = dict((acc,v) for acc,v in meta.items() if acc in seq_pep)
        
            ####################################
            #Obtain bags-of-words for each of the protein groups:
            opportunists, CUPPtides, bean_clusters, good_clusters = powderize(collection,original_seq_pep,lookup_acc_group,min_size,args.n_mer,args.ambiguous,args.gradient[current_iteration])
                        
            ####################################
            #Do not continue if no conserved peptides could be found:
            if not CUPPtides:
                sys.exit(" !  TERMINATED - Nothing can be grouped with the current collection!")
            args.cup_cut = args.cupp_cut_last if final_round and args.cupp_cut_last else args.cup_cut
            assign_cup,CUPPlist = relatedness(CUPPtides,cup_cut=args.cup_cut)
            
            ####################################
            #Locate accession numbers of ignored clusters:
            ignored_groups = {}
            for group,v in bean_clusters.items():
                if group not in good_clusters:
                    for a in meta_merge[group][args.header.index("Accession")]:
                        if len(v) < min_size:
                            args.ignored[a] = "%s of %s group_size" % (len(v),min_size)
                        else:
                            args.ignored[a] = "%s group with no_peptides" % (len(v))
                for acc in  meta_merge[group][args.header.index("Accession")]:
                    if acc in meta and acc in opportunists:
                        if opportunists[acc] < args.start_positions and final_round:
                            args.ignored[acc] = "%s peptides of %s as minimum CUPPgroup size" % (opportunists[acc],args.start_positions)
                    else:
                        args.ignored[acc] = "%s not in meta" % (acc)
                        
            ####################################
            #Rename the CUPP groups and remove protein which are considered outliers:
            long_assign, long_label = rename(CUPPlist,meta_merge,meta_cat,args.general_name,opportunists,assign_cup,args.ignored,args.start_positions,last=final_round)
            lookup_acc_group,meta_merge = pool_groups(long_assign,long_label,meta,meta_cat)
            opportunists2, CUPPtides, bean_clusters, good_clusters = powderize(collection,original_seq_pep,lookup_acc_group,min_size,args.n_mer,args.ambiguous,args.final_cut if final_round else args.gradient[current_iteration],last=True)
            print ("### Protein groups:\t%s" % len(meta_merge))
                    
            ####################################
            #Making a dendrogram of the sequences, create labels [THE RESTORED CUPP GROUP NUMBERS ARE NOT UPDATED IN DENDROGRAM]:
            meta_label_dict = dict((k,"") for k,v in labels.items())
            for k in labels:
            
                if k in lookup_acc_group:
                    #meta_label_dict[k] += lookup_acc_group[k][0].split(".")[0] if lookup_acc_group[k][0].split(".")[0]+".2" not in meta_merge else lookup_acc_group[k][0]
                    meta_label_dict[k] += lookup_acc_group[k][0] #WARNING
                if k in meta:
                    meta_label_dict[k] += "_"+("*" if "PDB" in meta[k] and meta[k]["PDB"] else "") + str(seq_length[k]) + "_" + str(opportunists[k] if k in opportunists else "0")
                    meta_label_dict[k] += "_"+k#[i for i,j in sorted(meta[k]["Accession"].items(),key=lambda x: x[1],reverse=True)][0].split(".")[0]
                    meta_label_dict[k] += "_"+[i for i,j in sorted(meta[k]["Classes"].items(),key=lambda x: x[1],reverse=True)][0][:5] if meta[k]["Classes"] else "_"
                    meta_label_dict[k] += "_"+"+".join([i for i,j in sorted(meta[k]["Subfam"].items(),key=lambda x: x[1],reverse=True)]) if "Subfam" in meta[k] else "_"
                    meta_label_dict[k] += "_"+"+".join([i for i,j in sorted(pool_pool(meta[k]["raw_function"]).items(),key=lambda x: x[1],reverse=True)]) if "raw_function" in meta[k] else "_"
                else:
                    meta_label_dict[k] += "_"+k
                    
            ####################################
            #Construction of the actual dendrogram:
            plot_fil = "%s/graphs%s/%s_%s_%s_acc.json" % (args.working_dir,args.database_version,args.general_name,args.setting,current_iteration)
            dendro(links,[meta_label_dict[l] for l,v in sorted(labels.items(),key=lambda x: x[1])],
                algo="ward",
                cc=1,
                name=plot_fil,
                fil=plot_fil,
                show=False,
                note=args.ignored)
            
            ####################################
            #Save the peptides and meta data as a json object:
            if final_round:

                ####################################
                #Save score for each accession number telling how much that proteins describe the group:
                meta_cat.append("Freq_sum")
                for group,v in sorted(meta_merge.items()):
                    for acc,score in v[args.header.index("Accession")].items():
                        meta_merge[group][args.header.index("Accession")][acc] = group_rep[redundant_proteins[acc][0]] if not score else opportunists[acc] if acc in opportunists and acc not in args.ignored else 0

                ####################################
                #Estimate the expected range of domain:
                for cupp,v in meta_merge.items():
                    ranges = []; scores = []
                    for acc,s in v[0].items():
                        if acc in seq_pep and seq_pep[acc] and str(s).isdigit():
                            ranges.append(max(seq_pep[acc].values())-min(seq_pep[acc].values()))
                            scores.append(s)
                    meta_merge[cupp].append({}) #Freq_sum
                    meta_merge[cupp][-1]["Minimum"] = round(min(ranges))
                    meta_merge[cupp][-1]["Average"] = round(sum(ranges)/len(ranges))
                    meta_merge[cupp][-1]["Maximum"] = round(max(ranges))
                    meta_merge[cupp][-1]["Coverage"] = round(sum(scores)/len(scores))
                    
                ####################################
                #Save the CUPP groupings:
                powder = {args.setting:{"peptides":CUPPtides,"meta":meta_merge,"meta_categories":meta_cat}}
                handle_out = open(CUPPpool_path,'w')
                json.dump(powder,handle_out,indent=3)
                handle_out.close()
                CUPPpool_path_list = [CUPPpool_path]
                print ("### The CUPPgroups have been saved: %s" % CUPPpool_path)
                
                ####################################
                #Save den dendrogram as Newick format for handling e.g. in iTOL:
                custom_colors = {}
                newick_path = "%s/itol%s/%s_%s.tree" % (args.working_dir,args.database_version,args.general_name,args.setting)            
                itol(args.itol_out,meta_merge,args.general_name,args.header,meta,newick_path,labels=[l.split(":")[0] for l,v in sorted(labels.items(),key=lambda x: x[1])],links=links,custom_colors = custom_colors,path=args.working_dir)
            else:  
                all_peptides = dict((pep,1) for pep in CUPPtides)
            print ("### The clustering run time:\t%s" % round(time.time() - round_start,3))

        print ("### Clustering COMPLETED:\t%s [s]" % str(round(time.time() - start,3)))
        print ("#"*60+"\n"+"#"*60)

    ####################################
    #Identify location of all relevant families:
    if args.recompile and not args.clustering:
        if not args.mix:
            for common in args.common:
                CUPPpool_path_list.append(os.path.join(args.working_dir,"CUPPpools%s" % args.database_version,common+"_CUPPpool.json"))
        else:
            CUPPpool_path_list.append(os.path.join(args.working_dir,"CUPPpools%s" % args.database_version,args.general_name+"_CUPPpool.json"))
    
    ####################################
    #Do relatedness without rerun:
    if "relatedness_alone" in args.beta_options and args.clustering:
        CUPPpool_path_list = [CUPPpool_path]
        
    ####################################
    #Specify a path to a single pool: NEW FEATURE, LOAD A FILE OF PATHS!!
    if args.query:
        CUPPpool_path_list = [args.query]
        
    ####################################
    #Specify a path to a single pool: NEW FEATURE, LOAD A FOLDER from PATHS!!
    if args.dir_query:
        CUPPpool_path_list = [os.path.join(args.dir_query,f) for f in os.listdir(args.dir_query)]
        
    ####################################
    #Check if the powder list is not emtry:
    if (not CUPPpool_path_list and (args.clustering or args.recompile)):
        args.clustering, args.recompile, args.predict = False, False, False
        
    ####################################
    ###    COMPILING CUPP Library    ###
    ####################################
    if args.recompile:

        ####################################
        #Merge one or more protein families into one large dictionary:
        print ("### The individual CUPPgroups are being merged for: %s (Families: %s)" % (args.setting,len(CUPPpool_path_list)))
        all_peptides = {}; all_meta = {}; meta_cat = []; experimentals = {}
        for js in CUPPpool_path_list:
            common = os.path.basename(js).replace("_CUPPpool.json","")
            if check_existence(js,args.setting):
                load_powder(js,all_peptides,all_meta,meta_cat,args.setting,experimentals,mimimum_group=args.final_group_size,rel_cut=args.rel_cut,final_cut=args.final_cut,mono_family=args.mono_family,limit=args.limit,aa=args.n_mer-args.ambiguous,special=args.beta_options)
            
            ####################################
            #Do not start to cluster during CUPPlibrary building:
            else:
                sys.stderr.write(" !  The family %s has not yet been clustered and disregarded in the CUPPlibrary\n" % common)
        
        ####################################
        #Inspect header:
        if [f for f in meta_cat if f not in ["Family","Classes","Taxonomy","Add","Freq_sum"]] != [f for f in args.header if f not in ["Family"]]:
            sys.stderr.write("!!! The headers of the individual protein families is not the same as the one in arguments\n")

        ####################################
        #Determine if peptides are shared across several families and down rank them:
        print ("#"*60+"\n### Adding important stats to the CUPPgroups (peptides %s)" % len(all_peptides))
        total_families = len(set((k.split(":")[0]) for k in all_meta))
        CUPPgroups = set((k) for k in all_meta)
        banned = {}; general = {}
        for pep,groups in all_peptides.items():
            
            ####################################
            #Family level:
            current_families = len(set((k.split(":")[0]) for k in groups))
            
            ####################################
            #Divide the score if two families share the same peptides:
            before = sum(all_peptides[pep].values())
            penalty = ((total_families-current_families+1)/total_families)**2
            for g,score in groups.items():
                corrected_score = score * penalty
                all_peptides[pep][g] = corrected_score
                if corrected_score <= args.final_cut**2:
                    if pep not in banned:
                        banned[pep] = {}
                    banned[pep][g] = corrected_score
            
        ####################################
        #If the down ranking have caused some peptides below the lower limit, remove them:
        args.compress = True 
        if args.compress:
            temp_peptides = {}
            for pep,v in all_peptides.items():
                groups = dict((group,score) for group,score in v.items() if pep not in banned or group not in banned[pep]) 
                temp_peptides[pep] = v
            all_peptides = temp_peptides
                        
        ####################################
        #Stats for CUPPlibrary
        print (("### The merged CUPPlibrary is being saved (%s) version: %s" % (args.setting,args.version)) if not args.clustering else "###",flush=True,end="")
        included = dict((os.path.basename(f).replace("_CUPPpool.json",""),1) for f in CUPPpool_path_list)
        stats = {"Version":args.version,"Created":int(psutil.Process(os.getpid()).create_time()),"Included":included,"Experimentals":experimentals}
        handle_out = open(args.compiled_json,'w')
        json.dump({"peptides":all_peptides,"meta":all_meta,"stats":stats,"meta_categories":meta_cat+["Function"]},handle_out,indent=3)
        handle_out.close()
        time.sleep(1)
        print (" DONE!")
       
    ####################################
    ###  CUPPprediction starts here  ###
    ####################################
    if args.predict or args.clustering:

        ####################################
        #Initiation of prediction and view of user arguments:
        show_arguments(args,predict=True)
        
        ####################################
        #Load the peptides for prediction:
        if os.path.exists(args.compiled_json):
            print ("### The precompiled CUPPlibrary is being loaded: %s" % args.setting,flush=True,end="")
            with open(args.compiled_json) as json_data:
                powder_merged = json.load(json_data)
            all_peptides = powder_merged["peptides"]; CUPPmeta = powder_merged["meta"]; meta_cat = powder_merged["meta_categories"]
            print ("_%s DONE" % powder_merged["stats"]["Version"])
            
            ####################################
            #Report if some models are missing in the library:
            missing = set(args.common) - (set(args.common) & set(powder_merged["stats"]["Included"]))
            for js in CUPPpool_path_list:
                family = "_".join(os.path.basename(js).split("_")[:-1])
                if family in missing:
                    print (' !  The powder of "%s" are not found in the precompiled CUPPlibrary (%s)' % (family,args.setting))
        else:
            sys.exit(" !  You need to compile a CUPPlibrary for the current settings %s" % args.setting)
        if not all_peptides:
            sys.exit(" !  The selected CUPPlibrary contained no peptides: %s" % args.compiled_json)
        
        ####################################
        #Make the meta data lists into a dict of dicts:
        CUPPmeta = create_table(CUPPmeta,targets=args.common,meta_cat=meta_cat)
        print ("### Total number of CUPP groups in CUPPlibrary: %s" % len(CUPPmeta))
        print ("### Total number of peptides in CUPPlibrary: %s" % len(all_peptides))
        
        ####################################
        #To be included during recompiling instead: 
        for cupp,v in CUPPmeta.items():
            if "Average" not in CUPPmeta[cupp]["Freq_sum"]:
                CUPPmeta[cupp]["Freq_sum"]["Average"] = 201
            else:
                CUPPmeta[cupp]["Freq_sum"]["Average"] = 201 if 200 < v["Freq_sum"]["Average"] else v["Freq_sum"]["Average"]+1 if 40 < v["Freq_sum"]["Average"] else 41

        
        ####################################
        #Determine peptide uniqueness:
        if 0:
            view = {} 
            dual = dict((k,0) for k in CUPPmeta)
            for pep,v in all_peptides.items():
                for cupp,s in v.items():
                    if 0 < s:
                        if cupp not in view:
                            view[cupp] = [0,0]
                        view[cupp][0 if len(v) == 1 else 1] += 1
                        
            for cupp,l in sorted(view.items(),key=lambda x: family_sort(x[0])):
                print ("%s\t%s\t%s\t%s\t%s" % (cupp,round(100*l[0]/sum(l),1),"%",view[cupp][0],view[cupp][1]))
            
            sys.exit()

        ####################################
        #Obtaining the collection of sequences:
        collection, fasta_meta = obtain_collection(original_list if not args.domains_only else domain_list, sep = args.sep, common = args.general_name, meta_cat = args.header,beta=args.beta_options)
        locate_accession = origin(original_list,cdhit_list,CUPPmeta,sorted(args.common),args.header,sep=args.sep)

        ####################################
        #Identify proteins not directly included in the clustering but represented by another protein in a high similarity group 
        for acc in fasta_meta:
            fasta_meta[acc]["cupp"] = locate_accession[acc][0] if acc in locate_accession else dict()

        ####################################
        #Define user options relevant for prediction:
        if "Experimentals" not in powder_merged["stats"]:
            args.exclude_family.update({"GH0","CE0","GT0","AA0","PL0","PL28","GH145","GH146","GH147","GH148","GH149","GH150","GH151","GH152","GH153","AA14","AA15","GT105"})
            args.exclude_family = args.exclude_family if "experimental_families" not in args.beta_options else {}
        opts = argparse.Namespace(all_peptides=all_peptides, meta=CUPPmeta, overlap=args.overlap, n_mer=args.n_mer, ambiguous=args.ambiguous,
            line_cut=args.position_cut, precision_ratio = args.best_hit,occurence = args.occurence,  type = args.type, cupp_minimum_score=args.cupp_minimum_score,
            minimum_cup = args.minimum_cup, evidence = args.evidence, domain_min = args.domain_min, relatedness = args.clustering, reduce_peptides = False,
             exclude_family = args.exclude_family, common = args.general_name, beta_options = args.beta_options, keep_only = args.keep_only)
                    
        ####################################
        #Determine relatedness with altered parameters:
        if args.clustering:
            opts = argparse.Namespace(all_peptides=all_peptides, meta=CUPPmeta, n_mer=args.n_mer, ambiguous=args.ambiguous, overlap=1,
                line_cut=0.4, precision_ratio = 0.5,occurence = False,  type = "none", keep_only = {},
                minimum_cup = 5, evidence = 0.01, domain_min = 20, relatedness = args.clustering, cupp_minimum_score=1,
                exclude_family = {}, common = args.general_name, beta_options = args.beta_options, reduce_peptides = False)
        
        ####################################
        #Predict the collection using one or several cores:
        predictions = {}; start = time.time()
        if args.jobs == 1:
        
            ####################################
            #Predict the entries one by one:
            for seq,acc in collection.items():
                for a in fasta_meta[acc]["Accession"]:
                    predictions[a] = predict(seq,opts)
        else:
            ####################################
            #Split the collection into smaller parts for parallel processing:
            parts = dict((n,{}) for n in range(args.jobs)); counter = 0
            for seq,acc in collection.items():
                parts[int(round(args.jobs * (counter/args.jobs-int(counter/args.jobs))))][seq] = acc
                counter += 1
                
            ####################################
            #Process  the collection using several processors:
            with Pool(processes=args.jobs) as pool:
                result_part = pool.map(parallel_predict, ([parts[q],opts] for q in range(args.jobs)))
                
            ####################################
            #Merge parts into single dict: 
            for q in result_part:
                predictions.update(q)
        print ("### Prediction COMPLETED: %s [s]" % (round(time.time() - start,3)))
        
        ####################################
        #Obtain summary of the predictions:
        eva = dict((tag,[0,0,0]) for tag in ["cupp","raw_function","Family","Subfam"])
        handle_fas = open(args.fasta_path,"w") if opts.type in ["both","fasta","quick"] else False
        handle_raw = open(args.raw_path,"w") if opts.type in ["both","raw"] else False
        itol_index = {}; custom_colors = {}

        for acc,result in predictions.items():
        
            ####################################
            #Output the result:
            seq = result["additional"]["seq"].replace("_","x") if "seq" in result["additional"] else "X" if "keep" in args.beta_options else ""
            out(handle_fas,handle_raw,acc,result,seq,args.beta_options,CUPPmeta,sep=args.sep)
        
            ####################################
            #Save den dendrogram as Newick format for handling e.g. in iTOL:
            if args.itol_out and not args.clustering:
                for fam,v in result["cupp"].items():
                    for cupp,score in v.items():
                        if cupp not in itol_index:
                            itol_index[cupp] = [{}]
                        itol_index[cupp][0][acc] = result["additional"][cupp][0] if cupp in result["additional"] else ""
    
            ####################################
            #For each category, determine stats:
            for tag in eva:
            
                ####################################
                #The known meta data of protein:
                if tag not in fasta_meta[acc] or not fasta_meta[acc][tag]:
                    continue
                eva[tag][2] += 1
                
                ####################################
                #Handle dual functions in query:
                parts = set()
                for obj in fasta_meta[acc][tag]:
                    obj = obj.split(":",1)[-1].replace("+","&").replace("-","&")
                    for p in obj.split("&"):
                        if tag != "cupp" or p.split(":")[-1][0] != "0":
                            parts.add(p)
                rec = parts
                
                ####################################
                #Handle dual functions in groups:
                parts = set()
                for fam,v in result[tag].items():
                    for obj in v:
                        obj = obj.split(":",1)[-1].replace("+","&").replace("-","&")
                        for p in obj.split("&"):
                            if tag != "cupp" or p.split(":")[-1][0] != "0":
                                parts.add(p)
                res = parts
                    
                ####################################
                #Identify incorrect predictions:
                score_preci = 1/(1+len(res - rec))
                eva[tag][1] += score_preci
                        
                ####################################
                #Identify correct predictions:
                score_sensi = len(rec & res)/len(rec)
                eva[tag][0] += score_sensi
                
                ####################################
                #Debug:
                #inspect = True if "&" in "".join(list(rec)) else False
                if "prediction_fails" in args.beta_options:
                    if score_preci != 1:
                        print ("-->",tag,"preci",acc,rec,result[tag],score_preci)
                    if score_sensi != 1:
                        print ("-->",tag,"sensi",acc,rec,result[tag],score_sensi)
        handle_fas.close() if handle_fas else None
        handle_raw.close() if handle_raw else None
        
        ####################################
        #Output labels for itol:
        if not args.clustering:
            itol(args.itol_out,itol_index,args.general_name+"_predict",args.header,fasta_meta,path=args.working_dir)
        
        ####################################
        #Output the overall stats:
        if eva["Family"][2]:
            if not args.clustering:
                for tag,report in sorted(eva.items()):
                    print ("### Validation: %s    \t%s (Sensitivity)\t  %s (Precision)\t%s Total" % (tag.replace("raw_f","F"),round(100*report[0]/report[2],2) if report[2] else "N/A",round(100*report[1]/report[2],2) if report[2] else "N/A",report[2]))
        else:
            sys.stderr.write("!!! No hits found during prediction (with family info)\n!")
        
        ####################################
        #Determine the connectivity among the groups by co-presence:
        if args.clustering:
            dendrogram_path = "%s/graphs%s/%s_%s_related.json" % (args.working_dir,args.database_version,args.general_name,args.setting)
            determine_relatedness(predictions,CUPPmeta,CUPPpool_path,dendrogram_path,rel_cut=args.rel_cut)
        
    ####################################
    #View the saved plot files:
    if not args.pipe and args.dendrogram:

        ####################################
        #The saved linkages can be inspected:
        print ("### Loading plots for: %s" % args.setting)
        plot_files_list = []
        plot_files_list.append("%s_%s_%s_acc.json" % (args.general_name,args.setting,len(args.gradient)-1))
        #plot_files_list.append("%s_%s_%s_cup.json" % (args.general_name,args.setting,len(args.gradient)-1))
        plot_files_list.append("%s_%s_related.json" % (args.general_name,args.setting))
        for graph_fil in plot_files_list+["test.json"]:
            path = os.path.join(args.working_dir,"graphs%s" % (args.database_version),graph_fil)
            if os.path.exists(path):
           
                with open(path) as json_data:
                    g = json.load(json_data)
                    
                dendro(np.array([np.array(e) for e in g["links"]]),g["labels"],algo=g["algo"],cc=g["cc"],name=g["name"],fil="",horisontal=g["horisontal"],note=g["note"] if "note" in g else "")
            elif graph_fil != "test.json":
                sys.stderr.write(" !  Could not find graph: %s\n" % (path))
        import matplotlib.pyplot as plt
        plt.show()
        
    ####################################
    #Output the complete run-time for successful clustering or recompile:
    if CUPPpool_path_list:  
        print ("### Full run-time:\t%s [s]" % str(round(time.time() - psutil.Process(os.getpid()).create_time(),3)))
