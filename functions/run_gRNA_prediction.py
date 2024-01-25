
import sys
import cPickle
import subprocess
from optparse import OptionParser

sys.path.insert(0, '/data/share/htp/perturb-seq/TF_selection_gRNA_design/Horlbeck_design/ScreenProcessing/')

execfile('/data/share/htp/perturb-seq/gRNA_design_workflow/functions/sgRNA_learning.py')

parser = OptionParser()
parser.add_option("-t", "--tss", dest="tss_file",  help="tss file", metavar="FILE")
#parser.add_option("-p", "--p1p2",dest="p1p2_file", help="p1p2 file", metavar="FILE" )
parser.add_option("-g", "--genome",dest="dict_file", help="picard tools dict file", metavar="FILE" )
parser.add_option("-a", "--atac",dest="atac_file", help="open chromatin bw file", metavar="FILE" )
parser.add_option("-m", "--model",dest="model_file", help="pickle with trained model", metavar="FILE" )
parser.add_option("--btf", dest="bt_folder", help="folder for all the bowtiefiles", metavar="FILE" )
parser.add_option("--btp", dest="bt_index_prom", help="bowtie indexed promoters", metavar="FILE" )
parser.add_option("--btg", dest="bt_index_genome", help="bowtie indexed genome", metavar="FILE" )
parser.add_option("-o", "--outbase", dest="outbase", help="name of csv files to write", metavar="FILE" )
parser.add_option("-n", "--ngRNAs", dest="ngRNAs", help="number of gRNAs to pick", type="int", default=40 )


(options, args) = parser.parse_args()


#reading all files
genomeDict = loadGenomeAsDict(options.dict_file) 
print("Genome loading done!")
tss = pd.read_csv(options.tss_file,sep='\t', index_col=range(2), header = 0, dtype = {"chromosome": "string"})
tss['primary TSS'] = tss['primary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]), int(tupString.strip('()').split(', ')[1].split('.')[0])))

#p1p2 = pd.read_csv(options.p1p2_file,sep='\t', header=0, index_col=range(2), dtype = {"chromosome": "string"})
#p1p2['primary TSS'] = p1p2['primary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]), int(tupString.strip('()').split(', ')[1].split('.')[0])))
#p1p2['secondary TSS'] = p1p2['secondary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]),int(tupString.strip('()').split(', ')[1].split('.')[0])))

bwhandleDict = {'atac':BigWigFile(open(options.atac_file))}

transformedParams_train ,fitTable, estimators, scaler, reg = pd.read_pickle(options.model_file)

#get initial gRNAs
#libraryTable, sgInfo = findAllGuides2(p1p2, genomeDict, (-25,500)) 
libraryTable, sgInfo = findAllGuides2(tss, genomeDict, (-25,500)) #changed
print("gRNA candidate selection done!")
#prepare predictors
paramTable = generateParamTableATAC(libraryTable, sgInfo , tss, genomeDict, bwhandleDict)
#paramTable.to_csv(options.outbase + "paramTable" + ".csv" ,sep="\t")
print("Predictor variables extracted.")


# scoring
transformedParams_new = transformParams(paramTable, fitTable, estimators)
colTups = []
for (l1, l2), col in transformedParams_new.iteritems():
    colTups.append((l1,str(l2)))
transformedParams_new.columns = pd.MultiIndex.from_tuples(colTups)
predictedScores_new = pd.Series(reg.predict(scaler.transform(transformedParams_new.loc[:, transformedParams_train.columns].fillna(0).values)), index=transformedParams_new.index)

print("scoring done.")

# off-targets
if not os.path.exists(options.bt_folder):
    os.mkdir( options.bt_folder )
    
fqFile=options.bt_folder+'/gRNAs.fq'
outputTempBowtieFastq(libraryTable, fqFile)

#specifying a list of parameters to run bowtie with
#each tuple contains
# *the mismatch threshold below which a site is considered a potential off-target (higher is more stringent)
# *the number of sites allowed (1 is minimum since each sgRNA should have one true site in genome)
# *the genome index against which to align the sgRNA sequences; these can be custom built to only consider sites near TSSs
# *a name for the bowtie run to create appropriately named output files
alignmentList = [(39,1, options.bt_index_prom,'39_nearTSS'),
                ( 31,1, options.bt_index_prom,'31_nearTSS'),
                ( 21,1, options.bt_index_genome,'21_genome'),
                ( 31,2, options.bt_index_prom,'31_2_nearTSS'),
                ( 31,3, options.bt_index_prom,'31_3_nearTSS')]

alignmentColumns = []
for btThreshold, mflag, bowtieIndex, runname in alignmentList:

    alignedFile = options.bt_folder + '/' + runname + '_aligned.txt'
    unalignedFile = options.bt_folder + '/' + runname + '_unaligned.fq'
    maxFile = options.bt_folder + '/' + runname + '_max.fq'
    
    bowtieString = 'bowtie -n 3 -l 15 -e '+str(btThreshold)+' -m ' + str(mflag) + ' --nomaqround -a --tryhard -p 16 --chunkmbs 256 ' + bowtieIndex + ' --suppress 5,6,7 --un ' + unalignedFile + ' --max ' + maxFile + ' '+ ' -q '+fqFile+' '+ alignedFile
    print (bowtieString)
    print (subprocess.call(bowtieString, shell=True)) #0 means finished without errors

    #parse through the file of sgRNAs that exceeded "m", the maximum allowable alignments, and mark "True" any that are found
    try:
        with open(maxFile) as infile:
            sgsAligning = set()
            for i, line in enumerate(infile):
                if i%4 == 0: #id line
                    sgsAligning.add(line.strip()[1:])
    except IOError: #no sgRNAs exceeded m, so no maxFile created
        sgsAligning = set()
                    
    alignmentColumns.append(libraryTable.apply(lambda row: row.name in sgsAligning, axis=1))
    
#collate results into a table, and flip the boolean values to yield the sgRNAs that passed filter as True
alignmentTable = pd.concat(alignmentColumns,axis=1, keys=zip(*alignmentList)[3]).ne(True)

print("Bowtie search done.")

predictedScores_new.name = 'predicted score'
v2Table = pd.concat((libraryTable, predictedScores_new, alignmentTable, sgInfo),sort=True, 
                    axis=1, keys=['library table v2', 'predicted score', 'off-target filters', 'sgRNA info'])
#v2Table.to_csv("new_data_files/mf6_all_gRNAs_mergedATAC.txt",sep="\t")


#minimum overlap between two sgRNAs targeting the same TSS
nonoverlapMin = 3

#number of sgRNAs to pick per gene/TSS
sgRNAsToPick=options.ngRNAs

#list of off-target filter (or combinations of filters) levels, matching the names in the alignment table above
offTargetLevels = [['31_nearTSS', '21_genome'],
                  ['31_nearTSS'],
                  ['21_genome'],
                  ['31_2_nearTSS'],
                  ['31_3_nearTSS']]

#for each gene/TSS, go through each sgRNA in descending order of predicted score
#if an sgRNA passes the restriction site, overlap, and off-target filters, accept it into the library
#if the number of sgRNAs accepted is less than sgRNAsToPick, reduce off-target stringency by one and continue
v2Groups = v2Table.groupby([('library table v2','gene'),('library table v2','transcripts')])
newSgIds = []
unfinishedTss = []
v2Table.to_csv(options.outbase + "_unfiltered.csv" ,sep="\t")


for (gene, transcript), group in v2Groups:
    geneSgIds = []
    geneLeftPositions = []
    empiricalSgIds = dict()
    
    stringency = 0
    
    while len(geneSgIds) < sgRNAsToPick and stringency < len(offTargetLevels):
        for sgId_v2, row in group.sort_values(('predicted score','predicted score'), ascending=False).iterrows():
            leftPos = row[('sgRNA info', 'position')] - (23 if row[('sgRNA info', 'strand')] == '-' else 0)
            if len(geneSgIds) < sgRNAsToPick and row['off-target filters'].loc[offTargetLevels[stringency]].all() \
                and checkOverlaps(leftPos, geneLeftPositions, nonoverlapMin):
                geneSgIds.append((sgId_v2,
                                  row[('sgRNA info', 'chromosome')],
                                  gene,transcript,
                                  row[('sgRNA info', 'position')],
                                  row[('sgRNA info', 'strand')],
                                  row[('library table v2','sequence')], 
                                  row[('library table v2','genomic sequence')], 
                                  row[('predicted score','predicted score')],
                                 stringency))
                geneLeftPositions.append(leftPos)
                
        stringency += 1
        
    newSgIds.extend(geneSgIds)
        
libraryTable_complete = pd.DataFrame(newSgIds, columns = ['sgID', 'chromosome', 'gene', 'transcript','position','strand', 'gRNA_sequence', 'genomic_sequence',
 'predicted_score', 'off_target_stringency']).set_index('sgID')

libraryTable_complete.to_csv(options.outbase + "_top" + str(sgRNAsToPick) + ".csv" ,sep="\t")

sys.exit("gRNA selection and scoring completed")
