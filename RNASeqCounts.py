#!/usr/bin/python

import sys
import glob
import os
import string
import pysam
from collections import defaultdict

referenceIDToCounts = {}

#Read is Paired
#Read Mapped in Proper Pair
#Mate Reverse Strand
# R1
#---->
#--------------------------
#					<----
#					  R2

def isSenseAlignment(record): 
	return True if (record.is_paired and record.is_proper_pair and record.mate_is_reverse) else False

#Read is Paired
#Read Mapped in Proper Pair
#Read Reverse Strand
#					  R1
#					<----
#--------------------------
#---->
# R2
def isAntiSenseAlignment(record): 
	return True if (record.is_paired and record.is_proper_pair and record.is_reverse) else False

#Used because when STAR aligner says there is a reverse match, it provides the sequence that maps to the transcript.
# However, this sequence is NOT the original read. The original read is the reverse complementary sequence. Therefore,
# If it is a sense alignment, R2 will be flipped, and antisense, R1 will be flipped. We use this function to flip
# R1 and R2 back so we can keep a standard way of viewing all our reads to determine how many of the SAME reads we have
# during the counting step.
def reverseComplementaryStrand(dna):
	return dna.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))[::-1]

#Process Endogenous Alignments
def assignRead(read, ReferenceIDDictionary, ReferenceIDCount, ReadAlignments, outputDir):
	OverallReads = defaultdict(list)
	keptAlignments = list()
	senseAlignment=False
	noUsableRecords=True
	for record in read:
		if record.is_secondary is False:
			if record.is_proper_pair:
				if record.is_read1:	
					if isSenseAlignment(record): 			#Either the specific record is a Sense Alignment
						noUsableRecords=False
						senseAlignment=True
						thisLib = "Sense"
						#Find the associated R2 to this R1
						read_one_reference = record.reference_name
						for look_for_read2 in read:
							if look_for_read2.is_read2 and look_for_read2.reference_name == read_one_reference:
								read_two_record = look_for_read2
								break
					elif isAntiSenseAlignment(record):		#Or the record is an Anti-Sense Alignment
						noUsableRecords=False
						antiSenseAlignment=True
						thisLib = "AntiSense"
						#Find the associated R2 to this R1
						read_one_reference = record.reference_name
						for look_for_read2 in read:
							if look_for_read2.is_read2 and look_for_read2.reference_name == read_one_reference:
								read_two_record = look_for_read2
								break

					#Group the Reads by their Sense or Antisense alignments
					#OverallReads[thisLib].append(record)

					#Our Library detects appears to be Reversely Stranded so we should just keep them all...
					keptAlignments.append(record)

					#print record.query_name
					#print reference
				elif record.is_read2:
					continue

	for record in keptAlignments:
		print record.reference_name
	print "HELLO"

	if noUsableRecords:
		return ReferenceIDDictionary, ReferenceIDCount, ReadAlignments	#If it isn't a proper pair, that means only one Read mapped, we need both

	#We prioritize Sense Alignments. So unless there are no sense alignments, it will usually be Sense
	# if (senseAlignment):
	# 	keptAlignments = OverallReads["Sense"]
	# else:
	# 	keptAlignments = OverallReads["AntiSense"]

	counter=0
	referenceID = list()
	readName = ""
	for record in keptAlignments:
		counter+=1 							#Used for outputting the single Read Information
		mapsTo = record.reference_name		#For the this specific read, we want to know what the records inside
											# map to. So if a read has 3 records, then it will map to 3 things
		if counter is 1:					#This stores the read information
			# This file is really big. We are opting not to output it anymore.
			# with open(os.path.join(outputDir+"ReadAlignments.dict"), "a") as readAlignments:
			if senseAlignment:
				ReadAlignments[record.query_name].append(record.query_sequence + "---" + reverseComplementaryStrand(read_two_record.query_sequence))
				ReadAlignments[record.query_name].append("Sense")
			# 		readAlignments.write(record.query_name + "\t" + record.query_sequence + "---" + reverseComplementaryStrand(read_two_record.query_sequence) + "\tSense")
			else:
				ReadAlignments[record.query_name].append(reverseComplementaryStrand(record.query_sequence) + "---" + read_two_record.query_sequence)
				ReadAlignments[record.query_name].append("AntiSense")
			# 		readAlignments.write(record.query_name + "\t" + reverseComplementaryStrand(record.query_sequence) + "---" + read_two_record.query_sequence + "\tAntiSense")
			readName = record.query_name

		if mapsTo not in referenceID:		#Combines all the different genes that the this specific read maps to
			referenceID.append(mapsTo)

	referenceIDs=""
	first=True
	for ID in referenceID: 						#Creates the individual references, so basically... if a read maps
		if first:								# to gene1, gene2, and gene3, then the first reference (1) would be
			referenceIDs = referenceIDs + ID 	# 1	gene1 | gene2 | gene3		Whereas another read that only maps
		else:									# to 2 genes would be: 		2	gene1 | gene3
			referenceIDs = referenceIDs + "|" + ID
		first=False

	if referenceIDs not in ReferenceIDDictionary:				#This part of the program is used to consolidate the information
		ReferenceIDDictionary[referenceIDs]=ReferenceIDCount		# between the individual Reads and their corresponding reference
		ReadAlignments[readName].append(ReferenceIDCount)	# information. Each read will have an associated reference ID
															# that reference ID will be linked to genes that the reads
															# have been aligned to
		# with open(os.path.join(outputDir+"AligmentReference.dict"), "a") as alignmentReference:
		# 	alignmentReference.write(str(ReferenceIDCount) + "\t" + referenceIDs + "\n")
		# with open(os.path.join(outputDir+"ReadAlignments.dict"), "a") as readAlignments:
		# 	readAlignments.write("\t" + str(ReferenceIDCount) + "\n")
		ReferenceIDCount+=1	
	else:
		ReadAlignments[readName].append(ReferenceIDDictionary[referenceIDs])
		# with open(os.path.join(outputDir+"ReadAlignments.dict"), "a") as readAlignments:
		# 	readAlignments.write("\t" + str(ReferenceIDDictionary[referenceIDs]) + "\n")

	return ReferenceIDDictionary, ReferenceIDCount, ReadAlignments

	#Basically...
	#Read Information 							#Reference ID Information
	#Read1	SEQUENCE 	#ReferenceID 			#ReferenceID 			#List of Genes
	#Read2	SEQUENCE 	#ReferenceID 			#AnotherReferenceID		#List of Genes

	#In this case, Read 1 and Read 2 map to the same list of genes, so instead of storing the information twice.
	# We only need to store the reference one and just use the referenceID as the link between the GeneList and Reads
def addInsert(insert, AlignmentDictionary):
	global referenceIDToCounts 		#Global Gene Counts List
	referenceIDs = list()
	readIDs = list()
	SenseOrAnti = insert[0][1][1]

	#For the group of inserts, we look for the associated reference ID value and pull the reference IDs
	# Aka if the group of inserts match up with ReferenceIDValue or 2, it would pull the 2nd Group of IDs
	# from the dictionary
	for referenceID in AlignmentDictionary[insert[0][1][2]].split("|"):
		referenceIDs.append(referenceID)
	for readID in insert:
		readIDs.append(readID[0])
	#We also want to know which Reads were associated with that specific group of inserts, aka which reads
	# have the same sequence

	for referenceID in referenceIDs:
		if SenseOrAnti not in referenceIDToCounts:
			referenceIDToCounts[SenseOrAnti] = {}

		if referenceID not in referenceIDToCounts[SenseOrAnti]:
			alignmentCount = [0,0,0,0]
		else:
			alignmentCount = referenceIDToCounts[SenseOrAnti][referenceID]

		#Unique Read - If there were 4 reads in this insert group, it would still only count it as 1 sequence for 1 Gene
		alignmentCount[0] += 1
		#TotalReadCount - If there were 4 reads with the same sequence, then this would add 4 to the Gene Count
		alignmentCount[1] += len(readIDs)
		#MultiMappingReadCount - If there were 4 reads, and they mapped to 3 genes, this would be 4/3 counts per gene
		alignmentCount[2] += 1.*len(readIDs)/1.*len(referenceIDs)
		#UniqueMappingReads - Taha Equation, where we only consider a count if there were no multi-mapping reads
		if (len(referenceIDs) == 1):
			alignmentCount[3] += len(readIDs)
		else:
			alignmentCount[3] += 0

		referenceIDToCounts[SenseOrAnti][referenceID] = alignmentCount

#Quantify Alignments
def readAndCountInserts(ReferenceIDDictionary, ReadAlignments):
	previousInsert = None
	thisInsert = list()
	#Converts ReferenceIDDictionary from:
	# Gene 1 | Gene 2 | Gene 3 		1
	# Gene 1 | Gene 5 | 			2

	# To:
	# 1		Gene 1 | Gene 2 | Gene 3
	# 2 	Gene 1 | Gene 5
	AlignmentDictionary=dict((v,k) for k,v in ReferenceIDDictionary.iteritems())
	#Sort alphabetically by the sequence in order to identify potential reads with repeating sequences
	#print sorted(ReadAlignments.iteritems(), key=lambda x: x[1])
	for read in sorted(ReadAlignments.iteritems(), key=lambda x: x[1]):
		#Read 			= Read[0]
		#Read Sequence 	= Read[1][0]
		#Sense or Anti  = Read[1][1]
		#Read Reference = Read[1][2]
		currentInsert = read[1][0]
		#There can be Reads with repetitive Sequences, we want to find all of them
		if (previousInsert == None or currentInsert == previousInsert):
			thisInsert.append(read)
			previousInsert = currentInsert
		else: #If it's not the same as the previous, then this is the whole insert
			addInsert(thisInsert, AlignmentDictionary)
			thisInsert = list()
			thisInsert.append(read)
			previousInsert = currentInsert
	addInsert(thisInsert, AlignmentDictionary)

def writeCounts(outputDir):
	global referenceIDToCounts
	for SenseOrAnti in referenceIDToCounts.keys():
		with open(os.path.join(outputDir+"ReadCounts" + SenseOrAnti + ".txt"), "a") as countsFile:
			countsFile.write("ReferenceID\tUniqueFragmentCount\tTotalReadCount\tMultimapAdjustedReadCount\tUniqueMappingReads\n")
			for referenceID in sorted(referenceIDToCounts[SenseOrAnti].iteritems(), key=lambda x:x[1], reverse = True):
				countsFile.write(referenceID[0]+"\t"+str(referenceID[1][0])+"\t"+str(referenceID[1][1])+"\t"+str(referenceID[1][2])+"\t"+str(referenceID[1][3]) + "\n")


if __name__ == "__main__":
	SAM_FILE = sys.argv[1]
	outputDir = sys.argv[2]
	samfile = pysam.AlignmentFile(SAM_FILE, "rb")		#Reads SAM or BAM file. Has to be sorted using samtools.
	lastReadID = None									#In order to know when we can assign a Read, we need to know when the Read changes
	thisRead = {}										#Dictionary for Read
														#ReadID :=	Alignment_1
														#			Alignment_2
														#			Alignment_3
	ReferenceIDDictionary = {}
	ReadAlignments = defaultdict(list)
	ReferenceIDCount = 1
	#Reading in the Sorted SAM file generated from STAR Aligner
	for record in samfile.fetch():
		#A read can align to multiple different references, each alignment is called a "record"
		#Once all of the read alignments have been processed, then we can assign the Read proper counts
		if (record.query_name != lastReadID and lastReadID != None):
			#Assign the Read
			ReferenceIDDictionary, ReferenceIDCount, ReadAlignments = assignRead(thisRead, ReferenceIDDictionary, ReferenceIDCount, ReadAlignments, outputDir)
			thisRead = {}
		thisRead[record]=record.reference_name 		#Add each alignment to the Read
		lastReadID = record.query_name				#To check when the Read changes

	ReferenceIDDictionary, ReferenceIDCount, ReadAlignments = assignRead(thisRead, ReferenceIDDictionary, ReferenceIDCount, ReadAlignments, outputDir)
	#For Debugging
	#print ReferenceIDDictionary
	#print ReadAlignments
	readAndCountInserts(ReferenceIDDictionary, ReadAlignments)
	writeCounts(outputDir)
	
	#alignmentDictionary=readAlignmentReference() 	#Global Reference(?)
	#readDictionary=readAlignments()
		#print record.query_name
		#print record.reference_name

	samfile.close()