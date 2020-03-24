/*
 * Code for phasing SV calls from their supporting reads and phased small variants along those reads.
 * After phasing these calls, it integrates them into the one large VCF file alongside the small variant calls.
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;

public class PhaseSVs
{
	// The factor by which the votes for one haplotype must be higher than the other to change a homozygous call to heterozygous
	static int OVERRULE_HOMOZYGOUS_FACTOR = 5;
	
	// The minimum read support needed for a homozygous variant to switch to heterozygous
	static int OVERRULE_HOMOZYGOUS_MINREADS = 5;
	
	// Whether or not to phase SVs
	static boolean PHASE_SVS = true;
	
	// Not currently used
	static double SVTOOLONG = 1.2;
	static double SVTOOSHORT = 0.8;
	static boolean REPORT_INVALID_LEN_INSERTIONS = true;
	
	// The number of phased SNVs a read must span for it to be assigned a haplotype
	static int MIN_PHASED_SNV = 5;
	
	// Maximum length of an SV to include
	static int MAX_SV_LEN = 100000;
	
	// Usage message if running the program incorrectly
	static String USAGE = "splicephase.pl phased.vcf "
			+ "sniffles.vcf loadreads.hairs spliced.vcf ref.fa\n";
	
	// The path to the samtools executable
	static String SAMTOOLS_PATH = "samtools";
	
	public static void main(String[] args) throws Exception
	{
		if(args.length < 5)
		{
			System.out.println(USAGE);
			System.exit(1);
		}
		else
		{
			// Check that the input files all exist
			
			if(!new File(args[0]).exists())
			{
				System.out.println("Invalid phased SNP file: " + args[0]);
				System.exit(1);
			}
			
			if(!new File(args[1]).exists())
			{
				System.out.println("Invalid Sniffles call file: " + args[1]);
				System.exit(1);
			}
			
			if(!new File(args[2]).exists())
			{
				System.out.println("Invalid HapCUT2 fragment file: " + args[1]);
				System.exit(1);
			}
			
			if(!new File(args[4]).exists())
			{
				System.out.println("Invalid reference file: " + args[4]);
				System.exit(1);
			}
			
			runSplicing(args[0], args[1], args[2], args[3], args[4]);
		}
	}
	
	/*
	 * Reverse complement of a String
	 */
	static char[] froms = new char[] {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
	static char[] tos = new char[] {'T', 'G', 'C', 'A', 't', 'g', 'c', 'a'};
	static String revComp(String s)
	{
		int n = s.length();
		char[] res = new char[n];
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(i);
			char to = c;
			for(int j = 0; j<froms.length; j++)
			{
				if(c == froms[j])
				{
					to = tos[j];
					break;
				}
			}
			res[n-1-i] = to;
		}
		return new String(res);
	}
	
	/*
	 * Queries a genomic substring - runs samtools faidx <genomeFile> chr:startPos-endPos
	 */
	static String genomeSubstring(String refFn, String chr, long startPos, long endPos) throws Exception
	{
		if(startPos > endPos)
		{
			return "";
		}
		String faidxCommand = String.format("%s faidx %s %s:%d-%d", SAMTOOLS_PATH, refFn, chr, startPos, endPos);
		Process child = Runtime.getRuntime().exec(faidxCommand);
		InputStream seqStream = child.getInputStream();
		Scanner seqInput = new Scanner(seqStream);
		
		// Make sure it produced an actual output
		if(!seqInput.hasNext())
		{
			seqInput.close();
			throw new Exception("samtools faidx did not produce an output: " + faidxCommand);
		}
		// Read in and ignore sequence name
		seqInput.next();
		
		// Make sure there's a sequence
		if(!seqInput.hasNext())
		{
			seqInput.close();
			throw new Exception("samtools faidx produced a sequence name but not an actual sequence: " + faidxCommand);
		}
		
		// Concatenate all lines of the output sequence
		StringBuilder res = new StringBuilder("");
		while(seqInput.hasNext())
		{
			res.append(seqInput.next());
		}
		seqInput.close();

		return res.toString();
	}
	
	/*
	 * Gets the reference sequence for a variant using samtools
	 */
	static String getVariantSeq(String chr, long startPos, long svLen, String refFn) throws Exception
	{
		return genomeSubstring(refFn, chr, startPos, startPos + svLen);
	}
	
	static void runSplicing(String snpsFn, String svsFn, 
			String hairsFn, String vcfOfn, String refFn) throws Exception
	{
		// Generate filenames for phasing information files based on output VCF
		String readPhaseOfn = vcfOfn + ".readphase";
		String svPhaseOfn = vcfOfn + ".svphase";
		String svPhaseDetailsOfn = vcfOfn + ".svphase.details";
		
		// Initialize PrintWriters for output files
		PrintWriter out = new PrintWriter(new File(vcfOfn));
		PrintWriter readPhaseOut = new PrintWriter(new File(readPhaseOfn));
		PrintWriter svPhaseOut = new PrintWriter(new File(svPhaseOfn));
		PrintWriter svPhaseDetailsOut = new PrintWriter(new File(svPhaseDetailsOfn));
		
		// Initialize Scanners for input files
		Scanner snpsInput = new Scanner(new FileInputStream(new File(snpsFn)));
		Scanner svsInput = new Scanner(new FileInputStream(new File(svsFn)));
		Scanner hairsInput = new Scanner(new FileInputStream(new File(hairsFn)));
		
		// Read in the phased SNPs
		System.err.println("Reading phased SNPs");
		ArrayList<String> snpsHeader = new ArrayList<String>();
		
		// Mapping from (chromosome, position) pair to SNP
		TreeMap<String, TreeMap<Integer, Variant>> snpData = new TreeMap<String, TreeMap<Integer, Variant>>();
		
		// List of all SNPs
		ArrayList<Variant> snpList = new ArrayList<Variant>();
		
		// Number of SNPs seen so far
		int numSnps = 0;
		
		// Read SNPs one at a time
		while(snpsInput.hasNext())
		{
			String line = snpsInput.nextLine();
			if(line.startsWith("#"))
			{
				snpsHeader.add(line);
			}
			else
			{
				Variant var = new Variant(line);
				numSnps++;
				if(!snpData.containsKey(var.chromosome))
				{
					snpData.put(var.chromosome, new TreeMap<Integer, Variant>());
				}
				snpData.get(var.chromosome).put(var.pos, var);
				snpList.add(var);
			}
		}
		
		snpsInput.close();
		System.err.printf("Loaded SNP file: %d header lines and %d variants\n", snpsHeader.size(), numSnps);
		
		// Now load the sniffles calls
		System.err.println("Loading Sniffles calls");
		int svCount = 0;
		
		// Information to report about each type
		HashMap<String, Report> typeReports = new HashMap<String, Report>();
		HashMap<String, Read> readsToPhase = new HashMap<String, Read>();
		TreeMap<String, TreeMap<Integer, Variant>> svData = new TreeMap<String, TreeMap<Integer, Variant>>();
		
		// Read SVs one at a time
		while(svsInput.hasNext())
		{
			String line = svsInput.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			
			svCount++;
			
			Variant v = new Variant(line);
			v.fillSvFields();
			
			// Update counter for the type
			if(!typeReports.containsKey(v.type))
			{
				typeReports.put(v.type, new Report());
			}
			typeReports.get(v.type).total++;
			
			// Update counters for each read supporting this SV
			for(String s : v.readNames)
			{
				if(!readsToPhase.containsKey(s))
				{
					readsToPhase.put(s, new Read());
				}
				readsToPhase.get(s).count++;
			}
			
			// Update map of (chromosome, position) pair to SV
			if(!svData.containsKey(v.chromosome))
			{
				svData.put(v.chromosome, new TreeMap<Integer, Variant>());
			}
			svData.get(v.chromosome).put(v.pos, v);
		}
		
		svsInput.close();
		System.err.printf("Loaded Sniffles calls: %d variants involving %d reads\n", svCount, readsToPhase.size());
		for(String s : typeReports.keySet())
		{
			System.err.println(s + ": " + typeReports.get(s).total);
		}
		
		// Determine read phasing information
		System.err.println("Determining read phasing");
		int svHairs = 0;
		int hairsLines = 0;
		
		// Go over each line of the fragment file from HapCUT2
		while(hairsInput.hasNext())
		{
			String line = hairsInput.nextLine();
			hairsLines++;
			
			String[] tokens = line.split(" ");
			String readId = tokens[1];
			
			// Ignore fragments not involved in SVs
			if(!readsToPhase.containsKey(readId))
			{
				continue;
			}
			
			svHairs++;
			
			// Get the number of variant blocks involved in the read
			int numBlocks = Integer.parseInt(tokens[0]);
			for(int i = 0; i<numBlocks; i++)
			{
				// Total number of variants in this block supporting each haplotype
				int hap1Vars = 0, hap2Vars = 0;
				
				// The index (within the VCF file) of the first variant in the consecutive block
				int firstVariantId = Integer.parseInt(tokens[i*2 + 2]);
				
				// String of 0 and 1 for whether it contains the ALT allele for each variant in the block
				String alleleString = tokens[2*i+3];
				
				// Go over each variant in the block, and let them vote on the phase of this read
				for(int j = 0; j<alleleString.length(); j++)
				{
					// Load the appropriate SNP
					int variantId = firstVariantId + j - 1;
					Variant snp = snpList.get(variantId);
					char allele = alleleString.charAt(j);
					
					char hap = 'A';
					
					// Consider the genotype of the SNP and whether or not the read has the ALT allele
					if(snp.genotype.equals("0|1"))
					{
						// 0|1 and absent means first haplotype
						if(allele == '0')
						{
							hap1Vars++;
							hap = 'A';
						}
						
						//  0|1 and present means second haplotype
						else
						{
							hap2Vars++;
							hap = 'B';
						}
					}
					else if(snp.genotype.equals("1|0"))
					{
						//  1|0 and present means first haplotype
						if(allele == '1')
						{
							hap1Vars++;
							hap = 'A';
						}
						
						//  1|0 and absent means second haplotype
						else
						{
							hap2Vars++;
							hap = 'B';
						}
					}
					
					// Store information about this SNP/read pair for logging
					readsToPhase.get(readId).snps.append(" " + snp.chromosome + ":" + snp.pos + ":" + allele + ":" + hap);
				}
				
				// Add the votes from this block to the read's total
				readsToPhase.get(readId).hap1 += hap1Vars;
				readsToPhase.get(readId).hap2 += hap2Vars;
			}
		}
		hairsInput.close();
		
		// Log some information about read phasing
		readPhaseOut.println("#READID\tNUMSV\t|\tHAP1\tHAP2\t| HAP HAPR | SNPS\n");
		ArrayList<String> sortedReadNames = new ArrayList<String>();
		sortedReadNames.addAll(readsToPhase.keySet());
		Collections.sort(sortedReadNames);
		for(String s : sortedReadNames)
		{
			Read r = readsToPhase.get(s);
			double hap1Prop = (r.hap1 + r.hap2 > 0) ? 100.0 * r.hap1 / (r.hap1 + r.hap2) : 0;
			String hap = r.hap1 >= r.hap2 ? "hapA" : "hapB";
			readPhaseOut.printf("%s\t%d\t|\t%d\t%d\t| %s %7.02f  |%s\n",
					s, r.count, r.hap1, r.hap2, hap, hap1Prop, r.snps.toString());
		}
		readPhaseOut.close();
		
		System.err.printf("Scanned fragments: %d lines with %d involved in SVs\n", hairsLines, svHairs);
		
		// Process Sniffles SVs and phase them based on their reads
		svPhaseOut.println("chr:pos:genotype\ttype\tsvlen\tseqlen\t|\tnumreads\thap1\thap2\t| hap\thap1r\t|\tnewgenotype\tincludesv\toverrulehomo");
		svPhaseDetailsOut.println("chr:pos:genotype\ttype\tsvlen\tseqlen\t|\tnumreads\thap1\thap2\t| hap\thap1r\t|\tnewgenotype\tincludesv\toverrulehomo");
		int snifflesCount = 0, reportedSvs = 0, phasedSvs = 0, unphasedSvs = 0, svLenErr = 0, genotypeErr = 0;
		
		// Iterate over the SVs one chromosome at a time
		for(String chr : svData.keySet())
		{
			TreeMap<Integer, Variant> svChrData = svData.get(chr);
			
			// Loop over the positions with variants on this chromosome
			for(int pos : svChrData.keySet())
			{
				snifflesCount++;
				Variant v = svChrData.get(pos);
				
				// The number of votes for each haplotype based on all of the reads supporting the SV
				int hap1Votes = 0, hap2Votes = 0;
				
				// Add up all the votes
				for(String readName : v.readNames)
				{
					Read r = readsToPhase.get(readName);
					hap1Votes += r.hap1;
					hap2Votes += r.hap2;
				}
				
				// Get the proportion of votes which are for haplotype1 and phase the SV based on that
				double hap1Prop = (hap1Votes + hap2Votes > 0) ? 100.0 * hap1Votes / (hap1Votes + hap2Votes) : 0;
				String hap = hap1Votes >= hap2Votes ? "hapA" : "hapB";
				
				// If it used to be homozygous, set the haplotype value to "hom"
				String oldGenotype = v.genotype;
				if(oldGenotype.equals("1/1"))
				{
					hap = "hom";
				}
				
				System.err.printf("Analyzing %s:%d:%s\t%s\t%d\t|\t%d\t%d\t%d\t| %s\t%7.02f  \n", 
						v.chromosome, v.pos, oldGenotype, v.type, v.svlen, v.readNames.length, hap1Votes, hap2Votes, hap, hap1Prop);
				
				if(v.type.equals("INS") || v.type.equals("DEL") || v.type.equals("INV"))
				{
					if(Math.abs(v.svlen) <= MAX_SV_LEN)
					{
						// Ignore SVs that are already phased
						if(oldGenotype.indexOf('|') == -1)
						{
							String newGenotype = oldGenotype;
							boolean overruleHomozygous = false;
							boolean isPhased = false;
							
							// Handle homozygous variants by updating their genotype to 1|1
							if(oldGenotype.equals("1/1"))
							{
								isPhased = true;
								newGenotype = "1|1";
								
								// Possibly overrule this genotype if the haplotype votes are one-sided enough
								if(v.readNames.length >= OVERRULE_HOMOZYGOUS_MINREADS)
								{
									if(hap2Votes >= hap1Votes * OVERRULE_HOMOZYGOUS_FACTOR)
									{
										newGenotype = "0|1";
										overruleHomozygous = true;
									}
									else if(hap1Votes >= hap2Votes * OVERRULE_HOMOZYGOUS_FACTOR)
									{
										newGenotype = "1|0";
										overruleHomozygous = true;
									}
								}
							}
							else
							{
								if(hap1Votes + hap2Votes < MIN_PHASED_SNV)
								{
									// Too few small variant calls to get a confident genotype, so leave unphased
									isPhased = false;
									if(hap.equals("hapA"))
									{
										newGenotype = "1/0";
									}
									else if(hap.equals("hapB"))
									{
										newGenotype = "0/1";
									}
								}
								else
								{
									// Phase according to which haplotype won the vote
									isPhased = true;
									if(hap.equals("hapA"))
									{
										newGenotype = "1|0";
									}
									else if(hap.equals("hapB"))
									{
										newGenotype = "0|1";
									}
								}
							}
								
							// This is the sequence length of the SV and is used for logging
							int seqLength = 0;
							if(v.seq != null)
							{
								seqLength = v.seq.length();
							}
							
							for(String readName : v.readNames)
							{
								svPhaseDetailsOut.println("== " + readName);
							}
							
							// Whether or not to include this SV - can be set to false based on certain criteria
							boolean includeSv = true;
							
							if(v.type.equals("DEL"))
							{
							}
							else if(v.type.equals("INS"))
							{
								// There was some logic about comparing the insertion sequence length and the SVLEN field that isn't working right
							}
							else if(v.type.equals("INV"))
							{
								StringBuilder refBuilder = new StringBuilder("");
								for(int i = 0; i<v.svlen; i++)
								{
									refBuilder.append("X");
								}
								v.ref = refBuilder.toString();
								v.alt = revComp(getVariantSeq(chr, (long)pos, (long)v.svlen, refFn));
							}
							
							// Log this SV to both logging files
							svPhaseOut.printf("%s:%d:%s\t<%s>\t%d\t%d\t|\t%d\t%d\t%d\t| %s\t%7.02f  \t|\t%s\t%d\t%d\n", 
									chr, pos, oldGenotype, v.type, v.svlen, seqLength, v.readNames.length, hap1Votes,
									hap2Votes, hap, hap1Prop, newGenotype, includeSv ? 1 : 0, overruleHomozygous ? 1 : 0);
							
							svPhaseDetailsOut.printf("%s:%d:%s\t<%s>\t%d\t%d\t|\t%d\t%d\t%d\t| %s\t%7.02f  \t|\t%s\t%d\t%d\n", 
									chr, pos, oldGenotype, v.type, v.svlen, seqLength, v.readNames.length, hap1Votes,
									hap2Votes, hap, hap1Prop, newGenotype, includeSv ? 1 : 0, overruleHomozygous ? 1 : 0);
							
							// If we include this in our phasing, mark it as phased and update the genotype
							if(includeSv && PHASE_SVS)
							{
								v.sample = newGenotype + v.sample.substring(3);
								snpData.get(chr).put(pos, v);
								reportedSvs++;
								
								typeReports.get(v.type).reported++;
								
								if(isPhased)
								{
									typeReports.get(v.type).phased++;
									phasedSvs++;
								}
								else
								{
									typeReports.get(v.type).unphased++;
									unphasedSvs++;
								}
							}
						}
						else
						{
							System.err.println("ERROR: Weird genotype call: " + v.type + " " + oldGenotype);
						}
					}
					else
					{
						System.err.println("ERROR: extreme SV length reported: " + v.type + " " + v.svlen);
						svLenErr++;
					}
				}
			}
		}
		
		svPhaseOut.close();
		svPhaseDetailsOut.close();
		
		// Output some SV phasing statistics
		System.err.printf("Reported %d, phased %d, unphased %d of %d attempted.  svlenerr: %d, genotypeerr: %d\n",
				reportedSvs, phasedSvs, unphasedSvs, snifflesCount, svLenErr, genotypeErr);
		System.err.println("type all reported phased unphased:");
		
		for(String s : typeReports.keySet())
		{
			Report r = typeReports.get(s);
			System.err.println(s + " " + r.total + " " + r.reported + " " + r.phased + " " + r.unphased);
		}
		
		System.err.println();
		
		// Print header to output file
		for(String headerLine : snpsHeader)
		{
			out.println(headerLine);
		}
		
		// Print final variants to output file
		for(String chr : snpData.keySet())
		{
			for(int pos : snpData.get(chr).keySet())
			{
				Variant v = snpData.get(chr).get(pos);
				out.println(v);
			}
		}
		
		out.close();
		
	}
	
	/*
	 * Information for each individual read
	 */
	static class Read
	{
		// The number of SVs it's included in
		int count;
		
		// The number of SNPs supporting each haplotype call
		int hap1, hap2;
		
		// A list of SNPs involved in this read
		StringBuilder snps;
		Read()
		{
			count = 0;
			hap1 = 0;
			hap2 = 0;
			snps = new StringBuilder("");
		}
	}
	
	/*
	 * Information to be reported for each type of variant
	 */
	static class Report
	{
		int total;
		int reported;
		int phased;
		int unphased;
		Report()
		{
			total = reported = phased = unphased = 0;
		}
	}
	
	/*
	 * Information for a variant
	 */
	static class Variant
	{
		String chromosome;
		int pos;
		String id, ref, alt;
		String qual, filter, info;
		String format, sample;
		String genotype;
		String[] readNames;
		String type;
		String seq;
		int svlen;
		Variant(String line)
		{
			String[] tokens = line.split("\t");
			chromosome = tokens[0];
			pos = Integer.parseInt(tokens[1]);
			id = tokens[2];
			ref = tokens[3];
			alt = tokens[4];
			qual = tokens[5];
			filter = tokens[6];
			info = tokens[7];
			format = tokens[8];
			sample = tokens.length > 9 ? tokens[9] : "";
			if(sample.length() > 0)
			{
				genotype = sample.substring(0, sample.indexOf(":"));
			}
		}
		
		/*
		 * Populates list of supporting reads from the RNAMES INFO field
		 */
		void fillSvFields()
		{
			String[] tokens = info.split(";");
			
			for(String token : tokens)
			{
				int equalsIdx = token.indexOf('=');
				if(equalsIdx != -1)
				{
					String key = token.substring(0, equalsIdx);
					String val = token.substring(1 + equalsIdx);
					if(key.equals("RNAMES"))
					{
						readNames = val.split(",");
					}
					else if(key.equals("SEQ"))
					{
						seq = val;
					}
					else if(key.equals("SVLEN"))
					{
						svlen = Math.abs(Integer.parseInt(val));
					}
					else if(key.equals("SVTYPE"))
					{
						type = val;
					}
				}
			}
			if(type == null && alt.contains("<"))
			{
				type = alt.substring(1, alt.length()-1);
			}
		}
		
		/*
		 * VCF formatted line (with no new line at the end)
		 */
		public String toString()
		{
			return chromosome + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual
					+ "\t" + filter + "\t" + info + "\t" + format + (sample.length() == 0 ? "" : ("\t" + sample));
		}
	}
}
