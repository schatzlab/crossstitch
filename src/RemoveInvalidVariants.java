/*
 * Removes invalid or ignored variants from a VCF file for CrossStitch processing
 * It has two modes:
 *   The first mode removes any types other than INS/DEL/INV, adds the type to the alt field, and removes variants not on chr1-22, X, Y
 *   The second mode (-o) removes overlapping variants and filters to only variants on chr1-22, X, and Y (if the karyotype is XY).
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Scanner;
import java.util.TreeMap;

public class RemoveInvalidVariants {
	
	static void usage()
	{
		System.out.println();
		System.out.println("USAGE: java RemoveInvalidVariants [-o XY|XX] orig.vcf new.vcf");
		System.out.println();
		System.out.println(" By default, remove variants that arent on chr1-22, chrX, chrY or are not INS, DEL, INV");
		System.out.println();
		System.out.println(" In -o mode, remove any overlapping variants or variants not on chr1-22, chrX, chrY");
		System.out.println("             but DONT filter by or update any variant types");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		if(args.length == 0 || args[0].startsWith("-h") || args[0].equalsIgnoreCase("--help"))
		{
			usage();
			System.exit(1);
		}
		if(args[0].equalsIgnoreCase("-o"))
		{
			if(args.length != 4)
			{
				usage();
				System.exit(1);
			}
			String karyotype = args[1];
			boolean removeY = !karyotype.equalsIgnoreCase("XY");
			
			String fn = args[2], ofn = args[3];
			
			removeInvalidVariants(fn, ofn, true, removeY);
		}
		else
		{
			if(args.length != 2)
			{
				usage();
				System.exit(1);
			}
			
			String fn = args[0], ofn = args[1];
			removeInvalidVariants(fn, ofn, false, false);
		}
	}
	
	/*
	 * Removes invalid or ignored variants from an input VCF file and writes filtered set to an output file.
	 * removeOverlaps is whether or not to remove overlapping variants, and if it is set to true, removeY is
	 * whether or not remove variants on chromosome Y.
	 */
	static void removeInvalidVariants(String fn, String ofn, boolean removeOverlaps, boolean removeY) throws Exception
	{
		// Make a list of allowed chromosome names
		HashSet<String> goodChrNames = new HashSet<String>();
		for(int i = 1; i<=22; i++)
		{
			goodChrNames.add(i+"");
			goodChrNames.add("chr" + i);
		}
		goodChrNames.add("X");
		goodChrNames.add("chrX");
		
		if(!removeY)
		{
			goodChrNames.add("Y");
			goodChrNames.add("chrY");
		}
		
		HashSet<String> goodTypes = new HashSet<String>();
		goodTypes.add("<INS>");
		goodTypes.add("<DEL>");
		goodTypes.add("<INV>");
		
		// The number of input and output variants
		int totalVariants = 0;
		int reportedVariants = 0;
		
		// Keep track of the last two variants for overlapping
		String lastChr = "";
		
		// The last variant seen on each haplotype
		Variant last0 = new Variant(), last1 = new Variant();
		
		// Whether or not the last variant seen was homozygous, meaning last0 and last1 are the same variant
		boolean lastHomo = false;
		
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		PrintWriter errOut = new PrintWriter(System.err);
		
		// The number of reported variants of each type
		TreeMap<String, Integer> typeCounts = new TreeMap<String, Integer>();
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				out.println(line);
				continue;
			}
			
			totalVariants++;
			
			// Get VCF fields from the line
			String[] fields = line.split("\t");
			String chr = fields[0];
			int pos = Integer.parseInt(fields[1]);
			String ref = fields[3];
			String alt = fields[4];
			String filter = fields[6];
			String info = fields[7];
			String genotype = fields[9].substring(0,  3);
			
			if(removeOverlaps)
			{
				// Remove lines without the PASS FILTER value
				if(!filter.equals("PASS"))
				{
					typeCounts.put("FAIL_PASS", 1 + typeCounts.getOrDefault("FAIL_PASS", 0));
					continue;
				}
				
				// Remove variants on disallowed chromosomes
				if(!goodChrNames.contains(chr))
				{
					typeCounts.put("FAIL_CHROM", 1 + typeCounts.getOrDefault("FAIL_CHROM", 0));
					continue;
				}
				
				// Print any previous variants whose ends we've already passed
				if(lastChr.length() > 0)
				{
					// Always print if we're on a different chromosome now
					if(!lastChr.equals(chr))
					{
						if(last0.nonNull())
						{
							out.println(last0.line);
							if(lastHomo)
							{
								last1.reset();
								lastHomo = false;
							}
						}
						last0.reset();
						
						if(last1.nonNull())
						{
							out.println(last1.line);
						}
						last1.reset();
					}
					
					// Otherwise, check if this variant's start is past the end of the previous variants
					else
					{
						if(last0.nonNull() && pos > last0.pos)
						{
							out.println(last0.line);
							last0.reset();
							if(lastHomo)
							{
								last1.reset();
								lastHomo = false;
							}
						}
						if(last1.nonNull() && pos > last1.pos)
						{
							out.println(last1.line);
							last1.reset();
						}
					}
				}
				
				// Whether or not to keep this variant for printing
				boolean printVar = true;
				
				// Get the type of the variant
				String type = "SUB";
				if(ref.length() < alt.length())
				{
					type = "INS";
				}
				if(ref.length() > alt.length())
				{
					type = "DEL";
				}
				
				// Which haplotype the variant is on (0 or 1, or 2 for both)
				int hap = 2;
				if(genotype.equals("0|1"))
				{
					hap = 1;
				}
				if(genotype.equals("1|0"))
				{
					hap = 0;
				}
				
				// Whether or not this variant is a sniffles SV call
				boolean sniffles = false;
				if(info.contains("Sniffles"))
				{
					sniffles = true;
				}
				
				// Check if we already found a variant with a greater end position
				if(lastChr.length() > 0 && chr.equals(lastChr))
				{
					// Check both haplotypes for overlapping calls
					for(int hapCheck = 0; hapCheck < 2; hapCheck++)
					{
						// Make sure that we haven't already found an overlapping variant and that the current variant is on this haplotype
						if(printVar && (hap == hapCheck || hap == 2))
						{
							// The last variant on this haplotype
							Variant v = (hapCheck == 0 ? last0 : last1);
							
							// The last variant on the other haplotype
							Variant other = (hapCheck == 0 ? last1 : last0);
							if(pos <= v.pos)
							{
								// If the last call was a SNP, throw out the SNP instead
								if(sniffles && !v.sniffles)
								{
									String key = "OVERRULE_" + hap;
									typeCounts.put(key, 1 + typeCounts.getOrDefault(key, 0));
									v.reset();
									if(lastHomo)
									{
										other.reset();
									}
									errOut.printf(" overrule overlap detected %d on hap %d (lastpos: %d @ %s) pos: %d @ %s, type: %s\n",
											hapCheck + 1, hap, v.pos, lastChr, pos, chr, type);
								}
								else
								{
									// Otherwise, we ignore this call because we already have a variant that overlaps it
									printVar = false;
									String key = "FAIL_" + hap;
									typeCounts.put(key, 1 + typeCounts.getOrDefault(key, 0));
									if(sniffles)
									{
										key = "FAIL_SNIFFLES_" + hap;
										typeCounts.put(key, 1 + typeCounts.getOrDefault(key, 0));
									}
									errOut.printf(" overlap detected %d on hap %d (lastpos: %d @ %s) pos: %d @ %s, type: %s\n",
											hapCheck + 1, hap, v.pos, lastChr, pos, chr, type);
								}
							}
						}
					}
				}
				
				if(printVar)
				{
					reportedVariants++;
					typeCounts.put(type, 1 + typeCounts.getOrDefault(type, 0));
					
					lastChr = chr;
					
					// The last variants will have their end positions stored for checking overlaps with later variants
					int newPos = pos + ref.length() - 1;
					
					// Update the last variant for any haplotypes this variant appears on
					for(int hapCheck = 0; hapCheck < 2; hapCheck++)
					{
						if(hap == hapCheck || hap == 2)
						{
							Variant v = (hapCheck == 0 ? last0 : last1);
							v.pos = newPos;
							v.ref = ref;
							v.alt = alt;
							v.line = line;
							v.sniffles = sniffles;
						}
					}
					
					lastHomo = (hap == 2);
				}
			}
			else
			{
				if(info.contains("SVTYPE=INS"))
				{
					alt = "<INS>";
				}
				if(info.contains("SVTYPE=DEL"))
				{
					alt = "<DEL>";
				}
				if(goodChrNames.contains(chr) && goodTypes.contains(alt))
				{
					reportedVariants++;
					typeCounts.put(alt, 1 + typeCounts.getOrDefault(alt, 0));
					out.println(line);
				}
			}
			
		}
		
		// Print the final variant from each haplotype
		if(last0.nonNull())
		{
			out.println(last0.line);
			last0.reset();
			if(lastHomo)
			{
				last1.reset();
			}
		}
		
		// Print some statistics to STDERR
		errOut.println("## Reported " + reportedVariants + " of " + totalVariants + " variants:");
		for(String s : typeCounts.keySet())
		{
			errOut.println(s + " " + typeCounts.get(s));
		}
		
		if(last1.nonNull())
		{
			out.println(last1.line);
			last1.reset();
		}
		
		input.close();
		out.close();
		errOut.close();
	}
	
	/*
	 * The information needed to store the last variant on each haplotype for checking overlaps
	 */
	static class Variant
	{
		int pos;
		String ref, alt, line;
		boolean sniffles;
		
		Variant()
		{
			reset();
		}
		
		void reset()
		{
			pos = -1;
			ref = "";
			alt = "";
			line = "";
		}
		boolean nonNull()
		{
			return line.length() > 0;
		}
	}
}
