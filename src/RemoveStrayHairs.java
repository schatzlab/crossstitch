/*
 * Module of CrossStitch - preprocessing for ExtractHairs
 * Replaces haploid genotypes in a VCF file of SNP data with diploid genotypes
 */
import java.util.*;
import java.io.*;
public class RemoveStrayHairs {
	static int GTFIELD = 9;
public static void main(String[] args) throws IOException
{
	String fn = args[0];
	String ofn = args.length > 1 ? args[1] : (fn + ".prehairs");
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.charAt(0) == '#') out.println(line);
		else
		{
			// Variant line - check for haplotypical genotype
			String[] vcfFields = line.split("\t");
			if(vcfFields.length < 10)
			{
				// Less than ten fields so not true variant line - just print it
				out.println(line);
				continue;
			}
			String gtField = vcfFields[GTFIELD];
			int gtLength = gtField.indexOf(':');
			if(gtLength == 1)
			{
				// Genotype is single character - duplicate it
				gtField = gtField.charAt(0) + "/" + gtField;
				
				// Print out fields one at a time and replace genotype field with updated string
				for(int i = 0 ; i<vcfFields.length; i++)
				{
					out.print((i == GTFIELD ? gtField : vcfFields[i]) + (i == vcfFields.length - 1 ? "\n" : "\t"));
				}
			}
            else if(gtLength == 3 || gtField.length() == 3)
            {
                if(gtField.substring(0, 3).equals("./."))
                {
                    // Genotype is single character - duplicate it
				    gtField = "0/0" + gtField.substring(3);
				
    				// Print out fields one at a time and replace genotype field with updated string
	    			for(int i = 0 ; i<vcfFields.length; i++)
		    		{
			    		out.print((i == GTFIELD ? gtField : vcfFields[i]) + (i == vcfFields.length - 1 ? "\n" : "\t"));
				    }

                }
                else
                {
                    out.println(line);
                }
            }
			else out.println(line);
		}
	}
	out.close();
}
}
