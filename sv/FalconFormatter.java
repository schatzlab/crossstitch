import java.util.*;
import java.io.*;
public class FalconFormatter {
public static void main(String[] args) throws IOException
{
	if(args.length == 3)
	{
		String assemblyFile = args[0];
		String readsFile = args[1];
		String outputFile = args[2];
		Scanner assemblyInput = new Scanner(new FileInputStream(new File(assemblyFile)));
		String seqName = assemblyInput.next().substring(1);
		String seq = assemblyInput.next();
		Scanner readInput = new Scanner(new FileInputStream(new File(readsFile)));
		PrintWriter out = new PrintWriter(new File(outputFile));
		out.println(seqName+" "+seq);
		while(readInput.hasNext())
		{
			String readSeqName = readInput.next().substring(1);
			String readSeq = readInput.next();
			out.println(readSeqName+" "+readSeq);
		}
		out.println("+ +");
		out.println("- -");
		
		out.close();
	}
	else if(args.length == 2)
	{
		String readsFile = args[0];
		String outputFile = args[1];
		Scanner readInput = new Scanner(new FileInputStream(new File(readsFile)));
		PrintWriter out = new PrintWriter(new File(outputFile));
		ArrayList<String> names = new ArrayList<String>();
		ArrayList<String> seqs = new ArrayList<String>();
		String name = "";
		StringBuilder sb = new StringBuilder("");
		while(readInput.hasNext())
		{
			String str = readInput.nextLine();
			if(str.startsWith(">"))
			{
				if(name.length() != 0)
				{
					names.add(name);
					seqs.add(sb.toString());
					sb = new StringBuilder("");
					name = str.substring(1);
				}
				else name = str.substring(1);
				continue;
			}
			sb.append(str);
		}
		if(name.length() != 0)
		{
			names.add(name);
			seqs.add(sb.toString());
		}
		int n = names.size();
		for(int i = 0; i<n; i++)
		{
			out.println(names.get(i)+" "+seqs.get(i));
			for(int j = 0; j<n; j++)
			{
				if(j != i) out.println(names.get(j)+" "+seqs.get(j));
			}
			out.println("+ +");
		}
		out.println("- -");
		out.close();
	}
}
}
