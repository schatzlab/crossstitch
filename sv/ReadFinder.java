import java.util.*;
import java.io.*;
public class ReadFinder {
public static void main(String[] args) throws IOException
{
	if(args.length < 2)
	{
		System.out.println("Usage: java ReadFinder vcfFileName outDir");
		return;
	}
	String vcfFn = args[0];
	String outDir = args[1];
	BufferedReader vcfScanner = new BufferedReader(new InputStreamReader(new FileInputStream(new File(vcfFn))));
	String find = "<INS>";
	String line = "";
	while(true)
	{
		line = vcfScanner.readLine();
		if(line == null) break;
		if(line.contains(find))
		{
			String[] tokens = line.split("\t");
			int index = 0;
			int curPos = -1;
			String chr = tokens[0];
			curPos = Integer.parseInt(tokens[1]);
			tokens = line.split(";");
			for(String t : tokens)
			{
				if(t.startsWith("RNAMES="))
				{
				    TreeSet<String> readList = new TreeSet<String>();
					String list = t.substring(7);
					StringTokenizer str = new StringTokenizer(list, ",");
					while(str.hasMoreTokens()) readList.add(str.nextToken());
					PrintWriter out = new PrintWriter(new File(outDir + "/" + curPos + ".txt"+"."+chr));
					for(String s : readList) out.println(s);
					out.close();
				}
				index++;
			}
		}
	}
}
}
